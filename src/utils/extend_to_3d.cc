/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/


#include "mesh.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "extend_to_3d.h"
#include "parallel.h"
#include "mesh_tools.h"
#include "boundary_info.h"

using PhysicalUnit::um;

void ExtendTo3D::operator() ()
{
  Mesh & mesh = _system.mesh();
  if(mesh.mesh_dimension()==3) return;

  this->save_old_system();

  // clear the system
  _system.clear();

  this->set_new_system();

  //setting doping/mole information and build simulation region
  this->set_new_regions();
}



void ExtendTo3D::save_old_system()
{
  Mesh & mesh = _system.mesh();

  //save 2d mesh
  magic_num    = mesh.magic_num();
  n_nodes      = mesh.n_nodes();
  n_elem       = mesh.n_elem();
  n_subs       = mesh.n_subdomains();

  // save point location
  MeshBase::node_iterator       node_it     = mesh.nodes_begin();
  const MeshBase::node_iterator node_it_end = mesh.nodes_end();
  for (; node_it != node_it_end; ++node_it)
  {
    assert (*node_it != NULL);
    assert ((*node_it)->id() == mesh_points.size());
    mesh_points.push_back(**node_it);
  }

  // save elem info
  MeshBase::element_iterator elem_it = mesh.elements_begin();
  const MeshBase::element_iterator elem_it_end = mesh.elements_end();
  for (; elem_it != elem_it_end; ++elem_it)
  {
    assert (*elem_it);
    const Elem* elem = *elem_it;
    pack_element (mesh_conn, elem);
  }

  // save subdomain info
  for(unsigned int n=0; n<mesh.n_subdomains(); n++)
  {
    subdomain_label.push_back(mesh.subdomain_label_by_id(n));
    subdomain_material.push_back(mesh.subdomain_material(n));
  }

  // get all the boundary element
  mesh.boundary_info->build_side_list(bd_elems, bd_sides, bd_ids);

  std::set<short int> boundary_ids;
  boundary_ids = mesh.boundary_info->get_boundary_ids();
  for(std::set<short int>::iterator it=boundary_ids.begin(); it!=boundary_ids.end(); ++it)
    {
      std::string label = mesh.boundary_info->get_label_by_id(*it);
      bd_map.insert(std::make_pair(label, *it));
    }


  // save solutions
  variables["T"     ].resize(_system.n_regions());
  variables["Na"    ].resize(_system.n_regions());
  variables["Nd"    ].resize(_system.n_regions());
  variables["psi"   ].resize(_system.n_regions());
  variables["n"     ].resize(_system.n_regions());
  variables["p"     ].resize(_system.n_regions());
  variables["mole_x"].resize(_system.n_regions());
  variables["mole_y"].resize(_system.n_regions());

  for(unsigned int r=0; r<_system.n_regions(); r++)
  {
    const SimulationRegion * region = _system.region(r);
    SimulationRegion::const_node_iterator it = region->nodes_begin();
    SimulationRegion::const_node_iterator it_end = region->nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * node = (*it).second;
      //if this node NOT belongs to this processor, continue
      if( !node->on_processor() ) continue;
      const FVM_NodeData * node_data = node->node_data();
      variables["T"     ][r][node->root_node()->id()] = node_data->T();
      variables["Na"    ][r][node->root_node()->id()] = node_data->Total_Na();
      variables["Nd"    ][r][node->root_node()->id()] = node_data->Total_Nd();
      variables["psi"   ][r][node->root_node()->id()] = node_data->psi();
      variables["n"     ][r][node->root_node()->id()] = node_data->n();
      variables["p"     ][r][node->root_node()->id()] = node_data->p();
      variables["mole_x"][r][node->root_node()->id()] = node_data->mole_x();
      variables["mole_y"][r][node->root_node()->id()] = node_data->mole_y();
    }
  }

  sync_solution("T");
  sync_solution("Na");
  sync_solution("Nd");
  sync_solution("psi");
  sync_solution("n");
  sync_solution("p");
  sync_solution("mole_x");
  sync_solution("mole_y");

}


void ExtendTo3D::pack_element (std::vector<int> &conn, const Elem* elem) const
{
  assert (elem != NULL);
  conn.push_back (static_cast<int>(elem->type()));
  conn.push_back (static_cast<int>(elem->subdomain_id()));
  for (unsigned int n=0; n<elem->n_nodes(); n++)
    conn.push_back (elem->node(n));
}


void ExtendTo3D::sync_solution(const std::string &sol)
{
  for(unsigned int r=0; r<_system.n_regions(); r++)
  {
    genius_assert(variables.find(sol)!=variables.end());
    variable_map & variabel = (variables[sol])[r];
    Parallel::allgather(variabel);
  }
}


void ExtendTo3D::set_new_system()
{
  Mesh & mesh = _system.mesh();

  std::map<const Node *, unsigned int> node_to_old_node_id;

  // set new mesh
  {
    mesh.magic_num() = magic_num + 2008;

    // z locations
    double zmin = _card.get_real("z.min", 0.0, "z.front")*um;
    double zmax = _card.get_real("z.max", 1.0, "z.back")*um;
    int    zdiv = _card.get_int("n.spaces", 1);

    // set new points
    for(int z=0; z<=zdiv; ++z)
    {
      double zloc = zmin + z*(zmax - zmin)/zdiv;
      for(unsigned n=0; n<n_nodes; ++n)
      {
        Node * node = mesh.add_point(Point(mesh_points[n].x(), mesh_points[n].y(), zloc));
        node_to_old_node_id[node] = n;
      }
    }

    // set new elem
    for(int z=0; z<zdiv; ++z)
    {
      unsigned int cnt = 0;

      for(unsigned n=0; n<n_elem; ++n)
      {
        const ElemType elem_type  = static_cast<ElemType>(mesh_conn[cnt++]);
        unsigned int subdomain_ID = mesh_conn[cnt++];
        Elem * elem;

        switch(elem_type)
        {
        case TRI3:
        case TRI3_FVM:
          {
            elem = Elem::build(PRISM6, NULL).release();
            elem->subdomain_id() = subdomain_ID;
            unsigned int A = mesh_conn[cnt++];
            unsigned int B = mesh_conn[cnt++];
            unsigned int C = mesh_conn[cnt++];
            elem->set_node(0) = mesh.node_ptr (A + z*n_nodes);
            elem->set_node(1) = mesh.node_ptr (B + z*n_nodes);
            elem->set_node(2) = mesh.node_ptr (C + z*n_nodes);
            elem->set_node(3) = mesh.node_ptr (A + (z+1)*n_nodes);
            elem->set_node(4) = mesh.node_ptr (B + (z+1)*n_nodes);
            elem->set_node(5) = mesh.node_ptr (C + (z+1)*n_nodes);
          }
          break;
        case QUAD4:
        case QUAD4_FVM:
          {
            elem = Elem::build(HEX8, NULL).release();
            elem->subdomain_id() = subdomain_ID;
            unsigned int A = mesh_conn[cnt++];
            unsigned int B = mesh_conn[cnt++];
            unsigned int C = mesh_conn[cnt++];
            unsigned int D = mesh_conn[cnt++];
            elem->set_node(0) = mesh.node_ptr (A + z*n_nodes);
            elem->set_node(1) = mesh.node_ptr (B + z*n_nodes);
            elem->set_node(2) = mesh.node_ptr (C + z*n_nodes);
            elem->set_node(3) = mesh.node_ptr (D + z*n_nodes);
            elem->set_node(4) = mesh.node_ptr (A + (z+1)*n_nodes);
            elem->set_node(5) = mesh.node_ptr (B + (z+1)*n_nodes);
            elem->set_node(6) = mesh.node_ptr (C + (z+1)*n_nodes);
            elem->set_node(7) = mesh.node_ptr (D + (z+1)*n_nodes);
          }
          break;
        default: genius_error();
        }
        mesh.add_elem(elem);
      }
    }

    // set new subdomains
    mesh.set_n_subdomains () = n_subs;
    for(unsigned int n=0; n<n_subs; ++n)
    {
      mesh.set_subdomain_label(n, subdomain_label[n]);
      mesh.set_subdomain_material(n, subdomain_material[n]);
    }

    // set boundary
    for(int z=0; z<zdiv; ++z)
      for(unsigned int n=0; n<bd_elems.size(); ++n)
      {
        unsigned int old_elem_id = bd_elems[n];
        unsigned int old_side    = bd_sides[n];

        unsigned int elem_id = old_elem_id + z*n_elem;
        const Elem* elem = mesh.elem(elem_id);
        unsigned int side = old_side + 1;

        mesh.boundary_info->add_side (elem, side, bd_ids[n]);
      }

    //set boundary on frond/back face
    std::map<unsigned int, short int>    boundary_id_map;
    for(unsigned int n=0; n<n_subs; ++n)
    {
      std::string bd_label = mesh.subdomain_label_by_id(n) + "_Neumann";
      if(bd_map.find(bd_label)==bd_map.end())
      {
        short int bd_id = 0;
        for(Bd_It it = bd_map.begin(); it!=bd_map.end(); ++it)
          bd_id = std::max(bd_id, it->second);
        bd_map[bd_label] = bd_id+1;
      }
      boundary_id_map[n] = bd_map[bd_label];
    }

    //frond side
    for(unsigned int n=0; n<n_elem; ++n)
    {
      const Elem* elem = mesh.elem(n);
      mesh.boundary_info->add_side (elem, 0, boundary_id_map[elem->subdomain_id()]);
    }

    // back side
    for(unsigned int n=0; n<n_elem; ++n)
    {
      const Elem* elem = mesh.elem(n + (zdiv-1)*n_elem);
      mesh.boundary_info->add_side (elem, elem->n_sides()-1, boundary_id_map[elem->subdomain_id()]);
    }

    // set boundary label
    for(Bd_It it = bd_map.begin(); it!=bd_map.end(); ++it)
      mesh.boundary_info->set_label_to_id(it->second, it->first);
  }

  // build simulation system
  _system.build_simulation_system();
  _system.sync_print_info();


  // set renumbered node id to original node id map
  {
    std::map<const Node *, unsigned int>::iterator  it=node_to_old_node_id.begin();
    for(; it!=node_to_old_node_id.end(); ++it)
      _node_to_old_node_id_map[it->first->id()] = it->second;
  }

}



void ExtendTo3D::set_new_regions()
{

  for(unsigned int r=0; r<_system.n_regions(); r++)
  {
    SimulationRegion * region = _system.region(r);

    SimulationRegion::const_node_iterator it = region->nodes_begin();
    SimulationRegion::const_node_iterator it_end = region->nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * node = (*it).second;

      //if this node NOT on local, continue
      if( !node->on_local() ) continue;

      FVM_NodeData * node_data = node->node_data();
      genius_assert(node_data);

      unsigned int old_node_id = _node_to_old_node_id_map[node->root_node()->id()];
      node_data->T()        =  variables["T"     ][r][old_node_id];
      node_data->Na()       =  variables["Na"    ][r][old_node_id];
      node_data->Nd()       =  variables["Nd"    ][r][old_node_id];
      node_data->psi()      =  variables["psi"   ][r][old_node_id];
      node_data->n()        =  variables["n"     ][r][old_node_id];
      node_data->p()        =  variables["p"     ][r][old_node_id];
      node_data->mole_x()   =  variables["mole_x"][r][old_node_id];
      node_data->mole_x()   =  variables["mole_x"][r][old_node_id];
    }
    region->reinit_after_import();
  }

}


