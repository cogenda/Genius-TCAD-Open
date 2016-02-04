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


#include "mesh_base.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "rotate_to_3d.h"
#include "mesh_communication.h"
#include "parallel.h"
#include "mesh_tools.h"
#include "boundary_info.h"

using PhysicalUnit::um;
using PhysicalUnit::nm;

void RotateTo3D::operator() ()
{
  MeshBase & mesh = _system.mesh();
  if(mesh.mesh_dimension()==3) return;

  this->save_old_system();

  // clear the system
  _system.clear();

  this->set_new_system();

  //setting doping/mole information and build simulation region
  this->set_new_regions();
}



void RotateTo3D::save_old_system()
{
  //save 2d mesh
  if(Genius::processor_id() == 0)
  {
    MeshBase & mesh = _system.mesh();

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
    
  
    for(unsigned n=0; n<n_nodes; ++n)
    {
      const Point & p = mesh_points[n];
      if( std::abs(p.x()) < 1*nm ) 
        on_axis_node.insert(n);
      else
        off_axis_node.insert(n);
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
      bool user_define = mesh.boundary_info->boundary_id_has_user_defined_label(*it);
      bd_map.insert(std::make_pair(label, std::make_pair(*it, user_define)));
    }
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
    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * node = *it;
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


void RotateTo3D::pack_element (std::vector<int> &conn, const Elem* elem) const
{
  assert (elem != NULL);
  conn.push_back (static_cast<int>(elem->type()));
  conn.push_back (static_cast<int>(elem->subdomain_id()));
  for (unsigned int n=0; n<elem->n_nodes(); n++)
    conn.push_back (elem->node(n));
}


void RotateTo3D::sync_solution(const std::string &sol)
{
  for(unsigned int r=0; r<_system.n_regions(); r++)
  {
    genius_assert(variables.find(sol)!=variables.end());
    variable_map & variabel = (variables[sol])[r];
    Parallel::allgather(variabel);
  }
}


void RotateTo3D::set_new_system()
{
  const double PI = 3.14159265358979323846;
  MeshBase & mesh = _system.mesh();

  
  std::map<const Node *, unsigned int> node_to_old_node_id;

  // set new mesh
  if(Genius::processor_id() == 0)
  {
    mesh.magic_num() = magic_num + 4096;
    
    // sectors in theta direction
    int theta_div = _card.get_int("n.sectors", 24);
    genius_assert(theta_div>=3);
    double dtheta = 2*PI/theta_div;
    
    _node_chain.resize(n_nodes);
    
    // set new points
    for(unsigned n=0; n<n_nodes; ++n)
    {
      const Point & p = mesh_points[n];
      for(int i=0; i<theta_div; i++)
      {
        double theta = i* dtheta;
        double x0 = p.x()*cos(theta);    
        double z0 = p.x()*sin(theta);
        Point new_p = Point(x0, p.y(), z0);
        
        if(on_axis_node.find(n) != on_axis_node.end())
        {
          if(i==0)
          {
            Node * node = mesh.add_point(new_p);
            node_to_old_node_id[node] = n;
            _node_chain[n].push_back(node);
          }
          else
          {
            _node_chain[n].push_back(_node_chain[n][0]);
          }
        }
        
        if(off_axis_node.find(n) != off_axis_node.end())
        {
          Node * node = mesh.add_point(new_p);
          node_to_old_node_id[node] = n;
          _node_chain[n].push_back(node);
        }
      }
      
      _node_chain[n].push_back(_node_chain[n][0]);
    }

    // set new elem
    {
      unsigned int cnt = 0;

      for(unsigned n=0; n<n_elem; ++n)
      {
        const ElemType elem_type  = static_cast<ElemType>(mesh_conn[cnt++]);
        unsigned int subdomain_ID = mesh_conn[cnt++];
        std::vector<Elem *> elems;

        switch(elem_type)
        {
            case TRI3:
            case TRI3_FVM:
            {
              unsigned int A = mesh_conn[cnt++];
              unsigned int B = mesh_conn[cnt++];
              unsigned int C = mesh_conn[cnt++];
              elems = rotate_elem(A, B, C);
            }
            break;
            case QUAD4:
            case QUAD4_FVM:
            {
              unsigned int A = mesh_conn[cnt++];
              unsigned int B = mesh_conn[cnt++];
              unsigned int C = mesh_conn[cnt++];
              unsigned int D = mesh_conn[cnt++];
              elems = rotate_elem(A, B, C, D);
            }
            break;
            default: genius_error();
        }
        for(unsigned int n=0; n<elems.size(); ++n)
        {
          Elem * elem = elems[n];
          elem->subdomain_id() = subdomain_ID;
          mesh.add_elem(elem);
        }
      }
    }

    // set new subdomains
    mesh.set_n_subdomains () = n_subs;
    for(unsigned int n=0; n<n_subs; ++n)
    {
      mesh.set_subdomain_label(n, subdomain_label[n]);
      mesh.set_subdomain_material(n, subdomain_material[n]);
    }

  }

  // broadcast mesh to all the processor
  MeshCommunication mesh_comm;
  mesh_comm.broadcast(mesh);

  // build simulation system
  _system.build_simulation_system();
  _system.sync_print_info();


  // set renumbered node id to original node id map
  if(Genius::processor_id() == 0)
  {
    std::map<const Node *, unsigned int>::iterator  it=node_to_old_node_id.begin();
    for(; it!=node_to_old_node_id.end(); ++it)
      _node_to_old_node_id_map[it->first->id()] = it->second;
  }

  Parallel::broadcast(_node_to_old_node_id_map);

}



std::vector<Elem *> RotateTo3D::rotate_elem(unsigned int A, unsigned int B, unsigned int C)
{
  bool A_on_axis = (on_axis_node.find(A) != on_axis_node.end());
  bool B_on_axis = (on_axis_node.find(B) != on_axis_node.end());
  bool C_on_axis = (on_axis_node.find(C) != on_axis_node.end());
  
  unsigned int AA=A,BB=B,CC=C;
  
  unsigned int n_on_axis = (unsigned int)A_on_axis + (unsigned int)B_on_axis + (unsigned int)C_on_axis;
  if( n_on_axis == 1 )
  {
    if(A_on_axis)
    { AA=A; BB=B; CC=C; }
    if(B_on_axis)
    { AA=B; BB=C; CC=A; }
    if(C_on_axis)
    { AA=C; BB=A; CC=B; }
  }
  if( n_on_axis == 2 )
  {
    if(A_on_axis && B_on_axis)
    { AA=A; BB=B; CC=C; }
    if(B_on_axis && C_on_axis)
    { AA=B; BB=C; CC=A; }
    if(C_on_axis && A_on_axis)
    { AA=C; BB=A; CC=B; }
  }
  
  const std::vector<Node *> & A_chian =  _node_chain[AA];
  const std::vector<Node *> & B_chian =  _node_chain[BB];
  const std::vector<Node *> & C_chian =  _node_chain[CC];
  
  genius_assert( A_chian.size() == B_chian.size() );
  genius_assert( A_chian.size() == C_chian.size() );
  
  std::vector<Elem *> res;
  
  if(n_on_axis == 0)
  {
    for(unsigned int i=0; i<A_chian.size()-1; i++)
    {
      Elem * elem = Elem::build(PRISM6, NULL).release();
      elem->set_node(3) = A_chian[i+1];
      elem->set_node(4) = B_chian[i+1];
      elem->set_node(5) = C_chian[i+1];
      elem->set_node(0) = A_chian[i];
      elem->set_node(1) = B_chian[i];
      elem->set_node(2) = C_chian[i];
      
      res.push_back(elem);
    }
  }
  
  if(n_on_axis == 1)
  {
    for(unsigned int i=0; i<A_chian.size()-1; i++)
    {
      Elem * elem = Elem::build(PYRAMID5, NULL).release();
      elem->set_node(4) = A_chian[i];
      elem->set_node(2) = B_chian[i];
      elem->set_node(1) = C_chian[i];
      elem->set_node(3) = B_chian[i+1];
      elem->set_node(0) = C_chian[i+1];
      res.push_back(elem);
    }
  }
  
  if(n_on_axis == 2)
  {
    for(unsigned int i=0; i<A_chian.size()-1; i++)
    {
      Elem * elem = Elem::build(TET4, NULL).release();
      elem->set_node(0) = A_chian[i];
      elem->set_node(1) = B_chian[i];
      elem->set_node(2) = C_chian[i];
      elem->set_node(3) = C_chian[i+1];
      res.push_back(elem);
    }
  }
  
  return res;
}


#if 1
std::vector<Elem *> RotateTo3D::rotate_elem(unsigned int A, unsigned int B, unsigned int C, unsigned int D)
{
  bool A_on_axis = (on_axis_node.find(A) != on_axis_node.end());
  bool B_on_axis = (on_axis_node.find(B) != on_axis_node.end());
  bool C_on_axis = (on_axis_node.find(C) != on_axis_node.end());
  bool D_on_axis = (on_axis_node.find(D) != on_axis_node.end());
  
  unsigned int AA=A,BB=B,CC=C,DD=D;
  
  unsigned int n_on_axis = (unsigned int)A_on_axis + (unsigned int)B_on_axis + (unsigned int)C_on_axis + (unsigned int)D_on_axis;
  genius_assert( n_on_axis==0 || n_on_axis==2 );

  if( n_on_axis == 2 )
  {
    if(A_on_axis && B_on_axis)
    { AA=A; BB=B; CC=C; DD=D;}
    if(B_on_axis && C_on_axis)
    { AA=B; BB=C; CC=D; DD=A;}
    if(C_on_axis && D_on_axis)
    { AA=C; BB=D; CC=A; DD=B;}
    if(D_on_axis && A_on_axis)
    { AA=D; BB=A; CC=B; DD=C;}
  }
  
  const std::vector<Node *> & A_chian =  _node_chain[AA];
  const std::vector<Node *> & B_chian =  _node_chain[BB];
  const std::vector<Node *> & C_chian =  _node_chain[CC];
  const std::vector<Node *> & D_chian =  _node_chain[DD];
  
  genius_assert( A_chian.size() == B_chian.size() );
  genius_assert( A_chian.size() == C_chian.size() );
  genius_assert( A_chian.size() == D_chian.size() );
  
  std::vector<Elem *> res;
  
  if(n_on_axis == 0)
  {
    for(unsigned int i=0; i<A_chian.size()-1; i++)
    {
      Elem * elem = Elem::build(HEX8, NULL).release();
      elem->set_node(0) = A_chian[i];
      elem->set_node(1) = B_chian[i];
      elem->set_node(2) = C_chian[i];
      elem->set_node(3) = D_chian[i];
      elem->set_node(4) = A_chian[i+1];
      elem->set_node(5) = B_chian[i+1];
      elem->set_node(6) = C_chian[i+1];
      elem->set_node(7) = D_chian[i+1];
      res.push_back(elem);
    }
  }
  
  
  if(n_on_axis == 2)
  {
    for(unsigned int i=0; i<A_chian.size()-1; i++)
    {
      Elem * elem = Elem::build(PRISM6, NULL).release();
      elem->set_node(5) = A_chian[i];
      elem->set_node(2) = B_chian[i];
      elem->set_node(1) = C_chian[i];
      elem->set_node(4) = D_chian[i];
      elem->set_node(0) = C_chian[i+1];
      elem->set_node(3) = D_chian[i+1];
      res.push_back(elem);
    }
  }
  
  return res;
}
#endif


void RotateTo3D::set_new_regions()
{

  for(unsigned int r=0; r<_system.n_regions(); r++)
  {
    SimulationRegion * region = _system.region(r);

    SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
    SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = (*node_it);

      FVM_NodeData * node_data = fvm_node->node_data();
      genius_assert(node_data);

      unsigned int old_node_id = _node_to_old_node_id_map[fvm_node->root_node()->id()];
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



