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

// C++ includes
#include <numeric>

// cgns lib include
#include <cgnslib.h>

// Local includes
#include "cgns_io.h"
#include "mesh_base.h"
#include "mesh.h"
#include "boundary_info.h"
#include "parallel.h"
#include "simulation_system.h"
#include "boundary_condition_collector.h"
#include "mesh_communication.h"

using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::V;
using PhysicalUnit::K;
using PhysicalUnit::eV;


inline std::vector<double> sort_it (const std::vector<double> & x, const std::vector<unsigned int > &id)
{
  std::vector<double> xx(id.size());

  for(unsigned int i=0; i<id.size(); i++)
    xx[id[i]] = x[i];

  return xx;
}


// ------------------------------------------------------------
// CGNSIO class members
// ------------------------------------------------------------


void CGNSIO::read (const std::string& filename)
{
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  Mesh & mesh = system.mesh();

  // clear the system
  system.clear();

  // map the global id stored in cgns file to node *
  // since Node * can be only used on local processor,
  // we will change it to unique id() for all the processor later.
  std::map<int, Node*>   global_id_to_node;

  if( Genius::processor_id() == 0)
  {
    // open CGNS file for read
    genius_assert(!cg_open(filename.c_str(), MODE_READ, &fn));

    // get the base number
    genius_assert(!cg_nbases(fn, &B));

    //cgns file should have only one base
    genius_assert(B==1);

    int  cell_dim, physical_dim;
    char basename[32];
    //read base information
    genius_assert(!cg_base_read(fn, B, basename, &cell_dim, &physical_dim));
    genius_assert(physical_dim==3);

    // get mesh magic number from base name
    {
      std::string basename_string = basename;
      //since the basename has the format of "GENIUS_Mesh_magicnum" we first get the substr of magicnum
      std::string basename_prefix = "GENIUS_Mesh_";
      std::string magic_num_string = basename_string.substr(basename_prefix.length());

      // translate string to int
      std::stringstream   ss;
      ss << magic_num_string;
      ss >> mesh.magic_num();
    }

    //cgns file may have several zones
    genius_assert(!cg_nzones(fn, B, &Z));
    // tell mesh structure how many zones
    mesh.set_n_subdomains() = Z;

    // search for each zone
    for(int z_id=1; z_id<=Z; z_id++)
    {
      // temporary variable
      int            isize[3];

      // assert the zone type is Unstructured
      {
        ZoneType_t     zonetype;
        genius_assert(!cg_zone_type(fn, B, z_id, &zonetype));
        genius_assert(zonetype==Unstructured);
      }

      // zone basic information
      {
        char           zone_name[32];
        genius_assert(!cg_zone_read(fn, B, z_id, zone_name, isize));
        mesh.set_subdomain_label(z_id -1, zone_name );
      }

      // read material of this zone
      {
        int            ndescriptor;
        genius_assert(!cg_goto(fn , B, "Zone_t", z_id, "end"));
        genius_assert(!cg_ndescriptors(&ndescriptor));
        if(ndescriptor)
        {
          char descriptorname[32];
          char *material;
          genius_assert(!cg_descriptor_read(ndescriptor, descriptorname, &material));
          mesh.set_subdomain_material(z_id -1, material);
          free(material);
        }
      }

      //alloc double array for x,y,z coordinate
      std::vector<double> x(isize[0]);
      std::vector<double> y(isize[0]);
      std::vector<double> z(isize[0]);

      //read the coordinate
      {
        int  range_min = 1;
        genius_assert(!cg_coord_read(fn, B, z_id, "CoordinateX", RealDouble, &range_min, &isize[0], &x[0]));
        genius_assert(!cg_coord_read(fn, B, z_id, "CoordinateY", RealDouble, &range_min, &isize[0], &y[0]));
        genius_assert(!cg_coord_read(fn, B, z_id, "CoordinateZ", RealDouble, &range_min, &isize[0], &z[0]));
      }

      // read the global index of region node in the mesh
      std::vector<int> global_node_id;
      {
        genius_assert(!cg_goto(fn, B, "Zone_t", z_id , "UserDefinedData_t", 1, "end"));
        int A;
        genius_assert(!cg_narrays( &A ));
        genius_assert(A==1);
        char       ArrayName[32];
        DataType_t DataType;
        int        DataDimension;
        int        DimensionVector;
        genius_assert(!cg_array_info(A, ArrayName , &DataType , &DataDimension , &DimensionVector ));
        genius_assert(DataType == Integer);
        genius_assert(DimensionVector==isize[0]);

        global_node_id.resize(DimensionVector);
        genius_assert(!cg_array_read_as(A, DataType, &global_node_id[0] ));
      }
      region_global_id.push_back(global_node_id);

      // fill point location into mesh
      for(int i=0; i<isize[0]; i++)
        if( global_id_to_node.find(global_node_id[i]) == global_id_to_node.end() )
        {
          Node * node = mesh.add_point(Point(x[i]*cm, y[i]*cm, z[i]*cm), global_node_id[i]);
          genius_assert( node->id() == static_cast<unsigned int>(global_node_id[i]) );
          global_id_to_node.insert(std::pair<int, Node*>(global_node_id[i], node) );
        }

      // find out how many sections
      genius_assert(!cg_nsections(fn, B, z_id, &S));
      // only have one section here
      genius_assert(S==1);

      //read section information
      {
        char           sectionname[32];
        ElementType_t  elemtype;
        int            range_min = 1;
        int            nbndry;
        int            iparent_flag;
        genius_assert(!cg_section_read(fn, B, z_id, S, sectionname, &elemtype, &range_min, &isize[1], &nbndry, &iparent_flag));
        genius_assert(elemtype==MIXED);
      }

      //read elements
      std::vector<Elem *>  elem_list;
      {
        int ElementDataSize;
        genius_assert(!cg_ElementDataSize(fn, B, z_id, S, &ElementDataSize));
        std::vector<int> elem(ElementDataSize);

        int  iparentdata;
        genius_assert(!cg_elements_read(fn, B, z_id, S, &elem[0], &iparentdata));
#ifdef ENABLE_AMR
        // read the element id and attribute
        std::vector<int> elem_id;
        std::vector<int> elem_attribute;
        {
          genius_assert(!cg_goto(fn, B, "Zone_t", z_id , "UserDefinedData_t", 2, "end"));
          int A;
          genius_assert(!cg_narrays( &A ));
          genius_assert(A==2);
          char       ArrayName[32];
          DataType_t DataType;
          int        DataDimension;
          int        DimensionVector;

          genius_assert(!cg_array_info(1, ArrayName , &DataType , &DataDimension , &DimensionVector ));
          genius_assert(DataType == Integer);
          elem_id.resize(DimensionVector);
          genius_assert(!cg_array_read_as(1, DataType, &elem_id[0] ));

          genius_assert(!cg_array_info(2, ArrayName , &DataType , &DataDimension , &DimensionVector ));
          genius_assert(DataType == Integer);
          elem_attribute.resize(DimensionVector);
          genius_assert(!cg_array_read_as(2, DataType, &elem_attribute[0] ));
        }
#endif
        std::vector<int>::iterator it = elem.begin();
        for(; it != elem.end(); )
          switch(*it++)
          {
            // here TRI_3, QUAD_4, TETRA_4, PENTA_6 and HEXA_8 are cgns defined element types!
          case  TRI_3:
            {
              Elem* elem = mesh.add_elem(Elem::build(TRI3).release());
              elem->set_node(0) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(1) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(2) = global_id_to_node[global_node_id[*it++ -1]];
              elem->subdomain_id() = z_id-1;
              elem_list.push_back(elem);
              break;
            }
          case  QUAD_4:
            {
              Elem* elem = mesh.add_elem(Elem::build(QUAD4).release());
              elem->set_node(0) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(1) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(2) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(3) = global_id_to_node[global_node_id[*it++ -1]];
              elem->subdomain_id() = z_id-1;
              elem_list.push_back(elem);
              break;
            }
          case  TETRA_4:
            {
              Elem* elem = mesh.add_elem(Elem::build(TET4).release());
              elem->set_node(0) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(1) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(2) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(3) = global_id_to_node[global_node_id[*it++ -1]];
              elem->subdomain_id() = z_id-1;
              elem_list.push_back(elem);
              break;
            }
          case  PYRA_5 :
          {
              Elem* elem = mesh.add_elem(Elem::build(PYRAMID5).release());
              elem->set_node(0) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(1) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(2) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(3) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(4) = global_id_to_node[global_node_id[*it++ -1]];
              elem->subdomain_id() = z_id-1;
              elem_list.push_back(elem);
              break;
          }
          case  PENTA_6:
            {
              Elem* elem = mesh.add_elem(Elem::build(PRISM6).release());
              elem->set_node(0) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(1) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(2) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(3) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(4) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(5) = global_id_to_node[global_node_id[*it++ -1]];
              elem->subdomain_id() = z_id-1;
              elem_list.push_back(elem);
              break;
            }
          case  HEXA_8:
            {
              Elem* elem = mesh.add_elem(Elem::build(HEX8).release());
              elem->set_node(0) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(1) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(2) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(3) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(4) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(5) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(6) = global_id_to_node[global_node_id[*it++ -1]];
              elem->set_node(7) = global_id_to_node[global_node_id[*it++ -1]];
              elem->subdomain_id() = z_id-1;
              elem_list.push_back(elem);
              break;
            }
          default :
            {
              MESSAGE << "ERROR: Unsupported element type found in this CGNS file!"<<std::endl; RECORD();
              genius_error();
            }

          }

#ifdef ENABLE_AMR
        // set element attribute if AMR is ensbled
        genius_assert(elem_id.size() == elem_list.size());
        genius_assert(elem_attribute.size() == 6*elem_list.size());

        std::map<int, Elem*> parents;
        for(unsigned int n=0; n<elem_list.size(); ++n)
          parents.insert(std::make_pair(elem_id[n], elem_list[n]));

        for(unsigned int n=0; n<elem_list.size(); ++n)
        {
          Elem* elem = elem_list[n];

          const int parent_ID         = elem_attribute[6*n+0];
          const int which_child       = elem_attribute[6*n+1];

          const int level             = elem_attribute[6*n+2];
          const int p_level           = elem_attribute[6*n+3];
          const Elem::RefinementState refinement_flag   = static_cast<Elem::RefinementState>(elem_attribute[6*n+4]);
          const Elem::RefinementState p_refinement_flag = static_cast<Elem::RefinementState>(elem_attribute[6*n+5]);

          if (parent_ID != -1) // Do a log(n) search for the parent
          {
            Elem* my_parent = parents.count(parent_ID) ? parents[parent_ID] : NULL;

            // If the parent was not previously added, we cannot continue.
            if (my_parent == NULL)
            {
              std::cerr << "Parent element with ID " << parent_ID
              << " not found." << std::endl;
              genius_error();
            }

            assert (my_parent->refinement_flag() == Elem::INACTIVE);

            my_parent->add_child(elem, which_child);
            elem->set_parent(my_parent);

            assert (my_parent->type() == elem->type());
            assert (my_parent->child(which_child) == elem);
          }

          else // level 0 element has no parent
          {
            assert (level == 0);
          }

          // Assign the IDs
          assert (elem->level() == static_cast<unsigned int>(level));
          elem->set_refinement_flag(refinement_flag);
          elem->set_p_refinement_flag(p_refinement_flag);
          elem->set_p_level(p_level);

        }
#endif

      }


      // read boundarys
      genius_assert(!cg_nbocos(fn, B, z_id, &BC));
      for(int bc_id=1; bc_id<=BC; bc_id++)
      {
        char           boconame[32];
        BCType_t       bocotype ;
        PointSetType_t ptset_type ;
        int            npnts ;
        int            normalindex ;
        int            normallistflag;
        DataType_t     normaldatatype;
        int            ndataset;
        int            nuserdata;

        genius_assert(!cg_boco_info(fn, B, z_id, bc_id, boconame, &bocotype, &ptset_type, &npnts, &normalindex, &normallistflag, &normaldatatype, &ndataset) );
        genius_assert(ptset_type==PointList);
        // for old cgns file, the bc_label is boconame
        // however, new cgns file will save bc label in descriptor, thus support long (>32 chars) bc label
        std::string bc_label(boconame);
        // new cgns file will contain string for boundary condition as well
        std::string bc_settings;
        {
          genius_assert(!cg_goto(fn, B, "Zone_t", z_id, "ZoneBC_t", 1, "BC_t", bc_id, "end"));
          int ndescriptor;
          genius_assert(!cg_ndescriptors(&ndescriptor));
          if(ndescriptor)
          {
            assert(ndescriptor == 2);

            // read first descriptor -- bc label
            char descriptorname[32];
            char *_bc_label;
            assert(!cg_descriptor_read(1, descriptorname, &_bc_label));
            assert(std::string(descriptorname) == "bc_label");
            bc_label = _bc_label;
            free(_bc_label);

            // read boundary condition
            char *_bc_settings;
            assert(!cg_descriptor_read(2, descriptorname, &_bc_settings));
            assert(std::string(descriptorname) == "bc_settings");
            bc_settings = _bc_settings;
            free(_bc_settings);
          }
        }

        // get the number of user defined data structure
        genius_assert(!cg_goto(fn, B, "Zone_t", z_id, "ZoneBC_t", 1, "BC_t", bc_id, "end"));
        genius_assert(!cg_nuser_data(&nuserdata));

        genius_assert(!cg_goto(fn, B, "Zone_t", z_id, "ZoneBC_t", 1, "BC_t", bc_id, "UserDefinedData_t", 1, "end"));

        int A;
        genius_assert(!cg_narrays( &A ));
        genius_assert(A==2);

        std::vector<int> bd_elem;
        std::vector<int> bd_side;
        // read the element-side array
        {
          char       ArrayName[32];
          DataType_t DataType;
          int        DataDimension;
          int        DimensionVector;

          genius_assert(!cg_array_info(1, ArrayName , &DataType , &DataDimension , &DimensionVector ));
          genius_assert(DataType == Integer);
          bd_elem.resize(DimensionVector);
          genius_assert(!cg_array_read_as(1, Integer, &bd_elem[0] ));

          genius_assert(!cg_array_info(2, ArrayName , &DataType , &DataDimension , &DimensionVector ));
          genius_assert(DataType == Integer);
          bd_side.resize(DimensionVector);
          genius_assert(!cg_array_read_as(2, Integer, &bd_side[0] ));
        }

        short int bd_id = mesh.boundary_info->get_id_by_label(bc_label);
        if(bd_id == BoundaryInfo::invalid_id)
          bd_id = mesh.boundary_info->n_boundary_ids() + 1;

        for(unsigned int n=0; n <bd_side.size(); n++)
          mesh.boundary_info->add_side(elem_list[bd_elem[n]-1], bd_side[n], bd_id);

        mesh.boundary_info->set_label_to_id( bd_id, bc_label );
        if(!bc_settings.empty())
          mesh.boundary_info->set_description_to_id(bd_id, bc_settings);

        // read extra boundary data: potential for electrode
        if( nuserdata == 2)
        {
          double potential=0;
          genius_assert(!cg_goto(fn, B, "Zone_t", z_id, "ZoneBC_t", 1, "BC_t", bc_id, "UserDefinedData_t", 2, "end"));
          genius_assert(!cg_array_read_as(1,RealDouble, &potential));
          electrode_potential[bd_id] = potential;
        }
      }

      //read all the solutions

      // the solution name to solution value for this region map
      std::map<std::pair<std::string,std::string> , std::vector<double> > sol_map;

      //how many fieldsol node?
      genius_assert(!cg_nsols(fn, B, z_id, &SOL));
      for(int sol_id=1; sol_id<=SOL; sol_id++)
      {
        genius_assert(!cg_nfields(fn, B, z_id, sol_id, &F));
        char           solutionname[32];
        GridLocation_t locationtype;
        genius_assert(!cg_sol_info(fn,B,z_id,sol_id,solutionname, &locationtype));
        for(int f_id=1; f_id<=F; f_id++)
        {
          int        imin=1;
          DataType_t datatype;
          char       fieldname[32];
          cg_field_info(fn, B, z_id, sol_id, f_id, &datatype, fieldname);
          genius_assert(datatype==RealDouble);

          std::pair<std::string,std::string> key(solutionname,fieldname);
          sol_map[key].resize(isize[0]);
          genius_assert(!cg_field_read(fn, B, z_id, sol_id, fieldname, datatype, &imin, &isize[0], &((sol_map[key])[0])));
        }
      }
      // save solution name to solution value map for later usage
      region_solutions.push_back(sol_map);

    }// search for each zone

    cg_close(fn);

  } //if( Genius::processor_id() == 0)



  // broadcast mesh to all the processor
  MeshCommunication mesh_comm;
  mesh_comm.broadcast(mesh);


  // build simulation system
  system.build_simulation_system();
  system.sync_print_info();

  // after mesh setup, the node id() may be changed. we get the new node id() by Node * we record before
  // this is done only for processor 0.
  if( Genius::processor_id() == 0)
  {
    for(unsigned int r=0; r<region_solutions.size(); r++)
    {
      std::map<int, unsigned int > global_id_to_node_id;
      std::vector<int>::iterator it = region_global_id[r].begin();
      for(; it!=region_global_id[r].end(); ++it )
        global_id_to_node_id[*it] = global_id_to_node[*it]->id();

      region_global_id_to_node_id.push_back(global_id_to_node_id);
    }
  }

  //distribute (solution's) global id to all the processor
  if(Genius::processor_id() != 0)
    region_global_id.resize(system.n_regions());
  for(unsigned int r=0; r<system.n_regions(); r++)
    Parallel::broadcast(region_global_id[r] , 0);


  //distribute (solution's) global id to node id information to all the processor
  if(Genius::processor_id() != 0)
    region_global_id_to_node_id.resize(system.n_regions());
  for(unsigned int r=0; r<system.n_regions(); r++)
    Parallel::broadcast(region_global_id_to_node_id[r] , 0);

  // distribute solution data to all the processor
  if(Genius::processor_id() != 0)
    region_solutions.resize(system.n_regions());

  for(unsigned int r=0; r<region_solutions.size(); r++)
  {
    unsigned int n_solutions = region_solutions[r].size();
    Parallel::broadcast(n_solutions , 0);

    std::string sol_name;
    std::string field_name;
    std::vector<double> sol_array;

    std::map<std::pair<std::string,std::string>, std::vector<double> >::iterator it = region_solutions[r].begin();

    for(unsigned int n=0; n<n_solutions; n++)
    {
      if(Genius::processor_id() == 0)
      {
        sol_name  = (*it).first.first;
        field_name  = (*it).first.second;
        sol_array = (*it).second;
        ++it;
      }

      Parallel::broadcast(sol_name   , 0);
      Parallel::broadcast(field_name , 0);
      Parallel::broadcast(sol_array  , 0);

      if(Genius::processor_id() != 0)
        (region_solutions[r])[std::pair<std::string,std::string>(sol_name,field_name)] = sol_array;
    }
  }

  // now all the processor have enough information for
  // updating region solution data
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);

    std::vector<int> & global_id = region_global_id[r];
    std::map<int, unsigned int > & global_id_to_node_id = region_global_id_to_node_id[r];
    std::map<std::pair<std::string,std::string>, std::vector<double> > & solution = region_solutions[r];

    switch ( region->type() )
    {
    case SemiconductorRegion :
      {

        //bool sigle   = Material::IsSingleCompSemiconductor(region->material());
        //bool complex = Material::IsComplexCompSemiconductor(region->material());

        for(unsigned int n=0; n<global_id.size(); n++)
        {
          unsigned int node_id = global_id_to_node_id[global_id[n]];
          FVM_Node * fvm_node = region->region_fvm_node(node_id);  genius_assert(fvm_node);
          if( !fvm_node->root_node()->on_local() ) continue;

          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

          for(std::map< std::pair<std::string,std::string>, std::vector<double> >::const_iterator var_it = solution.begin();
              var_it!=solution.end(); var_it++)
          {
            if(var_it->first.first == "Solution")
            {
              // field solution
              if(var_it->first.second == "elec_density")
              {
                node_data->n() = var_it->second[n] * pow(cm,-3);
                continue;
              }
              if(var_it->first.second == "hole_density")
              {
                node_data->p() = var_it->second[n] * pow(cm,-3);
                continue;
              }
              if(var_it->first.second == "potential")
              {
                node_data->psi() = var_it->second[n] * V;
                continue;
              }
              if(var_it->first.second == "lattice_temperature")
              {
                node_data->T() = var_it->second[n] * K;
                continue;
              }
              if(var_it->first.second == "elec_temperature")
              {
                node_data->Tn() = var_it->second[n] * K;
                continue;
              }
              if(var_it->first.second == "hole_temperature")
              {
                node_data->Tp() = var_it->second[n] * K;
                continue;
              }
              if(var_it->first.second == "elec_quantum_potential")
              {
                node_data->Eqc() = var_it->second[n] * eV;
                continue;
              }
              if(var_it->first.second == "hole_quantum_potential")
              {
                node_data->Eqv() = var_it->second[n] * eV;
                continue;
              }
            }
            if(var_it->first.first == "Doping")
            {
              // doping
              if(var_it->first.second == "Na")
              {
                node_data->Na()  = var_it->second[n] * pow(cm, -3);
                continue;
              }
              if(var_it->first.second == "Nd")
              {
                node_data->Nd()  = var_it->second[n] * pow(cm, -3);
                continue;
              }
              if(var_it->first.second == "P")
              {
                node_data->P()   = var_it->second[n] * pow(cm, -3);
                continue;
              }
              if(var_it->first.second == "As")
              {
                node_data->As()  = var_it->second[n] * pow(cm, -3);
                continue;
              }
              if(var_it->first.second == "Sb")
              {
                node_data->Sb()  = var_it->second[n] * pow(cm, -3);
                continue;
              }
              if(var_it->first.second == "B")
              {
                node_data->B()   = var_it->second[n] * pow(cm, -3);
                continue;
              }
            }
            if(var_it->first.first == "Mole")
            {
              // mole fraction
              if(var_it->first.second =="mole_x")
              {
                node_data->mole_x() = var_it->second[n];
                genius_assert(node_data->mole_x() <= 1.0);
                continue;
              }
              if(var_it->first.second =="mole_y")
              {
                node_data->mole_y() = var_it->second[n];
                genius_assert(node_data->mole_y() <= 1.0);
                continue;
              }
            }
            if(var_it->first.first == "Custom")
            {
              node_data->CreateUserScalarValue(var_it->first.second);
              node_data->UserScalarValue(var_it->first.second) = var_it->second[n];
            }
          }
        }

        // after import previous solutions, we re-init region here
        region->reinit_after_import();

        break;
      }
    case InsulatorRegion     :
      {
        for(unsigned int n=0; n<global_id.size(); n++)
        {
          unsigned int node_id = global_id_to_node_id[global_id[n]];
          FVM_Node * fvm_node = region->region_fvm_node(node_id);  genius_assert(fvm_node);
          if( !fvm_node->root_node()->on_local() ) continue;

          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

          for(std::map< std::pair<std::string,std::string>, std::vector<double> >::const_iterator var_it = solution.begin();
              var_it!=solution.end(); var_it++)
          {
            if(var_it->first.first == "Solution")
            {
              // field solution
              if(var_it->first.second == "potential")
              {
                node_data->psi() = var_it->second[n] * V;
                continue;
              }
              if(var_it->first.second == "lattice_temperature")
              {
                node_data->T() = var_it->second[n] * K;
                continue;
              }
            }
            if(var_it->first.first == "Custom")
            {
              node_data->CreateUserScalarValue(var_it->first.second);
              node_data->UserScalarValue(var_it->first.second) = var_it->second[n];
            }
          }
        }

        // after import previous solutions, we re-init region here
        region->reinit_after_import();

        break;
      }
    case ConductorRegion     :
      {
        for(unsigned int n=0; n<global_id.size(); n++)
        {
          unsigned int node_id = global_id_to_node_id[global_id[n]];
          FVM_Node * fvm_node = region->region_fvm_node(node_id);  genius_assert(fvm_node);
          if( !fvm_node->root_node()->on_local() ) continue;

          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

          for(std::map< std::pair<std::string,std::string>, std::vector<double> >::const_iterator var_it = solution.begin();
              var_it!=solution.end(); var_it++)
          {
            if(var_it->first.first == "Solution")
            {
              // field solution
              if(var_it->first.second == "potential")
              {
                node_data->psi() = var_it->second[n] * V;
                continue;
              }
              if(var_it->first.second == "lattice_temperature")
              {
                node_data->T() = var_it->second[n] * K;
                continue;
              }
            }
            if(var_it->first.first == "Custom")
            {
              node_data->CreateUserScalarValue(var_it->first.second);
              node_data->UserScalarValue(var_it->first.second) = var_it->second[n];
            }
          }
        }

        // after import previous solutions, we re-init region here
        region->reinit_after_import();

        break;
      }

    // no solution data in vacuum region?
    case VacuumRegion     :
      {
        region->reinit_after_import();
        break;
      }
    case PMLRegion     :
      {
        region->reinit_after_import();
        break;
      }
    default:
      {
        MESSAGE<<"ERROR: Unsupported region type found during CGNS import."<<std::endl; RECORD();
        genius_error();
      }
    }
  }

  // set potential of electrode
  Parallel::broadcast(electrode_potential, 0);
  std::map< short int, double >::iterator el_it = electrode_potential.begin();
  for(; el_it!=electrode_potential.end(); ++el_it )
  {
    std::string bc_label= mesh.boundary_info->get_label_by_id((*el_it).first);
    BoundaryCondition * bc = system.get_bcs()->get_bc(bc_label);
    if(bc && bc->is_electrode())
      bc->ext_circuit()->potential() = (*el_it).second;
  }
}




void CGNSIO::write (const std::string& filename)
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  const Mesh & mesh = system.mesh();

  if( Genius::processor_id() == 0)
  {
    // remove old file if exist
    remove(filename.c_str());

    // open CGNS file for write
    genius_assert(!cg_open(filename.c_str(), MODE_WRITE, &fn));

    // create base of three dimensional mesh, here mesh takes its magic number as postfix
    std::string cgns_mesh_label;
    std::stringstream   ss;
    ss << "GENIUS_Mesh_" << system.mesh().magic_num();
    ss >> cgns_mesh_label;

    genius_assert(!cg_base_write(fn, cgns_mesh_label.c_str(), 3, 3, &B));
  }

  // create zone
  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);


    int size[3];
    // region node number
    size[0] = region->n_node();

    //region cell number, since region only has local cell, we need to gather from all the processor
    size[1] = 0;

    //boundary cell number
    size[2] = 0;

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    x.reserve(size[0]);
    y.reserve(size[0]);
    z.reserve(size[0]);

    // addtional information about global id of region node
    std::vector<int> global_node_id;
    global_node_id.reserve(size[0]);

    if( Genius::processor_id() == 0)
    {
      // write region name
      genius_assert(!cg_zone_write(fn, B, region->name().c_str(), size, Unstructured, &Z));

      // goto the current region
      genius_assert(!cg_goto(fn,B,"Zone_t",Z,"end"));

      // write down material
      genius_assert(!cg_descriptor_write("Material", region->material().c_str() ));
    }

    std::map<const Node *, unsigned int> node_pointer_to_id;
    unsigned int local_id = 1;

    // set coordinates of each node in this zone
    SimulationRegion::const_node_iterator node_it = region->nodes_begin();
    for(; node_it!=region->nodes_end(); ++node_it, ++local_id)
    {
      //convert the unit to cm
      const Node * node = (*node_it).second->root_node();
      x.push_back( (*node)(0)/cm );
      y.push_back( (*node)(1)/cm );
      z.push_back( (*node)(2)/cm );
      // save global index of the region node
      global_node_id.push_back(node->id());
      // insert Node * to local node index map
      node_pointer_to_id.insert( std::pair<const Node *, unsigned int>(node, local_id) );
    }

    if( Genius::processor_id() == 0)
    {
      // write down coordinates
      cg_coord_write(fn, B, Z, RealDouble, "CoordinateX", &x[0], &C);
      cg_coord_write(fn, B, Z, RealDouble, "CoordinateY", &y[0], &C);
      cg_coord_write(fn, B, Z, RealDouble, "CoordinateZ", &z[0], &C);

      // write the global index of region node into cgns file
      genius_assert(!cg_goto(fn, B, "Zone_t", Z , "end"));
      genius_assert(!cg_user_data_write ("Global_Node_Index"));
      genius_assert(!cg_goto(fn, B, "Zone_t", Z, "UserDefinedData_t", 1, "end"));
      genius_assert(!cg_array_write("index", Integer, 1, &size[0], &global_node_id[0]));
    }


    // set element connectivity here

    std::vector<int> elem;           // for element connectivity
    std::vector<int> elem_attribute; // for element attribute
    std::vector<unsigned int> elem_id;

    Mesh::const_element_iterator elem_it = mesh.subdomain_elements_begin(r);
    Mesh::const_element_iterator elem_it_end = mesh.subdomain_elements_end(r);
    for(; elem_it != elem_it_end; ++elem_it)
    {
#ifdef ENABLE_AMR
      // use parent_ID of -1 to indicate a level 0 element
      if ((*elem_it)->level() == 0)
      {
        elem_attribute.push_back(-1);
        elem_attribute.push_back(-1);
      }
      else
      {
        elem_attribute.push_back((*elem_it)->parent()->id());
        elem_attribute.push_back((*elem_it)->parent()->which_child_am_i((*elem_it)));
      }
#endif

#ifdef ENABLE_AMR
      elem_attribute.push_back (static_cast<int>((*elem_it)->level()));
      elem_attribute.push_back (static_cast<int>((*elem_it)->p_level()));
      elem_attribute.push_back (static_cast<int>((*elem_it)->refinement_flag()));
      elem_attribute.push_back (static_cast<int>((*elem_it)->p_refinement_flag()));
#endif

      elem_id.push_back( (*elem_it)->id() );

      switch( (*elem_it)->type() )
      {
      case TRI3        :
      case TRI3_FVM    :
        {
          elem.push_back( TRI_3 );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(0)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(1)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(2)] );
          break;
        }
      case QUAD4       :
      case QUAD4_FVM   :
        {
          elem.push_back( QUAD_4 );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(0)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(1)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(2)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(3)] );
          break;
        }
      case TET4        :
      case TET4_FVM    :
        {
          elem.push_back( TETRA_4 );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(0)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(1)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(2)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(3)] );
          break;
        }
      case PYRAMID5      :
      case PYRAMID5_FVM  :
        {
          elem.push_back( PYRA_5 );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(0)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(1)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(2)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(3)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(4)] );
          break;
        }
      case PRISM6      :
      case PRISM6_FVM  :
        {
          elem.push_back( PENTA_6 );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(0)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(1)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(2)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(3)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(4)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(5)] );
          break;
        }
      case HEX8        :
      case HEX8_FVM    :
        {
          elem.push_back( HEXA_8 );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(0)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(1)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(2)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(3)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(4)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(5)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(6)] );
          elem.push_back( node_pointer_to_id[(*elem_it)->get_node(7)] );
          break;
        }
      default:
        {
          MESSAGE<<"ERROR: Unsupported element type found during CGNS export."<<std::endl; RECORD();
          genius_error();
        }
      }

    }


    // map elem->id() to the cell's index for cgns writing in this region
    std::map<unsigned int, int> elem_id_to_region_cell_index ;
    for(unsigned int n=0; n<elem_id.size(); n++)
      elem_id_to_region_cell_index[elem_id[n]] = n+1;

    if( Genius::processor_id() == 0)
    {
      genius_assert(!cg_section_write(fn, B ,Z, "GridElements", MIXED, 1, elem_id.size(), 0, &elem[0], &S));
#ifdef ENABLE_AMR
      // write element amr attribute
      int  elem_id_size = static_cast<int>(elem_id.size());
      int  elem_attribute_size = static_cast<int>(elem_attribute.size());
      genius_assert(!cg_goto(fn, B, "Zone_t", Z , "end"));
      genius_assert(!cg_user_data_write ("Element_Attribute"));
      genius_assert(!cg_goto(fn, B, "Zone_t", Z, "UserDefinedData_t", 2, "end"));
      genius_assert(!cg_array_write("elem_id_array", Integer, 1, &elem_id_size, &elem_id[0]));
      genius_assert(!cg_array_write("elem_attribute_array", Integer, 1, &elem_attribute_size, &elem_attribute[0]));
#endif

    }


    // write boundary condition into cgns file

    // the element-face boundary information
    std::map<short int, std::pair<std::vector<int>, std::vector<int> > > bd_info;
    {
      std::vector<unsigned int>       el;
      std::vector<unsigned short int> sl;
      std::vector<short int>          il;
      system.mesh().boundary_info->build_side_list (el, sl, il);
      for(unsigned int n=0; n<sl.size(); n++)
      {
        const Elem * elem = mesh.elem(el[n]);
        if( elem->subdomain_id() != r) continue;

        genius_assert( elem_id_to_region_cell_index.find(el[n]) != elem_id_to_region_cell_index.end() );
        bd_info[il[n]].first.push_back( elem_id_to_region_cell_index[el[n]] );
        bd_info[il[n]].second.push_back( sl[n] );
      }
    }


    // the point type boundary
    const BoundaryConditionCollector * bcs = system.get_bcs();
    for(unsigned int n=0; n< bcs->n_bcs(); n++)
    {
      const BoundaryCondition * bc = bcs->get_bc(n);
      //skip inter connector here
      if(bc->is_inter_connect_hub()) continue;

      short int bd_id = bcs->get_bd_id_by_bc_index(n);

      // collect boundary nodes
      std::vector<int> bd_point;
      std::vector<const Node *> bd_nodes = bc->nodes();
      for(unsigned int i=0; i<bd_nodes.size(); i++)
        if( node_pointer_to_id.find(bd_nodes[i]) != node_pointer_to_id.end() )
          bd_point.push_back(node_pointer_to_id[bd_nodes[i]]);

      // write down boundary information
      if( Genius::processor_id() == 0  && bd_info[bd_id].second.size() && bd_point.size() )
      {
        // write a dummy bc label
        char bc_label[32];
        sprintf( bc_label, "boundary_%d", bd_id );

        // boundary as point list
        cg_boco_write(fn, B, Z, bc_label, BCTypeUserDefined, PointList, bd_point.size(), &bd_point[0], &BC);

        // actual bc label is here
        assert(!cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, "end"));
        assert(!cg_descriptor_write("bc_label", bc->label().c_str()));

        // write boundary condition
        assert(!cg_descriptor_write("bc_settings", bc->boundary_condition_in_string().c_str()));

        // extra information as element-side list
        std::vector<int> & elems =  bd_info[bd_id].first;
        std::vector<int> & sides =  bd_info[bd_id].second;

        int n_side =  sides.size();
        assert(!cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, "end"));
        assert(!cg_user_data_write ("element_side_information"));
        assert(!cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, "UserDefinedData_t", 1, "end"));
        assert(!cg_array_write("element_list", Integer, 1, &n_side, &elems[0]));
        assert(!cg_array_write("side_list", Integer, 1, &n_side, &sides[0]));

        // write electrode potential
        if( bc->is_electrode() )
        {
          double potential = bc->ext_circuit()->potential();
          int    DimensionVector = 1;

          assert(!cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, "end"));
          assert(!cg_user_data_write ("Extra_data_for_electrode"));
          assert(!cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, "UserDefinedData_t", 2, "end"));
          assert(!cg_array_write("electrode_potential", RealDouble, 1, &DimensionVector, &potential));
        }
      }
    }


    // write zone 1-to-1 connect information
    // does this really needed?


    // write solution data to cgns file
    switch ( region->type() )
    {

    case SemiconductorRegion :
      {
        std::vector<unsigned int> region_id;
        std::vector<double> Na, Nd, P, As, Sb, Bo;
        std::vector<double> psi, p, n, T, Tn, Tp, Eqc, Eqv;
        std::vector<double> mole_x, mole_y;
        std::map< std::string, std::vector<double> > custom_variable;
        std::map< std::string, int > custom_variable_dummy;

        bool sigle   = Material::IsSingleCompSemiconductor(region->material());
        bool complex = Material::IsComplexCompSemiconductor(region->material());

        SimulationRegion::const_node_iterator node_it = region->nodes_begin();
        local_id = 0;
        for(; node_it!=region->nodes_end(); ++node_it, ++local_id)
        {
          if( (*node_it).second->root_node()->processor_id() != Genius::processor_id()) continue;

          const FVM_NodeData * node_data = (*node_it).second->node_data();
          genius_assert(node_data);

          region_id.push_back(local_id);

          // doping information
          Na.push_back( node_data->Na()/ pow(cm, -3) );
          Nd.push_back( node_data->Nd()/ pow(cm, -3) );
          P.push_back ( node_data->P() / pow(cm, -3) );
          As.push_back( node_data->As()/ pow(cm, -3) );
          Sb.push_back( node_data->Sb()/ pow(cm, -3) );
          Bo.push_back( node_data->B() /  pow(cm, -3) );

          // field solution
          psi.push_back( node_data->psi()/ V );
          n.push_back  ( node_data->n()  / pow(cm, -3) );
          p.push_back  ( node_data->p()  / pow(cm, -3) );
          T.push_back  ( node_data->T()  / K );
          Tn.push_back ( node_data->Tn() / K );
          Tp.push_back ( node_data->Tp() / K );
          Eqc.push_back( node_data->Eqc()/ eV );
          Eqv.push_back( node_data->Eqv()/ eV );

          if(sigle)
            mole_x.push_back(node_data->mole_x());

          if(complex)
          {
            mole_x.push_back(node_data->mole_x());
            mole_y.push_back(node_data->mole_y());
          }

          // user defined values.
          // some of which are defined only on a subset of nodes, we collect all names here first
          std::vector<std::string> user_var_list;
          node_data->GetUserScalarList(user_var_list);
          for(std::vector<std::string>::iterator name_it=user_var_list.begin();
              name_it!=user_var_list.end(); name_it++)
          {
            if(custom_variable.find(*name_it) == custom_variable.end())
            {
              custom_variable_dummy.insert(std::pair<std::string, int >(*name_it,0));
            }
          }
        }

        // actually collecting user defined values, and use zero as default.
        {
          std::vector<unsigned int> n_var;
          n_var.push_back(custom_variable_dummy.size());
          Parallel::allgather(n_var);
          for(unsigned int i=0; i<Genius::n_processors(); i++)
          {
            if (Genius::processor_id() == i)
            {
              for(std::map< std::string, int >::iterator var_it=custom_variable_dummy.begin();
                  var_it!=custom_variable_dummy.end();var_it++)
              {
                std::string key;
                key = var_it->first;
                Parallel::broadcast(key,i);
                custom_variable.insert(std::pair<std::string,std::vector<double> >(key,std::vector<double>()));
              }
            }
            else
            {
              for(unsigned int j=0; j<n_var[i]; j++)
              {
                std::string key;
                Parallel::broadcast(key,i);
                custom_variable.insert(std::pair<std::string,std::vector<double> >(key,std::vector<double>()));
              }
            }
          }

          if (custom_variable.size()>0)
          {
            local_id=0;
            for(node_it=region->nodes_begin(); node_it!=region->nodes_end(); ++node_it, ++local_id)
            {
              if( (*node_it).second->root_node()->processor_id() != Genius::processor_id()) continue;

              const FVM_NodeData * node_data = (*node_it).second->node_data();
              genius_assert(node_data);

              for(std::map< std::string, std::vector<double> >::iterator var_it=custom_variable.begin();
                  var_it!=custom_variable.end();var_it++)
              {
                var_it->second.push_back(0.0);
              }

              std::vector<std::string> user_var_list;
              node_data->GetUserScalarList(user_var_list);
              for(std::vector<std::string>::const_iterator name_it=user_var_list.begin();
                  name_it!=user_var_list.end(); name_it++)
              {
                double v = node_data->UserScalarValue(*name_it);
                std::map< std::string, std::vector<double> >::iterator var_it = custom_variable.find(*name_it);
                var_it->second.back() = v;
              }
            }
          }
        }

        // synchronization data with other processor
        Parallel::gather(0, region_id);
        Parallel::gather(0, Na);
        Parallel::gather(0, Nd);
        Parallel::gather(0, P);
        Parallel::gather(0, As);
        Parallel::gather(0, Sb);
        Parallel::gather(0, Bo);
        Parallel::gather(0, psi);
        Parallel::gather(0, n);
        Parallel::gather(0, p);
        Parallel::gather(0, T);
        Parallel::gather(0, Tn);
        Parallel::gather(0, Tp);
        Parallel::gather(0, Eqc);
        Parallel::gather(0, Eqv);

        if(sigle)
          Parallel::gather(0, mole_x);

        if(complex)
        {
          Parallel::gather(0, mole_x);
          Parallel::gather(0, mole_y);
        }

        if(custom_variable.size()>0)
        {
          for(std::map< std::string, std::vector<double> >::iterator var_it=custom_variable.begin();
              var_it!=custom_variable.end();var_it++)
          {
            Parallel::gather(0,var_it->second);
          }
        }

        if( Genius::processor_id() == 0)
        {
          //write solution
          genius_assert(!cg_sol_write  (fn, B, Z, "Solution", Vertex, &SOL));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "elec_density",           &(sort_it(n,  region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "hole_density",           &(sort_it(p,  region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "potential",              &(sort_it(psi,region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "lattice_temperature",    &(sort_it(T,  region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "elec_temperature",       &(sort_it(Tn, region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "hole_temperature",       &(sort_it(Tp, region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "elec_quantum_potential", &(sort_it(Eqc,region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "hole_quantum_potential", &(sort_it(Eqv,region_id)[0]), &F));

          //write doping
          genius_assert(!cg_sol_write  (fn, B, Z, "Doping", Vertex, &SOL));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "Na", &(sort_it(Na, region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "Nd", &(sort_it(Nd, region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "P",  &(sort_it(P,  region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "As", &(sort_it(As, region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "Sb", &(sort_it(Sb, region_id)[0]), &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "B",  &(sort_it(Bo, region_id)[0]), &F));

          //write mole fraction
          if(sigle)
          {
            genius_assert(!cg_sol_write  (fn, B, Z, "Mole", Vertex, &SOL));
            genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "mole_x", &(sort_it(mole_x, region_id)[0]), &F));
          }

          if(complex)
          {
            genius_assert(!cg_sol_write  (fn, B, Z, "Mole", Vertex, &SOL));
            genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "mole_x", &(sort_it(mole_x, region_id)[0]), &F));
            genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "mole_y", &(sort_it(mole_y, region_id)[0]), &F));
          }

          //write user defined values
          if(custom_variable.size()>0)
          {
            genius_assert(!cg_sol_write  (fn, B, Z, "Custom", Vertex, &SOL));
            for(std::map< std::string, std::vector<double> >::iterator var_it=custom_variable.begin();
                var_it!=custom_variable.end();var_it++)
            {
              genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, var_it->first.c_str(), &(sort_it(var_it->second, region_id)[0]), &F));
            }
          }
        }

        break;
      }


    case InsulatorRegion     :
      {
        std::vector<unsigned int> region_id;
        std::vector<double> psi, T;
        std::map< std::string, std::vector<double> > custom_variable;
        std::map< std::string, int > custom_variable_dummy;

        SimulationRegion::const_node_iterator node_it = region->nodes_begin();
        local_id = 0;
        for(; node_it!=region->nodes_end(); ++node_it, ++local_id)
        {

          if( (*node_it).second->root_node()->processor_id() != Genius::processor_id()) continue;

          const FVM_NodeData * node_data = (*node_it).second->node_data();
          genius_assert(node_data);

          region_id.push_back(local_id);

          // field solution
          psi.push_back( node_data->psi()/ V );
          T.push_back  ( node_data->T()  / K );

          // user defined values.
          // some of which are defined only on a subset of nodes, we collect all names here first
          std::vector<std::string> user_var_list;
          node_data->GetUserScalarList(user_var_list);
          for(std::vector<std::string>::iterator name_it=user_var_list.begin();
              name_it!=user_var_list.end(); name_it++)
          {
            if(custom_variable.find(*name_it) == custom_variable.end())
            {
              custom_variable_dummy.insert(std::pair<std::string, int >(*name_it,0));
            }
          }

        }

        // actually collecting user defined values, and use zero as default.
        {
          std::vector<unsigned int> n_var;
          n_var.push_back(custom_variable_dummy.size());
          Parallel::allgather(n_var);
          for(unsigned int i=0; i<Genius::n_processors(); i++)
          {
            if (Genius::processor_id() == i)
            {
              for(std::map< std::string, int >::iterator var_it=custom_variable_dummy.begin();
                  var_it!=custom_variable_dummy.end();var_it++)
              {
                std::string key;
                key = var_it->first;
                Parallel::broadcast(key,i);
                custom_variable.insert(std::pair<std::string,std::vector<double> >(key,std::vector<double>()));
              }
            }
            else
            {
              for(unsigned int j=0; j<n_var[i]; j++)
              {
                std::string key;
                Parallel::broadcast(key,i);
                custom_variable.insert(std::pair<std::string,std::vector<double> >(key,std::vector<double>()));
              }
            }
          }

          if (custom_variable.size()>0)
          {
            local_id=0;
            for(node_it=region->nodes_begin(); node_it!=region->nodes_end(); ++node_it, ++local_id)
            {
              if( (*node_it).second->root_node()->processor_id() != Genius::processor_id()) continue;

              const FVM_NodeData * node_data = (*node_it).second->node_data();
              genius_assert(node_data);

              for(std::map< std::string, std::vector<double> >::iterator var_it=custom_variable.begin();
                  var_it!=custom_variable.end();var_it++)
              {
                var_it->second.push_back(0.0);
              }

              std::vector<std::string> user_var_list;
              node_data->GetUserScalarList(user_var_list);
              for(std::vector<std::string>::const_iterator name_it=user_var_list.begin();
                  name_it!=user_var_list.end(); name_it++)
              {
                double v = node_data->UserScalarValue(*name_it);
                std::map< std::string, std::vector<double> >::iterator var_it = custom_variable.find(*name_it);
                var_it->second.back() = v;
              }
            }
          }
        }


        // synchronization data with other processor
        Parallel::gather(0, region_id);
        Parallel::gather(0, psi);
        Parallel::gather(0, T);

        if(custom_variable.size()>0)
        {
          for(std::map< std::string, std::vector<double> >::iterator var_it=custom_variable.begin();
              var_it!=custom_variable.end();var_it++)
          {
            Parallel::gather(0,var_it->second);
          }
        }


        if( Genius::processor_id() == 0)
        {
          genius_assert(!cg_sol_write  (fn, B, Z, "Solution", Vertex, &SOL));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "potential",              &(sort_it(psi,region_id)[0]),   &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "lattice_temperature",    &(sort_it(T,  region_id)[0]),   &F));

          //write user defined values
          if(custom_variable.size()>0)
          {
            genius_assert(!cg_sol_write  (fn, B, Z, "Custom", Vertex, &SOL));
            for(std::map< std::string, std::vector<double> >::iterator var_it=custom_variable.begin();
                var_it!=custom_variable.end();var_it++)
            {
              genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, var_it->first.c_str(), &(sort_it(var_it->second, region_id)[0]), &F));
            }
          }
        }
        break;
      }


    case ConductorRegion     :
      {
        std::vector<unsigned int> region_id;
        std::vector<double> psi, T;
        std::map< std::string, std::vector<double> > custom_variable;
        std::map< std::string, int > custom_variable_dummy;

        SimulationRegion::const_node_iterator node_it = region->nodes_begin();
        local_id = 0;
        for(; node_it!=region->nodes_end(); ++node_it, ++local_id)
        {

          if( (*node_it).second->root_node()->processor_id() != Genius::processor_id()) continue;

          const FVM_NodeData * node_data = (*node_it).second->node_data();
          genius_assert(node_data);

          region_id.push_back(local_id);

          // field solution
          psi.push_back( node_data->psi()/ V );
          T.push_back  ( node_data->T()  / K );

          // user defined values.
          // some of which are defined only on a subset of nodes, we collect all names here first
          std::vector<std::string> user_var_list;
          node_data->GetUserScalarList(user_var_list);
          for(std::vector<std::string>::iterator name_it=user_var_list.begin();
              name_it!=user_var_list.end(); name_it++)
          {
            if(custom_variable.find(*name_it) == custom_variable.end())
            {
              custom_variable_dummy.insert(std::pair<std::string, int >(*name_it,0));
            }
          }

        }

        // actually collecting user defined values, and use zero as default.
        {
          std::vector<unsigned int> n_var;
          n_var.push_back(custom_variable_dummy.size());
          Parallel::allgather(n_var);
          for(unsigned int i=0; i<Genius::n_processors(); i++)
          {
            if (Genius::processor_id() == i)
            {
              for(std::map< std::string, int >::iterator var_it=custom_variable_dummy.begin();
                  var_it!=custom_variable_dummy.end();var_it++)
              {
                std::string key;
                key = var_it->first;
                Parallel::broadcast(key,i);
                custom_variable.insert(std::pair<std::string,std::vector<double> >(key,std::vector<double>()));
              }
            }
            else
            {
              for(unsigned int j=0; j<n_var[i]; j++)
              {
                std::string key;
                Parallel::broadcast(key,i);
                custom_variable.insert(std::pair<std::string,std::vector<double> >(key,std::vector<double>()));
              }
            }
          }

          if (custom_variable.size()>0)
          {
            local_id=0;
            for(node_it=region->nodes_begin(); node_it!=region->nodes_end(); ++node_it, ++local_id)
            {
              if( (*node_it).second->root_node()->processor_id() != Genius::processor_id()) continue;

              const FVM_NodeData * node_data = (*node_it).second->node_data();
              genius_assert(node_data);

              for(std::map< std::string, std::vector<double> >::iterator var_it=custom_variable.begin();
                  var_it!=custom_variable.end();var_it++)
              {
                var_it->second.push_back(0.0);
              }

              std::vector<std::string> user_var_list;
              node_data->GetUserScalarList(user_var_list);
              for(std::vector<std::string>::const_iterator name_it=user_var_list.begin();
                  name_it!=user_var_list.end(); name_it++)
              {
                double v = node_data->UserScalarValue(*name_it);
                std::map< std::string, std::vector<double> >::iterator var_it = custom_variable.find(*name_it);
                var_it->second.back() = v;
              }
            }
          }
        }


        // synchronization data with other processor
        Parallel::gather(0, region_id);
        Parallel::gather(0, psi);
        Parallel::gather(0, T);

        if(custom_variable.size()>0)
        {
          for(std::map< std::string, std::vector<double> >::iterator var_it=custom_variable.begin();
              var_it!=custom_variable.end();var_it++)
          {
            Parallel::gather(0,var_it->second);
          }
        }


        if( Genius::processor_id() == 0)
        {
          genius_assert(!cg_sol_write  (fn, B, Z, "Solution", Vertex, &SOL));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "potential",              &(sort_it(psi,region_id)[0]),   &F));
          genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, "lattice_temperature",    &(sort_it(T,  region_id)[0]),   &F));

          //write user defined values
          if(custom_variable.size()>0)
          {
            genius_assert(!cg_sol_write  (fn, B, Z, "Custom", Vertex, &SOL));
            for(std::map< std::string, std::vector<double> >::iterator var_it=custom_variable.begin();
                var_it!=custom_variable.end();var_it++)
            {
              genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, var_it->first.c_str(), &(sort_it(var_it->second, region_id)[0]), &F));
            }
          }
        }
        break;
      }

    // no solution data in vacuum region?
    case VacuumRegion     :
      break;

    case PMLRegion     :
      break;

    default :
      {
        MESSAGE<<"ERROR: Unsupported region type found during CGNS export."<<std::endl; RECORD();
        genius_error();
      }

    }


  }


  // close CGNS file
  if( Genius::processor_id() == 0)
    cg_close(fn);

}


