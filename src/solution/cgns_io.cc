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
#include "boundary_info.h"
#include "material_define.h"
#include "simulation_system.h"
#include "boundary_condition_collector.h"

#include "parallel.h"
#include "mesh_communication.h"

using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::K;
using PhysicalUnit::eV;



// ------------------------------------------------------------
// CGNSIO class members
// ------------------------------------------------------------


void CGNSIO::read (const std::string& filename)
{
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

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
    char base_name[32];
    //read base information
    genius_assert(!cg_base_read(fn, B, base_name, &cell_dim, &physical_dim));
    genius_assert(physical_dim==3);

    // get mesh magic number from base name
    {
      std::string base_name_string = base_name;
      //since the base name has the format of "GENIUS_Mesh_magicnum" we first get the substr of magicnum
      std::string base_name_prefix = "GENIUS_Mesh_";
      std::string magic_num_string = base_name_string.substr(base_name_prefix.length());

      // translate string to int
      std::stringstream   ss;
      ss << magic_num_string;
      ss >> mesh.magic_num();
    }

    {
      genius_assert(!cg_goto(fn, B, "end"));
      int nd;
      genius_assert(!cg_ndescriptors(&nd));
      for(int d=1; d<=nd; ++d)
      {
        char descriptor_name[32];
        char *_description;
        genius_assert(!cg_descriptor_read(d, descriptor_name, &_description));

        if(std::string(descriptor_name) == "ResistiveMetal")
          system.set_resistive_metal_mode(true);
        cg_free(_description);
      }
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
      char           zone_name[32];
      {
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
          cg_free(material);
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
        std::string path = std::string("/") + base_name + "/" + zone_name + "/" + "Global_Node_Index";
        genius_assert(!cg_gopath(fn, path.c_str()));
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
          std::string path = std::string("/") + base_name + "/" + zone_name + "/" + "Element_Attribute";
          genius_assert(!cg_gopath(fn, path.c_str()));
          int A;
          genius_assert(!cg_narrays( &A ));
          genius_assert(A==2);

          for(int a=1; a<=A; a++)
          {
            char       ArrayName[32];
            DataType_t DataType;
            int        DataDimension;
            int        DimensionVector;

            genius_assert(!cg_array_info(a, ArrayName , &DataType , &DataDimension , &DimensionVector ));

            if( std::string(ArrayName) == "elem_id_array" )
            {
              genius_assert(DataType == Integer);
              elem_id.resize(DimensionVector);
              genius_assert(!cg_array_read_as(a, DataType, &elem_id[0] ));
            }

            if( std::string(ArrayName) == "elem_attribute_array" )
            {
              genius_assert(DataType == Integer);
              elem_attribute.resize(DimensionVector);
              genius_assert(!cg_array_read_as(a, DataType, &elem_attribute[0] ));
            }
          }
        }
#endif
        std::vector<int>::const_iterator it = elem.begin();
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
                elem->prepare_for_fvm();
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
                elem->prepare_for_fvm();
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
                elem->prepare_for_fvm();
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
                elem->prepare_for_fvm();
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
                elem->prepare_for_fvm();
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
                elem->prepare_for_fvm();
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
        char           boco_name[32];
        BCType_t       boco_type ;
        PointSetType_t ptset_type ;
        int            npnts ;
        int            normalindex ;
        int            normallistflag;
        DataType_t     normaldatatype;
        int            ndataset;
        int            nuserdata;

        genius_assert(!cg_boco_info(fn, B, z_id, bc_id, boco_name, &boco_type, &ptset_type, &npnts, &normalindex, &normallistflag, &normaldatatype, &ndataset) );
        genius_assert(ptset_type==PointList);
        // for old cgns file, the bc_label is boco_name
        // however, new cgns file will save bc label in descriptor, thus support long (>32 chars) bc label
        std::string bc_label(boco_name);
        // new cgns file will contain string for boundary condition as well
        std::string bc_settings;
        {
          genius_assert(!cg_goto(fn, B, "Zone_t", z_id, "ZoneBC_t", 1, "BC_t", bc_id, "end"));
          int ndescriptor;
          genius_assert(!cg_ndescriptors(&ndescriptor));

          char descriptorname[32];
          for(int nd=1; nd<=ndescriptor; ++nd)
          {
            char *_text;
            assert(!cg_descriptor_read(nd, descriptorname, &_text));
            if(std::string(descriptorname) == "bc_label")
            {
              bc_label = _text;
            }
            if(std::string(descriptorname) == "bc_settings")
            {
              bc_settings = _text;
            }
            cg_free(_text);
          }

        }

        // get the number of user defined data structure
        genius_assert(!cg_goto(fn, B, "Zone_t", z_id, "ZoneBC_t", 1, "BC_t", bc_id, "end"));
        genius_assert(!cg_nuser_data(&nuserdata));



        std::vector<int> bd_elem;
        std::vector<int> bd_side;
        // read the element-side array
        {
          genius_assert(!cg_gopath(fn, "element_side_information"));

          int A;
          genius_assert(!cg_narrays( &A ));
          genius_assert(A==2);

          for(int a=1; a<=A; a++)
          {
            char       ArrayName[32];
            DataType_t DataType;
            int        DataDimension;
            int        DimensionVector;

            genius_assert(!cg_array_info(a, ArrayName , &DataType , &DataDimension , &DimensionVector ));

            if( std::string(ArrayName) == "element_list" )
            {
              genius_assert(DataType == Integer);
              bd_elem.resize(DimensionVector);
              genius_assert(!cg_array_read_as(a, Integer, &bd_elem[0] ));
            }

            if( std::string(ArrayName) == "side_list" )
            {
              genius_assert(DataType == Integer);
              bd_side.resize(DimensionVector);
              genius_assert(!cg_array_read_as(a, Integer, &bd_side[0] ));
            }
          }
          genius_assert(!cg_gopath(fn, ".."));
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
          if(!cg_gopath(fn, "extra_data_for_electrode"))
          {
            int narray;
            genius_assert(!cg_narrays(&narray));
            for(int a=1; a<=narray; a++)
            {
              char array[32];
              DataType_t type;
              int dimension;
              int vector;
              double data;
              genius_assert(!cg_array_info(a, array, &type, &dimension, &vector));
              genius_assert(!cg_array_read_as(a, type, &data));
              if( std::string(array) == "electrode_potential" )
                electrode_potential[bd_id] = data;
              if( std::string(array) == "electrode_potential_old" )
                electrode_potential_old[bd_id] = data;
              if( std::string(array) == "electrode_vapp" )
                electrode_vapp[bd_id] = data;
              if( std::string(array) == "electrode_iapp" )
                electrode_iapp[bd_id] = data;
            }
            genius_assert(!cg_gopath(fn, ".."));
          }
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

          std::pair<std::string, std::string> key(solutionname,fieldname);
          sol_map[key].resize(isize[0]);
          genius_assert(!cg_field_read(fn, B, z_id, sol_id, fieldname, datatype, &imin, &isize[0], &((sol_map[key])[0])));
        }

        std::string path = std::string("/") + base_name + "/" + zone_name + "/" + solutionname;
        genius_assert(!cg_gopath(fn, path.c_str()));

        int ndescriptors;
        genius_assert(!cg_ndescriptors(&ndescriptors));
        for(int d_it=1; d_it<=ndescriptors; d_it++)
        {
          char descriptorname[32];
          char *_description;
          genius_assert(!cg_descriptor_read(d_it, descriptorname, &_description));
          region_solution_units.insert( std::make_pair( std::string(descriptorname),  std::string(_description)) );
          cg_free(_description);
        }
      }

      // save solution name to solution value map for later usage
      region_solutions.push_back(sol_map);

    }// search for each zone

    // read extra boundary info as interconnect or change boundary
    {
      genius_assert(!cg_goto(fn, B, "end"));
      int nuserdata; genius_assert(!cg_nuser_data(&nuserdata));
      for(int u=1; u<=nuserdata; ++u)
      {
        char userdataname[32];
        genius_assert(!cg_user_data_read(u, userdataname));
        if( std::string(userdataname) != "ExtraBoundaryInfo" ) continue;

        genius_assert(!cg_goto(fn, B, "UserDefinedData_t", u, "end"));
        int ndescriptor;
        genius_assert(!cg_ndescriptors(&ndescriptor));
        for(int i=1; i<=ndescriptor; ++i)
        {
          // read first descriptor -- bc label
          char descriptorname[32];
          char *_description;
          assert(!cg_descriptor_read(i, descriptorname, &_description));
          mesh.boundary_info->add_extra_description(std::string(_description));
          cg_free(_description);
        }
        genius_assert(!cg_goto(fn, B, "end"));
      }
    }

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

  Parallel::broadcast(region_solution_units   , 0);


  // now all the processor have enough information for
  // updating region solution data
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);

    const std::vector<int> & global_id = region_global_id[r];
    const std::map<int, unsigned int > & global_id_to_node_id = region_global_id_to_node_id[r];
    const std::map<std::pair<std::string,std::string>, std::vector<double> > & solution = region_solutions[r];

    switch ( region->type() )
    {
        case SemiconductorRegion :
        {

          for(std::map< std::pair<std::string, std::string>, std::vector<double> >::const_iterator var_it = solution.begin();
              var_it!=solution.end(); var_it++)
          {
            if(var_it->first.first == "Custom")
            {
              const std::string & variable = var_it->first.second;
              const std::string & variable_unit = region_solution_units[variable+".unit"];
              region->add_variable( SimulationVariable(variable, SCALAR, POINT_CENTER, variable_unit, invalid_uint, true, true)  );
            }
          }

          for(unsigned int n=0; n<global_id.size(); n++)
          {
            unsigned int node_id = global_id_to_node_id.find(global_id[n])->second;
            FVM_Node * fvm_node = region->region_fvm_node(node_id);
            if( !fvm_node || !fvm_node->root_node()->on_local() ) continue;

            FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

            for(std::map< std::pair<std::string,std::string>, std::vector<double> >::const_iterator var_it = solution.begin();
                var_it!=solution.end(); var_it++)
            {
              if(var_it->first.first == "Solution")
              {
                // field solution
                if(var_it->first.second == "electron" || var_it->first.second == "elec_density")
                {
                  node_data->n() = var_it->second[n] * pow(cm,-3);
                  continue;
                }
                if(var_it->first.second == "hole" || var_it->first.second == "hole_density")
                {
                  node_data->p() = var_it->second[n] * pow(cm,-3);
                  continue;
                }
                if(var_it->first.second == "potential")
                {
                  node_data->psi() = var_it->second[n] * V;
                  continue;
                }
                if(var_it->first.second == "temperature" || var_it->first.second == "lattice_temperature" )
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
              }
              if(var_it->first.first == "Doping")
              {
                // doping
                if(var_it->first.second == "Na" || var_it->first.second == "na")
                {
                  node_data->Na()  = var_it->second[n] * pow(cm, -3);
                  continue;
                }
                if(var_it->first.second == "Nd" || var_it->first.second == "nd")
                {
                  node_data->Nd()  = var_it->second[n] * pow(cm, -3);
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
                const SimulationVariable & variable = region->get_variable(var_it->first.second, POINT_CENTER);
                node_data->data<Real>(variable.variable_index) = var_it->second[n]*variable.variable_unit;
              }

            }
          }

          // after import previous solutions, we re-init region here
          region->reinit_after_import();

          break;
        }
        case InsulatorRegion     :
        {
          for(std::map< std::pair<std::string, std::string>, std::vector<double> >::const_iterator var_it = solution.begin();
              var_it!=solution.end(); var_it++)
          {
            if(var_it->first.first == "Custom")
            {
              const std::string & variable = var_it->first.second;
              const std::string & variable_unit = region_solution_units[variable+".unit"];
              region->add_variable( SimulationVariable(variable, SCALAR, POINT_CENTER, variable_unit, invalid_uint, true, true)  );
            }
          }

          for(unsigned int n=0; n<global_id.size(); n++)
          {
            unsigned int node_id = global_id_to_node_id.find(global_id[n])->second;
            FVM_Node * fvm_node = region->region_fvm_node(node_id);
            if( !fvm_node || !fvm_node->root_node()->on_local() ) continue;

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
                // field solution
                if(var_it->first.second == "electron" || var_it->first.second == "elec_density")
                {
                  node_data->n() = var_it->second[n] * pow(cm,-3);
                  continue;
                }
                if(var_it->first.second == "hole" || var_it->first.second == "hole_density")
                {
                  node_data->p() = var_it->second[n] * pow(cm,-3);
                  continue;
                }
                if(var_it->first.second == "temperature" || var_it->first.second == "lattice_temperature" )
                {
                  node_data->T() = var_it->second[n] * K;
                  continue;
                }
              }

              if(var_it->first.first == "Custom")
              {
                const SimulationVariable & variable = region->get_variable(var_it->first.second, POINT_CENTER);
                node_data->data<Real>(variable.variable_index) = var_it->second[n]*variable.variable_unit;
              }
            }
          }

          // after import previous solutions, we re-init region here
          region->reinit_after_import();

          break;
        }

        case ElectrodeRegion     :
        case MetalRegion         :
        {
          for(std::map< std::pair<std::string, std::string>, std::vector<double> >::const_iterator var_it = solution.begin();
              var_it!=solution.end(); var_it++)
          {
            if(var_it->first.first == "Custom")
            {
              const std::string & variable = var_it->first.second;
              const std::string & variable_unit = region_solution_units[variable+".unit"];
              region->add_variable( SimulationVariable(variable, SCALAR, POINT_CENTER, variable_unit, invalid_uint, true, true)  );
            }
          }

          for(unsigned int n=0; n<global_id.size(); n++)
          {
            unsigned int node_id = global_id_to_node_id.find(global_id[n])->second;
            FVM_Node * fvm_node = region->region_fvm_node(node_id);
            if( !fvm_node || !fvm_node->root_node()->on_local() ) continue;

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
                if(var_it->first.second == "electron" || var_it->first.second == "elec_density")
                {
                  node_data->n() = var_it->second[n] * pow(cm,-3);
                  continue;
                }
                if(var_it->first.second == "temperature" || var_it->first.second == "lattice_temperature" )
                {
                  node_data->T() = var_it->second[n] * K;
                  continue;
                }
              }

              if(var_it->first.first == "Custom")
              {
                const SimulationVariable & variable = region->get_variable(var_it->first.second, POINT_CENTER);
                node_data->data<Real>(variable.variable_index) = var_it->second[n]*variable.variable_unit;
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

  // set potential/Vapp/Iapp of electrode
  std::map< short int, double >::iterator el_it;

  Parallel::broadcast(electrode_potential, 0);
  for(el_it = electrode_potential.begin(); el_it!=electrode_potential.end(); ++el_it )
  {
    std::string bc_label= mesh.boundary_info->get_label_by_id(el_it->first);
    BoundaryCondition * bc = system.get_bcs()->get_bc(bc_label);
    if(bc && bc->is_electrode())
      bc->ext_circuit()->potential() = el_it->second * V;
  }

  Parallel::broadcast(electrode_potential_old, 0);
  for(el_it = electrode_potential_old.begin(); el_it!=electrode_potential_old.end(); ++el_it )
  {
    std::string bc_label= mesh.boundary_info->get_label_by_id(el_it->first);
    BoundaryCondition * bc = system.get_bcs()->get_bc(bc_label);
    if(bc && bc->is_electrode())
      bc->ext_circuit()->potential_old() = el_it->second * V;
  }

  Parallel::broadcast(electrode_vapp, 0);
  for(el_it = electrode_vapp.begin(); el_it!=electrode_vapp.end(); ++el_it )
  {
    std::string bc_label= mesh.boundary_info->get_label_by_id(el_it->first);
    BoundaryCondition * bc = system.get_bcs()->get_bc(bc_label);
    if(bc && bc->is_electrode())
      bc->ext_circuit()->Vapp() = el_it->second * V;
  }

  Parallel::broadcast(electrode_iapp, 0);
  for(el_it = electrode_iapp.begin(); el_it!=electrode_iapp.end(); ++el_it )
  {
    std::string bc_label= mesh.boundary_info->get_label_by_id(el_it->first);
    BoundaryCondition * bc = system.get_bcs()->get_bc(bc_label);
    if(bc && bc->is_electrode())
      bc->ext_circuit()->Iapp() = el_it->second * A;
  }

  system.init_region_post_process();
}




void CGNSIO::write (const std::string& filename)
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  const BoundaryConditionCollector * bcs = system.get_bcs();
  const MeshBase & mesh = system.mesh();

  if( Genius::processor_id() == 0)
  {
    // remove old file if exist
    remove(filename.c_str());
  }

  // classify node to region
  std::vector< std::vector<const Node *> > region_node_array(system.n_regions());
  if( Genius::processor_id() == 0)
  {
    std::vector< std::set<unsigned int> > region_node_id_array(system.n_regions());
    MeshBase::const_element_iterator elem_it = mesh.elements_begin();
    MeshBase::const_element_iterator elem_it_end = mesh.elements_end();
    for(; elem_it != elem_it_end; ++elem_it)
    {
      unsigned int region = (*elem_it)->subdomain_id();
      for(unsigned int n=0; n<(*elem_it)->n_nodes(); ++n)
      {
        region_node_id_array[region].insert((*elem_it)->get_node(n)->id());
      }
    }
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      std::set<unsigned int>::const_iterator it = region_node_id_array[r].begin();
      std::set<unsigned int>::const_iterator it_end = region_node_id_array[r].end();
      for(; it!=it_end; ++it)
        region_node_array[r].push_back( mesh.node_ptr(*it) );
    }
  }


  // classify elem to region
  std::vector< std::vector<const Elem *> > region_elem_array(system.n_regions());
  if( Genius::processor_id() == 0)
  {
    MeshBase::const_element_iterator elem_it = mesh.elements_begin();
    MeshBase::const_element_iterator elem_it_end = mesh.elements_end();
    for(elem_it = mesh.elements_begin(); elem_it != elem_it_end; ++elem_it)
    {
      region_elem_array[(*elem_it)->subdomain_id()].push_back(*elem_it);
    }
  }

  // get all the boundary elem-side-id trip
  std::vector<unsigned int>       boundary_el;
  std::vector<unsigned short int> boundary_sl;
  std::vector<short int>          boundary_il;
  // get nodes on each boundary
  std::map<short int, std::set<const Node *> > boundary_side_nodes_id_map;
  // classify boundaries to each region
  std::vector< std::vector<unsigned int> >region_boundary_array(system.n_regions());
  if( Genius::processor_id() == 0)
  {
    mesh.boundary_info->build_side_list (boundary_el, boundary_sl, boundary_il);
    mesh.boundary_info->boundary_side_nodes_with_id(boundary_side_nodes_id_map );
    for(unsigned int n=0; n<boundary_el.size(); n++)
    {
      const Elem * elem = mesh.elem(boundary_el[n]);
      region_boundary_array[elem->subdomain_id()].push_back(n);
    }
  }

  std::string base_name;

  // ok, create cgns here
  if( Genius::processor_id() == 0)
  {
    // open CGNS file for write
    genius_assert(!cg_open(filename.c_str(), MODE_WRITE, &fn));

    // create base of three dimensional mesh, here mesh takes its magic number as postfix
    std::stringstream   ss;
    ss << "GENIUS_Mesh_" << system.mesh().magic_num();
    ss >> base_name;

    genius_assert(!cg_base_write(fn, base_name.c_str(), 3, 3, &B));

    genius_assert(!cg_goto(fn,B,"end"));

    if(system.resistive_metal_mode())
      genius_assert(!cg_descriptor_write("ResistiveMetal", "true"));
  }

  // create cgns zone
  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    int size[3];
    // region node number
    size[0] = region_node_array[r].size();

    //region cell number
    size[1] = region_elem_array[r].size();

    //boundary cell number
    size[2] = 0;

    std::string zone_name = region->name();

    if( Genius::processor_id() == 0)
    {
      // write region name
      genius_assert(!cg_zone_write(fn, B, zone_name.c_str(), size, Unstructured, &Z));

      // goto the current region
      genius_assert(!cg_goto(fn,B,"Zone_t",Z,"end"));

      // write down material
      genius_assert(!cg_descriptor_write("Material", region->material().c_str() ));
    }


    // map node id to local node index in this region, alloc max_node_id with invalid_uint
    // it is a bit overkill in memory. however, i think vector is faster than map<id, local_id>
    std::vector<unsigned int> node_id_to_region_node_id(mesh.max_node_id (), invalid_uint);

    if( Genius::processor_id() == 0)
    {
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      x.reserve(size[0]);
      y.reserve(size[0]);
      z.reserve(size[0]);

      // addtional information about global id of region node
      std::vector<int> global_node_id;
      global_node_id.reserve(size[0]);

      // set coordinates of each node in this zone
      unsigned int local_id = 1;
      for(unsigned int n=0; n<region_node_array[r].size(); ++n, ++local_id)
      {
        //convert the unit to cm
        const Node * node = region_node_array[r][n];
        x.push_back( (*node)(0)/cm );
        y.push_back( (*node)(1)/cm );
        z.push_back( (*node)(2)/cm );
        // save global index of the region node
        global_node_id.push_back(node->id());
        // save local node index, the we can find local node index by node id
        node_id_to_region_node_id[node->id()] = local_id;
      }

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

    // map elem->id() to the cell's index for cgns writing in this region
    std::map<unsigned int, int> elem_id_to_region_cell_index ;

    if( Genius::processor_id() == 0)
    {
      const std::vector<const Elem *> & region_elem = region_elem_array[r];

      std::vector<int> elem_package;           // for element connectivity
      std::vector<int> elem_attribute;         // for element attribute
      std::vector<unsigned int> elem_id;
      elem_package.reserve(10*region_elem.size());// a bit overkill
      elem_attribute.reserve(6*region_elem.size());
      elem_id.reserve(region_elem.size());

      for(unsigned int n=0; n < region_elem.size(); ++n)
      {
        const Elem * elem = region_elem[n];
#ifdef ENABLE_AMR
        // use parent_ID of -1 to indicate a level 0 element
        if (elem->level() == 0)
        {
          elem_attribute.push_back(-1);
          elem_attribute.push_back(-1);
        }
        else
        {
          elem_attribute.push_back(elem->parent()->id());
          elem_attribute.push_back(elem->parent()->which_child_am_i(elem));
        }
#endif

#ifdef ENABLE_AMR
        elem_attribute.push_back (static_cast<int>(elem->level()));
        elem_attribute.push_back (static_cast<int>(elem->p_level()));
        elem_attribute.push_back (static_cast<int>(elem->refinement_flag()));
        elem_attribute.push_back (static_cast<int>(elem->p_refinement_flag()));
#endif

        elem_id.push_back( elem->id() );

        switch( elem->type() )
        {
            case TRI3        :
            case TRI3_FVM    :
            case TRI3_CY_FVM :
            {
              elem_package.push_back( TRI_3 );
              break;
            }
            case QUAD4       :
            case QUAD4_FVM   :
            case QUAD4_CY_FVM:
            {
              elem_package.push_back( QUAD_4 );
              break;
            }
            case TET4        :
            case TET4_FVM    :
            {
              elem_package.push_back( TETRA_4 );
              break;
            }
            case PYRAMID5      :
            case PYRAMID5_FVM  :
            {
              elem_package.push_back( PYRA_5 );
              break;
            }
            case PRISM6      :
            case PRISM6_FVM  :
            {
              elem_package.push_back( PENTA_6 );
              break;
            }
            case HEX8        :
            case HEX8_FVM    :
            {
              elem_package.push_back( HEXA_8 );
              break;
            }
            default:
            {
              MESSAGE<<"ERROR: Unsupported element type found during CGNS export."<<std::endl; RECORD();
              genius_error();
            }
        }

        for(unsigned int n=0; n<elem->n_nodes(); ++n)
          elem_package.push_back( node_id_to_region_node_id[elem->get_node(n)->id()] );

      }


      for(unsigned int n=0; n<elem_id.size(); n++)
        elem_id_to_region_cell_index[elem_id[n]] = n+1;

      genius_assert(!cg_section_write(fn, B ,Z, "GridElements", MIXED, 1, elem_id.size(), 0, &elem_package[0], &S));
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
    if( Genius::processor_id() == 0 )
    {
      const std::vector<unsigned int> & region_boundary = region_boundary_array[r];

      // the element-face boundary information
      std::map<short int, std::pair<std::vector<int>, std::vector<int> > > bd_info;
      //std::map<short int, std::set< const Node * > > bd_node_info;
      {
        for(unsigned int n=0; n<region_boundary.size(); n++)
        {
          unsigned int index = region_boundary[n];
          const Elem * elem = mesh.elem(boundary_el[index]);
          genius_assert( elem->subdomain_id() == r);

          bd_info[boundary_il[index]].first.push_back( elem_id_to_region_cell_index[boundary_el[index]] );
          bd_info[boundary_il[index]].second.push_back( boundary_sl[index] );
        }
      }

      // the point type boundary
      std::map<short int, std::pair<std::vector<int>, std::vector<int> > > ::const_iterator bd_info_it = bd_info.begin();
      for( ; bd_info_it != bd_info.end(); ++bd_info_it)
      {
        short int bd_id = bd_info_it->first;
        const std::vector<int> & elems =  bd_info_it->second.first;
        const std::vector<int> & sides =  bd_info_it->second.second;

        if( elems.empty() ) continue;

        const BoundaryCondition * bc = bcs->get_bc_by_bd_id(bd_id);

        // collect boundary nodes
        std::vector<int> bd_point;
        {
          const std::set<const Node *> & bd_nodes = boundary_side_nodes_id_map.find(bd_id)->second;
          for(std::set<const Node *>::const_iterator it=bd_nodes.begin(); it!=bd_nodes.end(); ++it)
              if( node_id_to_region_node_id[(*it)->id()] != invalid_uint )
                bd_point.push_back(node_id_to_region_node_id[(*it)->id()]);
        }

        // now write down boundary information

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
        int n_side =  sides.size();
        assert(!cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, "end"));
        assert(!cg_user_data_write ("element_side_information"));
        assert(!cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, "UserDefinedData_t", 1, "end"));
        assert(!cg_array_write("element_list", Integer, 1, &n_side, &elems[0]));
        assert(!cg_array_write("side_list", Integer, 1, &n_side, &sides[0]));

        // write electrode potential
        if( bc->is_electrode() )
        {
          double potential = bc->ext_circuit()->potential()/V;
          double potential_old = bc->ext_circuit()->potential_old()/V;
          double vapp = bc->ext_circuit()->Vapp()/V;
          double iapp = bc->ext_circuit()->Iapp()/A;
          int    DimensionVector = 1;
          assert(!cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, "end"));
          assert(!cg_user_data_write ("extra_data_for_electrode"));
          assert(!cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, "UserDefinedData_t", 2, "end"));
          assert(!cg_array_write("electrode_potential", RealDouble, 1, &DimensionVector, &potential));
          assert(!cg_array_write("electrode_potential_old", RealDouble, 1, &DimensionVector, &potential_old));
          assert(!cg_array_write("electrode_vapp", RealDouble, 1, &DimensionVector, &vapp));
          assert(!cg_array_write("electrode_iapp", RealDouble, 1, &DimensionVector, &iapp));
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
          bool sigle   = Material::IsSingleCompSemiconductor(region->material());
          bool complex = Material::IsComplexCompSemiconductor(region->material());

          std::vector<unsigned int> region_node_id;

          std::multimap< std::string, std::pair<SimulationVariable, std::vector<double> > >  region_data;
          typedef std::multimap< std::string, std::pair<SimulationVariable, std::vector<double> > > region_data_map;
          region_data_map::iterator region_data_it;

          region_data.insert( std::make_pair("Doping", std::make_pair(region->get_variable("na", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Doping", std::make_pair(region->get_variable("nd", POINT_CENTER), std::vector<double>())) );

          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("electron", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("hole", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("potential", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("temperature", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("elec_temperature", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("hole_temperature", POINT_CENTER), std::vector<double>())) );

          if(sigle)
            region_data.insert( std::make_pair("Mole", std::make_pair(region->get_variable("mole_x", POINT_CENTER), std::vector<double>())) );

          if(complex)
          {
            region_data.insert( std::make_pair("Mole", std::make_pair(region->get_variable("mole_x", POINT_CENTER), std::vector<double>())) );
            region_data.insert( std::make_pair("Mole", std::make_pair(region->get_variable("mole_y", POINT_CENTER), std::vector<double>())) );
          }

          std::vector<SimulationVariable> custom_variable;
          region->get_user_defined_variable(POINT_CENTER, SCALAR, custom_variable);
          for(unsigned int n=0; n<custom_variable.size(); ++n)
            region_data.insert( std::make_pair("Custom", std::make_pair(custom_variable[n], std::vector<double>())) );


          SimulationRegion::const_processor_node_iterator node_it = region->on_processor_nodes_begin();
          SimulationRegion::const_processor_node_iterator node_it_end = region->on_processor_nodes_end();
          for(; node_it!=node_it_end; ++node_it)
          {
            const FVM_Node * fvm_node = (*node_it);
            const FVM_NodeData * node_data = fvm_node->node_data();
            genius_assert(node_data);

            region_node_id.push_back(fvm_node->root_node()->id());

            for( region_data_it = region_data.begin(); region_data_it != region_data.end(); ++region_data_it)
            {
              const SimulationVariable & variable = region_data_it->second.first;
              region_data_it->second.second.push_back( node_data->data<double>( variable.variable_index ) / variable.variable_unit );
            }
          }

          // synchronization data with other processor
          Parallel::gather(0, region_node_id);

          for( region_data_it = region_data.begin(); region_data_it != region_data.end(); ++region_data_it)
            Parallel::gather(0, region_data_it->second.second);

          if( Genius::processor_id() == 0)
          {
            std::vector<unsigned int> region_node_local_id;
            for(unsigned int n=0; n<region_node_id.size(); ++n)
              region_node_local_id.push_back( node_id_to_region_node_id[region_node_id[n]]-1 );

            for( region_data_it = region_data.begin(); region_data_it != region_data.end(); )
            {
              const std::string & sol =  region_data_it->first;
              genius_assert(!cg_sol_write  (fn, B, Z, sol.c_str(), Vertex, &SOL));
              std::string path = std::string("/") + base_name + "/" + zone_name + "/" + sol;
              std::pair <region_data_map::iterator, region_data_map::iterator>  bounds = region_data.equal_range(sol);
              while (bounds.first != bounds.second)
              {
                const std::string & variable =  bounds.first->second.first.variable_name;
                std::string variable_unit = variable + ".unit";
                std::string variable_unit_string = bounds.first->second.first.variable_unit_string;
                genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, variable.c_str(),  &(_sort_it(bounds.first->second.second, region_node_local_id)[0]),   &F));
                genius_assert(!cg_gopath(fn, path.c_str()));
                genius_assert(!cg_descriptor_write(variable_unit.c_str(), variable_unit_string.c_str()));
                ++bounds.first;
              }
              region_data_it = bounds.second;
            }
          }

          break;
        }


        case InsulatorRegion     :
        {
          std::vector<unsigned int> region_node_id;

          std::multimap< std::string, std::pair<SimulationVariable, std::vector<double> > >  region_data;
          typedef std::multimap< std::string, std::pair<SimulationVariable, std::vector<double> > > region_data_map;
          region_data_map::iterator region_data_it;

          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("potential", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("electron", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("hole", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("temperature", POINT_CENTER), std::vector<double>())) );


          std::vector<SimulationVariable> custom_variable;
          region->get_user_defined_variable(POINT_CENTER, SCALAR, custom_variable);
          for(unsigned int n=0; n<custom_variable.size(); ++n)
            region_data.insert( std::make_pair("Custom", std::make_pair(custom_variable[n], std::vector<double>())) );

          SimulationRegion::const_processor_node_iterator node_it = region->on_processor_nodes_begin();
          SimulationRegion::const_processor_node_iterator node_it_end = region->on_processor_nodes_end();
          for(; node_it!=node_it_end; ++node_it)
          {
            const FVM_Node * fvm_node = (*node_it);
            const FVM_NodeData * node_data = fvm_node->node_data();

            region_node_id.push_back(fvm_node->root_node()->id());

            for( region_data_it = region_data.begin(); region_data_it != region_data.end(); ++region_data_it)
            {
              const SimulationVariable & variable = region_data_it->second.first;
              region_data_it->second.second.push_back( node_data->data<double>( variable.variable_index ) / variable.variable_unit );
            }
          }

          // synchronization data with other processor
          Parallel::gather(0, region_node_id);
          for( region_data_it = region_data.begin(); region_data_it != region_data.end(); ++region_data_it)
            Parallel::gather(0, region_data_it->second.second);

          if( Genius::processor_id() == 0)
          {
            std::vector<unsigned int> region_node_local_id;
            for(unsigned int n=0; n<region_node_id.size(); ++n)
              region_node_local_id.push_back( node_id_to_region_node_id[region_node_id[n]]-1 );

            for( region_data_it = region_data.begin(); region_data_it != region_data.end(); )
            {
              const std::string & sol =  region_data_it->first;
              genius_assert(!cg_sol_write  (fn, B, Z, sol.c_str(), Vertex, &SOL));
              std::string path = std::string("/") + base_name + "/" + zone_name + "/" + sol;
              std::pair <region_data_map::iterator, region_data_map::iterator>  bounds = region_data.equal_range(sol);
              while (bounds.first != bounds.second)
              {
                const std::string & variable =  bounds.first->second.first.variable_name;
                std::string variable_unit = variable + ".unit";
                std::string variable_unit_string = bounds.first->second.first.variable_unit_string;
                genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, variable.c_str(),  &(_sort_it(bounds.first->second.second, region_node_local_id)[0]),   &F));
                genius_assert(!cg_gopath(fn, path.c_str()));
                genius_assert(!cg_descriptor_write(variable_unit.c_str(), variable_unit_string.c_str()));
                ++bounds.first;
              }
              region_data_it = bounds.second;
            }
          }
          break;
        }
        case ElectrodeRegion     :
        {
          std::vector<unsigned int> region_node_id;

          std::multimap< std::string, std::pair<SimulationVariable, std::vector<double> > >  region_data;
          typedef std::multimap< std::string, std::pair<SimulationVariable, std::vector<double> > > region_data_map;
          region_data_map::iterator region_data_it;

          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("potential", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("temperature", POINT_CENTER), std::vector<double>())) );


          std::vector<SimulationVariable> custom_variable;
          region->get_user_defined_variable(POINT_CENTER, SCALAR, custom_variable);
          for(unsigned int n=0; n<custom_variable.size(); ++n)
            region_data.insert( std::make_pair("Custom", std::make_pair(custom_variable[n], std::vector<double>())) );

          SimulationRegion::const_processor_node_iterator node_it = region->on_processor_nodes_begin();
          SimulationRegion::const_processor_node_iterator node_it_end = region->on_processor_nodes_end();
          for(; node_it!=node_it_end; ++node_it)
          {
            const FVM_Node * fvm_node = (*node_it);
            const FVM_NodeData * node_data = fvm_node->node_data();

            region_node_id.push_back(fvm_node->root_node()->id());

            for( region_data_it = region_data.begin(); region_data_it != region_data.end(); ++region_data_it)
            {
              const SimulationVariable & variable = region_data_it->second.first;
              region_data_it->second.second.push_back( node_data->data<double>( variable.variable_index ) / variable.variable_unit );
            }
          }

          // synchronization data with other processor
          Parallel::gather(0, region_node_id);
          for( region_data_it = region_data.begin(); region_data_it != region_data.end(); ++region_data_it)
            Parallel::gather(0, region_data_it->second.second);

          if( Genius::processor_id() == 0)
          {
            std::vector<unsigned int> region_node_local_id;
            for(unsigned int n=0; n<region_node_id.size(); ++n)
              region_node_local_id.push_back( node_id_to_region_node_id[region_node_id[n]]-1 );

            for( region_data_it = region_data.begin(); region_data_it != region_data.end(); )
            {
              const std::string & sol =  region_data_it->first;
              genius_assert(!cg_sol_write  (fn, B, Z, sol.c_str(), Vertex, &SOL));
              std::string path = std::string("/") + base_name + "/" + zone_name + "/" + sol;
              std::pair <region_data_map::iterator, region_data_map::iterator>  bounds = region_data.equal_range(sol);
              while (bounds.first != bounds.second)
              {
                const std::string & variable =  bounds.first->second.first.variable_name;
                std::string variable_unit = variable + ".unit";
                std::string variable_unit_string = bounds.first->second.first.variable_unit_string;
                genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, variable.c_str(),  &(_sort_it(bounds.first->second.second, region_node_local_id)[0]),   &F));
                genius_assert(!cg_gopath(fn, path.c_str()));
                genius_assert(!cg_descriptor_write(variable_unit.c_str(), variable_unit_string.c_str()));
                ++bounds.first;
              }
              region_data_it = bounds.second;
            }
          }
          break;
        }
        case MetalRegion         :
        {
          std::vector<unsigned int> region_node_id;

          std::multimap< std::string, std::pair<SimulationVariable, std::vector<double> > >  region_data;
          typedef std::multimap< std::string, std::pair<SimulationVariable, std::vector<double> > > region_data_map;
          region_data_map::iterator region_data_it;

          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("potential", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("electron", POINT_CENTER), std::vector<double>())) );
          region_data.insert( std::make_pair("Solution", std::make_pair(region->get_variable("temperature", POINT_CENTER), std::vector<double>())) );


          std::vector<SimulationVariable> custom_variable;
          region->get_user_defined_variable(POINT_CENTER, SCALAR, custom_variable);
          for(unsigned int n=0; n<custom_variable.size(); ++n)
            region_data.insert( std::make_pair("Custom", std::make_pair(custom_variable[n], std::vector<double>())) );

          SimulationRegion::const_processor_node_iterator node_it = region->on_processor_nodes_begin();
          SimulationRegion::const_processor_node_iterator node_it_end = region->on_processor_nodes_end();
          for(; node_it!=node_it_end; ++node_it)
          {
            const FVM_Node * fvm_node = (*node_it);
            const FVM_NodeData * node_data = fvm_node->node_data();

            region_node_id.push_back(fvm_node->root_node()->id());

            for( region_data_it = region_data.begin(); region_data_it != region_data.end(); ++region_data_it)
            {
              const SimulationVariable & variable = region_data_it->second.first;
              region_data_it->second.second.push_back( node_data->data<double>( variable.variable_index ) / variable.variable_unit );
            }
          }

          // synchronization data with other processor
          Parallel::gather(0, region_node_id);
          for( region_data_it = region_data.begin(); region_data_it != region_data.end(); ++region_data_it)
            Parallel::gather(0, region_data_it->second.second);

          if( Genius::processor_id() == 0)
          {
            std::vector<unsigned int> region_node_local_id;
            for(unsigned int n=0; n<region_node_id.size(); ++n)
              region_node_local_id.push_back( node_id_to_region_node_id[region_node_id[n]]-1 );

            for( region_data_it = region_data.begin(); region_data_it != region_data.end(); )
            {
              const std::string & sol =  region_data_it->first;
              genius_assert(!cg_sol_write  (fn, B, Z, sol.c_str(), Vertex, &SOL));
              std::string path = std::string("/") + base_name + "/" + zone_name + "/" + sol;
              std::pair <region_data_map::iterator, region_data_map::iterator>  bounds = region_data.equal_range(sol);
              while (bounds.first != bounds.second)
              {
                const std::string & variable =  bounds.first->second.first.variable_name;
                std::string variable_unit = variable + ".unit";
                std::string variable_unit_string = bounds.first->second.first.variable_unit_string;
                genius_assert(!cg_field_write(fn, B, Z, SOL, RealDouble, variable.c_str(),  &(_sort_it(bounds.first->second.second, region_node_local_id)[0]),   &F));
                genius_assert(!cg_gopath(fn, path.c_str()));
                genius_assert(!cg_descriptor_write(variable_unit.c_str(), variable_unit_string.c_str()));
                ++bounds.first;
              }
              region_data_it = bounds.second;
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


  // write global boundary conditions, including interconnect and charge boundary
  if( Genius::processor_id() == 0)
  {
    genius_assert(!cg_goto(fn,B,"end"));
    assert(!cg_user_data_write ("ExtraBoundaryInfo"));
    assert(!cg_goto(fn, B, "UserDefinedData_t", 1, "end"));
    const BoundaryConditionCollector * bcs = system.get_bcs();
    for(unsigned int n=0; n<bcs->n_bcs(); ++n)
    {
      const BoundaryCondition * bc = bcs->get_bc(n);

      // interconnect bc
      if( bc->bc_type() == InterConnect || bc->bc_type() == ChargeIntegral)
      {
        genius_assert(!cg_descriptor_write(bc->label().c_str(), bc->boundary_condition_in_string().c_str() ));
      }
    }
  }

  // close CGNS file
  if( Genius::processor_id() == 0)
    cg_close(fn);

}


std::vector<double> & CGNSIO::_sort_it (std::vector<double> & x, const std::vector<unsigned int > &id)
{
  std::vector<double> xx(x) ;

  for(unsigned int i=0; i<id.size(); i++)
    x[id[i]] = xx[i];

  return x;
}



