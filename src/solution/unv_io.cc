// $Id: unv_io.C 4951 2011-11-10 21:14:30Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// C++ includes
#include <cassert>
#include <iomanip>
#include <cstdio>   // for std::sprintf
#include <algorithm> // for std::sort
#include <fstream>

// Local includes
#include "unv_io.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "mesh_communication.h"
#include "simulation_region.h"
#include "parallel.h"


using PhysicalUnit::mm;


//-----------------------------------------------------------------------------
// UNVIO class static members
const std::string UNVIO::_label_dataset_nodes    = "2411";
const std::string UNVIO::_label_dataset_elements = "2412";
const std::string UNVIO::_label_dataset_groups = "2467";


// ------------------------------------------------------------
// UNVIO class members
void UNVIO::clear ()
{
  /*
   * Initialize these to dummy values
   */
  this->_n_nodes     = 0;
  this->_n_elements  = 0;
  this->_need_D_to_e = true;

  this->_assign_nodes.clear();
  this->_ds_position.clear();

  this->_groups.clear();
  this->_element_in_group.clear();
  this->_regions.clear();
}


void UNVIO::read (const std::string& file_name)
{
  std::ifstream in_stream (file_name.c_str());
  this->read_implementation (in_stream);

  this->read_post_process();
}


void UNVIO::read_implementation (std::istream& in_stream)
{
  // clear everything, so that
  // we can start from scratch
  this->clear ();


  // Note that we read this file
  // @e twice.  First time to
  // detect the number of nodes
  // and elements (and possible
  // conversion tasks like D_to_e)
  // and the order of datasets
  // (nodes first, then elements,
  // or the other way around),
  // and second to do the actual
  // read.
  std::vector<std::string> order_of_datasets;
  order_of_datasets.reserve(2);

  {
    // the first time we read the file,
    // merely to obtain overall info
    if ( !in_stream.good() )
    {
      std::cerr << "ERROR: Input file not good."
      << std::endl;
      genius_error();
    }


    // Count nodes and elements, then let
    // other methods read the element and
    // node data.  Also remember which
    // dataset comes first: nodes or elements


    //    bool reached_eof = false;
    bool found_node  = false;
    bool found_elem  = false;
    bool found_group  = false;

    std::string olds, news;

    while (in_stream.good())
    {
      in_stream >> olds >> news;


      // a "-1" followed by a number means the beginning of a dataset
      // stop combing at the end of the file
      while ( ((olds != "-1") || (news == "-1") ) && !in_stream.eof() )
      {
        olds = news;
        in_stream >> news;
      }



      // if beginning of dataset, buffer it in
      // temp_buffer, if desired
      if (news == _label_dataset_nodes)
      {
        found_node = true;
        order_of_datasets.push_back (_label_dataset_nodes);
        this->count_nodes (in_stream);
      }

      else if (news == _label_dataset_elements)
      {
        found_elem = true;
        order_of_datasets.push_back (_label_dataset_elements);
        this->count_elements (in_stream);
      }
      else if (news == _label_dataset_groups)
      {
        found_group  = true;
        this->read_groups (in_stream);
      }



      if(found_node && found_elem && found_group) break;
    }


    // Here we should better have found
    // the datasets for nodes and elements,
    // otherwise the unv files is bad!
    if (!found_elem)
    {
      std::cerr << "ERROR: Could not find unv element set 2412!" << std::endl;
      genius_error();
    }

    if (!found_node)
    {
      std::cerr << "ERROR: Could not find unv node set 2411!" << std::endl;
      genius_error();
    }

    if (!found_group)
    {
      std::cerr << "ERROR: Could not find unv group set 2467!" << std::endl;
      genius_error();
    }


    // Don't close, just seek to the beginning
    in_stream.seekg(0, std::ios::beg);

    if (!in_stream.good() )
    {
      std::cerr << "ERROR: Cannot re-read input file."
      << std::endl;
      genius_error();
    }
  }


  // We finished scanning the file,
  // and our member data
  // \p this->_n_nodes,
  // \p this->_n_elements,
  // \p this->_need_D_to_e
  // should be properly initialized.
  {
    // Read the datasets in the order that
    // we already know

    for (unsigned int ds=0; ds < order_of_datasets.size(); ds++)
    {
      if (order_of_datasets[ds] == _label_dataset_nodes)
        this->node_in    (in_stream);
      else if (order_of_datasets[ds] == _label_dataset_elements)
        this->element_in (in_stream);
      else
        genius_error();
    }

  }


  // save memory
  this->_assign_nodes.clear();
  this->_ds_position.clear();

}



void UNVIO::build_boundary()
{
  std::map<std::string, std::pair<short int, bool> > bd_map;
  typedef std::map<std::string, std::pair<short int, bool> >::iterator Bd_It;

  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  //build neighbor information for boundary element. then elem->neighbor() is functional
  mesh.find_neighbors();
  //generate all the bounrary side with id -1
  //because I don't know if GMSH has marked all the side on.
  mesh.generate_boundary_info(-1);

  //mesh.boundary_info->print_info();

  //classfy boundary label
  std::vector<unsigned int>       elems;
  std::vector<unsigned short int> sides;
  std::vector<short int>          bds;

  // get all the boundary element
  mesh.boundary_info->build_side_list (elems, sides, bds);


  for (size_t nbd=0; nbd<elems.size(); nbd++ )
  {
    // get the element which has boundary/interface side
    const Elem* elem = mesh.elem(elems[nbd]);
    short int bd_index = bds[nbd];

    // face already has label (bd_id = -1 for labelless boundary side)
    if(bd_index >= 0 )  continue;

    //is it an interface side
    if( elem->neighbor(sides[nbd])!=NULL )
    {
      // the element and its neighbor should in diffetent subdomain
      unsigned int sbd_id1 = elem->subdomain_id();
      unsigned int sbd_id2 = elem->neighbor(sides[nbd])->subdomain_id();

      // delete the overkilled boundary side
      if (sbd_id1 == sbd_id2)
      {
        mesh.boundary_info->remove(elem, sides[nbd]);
        continue;
      }

      // the side should be an interface side
      genius_assert(elem->on_interface());
      genius_assert(elem->neighbor(sides[nbd])->on_interface());

      //remove the pair-element from boundary
      mesh.boundary_info->remove(elem, sides[nbd]);
      mesh.boundary_info->remove(elem->neighbor(sides[nbd]),
                                 elem->neighbor(sides[nbd])->which_neighbor_am_i(elem));

      // build the label for the interface, which has the form of RegionLabel1_to_RegionLabel2,
      // the two region is alpha ordered.
      std::string bd_label;
      std::string region1 = mesh.subdomain_label_by_id(sbd_id1);
      std::string region2 = mesh.subdomain_label_by_id(sbd_id2);
      if( region1 < region2)
        bd_label = region1 + "_to_" + region2;
      else
        bd_label = region2 + "_to_" + region1;


      // if the label already exist
      if( bd_map.find(bd_label) != bd_map.end() )
        bd_index = (*bd_map.find(bd_label)).second.first;
      else
      {
        //else, increase bd_index, insert it into bd_map
        bd_index = bd_map.size() + 1024;
        bd_map[bd_label] = std::make_pair(bd_index, false);
      }

      // add pair-element to boundary with new bd_index
      mesh.boundary_info->add_side(elem, sides[nbd], bd_index);
      mesh.boundary_info->add_side(elem->neighbor(sides[nbd]),
                                   elem->neighbor(sides[nbd])->which_neighbor_am_i(elem),
                                   bd_index);
    }
    // a boundary side
    else
    {
      unsigned int sbd_id = elem->subdomain_id();
      mesh.boundary_info->remove(elem, sides[nbd]);
      std::string bd_label = mesh.subdomain_label_by_id(sbd_id) + "_Neumann";

      // if the label already exist
      if( bd_map.find(bd_label) != bd_map.end() )
        bd_index = (*bd_map.find(bd_label)).second.first;
      else
      {
        //else, increase bd_index, insert it into bd_map
        bd_index = bd_map.size() + 1024;
        bd_map[bd_label] = std::make_pair(bd_index, false);
      }

      // add pair-element to boundary with new bd_index
      mesh.boundary_info->add_side(elem, sides[nbd], bd_index);
    }

  }


  // after the previous remove and insert operation, the number of boundary id maybe changed. renumber it here.
  mesh.boundary_info->rebuild_ids();

  //write down bd labels
  Bd_It bd_it = bd_map.begin();
  for(; bd_it != bd_map.end(); ++bd_it)
  {
    mesh.boundary_info->set_label_to_id( (*bd_it).second.first, (*bd_it).first, (*bd_it).second.second);
  }


}




void UNVIO::read_post_process()
{
  //set mesh subdomain
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  mesh.set_n_subdomains() = _regions.size();
  std::map<unsigned int, unsigned int>::const_iterator it = _regions.begin();
  for(; it != _regions.end(); ++it)
  {
    unsigned int group_id = it->first;
    unsigned int subdomain_id = it->second;
    std::string region_name = _groups[group_id];
    mesh.set_subdomain_label(subdomain_id,  region_name);

    // material name should be contained in the region name as region__material
    std::string material="Air";
    size_t pos = region_name.find("__");
    if(pos != std::string::npos)
      material = region_name.substr(pos+2);

    mesh.set_subdomain_material(subdomain_id, material);
  }

  build_boundary();

  // magic number, for 3D mesh, should > 2008
  if(_n_dim == 2)
    mesh.magic_num() = 1712;
  if(_n_dim == 3)
    mesh.magic_num() = 3712;

  // broadcast mesh to all the processor
  MeshCommunication mesh_comm;
  mesh_comm.broadcast(mesh);

  // build simulation system
  system.build_simulation_system();
  system.sync_print_info();


  // init each region
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);
    region->init(system.T_external());
  }

  system.init_region_post_process();

}


void UNVIO::count_nodes (std::istream& in_file)
{
  // if this->_n_nodes is not 0 the dataset
  // has already been scanned
  if (this->_n_nodes != 0)
  {
    std::cerr << "Error: Trying to scan nodes twice!"
    << std::endl;
    genius_error();
  }


  // Read from file, count nodes,
  // check if floats have to be converted
  std::string data;

  in_file >> data; // read the first node label


  if (data == "-1")
  {
    std::cerr << "ERROR: Bad, already reached end of dataset before even starting to read nodes!"
    << std::endl;
    genius_error();
  }


  // ignore the misc data for this node
  in_file.ignore(256,'\n');



  // Now we are there to verify whether we need
  // to convert from D to e or not
  in_file >> data;

  // When this "data" contains a "D", then
  // we have to convert each and every float...
  // But also assume when _this_ specific
  // line does not contain a "D", then the
  // other lines won't, too.
  {
    // #ifdef __HP_aCC
    //     // Use an "int" instead of unsigned int,
    //     // otherwise HP aCC may crash!
    //     const int position          = data.find("D",6);
    // #else
    //     const unsigned int position = data.find("D",6);
    // #endif
    std::string::size_type position = data.find("D",6);

    if (position!=std::string::npos) // npos means no position
    {
      this->_need_D_to_e = true;
    }
    else
      this->_need_D_to_e = false;
  }

  // read the remaining two coordinates
  in_file >> data;
  in_file >> data;


  // this was our first node
  this->_n_nodes++;



  // proceed _counting_ the remaining
  // nodes.
  while (in_file.good())
  {
    // read the node label
    in_file >> data;

    if (data == "-1")
      // end of dataset is reached
      break;

    // ignore the remaining data (coord_sys_label, color etc)
    in_file.ignore (256, '\n');
    // ignore the coordinates
    in_file.ignore (256, '\n');

    this->_n_nodes++;
  }


  if (in_file.eof())
  {
    std::cerr << "ERROR: File ended before end of node dataset!"  << std::endl;
    genius_error();
  }

}






void UNVIO::count_elements (std::istream& in_file)
{

  if (this->_n_elements != 0)
  {
    std::cerr << "Error: Trying to scan elements twice!" << std::endl;
    genius_error();
  }


  // Simply read the element
  // dataset for the @e only
  // purpose to count nodes!

  std::string data;
  unsigned int fe_id;

  while (!in_file.eof())
  {
    // read element label
    in_file >> data;

    // end of dataset?
    if (data == "-1")
      break;

    // read fe_id
    in_file >> fe_id;

    if(fe_id < 40 )_n_dim = std::max(1, _n_dim);//1D
    else if(fe_id < 100 )  _n_dim = std::max(2, _n_dim);//2D
    else _n_dim = std::max(3, _n_dim); //3D

    // Skip related data,
    // and node number list
    in_file.ignore (256,'\n');
    in_file.ignore (256,'\n');

    // For some elements the node numbers
    // are given more than one record

    // Rod
    if(fe_id == 11)
      in_file.ignore (256,'\n');

    // TET10 or QUAD9
    if (fe_id == 118 || fe_id == 300)
      in_file.ignore (256,'\n');

    // HEX20
    if (fe_id == 116)
    {
      in_file.ignore (256,'\n');
      in_file.ignore (256,'\n');
    }

    this->_n_elements++;
  }



  if (in_file.eof())
  {
    std::cerr << "ERROR: File ended before end of element dataset!"    << std::endl;
    genius_error();
  }

}



void UNVIO::node_in (std::istream& in_file)
{
  // adjust the \p istream to our position
  const bool ok = this->beginning_of_dataset(in_file, _label_dataset_nodes);

  if (!ok)
  {
    std::cerr << "ERROR: Could not find node dataset!" << std::endl;
    genius_error();
  }

  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  unsigned int node_lab;           // label of the node
  unsigned int exp_coord_sys_num,  // export coordinate system number       (not supported yet)
  disp_coord_sys_num, // displacement coordinate system number (not supported yet)
  color;              // color                                 (not supported yet)

  // allocate the correct amount
  // of memory for the node vector
  this->_assign_nodes.reserve (this->_n_nodes);


  // always 3 coordinates in the UNV file, no matter
  // which dimensionality libMesh is in
  //std::vector<Real> xyz (3);
  Point xyz;

  // depending on whether we have to convert each
  // coordinate (float), we offer two versions.
  // Note that \p count_nodes() already verified
  // whether this file uses "D" of "e"
  if (this->_need_D_to_e)
  {
    // ok, convert...
    std::string num_buf;

    for(unsigned int i=0; i<this->_n_nodes; i++)
    {
      assert (!in_file.eof());

      in_file >> node_lab                // read the node label
      >> exp_coord_sys_num       // (not supported yet)
      >> disp_coord_sys_num      // (not supported yet)
      >> color;                  // (not supported yet)

      // take care of the
      // floating-point data
      for (unsigned int d=0; d<3; d++)
      {
        in_file >> num_buf;
        xyz(d) = this->D_to_e (num_buf);
      }

      // set up the id map
      this->_assign_nodes.push_back (node_lab);
      mesh.add_point(xyz*mm,i);
    }
  }

  else
  {
    // very well, no need to convert anything,
    // just plain import.
    for (unsigned int i=0;i<this->_n_nodes;i++)
    {
      assert (!in_file.eof());

      in_file >> node_lab                // read the node label
      >> exp_coord_sys_num       // (not supported yet)
      >> disp_coord_sys_num      // (not supported yet)
      >> color                   // (not supported yet)
      >> xyz(0)                  // read x-coordinate
      >> xyz(1)                  // read y-coordinate
      >> xyz(2);                 // read z-coordinate

      // set up the id map
      this->_assign_nodes.push_back (node_lab);
      mesh.add_point(xyz*mm,i);
    }
  }

  // now we need to sort the _assign_nodes vector so we can
  // search it efficiently like a map
  std::sort (this->_assign_nodes.begin(),
             this->_assign_nodes.end());

}




void UNVIO::read_groups (std::istream& in_file)
{
  // Simply read the element
  // dataset for the @e only
  // purpose to count nodes!

  std::string data;
  std::string group_id;
  std::string group_name;
  unsigned int dummy;
  unsigned int group_count;

  while (!in_file.eof())
  {
    // read element label
    in_file >> data;

    // end of dataset?
    if (data == "-1")
      break;

    // read fe_id
    group_id = data;

    // dummy data
    for(unsigned int i=0; i<6; i++) in_file >> dummy;

    // how many items
    in_file >> group_count;

    // group name
    in_file >> group_name;

    // items
    for(unsigned int c=0; c<group_count; c++)
    {
      // Field 1       -- entity type code
      // Field 2       -- entity tag
      // Field 3       -- entity node leaf id.
      // Field 4       -- entity component/ ham id.

      unsigned int type,tag,a,b;
      in_file >> type >> tag >> a >> b;
      _element_in_group[tag] =  this->_n_groups;
    }

    _groups.push_back(group_name);

    this->_n_groups++;
  }


  if (in_file.eof())
  {
    std::cerr << "ERROR: File ended before end of group dataset!"    << std::endl;
    genius_error();
  }


}








void UNVIO::element_in (std::istream& in_file)
{
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  // adjust the \p istream to our
  // position
  const bool ok = this->beginning_of_dataset(in_file, _label_dataset_elements);

  if (!ok)
  {
    std::cerr << "ERROR: Could not find element dataset!" << std::endl;
    genius_error();
  }


  unsigned int      element_lab,       // element label (not supported yet)
  n_nodes;           // number of nodes on element
  unsigned long int fe_descriptor_id,  // FE descriptor id
  phys_prop_tab_num, // physical property table number (not supported yet)
  mat_prop_tab_num,  // material property table number (not supported yet)
  color;             // color (not supported yet)


  // vector that temporarily holds the node labels defining element
  std::vector<unsigned int> node_labels (21);


  // vector that assigns element nodes to their correct position
  // for example:
  // 44:plane stress      | QUAD4
  // linear quadrilateral |
  // position in UNV-file | position in libmesh
  // assign_elem_node[1]   = 0
  // assign_elem_node[2]   = 3
  // assign_elem_node[3]   = 2
  // assign_elem_node[4]   = 1
  //
  // UNV is 1-based, we leave the 0th element of the vectors unused in order
  // to prevent confusion, this way we can store elements with up to 20 nodes
  unsigned int assign_elem_nodes[21];


  // Get the beginning and end of the _assign_nodes vector
  // to eliminate repeated function calls
  const std::vector<unsigned int>::const_iterator it_begin =
    this->_assign_nodes.begin();

  const std::vector<unsigned int>::const_iterator it_end   =
    this->_assign_nodes.end();



  // read from the virtual file
  for (unsigned int i=0; i<this->_n_elements; i++)
  {
    in_file >> element_lab             // read element label
    >> fe_descriptor_id        // read FE descriptor id
    >> phys_prop_tab_num       // (not supported yet)
    >> mat_prop_tab_num        // (not supported yet)
    >> color                   // (not supported yet)
    >> n_nodes;                // read number of nodes on element

    if(fe_descriptor_id > 40)
    {
      // 2d/3d elem
      for (unsigned int j=1; j<=n_nodes; j++)
        in_file >> node_labels[j];       // read node labels
    }
    else
    {
      // 1d elem
      int d1, d2, d3;
      in_file >> d1 >> d2 >> d3;
      for (unsigned int j=1; j<=n_nodes; j++)
        in_file >> node_labels[j];       // read node labels
    }

    Elem* elem = NULL;                 // element pointer

    switch (fe_descriptor_id)
    {
    case 11: // rod
      {
        elem = Elem::build(EDGE2).release();  // create new element

        assign_elem_nodes[1]=0;
        assign_elem_nodes[2]=1;
        break;
      }
    case 41: // Plane Stress Linear Triangle
    case 91: // Thin Shell   Linear Triangle
      {
        elem = Elem::build(TRI3).release();  // create new element

        assign_elem_nodes[1]=0;
        assign_elem_nodes[2]=2;
        assign_elem_nodes[3]=1;
        break;
      }

    case 42: // Plane Stress Quadratic Triangle
    case 92: // Thin Shell   Quadratic Triangle
      {
        std::cerr << "ERROR: UNV-element type 43: Quadratic Triangle"
        << " not supported."
        << std::endl;
        genius_error();
        break;
      }

    case 43: // Plane Stress Cubic Triangle
      {
        std::cerr << "ERROR: UNV-element type 43: Plane Stress Cubic Triangle"
        << " not supported."
        << std::endl;
        genius_error();
        break;
      }

    case 44: // Plane Stress Linear Quadrilateral
    case 94: // Thin Shell   Linear Quadrilateral
      {
        elem = Elem::build(QUAD4).release(); // create new element

        assign_elem_nodes[1]=0;
        assign_elem_nodes[2]=3;
        assign_elem_nodes[3]=2;
        assign_elem_nodes[4]=1;
        break;
      }

    case 45: // Plane Stress Quadratic Quadrilateral
    case 95: // Thin Shell   Quadratic Quadrilateral
      {
        std::cerr << "ERROR: UNV-element type 46: Plane Stress Quadratic Quadrilateral"
        << " not supported."
        << std::endl;
        genius_error();
        break;
      }

    case 300: // Thin Shell   Quadratic Quadrilateral (nine nodes)
      {
        std::cerr << "ERROR: UNV-element type 46: PQuadratic Quadrilateral"
        << " not supported."
        << std::endl;
        genius_error();
        break;
      }

    case 46: // Plane Stress Cubic Quadrilateral
      {
        std::cerr << "ERROR: UNV-element type 46: Plane Stress Cubic Quadrilateral"
        << " not supported."
        << std::endl;
        genius_error();
        break;
      }

    case 111: // Solid Linear Tetrahedron
      {
        elem = Elem::build(TET4).release();  // create new element

        assign_elem_nodes[1]=0;
        assign_elem_nodes[2]=1;
        assign_elem_nodes[3]=2;
        assign_elem_nodes[4]=3;
        break;
      }

    case 112: // Solid Linear Prism
      {
        elem = Elem::build(PRISM6).release();  // create new element

        assign_elem_nodes[1]=0;
        assign_elem_nodes[2]=1;
        assign_elem_nodes[3]=2;
        assign_elem_nodes[4]=3;
        assign_elem_nodes[5]=4;
        assign_elem_nodes[6]=5;
        break;
      }

    case 115: // Solid Linear Brick
      {
        elem = Elem::build(HEX8).release();  // create new element

        assign_elem_nodes[1]=0;
        assign_elem_nodes[2]=4;
        assign_elem_nodes[3]=5;
        assign_elem_nodes[4]=1;
        assign_elem_nodes[5]=3;
        assign_elem_nodes[6]=7;
        assign_elem_nodes[7]=6;
        assign_elem_nodes[8]=2;
        break;
      }

    case 116: // Solid Quadratic Brick
      {
        std::cerr << "Error: UNV-element type 117: Solid Quadratic Brick"
        << " not supported."
        << std::endl;
        genius_error();
      }

    case 117: // Solid Cubic Brick
      {
        std::cerr << "Error: UNV-element type 117: Solid Cubic Brick"
        << " not supported."
        << std::endl;
        genius_error();
        break;
      }

    case 118: // Solid Quadratic Tetrahedron
      {
        std::cerr << "Error: UNV-element type 117: Solid Quadratic Tetrahedron"
        << " not supported."
        << std::endl;
        genius_error();
        break;
      }

    default: // Unrecognized element type
      {
        std::cerr << "ERROR: UNV-element type "
        << fe_descriptor_id
        << " not supported."
        << std::endl;
        genius_error();
        break;
      }
    }


    if(!elem) continue;


    // nodes are being stored in element
    for (unsigned int j=1; j<=n_nodes; j++)
    {
      // Find the position of node_labels[j] in the _assign_nodes vector.
      const std::pair<std::vector<unsigned int>::const_iterator,
      std::vector<unsigned int>::const_iterator>
      it = std::equal_range (it_begin,
                             it_end,
                             node_labels[j]);

      // it better be there, so libmesh_assert that it was found.
      assert (it.first  != it.second);
      assert (*(it.first) == node_labels[j]);

      // Now, the distance between this UNV id and the beginning of
      // the _assign_nodes vector will give us a unique id in the
      // range [0,n_nodes) that we can use for defining a contiguous
      // connectivity.
      const unsigned int assigned_node = std::distance (it_begin,
                                         it.first);

      // Make sure we didn't get an out-of-bounds id
      assert (assigned_node < this->_n_nodes);

      elem->set_node(assign_elem_nodes[j]) =  mesh.node_ptr(assigned_node);
    }

    // mesh element
    if(elem->dim() == _n_dim)
    {
      unsigned int group_id = 0;
      if(_element_in_group.find(element_lab)!= _element_in_group.end())
        group_id = _element_in_group.find(element_lab)->second;

      unsigned int subdomain_id = 0;
      if(_regions.find(group_id) == _regions.end())
      {
        subdomain_id = _regions.size();
        _regions[group_id]=subdomain_id;
      }
      else
      {
        subdomain_id = _regions.find(group_id)->second;
      }

      elem->subdomain_id() = subdomain_id;
      // add elem to the Mesh &
      // tell the MeshData object the foreign elem id
      // (note that mesh.add_elem() returns a pointer to the new element)
      mesh.add_elem(elem);
    }

    // boundary element
    if(elem->dim() < _n_dim)
    {
      delete elem;
    }

  }

}


bool UNVIO::beginning_of_dataset (std::istream& in_file,
                                  const std::string& ds_name) const
{
  assert (in_file.good());
  assert (!ds_name.empty());

  std::string olds, news;

  while (true)
  {
    in_file >> olds >> news;

    /*
    * a "-1" followed by a number means the beginning of a dataset
    * stop combing at the end of the file
    */
    while( ((olds != "-1") || (news == "-1") ) && !in_file.eof() )
    {
      olds = news;
      in_file >> news;
    }

    if (in_file.eof())
      return false;

    if (news == ds_name)
      return true;
  }

  // should never end up here
  return false;
}



Real UNVIO::D_to_e (std::string& number) const
{
  /* find "D" in string, start looking at
  * 6th element, to improve speed.
  * We dont expect a "D" earlier
  */


  const std::string::size_type position = number.find("D",6);

  assert (position != std::string::npos);
  number.replace(position, 1, "e");

  return std::atof (number.c_str());
}




