// $Id: gmsh_io.C 4951 2011-11-10 21:14:30Z roystgnr $

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

// This file was massively overhauled and extended by Martin L�thi, mluthi@tnoo.net

// C++ includes
#include <fstream>
#include <set>
#include <cstring> // std::memcpy, std::strncmp

// Local includes
//#include "genius_config.h"
#include "gmsh_io.h"
#include "elem.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "mesh_communication.h"
#include "simulation_region.h"
#include "parallel.h"

using PhysicalUnit::mm;

// anonymous namespace to hold local data
namespace
{

  /**
   * Defines a structure to hold boundary element information.
   *
   * We use a set because it keeps the nodes unique and ordered, and can be
   * easily compared to another set of nodes (the ones on the element side)
   */
  struct boundaryElementInfo
  {
    std::set<unsigned int> nodes;
    unsigned int id;
  };

  /**
   * Defines mapping from libMesh element types to Gmsh element types.
   */
  struct elementDefinition
  {
    std::string label;
    std::vector<unsigned int> nodes;
    ElemType type;
    unsigned int exptype;
    unsigned int dim;
    unsigned int nnodes;
  };


  // maps from a libMesh element type to the proper
  // Gmsh elementDefinition.  Placing the data structure
  // here in this anonymous namespace gives us the
  // benefits of a global variable without the nasty
  // side-effects
  std::map<ElemType, elementDefinition> eletypes_exp;
  std::map<unsigned int, elementDefinition> eletypes_imp;



  // ------------------------------------------------------------
  // helper function to initialize the eletypes map
  void init_eletypes ()
  {
    if (eletypes_exp.empty() && eletypes_imp.empty())
    {
      // This should happen only once.  The first time this method
      // is called the eletypes data struture will be empty, and
      // we will fill it.  Any subsequent calls will find an initialized
      // eletypes map and will do nothing.

      //==============================
      // setup the element definitions
      elementDefinition eledef;

      // use "swap trick" from Scott Meyer's "Effective STL" to initialize
      // eledef.nodes vector

      // POINT (only Gmsh)
      {
        eledef.exptype = 15;
        eledef.dim     = 0;
        eledef.nnodes  = 1;
        eledef.nodes.clear();

        // import only
        eletypes_imp[15] = eledef;
      }

      // EDGE2
      {
        eledef.type    = EDGE2;
        eledef.dim     = 1;
        eledef.nnodes  = 2;
        eledef.exptype = 1;
        eledef.nodes.clear();

        eletypes_exp[EDGE2] = eledef;
        eletypes_imp[1]     = eledef;
      }

      // EDGE3
      {
        eledef.type    = EDGE3;
        eledef.dim     = 1;
        eledef.nnodes  = 3;
        eledef.exptype = 8;
        eledef.nodes.clear();

        eletypes_exp[EDGE3] = eledef;
        eletypes_imp[8]     = eledef;
      }

      // TRI3
      {
        eledef.type    = TRI3;
        eledef.dim     = 2;
        eledef.nnodes  = 3;
        eledef.exptype = 2;
        eledef.nodes.clear();

        eletypes_exp[TRI3] = eledef;
        eletypes_imp[2] = eledef;
      }

      // TRI6
      {
        eledef.type    = TRI6;
        eledef.dim     = 2;
        eledef.nnodes  = 6;
        eledef.exptype = 9;
        eledef.nodes.clear();

        eletypes_exp[TRI6] = eledef;
        eletypes_imp[9]    = eledef;
      }

      // QUAD4
      {
        eledef.type    = QUAD4;
        eledef.dim     = 2;
        eledef.nnodes  = 4;
        eledef.exptype = 3;
        eledef.nodes.clear();

        eletypes_exp[QUAD4] = eledef;
        eletypes_imp[3]     = eledef;
      }

      // QUAD8
      // TODO: what should be done with this on writing?
      {
        eledef.type    = QUAD8;
        eledef.dim     = 2;
        eledef.nnodes  = 8;
        eledef.exptype = 100;
        const unsigned int nodes[] = {1,2,3,4,5,6,7,8};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[QUAD8] = eledef;
        eletypes_imp[10]    = eledef;
      }

      // QUAD9
      {
        eledef.type    = QUAD9;
        eledef.dim     = 2;
        eledef.nnodes  = 9;
        eledef.exptype = 10;
        eledef.nodes.clear();

        eletypes_exp[QUAD9] = eledef;
        eletypes_imp[10]    = eledef;
      }

      // HEX8
      {
        eledef.type    = HEX8;
        eledef.dim     = 3;
        eledef.nnodes  = 8;
        eledef.exptype = 5;
        eledef.nodes.clear();

        eletypes_exp[HEX8] = eledef;
        eletypes_imp[5]    = eledef;
      }

      // HEX20
      // TODO: what should be done with this on writing?
      {
        eledef.type    = HEX20;
        eledef.dim     = 3;
        eledef.nnodes  = 20;
        eledef.exptype = 101;
        const unsigned int nodes[] = {1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,16};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[HEX20] = eledef;
        eletypes_imp[12]    = eledef;
      }

      // HEX27
      {
        eledef.type    = HEX27;
        eledef.dim     = 3;
        eledef.nnodes  = 27;
        eledef.exptype = 12;
        const unsigned int nodes[] =
          {
            0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,
            15,16,19,17,18,20,21,24,22,23,25,26
          };
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[HEX27] = eledef;
        eletypes_imp[12]    = eledef;
      }

      // TET4
      {
        eledef.type    = TET4;
        eledef.dim     = 3;
        eledef.nnodes  = 4;
        eledef.exptype = 4;
        eledef.nodes.clear();

        eletypes_exp[TET4] = eledef;
        eletypes_imp[4]    = eledef;
      }

      // TET10
      {
        eledef.type    = TET10;
        eledef.dim     = 3;
        eledef.nnodes  = 10;
        eledef.exptype = 11;
        const unsigned int nodes[] = {0,1,2,3,4,5,6,7,9,8};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);
        eletypes_exp[TET10] = eledef;
        eletypes_imp[11]    = eledef;
      }

      // PRISM6
      {
        eledef.type    = PRISM6;
        eledef.dim     = 3;
        eledef.nnodes  = 6;
        eledef.exptype = 6;
        eledef.nodes.clear();

        eletypes_exp[PRISM6] = eledef;
        eletypes_imp[6]      = eledef;
      }

      // PRISM15
      // TODO: what should be done with this on writing?
      {
        eledef.type    = PRISM15;
        eledef.dim     = 3;
        eledef.nnodes  = 15;
        eledef.exptype = 103;
        eledef.nodes.clear();

        eletypes_exp[PRISM15] = eledef;
        eletypes_imp[13] = eledef;
      }

      // PRISM18
      {
        eledef.type    = PRISM18;
        eledef.dim     = 3;
        eledef.nnodes  = 18;
        eledef.exptype = 13;
        const unsigned int nodes[] =
          {
            0,1,2,3,4,5,6,8,9,7,10,11,
            12,14,13,15,17,16
          };
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[PRISM18] = eledef;
        eletypes_imp[13]      = eledef;
      }

      // PYRAMID5
      {
        eledef.type    = PYRAMID5;
        eledef.dim     = 3;
        eledef.nnodes  = 5;
        eledef.exptype = 7;
        eledef.nodes.clear();

        eletypes_exp[PYRAMID5] = eledef;
        eletypes_imp[7]        = eledef;
      }

      //==============================
    }
  }

} // end anonymous namespace


// ------------------------------------------------------------
// GmshIO  members
void GmshIO::read (const std::string& name)
{
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();
  mesh.clear();

  if(Genius::processor_id() == 0)
  {
    this->read_info (name+".info");
    this->read_mesh (name+".msh");
  }

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


void GmshIO::read_info (const std::string& name)
{
  std::ifstream in (name.c_str());
  if(!in.good())
  {
    std::cerr << "Open GMSH info file " << name << " failed." << "\n";
    genius_error();
  }

  const int  bufLen = 256;
  char       buf[bufLen+1];
  while (!in.eof())
  {
    in >> buf;
    if (!std::strncmp(buf,"$RegionInfo",11))
    {
      int physical;
      std::string name, material;
      in >> physical >> name >> material;
      in >> buf;
      region_info[physical]=std::make_pair(name, material);
    }
    else if (!std::strncmp(buf,"$BoundaryInfo",13))
    {
      int physical;
      std::string name;
      in >> physical >> name;
      in >> buf;
      boundary_info[physical]=name;
    }
  }
 
}


void GmshIO::read_mesh(const std::string& name)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  genius_assert(Genius::processor_id() == 0);

  std::ifstream in (name.c_str());
  if(!in.good())
  {
    std::cerr << "Open GMSH file " << name << " failed." << "\n";
    genius_error();
  }

  // initialize the map with element types
  init_eletypes();

  // clear any data in the mesh
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  const unsigned int dim = 3;

  // some variables
  const int  bufLen = 256;
  char       buf[bufLen+1];
  int        format=0, size=0;
  Real       version = 1.0;

  // map to hold the node numbers for translation
  // note the the nodes can be non-consecutive
  std::map<unsigned int, unsigned int> nodetrans;

  std::map<int, unsigned int> elem_physical_map;
  std::map<int, unsigned int> boundary_physical_map;

  {
    while (!in.eof())
    {
      in >> buf;

      if (!std::strncmp(buf,"$MeshFormat",11))
      {
        in >> version >> format >> size;
        if ((version != 2.0) && (version != 2.1) && (version != 2.2))
        {
          // Some notes on gmsh mesh versions:
          //
          // Mesh version 2.0 goes back as far as I know.  It's not explicitly
          // mentioned here: http://www.geuz.org/gmsh/doc/VERSIONS.txt
          //
          // As of gmsh-2.4.0:
          // bumped mesh version format to 2.1 (small change in the $PhysicalNames
          // section, where the group dimension is now required);
          // [Since we don't even parse the PhysicalNames section at the time
          //  of this writing, I don't think this change affects us.]
          std::cerr << "Error: Wrong msh file version " << version << "\n";
          genius_error();
        }
        if(format)
        {
          std::cerr << "Error: Unknown data format for mesh\n";
          genius_error();
        }
      }

      // read the node block
      else if (!std::strncmp(buf,"$NOD",4) ||
               !std::strncmp(buf,"$NOE",4) ||
               !std::strncmp(buf,"$Nodes",6)
              )
      {
        unsigned int numNodes = 0;
        in >> numNodes;
        mesh.reserve_nodes (numNodes);

        // read in the nodal coordinates and form points.
        Real x, y, z;
        unsigned int id;

        // add the nodal coordinates to the mesh
        for (unsigned int i=0; i<numNodes; ++i)
        {
          in >> id >> x >> y >> z;
          mesh.add_point (Point(x, y, z)*mm, i);
          nodetrans[id] = i;
        }
        // read the $ENDNOD delimiter
        in >> buf;
        
      }

      /**
       * Read the element block
       *
       * If the element dimension is smaller than the mesh dimension, this is a
       * boundary element and will be added to mesh.boundary_info.
       *
       * Because the elements might not yet exist, the sides are put on hold
       * until the elements are created, and inserted once reading elements is
       * finished
       */
      else if (!std::strncmp(buf,"$ELM",4) ||
               !std::strncmp(buf,"$Elements",9)
              )
      {
        unsigned int numElem = 0;
        std::vector< boundaryElementInfo > boundary_elem;

        // read how many elements are there, and reserve space in the mesh
        in >> numElem;
        mesh.reserve_elem (numElem);

        // read the elements
        unsigned int elem_id_counter = 0;
        for (unsigned int iel=0; iel<numElem; ++iel)
        {
          unsigned int id, type, physical=0, elementary=0,
          /* partition = 1,*/ nnodes, ntags;
          // note - partition was assigned but never used - BSK
          if(version <= 1.0)
          {
            in >> id >> type >> physical >> elementary >> nnodes;
          }
          else
          {
            in >> id >> type >> ntags;
            elementary = physical = /* partition = */ 1;
            for(unsigned int j = 0; j < ntags; j++)
            {
              int tag;
              in >> tag;
              if(j == 0)
                physical = tag;
              else if(j == 1)
                elementary = tag;
              // else if(j == 2)
              //  partition = tag;
              // ignore any other tags for now
            }
          }


          // consult the import element table which element to build
          const elementDefinition& eletype = eletypes_imp[type];
          nnodes = eletype.nnodes;

          // only elements that match the mesh dimension are added
          // if the element dimension is one less than dim, the nodes and
          // sides are added to the mesh.boundary_info
          if (eletype.dim == dim)
          {
            // add the elements to the mesh
            Elem* elem = Elem::build(eletype.type).release();
            elem->set_id(elem_id_counter);
            mesh.add_elem(elem);

            // different to iel, lower dimensional elems aren't added
            elem_id_counter++;

            // check number of nodes. We cannot do that for version 2.0
            if (version <= 1.0)
            {
              if (elem->n_nodes() != nnodes)
              {
                std::cerr << "Number of nodes for element " << id
                << " of type " << eletypes_imp[type].type
                << " (Gmsh type " << type
                << ") does not match Libmesh definition. "
                << "I expected " << elem->n_nodes()
                << " nodes, but got " << nnodes << "\n";
                genius_error();
              }
            }

            // add node pointers to the elements
            int nod = 0;
            // if there is a node translation table, use it
            if (eletype.nodes.size() > 0)
              for (unsigned int i=0; i<nnodes; i++)
              {
                in >> nod;
                elem->set_node(eletype.nodes[i]) = mesh.node_ptr(nodetrans[nod]);
              }
            else
            {
              for (unsigned int i=0; i<nnodes; i++)
              {
                in >> nod;
                elem->set_node(i) = mesh.node_ptr(nodetrans[nod]);
              }
            }

            // Finally, set the subdomain ID to physical
            if(elem_physical_map.find(physical) == elem_physical_map.end())
              elem_physical_map[physical]=elem_physical_map.size();
            elem->subdomain_id() = elem_physical_map[physical];
          } // if element.dim == dim
          // if this is a boundary
          else if (eletype.dim == dim-1)
          {
            /**
             * add the boundary element nodes to the set of nodes
             */

            boundaryElementInfo binfo;
            std::set<unsigned int>::iterator iter = binfo.nodes.begin();
            int nod = 0;
            for (unsigned int i=0; i<nnodes; i++)
            {
              in >> nod;
              mesh.boundary_info->add_node(nodetrans[nod], physical);
              binfo.nodes.insert(iter, nodetrans[nod]);
            }
            if(boundary_physical_map.find(physical) == boundary_physical_map.end())
              boundary_physical_map[physical]=boundary_physical_map.size();
            binfo.id = boundary_physical_map[physical];
            boundary_elem.push_back(binfo);
          }
          /**
           * If the element yet another dimension, just read in the nodes
           * and throw them away
           */
          else
          {
            static bool seen_high_dim_element = false;
            if (!seen_high_dim_element)
            {
              std::cerr << "Warning: can't load an element of dimension "
              << eletype.dim << " into a mesh of dimension "
              << dim << std::endl;
              seen_high_dim_element = true;
            }
            int nod = 0;
            for (unsigned int i=0; i<nnodes; i++)
              in >> nod;
          }
        }//element loop


        // read the $ENDELM delimiter
        in >> buf;

        
        /**
         * If any lower dimensional elements have been found in the file,
         * try to add them to the mesh.boundary_info as sides and nodes with
         * the respecitve id's (called "physical" in Gmsh).
         */
        if (boundary_elem.size() > 0)
        {
          // create a index of the boundary nodes to easily locate which
          // element might have that boundary
          std::map<unsigned int, std::vector<unsigned int> > node_index;
          for (unsigned int i=0; i<boundary_elem.size(); i++)
          {
            boundaryElementInfo binfo = boundary_elem[i];
            std::set<unsigned int>::iterator iter = binfo.nodes.begin();
            for (;iter!= binfo.nodes.end(); iter++)
              node_index[*iter].push_back(i);
          }

          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          // iterate over all elements and see which boundary element has
          // the same set of nodes as on of the boundary elements previously read
          for ( ; it != end; ++it)
          {
            const Elem* elem = *it;
            for (unsigned int s=0; s<elem->n_sides(); s++)
              if (elem->neighbor(s) == NULL)
              {
                AutoPtr<Elem> side (elem->build_side(s));
                std::set<unsigned int> side_nodes;
                std::set<unsigned int>::iterator iter = side_nodes.begin();

                // make a set with all nodes from this side
                // this allows for easy comparison
                for (unsigned int ns=0; ns<side->n_nodes(); ns++)
                  side_nodes.insert(iter, side->node(ns));

                // See whether one of the side node occurs in the list
                // of tagged nodes. If we would loop over all side
                // nodes, we would just get multiple hits, so taking
                // node 0 is enough to do the job
                unsigned int sn = side->node(0);
                if (node_index.count(sn) > 0)
                {
                  // Loop over all tagged ("physical") "sides" which
                  // contain the node sn (typically just 1 to
                  // three). For each of these the set of nodes is
                  // compared to the current element's side nodes
                  for (unsigned int n=0; n<node_index[sn].size(); n++)
                  {
                    unsigned int bidx = node_index[sn][n];
                    if (boundary_elem[bidx].nodes == side_nodes)
                    {
                      mesh.boundary_info->add_side(elem, s, boundary_elem[bidx].id);
                    }
                  }
                }
              } // if elem->neighbor(s) == NULL
          } // element loop
        } // if boundary_elem.size() > 0
      } // if $ELM

    } // while !in.eof()

  }

  // set mesh subdomain info


  mesh.set_n_subdomains() = elem_physical_map.size();

  std::vector<std::string> regions(mesh.n_subdomains());
  {
    std::map<int, unsigned int>::const_iterator it = elem_physical_map.begin();
    for(; it!=elem_physical_map.end(); ++it)
    {
      if(region_info.find(it->first) == region_info.end())
      {
        std::cerr << "Region with physical tag " << it->first << " do not have name/material info in the inf file" << "\n";
        genius_error();
      }
      const std::string & name = region_info.find(it->first)->second.first;
      const std::string & material = region_info.find(it->first)->second.second;

      mesh.set_subdomain_label(it->second, name );
      mesh.set_subdomain_material(it->second, material);

      regions[it->second] = name;
    }
  }

  // set mesh boundary info

  std::map<std::string, std::pair<short int, bool> > bd_map;
  typedef std::map<std::string, std::pair<short int, bool> >::iterator Bd_It;
  {
    std::map<int, unsigned int>::const_iterator it = boundary_physical_map.begin();
    for(; it!=boundary_physical_map.end(); ++it)
    {
      std::string name;
      if(boundary_info.find(it->first) == boundary_info.end())
      {
        std::stringstream ss;
        ss << it->first;
        name =  ss.str();
      }
      else
       name = boundary_info.find(it->first)->second;

      bd_map[name] = std::make_pair(it->second, true);
    }
  }

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
      std::string region1 = regions[sbd_id1];
      std::string region2 = regions[sbd_id2];
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
      std::string bd_label = regions[sbd_id] + "_Neumann";

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

  // magic number, for 3D mesh, should > 2008
  mesh.magic_num() = 3412;

}




