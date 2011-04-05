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

#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <iostream>
#include <fstream>

#include "dfise_block.h"
#include "dfise.h"

#include "config.h"
#ifdef CYGWIN
  #include <io.h>      // for windows _access function
#else
  #include <unistd.h>  // for POSIX access function
#endif

namespace DFISE
{


  //---------------------------------------------------
  // functions for class DFISE_MESH


  // avoid isatty() problem of Bison 2.3
#define YY_NEVER_INTERACTIVE 1
#ifdef CYGWIN
  #define YY_NO_UNISTD_H 1
#endif

#include "dfise_lex.yy.c"
#include "dfise_parser.tab.c"

  int DFISE_MESH::parse_dfise(const std::string & file)
  {
    std::string grid_file = file + ".grd";
    std::string data_file = file + ".dat";

    //check if we can fild the grd file
#ifdef CYGWIN
    if ( _access( (char *)grid_file.c_str(),  04 ) == -1 )
#else
    if (  access( grid_file.c_str(),  R_OK ) == -1 )
#endif
    {
      // still no find?
      std::cout<<"No DFISE grid file " << std::string(file + ".grd") << " can be found."<< std::endl;
      return 1;
    }


    //check if we can fild the data file
#ifdef CYGWIN
    if ( _access( (char *)data_file.c_str(),  04 ) == -1 )
#else
    if (  access( data_file.c_str(),  R_OK ) == -1 )
#endif
    {
      std::cout<<"No DFISE data file " << data_file << " can be found."<< std::endl;
      return 1;
    }


    //parse grid file
    if(parse_dfise_grid_file(grid_file)) return 1;

    //parse dataset file
    if(parse_dfise_dataset_file(data_file)) return 1;

    return 0;
  }



  void DFISE_MESH::write_dfise(const std::string & file) const
  {
    std::string grid_file = file + ".grd";
    std::string data_file = file + ".dat";

    // write grid file
    std::cout<<"  Writing DF-ISE grid file " << grid_file << "..."<< std::endl;
    std::ofstream grid_out(grid_file.c_str());
    grid_out<<"DF-ISE text\n"<<std::endl;
    grid_info.print(grid_out);
    grid.print(grid_out);
    grid_out.close();

    // write data file
    std::cout<<"  Writing DF-ISE dataset file " << data_file << "..."<< std::endl;
    std::ofstream data_out(data_file.c_str());
    data_out<<"DF-ISE text\n"<<std::endl;
    data_info.print(data_out);
    data_out<<"Data {"<<std::endl;
    for(unsigned int n=0; n<data_sets.size(); ++n)
      data_sets[n]->print(data_out);
    data_out<<"}"<<std::endl;
    data_out.close();
  }



  int DFISE_MESH::parse_dfise_grid_file(const std::string & grid_file)
  {
    std::cout<<"  Reading DF-ISE grid file " << grid_file << "..."<< std::endl;

    // top block
    BLOCK *block= new BLOCK;

    yyin = fopen(grid_file.c_str(), "r");
    assert( yyin != NULL );
    assert(!yyparse(block));
    fclose(yyin);
    YY_FLUSH_BUFFER;

    for(unsigned int n=0; n<block->n_sub_blocks(); ++n)
    {
      BLOCK * sub = block->get_sub_block(n);
      if(sub->keyword()=="Info") read_grid_info(sub);
      if(sub->keyword()=="Data") read_grid_data(sub);
    }

    block->clear();
    delete block;

    return 0;
  }



  void DFISE_MESH::read_grid_info(BLOCK * block)
  {
    grid_info.version = (block->get_float_parameter("version", 0));

    std::string  type = (block->get_string_parameter("type", 0));
    if(type=="boundary")
      grid_info.type=INFO::boundary;
    if(type=="grid")
      grid_info.type=INFO::grid;

    grid_info.dimension    = (block->get_int_parameter("dimension", 0));
    grid_info.nb_vertices  = (block->get_int_parameter("nb_vertices", 0));
    grid_info.nb_edges     = (block->get_int_parameter("nb_edges", 0));
    grid_info.nb_faces     = (block->get_int_parameter("nb_faces", 0));
    grid_info.nb_elements  = (block->get_int_parameter("nb_elements", 0));
    grid_info.nb_regions   = (block->get_int_parameter("nb_regions", 0));

    unsigned int regions   = block->n_values_in_parameter("regions");
    unsigned int materials = block->n_values_in_parameter("materials");
    assert(regions == grid_info.nb_regions);
    assert(materials == grid_info.nb_regions);

    // save region name and materials here
    for(unsigned int n=0; n<grid_info.nb_regions; ++n)
    {
      // replace any black char ' ' in region name with '_'
      std::string region_name = fix_region_name(block->get_string_parameter("regions", n));
      grid_info.regions.push_back(region_name);
      grid_info.materials.push_back(block->get_string_parameter("materials", n));
    }

    // set field region and boundary region information
    int n_field_region=0;
    int n_boundary_region=0;
    for(unsigned int n=0; n<grid_info.nb_regions; ++n)
    {
      if(grid_info.is_boundary_region(n))
      {
        grid_info.set_boundary_region(n, n_boundary_region++);
      }
      else if(grid_info.is_interface_region(n))
      {
        // do nothing
      }
      else
      {
        grid_info.set_field_region(n, n_field_region++);
      }
    }

    //grid_info.print(std::cout);
  }



  void DFISE_MESH::read_grid_data(BLOCK *block)
  {
    grid.dimension = grid_info.dimension;

    // set CoordSystem
    {
      BLOCK * CoordSystem = block->get_sub_block("CoordSystem");
      assert(CoordSystem!=NULL);

      for(unsigned int n=0; n<3; ++n)
        grid.translate[n] = (CoordSystem->get_float_parameter("translate", n));

      for(unsigned int m=0; m<3; ++m)
        for(unsigned int n=0; n<3; ++n)
          grid.transform[m*3+n] = (CoordSystem->get_float_parameter("transform", m*3+n));
    }

    // read Vertices
    {
      BLOCK * Vertices = block->get_sub_block("Vertices");
      assert(Vertices!=NULL);
      assert(Vertices->index()==static_cast<int>(grid_info.nb_vertices));

      grid.Vertices.reserve(grid_info.nb_vertices);
      if(grid_info.dimension == 3)
        for(unsigned int n=0; n<grid_info.nb_vertices; ++n)
        {
          Point p;
          p[0] = Vertices->get_float_value(3*n);
          p[1] = Vertices->get_float_value(3*n+1);
          p[2] = Vertices->get_float_value(3*n+2);
          grid.Vertices.push_back(p);
        }
      if(grid_info.dimension == 2)
        for(unsigned int n=0; n<grid_info.nb_vertices; ++n)
        {
          Point p;
          p[0] = Vertices->get_float_value(2*n);
          p[1] = Vertices->get_float_value(2*n+1);
          p[2] = 0;
          grid.Vertices.push_back(p);
        }
    }

    // read edges
    {
      BLOCK * Edges = block->get_sub_block("Edges");
      assert(Edges!=NULL);
      assert(Edges->index()==static_cast<int>(grid_info.nb_edges));

      grid.Edges.reserve(grid_info.nb_edges);
      for(unsigned int n=0; n<grid_info.nb_edges; ++n)
      {
        grid.add_edge(std::make_pair((Edges->get_int_value(2*n)), (Edges->get_int_value(2*n+1))));
      }
    }

    // read faces
    if(grid_info.dimension == 3)
    {
      BLOCK * Faces = block->get_sub_block("Faces");
      assert(Faces!=NULL);
      assert(Faces->index()==static_cast<int>(grid_info.nb_faces));
      grid.Faces.resize(grid_info.nb_faces);

      unsigned int next=0;
      for(unsigned int n=0; n<grid_info.nb_faces; ++n)
      {
        int nb_edges = Faces->get_int_value(next++);
        for(int e=0; e<nb_edges; e++)
          grid.Faces[n].push_back((Faces->get_int_value(next++)));
      }
    }

    // read locations
    {
      BLOCK * Locations = block->get_sub_block("Locations");
      assert(Locations!=NULL);

      if(grid_info.dimension == 2)
        assert(Locations->index()==static_cast<int>(grid_info.nb_edges));

      if(grid_info.dimension == 3)
        assert(Locations->index()==static_cast<int>(grid_info.nb_faces));

      for(unsigned int n=0; n<Locations->n_values(); ++n)
      {
        std::string location_string = Locations->get_string_value(n);
        for(unsigned int c=0; c<location_string.size(); ++c )
          grid.Locations.push_back(location_string[c]);
      }

      assert( grid.Locations.size() == Locations->index() );

    }

    // read Elements
    {
      BLOCK * Elements = block->get_sub_block("Elements");
      assert(Elements!=NULL);
      assert(Elements->index()==static_cast<int>(grid_info.nb_elements));
      grid.Elements.resize(grid_info.nb_elements);

      unsigned int next=0;
      for(unsigned int n=0; n<grid_info.nb_elements; ++n)
      {
        int elem_code = Elements->get_int_value(next++);
        grid.Elements[n].elem_code = elem_code;
        switch(elem_code)
        {
        case 1: //Segment
          {
            for(unsigned int i=0; i<2; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            break;
          }
        case 2: //Triangle
          {
            for(unsigned int i=0; i<3; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            break;
          }
        case 3: //Rectangle
          {
            for(unsigned int i=0; i<4; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            break;
          }
        case 5: //Tetrahedron
          {
            for(unsigned int i=0; i<4; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            break;
          }
        case 6: //Pyramid
          {
            for(unsigned int i=0; i<5; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            break;
          }
        case 7: //Prism
          {
            for(unsigned int i=0; i<5; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            break;
          }
        case 8: //Brick
          {
            for(unsigned int i=0; i<6; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            break;
          }
        default :
          {
            //
            std::cout<<"  ERROR: DFISE with unsupported element."<< std::endl;
            exit(0);
          }
        }
        grid.build_node(grid.Elements[n]);
      }
    }

    // read region(material) information
    //NOTE! when material equals to "Interface" or "Contact", we need some special process
    int region_index=0;
    grid.region_elements.resize(grid_info.nb_regions);
    for(unsigned int n=0; n<block->n_sub_blocks(); ++n)
    {
      BLOCK * region = block->get_sub_block(n);

      if(region->keyword()=="Region")
      {
        region->set_label( fix_region_name(region->label()) );
        assert(region->label()==grid_info.regions[region_index]);
        std::string material = region->get_string_parameter("material", 0);
        assert(material==grid_info.materials[region_index]);

        grid.regions.push_back(region->label());
        grid.materials.push_back(material);

        BLOCK * Elements = region->get_sub_block("Elements");
        unsigned int n_elem = Elements->index();
        for(unsigned int i=0; i<n_elem; ++i)
        {
          unsigned int elem_index = (Elements->get_int_value(i));
          grid.Elements[elem_index].region_index = region_index;
          grid.region_elements[region_index].push_back(elem_index);
        }
        region_index++;
      }
    }


  }



  int DFISE_MESH::parse_dfise_dataset_file(const std::string & dataset_file)
  {
    std::cout<<"  Reading DF-ISE dataset file " << dataset_file << "..."<< std::endl;

    BLOCK *block = new BLOCK;

    yyin = fopen(dataset_file.c_str(), "r");
    assert( yyin != NULL );
    assert(!yyparse(block));
    fclose(yyin);
    YY_FLUSH_BUFFER;

    for(unsigned int n=0; n<block->n_sub_blocks(); ++n)
    {
      BLOCK * sub = block->get_sub_block(n);
      if(sub->keyword()=="Info") read_dataset_info(sub);
      if(sub->keyword()=="Data") read_dataset_data(sub);
    }

    block->clear();
    delete block;

    return 0;
  }


  void DFISE_MESH::read_dataset_info(BLOCK * block)
  {
    data_info.version = (block->get_float_parameter("version", 0));

    std::string  type = (block->get_string_parameter("type", 0));
    assert(type=="dataset");
    data_info.type=INFO::dataset;

    data_info.dimension    = (block->get_int_parameter("dimension", 0));
    data_info.nb_vertices  = (block->get_int_parameter("nb_vertices", 0));
    data_info.nb_edges     = (block->get_int_parameter("nb_edges", 0));
    data_info.nb_faces     = (block->get_int_parameter("nb_faces", 0));
    data_info.nb_elements  = (block->get_int_parameter("nb_elements", 0));
    data_info.nb_regions   = (block->get_int_parameter("nb_regions", 0));

    assert(data_info.dimension   == grid_info.dimension);
    assert(data_info.nb_vertices == grid_info.nb_vertices);
    assert(data_info.nb_edges    == grid_info.nb_edges);
    assert(data_info.nb_faces    == grid_info.nb_faces);
    assert(data_info.nb_elements == grid_info.nb_elements);
    assert(data_info.nb_regions  == grid_info.nb_regions);

    unsigned int datasets  = block->n_values_in_parameter("datasets");
    unsigned int functions = block->n_values_in_parameter("functions");
    assert(datasets == functions);
    for(unsigned int n=0; n<datasets; ++n)
    {
      data_info.datasets.push_back((block->get_string_parameter("datasets", n)));
      data_info.functions.push_back((block->get_string_parameter("functions", n)));
    }

    //data_info.print(std::cout);
  }


  void DFISE_MESH::read_dataset_data(BLOCK *block)
  {
    for(unsigned int n=0; n<block->n_sub_blocks(); ++n)
    {
      BLOCK * dataset_block = block->get_sub_block(n);
      assert(dataset_block->keyword()=="Dataset");

      DATASET * dataset = new DATASET;

      dataset->name = dataset_block->label();
      dataset->function = dataset_block->get_string_parameter("function",0);
      dataset->type = dataset_block->get_string_parameter("type",0)=="scalar" ? DATASET::scalar : DATASET::vector;
      dataset->dimension= dataset_block->get_int_parameter("dimension",0);
      dataset->location=DATASET::location_string_to_enum(dataset_block->get_string_parameter("location",0));

      int n_validity = dataset_block->n_values_in_parameter("validity");
      for(int i=0; i<n_validity; ++i)
      {
        std::string region = fix_region_name(dataset_block->get_string_parameter("validity", i));
        int region_index = grid_info.fieldregion_index_by_label(region);
        assert(region_index>=0);

        dataset->validity.push_back(region);
        dataset->Regions.push_back(static_cast<unsigned int>(region_index));
      }

      //read data
      if(dataset->type==DATASET::scalar)
      {
        BLOCK * Values = dataset_block->get_sub_block("Values");
        dataset->n_data = Values->index();
        assert(Values->n_values()==dataset->n_data);
        for(unsigned int i=0; i<dataset->n_data; ++i)
          dataset->Scalar_Values.push_back(Values->get_float_value(i));
      }

      //
      if(dataset->type==DATASET::vector)
      {
        int index=0;

        BLOCK * Values = dataset_block->get_sub_block("Values");
        dataset->n_data = Values->index();
        assert(Values->n_values()==dataset->n_data);

        int n_vector_data = dataset->n_data/dataset->dimension;
        for(int i=0; i<n_vector_data; ++i)
        {
          std::vector<double> vector_value;
          for(unsigned int n=0; n<dataset->dimension; ++n)
            vector_value.push_back(Values->get_float_value(index++));
          dataset->Vector_Values.push_back(vector_value);
        }
      }

      // build dataset value -> grid vertex map
      //dataset->node_to_value_index_vector.resize(dataset->n_data);
      std::set<unsigned int> node_set;
      for(unsigned int r=0; r<grid_info.nb_regions; ++r)
      {
        if(!dataset->is_valid(grid_info.region_label(r))) continue;

        for(unsigned int e=0; e<grid.region_elements[r].size(); ++e)
        {
          unsigned int elem_id = grid.region_elements[r][e];
          const  Element & elem = grid.Elements[elem_id];

          for(unsigned int m=0; m<elem.vertices.size(); ++m)
          {
            unsigned int id = elem.vertices[m];
            node_set.insert(id);
          }
        }
      }

      for(std::set<unsigned int>::iterator it=node_set.begin(); it!=node_set.end(); ++it)
          dataset->node_to_value_index_map.insert(std::make_pair(*it, dataset->node_to_value_index_map.size()));

      assert(dataset->node_to_value_index_map.size() == dataset->n_data);

      data_sets.push_back(dataset);
    }
  }


  void DFISE_MESH::add_dataset(DATASET * dataset, bool copy)
  {
    if(copy)
    {
      DATASET * new_dataset = new DATASET(*dataset);
      data_sets.push_back(new_dataset);
    }
    else
    {
      data_sets.push_back(dataset);
    }
  }


  std::string DFISE_MESH::fix_region_name(const std::string &name)
  {
    std::string fixed_name = name;
    // remove any blank char at the begin/end of string
    if( fixed_name.find_first_not_of(' ') != 0 )
      fixed_name  = fixed_name.substr(fixed_name.find_first_not_of(' '), fixed_name.size()-1);

    if( fixed_name.rfind(' ') != std::string::npos )
      fixed_name  = fixed_name.substr(0, fixed_name.find_last_not_of(' ')+1);

    // replace blank char with '_'
    std::string filt_elems(" \t");
    std::string::size_type pos = 0;
    while (( pos = fixed_name.find_first_of( filt_elems, pos )) != std::string::npos )
    {
      fixed_name.replace(pos, 1, "_");
    }
    return fixed_name;
  }



  void DFISE_MESH::export_vtk(const std::string & file)
  {
    std::ofstream out(file.c_str(), std::ofstream::trunc);

    out << "# vtk DataFile Version 3.0"      <<'\n';
    out << "Date Generated by Genius TCAD"   <<'\n';
    out << "ASCII"                           <<'\n';
    out << "DATASET UNSTRUCTURED_GRID"       <<'\n';
    out << "POINTS " << grid_info.nb_vertices       << " float" << '\n';

    // write out nodal data
    for(unsigned int n=0; n<grid.Vertices.size(); ++n)
    {
      out << grid.Vertices[n][0] << " " << grid.Vertices[n][1] << " " << grid.Vertices[n][2] << '\n';
    }

    out << std::endl;

    int cell_size = grid.Elements.size();

    for(unsigned int n=0; n<grid.Elements.size(); ++n)
    {
      cell_size   += grid.Elements[n].vertices.size();
    }

    out<<"CELLS "<<grid.Elements.size()<<" "<<cell_size<<'\n';

    for(unsigned int n=0; n<grid.Elements.size(); ++n)
    {
      const Element & elem = grid.Elements[n];
      out << elem.vertices.size()<<" ";
      switch(elem.elem_code)
      {
      case 1:
        out << elem.vertices[0]<<" "<< elem.vertices[1];  //VTK_LINE;
        break;
      case 2:
        out << elem.vertices[0]<<" "<< elem.vertices[1]<<" "<< elem.vertices[2];  //VTK_TRIANGLE;
        break;// 3
      case 3:
        out << elem.vertices[0]<<" "<< elem.vertices[1]<<" "<< elem.vertices[2]<<" "<< elem.vertices[3];  //VTK_QUAD;
        break;// 5
      case 5:
        out << elem.vertices[2]<<" "<< elem.vertices[1]<<" "<< elem.vertices[0]<<" "<< elem.vertices[3]; //VTK_TETRA;
        break;// 8
      case 6:
        out << elem.vertices[3]<<" "<< elem.vertices[2]<<" "<< elem.vertices[1]<<" "<< elem.vertices[0]<<" "<< elem.vertices[4]; //VTK_PYRAMID;
        break;// 16
      case 7:
        out << elem.vertices[2]<<" "<< elem.vertices[1]<<" "<< elem.vertices[0]<<" "
        << elem.vertices[3]<<" "<< elem.vertices[4]<<" "<< elem.vertices[5]; //VTK_WEDGE;
        break;// 13
      case 8:
        out << elem.vertices[3]<<" "<< elem.vertices[2]<<" "<< elem.vertices[1]<<" "<< elem.vertices[0]<<" "
            << elem.vertices[4]<<" "<< elem.vertices[5]<<" "<< elem.vertices[6]<<" "<< elem.vertices[7]; //VTK_HEXAHEDRON;
        break;// 10
      default:
        {
          std::cerr<<"element type "<<grid.Elements[n].elem_code<<" not implemented"<<std::endl;
        }
      }
      out << std::endl;
    }
    out << std::endl;

    out << "CELL_TYPES " << grid.Elements.size() << '\n';
    for(unsigned int n=0; n<grid.Elements.size(); ++n)
    {
      unsigned int celltype;
      switch(grid.Elements[n].elem_code)
      {
      case 1:
        celltype = 3;  //VTK_LINE;
        break;
      case 2:
        celltype = 5;  //VTK_TRIANGLE;
        break;// 3
      case 3:
        celltype = 9;  //VTK_QUAD;
        break;// 5
      case 5:
        celltype = 10; //VTK_TETRA;
        break;// 8
      case 6:
        celltype = 14; //VTK_PYRAMID;
        break;// 16
      case 7:
        celltype = 13; //VTK_WEDGE;
        break;// 13
      case 8:
        celltype = 12; //VTK_HEXAHEDRON;
        break;// 10
      default:
        {
          std::cerr<<"element type "<<grid.Elements[n].elem_code<<" not implemented"<<std::endl;
        }
      }
      out << celltype << std::endl;
    }

    out << std::endl;

    // region information
    out << "CELL_DATA "<<grid.Elements.size()     <<'\n';
    out << "SCALARS region float 1"  <<'\n';
    out << "LOOKUP_TABLE default"    <<'\n';
    for(unsigned int n=0; n<grid.Elements.size(); ++n)
    {
      out << static_cast<float>(grid.Elements[n].region_index)<<'\n';
    }
    out<<std::endl;

    /*
        //write node based boundary info to vtk
        out<<"POINT_DATA "<<grid.Vertices.size()      <<'\n';
        out<<std::endl;

        out<<"SCALARS node_value float 1" <<'\n';
        out<<"LOOKUP_TABLE default"       <<'\n';
        for(unsigned int n=0; n<grid.Vertices.size(); ++n)
        {
          double value = get_scaler_value("DopingConcentration", "body", n);
          out<<value<<'\n' ;
        }
        out<<std::endl;
    */
    out.close();
  }

}
