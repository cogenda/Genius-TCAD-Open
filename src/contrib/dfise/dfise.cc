#include <cstdio>
#include <iostream>
#include <fstream>

#include "config.h"
#include "genius_env.h"

#include "dfise_block.h"
#include "dfise.h"

#ifdef CYGWIN
  #include <io.h>      // for windows _access function
#else
  #include <unistd.h>  // for POSIX access function
#endif

namespace DFISE
{
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

    for(unsigned int n=0; n<grid_info.nb_regions; ++n)
    {
      grid_info.regions.push_back((block->get_string_parameter("regions", n)));
      grid_info.materials.push_back((block->get_string_parameter("materials", n)));
    }

    // set field region and boundary region information
    int n_field_region=0;
    int n_boundary_region=0;
    for(unsigned int n=0; n<grid_info.nb_regions; ++n)
    {
      if(grid_info.is_boundary_region(n))
      {
        grid_info.region_to_boundaryregion_map[n] = n_boundary_region;
        grid_info.boundaryregion_to_region_map[n_boundary_region] = n;
        ++n_boundary_region;
      }
      else if(grid_info.is_interface_region(n))
      {
        // do nothing
      }
      else
      {
        grid_info.region_to_fieldregion_map[n] = n_field_region;
        grid_info.fieldregion_to_region_map[n_field_region] = n;
        ++n_field_region;
      }
    }

    //grid_info.print();
  }



  void DFISE_MESH::read_grid_data(BLOCK *block)
  {
    // set CoordSystem
    {
      BLOCK * CoordSystem = block->get_sub_block("CoordSystem");
      assert(CoordSystem!=NULL);

      for(unsigned int n=0; n<3; ++n)
        grid.translate[n] = (CoordSystem->get_float_parameter("translate", n));

      for(unsigned int m=0; m<3; ++m)
        for(unsigned int n=0; n<3; ++n)
          grid.transform(m,n) = (CoordSystem->get_float_parameter("transform", m*3+n));
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
          VectorValue<double> point((Vertices->get_float_value(3*n)),
                                    (Vertices->get_float_value(3*n+1)),
                                    (Vertices->get_float_value(3*n+2)));
          point = grid.transform*point + grid.translate;
          grid.Vertices.push_back(point);
        }
      if(grid_info.dimension == 2)
        for(unsigned int n=0; n<grid_info.nb_vertices; ++n)
        {
          VectorValue<double> point((Vertices->get_float_value(2*n)),
                                    (Vertices->get_float_value(2*n+1)), 0.0);
          point = grid.transform*point + grid.translate;
          grid.Vertices.push_back(point);
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
        grid.Edges.push_back(std::make_pair((Edges->get_int_value(2*n)), (Edges->get_int_value(2*n+1))));
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
        int elem_id = Elements->get_int_value(next++);
        switch(elem_id)
        {
        case 1: //Segment
          {
            grid.Elements[n].elem_type = EDGE2;
            for(unsigned int i=0; i<2; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            grid.build_edge2_node(grid.Elements[n]);
            break;
          }
        case 2: //Triangle
          {
            grid.Elements[n].elem_type = TRI3;
            for(unsigned int i=0; i<3; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            grid.build_tri3_node(grid.Elements[n]);
            break;
          }
        case 3: //Rectangle
          {
            grid.Elements[n].elem_type = QUAD4;
            for(unsigned int i=0; i<4; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            grid.build_quad4_node(grid.Elements[n]);
            break;
          }
        case 5: //Tetrahedron
          {
            grid.Elements[n].elem_type = TET4;
            for(unsigned int i=0; i<4; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            grid.build_tet4_node(grid.Elements[n]);
            break;
          }
        case 6: //Pyramid
          {
            grid.Elements[n].elem_type = PYRAMID5;
            for(unsigned int i=0; i<5; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            grid.build_pyramid5_node(grid.Elements[n]);
            break;
          }
        case 7: //Prism
          {
            grid.Elements[n].elem_type = PRISM6;
            for(unsigned int i=0; i<5; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            grid.build_prism6_node(grid.Elements[n]);
            break;
          }
        case 8: //Brick
          {
            grid.Elements[n].elem_type = HEX8;
            for(unsigned int i=0; i<6; i++)
              grid.Elements[n].faces.push_back((Elements->get_int_value(next++)));
            grid.build_hex8_node(grid.Elements[n]);
            break;
          }
        default :
          {
            //
            std::cout<<"  ERROR: DFISE with unsupported element."<< std::endl;
            genius_error(); //not supported
          }
        }
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
        assert(region->label()==grid_info.regions[region_index]);
        std::string material = region->get_string_parameter("material", 0);
        assert(material==grid_info.materials[region_index]);

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

    //data_info.print();
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
        std::string region = dataset_block->get_string_parameter("validity", i);
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
      out << grid.Vertices[n](0) << " " << grid.Vertices[n](1) << " " << grid.Vertices[n](2) << '\n';
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
      switch(elem.elem_type)
      {
      case EDGE2:
        out << elem.vertices[0]<<" "<< elem.vertices[1];  //VTK_LINE;
        break;
      case TRI3:
        out << elem.vertices[0]<<" "<< elem.vertices[1]<<" "<< elem.vertices[2];  //VTK_TRIANGLE;
        break;// 3
      case QUAD4:
        out << elem.vertices[0]<<" "<< elem.vertices[1]<<" "<< elem.vertices[2]<<" "<< elem.vertices[3];  //VTK_QUAD;
        break;// 5
      case TET4:
        out << elem.vertices[2]<<" "<< elem.vertices[1]<<" "<< elem.vertices[0]<<" "<< elem.vertices[3]; //VTK_TETRA;
        break;// 8
      case HEX8:
        out << elem.vertices[3]<<" "<< elem.vertices[2]<<" "<< elem.vertices[1]<<" "<< elem.vertices[0]<<" "
        << elem.vertices[4]<<" "<< elem.vertices[5]<<" "<< elem.vertices[6]<<" "<< elem.vertices[7]; //VTK_HEXAHEDRON;
        break;// 10
      case PRISM6:
        out << elem.vertices[2]<<" "<< elem.vertices[1]<<" "<< elem.vertices[0]<<" "
        << elem.vertices[3]<<" "<< elem.vertices[4]<<" "<< elem.vertices[5]; //VTK_WEDGE;
        break;// 13
      case PYRAMID5:
        out << elem.vertices[3]<<" "<< elem.vertices[2]<<" "<< elem.vertices[1]<<" "<< elem.vertices[0]<<" "<< elem.vertices[4]; //VTK_PYRAMID;
        break;// 16
      default:
        {
          std::cerr<<"element type "<<grid.Elements[n].elem_type<<" not implemented"<<std::endl;
        }
      }
      out << std::endl;
    }
    out << std::endl;

    out << "CELL_TYPES " << grid.Elements.size() << '\n';
    for(unsigned int n=0; n<grid.Elements.size(); ++n)
    {
      unsigned int celltype;
      switch(grid.Elements[n].elem_type)
      {
      case EDGE2:
        celltype = 3;  //VTK_LINE;
        break;
      case TRI3:
        celltype = 5;  //VTK_TRIANGLE;
        break;// 3
      case QUAD4:
        celltype = 9;  //VTK_QUAD;
        break;// 5
      case TET4:
        celltype = 10; //VTK_TETRA;
        break;// 8
      case HEX8:
        celltype = 12; //VTK_HEXAHEDRON;
        break;// 10
      case PRISM6:
        celltype = 13; //VTK_WEDGE;
        break;// 13
      case PYRAMID5:
        celltype = 14; //VTK_PYRAMID;
        break;// 16
      default:
        {
          std::cerr<<"element type "<<grid.Elements[n].elem_type<<" not implemented"<<std::endl;
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
