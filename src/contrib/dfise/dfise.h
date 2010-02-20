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



#ifndef __dfise_h__
#define __dfise_h__


#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <cstring>

#include "enum_elem_type.h" //for libmesh compatible element type
#include "vector_value.h"
#include "tensor_value.h"


namespace DFISE
{

  class BLOCK;

  /**
   * info block
   */
  struct INFO
  {
    /**
     * version of DFISE file, it seems should be at present 1.0
     */
    double version;

    /**
    * different types of DF-ISE file in INFO block
    */
    enum TYPE{property, layout, cell, recursive_tensor, variable_tensor, boundary, grid, dataset};

    /**
     * the type of INFO block, different DFISE file has different TYPE
     */
    TYPE type;

    unsigned int dimension;  // 1 | 2 | 3

    unsigned int nb_vertices;

    unsigned int nb_edges;

    unsigned int nb_faces;

    unsigned int nb_elements;

    unsigned int nb_regions;

    /**
     * regions in df-ise file.
     * however, it contains field region (dim element) and boundary region (dim-1 element)
     * the "material" for boundary region is "Interface" and "Contact"
     */
    std::vector<std::string> regions;

    const std::string & region_label(const int region_index) const
    {
      assert(region_index<static_cast<int>(regions.size()));
      return regions[region_index];
    }

    /**
     * material info of the region
     */
    std::vector<std::string> materials;


    /**
     * region to field region map
     */
    std::map<int, int> region_to_fieldregion_map;

    /**
     * field region to region map
     */
    std::map<int, int> fieldregion_to_region_map;

    /**
     * @return the number of field regions
     */
    unsigned int n_field_regions() const
      { return region_to_fieldregion_map.size(); }

    /**
     * @return field region index by region index
     */
    int region_to_fieldregion(int region_index) const
    {
      if(region_to_fieldregion_map.find(region_index)!=region_to_fieldregion_map.end())
        return region_to_fieldregion_map.find(region_index)->second;
      return -1;
    }

    /**
     * @return region_label by fieldregion index
     */
    const std::string & fieldregion_label(int fieldregion_index) const
    {
      return region_label(fieldregion_to_region_map.find(fieldregion_index)->second);
    }

    int fieldregion_index_by_label(const std::string & label) const
    {
      for(unsigned int n=0; n<regions.size(); ++n)
        if(regions[n]==label)
          return region_to_fieldregion_map.find(n)->second;
      return -1;
    }

    /**
     * region to boundary map
     */
    std::map<int, int> region_to_boundaryregion_map;
    /**
     * boundary to region map
     */
    std::map<int, int> boundaryregion_to_region_map;

    /**
     * @return true if region is only a boundary region (all the element has dimension-1)
     */
    bool is_boundary_region(int region_index) const
    { return (materials[region_index]=="Contact"); }

    /**
     * @return true if region is only a interface region (all the element has dimension-1)
     */
    bool is_interface_region(int region_index) const
      { return (materials[region_index]=="Interface"); }

    /**
     * @return the number of boundary regions
     */
    unsigned int n_boundary_regions() const
      { return region_to_boundaryregion_map.size(); }

    /**
     * @return boundary region index by region index
     */
    int region_to_boundaryregion(int region_index) const
    {
      if(region_to_boundaryregion_map.find(region_index)!=region_to_boundaryregion_map.end())
        return region_to_boundaryregion_map.find(region_index)->second;
      return -1;
    }

    /**
     * @return boundary_label by boundary index
     */
    const std::string & boundary_label(int boundary_index) const
    {
      return region_label(boundaryregion_to_region_map.find(boundary_index)->second);
    }

    std::vector<std::string> datasets;

    std::vector<std::string> functions;


    /**
     * debug output
     */
    void print()
    {
      std::cout<<"DF-ISE Info"  << std::endl;
      std::cout<<" version "    << version     <<std::endl;

      switch(type)
      {
      case boundary : std::cout<<" type boundary"<<std::endl;  break;
      case grid     : std::cout<<" type grid"    <<std::endl;  break;
      case dataset  : std::cout<<" type dataset" <<std::endl;  break;
      default : break;
      }

      std::cout<<" dimension "  << dimension   <<std::endl;
      std::cout<<" nb_vertices "<< nb_vertices <<std::endl;
      std::cout<<" nb_edges "   << nb_edges    <<std::endl;
      std::cout<<" nb_faces "   << nb_faces    <<std::endl;
      std::cout<<" nb_elements "<< nb_elements <<std::endl;
      std::cout<<" nb_regions " << nb_regions  <<std::endl;

      if(regions.size())
      {
        for(unsigned int n=0; n<regions.size(); ++n)
          std::cout<<" region: "<<regions[n]<<", material: "<<materials[n]<<std::endl;
      }

      if(datasets.size())
      {
        for(unsigned int n=0; n<datasets.size(); ++n)
          std::cout<<" dataset: "<<datasets[n]<<", function: "<<functions[n]<<std::endl;
      }
    }
  };


  struct Element
  {
    ElemType elem_type;
    std::vector<int> faces;
    std::vector<int> vertices;
    int region_index;
    int bc_index;
  };

  /**
   * data block of grid/boundary file
   */
  struct GRID
  {
    VectorValue<double> translate;

    TensorValue<double> transform;

    std::vector< VectorValue<double> > Vertices;

    std::vector< std::pair<int, int> > Edges;

    std::vector< std::vector<int> >    Faces;

    /**
     * @return the nodes belongs to a face in special order
     */
    std::vector<int> get_face_nodes(int face_index)
    {
      bool face_inverse = face_index < 0;
      face_index = face_index < 0 ? -face_index-1 : face_index;
      assert(static_cast<unsigned int>(face_index) < Faces.size());

      std::vector<int> nodes;

      const std::vector<int> & face = Faces[face_index];
      for(unsigned int n=0; n<face.size(); ++n)
      {
        int edge_index = face[n];
        bool edge_inverse = edge_index < 0;
        edge_index = edge_index < 0 ? -edge_index-1 : edge_index;
        if(!edge_inverse)
        {
          nodes.push_back(Edges[edge_index].first);
          nodes.push_back(Edges[edge_index].second);
        }
        else
        {
          nodes.push_back(Edges[edge_index].second);
          nodes.push_back(Edges[edge_index].first);
        }
      }
      std::unique(nodes.begin(), nodes.end());
      nodes.resize(face.size());

      if(face_inverse) std::reverse(nodes.begin(), nodes.end());

      return nodes;
    }


    enum LOCATION{i, f, e, u};

    std::vector< LOCATION >            Locations;

    std::vector< Element >             Elements;

    void build_edge2_node(Element & elem)
    {
      std::vector<int> nodes;

      int edge_index = elem.faces[0];
      bool inverse = edge_index < 0;
      edge_index = edge_index < 0 ? -edge_index-1 : edge_index;
      if(!inverse)
      {
        nodes.push_back(Edges[edge_index].first);
        nodes.push_back(Edges[edge_index].second);
      }
      else
      {
        nodes.push_back(Edges[edge_index].second);
        nodes.push_back(Edges[edge_index].first);
      }
      elem.vertices = nodes;
      assert(elem.vertices.size()==2);
    }

    void build_tri3_node(Element & elem)
    {
      std::vector<int> nodes;
      for(unsigned int n=0; n<elem.faces.size(); ++n)
      {
        int edge_index = elem.faces[n];
        bool inverse = edge_index < 0;
        edge_index = edge_index < 0 ? -edge_index-1 : edge_index;
        if(!inverse)
        {
          nodes.push_back(Edges[edge_index].first);
          nodes.push_back(Edges[edge_index].second);
        }
        else
        {
          nodes.push_back(Edges[edge_index].second);
          nodes.push_back(Edges[edge_index].first);
        }
      }
      std::unique(nodes.begin(), nodes.end());
      nodes.resize(elem.faces.size());

      // check for ccw
      {
        VectorValue<double> & p0 = Vertices[nodes[0]];
        VectorValue<double> & p1 = Vertices[nodes[1]];
        VectorValue<double> & p2 = Vertices[nodes[2]];

        if( ((p1-p0).cross(p2-p0))(2) < 0)
          std::reverse(nodes.begin(), nodes.end());
      }
      elem.vertices = nodes;
      assert(elem.vertices.size()==3);
    }


    void build_quad4_node(Element & elem)
    {
      std::vector<int> nodes;
      for(unsigned int n=0; n<elem.faces.size(); ++n)
      {
        int edge_index = elem.faces[n];
        bool inverse = edge_index < 0;
        edge_index = edge_index < 0 ? -edge_index-1 : edge_index;
        if(!inverse)
        {
          nodes.push_back(Edges[edge_index].first);
          nodes.push_back(Edges[edge_index].second);
        }
        else
        {
          nodes.push_back(Edges[edge_index].second);
          nodes.push_back(Edges[edge_index].first);
        }
      }
      std::unique(nodes.begin(), nodes.end());
      nodes.resize(elem.faces.size());

      // check for ccw
      {
        VectorValue<double> & p0 = Vertices[nodes[0]];
        VectorValue<double> & p1 = Vertices[nodes[1]];
        VectorValue<double> & p2 = Vertices[nodes[2]];

        if( ((p1-p0).cross(p2-p0))(2) < 0)
          std::reverse(nodes.begin(), nodes.end());
      }

      elem.vertices = nodes;
      assert(elem.vertices.size()==4);
    }

    void build_tet4_node(Element & elem)
    {
      std::vector<int> nodes3 = get_face_nodes(elem.faces[3]);
      std::vector<int> nodes1 = get_face_nodes(elem.faces[1]);
      elem.vertices = nodes3;
      for(unsigned int m=0; m<nodes1.size(); ++m)
      {
        int p4 = nodes1[m];
        bool skip=false;
        for(unsigned int n=0; n<nodes3.size(); ++n)
          if(p4==nodes3[n])
          {
            skip=true;
            continue;
          }
        if(!skip) {elem.vertices.push_back(p4); break;}
      }
      assert(elem.vertices.size()==4);
    }

    void build_pyramid5_node(Element & elem)
    {
      std::vector<int> nodes4 = get_face_nodes(elem.faces[4]);
      std::vector<int> nodes1 = get_face_nodes(elem.faces[1]);
      elem.vertices = nodes4;
      for(unsigned int m=0; m<nodes1.size(); ++m)
      {
        int p5 = nodes1[m];
        bool skip=false;
        for(unsigned int n=0; n<nodes4.size(); ++n)
          if(p5==nodes4[n])
          {
            skip=true;
            continue;
          }
        if(!skip) {elem.vertices.push_back(p5); break;}
      }
      assert(elem.vertices.size()==5);
    }

    void build_prism6_node(Element & elem)
    {
      std::vector<int> nodes3 = get_face_nodes(elem.faces[3]);
      std::vector<int> nodes4 = get_face_nodes(elem.faces[4]);
      elem.vertices = nodes3;
      for(unsigned int m=0; m<nodes4.size(); ++m)
        elem.vertices.push_back(nodes4[m]);
      assert(elem.vertices.size()==6);
    }

    void build_hex8_node(Element & elem)
    {
      std::vector<int> nodes4 = get_face_nodes(elem.faces[4]);
      std::vector<int> nodes5 = get_face_nodes(elem.faces[5]);
      elem.vertices = nodes4;
      for(unsigned int m=0; m<nodes5.size(); ++m)
        elem.vertices.push_back(nodes5[m]);
      assert(elem.vertices.size()==8);
    }


    std::vector< std::vector<unsigned int> > region_elements;

  };


  /**
   * data block of dataset file
   */
  struct DATASET
  {
    std::string name;

    std::string function;

    enum DATATYPE{scalar, vector};

    DATATYPE type;

    /**
     * (1 for scalar, > 1 for vectors or arrays)
     */
    unsigned int  dimension;

    enum LOCATION{vertex, edge, face, element, region, invalid_location};

    LOCATION location;

    static LOCATION location_string_to_enum(const std::string & location)
    {
      if(location=="vertex")      return vertex;
      else if(location=="edge")   return edge;
      else if(location=="face")   return face;
      else if(location=="element")return element;
      else if(location=="region") return region;
      return invalid_location;
    }

    unsigned int n_data;

    /**
     * which regions this dataset is valid
     * ["<region0>" "<region1>" ...]
     */
    std::vector<std::string> validity;//

    /**
     * convert validity string to region index
     */
    std::vector<unsigned int> Regions;

    /**
     * @return true if r_index in the Regions record
     */
    bool is_valid(unsigned int r_index) const
    {
      for(unsigned int n=0; n<Regions.size(); ++n)
        if(Regions[n]==r_index) return true;
      return false;
    }

    /**
     * @return true if region in the validity vector
     */
    bool is_valid(const std::string & region) const
    {
      for(unsigned int n=0; n<validity.size(); ++n)
        if(validity[n]==region) return true;
      return false;
    }

    /**
     * map the node id to value
     */
    std::map<unsigned int, unsigned int> node_to_value_index_map;

    /**
     * scaler value
     */
    std::vector<double> Scalar_Values;

    double get_scaler_value_by_node_index(unsigned int node_index)
    {
      assert(node_to_value_index_map.find(node_index)!=node_to_value_index_map.end());
      //return Scalar_Values[node_to_value_index_vector[node_id]];
      return Scalar_Values[node_to_value_index_map.find(node_index)->second];
    }

    /**
     * vector value
     */
    std::vector< std::vector<double> > Vector_Values;
  };



  /**
   * a dfise mesh is consisted by grid/boundary file and data file
   */
  class DFISE_MESH
  {
  public:

    DFISE_MESH() {}

    /**
     * free the datasets
     */
    ~DFISE_MESH()
    {
      for(unsigned int n=0; n<data_sets.size(); ++n)
        delete data_sets[n];
    }

    /**
     * read df-ise data into internal data structure
     */
    int parse_dfise(const std::string & file);

    void add_dataset(DATASET * dataset)
  { data_sets.push_back(dataset); }

    const INFO & get_grid_info() const
      { return grid_info;}

    const GRID & get_grid() const
      { return grid; }

    const INFO & get_dataset_info() const
      { return grid_info;}

    unsigned int n_datasets() const
      { return data_sets.size(); }

    DATASET * get_dataset(unsigned int n)
    { return data_sets[n]; }

    /**
     * @return true when value with name value_name exist and valid in region
     */
    bool is_value_exist(const std::string & value_name, unsigned int region) const
    {
      for(unsigned int n=0; n<data_sets.size(); ++n)
        if(data_sets[n]->name == value_name && data_sets[n]->is_valid(region))
          return true;
      return false;
    }

    /**
     * @return the node scalar value with value_name, and the value should be valid in region
     */
    double get_scaler_value(const std::string & value_name, unsigned int region, unsigned int node_index) const
    {
      double value =0;
      for(unsigned int n=0; n<data_sets.size(); ++n)
        if(data_sets[n]->name == value_name && data_sets[n]->is_valid(region))
        {
          value += data_sets[n]->get_scaler_value_by_node_index(node_index);
        }
      return value;
    }

    /**
     * @return the fuzzy matching node scalar value (the name of the value contains value_name),
     * and the value should be valid in region
     */
    double get_scaler_value_by_fuzzy(const std::string & value_name, unsigned int region, unsigned int node_index) const
    {
      double value =0;
      for(unsigned int n=0; n<data_sets.size(); ++n)
        if(data_sets[n]->name.find(value_name)!=std::string::npos && data_sets[n]->is_valid(region))
        {
          value += data_sets[n]->get_scaler_value_by_node_index(node_index);
        }
      return value;
    }

    /**
     * export df-ise mesh in vtk format, debug only
     */
    void export_vtk(const std::string & file);

  private:

    int parse_dfise_grid_file(const std::string & grid_file);

    int parse_dfise_dataset_file(const std::string & dataset_file);

    void read_grid_info(BLOCK *);

    void read_grid_data(BLOCK *);

    void read_dataset_info(BLOCK *);

    void read_dataset_data(BLOCK *);

    /**
     * DFISE file information, should be compatible in grid/dataset file
     */
    INFO grid_info;

    /**
     * DFISE grid information
     */
    GRID grid;

    /**
     * DFISE file information, should be compatible in grid/dataset file
     */
    INFO data_info;

    /**
     * DFISE data information
     */
    std::vector<DATASET *> data_sets;

  };



}








#endif //__dfise_h__
