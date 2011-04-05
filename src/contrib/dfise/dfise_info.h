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

#include <cassert>
#include <vector>
#include <string>
#include <map>



namespace DFISE
{

  /**
   * info block
   */
  class INFO
  {
  public:

    INFO():version(1.0), type(invalid_type)
    {
      dimension = 0;
      nb_vertices = 0;
      nb_edges = 0;
      nb_faces = 0;
      nb_elements = 0;
      nb_regions = 0;
    }

    ~INFO()
    { this->clear(); }


    /**
     * different types of DF-ISE file in INFO block
     */
    enum TYPE{property, layout, cell, recursive_tensor, variable_tensor, boundary, grid, dataset, invalid_type};

    /**
     * clear INFO
     */
    void clear()
    {
      version = 1.0;
      type = invalid_type;
      dimension = 0;
      nb_vertices = 0;
      nb_edges = 0;
      nb_faces = 0;
      nb_elements = 0;
      nb_regions = 0;
      regions.clear();
      materials.clear();
      _region_to_fieldregion_map.clear();
      _fieldregion_to_region_map.clear();
      _region_to_boundaryregion_map.clear();
      _boundaryregion_to_region_map.clear();
      datasets.clear();
      functions.clear();
    }

    /**
     * @return the region label of ith region
     */
    const std::string & region_label(const int region_index) const
    {
      assert(region_index<static_cast<int>(regions.size()));
      return regions[region_index];
    }

    /**
     * set field region to region map
     */
    void set_field_region(int region_index, int fieldregion_index)
    {
      _region_to_fieldregion_map[region_index] = fieldregion_index;
      _fieldregion_to_region_map[fieldregion_index] = region_index;
    }

    /**
     * set boundary region to region map
     */
    void set_boundary_region(int region_index, int boundaryregion_index)
    {
      _region_to_boundaryregion_map[region_index] = boundaryregion_index;
      _boundaryregion_to_region_map[boundaryregion_index] = region_index;
    }


    /**
     * @return field region index by region index
     */
    int region_to_fieldregion(int region_index) const
    {
      if(_region_to_fieldregion_map.find(region_index) != _region_to_fieldregion_map.end())
        return _region_to_fieldregion_map.find(region_index)->second;
      return -1;
    }

    /**
     * @return region_label by fieldregion index
     */
    const std::string & fieldregion_label(int fieldregion_index) const
    {
      return region_label(_fieldregion_to_region_map.find(fieldregion_index)->second);
    }

    /**
     * @return fieldregion index by region label
     */
    int fieldregion_index_by_label(const std::string & label) const
    {
      for(unsigned int n=0; n<regions.size(); ++n)
        if(regions[n]==label)
          return _region_to_fieldregion_map.find(n)->second;
      return -1;
    }

    /**
     * @return the number of field regions
     */
    unsigned int n_field_regions() const
    { return _region_to_fieldregion_map.size(); }


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
      { return _region_to_boundaryregion_map.size(); }

    /**
     * @return boundary region index by region index
     */
    int region_to_boundaryregion(int region_index) const
    {
      if(_region_to_boundaryregion_map.find(region_index) != _region_to_boundaryregion_map.end())
        return _region_to_boundaryregion_map.find(region_index)->second;
      return -1;
    }

    /**
     * @return boundary_label by boundary index
     */
    const std::string & boundary_label(int boundary_index) const
    {
      return region_label(_boundaryregion_to_region_map.find(boundary_index)->second);
    }

    /**
     * output
     */
    void print( std::ostream & out ) const
    {
      out<<"Info  {"<< std::endl;
      out<<"  version = "    << std::fixed << std::setprecision(1) << version     <<std::endl;

      switch(type)
      {
          case boundary : out<<"  type = boundary"<<std::endl;  break;
          case grid     : out<<"  type = grid"    <<std::endl;  break;
          case dataset  : out<<"  type = dataset" <<std::endl;  break;
          default : break;
      }

      out<<"  dimension   = " << dimension   <<std::endl;
      out<<"  nb_vertices = " << nb_vertices <<std::endl;
      out<<"  nb_edges    = " << nb_edges    <<std::endl;
      out<<"  nb_faces    = " << nb_faces    <<std::endl;
      out<<"  nb_elements = " << nb_elements <<std::endl;
      out<<"  nb_regions  = " << nb_regions  <<std::endl;

      if(regions.size())
      {
        out<<"  regions  = [";
        for(unsigned int n=0; n<regions.size(); ++n)
          out<<" \""<<regions[n]<<"\" ";
        out<< "]" << std::endl;

        out<<"  materials  = [";
        for(unsigned int n=0; n<materials.size(); ++n)
          out<<" "<<materials[n]<<" ";
        out<< "]" << std::endl;
      }

      if(datasets.size())
      {
        out<<"  datasets  = [";
        for(unsigned int n=0; n<datasets.size(); ++n)
          out<<" \""<<datasets[n]<<"\" ";
        out<< "]" << std::endl;

        out<<"  functions  = [";
        for(unsigned int n=0; n<functions.size(); ++n)
          out<<" "<<functions[n]<<" ";
        out<< "]" << std::endl;
      }

      out<<"}"<< std::endl;
      out<<std::endl;
    }


  public:

    /**
     * version of DFISE file, it seems should be at present 1.0
     */
    double version;

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

    /**
     * material info of the region
     */
    std::vector<std::string> materials;

    /**
     * label of all the data set
     */
    std::vector<std::string> datasets;

    /**
     * label of all the functions
     */
    std::vector<std::string> functions;

  private:

    /**
     * region to field region map
     */
    std::map<int, int> _region_to_fieldregion_map;

    /**
       * field region to region map
     */
    std::map<int, int> _fieldregion_to_region_map;

    /**
       * region to boundary map
     */
    std::map<int, int> _region_to_boundaryregion_map;

    /**
       * boundary to region map
     */
    std::map<int, int> _boundaryregion_to_region_map;

  };

}

