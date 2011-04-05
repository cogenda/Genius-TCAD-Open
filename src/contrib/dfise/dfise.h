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

#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>


#include "dfise_info.h"
#include "dfise_grid.h"
#include "dfise_dataset.h"

namespace DFISE
{

  class BLOCK;

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

    /**
     * write dfise file
     */
    void write_dfise(const std::string & file) const;



    /**
     * @return const reference to grid info
     */
    const INFO & get_grid_info() const
      { return grid_info;}

    /**
     * @return reference to grid info
     */
    INFO & get_grid_info()
    { return grid_info;}

    /**
     * @return const reference to grid data
     */
    const GRID & get_grid() const
      { return grid; }

    /**
     * @return reference to grid data
     */
    GRID & get_grid()
    { return grid; }

    /**
     * @return const reference to dataset info
     */
    const INFO & get_dataset_info() const
    { return data_info;}

    /**
     * @return reference to dataset info
     */
    INFO & get_dataset_info()
    { return data_info;}



    /**
     * @return total numbers of dataset
     */
    unsigned int n_datasets() const
      { return data_sets.size(); }


    /**
     * add a new data set
     */
    void add_dataset(DATASET * dataset, bool copy=false);

    /**
     * get ith dataset
     */
    DATASET * get_dataset(unsigned int n)
    { return data_sets[n]; }


    /**
     * @return true when value with name value_name exist
     */
    bool is_value_exist(const std::string & value_name) const
    {
      for(unsigned int n=0; n<data_sets.size(); ++n)
        if(data_sets[n]->name == value_name)
          return true;
      return false;
    }

    /**
     * @return true when value with fuzzy name value_name exist
     */
    bool is_value_exist_fuzzy(const std::string & value_name) const
    {
      for(unsigned int n=0; n<data_sets.size(); ++n)
        if(data_sets[n]->name.find(value_name)!=std::string::npos)
          return true;
      return false;
    }

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
     * @return true when value with fuzzy name value_name exist and valid in region
     */
    bool is_value_exist_fuzzy(const std::string & value_name, unsigned int region) const
    {
      for(unsigned int n=0; n<data_sets.size(); ++n)
        if(data_sets[n]->name.find(value_name)!=std::string::npos && data_sets[n]->is_valid(region))
          return true;
      return false;
    }

    /**
     * @return the node scalar value with value_name, and the value should be valid in region
     */
    double get_scaler_value(const std::string & value_name, unsigned int region, unsigned int node_index) const
    {
      double value = 0.0;
      for(unsigned int n=0; n<data_sets.size(); ++n)
        if(data_sets[n]->name == value_name && data_sets[n]->is_valid(region))
        {
          value += data_sets[n]->get_scaler_value_by_node_index(node_index);
        }
      return value;
    }


    /**
     * @return the node scalar value with value_name, and the value should be valid in region
     */
    double get_scaler_value(int data_set_index, unsigned int region, unsigned int node_index) const
    {
      double value =0;
      if(data_sets[data_set_index]->is_valid(region))
      {
        value = data_sets[data_set_index]->get_scaler_value_by_node_index(node_index);
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
     * DFISE region may contain black char ' ', replace it with '_'
     */
    std::string fix_region_name(const std::string &s);

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
