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


namespace DFISE
{

  /**
   * data block of dataset file
   */
  class DATASET
  {
  public:

    DATASET() {}

    ~DATASET() { this->clear(); }

    void clear()
    {
      dimension = 0;
      location  = invalid_location;
      n_data = 0;
      validity.clear();
      Regions.clear();
      node_to_value_index_map.clear();
      Scalar_Values.clear();
      Vector_Values.clear();
    }

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

    double get_scaler_value_by_node_index(unsigned int node_index) const
    {
      assert(node_to_value_index_map.find(node_index)!=node_to_value_index_map.end());
      //return Scalar_Values[node_to_value_index_vector[node_id]];
      return Scalar_Values[node_to_value_index_map.find(node_index)->second];
    }

    /**
     * output
     */
    void print( std::ostream & out ) const
    {
      out << "  Dataset (" << '"' << name << '"' << ") {" << std::endl;
      out << "    function  = " << function << std::endl;
      if(type == scalar)
        out << "    type      = scalar" << std::endl;
      if(type == vector)
        out << "    type      = vector" << std::endl;
      out<<  "    dimension = " << dimension << std::endl;
      switch(location)
      {
        case vertex  : out<<  "    location  = vertex" << std::endl;  break;
        case edge    : out<<  "    location  = edge" << std::endl;break;
        case face    : out<<  "    location  = face" << std::endl;break;
        case element : out<<  "    location  = element" << std::endl;break;
        case region  : out<<  "    location  = region" << std::endl;break;
      }
      out<<  "    validity  = [";
      for(unsigned int n=0; n<validity.size(); ++n )
        out << " \"" << validity[n]<<"\" ";
      out << "]" << std::endl;

      out<<  "    Values (" << n_data << ") {" << std::endl;
      out<< std::scientific << std::setprecision(15);
      if(type == scalar)
      {
        for(unsigned int n=0; n<Scalar_Values.size(); ++n)
          out <<  "      "  << Scalar_Values[n] << std::endl;
      }
      if(type == vector)
      {
        for(unsigned int n=0; n<Vector_Values.size(); ++n)
        {
          out <<  "      ";
          for(unsigned int d=0; d<Vector_Values[n].size(); ++d)
            out << Vector_Values[n][d] << " ";
          out<<std::endl;
        }
      }

      out<<"    }"<< std::endl;

      out<<"  }"<< std::endl;
    }

  public:

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
     * map the node id to value.
     * please note: the order of dataset is sort as node index
     * as a result, although dataset may have multi valid regions,
     * the point on interface of adjacent regions will only count once
     */
    std::map<unsigned int, unsigned int> node_to_value_index_map;

    /**
     * scaler value
     */
    std::vector<double> Scalar_Values;

    /**
     * vector value
     */
    std::vector< std::vector<double> > Vector_Values;
  };

}
