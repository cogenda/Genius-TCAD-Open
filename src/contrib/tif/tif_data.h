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



#ifndef __tif_data_h__
#define __tif_data_h__

#include <list>
#include <vector>
#include <string>
#include <iostream>

namespace TIF
{

  /**
   * Nodal data structure
   */
  struct Node_t
  {
    /**
     * node index, (sequential, starts at 0)
     */
    int    index;

    /**
     * X-coordinate of node
     */
    double x;

    /**
     * Y-coordinate of node
     */
    double y;

    /**
     * mesh spacing parameter associated with the node
     */
    double h;

  } ;


  /**
   * Edge structure
   */
  struct Edge_t
  {
    /**
     * edge index, (sequential, starts at 0)
     */
    int index;

    /**
     * index of the starting coordinate node.
     */
    int point1;

    /**
     * index of the ending coordinate node.
     */
    int point2;

    /**
     * has the value 2 for edges on exposed boundaries, and zero otherwise
     */
    int bcode;

    /**
     * bc index
     */
    short int bc_index;

  };


  /**
   * the '<' operator used in std::map
   */
  struct lt_edge
  {
    bool operator()(const Edge_t &e1, const Edge_t &e2) const
    {
      if(e1.point2 != e2.point2)
        return e1.point2 < e2.point2;

      return e1.point1 < e2.point1;
    }
  };


  /**
   * Region information structure
   */
  struct Region_t
  {

    /**
     * region index (sequential, starts at 0)
     */
    int  index;

    /**
     * the node num of this region
     */
    int  node_num;

    /**
     * the triangle num of this region
     */
    int  tri_num;

    /**
     *  material type
     */
    std::string  type;

    /**
     * name of the region
     */
    std::string name;

    /**
     * edge index for edge on boundary of a region
     */
    std::vector<int> boundary;

  };


  /**
   * Interface or Electrode structure
   */
  struct Interface_t
  {
    /**
     * interface index number (sequential, starts at 0).
     */
    int  index;

    /**
     * Interface type (refers to the interface symbol from the material
     * database file). Can specify an electrode.
     */
    std::string  type;

    /**
     * name of the interface.
     */
    std::string  name;

    /**
     * Region index if an entire region is the electrode.
     * -1 means it is only a boundary
     */
    int  region;

    /**
     * edge index for edge on an interface
     */
    std::vector<int> boundary;
  };


  /**
   * Triangle structure
   */
  struct Tri_t
  {
    /**
     * triangle index (sequential, starts at 0)
     */
    int  index;

    /**
     * region index of the region the triangle is part of
     */
    int  region;

    /**
     * coordinate index of the triangle node 1
     */
    int c1;

    /**
     * coordinate index of the triangle node 2
     */
    int c2;

    /**
     * coordinate index of the triangle node 3
     */
    int c3;

    /**
     * triangle index of neighbor triangle opposite node c1
     * A code of -1024 is used instead of a neighbor triangle index for nodes opposite a
     * reflecting boundary, and -1022 is used instead of a neighbor triangle index for
     * nodes opposite an exposed boundary.
     */
    int  t1;

    /**
     * triangle index of neighbor triangle opposite node c1
     */
    int  t2;

    /**
     * triangle index of neighbor triangle opposite node c1
     */
    int  t3;
  };


  /**
   * material parameters
   */
  struct Parameter_t
  {
    /**
     * region index, the material parameters belongs to which region
     */
    int  region;

    /**
     * a vector for all the parameter's name
     */
    std::vector<std::string> parameter_name_array;

    /**
     * a vector for all the parameter's value
     */
    std::vector<double> parameter_value_array;
  };


  /**
   * solution head structure
   */
  struct SolHead_t
  {
    /**
     * Number of solution variables
     */
    int sol_num;

    /**
     * the name of each solution variables
     */
    std::vector<std::string> sol_name_array;

    /**
     * clear solution head structure
     */
    void clear()
    {
      sol_num = 0;
      sol_name_array.clear();
    }

    /**
     * @return the index of sol_name in sol_name_array
     */
    unsigned int solution_index(const std::string & sol_name) const
    {
       for(unsigned int i=0; i<sol_name_array.size(); i++)
         if( sol_name_array[i] == sol_name )
           return i;
       return static_cast<unsigned int>(-1);
    }
  };

  /**
   * solution data structure
   */
  struct SolData_t
  {
    /**
     * the node index of this solution data belongs to
     */
    int index;

    /**
     * the material of this node
     */
    std::string  material;

    /**
     * solution data array. the data has the same order as sol_name_array in SolHead_t
     */
    std::vector<double> data_array;
  };


  /**
   * dc or transient simulation result
   */
  struct TranSol_t
  {
    /**
     * a string for solution name
     */
    char sol_name[32];

    /**
     * vector for data name
     */
    std::vector<std::string> data_name_array;

    /**
     * vector for data unit
     */
    std::vector<std::string> data_unit_array;

    /**
     * vector for data value
     */
    std::vector<double> data_value_array;

    /**
     * clear simulation result
     */
    void clear()
    {
      data_name_array.clear();
      data_unit_array.clear();
      data_value_array.clear();
    }
  };


  /**
   * structure for compound semiconductor material
   */
  class Component_t
  {
  public:
    int region;
    double xmin,xmax;
    double ymin,ymax;
    int    direction;
    double mole_begin;
    double mole_ratio;
  public:
    void mole(double x, double y, double & mole_x)
    {
      if(x>=xmin-1e-10 && x <=xmax+1e-10 && y>=ymin-1e-10 && y <=ymax+1e-10)
      {
        if(direction == 0) 	 mole_x = mole_begin;
        if(direction == 1)	 mole_x = mole_begin + mole_ratio*(x-xmin);
        if(direction == 2)	 mole_x = mole_begin - mole_ratio*(y-ymax);
        if(mole_x<0) 	         mole_x = 0;
      }
    }
  };

  // parser
  extern FILE * yyin;
  extern int yyparse();

  //  the parser will fill data into these structure
  extern std::vector<Node_t>        node_array;
  extern std::vector<Edge_t>        edge_array;
  extern std::vector<Tri_t>         tri_array;
  extern std::vector<Region_t>      region_array;
  extern std::vector<Interface_t>   interface_array;
  extern std::vector<Component_t>   component_array;
  extern std::vector<Parameter_t>   parameter_array;
  extern SolHead_t                  sol_head;
  extern std::vector<SolData_t>     sol_data;
  extern TranSol_t                  tran_sol;

  // clear all the data
  extern void clear();

}

#endif

