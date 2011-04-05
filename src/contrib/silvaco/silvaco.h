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

#ifndef __silvaco_h__
#define __silvaco_h__

#include <string>
#include <vector>
#include <map>


/**
 * read silvaco TIF file
 */
class SilvacoTIF
{
public:

  /**
   * constructor
   */
  SilvacoTIF(const std::string & file);

  /**
   * destroy internal data?
   */
  ~SilvacoTIF() {}

  /**
   * read the silvaco tif file
   */
  int read();

  //void write_vtk();

public:
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
     * Z-coordinate of node
     */
    double z;
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
   * Region information structure,
   * NOTE: when tri_num==0, region is a boundary segment
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
     * the region is an electrode
     */
    bool electrode;

    /**
     * the region is segment
     */
    bool segment;

    /**
     *  material
     */
    std::string  material;

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
   * solution head structure
   */
  struct SolHead_t
  {
    SolHead_t():sol_num(0) {}

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

public:

  const std::vector<Node_t> & tif_nodes() const
  { return  _nodes; }

  std::vector<Node_t> & tif_nodes()
  { return  _nodes; }

  const std::vector<Edge_t> & tif_edges() const
  { return  _edges; }

  std::vector<Edge_t> & tif_edges()
  { return  _edges; }

  const std::vector<Tri_t> & tif_tris() const
  { return  _tris; }

  std::vector<Tri_t> & tif_tris()
  { return  _tris; }

  const std::vector<Region_t> & region_array() const
  { return _regions; }

  std::vector<Region_t> & region_array()
  { return _regions; }

  const std::vector<SolData_t> & sol_data_array() const
  { return _sol_data; }

  std::vector<SolData_t> & sol_data_array()
  { return _sol_data; }

public:
  const Edge_t & edge(unsigned int n) const
  { return _edges[n]; }

  Edge_t & edge(unsigned int n)
  { return _edges[n]; }

  const Tri_t & tri(unsigned int n) const
  { return _tris[n]; }

  Tri_t & tri(unsigned int n)
  { return _tris[n]; }

  const Region_t & region(unsigned int n) const
  { return _regions[n]; }

  Region_t & region(unsigned int n)
  { return _regions[n]; }

  const SolHead_t & sol_head() const
  { return _sol_head; }

  SolHead_t & sol_head()
  { return _sol_head; }

  const SolData_t & sol_data(unsigned int n) const
  { return _sol_data[n]; }

  SolData_t & sol_data(unsigned int n)
  { return _sol_data[n]; }

private:

  const std::string      _file;

  std::vector<Node_t>    _nodes;

  std::vector<Edge_t>    _edges;

  std::vector<Tri_t>     _tris;

  std::vector<Region_t>  _regions;

  SolHead_t              _sol_head;

  std::vector<SolData_t> _sol_data;

  std::vector<int>       _electrode_info;

  /**
   * implant map
   */
  std::map<int, std::string> SilImp;

  /**
   * material map
   */
  std::map<int, std::string> SilMat;

};

#endif // #define __silvaco_h__

