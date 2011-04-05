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

#ifndef __tif3d_h__
#define __tif3d_h__

#include <string>
#include <vector>
#include <map>


/**
 * read 3D TIF file
 */
class TIF3D
{
public:

  /**
   * constructor
   */
  TIF3D(const std::string & file);

  /**
   * destroy internal data?
   */
  ~TIF3D() {}

  /**
   * read the silvaco tif file
   */
  int read();


public:
  /** Nodal data structure */
  struct Node_t
  {
    /** node index, (sequential, starts at 0) */
    int    index;

    /** X-coordinate of node  */
    double x;

    /** Y-coordinate of node  */
    double y;

    /** Z-coordinate of node  */
    double z;

    /** mesh spacing parameter associated with the node  */
    double h;
  };


  /**
   * Face structure
   */
  struct Face_t
  {
    /**
     * face index, (sequential, starts at 0)
     */
    int index;

    /**
     * index of the first coordinate node.
     */
    int point1;

    /**
     * index of the second coordinate node.
     */
    int point2;

    /**
     * index of the third coordinate node.
     */
    int point3;

    /**
     * has the value 2 for face on exposed boundaries, and zero otherwise
     */
    int bcode;

    /**
     * bc index
     */
    int bc_index;

  };


  /**
   * the '<' operator used in std::map
   */
  struct lt_face
  {
    bool operator()(const Face_t &_f1, const Face_t &_f2) const
    {
      Face_t f1 = _f1;
      Face_t f2 = _f2;

      if ( f1.point1 > f1.point2 ) std::swap ( f1.point1, f1.point2 );
      if ( f1.point2 > f1.point3 )
      {
        if ( f1.point3 < f1.point1 ) std::swap ( f1.point1, f1.point3 );
        std::swap ( f1.point2 , f1.point3 );
      }

      if ( f2.point1 > f2.point2 ) std::swap ( f2.point1, f2.point2 );
      if ( f2.point2 > f2.point3 )
      {
        if ( f2.point3 < f2.point1 ) std::swap ( f2.point1, f2.point3 );
        std::swap ( f2.point2 , f2.point3 );
      }

      if(f1.point1 != f2.point1)
        return f1.point1 < f2.point1;
      else if(f1.point2 != f2.point2)
        return f1.point2 < f2.point2;
      return f1.point3 < f2.point3;
    }
  };


  /** Tet structure */
  struct Tet_t
  {
    /** tetrahedron index (sequential, starts at 0) */
    int  index;

    /** region index of the region the tetrahedron is part of */
    int  region;

    /** coordinate index of the tetrahedron node 1 */
    int c1;

    /** coordinate index of the tetrahedron node 2 */
    int c2;

    /** coordinate index of the tetrahedron node 3 */
    int c3;

    /** coordinate index of the tetrahedron node 4 */
    int c4;

    /**
     * tetrahedron index of neighbor triangle opposite node c1
     * A negative code is used instead of a neighbor tetrahedron index for nodes opposite an exposed boundary.
     */
    int  t1;

    /** tetrahedron index of neighbor triangle opposite node c2 */
    int  t2;

    /** tetrahedron index of neighbor triangle opposite node c3 */
    int  t3;

    /** tetrahedron index of neighbor triangle opposite node c4 */
    int  t4;
  };


  /** Region information structure */
  struct Region_t
  {
    /** region index (sequential, starts at 0) */
    int  index;

    /** the node num of this region */
    int  node_num;

    /** the tetrahedron num of this region */
    int  tet_num;

    /** material type */
    std::string  material;

    /** name of the region */
    std::string name;

    /** group of the region */
    std::string group;

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
     * the region label of this node
     */
    int  region_index;

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

  const std::vector<Face_t> & tif_faces() const
  { return  _faces; }

  std::vector<Face_t> & tif_faces()
  { return  _faces; }

  const std::vector<Tet_t> & tif_tets() const
    { return  _tets; }

  std::vector<Tet_t> & tif_tets()
  { return  _tets; }

  const std::vector<Region_t> & region_array() const
    { return _regions; }

  std::vector<Region_t> & region_array()
  { return _regions; }

  const std::vector<SolData_t> & sol_data_array() const
    { return _sol_data; }

  std::vector<SolData_t> & sol_data_array()
  { return _sol_data; }

public:

  const Face_t & face(unsigned int n) const
  { return _faces[n]; }

  Face_t & face(unsigned int n)
  { return _faces[n]; }

  const Tet_t & tet(unsigned int n) const
    { return _tets[n]; }

  Tet_t & tet(unsigned int n)
  { return _tets[n]; }

  const Region_t & region(unsigned int n) const
    { return _regions[n]; }

  bool face_has_label(int bc_index) const
  { return _face_labels.find(bc_index) != _face_labels.end(); }

  const std::string & face_label(int bc_index) const
  { return _face_labels.find(bc_index)->second; }

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

  std::vector<Face_t>    _faces;

  std::map<int, std::string> _face_labels;

  std::vector<Tet_t>     _tets;

  std::vector<Region_t>  _regions;

  SolHead_t              _sol_head;

  std::vector<SolData_t> _sol_data;

  std::vector<int>       _electrode_info;

};

#endif // #define __tif3d_h__

