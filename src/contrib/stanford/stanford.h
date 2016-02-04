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

#ifndef __stanford_tif_h__
#define __stanford_tif_h__

#include <string>
#include <vector>
#include <map>


/**
 * base class
 * read structure file from stanford TCAD or its derived code (medici, silvaco, etc)
 */
class StanfordTIF
{
public:


  /**
   * constructor
   */
  StanfordTIF() {}

  /**
   * constructor
   */
  StanfordTIF(const std::string & file):_file(file) {}

  /**
   * destroy internal data?
   */
  virtual ~StanfordTIF() {}

  /**
   * read file
   */
  virtual bool read(std::string &err) {}

  /**
   * DIM, most TIF is 2D, except silvaco has 3D extersion
   */
  virtual unsigned int dim() { return 2; }

  /**
   * sync mesh
   */
  //void broadcast_mesh(unsigned int root = 0);

  /**
   * sync  data
   */
  void broadcast(unsigned int root = 0);


  /// debug only
  void export_vtk(const std::string & file) const;

  
  virtual void export_tif(const std::string & file) const {}
  
  
  virtual void export_sup(const std::string & file) const {}

  /**
   * merge several tif files to one
   */
  static StanfordTIF * merge( const std::vector<const StanfordTIF *> &);

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
     * coordinate index of the triangle node 4 (quad element), or -1
     */
    int c4;

    /**
     * triangle index of neighbor triangle opposite node c1
     * A code of -1024 is used instead of a neighbor triangle index for nodes opposite a
     * reflecting boundary, and -1022 is used instead of a neighbor triangle index for
     * nodes opposite an exposed boundary.
     */
    int  t1;

    /**
     * triangle index of neighbor triangle opposite node c2
     */
    int  t2;

    /**
     * triangle index of neighbor triangle opposite node c3
     */
    int  t3;

    /**
     * triangle index of neighbor triangle opposite node c4 (quad element), or -1
     */
    int  t4;


    /** ? */
    int father;

    /** ? */
    int offspr;
  };


  /**
   * Prism structure
   */
  struct Prism_t
  {

    /**
     * Prism index (sequential, starts at 0)
     */
    int  index;

    /**
     * region index
     */
    int  region;

    /**
     * tri index, the Prism made of
     */
    int tri;

    /**
     * top plane
     */
    int z1;

    /**
     * bot plane
     */
    int z2;

  };


  /**
   * Region information structure,
   * NOTE: when tri_num==0, region is a boundary segment
   */
  struct Region_t
  {
  public:
    /**
     * if an entire region is the electrode.
     * -1 means it is only a boundary
     */
    bool is_interface() const
    {
      if (-1==region)
        return false;
      return true;
    }

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
     * the region is an electrode, silvaco tif
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
     * Region index if an entire region is the electrode.
     * -1 means it is only a boundary, medici tif
     */
    int  region;

    /**
     * edge index for edge on boundary of a region
     */
    std::vector<int> boundary;

    int z1;

    int z2;
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
     * the unit of each solution variables
     */
    std::vector<std::string> sol_unit_array;

    /**
     * clear solution head structure
     */
    void clear()
    {
      sol_num = 0;
      sol_name_array.clear();
      sol_unit_array.clear();
    }

    /**
     * @return true when solution exist
     */
    bool has_solution(const std::string & sol_name) const
    {
      for(unsigned int i=0; i<sol_name_array.size(); i++)
        if( sol_name_array[i] == sol_name )
          return true;
      return false;
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

    /**
     * @return the unit of sol_name in sol_name_array
     */
    std::string solution_unit(const std::string & sol_name) const
    {
      unsigned int index = solution_index(sol_name);
      if( index != static_cast<unsigned int>(-1) && index < sol_unit_array.size() )
        return  sol_unit_array[index];
      return std::string("1.0");
    }
    
    std::vector<std::string> extra_infos;
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
     * convert material to region index
     */
    int region_index;

    /**
     * z plane index
     */
    int zplane;

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

  const std::vector<Prism_t> & tif_prisms() const
    { return  _prisms; }

  const std::vector<Region_t> & region_array() const
    { return _regions; }

  std::vector<Region_t> & region_array()
    { return _regions; }

  const SolHead_t & sol_head() const
    { return _sol_head; }

  SolHead_t & sol_head()
  { return _sol_head; }

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

  const Prism_t & prism(unsigned int n) const
  { return _prisms[n]; }

  Prism_t & prism(unsigned int n)
  { return _prisms[n]; }

  int nz() const
  { return _z_slice.size(); }

  double z(unsigned int n) const
  { return _z_slice[n]; }

  const Region_t & region(unsigned int n) const
    { return _regions[n]; }

  Region_t & region(unsigned int n)
  { return _regions[n]; }

  const SolData_t & sol_data(unsigned int n) const
    { return _sol_data[n]; }

  SolData_t & sol_data(unsigned int n)
  { return _sol_data[n]; }


public:

  double solution(unsigned int sol_index, unsigned int index) const;
  
  double solution(unsigned int sol_index, unsigned int r, unsigned int n) const;

  /// return the acceptor in sol_data
  virtual double acceptor(unsigned int index) const;
  /// return the donor in sol_data
  virtual double donor(unsigned int index) const;

  /// return the mole_x in sol_data
  virtual double mole_x(unsigned int index) const;
  /// return the mole_y in sol_data
  virtual double mole_y(unsigned int index) const;


  /// return the acceptor in sol_data by region index and node index
  virtual double acceptor(unsigned int r, unsigned int n) const;
  /// return the donor in sol_data by region index and node index
  virtual double donor(unsigned int r, unsigned int n) const;

  /// return the mole_x in sol_data by region index and node index
  virtual double mole_x(unsigned int r, unsigned int n) const;
  /// return the mole_y in sol_data by region index and node index
  virtual double mole_y(unsigned int r, unsigned int n) const;

protected:

  const std::string      _file;

  std::vector<Node_t>    _nodes;

  std::vector<Edge_t>    _edges;

  std::vector<Tri_t>     _tris;

  std::vector<double>    _z_slice;

  std::vector<Prism_t>   _prisms;

  std::vector<Region_t>  _regions;

  SolHead_t              _sol_head;

  std::vector<SolData_t> _sol_data;

  std::map<std::string, unsigned int> _solution_index_map;
  
  std::map<int, std::string> _region_material; 

  /// map (node_index, material) to data_index
  std::map<std::pair<unsigned int, std::string>, unsigned int> _solution_map;

  /**
   * make sure each region/boundary has unique name
   */
  bool _check_duplicate_name() const;

  /**
   * solution has a material string to indicate which region it belongs to
   */
  void _find_solution_region_by_material();

};

#endif // #define __stanford_tif_h__

