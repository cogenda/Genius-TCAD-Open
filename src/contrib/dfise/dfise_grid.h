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


#include <vector>
#include <set>



namespace DFISE
{

  /**
   * DFISE Point
   */
  struct Point
  {
    Point(double x=0.0, double y=0.0, double z=0.0 )
    {
      coords[0] = x;
      coords[1] = y;
      coords[2] = z;
    }

    double operator [] (const unsigned int i) const
    { return coords[i]; }

    double & operator [] (const unsigned int i)
    { return coords[i]; }

    Point operator + (const Point &p) const
    {
      return Point(coords[0] + p.coords[0], coords[1] + p.coords[1], coords[2] + p.coords[2]);
    }

    Point operator - (const Point &p) const
    {
      return Point(coords[0] - p.coords[0], coords[1] - p.coords[1], coords[2] - p.coords[2]);
    }

    Point operator - () const
    {
      return Point(-coords[0], -coords[1], -coords[2]);
    }

    double  dot (const Point &p) const
    { return (coords[0]*p.coords[0] +  coords[1]*p.coords[1] + coords[2]*p.coords[2]);  }

    Point  cross (const Point &p) const
    {
      return Point(coords[1]*p.coords[2] - coords[2]*p.coords[1],
                  -coords[0]*p.coords[2] + coords[2]*p.coords[0],
                   coords[0]*p.coords[1] - coords[1]*p.coords[0]);
    }

    double coords[3];

  };

  /**
   * @return ture when 2D point p1-p2-p3 in ccw
   */
  inline bool is_ccw_2d(const Point &p1,  const Point &p2, const Point &p3)
  {
    double ax = p2.coords[0] - p1.coords[0];
    double ay = p2.coords[1] - p1.coords[1];

    double bx = p3.coords[0] - p1.coords[0];
    double by = p3.coords[1] - p1.coords[1];

    return ( ax*by - ay*bx > 0);
  }

  /**
   * @return true when p above the plane made by p1-p2-p3
   */
  inline bool is_above(const Point &p1,  const Point &p2, const Point &p3, const Point &p)
  {
    Point n  = (p2 - p1).cross(p3 - p1);
    double proj = (p - p1).dot(n);
    return (proj > 0.);
  }


  /**
   * DFISE element
   */
  struct Element
  {
    int  elem_code;            // DFISE element type
    std::vector<int> faces;    // face index
    std::vector<int> vertices; // vertex index
    int region_index;
    int bc_index;
  };

  /**
   * data block of grid/boundary file
   */
  class GRID
  {
  public:
    /**
     * constructor
     */
    GRID();

    ~GRID() { this->clear(); }

    /**
     * clear the struct
     */
    void clear();

    /**
     * add edge
     */
    void add_edge(const std::pair<int, int> & edge)
    {
      Edges.push_back(edge);
      _edge_set.insert(edge);
    }


    /**
     * @return the nodes belongs to a face in special order
     */
    std::vector<int> get_face_nodes(int face_index);

    /**
     * build element node index
     */
    void build_node(Element & elem);


    /**
     * output
     */
    void print( std::ostream & out ) const;

  public:

    unsigned int dimension;

    double translate[3];

    double transform[9];

    std::vector< Point > Vertices;

    std::vector< std::pair<int, int> > Edges;

    std::vector< std::vector<int> >    Faces;

    std::vector< char >                Locations;

    std::vector< Element >             Elements;

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
     * elements of each region
     */
    std::vector< std::vector<unsigned int> > region_elements;


  private:

    /**
     * the '<' operator used in set/map key
     */
    struct _lt_pair_int
    {
      bool operator()(const std::pair<int, int> &e1, const std::pair<int, int> &e2) const
      {
        int e1_min = std::min(e1.first, e1.second);
        int e1_max = std::max(e1.first, e1.second);
        int e2_min = std::min(e2.first, e2.second);
        int e2_max = std::max(e2.first, e2.second);

        if( e1_min != e2_min ) return e1_min < e2_min;
        return e1_max < e2_max;
      }
    };

    std::set< std::pair<int, int>,  _lt_pair_int> _edge_set;

    /**
     * build each element
     */
    void build_edge2_node(Element & elem);
    void build_tri3_node(Element & elem);
    void build_quad4_node(Element & elem);
    void build_tet4_node(Element & elem);
    void build_pyramid5_node(Element & elem);
    void build_prism6_node(Element & elem);
    void build_hex8_node(Element & elem);

  };

}

