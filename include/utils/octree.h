#ifndef __octree_h__
#define __octree_h__

#include <vector>
#include <map>
#include <set>
#include <string>

#include "point.h"
#include "stl_tree.h"


/*	      y
              |top
              |
              |bot
              |_______________ x
             /  left  right
            / back
           /
        z / front
*/


class OcTreeLocation
{
public:

  enum Location{L_Top, L_Bottom, L_Right, L_Left, L_Front, L_Back, ROOT};

  /**
   * constructor
   */
  OcTreeLocation ( Location x=ROOT, Location y=ROOT, Location z=ROOT )
      :_x ( x ), _y ( y ), _z ( z )
  {
    assert(_x == L_Right || _x == L_Left    || _x == ROOT);
    assert(_y == L_Top   || _y == L_Bottom  || _y == ROOT);
    assert(_z == L_Front || _z == L_Back    || _z == ROOT);
  }

  /**
   * copy constrctor
   */
  OcTreeLocation ( const OcTreeLocation & other )
      :_x ( other._x ), _y ( other._y ), _z ( other._z )
  {}

  /**
   * constrctor by key
   */
  OcTreeLocation ( int k )
      :_x ( ROOT ), _y ( ROOT ), _z ( ROOT )
  {
    _x = k%2 ? L_Left:L_Right;
    _y = (k>>1)%2 ? L_Bottom:L_Top;
    _z = (k>>2)%2 ? L_Back:L_Front;
    assert(this->key()==k);
  }


  void set_x(Location loc) { _x = loc; }
  void set_y(Location loc) { _y = loc; }
  void set_z(Location loc) { _z = loc; }

  static OcTreeLocation x_symmetry ( const OcTreeLocation & location )
  {
    OcTreeLocation symmetry = location;
    if ( location._x == L_Right )
      symmetry._x = L_Left;
    else if ( location._x == L_Left )
      symmetry._x = L_Right;
    return symmetry;
  }

  static OcTreeLocation y_symmetry ( const OcTreeLocation & location )
  {
    OcTreeLocation symmetry = location;
    if ( location._y == L_Top )
      symmetry._y = L_Bottom;
    else if ( location._y == L_Bottom )
      symmetry._y = L_Top;
    return symmetry;
  }

  static OcTreeLocation z_symmetry ( const OcTreeLocation & location )
  {
    OcTreeLocation symmetry = location;
    if ( location._z == L_Front )
      symmetry._z = L_Back;
    else if ( location._z == L_Back )
      symmetry._z = L_Front;
    return symmetry;
  }

  /**
   * @return true if this node has a location
   */
  bool has_location ( Location loc ) const  { return ( _x==loc || _y==loc || _z==loc); }

  /**
   * @return true when 2 OcTreeLocation equals
   */
  bool operator == ( const OcTreeLocation & other ) const     { return ( _x == other._x && _y == other._y && _z == other._z ); }

  /**
   * @return true when OcTreeLocation is valid
   */
  bool is_valid() const    {return ( _x != ROOT && _y != ROOT && _z != ROOT );}

  /**
   * OcTreeLocation has a uniform key. from valid:0-7, invalid -1
   */
  int key() const
  {
    if( !is_valid() ) return -1;
    if(_x==L_Left && _y==L_Bottom && _z==L_Front)  return 0;
    if(_x==L_Right && _y==L_Bottom && _z==L_Front) return 1;
    if(_x==L_Right && _y==L_Bottom && _z==L_Back)  return 2;
    if(_x==L_Left && _y==L_Bottom && _z==L_Back)   return 3;
    if(_x==L_Left && _y==L_Top && _z==L_Front)     return 4;
    if(_x==L_Right && _y==L_Top && _z==L_Front)    return 5;
    if(_x==L_Right && _y==L_Top && _z==L_Back)     return 6;
    if(_x==L_Left && _y==L_Top && _z==L_Back)      return 7;

    return -1;
  }

  void print ( std::ostream& os ) const
  {
    os <<"(";
    if ( _x == L_Right )  os << "Right ";
    if ( _x == L_Left )   os << "Left ";
    if ( _y == L_Top )    os << "Top ";
    if ( _y == L_Bottom ) os << "Bottom ";
    if ( _z == L_Front )  os << "Front";
    if ( _z == L_Back )   os << "Back";
    os<<")";
    //os<<" key="<<key()<<std::endl;
  }


  friend std::ostream& operator << (std::ostream& os, const OcTreeLocation& t)
  {
    t.print(os);
    return os;
  }

private:
  Location _x;
  Location _y;
  Location _z;
};


class OcTreeNode;


/**
 * empty octree data
 */
class OcTreeDataBase
{
public:
  OcTreeDataBase() {}

  virtual ~OcTreeDataBase() {}

  virtual std::vector<OcTreeDataBase *> subdivide(const OcTreeNode &) = 0;

  virtual bool refine(const OcTreeNode & , const std::vector<OcTreeNode> &) const { return false; }

  virtual double value(const std::string &) const { return 0.0; }
};


/**
 * octree node, light weighted class, can be copy constructed safely
 */
class OcTreeNode
{
public:

  /**
   * empty constructor
   */
  OcTreeNode (const OcTreeLocation &location=OcTreeLocation() )
      :_location ( location ), _data(0), _divide_flag ( false )
  {
    for(int i=0; i<8; ++i)
      _p[i] = 0;
  }

  /**
   * copy constructor
   */
  OcTreeNode ( const OcTreeNode & other )
  {
    for(int i=0; i<8; ++i)
      _p[i] = other._p[i];
    _location = other._location;
    _divide_flag = other._divide_flag;
    _data = other._data;
  }


  /**
   * set octree node data
   */
  void set_data(OcTreeDataBase * data) { _data = data; }


  /**
   * @return octree node data
   */
  OcTreeDataBase * data() const { return _data; }


  /**
   * set the point at corner of octree node by OcTreeLocation key
   */
  void set_corner_point(const Point *p, const OcTreeLocation &loc )  { _p[loc.key()] = p; }


  /**
   * get the corner point at given location
   */
  const Point *  get_point(const OcTreeLocation &loc) const     { return _p[loc.key()]; }


  /**
   * get the center point at given location
   */
  Point  get_center_point(const OcTreeLocation &loc1, const OcTreeLocation &loc2) const
    { return 0.5*(*_p[loc1.key()] + *_p[loc2.key()]); }


  /**
   * get the center point at given location
   */
  Point  get_center_point() const
  {
    Point p1 = (*get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top, OcTreeLocation::L_Back)));
    Point p2 = (*get_point(OcTreeLocation(OcTreeLocation::L_Left, OcTreeLocation::L_Bottom, OcTreeLocation::L_Front)));
    return 0.5*(p1 + p2);
  }


  /**
   * @return the location of octree node
   */
  const OcTreeLocation & get_location() const    { return _location; }

  /**
   * size in x direction
   */
  Real x_size() const
  {
    Point p1 = (*get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top, OcTreeLocation::L_Front)));
    Point p2 = (*get_point(OcTreeLocation(OcTreeLocation::L_Left, OcTreeLocation::L_Top, OcTreeLocation::L_Front)));
    return (p1-p2).coord(0);
  }


  /**
   * size in y direction
   */
  Real y_size() const
  {
    Point p1 = (*get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top, OcTreeLocation::L_Front)));
    Point p2 = (*get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Bottom, OcTreeLocation::L_Front)));
    return (p1-p2).coord(1);
  }


  /**
   * size in z direction
   */
  Real z_size() const
  {
    Point p1 = (*get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top, OcTreeLocation::L_Front)));
    Point p2 = (*get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Top, OcTreeLocation::L_Back)));
    return (p1-p2).coord(2);
  }


  /**
   * volume of octree leaf
   */
  Real volume() const
  {
    return x_size()*y_size()*z_size();
  }


  bool has_point ( const Point & point ) const
  {
    const Point * b = get_point(OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
    const Point * t = get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top, OcTreeLocation::L_Front));

    return ( point.x() >= b->x()  &&
             point.x() <= t->x()  &&
             point.y() >= b->y()  &&
             point.y() <= t->y()  &&
             point.z() >= b->z()  &&
             point.z() <= t->z() );
  }


  /**
   * @return point location in its sub node 
   */
  OcTreeLocation get_sub_location(const Point & point ) const
  {
    if(!has_point(point)) return OcTreeLocation();

    const Point c = this->get_center_point();
    OcTreeLocation loc;
    loc.set_x( point.x() < c.x() ? OcTreeLocation::L_Left : OcTreeLocation::L_Right );
    loc.set_y( point.y() < c.y() ? OcTreeLocation::L_Bottom : OcTreeLocation::L_Top );
    loc.set_z( point.z() < c.z() ? OcTreeLocation::L_Back : OcTreeLocation::L_Front );
    return loc;
  }

  // FIXME http://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms
  bool hit(const Point &p, const Point &dir, std::pair<double, double>& t)
  {
    // if point p inside the bounding_box
    bool inside = has_point(p);

    Point upp = (*get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top, OcTreeLocation::L_Front)));
    Point low = (*get_point(OcTreeLocation(OcTreeLocation::L_Left, OcTreeLocation::L_Bottom, OcTreeLocation::L_Back)));

    // point p outside the bounding_box
    double tmin, tmax;
    double txmin = -std::numeric_limits<double>::infinity();
    double tymin = -std::numeric_limits<double>::infinity();
    double tzmin = -std::numeric_limits<double>::infinity();
    double txmax = std::numeric_limits<double>::infinity();
    double tymax = std::numeric_limits<double>::infinity();
    double tzmax = std::numeric_limits<double>::infinity();

    double eps = 1e-10*(x_size() + y_size() + z_size());

    // assume p outside the bbox
    if( dir.x() != 0.0 )
    {
      double divx = 1 / dir.x();
      if (divx >= 0)
      {
        txmin = (low.x()  - p.x()) * divx;
        txmax = (upp.x()  - p.x()) * divx;
      }
      else
      {
        txmin = (upp.x() - p.x()) * divx;
        txmax = (low.x() - p.x()) * divx;
      }
    }
    else
    {
      if( p.y() < low.y() || p.y() >upp.y() ) return false;
      if( p.z() < low.z() || p.z() >upp.z() ) return false;
    }

    if( dir.y() != 0.0 )
    {
      double divy = 1 / dir.y();
      if (divy >= 0)
      {

        tymin = (low.y() - p.y()) * divy;
        tymax = (upp.y() - p.y()) * divy;
      }
      else
      {
        tymin = (upp.y() - p.y()) * divy;
        tymax = (low.y() - p.y()) * divy;
      }
    }
    else
    {
      if( p.x() < low.x() || p.x() >upp.x() ) return false;
      if( p.z() < low.z() || p.z() >upp.z() ) return false;
    }

    if( dir.z() != 0.0 )
    {
      double divz = 1 / dir.z();
      if (divz >= 0)
      {
        tzmin = (low.z() - p.z()) * divz;
        tzmax = (upp.z() - p.z()) * divz;
      }
      else
      {
        tzmin = (upp.z() - p.z()) * divz;
        tzmax = (low.z() - p.z()) * divz;
      }
    }
    else
    {
      if( p.y() < low.y() || p.y() >upp.y() ) return false;
      if( p.x() < low.x() || p.x() >upp.x() ) return false;
    }

    //std::cout<< txmin << " " << txmax << std::endl;
    //std::cout<< tymin << " " << tymax << std::endl;
    //std::cout<< tzmin << " " << tzmax << std::endl;

    if(inside)
    {
      tmin = std::min(-txmin, -tymin);
      tmin = std::min(tmin,   -tzmin);

      tmax = std::max(-txmax, -tymax);
      tmax = std::max(tmax,   -tzmax);

      t.first  = -tmin;
      t.second = -tmax;
      return true;
    }
    else
    {
      if ( (txmin > tymax) || (tymin > txmax) )
        return false;

      tmin = std::max(txmin, tymin);
      tmax = std::min(txmax, tymax);

      if ( (tmin > tzmax) || (tzmin > tmax) )
        return false;

      tmin = std::max(tmin,  tzmin);
      tmax = std::min(tmax,  tzmax);

      if ( tmin >= tmax ) return false;
    }

    // the ray start point is on boundbox, and it is leaving the bounding box
    if(tmin<eps && tmax<0) return false;

    if ( (tmin >= -eps) && (tmax < 1e+30) )
    {
      t.first  = tmin;
      t.second = tmax;
      return true;
    }

    return false;
  }


  bool divide_flag() const { return _divide_flag; }



  bool & divide_flag()     { return _divide_flag; }

  void print ( std::ostream& os ) const
  {
    Point upp = (*get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top, OcTreeLocation::L_Front)));
    Point low = (*get_point(OcTreeLocation(OcTreeLocation::L_Left, OcTreeLocation::L_Bottom, OcTreeLocation::L_Back)));

    os << "octree leaf node begin" << std::endl;
    os << upp;
    os << low;
    os << "octree leaf node end " << std::endl;
  }


  friend std::ostream& operator << (std::ostream& os, const OcTreeNode& t)
  {
    t.print(os);
    return os;
  }


private:

  /**
   * all the points;
   */
  const Point * _p[8];


  /**
   * the location of octree leaf
   */
  OcTreeLocation _location;


  /**
   * base class of octree data
   */
  OcTreeDataBase * _data;


  /**
   * indicate if this leaf should be divide
   */
  bool _divide_flag;

};



class OcTree : public tree<OcTreeNode>
{
public:

  /**
  * construct with bound box
   */
  OcTree(const Point &b1, const Point &b2, OcTreeDataBase * data = 0, int depth=0);


  /**
  * constructor, set root node
   */
  OcTree ( const OcTreeNode & root_node );

  /**
  * free _points array
   */
  ~OcTree();

  /**
   * refine the octree, when OcTreeDataBase::refine is true
   */
  void refine();


  typedef typename tree< OcTreeNode >::iterator_base     tree_iterator_base;
  typedef typename tree< OcTreeNode >::leaf_iterator     tree_leaf_iterator;
  typedef typename tree< OcTreeNode >::sibling_iterator  tree_sibling_iterator;

  /**
  * add point to quadtree.
  * we will first search _points array to see if \p point already exist
  * or we will create new Point<T> by the \p point and strore the new created point in _points
  * for both situation, return the pointer point to  _points
   */
  const Point * add_point ( const Point & point );


  /**
  * divide leaf into 8 child leaf
   */
  void subdivide ( tree_iterator_base & leaf_it );

  /**
  * balance the quadtree (neighbor leaf only have one depth difference )
   */
  void balance();


  /**
  * return the nearest octree point to given point p
   */
  Point get_nearest_point ( const Point & p ) const;


  /**
   * write octree mesh in vtk format
   */
  void export_vtk ( const std::string &file );

  /**
   * return a unique key for octree leaf
   */
  std::string key( const tree_iterator_base & leaf_it) const;

  /**
   * print the absolute path of an iterator
   */
  void print_path ( const tree_iterator_base & it ) const;

  /**
   * find the leaf which has point p 
   */
  tree_iterator_base find_leaf_has_point ( const Point & p ) const;


  /**
   * get all the OcTreeNode intersect with edge (p1, p2)
   */
  void intersect(const Point &p1, const Point &p2, std::vector<std::pair<OcTreeNode, double> > & );


  /**
  * find the neighbor of leaf in the direction  x
  * return empty iterator if no find
   */
  tree_iterator_base find_neighbor ( const tree_iterator_base & it, const OcTreeLocation::Location x ) const;

  /**
  * travel to child by given location
   */
  tree_iterator_base goto_octree_child ( const tree_iterator_base & it, const OcTreeLocation & location ) const;

private:

  std::vector<const Point *> _points;

  struct Point_Less
  {
    bool operator() ( const Point *v1, const Point *v2 ) const
    {
      return ( *v1 ) < ( *v2 );
    }
  };

  std::map<const Point *, unsigned int, Point_Less> _point_to_id;


  void insert_point( const Point * p );


  std::set<OcTreeDataBase *> _data;


};

#endif
