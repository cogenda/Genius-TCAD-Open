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

//  $Id: skeleton.h,v 1.6 2008/07/09 05:58:16 gdiso Exp $

#ifndef _skeleton_h_
#define _skeleton_h_

#include <vector>
#include <map>
#include <string>
#include <math.h>

enum   mole_grad{MOLE_GRAD_X,MOLE_GRAD_Y,MOLE_GRAD_Z};


/**
 * The background (Skeleton) point data class
 */
class SkeletonPoint
{
public:
  /**
   *  the location of this point
   */
  double          x,y,z;

  /**
   * the index of this point in global point array
   */
  unsigned int    index;

  /**
   * indicate if this point has been eliminated
   */
  bool            eliminated;

  /**
   * set the location of point
   */
  void            set_location(double _x,double _y,double _z)
  {x=_x; y=_y; z=_z;}

  /**
   * default constructor
   */
  SkeletonPoint():x(0),y(0),z(0),eliminated(false) {}

  /**
   * constructor
   */
  SkeletonPoint(double _x,double _y, double _z):x(_x),y(_y),z(_z),eliminated(false) {}
};



/**
 * A Skeleton line is the collection of points lie on a line
 */
struct SkeletonLine
{
  /**
   * the array which contains the  Skeleton Point
   */
  std::vector<SkeletonPoint>  point_array;

  /**
   * used to insert point into point_array
   */
  void insert(const SkeletonPoint &p)
  { point_array.push_back(p); }

} ;



/**
 * A Skeleton edge is made up of two points
 * this is used in XY plan 2D mesh generator
 */
struct SkeletonEdge
{
  /**
   * the index of this edge in the global 3D array
   */
  unsigned int    IX,IY,IZ;

  // p1[index_x][index_y][index_z]
  unsigned int    p1[3];

  // p2[index_x][index_y][index_z]
  unsigned int    p2[3];

};



/**
 * the '<' operator used in std::map SkeletonEdge
 */
struct lt_edge
{
  bool operator()(const SkeletonEdge &e1, const SkeletonEdge &e2) const
  {
    unsigned long int space1 = e1.IX*e1.IY;
    unsigned long int space2 = e2.IX*e2.IY;
    return (e1.p1[1]*space1 + e1.p1[0] + e1.p2[1]*space1 + e1.p2[0] <
            e2.p1[1]*space2 + e2.p1[0] + e2.p2[1]*space2 + e2.p2[0]);
  }
};



/**
 * SkeletonSegment:  SkeletonSegment is the collection of SkeletonEdge
 * with same boundary mark (also boundary label)
 */
struct SkeletonSegment
{
  /**
   * the collection of SkeletonEdge
   */
  std::vector<SkeletonEdge>  edge_array;

  /**
   *  the mark of this face
   */
  int            segment_mark;

  /**
   * we should give the segment a label
   */
  std::string    segment_label;

};


/**
 * SkeletonSector:
 */
struct SkeletonSector
{
  
  /**
   * bound box of this face
   */
  double rmin,rmax;
  double thetamin,thetamax;
  double zmin,zmax;

  /**
   *  the mark of this sector
   */
  int    sector_mark;

 
  /**
   * we should give the sector a label
   */
  std::string    sector_label;
};



/**
 * SkeletonRegion2D:  SkeletonRegion2D contains the 2D region information
 *
 */
struct SkeletonRegion2D
{
  /**
   * the user defined label of this region
   */
  std::string    label;

  /**
   * the material type of this region
   */
  std::string    material;

  /**
   * we allow regtangle and ellipse region shape
   */
  int            shape;

  /**
   * half point x location within the region
   * used to specify region material
   */
  std::vector<double> px;
  /**
   * half point y location within the region
   * used to specify region material
   */
  std::vector<double> py;

  // for regtangle
  unsigned int ixmin,ixmax,iymin,iymax;  //bound box
  double xmin,xmax,ymin,ymax;      //bound box

  // for ellipse
  double centrex,centrey;          //the centre of the ellipse
  double major_radii,minor_radii;  //major and minor radii
  double theta;                    //the rotary angle
  int    division;                 //the division number of its boundary
  
  //for sector
  double rmin, rmax;
  double thetamin, thetamax;

  double    mole_x1;              // for compound materials
  double    mole_x1_slope;
  mole_grad mole_x1_grad;

  //  output value
  int              node_num;       //output total node num
  int              tri_num;        //output total triangle num
  std::vector<int> boundary;       //output total boundary

} ;



/**
 * SkeletonQuad:  This struct defines the quadrangle
 * used for model surface in 3D.
 */
struct  SkeletonQuad
{
  // the index of four points in global 1D array
  unsigned int    index1;
  unsigned int    index2;
  unsigned int    index3;
  unsigned int    index4;

  // the pointer of 4 points
  SkeletonPoint * pp1;
  SkeletonPoint * pp2;
  SkeletonPoint * pp3;
  SkeletonPoint * pp4;

  // the index of four points in global 3D array
  unsigned int    p1[3];    // p1[index_z][index_y][index_x]
  unsigned int    p2[3];    // p2[index_z][index_y][index_x]
  unsigned int    p3[3];    // p3[index_z][index_y][index_x]
  unsigned int    p4[3];    // p4[index_z][index_y][index_x]

  //// the center of quadrangle
  double xc,yc,zc;

  // functions for set each point

  void set_p1(unsigned int x, unsigned int y, unsigned int z)
  { p1[0]=x; p1[1]=y; p1[2]=z; }

  void set_p2(unsigned int x, unsigned int y, unsigned int z)
  { p2[0]=x; p2[1]=y; p2[2]=z; }

  void set_p3(unsigned int x, unsigned int y, unsigned int z)
  { p3[0]=x; p3[1]=y; p3[2]=z; }

  void set_p4(unsigned int x, unsigned int y, unsigned int z)
  { p4[0]=x; p4[1]=y; p4[2]=z; }

};



/**
 * lt_quadrangle: the structure contains '<' operator used for stl::map
 * we defined the '<' operator of two quadrangle by its center (X,Y,Z) index
 */
struct lt_quadrangle
{
  bool operator()(const SkeletonQuad &q1, const SkeletonQuad &q2) const
  {
    //the TOLERANCE is 1e-8, I think it is enough.
    if( fabs(q1.zc-q2.zc) > TOLERANCE )
      return  q1.zc<q2.zc;
    else if ( fabs(q1.yc-q2.yc) > TOLERANCE )
      return  q1.yc<q2.yc;
    else if ( fabs(q1.xc-q2.xc) > TOLERANCE)
      return  q1.xc<q2.xc;
    else return false; //meet strict week ordering
  }
};


/**
 * SkeletonFace:  SkeletonFace is the collection of SkeletonQuad
 * with same boundary mark (also boundary label)
 */
struct SkeletonFace
{

  //bound box of this face
  unsigned int ixmin,ixmax,iymin,iymax,izmin,izmax;

  /**
   * the collection of SkeletonQuad
   */
  std::vector<SkeletonQuad>  quad_array;

  /**
   *  the mark of this face
   */
  int    face_mark;

  /**
   * we should give the face a label
   */
  std::string face_label;

  /**
   * true for user defined face
   */
  bool   user_define;

};





/**
 * SkeletonRegion3D:  SkeletonRegion3D contains the 3D region information
 * It is useful to work with tetgen
 */
struct SkeletonRegion3D
{
  //bound box by SkeletonLine
  unsigned int ixmin,ixmax,iymin,iymax,izmin,izmax;

  //bound box by coordinate
  double xmin,xmax,ymin,ymax,zmin,zmax;

  //for sector
  double rmin, rmax;
  double thetamin, thetamax;

  std::vector<double> px;             //half point x location within the region
  std::vector<double> py;             //half point y location within the region
  std::vector<double> pz;             //half point z location within the region

  /**
   * the user defined label of this region
   */
  std::string label;

  /**
   * the material of this region
   */
  std::string material;

  // for compound materials
  double mole_x1;
  double mole_x1_slope;
  mole_grad mole_x1_grad;

  //  output value
  int    node_num;             //output
  int    tetra_num;            //output
  std::vector<int> boundary;        //output

};

#endif
