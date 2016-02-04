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

//  $Id: mesh_generation.h,v 1.9 2008/07/09 05:58:16 gdiso Exp $


#ifndef __mesh_generation_h__
#define __mesh_generation_h__



// C++ Includes   -----------------------------------
#include <vector>

// Local Includes -----------------------------------
#include "genius_common.h"

#include "enum_elem_type.h" // needed for ElemType enum
#include "point.h"
#include "parser.h"
#include "skeleton.h"
#include "mesh_generation_base.h"

class MeshGeneratorStruct : public MeshGeneratorBase
{
public:

  /**
   * Constructor.
   */
  MeshGeneratorStruct (MeshBase& mesh)
  : MeshGeneratorBase(mesh), dscale(1e2), point_num(0), point_array3d(0),
  IX(0), IY(0), IZ(0), h_min(1e30), h_max(0.0)
  {}


  /**
   * destructor.
   */
  virtual ~MeshGeneratorStruct() {}

protected:
  /**
   * scale the dimension of geometry
   * the unit of dimension in input file should be um
   * default scale is multiply 1e2 to um
   */
  Real dscale;

  /**
   * the total point number
   */
  int  point_num;

  /**
   * the 3D point array which makes the background point clouds
   */
  SkeletonPoint ***point_array3d;

  /**
   * the bound box of 3D points array by index
   */
  unsigned int  IX,IY,IZ;

  /**
   * the bound box of point clouds by coordinate
   */
  double xmin,xmax,ymin,ymax,zmin,zmax;

  /**
   * the min/max grid size
   */
  double h_min, h_max;

  /**
   * the SkeletonLine in x dim
   */
  SkeletonLine   skeleton_line_x;

  /**
   * the SkeletonLine in y dim
   */
  SkeletonLine   skeleton_line_y;

  /**
   * the SkeletonLine in z dim
   */
  SkeletonLine   skeleton_line_z;

  /**
   * set_x_line:  This function check and do X.MESH card
   * which init mesh line in x direction.
   */
  int set_x_line(const Parser::Card &c);

  /**
   * set_y_line:  This function check and do Y.MESH card
   * which init mesh line in y direction.
   */
  int set_y_line(const Parser::Card &c);

  /**
   * set_z_line:  This function check and do Z.MESH card
   * which init mesh line in z direction.
   */
  int set_z_line(const Parser::Card &c);

  /**
   * @return the index of x direction skeleton line by x coordinate
   */
  int  find_skeleton_line_x(double x);

  /**
   * @return the index of y direction skeleton line by y coordinate
   */
  int  find_skeleton_line_y(double y);

  /**
   * @return the index of z direction skeleton line by z coordinate
   */
  int  find_skeleton_line_z(double z);

};




#endif // #define __mesh_generation_h__
