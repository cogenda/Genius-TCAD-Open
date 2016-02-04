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



#ifndef __mesh_generation_cy_h__
#define __mesh_generation_cy_h__



// C++ Includes   -----------------------------------
#include <vector>

// Local Includes -----------------------------------
#include "genius_common.h"

#include "enum_elem_type.h" // needed for ElemType enum
#include "point.h"
#include "parser.h"
#include "skeleton.h"
#include "mesh_generation_base.h"



class MeshGeneratorCylinder : public MeshGeneratorBase
{
public:

  /**
   * Constructor.
   */
  MeshGeneratorCylinder (MeshBase& mesh)
      : MeshGeneratorBase(mesh), h_min(1e30), h_max(0.0)
  {}


  /**
   * destructor.
   */
  virtual ~MeshGeneratorCylinder() {}


protected:

  /**
   * the min/max grid size
   */
  double h_min, h_max;

  /**
   * points in r direction
   */
  std::vector<double> r_points;

  /**
   * find exact r point neraest to ra  
   */
  double find_r(double ra) const;

  /**
   *  theta division
   */
  std::vector<int> theta_points;

  /**
   * find exact theta neraest to theta  
   */
  double find_theta(double rmin, double thetaa) const;

  /**
   *  points in z direction
   */
  std::vector<double> z_points;


  /**
   * find exact z point neraest to za
   */
  double find_z(double za) const;

  /**
   * set_z_line:  This function check and do R.MESH card
   * which init mesh line in r direction.
   */
  int set_rt_line(const Parser::Card &c);

  /**
   * set_z_line:  This function check and do Z.MESH card
   * which init mesh line in z direction.
   */
  int set_z_line(const Parser::Card &c);

  /**
   * aux function to get points in circle
   */
  std::vector<Point>  divide_circle(Point center, double r, int division) const;


};




#endif // #define __mesh_generation_cy_h__
