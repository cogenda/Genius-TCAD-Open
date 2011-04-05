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

//  $Id: mesh_generation.cc,v 1.7 2008/07/09 05:58:16 gdiso Exp $


// C++ includes
#include <cmath> // for std::sqrt


// Local includes
#include "mesh_generation.h"


/* ----------------------------------------------------------------------------
 * set_x_line:  This function check and do X.MESH card
 * which init mesh line in x direction.
 */
int MeshGenerator::set_x_line ( const Parser::Card &c )
{
  // get parameter value from command structure
  double xmax = c.get_real("x.max",0.0,"x.right");
  double xmin = c.get_real("x.min",0.0,"x.left");
  double width = c.get_real("width",0.0);

  double x,dx;

  // if no previous X.MESH card (point_array is empty), must insert the first point.
  if(!skeleton_line_x.point_array.size())
  {
      SkeletonPoint  p(xmin,0,0);
      skeleton_line_x.insert(p);
  }
  else  // xmin is determined by previous X.MESH card.
      xmin = skeleton_line_x.point_array[skeleton_line_x.point_array.size()-1].x;

  if ( width > 0 )
      xmax = xmin+width;
  else
      width = xmax-xmin;

  double ratio = c.get_real("ratio",1.0) ;
  int  nspaces = c.get_int ("n.spaces",1);
  double min_space = c.get_real("min.space", 1e-4);

  //if both h1 and h2 defined
  if(c.is_parameter_exist("h1") && c.is_parameter_exist("h2"))
  {
        double h1 = c.get_real("h1",0.1);
        double h2 = c.get_real("h2",h1);
        if(h1==h2)
        {
                ratio = 1.0;
                nspaces = int(width/h1+0.5);
        }
        else
        {
                ratio = (width-h1)/(width-h2);
                nspaces = (int)(log(h2/h1)/log(ratio)+ 1 +0.5);
        }
  }
  //  h1 and ratio
  else if(c.is_parameter_exist("h1"))
  {
        double h1 = c.get_real("h1",0.1);
        ratio = c.get_real("ratio",1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(width/h1+0.5);
        }
        else
        {
		nspaces = (int)(log(1-(1-ratio)*width/h1)/log(ratio)+0.5);
        }
  }
  //  h2 and ratio
  else if(c.is_parameter_exist("h2"))
  {
        double h2 = c.get_real("h2",0.1);
        ratio = 1.0/c.get_real("ratio",1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(width/h2+0.5);
        }
        else
        {
		double h1 = width*(1-ratio)+h2*ratio;
		nspaces = int(log(h2/h1)/log(ratio)+1.0+0.5);
        }
  }

  if ( nspaces==1 )
  {
      SkeletonPoint  p(xmax,0,0);
      skeleton_line_x.insert(p);
      return 0;
  }

  if ( ratio==1.0 )
      dx = width/nspaces;
  else
      dx = width*(ratio-1)/(pow(ratio,nspaces)-1);

  // limit min_space
  while( dx < min_space && nspaces > 1)
  {
      nspaces--;
      if ( ratio==1.0 )
        dx = width/nspaces;
      else
        dx = width*(ratio-1)/(pow(ratio,nspaces)-1);
  }

  x = xmin;

  for(int i=1; i<=nspaces; i++)
  {
      // record min/max grid size
      if( dx < h_min) h_min = dx;
      if( dx > h_max) h_max = dx;

      // the x position of next point
      x += dx;
      // add the point to skeleton line
      SkeletonPoint  p(x,0,0);
      skeleton_line_x.insert(p);

      dx*=ratio;
  }

  return 0;
}



/* ----------------------------------------------------------------------------
 * set_y_line:  This function check and do Y.MESH card
 * which init mesh line in y direction.
 */
int MeshGenerator::set_y_line ( const Parser::Card &c )
{
  // get parameter value from command structure
  double ymax = c.get_real("y.max",0.0,"y.bottom");
  double ymin = c.get_real("y.min",0.0,"y.top");
  double depth = c.get_real("depth",0.0);

  double y, dy;

  // if no previous Y.MESH card (point_array is empty), must insert the first point.
  if(!skeleton_line_y.point_array.size())
  {
      SkeletonPoint  p(0,ymin,0);
      skeleton_line_y.insert(p);
  }
  else  // ymin is determined by previous Y.MESH card.
      ymin = skeleton_line_y.point_array[skeleton_line_y.point_array.size()-1].y;

  if ( depth > 0 )
      ymax = ymin+depth;
  else
      depth = ymax-ymin;

  double ratio = c.get_real("ratio",1.0) ;
  int  nspaces = c.get_int ("n.spaces",1);
  double min_space = c.get_real("min.space", 1e-4);

  //if both h1 and h2 defined
  if(c.is_parameter_exist("h1") && c.is_parameter_exist("h2"))
  {
        double h1 = c.get_real("h1",0.1);
        double h2 = c.get_real("h2",h1);
        if(h1==h2)
        {
                ratio = 1.0;
                nspaces = int(depth/h1+0.5);
        }
        else
        {
                ratio = (depth-h1)/(depth-h2);
                nspaces = (int)(log(h2/h1)/log(ratio)+ 1 +0.5);
        }
  }
  //  h1 and ratio
  else if(c.is_parameter_exist("h1"))
  {
        double h1 = c.get_real("h1",0.1);
        ratio = c.get_real("ratio",1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(depth/h1+0.5);
        }
        else
        {
		nspaces = (int)(log(1-(1-ratio)*depth/h1)/log(ratio)+0.5);
        }
  }
  //  h2 and ratio
  else if(c.is_parameter_exist("h2"))
  {
        double h2 = c.get_real("h2",0.1);
        ratio = 1.0/c.get_real("ratio",1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(depth/h2+0.5);
        }
        else
        {
		double h1 = depth*(1-ratio)+h2*ratio;
		nspaces = int(log(h2/h1)/log(ratio)+1.0+0.5);
        }
  }

  if ( nspaces==1 )
  {
      SkeletonPoint  p(0,ymax,0);
      skeleton_line_y.insert(p);
      return 0;
  }

  if ( ratio==1.0 )
      dy = depth/nspaces;
  else
      dy = depth*(ratio-1)/(pow(ratio,nspaces)-1);

  // limit min_space
  while( dy < min_space && nspaces > 1)
  {
    nspaces--;
    if ( ratio==1.0 )
      dy = depth/nspaces;
    else
      dy = depth*(ratio-1)/(pow(ratio,nspaces)-1);
  }

  y = ymin;

  for(int i=1; i<=nspaces; i++)
  {
      // record min/max grid size
      if( dy < h_min) h_min = dy;
      if( dy > h_max) h_max = dy;

      // the y position of next point
      y += dy;

      // add the point to skeleton line
      SkeletonPoint  p(0,y,0);
      skeleton_line_y.insert(p);

      dy*=ratio;
  }

  return 0;
}



/* ----------------------------------------------------------------------------
 * set_z_line:  This function check and do Z.MESH card
 * which init mesh line in z direction.
 */
int MeshGenerator::set_z_line ( const Parser::Card &c )
{
  // get parameter value from command structure
  double zmax = c.get_real("z.max",0.0,"z.back");
  double zmin = c.get_real("z.min",0.0,"z.front");
  double width = c.get_real("width",0.0);

  double z,dz;

  // if no previous X.MESH card (point_array is empty), must insert the first point.
  if(!skeleton_line_z.point_array.size())
  {
      SkeletonPoint  p(0,0,zmin);
      skeleton_line_z.insert(p);
  }
  else  // xmin is determined by previous X.MESH card.
      zmin = skeleton_line_z.point_array[skeleton_line_z.point_array.size()-1].z;

  if ( width > 0 )
      zmax = zmin+width;
  else
      width = zmax-zmin;

  double ratio = c.get_real("ratio",1.0) ;
  int  nspaces = c.get_int ("n.spaces",1);
  double min_space = c.get_real("min.space", 1e-4);

  //if both h1 and h2 defined
  if(c.is_parameter_exist("h1") && c.is_parameter_exist("h2"))
  {
        double h1 = c.get_real("h1",0.1);
        double h2 = c.get_real("h2",h1);
        if(h1==h2)
        {
                ratio = 1.0;
                nspaces = int(width/h1+0.5);
        }
        else
        {
                ratio = (width-h1)/(width-h2);
                nspaces = (int)(log(h2/h1)/log(ratio)+ 1 +0.5);
        }
  }
  //  h1 and ratio
  else if(c.is_parameter_exist("h1"))
  {
        double h1 = c.get_real("h1",0.1);
        ratio = c.get_real("ratio",1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(width/h1+0.5);
        }
        else
        {
		nspaces = (int)(log(1-(1-ratio)*width/h1)/log(ratio)+0.5);
        }
  }
  //  h2 and ratio
  else if(c.is_parameter_exist("h2"))
  {
        double h2 = c.get_real("h2",0.1);
        ratio = 1.0/c.get_real("ratio",1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(width/h2+0.5);
        }
        else
        {
		double h1 = width*(1-ratio)+h2*ratio;
		nspaces = int(log(h2/h1)/log(ratio)+1.0+0.5);
        }
  }

  if ( nspaces==1 )
  {
      SkeletonPoint  p(0,0,zmax);
      skeleton_line_z.insert(p);
      return 0;
  }

  if ( ratio==1.0 )
      dz = width/nspaces;
  else
      dz = width*(ratio-1)/(pow(ratio,nspaces)-1);

  // limit min_space
  while( dz < min_space && nspaces > 1)
  {
    nspaces--;
    if ( ratio==1.0 )
      dz = width/nspaces;
    else
      dz = width*(ratio-1)/(pow(ratio,nspaces)-1);
  }

  z = zmin;

  for(int i=1; i<=nspaces; i++)
  {
      // record min/max grid size
      if( dz < h_min) h_min = dz;
      if( dz > h_max) h_max = dz;

      // the z position of next point
      z += dz;

      // add the point to skeleton line
      SkeletonPoint  p(0,0,z);
      skeleton_line_z.insert(p);

      dz*=ratio;
  }

  return 0;
}


/* ----------------------------------------------------------------------------
 * find_skeleton_line_x:  find the nearest x skeleton line index by x location
 */
int MeshGenerator::find_skeleton_line_x(double x)
{
  int ix=0;
  double dx=1e10;
  for(unsigned int i=0;i<IX;i++)
    if(fabs(x-skeleton_line_x.point_array[i].x)<dx)
    {
      dx=fabs(x-skeleton_line_x.point_array[i].x);
      ix=i;
    }
  return ix;
}


/* ----------------------------------------------------------------------------
 * find_skeleton_line_y:  find the nearest y skeleton line index by y location
 */
int MeshGenerator::find_skeleton_line_y(double y)
{
  int iy=0;
  double dy=1e10;
  for(unsigned int i=0;i<IY;i++)
    if(fabs(y-skeleton_line_y.point_array[i].y)<dy)
    {
      dy=fabs(y-skeleton_line_y.point_array[i].y);
      iy=i;
    }
  return iy;
}


/* ----------------------------------------------------------------------------
 * find_skeleton_line_z:  find the nearest z skeleton line index by z location
 */
int MeshGenerator::find_skeleton_line_z(double z)
{
  int iz=0;
  double dz=1e10;
  for(unsigned int i=0;i<IZ;i++)
    if(fabs(z-skeleton_line_z.point_array[i].z)<dz)
    {
      dz=fabs(z-skeleton_line_z.point_array[i].z);
      iz=i;
    }
  return iz;
}



