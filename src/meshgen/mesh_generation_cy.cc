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


// Local includes
#include "mesh_generation_cy.h"




/* ----------------------------------------------------------------------------
 * set_r_line:  This function check and do Z.MESH card
 * which init mesh line in z direction.
 */
int MeshGeneratorCylinder::set_rt_line ( const Parser::Card &c )
{
  // get parameter value from command structure
  double rmax = c.get_real("r.max",0.0,"r.back");
  double rmin = c.get_real("r.min",0.0,"r.front");
  double width = c.get_real("width",0.0,"r.width");
  int theta    = c.get_int("theta.division", theta_points.empty() ? 12 : theta_points.back());
  bool theta_refine = c.get_bool("theta.refine", false);
  if(theta_refine && !theta_points.empty())
    theta = theta_points.back()*2;

  double r,dr;

  // if no previous R.MESH card (point_array is empty), must insert the first point.
  if(!r_points.size())
  {
    r_points.push_back(rmin);
    theta_points.push_back(0);
  }
  else  // rmin is determined by previous R.MESH card.
    rmin = r_points.back();

  if ( width > 0 )
    rmax = rmin+width;
  else
    width = rmax-rmin;

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
    r_points.push_back(rmax);
    theta_points.push_back(theta);
    return 0;
  }

  if ( ratio==1.0 )
    dr = width/nspaces;
  else
    dr = width*(ratio-1)/(pow(ratio,nspaces)-1);

  // limit min_space
  while( dr < min_space && nspaces > 1)
  {
    nspaces--;
    if ( ratio==1.0 )
      dr = width/nspaces;
    else
      dr = width*(ratio-1)/(pow(ratio,nspaces)-1);
  }

  r = rmin;

  for(int i=1; i<=nspaces; i++)
  {
    // record min/max grid size
    if( dr < h_min) h_min = dr;
    if( dr > h_max) h_max = dr;

    // the r position of next point
    r += dr;

    // add the point to skeleton line
    r_points.push_back(r);
    theta_points.push_back(theta);

    dr*=ratio;
  }

  return 0;
}



/* ----------------------------------------------------------------------------
 * set_z_line:  This function check and do Z.MESH card
 * which init mesh line in z direction.
 */
int MeshGeneratorCylinder::set_z_line ( const Parser::Card &c )
{
  // get parameter value from command structure
  double zmax = c.get_real("z.max",0.0,"z.back");
  double zmin = c.get_real("z.min",0.0,"z.front");
  double width = c.get_real("width",0.0);

  double z,dz;

  // if no previous X.MESH card (point_array is empty), must insert the first point.
  if(!z_points.size())
  {
    z_points.push_back(zmin);
  }
  else  // xmin is determined by previous X.MESH card.
    zmin = z_points.back();

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
    z_points.push_back(zmax);
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
    z_points.push_back(z);

    dz*=ratio;
  }

  return 0;
}


double MeshGeneratorCylinder::find_r(double ra) const
{
  int ir=0;
  double dr=1e10;
  for(unsigned int i=0; i<r_points.size(); i++)
    if(std::abs(ra-r_points[i])<dr)
    {
      dr=std::abs(ra-r_points[i]);
      ir=i;
    }
  return r_points[ir];
}

double MeshGeneratorCylinder::find_theta(double rmin, double thetaa) const
{
  if(rmin == 0.0) return thetaa;

  int ir=0;
  double dr=1e10;
  for(unsigned int i=0; i<r_points.size(); i++)
    if(std::abs(rmin-r_points[i])<dr)
    {
      dr=std::abs(rmin-r_points[i]);
      ir=i;
    }

  double theta_div =  360.0/theta_points[ir];
  int loc = static_cast<int>(thetaa/theta_div+0.5);

  return theta_div*loc;
}


double MeshGeneratorCylinder::find_z(double za) const
{
  int iz=0;
  double dz=1e10;
  for(unsigned int i=0; i<z_points.size(); i++)
    if(std::abs(za-z_points[i])<dz)
    {
      dz=std::abs(za-z_points[i]);
      iz=i;
    }
  return z_points[iz];

}



std::vector<Point> MeshGeneratorCylinder::divide_circle(Point center, double r, int division) const
{
  double angle = 360.0 / division;
  const double pi = 3.14159265358979323846;
  std::vector<Point>  pts;
  for(int i=0; i<division; i++)
  {
    double x = center.x() +  r * cos(i*angle / 180.0 * pi);
    double y = center.y() +  r * sin(i*angle / 180.0 * pi);
    pts.push_back( Point(x, y, 0.0) );
  }

  pts.push_back( pts.front() );
  return pts;
}


