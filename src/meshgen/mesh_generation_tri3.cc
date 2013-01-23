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

//  $Id: mesh_generation_Tri3.cc,v 1.32 2008/07/09 05:58:16 gdiso Exp $


// C++ includes
#include <cmath> // for std::sqrt
#include <cstdlib>

// Local includes
#include "genius_env.h"
#include "log.h"
#include "unstructured_mesh.h"
#include "mesh_generation_tri3.h"

#include "mathfunc.h" //for Erfc


#include "point.h"
#include "boundary_info.h"
#include "physical_unit.h"

#include "perf_log.h"


void MeshGeneratorTri3::triangulateio_init()
{
#ifdef __triangle_h__
  // init Triangle io structure.
  in.pointlist = (double *) NULL;
  in.pointattributelist = (double *) NULL;
  in.pointmarkerlist = (int *) NULL;
  in.segmentlist = (int *) NULL;
  in.segmentmarkerlist = (int *) NULL;
  in.regionlist = (double *)NULL;

  out.pointlist = (double *) NULL;
  out.pointattributelist = (double *) NULL;
  out.pointmarkerlist = (int *) NULL;
  out.trianglelist = (int *) NULL;
  out.triangleattributelist = (double *) NULL;
  out.segmentlist = (int *) NULL;
  out.segmentmarkerlist = (int *) NULL;
#else


#endif
}


void MeshGeneratorTri3::triangulateio_finalize()
{
#ifdef __triangle_h__
  free(in.pointlist);
  free(in.pointmarkerlist);
  free(in.pointattributelist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);
  free(in.regionlist);

  free(out.pointlist);
  free(out.pointmarkerlist);
  free(out.pointattributelist);
  free(out.trianglelist);
  free(out.triangleattributelist);
  free(out.segmentlist);
  free(out.segmentmarkerlist);
#else
  in.clear();
  out.clear();
#endif
}



/* ----------------------------------------------------------------------------
 * build the rectangle Skeleton mesh in XY plane
 */
void MeshGeneratorTri3::build_rectangle_mesh()
{

  // the point number in XY direction
  IX = skeleton_line_x.point_array.size();
  IY = skeleton_line_y.point_array.size();
  IZ = 0;

  point_num = IX*IY;

  // set up the 3d array for gird points
  // however, we only need xy plane information
  point_array3d = new SkeletonPoint **[1];

  point_array3d[0]= new SkeletonPoint *[IY];
  for(unsigned int j=0;j<IY;j++)
    point_array3d[0][j]= new SkeletonPoint [IX];

  for(unsigned int j=0;j<IY;j++)
    for(unsigned int i=0;i<IX;i++)
    {
      double x = skeleton_line_x.point_array[i].x;
      double y = skeleton_line_y.point_array[j].y;
      double z = 0.0;
      point_array3d[0][j][i].set_location(x,y,z);
    }

  // set the bound box
  xmin  =  point_array3d[0][0][0].x;
  xmax  =  point_array3d[0][0][IX-1].x;
  ymax  =  point_array3d[0][0][0].y;
  ymin  =  point_array3d[0][IY-1][0].y;
  zmin  =  0.0;
  zmax  =  0.0;
}



/* ----------------------------------------------------------------------------
 * eliminate unnecessary points in x direction
 */
int MeshGeneratorTri3::x_eliminate(int ixmin,int ixmax, int iymin, int iymax)
{
  for( int i=ixmin; i<=ixmax; i++ )
  {
    bool eliminate_flag = true;

    for( int j=std::max(1,iymin+1); j<std::min(iymax+1,static_cast<int>(IY-1)); j++ )
      if( !point_array3d[0][j][i].eliminated )
      {
        if( eliminate_flag == true )
        {
          point_array3d[0][j][i].eliminated = true;
          eliminate_flag = false;
        }
        else
          eliminate_flag = true;
      }
  }
  return 0;
}



/* ----------------------------------------------------------------------------
 * eliminate unnecessary points in y direction
 */
int MeshGeneratorTri3::y_eliminate(int iymin,int iymax, int ixmin, int ixmax)
{
  for( int j=iymin; j<=iymax; j++ )
  {
    bool eliminate_flag = true;
    for ( int i=std::max(1,ixmin+1); i<std::min(ixmax+1,static_cast<int>(IX-1)); i++)
      if( !point_array3d[0][j][i].eliminated )
      {
        if( eliminate_flag == true )
        {
          point_array3d[0][j][i].eliminated = true;
          eliminate_flag = false;
        }
        else
          eliminate_flag = true;
      }
  }
  return 0;
}



/* ----------------------------------------------------------------------------
 * renumber the nodes
 */
void MeshGeneratorTri3::make_node_index()
{
  point_num=0;

  for(unsigned int i=0; i<IX; i++)
    for(unsigned int j=0; j<IY; j++)
    {
      if(!point_array3d[0][j][i].eliminated)
        point_array3d[0][j][i].index = point_num++;
      else
        point_array3d[0][j][i].index = invalid_uint;
    }
}



/* ----------------------------------------------------------------------------
 * this function do 2D spread on xy plane. copied from Pisces code
 */
int MeshGeneratorTri3::make_spread(  const std::string &location, double width,
                                     unsigned int upperline,unsigned int lowerline,
                                     double yuploc,double yloloc,double encroach,double grading)
{
  double xloc,yupold,yloold,ybot;
  if(location == "left")
    xloc = point_array3d[0][0][0].x + width;
  else
    xloc = point_array3d[0][0][IX-1].x - width;

  double  erfarg,erfar2,erfval,erfvl2;

  //proc colume by colume
  for(unsigned int i=0;i<IX;i++)
  {
    double xco = point_array3d[0][0][i].x;

    //evaluate error function for this coord.
    erfarg=(xco-xloc)*encroach;
    if (location== "right")
      erfar2=1.5*(erfarg+0.6);
    else
      erfar2=1.5*(erfarg-0.6);
    if (location== "left")
    {
#ifdef WINDOWS
      erfval=Erfc(erfarg);
      erfvl2=Erfc(erfar2);
#else
      erfval=erfc(erfarg);
      erfvl2=erfc(erfar2);
#endif

    }
    else
    {
#ifdef WINDOWS
      erfval=Erfc(-erfarg);
      erfvl2=Erfc(-erfar2);
#else
      erfval=erfc(-erfarg);
      erfvl2=erfc(-erfar2);
#endif

    }
    erfval=erfval*0.5;
    erfvl2=erfvl2*0.5;

    //
    // compute new node locations on this column

    //  get upper, lower, bottom current loc.
    yupold = point_array3d[0][upperline][i].y;
    yloold = point_array3d[0][lowerline][i].y;
    ybot = point_array3d[0][IY-1][i].y;

    // compute upward shift and downward
    double deltup=erfval*(yupold-yuploc);
    double deltlo=erfval*(yloloc-yloold);

    //  compute new locations
    double yupnew=yupold-deltup;
    double ylonew=yloold+deltlo;

    // compute old and new spreads of middle, bottom
    double spmdol=yloold-yupold;
    double spmdnw=spmdol+deltup+deltlo;
    double spbtol=ybot-yloold;
    double spbtnw=spbtol-deltlo;

    // grading ratio
    double y,dy;
    if(grading!=1.0)
    {
      y  =  yupnew;
      dy = (ylonew-yupnew)*(grading-1)/(pow(grading,static_cast<double>(lowerline-upperline))-1);
    }

    // scan y nodes and move as necessary
    for(unsigned int j=0; j<IY; j++)
    {
      // get node number and y-coord.
      double yco=point_array3d[0][j][i].y;

      // consider top, middle, and bottom regions

      //  top region (j<upperline) shift all by deltup
      if(j<=upperline)
      {
        yco=yco-deltup;
      }
      // bottom region, spread proportionally
      else if (j>=lowerline)
      {
        double rat=(yco-yloold)/spbtol;
        yco=ylonew+rat*spbtnw;
      }
      // middle region spread proportionally unless new grading
      // is requested.
      else
      {
        if(grading==1.0)
        {
          double rat=(yco-yupold)/spmdol;
          yco=yupnew+rat*spmdnw;
        }
        else // new grading requested
        {
          double ycordg = y+dy;
          y += dy;
          dy*= grading;
          // vary from new grading to proportional grading
          // based on erfvl2 (2/3 the spread of erfval and
          // centered at 30% point of erfval instead of
          // 50% point.
          yco=yco+erfvl2*(ycordg-yco);
        }
      }
      point_array3d[0][j][i].y = yco;
    }
  }
  return 0;
}



/* ----------------------------------------------------------------------------
 * record the boundary segment of the region
 */
int  MeshGeneratorTri3::make_region_segment()
{
  // set inner point for each region;
  for(unsigned int i=0;i<IX-1;i++)
    for(unsigned int j=0;j<IY-1;j++)
    {
      for(int k=region_array1d.size()-1;k>=0;k--)
        if(//region_array1d[k].shape==Rectangle &&
          i>=region_array1d[k].ixmin && i+1<=region_array1d[k].ixmax &&
          j>=region_array1d[k].iymin && j+1<=region_array1d[k].iymax)
        {
          region_array1d[k].px.push_back((point_array3d[0][j][i].x+point_array3d[0][j][i+1].x)/2.0);
          region_array1d[k].py.push_back((point_array3d[0][j][i].y+point_array3d[0][j+1][i].y)/2.0);
          break;
        }
    }

  // set segment bounding
  std::map<SkeletonEdge, int, lt_edge>::iterator pt;

  for(size_t r=0;r<region_array1d.size();r++)
  {
    //if(region_array1d[r].shape! = Rectangle) continue;

    SkeletonFace  face;
    face.face_label=region_array1d[r].label+"_Neumann";
    face.face_mark=r+1;
    face.user_define = false;
    face_array1d.push_back(face);

    SkeletonEdge edge;
    edge.IX = IX;
    edge.IY = IY;

    //process bottom line
    for( unsigned int  i=region_array1d[r].ixmin; i<region_array1d[r].ixmax; )
    {
      edge.p1[0] = i++;
      edge.p1[1] = region_array1d[r].iymax;
      while( point_array3d[0][region_array1d[r].iymax][i].eliminated ) i++;
      edge.p2[0] = i;
      edge.p2[1] = region_array1d[r].iymax;
      edge_table[edge] =  (r+1);
    }
    //process top line
    for( unsigned int  i=region_array1d[r].ixmin; i<region_array1d[r].ixmax; )
    {
      edge.p1[0] = i++;
      edge.p1[1] = region_array1d[r].iymin;
      while( point_array3d[0][region_array1d[r].iymin][i].eliminated ) i++;
      edge.p2[0] = i;
      edge.p2[1] = region_array1d[r].iymin;
      edge_table[edge] =  (r+1);
    }
    //process left line
    for( unsigned int  i=region_array1d[r].iymin; i<region_array1d[r].iymax; )
    {
      edge.p1[0] = region_array1d[r].ixmin;
      edge.p1[1] = i++;
      while( point_array3d[0][i][region_array1d[r].ixmin].eliminated ) i++;
      edge.p2[0] = region_array1d[r].ixmin;
      edge.p2[1] = i;
      // is it the global left line
      edge_table[edge] =  (r+1);
    }
    //process right line
    for( unsigned int  i=region_array1d[r].iymin; i<region_array1d[r].iymax; )
    {
      edge.p1[0] = region_array1d[r].ixmax;
      edge.p1[1] = i++;
      while( point_array3d[0][i][region_array1d[r].ixmax].eliminated ) i++;
      edge.p2[0] = region_array1d[r].ixmax;
      edge.p2[1] = i;
      edge_table[edge] =  (r+1);
    }
  }
  return 0;
}



/* ----------------------------------------------------------------------------
 * make_face :  set face (segments in 2D) with user defined label
 */
int  MeshGeneratorTri3::make_face(int ixmin,int ixmax,int iymin,int iymax, const std::string &label)
{
  int face_mark = face_array1d.size() + 1;

  SkeletonFace face;
  face.face_label=label;
  face.face_mark=face_mark;
  face.user_define=true;
  face.ixmin = ixmin;
  face.ixmax = ixmax;
  face.iymin = iymin;
  face.iymax = iymax;
  face.izmin = 0;
  face.izmax = 0;
  face_array1d.push_back(face);

  SkeletonEdge edge;
  edge.IX = IX;
  edge.IY = IY;

  if(iymin==iymax) //top or bottom face
  {
    // set the projection segment of this face on xy plane
    for(int i=ixmin;i<ixmax;)
    {
      edge.p1[0]=i++;
      edge.p1[1]=iymax;
      while(point_array3d[0][iymax][i].eliminated) i++;
      edge.p2[0]=i;
      edge.p2[1]=iymax;
      edge_table[edge]=face_mark;
    }
  }
  else if(ixmin==ixmax) //left or right face
  {
    // set the projection segment of this face on xy plane
    for(int i=iymin;i<iymax;)
    {
      edge.p1[0]=ixmin;
      edge.p1[1]=i++;
      while(point_array3d[0][i][ixmin].eliminated) i++;
      edge.p2[0]=ixmin;
      edge.p2[1]=i;
      edge_table[edge]=face_mark;
    }
  }

  return 0;
}





/* ----------------------------------------------------------------------------
 * set_eliminate:  This function check and do ELIMINATE card
 * eliminate unnecessary lines in x or y direction.
 */
int MeshGeneratorTri3::set_eliminate(const Parser::Card &c)
{

  // get parameter value from command structure
  unsigned int ixmin = c.get_int("ix.min",0,"ix.left");
  unsigned int ixmax = c.get_int("ix.max",IX-1,"ix.right");
  unsigned int iymin = c.get_int("iy.min",0,"iy.top");
  unsigned int iymax = c.get_int("iy.max",IY-1,"iy.bottom");

  std::string dir = c.get_string("direction","");

  if(c.is_parameter_exist("x.max") || c.is_parameter_exist("x.right"))
  {
    double _xmax = c.get_real("x.max",xmax,"x.right");
    ixmax = find_skeleton_line_x(_xmax);
  }
  if(c.is_parameter_exist("x.min") || c.is_parameter_exist("x.left"))
  {
    double _xmin = c.get_real("x.min",xmin,"x.left");
    ixmin = find_skeleton_line_x(_xmin);
  }
  if(c.is_parameter_exist("y.max") || c.is_parameter_exist("y.bottom"))
  {
    double _ymax = c.get_real("y.max",ymax,"y.bottom");
    iymax = find_skeleton_line_y(_ymax);
  }
  if(c.is_parameter_exist("y.min") || c.is_parameter_exist("y.top"))
  {
    double _ymin = c.get_real("y.min",ymin,"y.top");
    iymin = find_skeleton_line_y(_ymin);
  }

  if( dir == "ynorm" )
    y_eliminate(iymin,iymax,ixmin,ixmax);
  else if( dir == "xnorm" )
    x_eliminate(ixmin,ixmax,iymin,iymax);
  else
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()
                           << " ELIMINATE: Tri3 mesh generator not support eliminate along Z direction!\n";
    RECORD();
    return 1;
  }

  return 0;
}




/* ----------------------------------------------------------------------------
 * set_spread:  This function check and do SPREAD card
 * I copyed some code from PISCES. please forgive me...
 */
int MeshGeneratorTri3::set_spread(const Parser::Card &c)
{

  // get parameter value from command structure

  // get x location of distorted region.
  int xpoint;

  if( ! c.is_parameter_exist("location") )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " SPREAD: you must give location of spread region!\n";
    RECORD();
    return 1;
  }

  std::string location = c.get_string("location","");

  if( location == "left" )
  { xpoint=0; }
  else if( location == "right" )
  { xpoint=IX-1; }
  else
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " SPREAD: you must give correct location of spread region!\n";
    RECORD();
    return 1;
  }
  double width = c.get_real("width",0.0);

  // the upper and lower grid line of distorted region
  int  upperline = c.get_int("upper",0);
  int  lowerline = c.get_int("lower",0);

  // current thickness
  double cthick = point_array3d[0][upperline][xpoint].y - point_array3d[0][lowerline][xpoint].y;

  // get the new location of upper and lower y grid line
  double yuploc,yloloc,thick;
  if(c.is_parameter_exist("y.lower"))
  {
    yuploc = point_array3d[0][upperline][xpoint].y;
    yloloc = c.get_real("y.lower",0.0);
    thick = yloloc - yuploc ;
  }
  else if(c.is_parameter_exist("thickness"))
  {
    thick = c.get_real("thickness",1.0);
    double vol_rat = c.get_real("vol.ratio",0.44);
    double dthick = vol_rat*(thick-cthick);
    double uthick = thick-dthick-cthick;
    yuploc = point_array3d[0][upperline][xpoint].y - uthick;
    yloloc = point_array3d[0][lowerline][xpoint].y + dthick;
  }

  // this parameter give the transition between distorted and undistorted grid. in x direction.
  // from PISCES:
  //   scale so that 80% of thickness change occurs in a
  //   distance equal to the thickness change (for an
  //   encroachment factor of 1)
  double encroach =  c.get_real("encroach",1.0);
  if(encroach<0.1) encroach = 0.1;

  const double erfc80 = 1.812386;
  encroach=fabs(erfc80/((thick-cthick)*encroach));

  // get grading parameter(s)
  double grading =  c.get_real("grading",1.0);

  make_spread(location,width,upperline,lowerline,yuploc,yloloc,encroach,grading);

  return 0;
}




/* ----------------------------------------------------------------------------
 * set_region:  This function check and do ellipse REGION card
 * it specify the material of a region
 */
int MeshGeneratorTri3::set_region_ellipse(const Parser::Card &c)
{

  // get parameter value from command structure
  SkeletonRegion2D region;

  region.shape = 1;

  region.division = c.get_int("division", 12);
  region.centrex = c.get_real("centrex", 0.0);
  region.centrey = c.get_real("centrey", 0.0);
  region.major_radii  = c.get_real("majorradii", 1.0);
  region.minor_radii  = c.get_real("minorradii", region.major_radii);
  region.theta        = c.get_real("theta", 0.0)/360*(2*3.14159265359);

  if(region.major_radii<=0 || region.minor_radii<=0)
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " REGION: Circle radius should be positive!\n";
    RECORD();
    return 1;
  }

  region.label    = c.get_string("label", "");
  region.material = c.get_string("material", "");

  // test if region label is unique
  for(unsigned int n=0; n<region_array1d.size(); ++n)
  {
    if(region_array1d[n].label == region.label)
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " REGION: Region label " << region.label << " has been defined before!\n";
      RECORD();
      return 1;
    }
  }

  region_array1d.push_back(region);

  return 0;
}




/* ----------------------------------------------------------------------------
 * set_region:  This function check and do rectangle REGION card
 * it specify the material of a region
 */
int MeshGeneratorTri3::set_region_rectangle(const Parser::Card &c)
{


  // get parameter value from command structure
  SkeletonRegion2D region;

  // some problem: ixxx is unsigned int, but get_int may return nagetive vlaue...
  // however, for this case, the unsigned int value will be very large, much large than
  // IX/IY/IZ
  region.ixmin = c.get_int("ix.min",0,"ix.left");
  region.ixmax = c.get_int("ix.max",IX-1,"ix.right");
  region.iymin = c.get_int("iy.min",0,"iy.top");
  region.iymax = c.get_int("iy.max",IY-1,"iy.bottom");
  region.label = c.get_string("label","") ;

  if(!c.is_parameter_exist("material"))
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " REGION: Material should be given for each region!\n";
    RECORD();
  }

  region.material= c.get_string("material","") ;

  region.shape = 0;

  if(c.is_parameter_exist("x.max") || c.is_parameter_exist("x.right"))
  {
    double _xmax = c.get_real("x.max",xmax,"x.right");
    region.ixmax = find_skeleton_line_x(_xmax);
  }
  if(c.is_parameter_exist("x.min") || c.is_parameter_exist("x.left"))
  {
    double _xmin = c.get_real("x.min",xmin,"x.left");
    region.ixmin = find_skeleton_line_x(_xmin);
  }
  if(c.is_parameter_exist("y.max") || c.is_parameter_exist("y.bottom"))
  {
    double _ymax = c.get_real("y.max",ymax,"y.bottom");
    region.iymax = find_skeleton_line_y(_ymax);
  }
  if(c.is_parameter_exist("y.min") || c.is_parameter_exist("y.top"))
  {
    double _ymin = c.get_real("y.min",ymin,"y.top");
    region.iymin = find_skeleton_line_y(_ymin);
  }


  if( static_cast<int>(region.ixmin) < 0  ||  region.ixmax > (IX-1) || region.ixmin > region.ixmax )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " REGION: Can't locate left/right boundary of region!\n";
    RECORD();
    return 1;
  }
  if( static_cast<int>(region.iymin) < 0  ||  region.iymax > (IY-1) || region.iymin > region.iymax )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " REGION: Can't locate top/bottom boundary of region!\n";
    RECORD();
    return 1;
  }

  // test if region label is unique
  for(unsigned int n=0; n<region_array1d.size(); ++n)
  {
    if(region_array1d[n].label == region.label)
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " REGION: Region label " << region.label << " has been defined before!\n";
      RECORD();
      return 1;
    }
  }

  region_array1d.push_back(region);

  return 0;


}




/* ----------------------------------------------------------------------------
 * set_region:  This function check and do REGION card
 * specify the material of a region
 */
int MeshGeneratorTri3::set_region(const Parser::Card &c)
{

  std::string shape = c.get_string("shape","rectangle");

  if( shape == "rectangle" ) return set_region_rectangle(c);

  if( shape == "ellipse" )   return set_region_ellipse(c);

  return 0; //prevent warning
}








/* ----------------------------------------------------------------------------
 * set_face:  This function check and do FACE card
 * specify the label and location of a face.
 * these segments may be specified as electrode later by BOUNDARY card.
 */
int MeshGeneratorTri3::set_face(const Parser::Card &c)
{
  unsigned int ixmin,ixmax,iymin,iymax;
  int ix,iy;

  // get parameter value from command structure
  ix    = c.get_int("ix",0);
  ixmin = c.get_int("ix.min",0,"ix.left");
  ixmax = c.get_int("ix.max",IX-1,"ix.right");
  iy    = c.get_int("iy",0);
  iymin = c.get_int("iy.min",0,"iy.top");
  iymax = c.get_int("iy.max",IY-1,"iy.bottom");

  std::string  label = c.get_string("label","");
  if(c.is_parameter_exist("x"))
  {
    double x = c.get_real("x",0.0);
    ix = find_skeleton_line_x(x);
  }
  if(c.is_parameter_exist("x.max") || c.is_parameter_exist("x.right"))
  {
    double _xmax = c.get_real("x.max",xmax,"x.right");
    ixmax = find_skeleton_line_x(_xmax);
  }
  if(c.is_parameter_exist("x.min") || c.is_parameter_exist("x.left"))
  {
    double _xmin = c.get_real("x.min",xmin,"x.left");
    ixmin = find_skeleton_line_x(_xmin);
  }
  if(c.is_parameter_exist("y"))
  {
    double y = c.get_real("y",0.0);
    iy = find_skeleton_line_y(y);
  }
  if(c.is_parameter_exist("y.max") || c.is_parameter_exist("y.bottom"))
  {
    double _ymax = c.get_real("y.max",ymax,"y.bottom");
    iymax = find_skeleton_line_y(_ymax);
  }
  if(c.is_parameter_exist("y.min") || c.is_parameter_exist("y.top"))
  {
    double _ymin = c.get_real("y.min",ymin,"y.top");
    iymin = find_skeleton_line_y(_ymin);
  }

  if(c.is_parameter_exist("location"))
  {
    if(c.is_enum_value("location","top"))
      iymin=iymax=0;
    else if(c.is_enum_value("location","bottom"))
      iymin=iymax=IY-1;
    else if(c.is_enum_value("location","left"))
      ixmin=ixmax=0;
    else if(c.is_enum_value("location","right"))
      ixmin=ixmax=IX-1;
    else if(c.is_enum_value("location","back") || c.is_enum_value("location","front"))
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " FACE: for Tet3 mesh, face can not locate on front/back!\n";
      RECORD();
      return 1;
    }
  }

  if(c.is_parameter_exist("direction"))
  {
    if(c.is_enum_value("direction","ynorm"))
      iymin=iymax=iy;
    else if(c.is_enum_value("direction","xnorm"))
      ixmin=ixmax=ix;
    else if(c.is_enum_value("direction","znorm"))
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " FACE: for Tet3 mesh, face norm can't along Z direction!\n";
      RECORD();
      return 1;
    }
  }

  if(ixmin!=ixmax && iymin!=iymax)
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " FACE: face must along grid lines!\n";
    RECORD();
    return 1;
  }

  if( ixmax>(IX-1) || ixmin > ixmax )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " FACE: face xline exceed range!\n";
    RECORD();
    return 1;
  }

  if( iymax>(IY-1) || iymin > iymax )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " FACE: face yline exceed range!\n";
    RECORD();
    return 1;
  }

  make_face(ixmin,ixmax,iymin,iymax,label);

  return 0;

}




/* ----------------------------------------------------------------------------
 * triangle_mesh:  call Triangle code (developed by Jonathan Richard Shewchuk)
 * to generate triangle mesh on xy plane
 */
int MeshGeneratorTri3::triangle_mesh()
{
  //set point
  in.numberofpoints = point_num;
  in.numberofpointattributes = 0;
  in.pointattributelist = (double *)NULL;
#ifdef __triangle_h__
  in.pointlist = (double *) calloc(in.numberofpoints*2, sizeof(double));
  in.pointmarkerlist = (int *) calloc(in.numberofpoints, sizeof(int));
#else
  in.pointlist =  new double[in.numberofpoints*2];
  in.pointmarkerlist = new int[in.numberofpoints];
#endif
  double *ppointlist=in.pointlist;
  int *ppointmarkerlist=in.pointmarkerlist;
  //the points belongs to rectangle region
  for(unsigned int i=0;i<IX;i++)
    for(unsigned int j=0;j<IY;j++)
      if(!point_array3d[0][j][i].eliminated)
      {
        *ppointlist++ = point_array3d[0][j][i].x;
        *ppointlist++ = point_array3d[0][j][i].y;
        *ppointmarkerlist++ = 0;
      }

  //do necessarily prepare for call triangulate
  in.numberoftriangles = 0;
  in.numberofcorners = 3;
  in.numberoftriangleattributes = 0;
  in.trianglelist =  (int *) NULL;
  in.trianglearealist = (double *) NULL;
  in.triangleattributelist = NULL;

  // set segment information
  in.numberofsegments = edge_table.size();
#ifdef __triangle_h__
  in.segmentlist =  (int *) calloc(in.numberofsegments*2,  sizeof(int));
  in.segmentmarkerlist = (int *) calloc(in.numberofsegments, sizeof(int));
#else
  in.segmentlist =  new int[in.numberofsegments*2];
  in.segmentmarkerlist = new int[in.numberofsegments];
#endif
  int *psegmentlist =  in.segmentlist;
  int *psegmentmarkerlist = in.segmentmarkerlist;
  std::map<SkeletonEdge, int, lt_edge>::iterator pt = edge_table.begin();
  for(size_t i=0;i<edge_table.size();i++)
  {
    *psegmentlist++ = point_array3d[0][pt->first.p1[1]][pt->first.p1[0]].index;
    *psegmentlist++ = point_array3d[0][pt->first.p2[1]][pt->first.p2[0]].index;
    *psegmentmarkerlist++ = pt->second;
    ++pt;
  }

  //set region information
  in.numberofholes = 0;
  in.numberofregions = region_array1d.size();
#ifdef __triangle_h__
  in.regionlist = (double *) calloc(in.numberofregions*4, sizeof(double));
#else
  in.regionlist = new double[in.numberofregions*4];
#endif
  double *pregionlist =  in.regionlist;
  for(int i=0;i<in.numberofregions;i++)
  {
    *pregionlist++ = region_array1d[i].px[0];
    *pregionlist++ = region_array1d[i].py[0];
    *pregionlist++ = double(i);
    *pregionlist++ = 0;
  }

  // call Triangle here
#ifdef __triangle_h__
  triangulate(const_cast<char *>(tri_cmd.c_str()), &in, &out, (struct triangulateio *) NULL);
#else
  ctri_triangulate(const_cast<char *>(tri_cmd.c_str()), &in, &out);
#endif
  // set boundary mark to output edge
  out_edge_table.clear();
  for(int i=0; i<out.numberofsegments; i++)
  {
    OutEdge edge;
    edge.pointno = out.numberofpoints;
    if(out.segmentlist[2*i+0] < out.segmentlist[2*i+1])
    {
      edge.p1 = out.segmentlist[2*i+0];
      edge.p2 = out.segmentlist[2*i+1];
    }
    else
    {
      edge.p1 = out.segmentlist[2*i+1];
      edge.p2 = out.segmentlist[2*i+0];
    }
    edge.mark = out.segmentmarkerlist[i];
    out_edge_table[edge] = edge.mark;

  }

  return 0;
}



/* ----------------------------------------------------------------------------
 * get the bc index of a Tet3 cell's side by cell's region_id,
 * and segment bc (get by p1 and p2)
 */
int  MeshGeneratorTri3::get_bc_id(const Elem *, int p1, int p2)
{
  int bc = 0;

  OutEdge edge;
  edge.pointno = out.numberofpoints;
  edge.p1 = p1 < p2 ? p1:p2;
  edge.p2 = p1 < p2 ? p2:p1;

  if ( out_edge_table.find(edge) != out_edge_table.end() )
  {
    bc = (*out_edge_table.find(edge)).second;
  }

  return bc;

}


/* ----------------------------------------------------------------------------
 * over write the do_mesh() virtual function in meshgen base class
 */
int MeshGeneratorTri3::do_mesh()
{

  START_LOG("do_mesh()", "MeshGenerator");

  MESSAGE<<"Generating Tri3 mesh by MeshGeneratorTri3...\n"; RECORD();

  dscale = PhysicalUnit::um;


  // build x-y grid skeleton
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if(c.key() == "MESH")      // It's a MESH card
      tri_cmd = c.get_string("triangle", "pzADq30Q");

    if(c.key() == "X.MESH")   // It's a X.MESH card
      if(set_x_line(c)) return 1;

    if(c.key() == "Y.MESH")   // It's a Y.MESH card
      if(set_y_line(c)) return 1;

    if(c.key() == "Z.MESH")   // It's a Z.MESH card, we do not like it.
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " Z.MESH: Tet3 mesh generator not support this statement!\n";
      RECORD();
      return 1;
    }
  }

  //after the initializtion of x and y mesh lines, we can build rectangle mesh now.
  build_rectangle_mesh();

  //use ELIMINATE and SPREAD card to do necessary modification to the rectangle mesh.
  //after that, use REGION card to set material region.
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if( c.key() == "ELIMINATE" )   // do eliminate
      if( set_eliminate(c) ) return 1;

    if( c.key() == "SPREAD" )      // do spread
      if( set_spread(c) ) return 1;

    if( c.key() == "SPREAD3D" )    //  spread3d?
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " SPREAD3D: Tet3 mesh generator not support this statement!\n";
      RECORD();
      return 1;
    }

    if( c.key() == "REGION" )      // set region here.
      if( set_region(c) ) return 1;
  }

  // pick up the boundary of xy plane
  // Thus, the boundary will be reserved during triangle mesh.
  make_region_segment();

  make_node_index();

  // set labeled faces(most of them are electrode) defined by user here.
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if(c.key() == "FACE")      // set face here.
      if( set_face(c) ) return 1;
  }

  triangulateio_init();

  // call triangle to generate 2D mesh on xy plane
  triangle_mesh();


  // fill the _mesh structure
  _mesh.magic_num() = this->magic_num();

  // Build the nodes. we do dimension scale here
  _mesh.reserve_nodes( out.numberofpoints );
  _mesh.reserve_elem( out.numberoftriangles );

  for(int i = 0; i < out.numberofpoints; i++)
  {
    _mesh.add_point(Point(out.pointlist[2*i+0]*dscale,
                          out.pointlist[2*i+1]*dscale,
                          0.0));
  }

  // Build the elements
  for(int i = 0; i < out.numberoftriangles; i++)
  {
    Elem* elem = _mesh.add_elem(Elem::build(TRI3).release());
    genius_assert(elem);

    elem->set_node(0) = _mesh.node_ptr( out.trianglelist[3*i+0] );
    elem->set_node(1) = _mesh.node_ptr( out.trianglelist[3*i+1] );
    elem->set_node(2) = _mesh.node_ptr( out.trianglelist[3*i+2] );
    elem->subdomain_id() = static_cast<int>(out.triangleattributelist[i]+0.5);

    int bc_index;

    if( (bc_index = get_bc_id(elem, out.trianglelist[3*i+0], out.trianglelist[3*i+1])) > 0 )
      _mesh.boundary_info->add_side(elem, 0, bc_index);

    if( (bc_index = get_bc_id(elem, out.trianglelist[3*i+1], out.trianglelist[3*i+2])) > 0 )
      _mesh.boundary_info->add_side(elem, 1, bc_index);

    if( (bc_index = get_bc_id(elem, out.trianglelist[3*i+2], out.trianglelist[3*i+0])) > 0 )
      _mesh.boundary_info->add_side(elem, 2, bc_index);

  }

  // write region label and material info to _mesh

  _mesh.set_n_subdomains() = region_array1d.size();

  for(size_t r=0; r < region_array1d.size(); r++)
  {
    _mesh.set_subdomain_label(r, region_array1d[r].label );
    _mesh.set_subdomain_material(r, region_array1d[r].material);
  }

  // write boundary label into mesh.boundary_info structure
  for(size_t f=0; f < face_array1d.size(); f++)
  {
    _mesh.boundary_info->set_label_to_id( face_array1d[f].face_mark, face_array1d[f].face_label, face_array1d[f].user_define );
  }


  //however, the interface information should be set here

  std::vector<unsigned int>       elems;
  std::vector<unsigned short int> sides;
  std::vector<short int>          bds;
  std::map<std::string, short int> bd_map;
  typedef std::map<std::string, short int>::iterator Bd_It;

  // get all the boundary element
  _mesh.boundary_info->build_side_list (elems, sides, bds);

  //build neighbor information for mesh. then elem->neighbor() is functional
  _mesh.boundary_info->find_neighbors();


  for (size_t nbd=0; nbd<elems.size(); nbd++ )
  {
    // get the element which has boundary/interface side
    const Elem* elem = _mesh.elem(elems[nbd]);
    genius_assert(elem->on_boundary() || elem->on_interface());

    //is it an interface side && not a label face
    if( elem->neighbor(sides[nbd])!=NULL && bds[nbd]<=static_cast<short int>(region_array1d.size()) )
    {
      // the element and its neighbor should in diffetent subdomain
      unsigned int sbd_id1 = elem->subdomain_id();
      unsigned int sbd_id2 = elem->neighbor(sides[nbd])->subdomain_id();

      // delete the overkilled boundary side
      if (sbd_id1 == sbd_id2)
      {
        _mesh.boundary_info->remove(elem, sides[nbd]);
        continue;
      }

      // the side should be an interface side
      genius_assert(elem->on_interface());
      genius_assert(elem->neighbor(sides[nbd])->on_interface());

      //remove the pair-element from boundary
      _mesh.boundary_info->remove(elem, sides[nbd]);
      _mesh.boundary_info->remove(elem->neighbor(sides[nbd]),
                                  elem->neighbor(sides[nbd])->which_neighbor_am_i(elem));

      // build the label for the interface, which has the form of RegionLabel1_to_RegionLabel2,
      // the two region is alpha ordered.
      std::string bd_label;
      if( region_array1d[sbd_id1].label < region_array1d[sbd_id2].label)
        bd_label = region_array1d[sbd_id1].label + "_to_" + region_array1d[sbd_id2].label;
      else
        bd_label = region_array1d[sbd_id2].label + "_to_" + region_array1d[sbd_id1].label;

      short int bd_index;

      // if the label already exist
      if( bd_map.find(bd_label) != bd_map.end() )
        bd_index = (*bd_map.find(bd_label)).second;
      else
      {
        //else, increase bd_index, insert it into bd_map
        bd_index = face_array1d.size() + bd_map.size() + 1;
        bd_map.insert(std::pair<const std::string, short int>(bd_label,bd_index));
      }

      // add pair-element to boundary with new bd_index
      _mesh.boundary_info->add_side(elem, sides[nbd], bd_index);
      _mesh.boundary_info->add_side(elem->neighbor(sides[nbd]),
                                    elem->neighbor(sides[nbd])->which_neighbor_am_i(elem),
                                    bd_index);
    }
  }

  // after the previous remove and insert operation, the number of boundary id maybe changed. renumber it here.
  _mesh.boundary_info->rebuild_ids();

  //write down new labels
  Bd_It bd_it = bd_map.begin();
  for(; bd_it != bd_map.end(); bd_it++)
    _mesh.boundary_info->set_label_to_id( (*bd_it).second, (*bd_it).first, false );


  // free tri_io structure
  triangulateio_finalize();

  MESSAGE<<"Tri3 mesh successfully generated.\n"<<std::endl; RECORD();

  STOP_LOG("do_mesh()", "MeshGenerator");

  return 0;

}





/* ----------------------------------------------------------------------------
 * over write the do_refine() virtual function in meshgen base class
 * this is a general mesh refine method for 2D mesh, both TRI and QUAD cells
 * the refine result is TRI mesh, of course.
 */
int MeshGeneratorTri3::do_refine(MeshRefinement & mesh_refinement)
{

  START_LOG("do_mesh()", "MeshGenerator");

  MESSAGE<<"Regriding Tri3 mesh by MeshGeneratorTri3...\n"; RECORD();

  // call MeshRefinement class to do FEM refine
  mesh_refinement.refine_and_coarsen_elements ();
  _mesh.find_neighbors();

  // we can simply flat refined mesh to level 0
  // since we do not need these information
  MeshTools::Modification::flatten(_mesh);

  // do necessarily prepare for call triangulate
  triangulateio_init();

  // prepare point information for triangle io
  {
    in.numberofpoints           = _mesh.max_node_id();
    in.numberofpointattributes  = 0;
    in.pointattributelist       = NULL;
#ifdef __triangle_h__
    in.pointlist                = (double *)calloc(in.numberofpoints*2, sizeof(double));
    in.pointmarkerlist          = (int *) calloc(in.numberofpoints, sizeof(int));
#else
    in.pointlist                = new double[in.numberofpoints*2];
    in.pointmarkerlist          = new int[in.numberofpoints];
#endif
    // push mesh nodes into triangle io
    MeshBase::node_iterator node_it = _mesh.active_nodes_begin();
    for(int i=0; node_it != _mesh.active_nodes_end() ; ++i, ++node_it)
    {
      in.pointlist[2*i+0] = (*(*node_it))(0);
      in.pointlist[2*i+1] = (*(*node_it))(1);
      in.pointmarkerlist[i] = _mesh.boundary_info->boundary_id(*node_it);
    }
  }

  // fill mesh elements into triangle io. in fact no need to do this.
  std::vector< std::vector<Point> > region_point;
  region_point.resize(_mesh.n_subdomains());
  {

    MeshBase::element_iterator elem_it = _mesh.active_elements_begin();
    for(; elem_it != _mesh.active_elements_end(); ++elem_it)
    {
      region_point[(*elem_it)->subdomain_id()].push_back( (*elem_it)->centroid () );
    }

  }


  // set region information
  std::map<unsigned int, std::pair<std::string, std::string> > region_info_map;
  {
    in.numberofholes   = 0;
    in.numberofregions = _mesh.n_subdomains();
#ifdef __triangle_h__
    in.regionlist = (double *) calloc(in.numberofregions*4, sizeof(double));
#else
    in.regionlist = new double[in.numberofregions*4];
#endif
    double *pregionlist =  in.regionlist;
    for(int i=0; i<in.numberofregions; i++)
    {
      *pregionlist++ = (region_point[i])[0](0);
      *pregionlist++ = (region_point[i])[0](1);
      *pregionlist++ = double(i);
      *pregionlist++ = 0;
    }

    for(unsigned int n_sub=0; n_sub<_mesh.n_subdomains(); n_sub++)
      region_info_map[n_sub] = std::make_pair( _mesh.subdomain_label_by_id(n_sub),_mesh.subdomain_material(n_sub) );

  }

  // set mesh boundary to triangle io
  std::map<short int, std::string> bd_label_map;
  std::set<short int>  bd_user_defined;
  {
    in.numberofsegments =  _mesh.boundary_info->n_boundary_conds();
#ifdef __triangle_h__
    in.segmentlist =  (int *) calloc(in.numberofsegments*2, sizeof(int));
    in.segmentmarkerlist = (int *) calloc(in.numberofsegments, sizeof(int));
#else
    in.segmentlist =  new int[in.numberofsegments*2];
    in.segmentmarkerlist = new int[in.numberofsegments];
#endif
    std::vector<unsigned int>       el;
    std::vector<unsigned short int> sl;
    std::vector<short int>          il;
    _mesh.boundary_info->build_side_list (el, sl, il);

    for(unsigned int i=0; i<el.size(); i++)
    {
      AutoPtr<Elem> side_elem = _mesh.elem(el[i])->build_side(sl[i]);
      in.segmentlist[2*i]   = side_elem->get_node(0)->id();
      in.segmentlist[2*i+1] = side_elem->get_node(1)->id();
      in.segmentmarkerlist[i]  = il[i];
    }

    const std::set<short int>& boundary_ids = _mesh.boundary_info->get_boundary_ids();
    std::set<short int>::const_iterator it = boundary_ids.begin();
    for(; it != boundary_ids.end(); ++it)
    {
      bd_label_map[*it] =  _mesh.boundary_info->get_label_by_id(*it);
      if( _mesh.boundary_info->boundary_id_has_user_defined_label(*it) )
        bd_user_defined.insert(*it);
    }
  }

  // rebuild the triangulation mesh
#ifdef __triangle_h__
  triangulate(const_cast<char *>(tri_cmd.c_str()), &in, &out, (struct triangulateio *) NULL);
#else
  ctri_triangulate(const_cast<char *>(tri_cmd.c_str()), &in, &out);
#endif

  // set boundary mark to output edge
  out_edge_table.clear();
  for(int i=0; i<out.numberofsegments; i++)
  {
    OutEdge edge;
    edge.pointno = out.numberofpoints;
    if(out.segmentlist[2*i+0] < out.segmentlist[2*i+1])
    {
      edge.p1 = out.segmentlist[2*i+0];
      edge.p2 = out.segmentlist[2*i+1];
    }
    else
    {
      edge.p1 = out.segmentlist[2*i+1];
      edge.p2 = out.segmentlist[2*i+0];
    }
    edge.mark = out.segmentmarkerlist[i];
    out_edge_table[edge] = edge.mark;
  }


  // clear old mesh structure
  _mesh.clear();

  // fill the _mesh structure
  _mesh.magic_num() = this->magic_num();

  // Build the nodes.
  _mesh.reserve_nodes( out.numberofpoints );
  _mesh.reserve_elem( out.numberoftriangles );

  for(int i = 0; i < out.numberofpoints; i++)
    _mesh.add_point(Point(out.pointlist[2*i+0], out.pointlist[2*i+1], 0.0));

  // build the elements
  for(int i = 0; i < out.numberoftriangles; i++)
  {
    Elem* elem = _mesh.add_elem(Elem::build(TRI3).release());

    elem->set_node(0) = _mesh.node_ptr( out.trianglelist[3*i+0] );
    elem->set_node(1) = _mesh.node_ptr( out.trianglelist[3*i+1] );
    elem->set_node(2) = _mesh.node_ptr( out.trianglelist[3*i+2] );

    elem->subdomain_id() = static_cast<int>(out.triangleattributelist[i]+0.5);

    int bc_index;

    if( (bc_index = get_bc_id(elem, out.trianglelist[3*i+0], out.trianglelist[3*i+1])) > 0 )
      _mesh.boundary_info->add_side(elem, 0, bc_index);

    if( (bc_index = get_bc_id(elem, out.trianglelist[3*i+1], out.trianglelist[3*i+2])) > 0 )
      _mesh.boundary_info->add_side(elem, 1, bc_index);

    if( (bc_index = get_bc_id(elem, out.trianglelist[3*i+2], out.trianglelist[3*i+0])) > 0 )
      _mesh.boundary_info->add_side(elem, 2, bc_index);

  }

  // write region label and material info to _mesh
  _mesh.set_n_subdomains() = region_info_map.size();
  std::map<unsigned int, std::pair<std::string, std::string> >::iterator region_it = region_info_map.begin();
  for(; region_it != region_info_map.end(); ++region_it)
  {
    _mesh.set_subdomain_label((*region_it).first, (*region_it).second.first );
    _mesh.set_subdomain_material((*region_it).first, (*region_it).second.second);
  }

  // write boundary label into mesh.boundary_info structure
  std::map<short int, std::string>::iterator bd_it = bd_label_map.begin();
  for(; bd_it != bd_label_map.end(); ++bd_it)
  {
    bool user_defined = (bd_user_defined.find((*bd_it).first) != bd_user_defined.end());
    _mesh.boundary_info->set_label_to_id( (*bd_it).first, (*bd_it).second, user_defined );
  }

  // free tri_io structure
  triangulateio_finalize();

  MESSAGE<<"Tri3 mesh successfully Regrided.\n"<<std::endl; RECORD();

  STOP_LOG("do_mesh()", "MeshGenerator");

  return 0;

}
