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



// C++ includes
#include <cmath> // for std::sqrt


// Local includes
#include "genius_env.h"
#include "genius_common.h"
#include "log.h"
#include "unstructured_mesh.h"
#include "mesh_generation_quad4.h"

#include "point.h"
#include "boundary_info.h"


#include "perf_log.h"


/* ----------------------------------------------------------------------------
 * build the rectangle Skeleton mesh in XY plane
 */
void MeshGeneratorQuad4::build_rectangle_mesh()
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
 * set_region:  This function check and do REGION card
 * specify the material of a region
 */
int MeshGeneratorQuad4::set_region(const Parser::Card &c)
{

  // get parameter value from command structure
  SkeletonRegion2D region;

  // some problem: ixxx is unsigned int, but get_int may return nagetive vlaue...
  // however, for this case, the unsigned int value will be very large, much large than
  // IX/IY
  region.ixmin = c.get_int("ix.min",0,"ix.left");
  region.ixmax = c.get_int("ix.max",IX-1,"ix.right");
  region.iymin = c.get_int("iy.min",0,"iy.top");
  region.iymax = c.get_int("iy.max",IY-1,"iy.bottom");
  region.label = c.get_string("label","") ;

  region.material= c.get_string("material","") ;

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
 * get the region index of a quad4 cell by its bound
 * if region is overlaped, return the last defined region id
 */
int  MeshGeneratorQuad4::get_cell_region_id(unsigned int ixmin, unsigned int  ixmax,
                                            unsigned int iymin, unsigned int  iymax)
{
  int id = 0;

  for(size_t i=0; i < region_array1d.size(); i++)
  {
    if( ixmin >= region_array1d[i].ixmin && ixmax <= region_array1d[i].ixmax &&
        iymin >= region_array1d[i].iymin && iymax <= region_array1d[i].iymax )
    {
      id = i;
    }
  }

  return id;

}



/* ----------------------------------------------------------------------------
 * set_face:  This function check and do FACE card
 * specify the label and location of a face.
 * these segments may be specified as electrode later by BOUNDARY card.
 */
int MeshGeneratorQuad4::set_face(const Parser::Card &c)
{
  int ixmin,ixmax,iymin,iymax;
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
  }

  if(c.is_parameter_exist("direction"))
  {
    if(c.is_enum_value("direction","YNORM"))
      iymin=iymax=iy;
    else if(c.is_enum_value("direction","XNORM"))
      ixmin=ixmax=ix;
  }

  if(ixmin!=ixmax && iymin!=iymax)
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " FACE: face must along grid lines!\n";
    RECORD();
    return 1;
  }

  if(ixmin<0 || ixmax>static_cast<int>(IX-1) || ixmin > ixmax )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " FACE: face xline exceed range!\n";
    RECORD();
    return 1;
  }

  if(iymin<0 || iymax>static_cast<int>(IY-1) || iymin > iymax )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " FACE: face yline exceed range!\n";
    RECORD();
    return 1;
  }


  make_face(ixmin,ixmax,iymin,iymax,label);

  return 0;

}





/* ----------------------------------------------------------------------------
 * make_face:  set labeled faces
 */
int  MeshGeneratorQuad4::make_face(unsigned int ixmin,unsigned int ixmax,
                                   unsigned int iymin,unsigned int iymax,
                                   const std::string &label)
{
  int face_mark = face_array1d.size() + 1;

  SkeletonFace face;
  face.face_label=label;
  face.face_mark=face_mark;
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
    for(unsigned int i=ixmin;i<ixmax;)
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
    for(unsigned int i=iymin;i<iymax;)
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
 * make_region_boundary:  set all the boundary face of regions to Neumann
 */
int  MeshGeneratorQuad4::make_region_boundary()
{
  // set segment bounding
  std::map<SkeletonEdge, int, lt_edge>::iterator pt;

  for(size_t r=0;r<region_array1d.size();r++)
  {
    //if(region_array1d[r].shape! = Rectangle) continue;

    SkeletonFace  face;
    face.face_label=region_array1d[r].label+"_Neumann";
    face.face_mark=r+1;
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
 * over write the virtual function in meshgen base class
 */
int MeshGeneratorQuad4::do_mesh()
{
  START_LOG("do_mesh()", "MeshGenerator");

  MESSAGE<<"Generating Quad4 mesh by MeshGeneratorQuad4...\n"; RECORD();

  // reset default dimension scale value if user provide it
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if(c.key() == "GLOBAL")
    {
      dscale = pow(c.get_real("dopingscale",1e18), 1.0/3)/1e4;
      _decks.breakloop();
      break;
    }
  }

  // build x-y-z grid skeleton
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if(c.key() == "X.MESH")   // It's a X.MESH card
      if(set_x_line(c)) return 1;
    if(c.key() == "Y.MESH")   // It's a Y.MESH card
      if(set_y_line(c)) return 1;
    if(c.key() == "Z.MESH")   // It's a Z.MESH card, we do not like it.
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " Z.MESH: Quad4 mesh generator not support this statement!\n";
      RECORD();
      return 1;
    }
  }

  //after the initializtion of x and y mesh lines, we can build rectangle now.
  build_rectangle_mesh();

  //use REGION card to set material of the region.
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if( c.key() == "ELIMINATE" )   //  eliminate?
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " ELIMINATE: Quad4 mesh generator not support this statement!\n";
      RECORD();
      return 1;
    }

    if( c.key() == "SPREAD" || c.key() == "SPREAD3D" )      //  spread?
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " SPREAD: Quad4 mesh generator not support this statement!\n";
      RECORD();
      return 1;
    }

    if(c.key() == "REGION")      // set region here.
      if( set_region(c) ) return 1;
  }

  // set all the boundary face to Neumann, this may be overwrite
  // by following "FACE" card
  make_region_boundary();

  // set labeled faces(most of them are electrode) defined by user here.
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if(c.key() == "FACE")      // set face here.
      if( set_face(c) ) return 1;
  }

  // fill the _mesh structure
  _mesh.magic_num() = this->magic_num();

  //write to mesh
  _mesh.reserve_nodes( IX*IY );
  _mesh.reserve_elem( (IX-1)*(IY-1));

  // Build the nodes. we do dimension scale here
  {
      for (unsigned int j=0; j<IY; j++)
        for (unsigned int i=0; i<IX; i++)
        {
          _mesh.add_point(Point(point_array3d[0][j][i].x*dscale,
                                point_array3d[0][j][i].y*dscale,
                                0));
        }

  }

  //set the elem , the region(subdomain) id is also set
    for (unsigned int j=0; j<IY-1; j++)
      for (unsigned int i=0; i<IX-1; i++)
      {
        Elem* elem = _mesh.add_elem(Elem::build(QUAD4).release());


        elem->set_node(0) = _mesh.node_ptr( i + IX*(j)             );
        elem->set_node(1) = _mesh.node_ptr( (i+1) + IX*(j)         );
        elem->set_node(2) = _mesh.node_ptr( (i+1) + IX*((j+1))     );
        elem->set_node(3) = _mesh.node_ptr( i + IX*((j+1))         );

        // set subdomain id
        elem->subdomain_id() = get_cell_region_id (i,i+1,j,j+1);

        // set boundary id, however, all the region faces are set to Neumann at present,
        // we should pick up the face lies on interface of two regions later.
        SkeletonEdge  edge;
        edge.IX = IX;
        edge.IY = IY;

        //top face
        {
          edge.p1[0] = i;
          edge.p1[1] = j+1;
          edge.p2[0] = i+1;
          edge.p2[1] = j+1;
          if ( edge_table.find(edge) != edge_table.end() )
            _mesh.boundary_info->add_side(elem, 2, static_cast<short int>(edge_table[edge]) );
        }

        //bottom face
        {
          edge.p1[0] = i;
          edge.p1[1] = j;
          edge.p2[0] = i+1;
          edge.p2[1] = j;
          if ( edge_table.find(edge) != edge_table.end() )
            _mesh.boundary_info->add_side(elem, 0, static_cast<short int>(edge_table[edge]) );
        }

        //left face
        {
          edge.p1[0] = i;
          edge.p1[1] = j;
          edge.p2[0] = i;
          edge.p2[1] = j+1;
          if ( edge_table.find(edge) != edge_table.end() )
            _mesh.boundary_info->add_side(elem, 3, static_cast<short int>(edge_table[edge]) );
        }

        //right face
        {
          edge.p1[0] = i+1;
          edge.p1[1] = j;
          edge.p2[0] = i+1;
          edge.p2[1] = j+1;
          if ( edge_table.find(edge) != edge_table.end() )
            _mesh.boundary_info->add_side(elem, 1, static_cast<short int>(edge_table[edge]) );
        }

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
    _mesh.boundary_info->set_label_to_id( face_array1d[f].face_mark, face_array1d[f].face_label );

  //the interface information should be set here

  std::vector<unsigned int>       elems;
  std::vector<unsigned short int> sides;
  std::vector<short int>          bds;
  std::map<const std::string, short int> bd_map;
  typedef std::map<const std::string, short int>::iterator Bd_It;

  // get all the boundary element
  _mesh.boundary_info->build_side_list (elems, sides, bds);

  //build neighbor information for mesh. then elem->neighbor() is functional
  _mesh.boundary_info->find_neighbors();


  for (size_t nbd=0; nbd<elems.size(); nbd++ )
  {
    // get the element which has boundary/interface side
    const Elem* elem = _mesh.elem(elems[nbd]);
    //is it an interface side && not a label face?
    if( elem->neighbor(sides[nbd])!=NULL && bds[nbd]<=static_cast<short int>(region_array1d.size()) )
    {
       // the element and its neighbor should in diffetent subdomain
       unsigned int sbd_id1 = elem->subdomain_id();
       unsigned int sbd_id2 = elem->neighbor(sides[nbd])->subdomain_id();
       assert (sbd_id1 != sbd_id2);

       //remove the side from boundary
       _mesh.boundary_info->remove(elem, sides[nbd]);

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

       // add element to boundary again with new bd_index
       _mesh.boundary_info->add_side(elem, sides[nbd], bd_index);
    }
  }

  // after the previous remove and insert operation, the number of boundary id maybe changed. renumber it here.
  _mesh.boundary_info->rebuild_ids();

  // write down new labels
  Bd_It bd_it = bd_map.begin();
  for(; bd_it != bd_map.end(); bd_it++)
    _mesh.boundary_info->set_label_to_id( (*bd_it).second, (*bd_it).first );

  // addtinal work: rotate the mesh as y positive point to top direction
  //MeshTools::Modification::rotate (_mesh, 180.0, 0.0, 0.0);

  MESSAGE<<"Quad4 mesh successfully generated.\n"<<std::endl; RECORD();

  STOP_LOG("do_mesh()", "MeshGenerator");

  return 0;
}


/* ----------------------------------------------------------------------------
 * over write the do_refine() virtual function in meshgen base class
 */
int MeshGeneratorQuad4::do_refine(MeshRefinement & )
{

  START_LOG("do_mesh()", "MeshGenerator");

  MESSAGE<<"Warning: Regriding Quad4 mesh is not supported at present, return now.\n"; RECORD();

  STOP_LOG("do_mesh()", "MeshGenerator");

  return 0;

}

