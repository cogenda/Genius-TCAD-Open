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
#include "genius_env.h"
#include "genius_common.h"
#include "log.h"
#include "unstructured_mesh.h"
#include "mesh_generation_cy3d.h"

#include "point.h"
#include "boundary_info.h"
#include "physical_unit.h"

#include "perf_log.h"


/* ----------------------------------------------------------------------------
 * set_region:  This function check and do REGION card
 * specify the material of a region
 */
int MeshGeneratorCylinder3D::set_region(const Parser::Card &c)
{
  const double pi = 3.14159265358979323846;
  double rmin = r_points.front();
  double rmax = r_points.back();
  double zmin = z_points.front();
  double zmax = z_points.back();

  // get parameter value from command structure
  SkeletonRegion3D region;

  // region label
  region.label = c.get_string("label","") ;

  if(!c.is_parameter_exist("material"))
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " REGION: Material should be given for each region!\n";
    RECORD();
  }

  // region material
  region.material= c.get_string("material","") ;

  // r range
  double rmax_ = c.get_real("r.max",rmax);
  double rmin_ = c.get_real("r.min",rmin);
  region.rmax = find_r(rmax_);
  region.rmin = find_r(rmin_);

  // theta range
  double thetamax_ = c.get_real("theta.max", 360.0);
  double thetamin_ = c.get_real("theta.min", 0.0);

  region.thetamax = find_theta(region.rmin, thetamax_);
  region.thetamin = find_theta(region.rmin, thetamin_);

  // z range
  double zmax_ = c.get_real("z.max", zmax);
  double zmin_ = c.get_real("z.min", zmin);

  region.zmax = find_z(zmax_);
  region.zmin = find_z(zmin_);

  if( region.rmin < 0  ||  region.rmax < 0 || region.rmin >= region.rmax )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " REGION " << region.label << ": Can't locate radial position of region!\n";
    RECORD();
    return 1;
  }

  if(region.rmin == region.rmax || region.thetamin == region.thetamax)
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " REGION " << region.label << ": Region area is zero!\n";
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

  if(region.thetamin != region.thetamax &&  region.thetamin+360.0 != region.thetamax )
  {
    {
      SkeletonSector face;
      face.rmin=region.rmin;
      face.rmax=region.rmax;
      face.thetamin=region.thetamin;
      face.thetamax=region.thetamin;
      face.zmin=region.zmin;
      face.zmax=region.zmax;
      face.sector_label=region.label+"_Neumann";
      face_theta_array1d.push_back(face);
    }

    {
      SkeletonSector face;
      face.rmin=region.rmin;
      face.rmax=region.rmax;
      face.thetamin=region.thetamax;
      face.thetamax=region.thetamax;
      face.zmin=region.zmin;
      face.zmax=region.zmax;
      face.sector_label=region.label+"_Neumann";
      face_theta_array1d.push_back(face);
    }
  }

  if(region.rmin > 0.0)
  {
    SkeletonSector face;
    face.rmin=region.rmin;
    face.rmax=region.rmin;
    face.thetamin=region.thetamin;
    face.thetamax=region.thetamax;
    face.zmin=region.zmin;
    face.zmax=region.zmax;
    face.sector_label=region.label+"_Neumann";
    face_r_array1d.push_back(face);
  }

  {
    SkeletonSector face;
    face.rmin=region.rmax;
    face.rmax=region.rmax;
    face.thetamin=region.thetamin;
    face.thetamax=region.thetamax;
    face.zmin=region.zmin;
    face.zmax=region.zmax;
    face.sector_label=region.label+"_Neumann";
    face_r_array1d.push_back(face);
  }

  {
    SkeletonSector face;
    face.rmin=region.rmin;
    face.rmax=region.rmax;
    face.thetamin=region.thetamin;
    face.thetamax=region.thetamax;
    face.zmin=region.zmin;
    face.zmax=region.zmin;
    face.sector_label=region.label+"_Neumann";
    face_z_array1d.push_back(face);
  }

  {
    SkeletonSector face;
    face.rmin=region.rmin;
    face.rmax=region.rmax;
    face.thetamin=region.thetamin;
    face.thetamax=region.thetamax;
    face.zmin=region.zmax;
    face.zmax=region.zmax;
    face.sector_label=region.label+"_Neumann";
    face_z_array1d.push_back(face);
  }


  face_id_map.insert(std::make_pair(region.label+"_Neumann", face_id_map.size()+1));

  return 0;

}



/* ----------------------------------------------------------------------------
 * get the region index of a cell by its bound
 * if region is overlaped, return the last defined region id
 */
int  MeshGeneratorCylinder3D::get_cell_region_id(const Point &p)
{
  const double pi = 3.14159265358979323846;

  double r = sqrt(p.x()*p.x() + p.y()*p.y());
  double theta = 0;
  if( p.x() == 0.0)
  {
    if(p.y() > 0.0) theta = 0.5*pi;
    if(p.y() < 0.0) theta = 1.5*pi;
  }
  else
  {
    theta=atan(p.y()/p.x());
    if(p.x()<0) theta = theta + pi;
  }

  theta *= 180/pi;
  double z = p.z();

  int id = 0;
  for(size_t i=0; i < region_array1d.size(); i++)
  {
    if( r >= region_array1d[i].rmin && r <= region_array1d[i].rmax &&
        ((theta >= region_array1d[i].thetamin && theta <= region_array1d[i].thetamax) ||
         (theta+360 >= region_array1d[i].thetamin && theta+360 <= region_array1d[i].thetamax)) &&
        z >= region_array1d[i].zmin && z<=region_array1d[i].zmax
      )
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
int MeshGeneratorCylinder3D::set_face(const Parser::Card &c)
{
  const double pi = 3.14159265358979323846;
  std::string  label = c.get_string("label","");

  double rmin = r_points.front();
  double rmax = r_points.back();
  double thetamin=0, thetamax=360;
  double zmin = z_points.front();
  double zmax = z_points.back();

  if(c.is_parameter_exist("r"))
  {
    double r = c.get_real("r", 0);
    rmin = rmax = r;
  }

  if(c.is_parameter_exist("r.max"))
  {
    rmax = c.get_real("r.max",rmax);
  }
  if(c.is_parameter_exist("r.min"))
  {
    rmin = c.get_real("r.min",rmin);
  }

  if(c.is_parameter_exist("theta"))
  {
    double theta = c.get_real("theta", 0);
    thetamin = thetamax = theta;
  }

  if(c.is_parameter_exist("theta.max"))
  {
    thetamax = c.get_real("theta.max",360.0);
  }

  if(c.is_parameter_exist("theta.min"))
  {
    thetamin = c.get_real("theta.min",0);
  }

  if(c.is_parameter_exist("z"))
  {
    double z = c.get_real("z", 0);
    zmin = zmax = z;
  }

  if(c.is_parameter_exist("z.max"))
  {
    zmax = c.get_real("z.max", zmax);
  }

  if(c.is_parameter_exist("z.min"))
  {
    zmin = c.get_real("z.min", zmin);
  }


  if(c.is_parameter_exist("location"))
  {
    if(c.is_enum_value("location","top"))
      zmin=zmax=z_points.back();
    else if(c.is_enum_value("location","bottom"))
      zmin=zmax=z_points.front();
  }

  SkeletonSector face;
  face.rmin= find_r(rmin);
  face.rmax= find_r(rmax);
  face.thetamin=find_theta(face.rmin, thetamin);
  face.thetamax=find_theta(face.rmin, thetamax);
  face.zmin = find_z(zmin);
  face.zmax = find_z(zmax);
  face.sector_label=label;

  if(rmin == rmax)
    face_r_array1d.push_back(face);

  if(thetamin == thetamax)
    face_theta_array1d.push_back(face);

  if(zmin == zmax)
    face_z_array1d.push_back(face);

  face_id_map.insert(std::make_pair(label, face_id_map.size() +1));

  return 0;

}

int  MeshGeneratorCylinder3D::get_face_id_r(const Point &p1, const Point &p2, const Point &p3, const Point &p4)
{
  const double pi = 3.14159265358979323846;
  Point mid = 0.25*(p1+p2+p3+p4);

  double r = sqrt(p1.x()*p1.x() + p1.y()*p1.y());
  double theta = 0;
  if( mid.x() == 0.0)
  {
    if(mid.y() > 0.0) theta = 0.5*pi;
    if(mid.y() < 0.0) theta = 1.5*pi;
  }
  else
  {
    theta=atan(mid.y()/mid.x());
    if(mid.x()<0) theta = theta + pi;
  }

  theta *= 180/pi;
  double z = mid.z();

  {
    int id = 0;
    for(unsigned int i=0; i<face_r_array1d.size(); ++i)
    {
      const SkeletonSector & face = face_r_array1d[i];
      assert(face.rmin == face.rmax);
      if( std::abs(r - face.rmin) < 1e-10 &&
          ((theta > face.thetamin && theta < face.thetamax) ||
           (theta+360 > face.thetamin && theta+360 < face.thetamax) ) &&
          z > face.zmin && z < face.zmax
        )
      {
        id = face_id_map.find(face.sector_label)->second;
      }
    }
    return id;
  }

  return 0;

}




int  MeshGeneratorCylinder3D::get_face_id_theta(double rmid, const Point &mid)
{
  const double pi = 3.14159265358979323846;
  double theta = 0;
  if( mid.x() == 0.0)
  {
    if(mid.y() > 0.0) theta = 0.5*pi;
    if(mid.y() < 0.0) theta = 1.5*pi;
  }
  else
  {
    theta=atan(mid.y()/mid.x());
    if(mid.x()<0) theta = theta + pi;
  }

  theta *= 180/pi;
  double z = mid.z();

  {
    int id = 0;
    for(unsigned int i=0; i<face_theta_array1d.size(); ++i)
    {
      const SkeletonSector & face = face_theta_array1d[i];
      assert(face.thetamin == face.thetamax);
      if( std::abs(theta - face.thetamin) < 1e-3 && rmid > face.rmin && rmid < face.rmax && z > face.zmin && z < face.zmax  )
        id = face_id_map.find(face.sector_label)->second;
    }
    return id;
  }


  return 0;
}



int  MeshGeneratorCylinder3D::get_face_id_z(double rmid, const Point &mid)
{
  const double pi = 3.14159265358979323846;
  double theta = 0;
  if( mid.x() == 0.0)
  {
    if(mid.y() > 0.0) theta = 0.5*pi;
    if(mid.y() < 0.0) theta = 1.5*pi;
  }
  else
  {
    theta=atan(mid.y()/mid.x());
    if(mid.x()<0) theta = theta + pi;
  }

  theta *= 180/pi;
  double z = mid.z();

  {
    int id = 0;
    for(unsigned int i=0; i<face_z_array1d.size(); ++i)
    {
      const SkeletonSector & face = face_z_array1d[i];
      assert(face.zmin == face.zmax);
      if( std::abs(z - face.zmin) < 1e-10 &&
          ((theta > face.thetamin && theta < face.thetamax) ||
           (theta+360 > face.thetamin && theta+360 < face.thetamax) ) &&
          rmid > face.rmin && rmid < face.rmax
        )
        id = face_id_map.find(face.sector_label)->second;
    }

    return id;
  }


  return 0;
}

/* ----------------------------------------------------------------------------
 * over write the virtual function in meshgen base class
 */
int MeshGeneratorCylinder3D::do_mesh()
{
  START_LOG("do_mesh()", "MeshGenerator");

  MESSAGE<<"Generating mesh by MeshGeneratorCylinder3D...\n"; RECORD();

  double dscale = PhysicalUnit::um;

  // build x-y-z grid skeleton
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if(c.key() == "R.MESH")   // It's a R.MESH card
      if(set_rt_line(c)) return 1;

    if(c.key() == "Z.MESH")   // It's a Z.MESH card
      if(set_z_line(c)) return 1;
  }

  // set region here.
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();
    if(c.key() == "REGION")
      if( set_region(c) ) return 1;
  }

  // set face here.
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();
    if(c.key() == "FACE")
      if( set_face(c) ) return 1;
  }

  // mesh points
  std::map< std::pair<unsigned int, unsigned int>,  std::vector<Node *> >  all_nodes;

  // fill the _mesh structure
  _mesh.magic_num() = this->magic_num();

  for(unsigned int i=0; i<r_points.size(); i++)
  {
    double radiuse = r_points[i];
    int div = theta_points[i];

    if(i==0)
    {
      for(unsigned int j=0; j<z_points.size(); ++j)
      {
        std::vector<Node *> & nodes = all_nodes[std::make_pair(i, j)];
        Node * node =  _mesh.add_point(Point(0, 0, z_points[j]));
        nodes.push_back(node);
      }
    }

    if(i==1)
    {
      double radiuse_in = r_points[0];
      double r_mid = 0.5*(radiuse_in+radiuse);

      for(unsigned int j=0; j<z_points.size(); ++j)
      {
        std::vector<Node *> & nodes = all_nodes[std::make_pair(i, j)];
        std::vector<Point> pts = divide_circle(Point(), radiuse, div);
        for(unsigned int n=0; n<pts.size()-1; ++n)
        {
          Point p(pts[n].x(), pts[n].y(), z_points[j]);
          Node * node =  _mesh.add_point(p);
          nodes.push_back(node);
        }
        nodes.push_back(nodes.front());
      }

      for(unsigned int j=0; j<z_points.size()-1; ++j)
      {
        Node * in_node_bot = all_nodes[std::make_pair(0, j)][0];
        Node * in_node_top = all_nodes[std::make_pair(0, j+1)][0];
        const std::vector<Node *> & nodes_bot = all_nodes[std::make_pair(i, j)];
        const std::vector<Node *> & nodes_top = all_nodes[std::make_pair(i, j+1)];
        for(unsigned int n=0; n<nodes_top.size()-1; ++n)
        {
          Elem* elem = _mesh.add_elem(Elem::build(PRISM6).release());
          elem->set_node(0) = in_node_bot;
          elem->set_node(1) = nodes_bot[n];
          elem->set_node(2) = nodes_bot[n+1];
          elem->set_node(3) = in_node_top;
          elem->set_node(4) = nodes_top[n];
          elem->set_node(5) = nodes_top[n+1];

          Point c = elem->centroid();
          Point cp(c.x(), c.y(), 0.0);
          cp = cp.unit()*r_mid;
          elem->subdomain_id() = get_cell_region_id(Point(cp.x(), cp.y(), c.z()));

          Point s1_mid = 0.25*(*in_node_bot + *nodes_bot[n] + *in_node_top + *nodes_top[n]);
          int s1_id = get_face_id_theta(r_mid, s1_mid);
          if ( s1_id > 0 )
            _mesh.boundary_info->add_side(elem, 1, static_cast<short int>(s1_id) );

          Point s3_mid = 0.25*(*in_node_bot + *nodes_bot[n+1] + *in_node_top + *nodes_top[n+1]);
          int s3_id = get_face_id_theta(r_mid, s3_mid);
          if ( s3_id > 0 )
            _mesh.boundary_info->add_side(elem, 3, static_cast<short int>(s3_id) );

          int s2_id = get_face_id_r(*nodes_bot[n], *nodes_bot[n+1], *nodes_top[n], *nodes_top[n+1]);
          if ( s2_id > 0 )
            _mesh.boundary_info->add_side(elem, 2, static_cast<short int>(s2_id) );

          Point s0_mid = (*in_node_bot + *nodes_bot[n] + *nodes_bot[n+1])/3.0;
          int s0_id = get_face_id_z(r_mid, s0_mid);
          if ( s0_id > 0 )
            _mesh.boundary_info->add_side(elem, 0, static_cast<short int>(s0_id) );

          Point s4_mid = (*in_node_top + *nodes_top[n] + *nodes_top[n+1])/3.0;
          int s4_id = get_face_id_z(r_mid, s4_mid);
          if ( s4_id > 0 )
            _mesh.boundary_info->add_side(elem, 4, static_cast<short int>(s4_id) );

        }
      }
    }

    if(i>1)
    {
      double radiuse_in = r_points[i-1];
      double r_mid = 0.5*(radiuse_in+radiuse);

      for(unsigned int j=0; j<z_points.size(); ++j)
      {
        std::vector<Node *> & nodes = all_nodes[std::make_pair(i, j)];

        std::vector<Point> pts = divide_circle(Point(), radiuse, div);
        for(unsigned int n=0; n<pts.size()-1; ++n)
        {
          Point p(pts[n].x(), pts[n].y(), z_points[j]);
          Node * node =  _mesh.add_point(p);
          nodes.push_back(node);
        }
        nodes.push_back(nodes.front());
      }

      for(unsigned int j=0; j<z_points.size()-1; ++j)
      {
        const std::vector<Node *> & in_nodes_bot = all_nodes[std::make_pair(i-1, j)];
        const std::vector<Node *> & in_nodes_top = all_nodes[std::make_pair(i-1, j+1)];

        const std::vector<Node *> & nodes_bot = all_nodes[std::make_pair(i, j)];
        const std::vector<Node *> & nodes_top = all_nodes[std::make_pair(i, j+1)];

        if(in_nodes_top.size() == nodes_top.size())
        {
          for(unsigned int n=0; n<in_nodes_top.size()-1; ++n)
          {
            Elem* elem = _mesh.add_elem(Elem::build(HEX8).release());

            elem->set_node(0) = in_nodes_bot[n+1];
            elem->set_node(1) = in_nodes_bot[n];
            elem->set_node(2) = nodes_bot[n];
            elem->set_node(3) = nodes_bot[n+1];
            elem->set_node(4) = in_nodes_top[n+1];
            elem->set_node(5) = in_nodes_top[n];
            elem->set_node(6) = nodes_top[n];
            elem->set_node(7) = nodes_top[n+1];

            Point c = elem->centroid();
            Point cp(c.x(), c.y(), 0.0);
            cp = cp.unit()*r_mid;
            elem->subdomain_id() = get_cell_region_id(Point(cp.x(), cp.y(), c.z()));

            int s1_id = get_face_id_r(*in_nodes_bot[n+1], *in_nodes_bot[n], *in_nodes_top[n+1], *in_nodes_top[n]);
            if ( s1_id > 0 )
              _mesh.boundary_info->add_side(elem, 1, static_cast<short int>(s1_id) );

            int s3_id = get_face_id_r(*nodes_bot[n], *nodes_bot[n+1], *nodes_top[n], *nodes_top[n+1]);
            if ( s3_id > 0 )
              _mesh.boundary_info->add_side(elem, 3, static_cast<short int>(s3_id) );

            Point s2_mid = 0.25*(*in_nodes_bot[n] + *nodes_bot[n] + *in_nodes_top[n] + *nodes_top[n]);
            int s2_id = get_face_id_theta(r_mid, s2_mid);
            if ( s2_id > 0 )
              _mesh.boundary_info->add_side(elem, 2, static_cast<short int>(s2_id) );

            Point s4_mid = 0.25*(*in_nodes_bot[n+1] + *nodes_bot[n+1], *in_nodes_top[n+1] + *nodes_top[n+1]);
            int s4_id = get_face_id_theta(r_mid, s4_mid);
            if ( s4_id > 0 )
              _mesh.boundary_info->add_side(elem, 4, static_cast<short int>(s4_id) );

            Point s0_mid = 0.25*(*in_nodes_bot[n] + *in_nodes_bot[n+1] + *nodes_bot[n] + *nodes_bot[n+1]);
            int s0_id = get_face_id_z(r_mid, s0_mid);
            if ( s0_id > 0 )
              _mesh.boundary_info->add_side(elem, 0, static_cast<short int>(s0_id) );

            Point s5_mid = 0.25*(*in_nodes_top[n] + *in_nodes_top[n+1] + *nodes_top[n] + *nodes_top[n+1]);
            int s5_id = get_face_id_z(r_mid, s5_mid);
            if ( s5_id > 0 )
              _mesh.boundary_info->add_side(elem, 5, static_cast<short int>(s5_id) );
          }
        }


        if((in_nodes_top.size()-1)*2 == nodes_top.size()-1)
        {
          for(unsigned int n=0; n<in_nodes_top.size()-1; ++n)
          {
            {
              Elem* elem1 = _mesh.add_elem(Elem::build(PRISM6).release());
              elem1->set_node(0) = in_nodes_bot[n];
              elem1->set_node(1) = nodes_bot[2*n];
              elem1->set_node(2) = nodes_bot[2*n+1];
              elem1->set_node(3) = in_nodes_top[n];
              elem1->set_node(4) = nodes_top[2*n];
              elem1->set_node(5) = nodes_top[2*n+1];

              Point c = elem1->centroid();
              Point cp(c.x(), c.y(), 0.0);
              cp = cp.unit()*r_mid;
              elem1->subdomain_id() = get_cell_region_id(Point(cp.x(), cp.y(), c.z()));

              Point s1_mid = 0.25*(*in_nodes_bot[n] + *nodes_bot[2*n] + *in_nodes_top[n] + *nodes_top[2*n]);
              int s1_id = get_face_id_theta(r_mid, s1_mid);
              if ( s1_id > 0 )
                _mesh.boundary_info->add_side(elem1, 1, static_cast<short int>(s1_id) );

              int s2_id = get_face_id_r(*nodes_bot[2*n], *nodes_top[2*n], *nodes_bot[2*n+1], *nodes_top[2*n+1]);
              if ( s2_id > 0 )
                _mesh.boundary_info->add_side(elem1, 2, static_cast<short int>(s2_id) );

              Point s0_mid = (*in_nodes_bot[n] + *nodes_bot[2*n] + *nodes_bot[2*n+1])/3.0;
              int s0_id = get_face_id_z(r_mid, s0_mid);
              if ( s0_id > 0 )
                _mesh.boundary_info->add_side(elem1, 0, static_cast<short int>(s0_id) );

              Point s4_mid = (*in_nodes_top[n] + *nodes_top[2*n] + *nodes_top[2*n+1])/3.0;
              int s4_id = get_face_id_z(r_mid, s4_mid);
              if ( s4_id > 0 )
                _mesh.boundary_info->add_side(elem1, 4, static_cast<short int>(s4_id) );
            }

            {
              Elem* elem2 = _mesh.add_elem(Elem::build(PRISM6).release());
              elem2->set_node(0) = in_nodes_bot[n+1];
              elem2->set_node(1) = nodes_bot[2*n+1];
              elem2->set_node(2) = nodes_bot[2*n+2];
              elem2->set_node(3) = in_nodes_top[n+1];
              elem2->set_node(4) = nodes_top[2*n+1];
              elem2->set_node(5) = nodes_top[2*n+2];

              Point c = elem2->centroid();
              Point cp(c.x(), c.y(), 0.0);
              cp = cp.unit()*r_mid;
              elem2->subdomain_id() = get_cell_region_id(Point(cp.x(), cp.y(), c.z()));

              int s2_id = get_face_id_r(*nodes_bot[2*n+1], *nodes_top[2*n+1], *nodes_bot[2*n+2], *nodes_top[2*n+2]);
              if ( s2_id > 0 )
                _mesh.boundary_info->add_side(elem2, 2, static_cast<short int>(s2_id) );

              Point s3_mid = 0.25*(*nodes_bot[2*n+2] + *nodes_top[2*n+2] + *in_nodes_bot[n+1] + *in_nodes_top[n+1]);
              int s3_id = get_face_id_theta(r_mid, s3_mid);
              if ( s3_id > 0 )
                _mesh.boundary_info->add_side(elem2, 3, static_cast<short int>(s3_id) );

              Point s0_mid = (*in_nodes_bot[n+1] + *nodes_bot[2*n+1] + *nodes_bot[2*n+2])/3.0;
              int s0_id = get_face_id_z(r_mid, s0_mid);
              if ( s0_id > 0 )
                _mesh.boundary_info->add_side(elem2, 0, static_cast<short int>(s0_id) );

              Point s4_mid = (*in_nodes_top[n+1] + *nodes_top[2*n+1] + *nodes_top[2*n+2])/3.0;
              int s4_id = get_face_id_z(r_mid, s4_mid);
              if ( s4_id > 0 )
                _mesh.boundary_info->add_side(elem2, 4, static_cast<short int>(s4_id) );
            }

            {
              Elem* elem3 = _mesh.add_elem(Elem::build(PRISM6).release());
              elem3->set_node(0) = in_nodes_bot[n+1];
              elem3->set_node(1) = in_nodes_bot[n];
              elem3->set_node(2) = nodes_bot[2*n+1];
              elem3->set_node(3) = in_nodes_top[n+1];
              elem3->set_node(4) = in_nodes_top[n];
              elem3->set_node(5) = nodes_top[2*n+1];

              Point c = elem3->centroid();
              Point cp(c.x(), c.y(), 0.0);
              cp = cp.unit()*r_mid;
              elem3->subdomain_id() = get_cell_region_id(Point(cp.x(), cp.y(), c.z()));

              int s1_id = get_face_id_r(*in_nodes_bot[n], *in_nodes_top[n], *in_nodes_bot[n+1], *in_nodes_top[n+1]);
              if ( s1_id > 0 )
                _mesh.boundary_info->add_side(elem3, 1, static_cast<short int>(s1_id) );

              Point s0_mid = (*in_nodes_bot[n] + *nodes_bot[2*n+1] + *in_nodes_bot[n+1])/3.0;
              int s0_id = get_face_id_z(r_mid, s0_mid);
              if ( s0_id > 0 )
                _mesh.boundary_info->add_side(elem3, 0, static_cast<short int>(s0_id) );

              Point s4_mid = (*in_nodes_top[n] + *nodes_top[2*n+1] + *in_nodes_top[n+1])/3.0;
              int s4_id = get_face_id_z(r_mid, s4_mid);
              if ( s4_id > 0 )
                _mesh.boundary_info->add_side(elem3, 4, static_cast<short int>(s4_id) );
            }
          }
        }
      }//for z_points
    }

  }// for r_points

  // scale

  for (unsigned int n=0; n<_mesh.n_nodes(); n++)
  {
    _mesh.node(n) *= dscale;
  }

  // write region label and material info to _mesh

  _mesh.set_n_subdomains() = region_array1d.size();

  for(size_t r=0; r < region_array1d.size(); r++)
  {
    _mesh.set_subdomain_label(r, region_array1d[r].label );
    _mesh.set_subdomain_material(r, region_array1d[r].material);
  }


  // write boundary label into mesh.boundary_info structure
  std::map<std::string, int>::const_iterator seg_it = face_id_map.begin();
  for(; seg_it!=face_id_map.end(); seg_it++)
  {
    _mesh.boundary_info->set_label_to_id( seg_it->second, seg_it->first, true);
  }


  //the interface information should be set here
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
        bd_index = face_id_map.size() + bd_map.size() + 1;
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
    _mesh.boundary_info->set_label_to_id( (*bd_it).second, (*bd_it).first, false );


  STOP_LOG("do_mesh()", "MeshGenerator");

  return 0;

}




