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
#include "mesh_generation_cy2d.h"

#include "point.h"
#include "boundary_info.h"
#include "physical_unit.h"

#include "perf_log.h"


/* ----------------------------------------------------------------------------
 * set_region:  This function check and do REGION card
 * specify the material of a region
 */
int MeshGeneratorCylinder2D::set_region(const Parser::Card &c)
{
  const double pi = 3.14159265358979323846;
  double rmin = r_points.front();
  double rmax = r_points.back();

  // get parameter value from command structure
  SkeletonRegion2D region;

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
      SkeletonSector segment;
      segment.rmin=region.rmin;
      segment.rmax=region.rmax;
      segment.thetamin=region.thetamin;
      segment.thetamax=region.thetamin;
      segment.sector_label=region.label+"_Neumann";
      segment_theta_array1d.push_back(segment);
    }

    {
      SkeletonSector segment;
      segment.rmin=region.rmin;
      segment.rmax=region.rmax;
      segment.thetamin=region.thetamax;
      segment.thetamax=region.thetamax;
      segment.sector_label=region.label+"_Neumann";
      segment_theta_array1d.push_back(segment);
    }
  }

  if(region.rmin > 0.0)
  {
    SkeletonSector segment;
    segment.rmin=region.rmin;
    segment.rmax=region.rmin;
    segment.thetamin=region.thetamin;
    segment.thetamax=region.thetamax;
    segment.sector_label=region.label+"_Neumann";
    segment_r_array1d.push_back(segment);
  }

  {
    SkeletonSector segment;
    segment.rmin=region.rmax;
    segment.rmax=region.rmax;
    segment.thetamin=region.thetamin;
    segment.thetamax=region.thetamax;
    segment.sector_label=region.label+"_Neumann";
    segment_r_array1d.push_back(segment);
  }

  segment_id_map.insert(std::make_pair(region.label+"_Neumann", segment_id_map.size()+1));

  return 0;

}



/* ----------------------------------------------------------------------------
 * get the region index of a cell by its bound
 * if region is overlaped, return the last defined region id
 */
int  MeshGeneratorCylinder2D::get_cell_region_id(const Point &p)
{
  const double pi = 3.14159265358979323846;
  double r = p.size();
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

  int id = -1;
  for(size_t i=0; i < region_array1d.size(); i++)
  {
    if( r >= region_array1d[i].rmin && r <= region_array1d[i].rmax &&
        ((theta >= region_array1d[i].thetamin && theta <= region_array1d[i].thetamax) ||
         (theta+360 >= region_array1d[i].thetamin && theta+360 <= region_array1d[i].thetamax))
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
int MeshGeneratorCylinder2D::set_face(const Parser::Card &c)
{
  const double pi = 3.14159265358979323846;
  std::string  label = c.get_string("label","");

  double rmin = r_points.front();
  double rmax = r_points.back();
  double thetamin=0, thetamax=360;

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


  SkeletonSector segment;
  segment.rmin= find_r(rmin);
  segment.rmax= find_r(rmax);
  segment.thetamin=find_theta(segment.rmin, thetamin);
  segment.thetamax=find_theta(segment.rmin, thetamax);
  segment.sector_label=label;

  if(rmin == rmax)
    segment_r_array1d.push_back(segment);

  if(thetamin == thetamax)
    segment_theta_array1d.push_back(segment);

  segment_id_map.insert(std::make_pair(label, segment_id_map.size() +1));

  return 0;

}

int  MeshGeneratorCylinder2D::get_segment_id_theta(double rmid, const Point &mid)
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

  {
    int id = 0;
    for(unsigned int i=0; i<segment_theta_array1d.size(); ++i)
    {
      const SkeletonSector & segment = segment_theta_array1d[i];
      assert(segment.thetamin == segment.thetamax);
      if( std::abs(theta - segment.thetamin) < 1e-3 && rmid > segment.rmin && rmid < segment.rmax )
        id = segment_id_map.find(segment.sector_label)->second;
    }
    return id;
  }


  return 0;
}




int  MeshGeneratorCylinder2D::get_segment_id_r(const Point &p1, const Point &p2)
{
  const double pi = 3.14159265358979323846;
  double r = 0.5*(p1.size() + p2.size());

  Point mid = 0.5*(p1+p2);
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

  {
    int id = 0;
    for(unsigned int i=0; i<segment_r_array1d.size(); ++i)
    {
      const SkeletonSector & segment = segment_r_array1d[i];
      assert(segment.rmin == segment.rmax);
      if( std::abs(r - segment.rmin) < 1e-10 &&
          ((theta > segment.thetamin && theta < segment.thetamax) ||
           (theta+360 > segment.thetamin && theta+360 < segment.thetamax) )
        )
      {
        id = segment_id_map.find(segment.sector_label)->second;
      }
    }
    return id;
  }

  return 0;
}





/* ----------------------------------------------------------------------------
 * over write the virtual function in meshgen base class
 */
int MeshGeneratorCylinder2D::do_mesh()
{
  START_LOG("do_mesh()", "MeshGenerator");

  MESSAGE<<"Generating mesh by MeshGeneratorCylinder2D...\n"; RECORD();

  double dscale = PhysicalUnit::um;

  // build x-y-z grid skeleton
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if(c.key() == "R.MESH")   // It's a R.MESH card
      if(set_rt_line(c)) return 1;
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
  std::vector< std::vector<Node *> > all_nodes(r_points.size());

  // fill the _mesh structure
  _mesh.magic_num() = this->magic_num();

  for(unsigned int i=0; i<r_points.size(); i++)
  {
    double radiuse = r_points[i];
    int div = theta_points[i];
    

    if(i==0)
    {
      std::vector<Node *> & nodes = all_nodes[i];
      Node * node =  _mesh.add_point(Point());
      nodes.push_back(node);
    }

    if(i==1)
    {
      double radiuse_in = r_points[0];
      double r_mid = 0.5*(radiuse_in+radiuse);

      std::vector<Node *> & nodes = all_nodes[i];

      std::vector<Point> pts = divide_circle(Point(), radiuse, div);
      for(unsigned int n=0; n<pts.size()-1; ++n)
      {
        Node * node =  _mesh.add_point(pts[n]);
        nodes.push_back(node);
      }
      nodes.push_back(nodes.front());

      Node * in_node = all_nodes[0][0];
      for(unsigned int n=0; n<nodes.size()-1; ++n)
      {
        Elem* elem = _mesh.add_elem(Elem::build(TRI3).release());
        elem->set_node(0) = in_node;
        elem->set_node(1) = nodes[n];
        elem->set_node(2) = nodes[n+1];
        elem->subdomain_id() = get_cell_region_id (elem->centroid().unit()*0.5*(radiuse_in+radiuse));

        Point e0_mid = 0.5*(*in_node + *nodes[n]);
        int e0_id = get_segment_id_theta(r_mid, e0_mid);
        if ( e0_id > 0 )
          _mesh.boundary_info->add_side(elem, 0, static_cast<short int>(e0_id) );

        Point e2_mid = 0.5*(*in_node + *nodes[n+1]);
        int e2_id = get_segment_id_theta(r_mid, e2_mid);
        if ( e2_id > 0 )
          _mesh.boundary_info->add_side(elem, 2, static_cast<short int>(e2_id) );

        int e1_id = get_segment_id_r(*nodes[n], *nodes[n+1]);
        if ( e1_id > 0 )
          _mesh.boundary_info->add_side(elem, 1, static_cast<short int>(e1_id) );
      }
    }

    if(i>1)
    {
      double radiuse_in = r_points[i-1];
      double r_mid = 0.5*(radiuse_in+radiuse);

      std::vector<Node *> & nodes = all_nodes[i];

      std::vector<Point> pts = divide_circle(Point(), radiuse, div);
      for(unsigned int n=0; n<pts.size()-1; ++n)
      {
        Node * node =  _mesh.add_point(pts[n]);
        nodes.push_back(node);
      }
      nodes.push_back(nodes.front());

      const std::vector<Node *> & in_nodes = all_nodes[i-1];
      if(in_nodes.size() == nodes.size())
      {
        for(unsigned int n=0; n<in_nodes.size()-1; ++n)
        {
          Elem* elem = _mesh.add_elem(Elem::build(QUAD4).release());

          elem->set_node(0) = in_nodes[n+1];
          elem->set_node(1) = in_nodes[n];
          elem->set_node(2) = nodes[n];
          elem->set_node(3) = nodes[n+1];
          elem->subdomain_id() = get_cell_region_id (elem->centroid().unit()*0.5*(radiuse_in+radiuse));

          int e0_id = get_segment_id_r(*in_nodes[n+1], *in_nodes[n]);
          if ( e0_id > 0 )
            _mesh.boundary_info->add_side(elem, 0, static_cast<short int>(e0_id) );

          int e2_id = get_segment_id_r(*nodes[n], *nodes[n+1]);
          if ( e2_id > 0 )
            _mesh.boundary_info->add_side(elem, 2, static_cast<short int>(e2_id) );

          Point e1_mid = 0.5*(*in_nodes[n] + *nodes[n]);
          int e1_id = get_segment_id_theta(r_mid, e1_mid);
          if ( e1_id > 0 )
            _mesh.boundary_info->add_side(elem, 1, static_cast<short int>(e1_id) );

          Point e3_mid = 0.5*(*in_nodes[n+1] + *nodes[n+1]);
          int e3_id = get_segment_id_theta(r_mid, e3_mid);
          if ( e3_id > 0 )
            _mesh.boundary_info->add_side(elem, 3, static_cast<short int>(e3_id) );
        }
      }

      if((in_nodes.size()-1)*2 == nodes.size()-1)
      {
        for(unsigned int n=0; n<in_nodes.size()-1; ++n)
        {
          {
            Elem* elem1 = _mesh.add_elem(Elem::build(TRI3).release());
            elem1->set_node(0) = in_nodes[n];
            elem1->set_node(1) = nodes[2*n];
            elem1->set_node(2) = nodes[2*n+1];
            elem1->subdomain_id() = get_cell_region_id (elem1->centroid().unit()*0.5*(radiuse_in+radiuse));

            Point e0_mid = 0.5*(*in_nodes[n] + *nodes[2*n]);
            int e0_id = get_segment_id_theta(r_mid, e0_mid);
            if ( e0_id > 0 )
              _mesh.boundary_info->add_side(elem1, 0, static_cast<short int>(e0_id) );

            int e1_id = get_segment_id_r(*nodes[2*n], *nodes[2*n+1]);
            if ( e1_id > 0 )
              _mesh.boundary_info->add_side(elem1, 1, static_cast<short int>(e1_id) );
          }

          {
            Elem* elem2 = _mesh.add_elem(Elem::build(TRI3).release());
            elem2->set_node(0) = in_nodes[n+1];
            elem2->set_node(1) = nodes[2*n+1];
            elem2->set_node(2) = nodes[2*n+2];
            elem2->subdomain_id() = get_cell_region_id (elem2->centroid().unit()*0.5*(radiuse_in+radiuse));

            int e1_id = get_segment_id_r(*nodes[2*n+1], *nodes[2*n+2]);
            if ( e1_id > 0 )
              _mesh.boundary_info->add_side(elem2, 1, static_cast<short int>(e1_id) );

            Point e2_mid = 0.5*(*nodes[2*n+2] + *in_nodes[n+1]);
            int e2_id = get_segment_id_theta(r_mid, e2_mid);
            if ( e2_id > 0 )
              _mesh.boundary_info->add_side(elem2, 2, static_cast<short int>(e2_id) );
          }

          {
            Elem* elem3 = _mesh.add_elem(Elem::build(TRI3).release());
            elem3->set_node(0) = in_nodes[n+1];
            elem3->set_node(1) = in_nodes[n];
            elem3->set_node(2) = nodes[2*n+1];
            elem3->subdomain_id() = get_cell_region_id (elem3->centroid().unit()*0.5*(radiuse_in+radiuse));

            int e0_id = get_segment_id_r(*in_nodes[n], *in_nodes[n+1]);
            if ( e0_id > 0 )
              _mesh.boundary_info->add_side(elem3, 0, static_cast<short int>(e0_id) );
          }
        }
      }
    }

  }

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
  std::map<std::string, int>::const_iterator seg_it = segment_id_map.begin();
  for(; seg_it!=segment_id_map.end(); seg_it++)
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
        bd_index = segment_id_map.size() + bd_map.size() + 1;
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




