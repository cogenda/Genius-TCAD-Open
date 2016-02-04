#include "doping_fun.h"

#include <cmath>
//win32 does not have erfc function
#ifdef WINDOWS
#include "mathfunc.h"
#endif
#include "physical_unit.h"

using PhysicalUnit::um;

//------------------------------------------------------------------

double UniformDopingFunction::profile(double x,double y,double z)
{
  if( x >= _xmin-1e-6*um && x <= _xmax+1e-6*um &&
      y >= _ymin-1e-6*um && y <= _ymax+1e-6*um &&
      z >= _zmin-1e-6*um && z <= _zmax+1e-6*um )
    return _ion*_peak;
  else
    return 0.0;
}


//------------------------------------------------------------------


double LinearDopingFunction::profile(double x,double y,double z)
{
    if( x >= _xmin-1e-6*um && x <= _xmax+1e-6*um &&
        y >= _ymin-1e-6*um && y <= _ymax+1e-6*um &&
        z >= _zmin-1e-6*um && z <= _zmax+1e-6*um )
    {
      switch(_dir)
      {
        case X_Direction:
         if(x >= _base_plane && x <= _end_plane)
           return  _ion*(_doping_begin + (x-_base_plane)*_doping_slope);
         return 0.0;
        case Y_Direction:
         if(y >= _base_plane && y <= _end_plane)
           return  _ion*(_doping_begin + (y-_base_plane)*_doping_slope);
         return 0.0;
        case Z_Direction:
         if(z >= _base_plane && z <= _end_plane)
           return  _ion*(_doping_begin + (z-_base_plane)*_doping_slope);
         return 0.0;
        default: genius_error();
      }
    }
    
    return 0;
  
}


//------------------------------------------------------------------

/**
 * compute doping concentration by point location
 */
double AnalyticDopingFunction::profile(double x, double y, double z)
{
  double dx, dy, dz;
  if ( _XERFC )
#ifdef WINDOWS
    dx = (Erfc((x-_xmax)/_XCHAR)-Erfc((x-_xmin)/_XCHAR))/2.0;
#else
    dx = (erfc((x-_xmax)/_XCHAR)-erfc((x-_xmin)/_XCHAR))/2.0;
#endif
  else
  {
    if(x<_xmin)
      dx = exp(-(x-_xmin)*(x-_xmin)/(_XCHAR*_XCHAR));
    else if(x>=_xmin&&x<=_xmax)
      dx = 1.0;
    else
      dx = exp(-(x-_xmax)*(x-_xmax)/(_XCHAR*_XCHAR));
  }

  if ( _YERFC )
#ifdef WINDOWS
    dy = (Erfc((y-_ymax)/_YCHAR)-Erfc((y-_ymin)/_YCHAR))/2.0;
#else
    dy = (erfc((y-_ymax)/_YCHAR)-erfc((y-_ymin)/_YCHAR))/2.0;
#endif
  else
  {
    if(y<_ymin)
      dy = exp(-(y-_ymin)*(y-_ymin)/(_YCHAR*_YCHAR));
    else if(y>=_ymin&&y<=_ymax)
      dy = 1.0;
    else
      dy = exp(-(y-_ymax)*(y-_ymax)/(_YCHAR*_YCHAR));
  }

  if ( _ZERFC )
#ifdef WINDOWS
    dz = (Erfc((z-_zmax)/_ZCHAR)-Erfc((z-_zmin)/_ZCHAR))/2.0;
#else
    dz = (erfc((z-_zmax)/_ZCHAR)-erfc((z-_zmin)/_ZCHAR))/2.0;
#endif
  else
  {
    if(z<_zmin)
      dz = exp(-(z-_zmin)*(z-_zmin)/(_ZCHAR*_ZCHAR));
    else if(z>=_zmin&&z<=_zmax)
      dz = 1.0;
    else
      dz = exp(-(z-_zmax)*(z-_zmax)/(_ZCHAR*_ZCHAR));
  }

  return _ion*_peak*dx*dy*dz;
}


//------------------------------------------------------------------

PolyMaskDopingFunction::PolyMaskDopingFunction(double ion, const std::vector<Point> &poly, double theta, double phi,
                       double rmin, double rmax, double npeak, double char_depth, double char_lateral,
                       double resolution_factor)
  : DopingFunction(ion),
    _rmin(rmin), _rmax(rmax), _peak(npeak), _char_depth(char_depth), _char_lateral(char_lateral), _resolution_factor(resolution_factor),
    _mask_mesh(poly, char_lateral/resolution_factor)
{
    // doping line direction
  double deg = 3.14159265359/180.0;
  _dir = Point( sin(theta*deg)*cos(phi*deg), cos(theta*deg), sin(theta*deg)*sin(phi*deg));
}



double PolyMaskDopingFunction::profile(double x, double y, double z)
{
  Point p(x, y, z); // doping point
  Plane mask_plane = _mask_mesh.plane(); // mask plane
  Point point_on_plane;
  {
    Point  plane_point = mask_plane.point();
    Point  plane_norm = mask_plane.normal();

    double t = (plane_point - p)*plane_norm/(plane_norm*_dir);
    point_on_plane = p + _dir*t;
  }

  double r = (p-point_on_plane)*_dir;
  if( r > _rmax + 5*_char_lateral || r < _rmin - 5*_char_lateral ) return 0.0;

  double profile = 0.0;
  const double A = prof_func_lnorm();
  std::vector<double> vecr, vecl;   // array of range/lateral coordinates

  // if point outside the doping box, skip
  typedef std::pair<unsigned int, unsigned int> Idx;
  std::pair<Idx, Idx> range;
  if( !_mask_mesh.index_range(point_on_plane, 5*_char_lateral, range) ) return 0.0;

    // or we have to integral all the doping lines
  for(unsigned int i=range.first.first; i<=range.second.first; ++i)
    for(unsigned int j=range.first.second; j<=range.second.second; ++j)
  {
    Idx idx = std::make_pair(i,j);
    if( !_mask_mesh.inside(idx) ) continue;

    Point loc = _mask_mesh.location(idx);
    Point project_p = loc + _dir*((p-loc)*_dir);

    double r = (project_p - loc)*_dir;
    double dist = (project_p - p).size();

    vecr.push_back(r);
    vecl.push_back(dist);
  }

  std::vector<double> vecfr = prof_func_r(vecr);
  std::vector<double> vecfl = prof_func_l(vecl);

  for (size_t i=0; i<vecfr.size(); ++i)
  {
    profile += A*vecfr[i]*vecfl[i];
  }

  return _ion*profile;
}


std::vector<double> PolyMaskDopingFunction::prof_func_r(const std::vector<double> &rs) const
{
  std::vector<double> res;
  for (std::vector<double>::const_iterator i=rs.begin(); i!=rs.end(); i++)
  {
    const double r = *i;

    double dr;
    if(r<_rmin)
      dr = exp(-(r-_rmin)*(r-_rmin)/(_char_depth*_char_depth));
    else if(r<=_rmax)
      dr = 1.0;
    else
      dr = exp(-(r-_rmax)*(r-_rmax)/(_char_depth*_char_depth));

    res.push_back(_peak*dr);
  }
  return res;
}

std::vector<double> PolyMaskDopingFunction::prof_func_l(const std::vector<double> &dists) const
{
  std::vector<double> res;
  for (std::vector<double>::const_iterator i=dists.begin(); i!=dists.end(); i++)
  {
    const double dist = *i;
    if (dist==0.0)
      res.push_back(1.0);
    else
      res.push_back( exp( -dist*dist/(_char_lateral*_char_lateral) ) );
  }
  return res;
}

double PolyMaskDopingFunction::prof_func_lnorm() const
{
  double density_of_doping_line = (_char_lateral/_resolution_factor);
  return (density_of_doping_line*density_of_doping_line)/(3.1415927*_char_lateral*_char_lateral);
}

//------------------------------------------------------------------


void RecMaskDopingFunction::_init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double theta, double phi)
{
  // mask rectangle
  _mask_min = Point(xmin,ymin,zmin);
  _mask_max = Point(xmax,ymax,zmax);

  // mask plane
  if(xmin == xmax)
    _mask.yz_plane(xmin);

  if(ymin == ymax)
    _mask.xz_plane(ymin);

  if(zmin == zmax)
    _mask.xy_plane(zmin);

  // doping line direction
  double deg = 3.14159265359/180.0;
  _dir = Point( sin(theta*deg)*cos(phi*deg), cos(theta*deg), sin(theta*deg)*sin(phi*deg));

  //doping bounding box
  double doping_xmin = _mask.point().x() + _rmin*_dir.x() - 5*_char_depth;
  double doping_xmax = _mask.point().x() + _rmax*_dir.x() + 5*_char_depth;
  double doping_ymin = _mask.point().y() + _rmin*_dir.y() - 5*_char_depth;
  double doping_ymax = _mask.point().y() + _rmax*_dir.y() + 5*_char_depth;
  double doping_zmin = _mask.point().z() + _rmin*_dir.z() - 5*_char_depth;
  double doping_zmax = _mask.point().z() + _rmax*_dir.z() + 5*_char_depth;
  _doping_min = Point(doping_xmin,doping_ymin,doping_zmin);
  _doping_max = Point(doping_xmax,doping_ymax,doping_zmax);
}



double RecMaskDopingFunction::profile(double x, double y, double z)
{
  Point p(x,y,z);

  if( _mask.is_xy_plane() )
    if( p.z() > _doping_max.z() || p.z() < _doping_min.z() ) return 0.0;
  if( _mask.is_xz_plane() )
    if( p.y() > _doping_max.y() || p.y() < _doping_min.y() ) return 0.0;
  if( _mask.is_yz_plane() )
    if( p.x() > _doping_max.x() || p.x() < _doping_min.x() ) return 0.0;

  Point point_on_plane = p + _dir*((_mask.point() - p)*_mask.normal()/(_mask.normal()*_dir));

  Point range_max = _mask.closest_point(point_on_plane + Point(5*_char_depth, 5*_char_depth, 5*_char_depth));
  Point range_min = _mask.closest_point(point_on_plane - Point(5*_char_depth, 5*_char_depth, 5*_char_depth));

  if(range_min > _mask_max || range_max < _mask_min) return 0.0;

  Point hot_area_min(std::max(range_min.x(), _mask_min.x()), std::max(range_min.y(), _mask_min.y()),std::max(range_min.z(), _mask_min.z()));
  Point hot_area_max(std::min(range_max.x(), _mask_max.x()), std::min(range_max.y(), _mask_max.y()),std::min(range_max.z(), _mask_max.z()));

  double profile = 0.0;
  double sample_distance = _char_depth/_resolution_factor;
  const double A = (sample_distance*sample_distance)/(3.1415927*_char_lateral*_char_lateral);
  for(double sx = hot_area_min.x(); sx<=hot_area_max.x(); sx+=sample_distance)
    for(double sy = hot_area_min.y(); sy<=hot_area_max.y(); sy+=sample_distance)
      for(double sz = hot_area_min.z(); sz<=hot_area_max.z(); sz+=sample_distance)
      {
        Point loc_start(sx,sy,sz);
        Point project_p = loc_start + _dir*((p-loc_start)*_dir);

        double r = (project_p - loc_start)*_dir;
        double dist = (project_p - p).size();
        profile += _profile_r(r)*A*exp( -dist*dist/(_char_lateral*_char_lateral) ); // NOTE: this line is time consuming
      }
  return _ion*profile;
}


double RecMaskDopingFunction::_profile_r(double r) const
{
  double dr;
  if(r<_rmin)
    dr = exp(-(r-_rmin)*(r-_rmin)/(_char_depth*_char_depth));
  else if(r<=_rmax)
    dr = 1.0;
  else
    dr = exp(-(r-_rmax)*(r-_rmax)/(_char_depth*_char_depth));
  return _peak*dr;
}

