#if 0

#include <cassert>
#include <cmath>

#include "genius_common.h"
#include "asinh.hpp"
#include "interpolation_3d_qshep.h"
#include "parallel.h"

#ifdef WINDOWS
  #include "qshep_win.h"
#else
  #include "qshep_linux.h"
#endif

#include "log.h"

Interpolation3D_qshep::Interpolation3D_qshep()
{}


Interpolation3D_qshep::~Interpolation3D_qshep()
{
  clear();
}


void Interpolation3D_qshep::clear()
{
  field.clear();
  field_limit.clear();
}



void Interpolation3D_qshep::add_scatter_data(const Point & point, int group, double value)
{
  field[group].x.push_back(point.x());
  field[group].y.push_back(point.y());
  field[group].z.push_back(point.z());

  if(field_limit.find(group)==field_limit.end())
    field_limit[group] = std::make_pair(value, value);
  else
  {
    field_limit[group].first  = std::min(field_limit[group].first, value);
    field_limit[group].second = std::max(field_limit[group].second, value);
  }

  InterpolationType type = _interpolation_type[group];
  switch(type)
  {
  case Linear    : break;
  case SignedLog : value = (value>0 ? 1.0 : -1.0)*std::log(1.0+std::abs(value)); break;
  case Asinh     : value = asinh(value); break;
  }

  field[group].f.push_back(value);
  //std::cout<<field[group].f.size()<<std::endl;
}



void Interpolation3D_qshep::setup(int group)
{
  QSHEP & qshep = field[group];

  qshep.n  = qshep.f.size();
  qshep.nq = 17;
  qshep.nw = 32;
  qshep.nr = int(std::pow(qshep.n/3.0, 1.0/3.0));
  qshep.lcell.resize(qshep.nr*qshep.nr*qshep.nr);
  qshep.lnext.resize(qshep.n);
  qshep.rsq.resize(qshep.n);
  qshep.A.resize(9*qshep.n);

  QSHEP3 ( qshep.n, &(qshep.x[0]), &(qshep.y[0]), &(qshep.z[0]), &(qshep.f[0]),
           qshep.nq, qshep.nw, qshep.nr, &(qshep.lcell[0]), &(qshep.lnext[0]), &(qshep.xyzmin[0]), &(qshep.xyzdel[0]),
           qshep.rmax, &(qshep.rsq[0]), &(qshep.A[0]), qshep.ier);

  if(qshep.ier)
  {
     MESSAGE<<"\nERROR: 3D interpolation failed."<<std::endl; RECORD();
     switch(qshep.ier)
     {
       case 1 : MESSAGE<<"N, NQ, NW, or NR is out of range."<<std::endl; RECORD(); break;
       case 2 : MESSAGE<<"Duplicate nodes were encountered."<<std::endl; RECORD(); break;
       case 3 : MESSAGE<<"No unique solution due to collinear nodes."<<std::endl; RECORD(); break;
     }
     genius_error();
  }
}


void Interpolation3D_qshep::broadcast(unsigned int root)
{
  std::vector<int> groups;
  std::map<int, QSHEP>::iterator it=field.begin();
  for(; it!=field.end(); ++it)
    groups.push_back(it->first);
  Parallel::broadcast(groups, root);

  for(unsigned int n=0; n<groups.size(); ++n)
  {
     QSHEP & qshep = field[groups[n]];
     Parallel::broadcast(qshep.x, root);
     Parallel::broadcast(qshep.y, root);
     Parallel::broadcast(qshep.z, root);
     Parallel::broadcast(qshep.f, root);

     std::pair<double, double> & limits = field_limit[groups[n]];
     Parallel::broadcast(limits.first,  root);
     Parallel::broadcast(limits.second, root);
  }
}


double Interpolation3D_qshep::get_interpolated_value(const Point & point, int group) const
{
  QSHEP & qshep = const_cast<QSHEP &>(field.find(group)->second);

  double px = point[0];
  double py = point[1];
  double pz = point[2];
  double value = QS3VAL (px, py, pz, qshep.n, &(qshep.x[0]), &(qshep.y[0]), &(qshep.z[0]), &(qshep.f[0]),
                         qshep.nr, &(qshep.lcell[0]), &(qshep.lnext[0]), &(qshep.xyzmin[0]), &(qshep.xyzdel[0]),
                         qshep.rmax, &(qshep.rsq[0]), &(qshep.A[0]));

  InterpolationType type = _interpolation_type.find(group)->second;

  switch(type)
  {
  case Linear    : break;
  case SignedLog : value = value>0 ? exp(value)-1 : 1-exp(-value) ; break;
  case Asinh     : value = sinh(value); break;
  }

  double vmin = field_limit.find(group)->second.first;
  double vmax = field_limit.find(group)->second.second;
  if(value<vmin) value = vmin;
  if(value>vmax) value = vmax;

  return value;
}


#endif

