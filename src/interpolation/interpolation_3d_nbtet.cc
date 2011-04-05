/*
 * interpolation_3d_nbtet.cc
 *
 *  Created on: May 18, 2009
 *      Author: hash
 */

#include <cassert>
#include <cmath>

#include "genius_common.h"
#include "asinh.hpp"
#include "interpolation_3d_qshep.h"
#include "parallel.h"

#include "ANN/ANN.h"
#include "interpolation_3d_nbtet.h"

#include "log.h"

ANNSession::ANNSession(int dim) :
  _dim(dim), _tree_built(false), _datapts(NULL), _databuf(NULL), _kdTree(NULL)
{ }


void ANNSession::clear()
{
  if (_datapts!=NULL)
    delete _datapts;
  if (_databuf!=NULL)
    delete _databuf;
  if (_kdTree!=NULL)
    delete _kdTree;

  _databuf = NULL;
  _datapts = NULL;
  _kdTree  = NULL;
  _tree_built = false;
}

ANNSession::~ANNSession()
{
  this->clear();
}

bool ANNSession::addPoint(const Point &pt)
{
  if (_tree_built)
    return false;

  _pts.push_back(pt);

  return true;
}

void ANNSession::broadcast(unsigned int root)
{
  // data buffer must not have been setup
  assert(_databuf==NULL);
  // kd-tree not built yet.
  assert(_tree_built==false);

  unsigned int size = _pts.size();
  Parallel::broadcast(size, root);

  std::vector<double> x,y,z;
  x.resize(size);
  y.resize(size);
  z.resize(size);

  if(Genius::processor_id()==root)
  {
    for (unsigned int i=0; i<size; i++){
      x[i] = _pts[i](0);
      y[i] = _pts[i](1);
      z[i] = _pts[i](2);
    }
  }
  else
  {
    _pts.resize(size);
  }
  Parallel::broadcast(x, root);
  Parallel::broadcast(y, root);
  Parallel::broadcast(z, root);
  if(Genius::processor_id()!=root)
  {
    for (unsigned int i=0; i<size; i++){
      _pts[i](0) = x[i];
      _pts[i](1) = y[i];
      _pts[i](2) = z[i];
    }
  }
}

void ANNSession::setup()
{
  assert(_tree_built==false);
  assert(_kdTree==NULL);

  _databuf = new ANNcoord[_dim*_pts.size()];
  _datapts = new ANNpoint[_pts.size()];

  for (unsigned int i = 0; i < _pts.size(); i++) {
    _databuf[i*_dim+0] = _pts[i](0);
    _databuf[i*_dim+1] = _pts[i](1);
    if (_dim==3)
      _databuf[i*_dim+2] = _pts[i](2);

    _datapts[i] = &(_databuf[i*_dim]);
  }
  _kdTree = new ANNkd_tree(_datapts,_pts.size(),_dim);
  _tree_built = true;
}

ANNResultType ANNSession::search(const Point &pt, const unsigned int k) const
{
  assert(_tree_built);
  assert(_kdTree);

  ANNpoint ann_pt = new ANNcoord[_dim];  // allocate search point
  ann_pt[0] = pt(0);
  ann_pt[1] = pt(1);
  if (_dim==3)
    ann_pt[2] = pt(2);

  ANNResultType res;
  ANNidxArray  _residx;       // indices of result points
  ANNdistArray _resdist;      // distances of result points

  if (_pts.size()<k)
    return res;

  _residx  = new ANNidx[k];            // allocate near neighbor indices
  _resdist = new ANNdist[k];           // allocate near neighbor dists
  _kdTree->annkSearch(ann_pt, k, _residx, _resdist, 0);

  for(unsigned int i=0; i<k; i++)
  {
    res.push_back(std::pair<int, double>(_residx[i],_resdist[i]));
  }

  delete _residx;
  delete _resdist;

  return res;
}

Point ANNSession::getPointCoord(unsigned int i) const
{
  if (i<_pts.size())
    return _pts[i];
  else
    return Point();
}

long ANNSession::size() const
{
  return _pts.size();
}

Interpolation3D_nbtet::Interpolation3D_nbtet() :
  _ann(3)
{
}

void Interpolation3D_nbtet::clear()
{
  _field.clear();
  _ann.clear();
}

void Interpolation3D_nbtet::broadcast(unsigned int root)
{
  _ann.broadcast(root);
  Parallel::broadcast(_field, root);
}

void Interpolation3D_nbtet::add_scatter_data(const Point & pt, int group, double value)
{
  assert(group==0); // we currently support 1 group only
  _ann.addPoint(pt);

  InterpolationType type = _interpolation_type[group];
  _field.push_back(scaleValue(type, value));
}

void Interpolation3D_nbtet::setup(int group)
{
  assert(group==0);
  _ann.setup();
}

double Interpolation3D_nbtet::get_interpolated_value(const Point & pt, int group) const
{
  assert(group==0); // we currently support 1 group only

  double toler;
  unsigned int maxpt = 20;
  int ia, ib, ic, id;
  Point tetp[4];
  ANNResultType ann_res = _ann.search(pt,maxpt);

  assert(ann_res.size()==maxpt);

  ia = ann_res[0].first;
  ib = ann_res[1].first;
  tetp[0] = _ann.getPointCoord(ia); // coord of 1st point
  tetp[1] = _ann.getPointCoord(ib); // coord of 2nd point

  InterpolationType type = _interpolation_type.find(group)->second;

  for (int pass=0; pass<2; pass++)
  {
    for (unsigned int i=2; i<maxpt; i++)
    {
      // the 3rd point must NOT fall on the same line
      unsigned int j;
      for (j=i; j<maxpt; j++)
      {
        Point vab, vac, r, vtmp;

        ic = ann_res[j].first;
        tetp[2] = _ann.getPointCoord(ic); // coord of 3rd point
        vab = tetp[0]-tetp[1];
        vac = tetp[0]-tetp[2];
        r = vab.cross(vac);

        vtmp = tetp[1]-tetp[2];
        toler = vab.size();
        toler = vac.size()<toler ? vac.size() : toler;
        toler = vtmp.size()<toler ? vtmp.size() : toler;

        if (r.size()>toler*toler/100)
          break; // area > 0, acceptable point!
      }
      j++;
      // the 4th point must NOT fall on the same plane
      double vol, orient;
      for (; j<maxpt; j++)
      {
        Point vab, vac, vad, vtmp;

        id = ann_res[j].first;
        tetp[3] = _ann.getPointCoord(id); // coord of 4th point
        vab = tetp[1]-tetp[0];
        vac = tetp[2]-tetp[0];
        vad = tetp[3]-tetp[0];
        vol = vab.dot(vad.cross(vac));
        orient = vol>0 ? 1.0 : -1.0;
        vol = std::abs(vol);

        toler = vab.size();
        toler = vac.size()<toler ? vac.size() : toler;
        toler = vad.size()<toler ? vac.size() : toler;
        vtmp = tetp[1]-tetp[2];
        toler = vtmp.size()<toler ? vtmp.size() : toler;
        vtmp = tetp[1]-tetp[3];
        toler = vtmp.size()<toler ? vtmp.size() : toler;
        vtmp = tetp[2]-tetp[3];
        toler = vtmp.size()<toler ? vtmp.size() : toler;

        if ( vol > toler*toler*toler/1000 )
          break; // volume > 0, acceptable 4th point.
      }
      if (j<maxpt)
      {
        // we have a tetrahedron now, check if it contains the point
        double vol_a, vol_b, vol_c, vol_d;
        Point vpa, vpb, vpc, vpd;
        vpa = tetp[0] - pt;
        vpb = tetp[1] - pt;
        vpc = tetp[2] - pt;
        vpd = tetp[3] - pt;
        vol_a = orient*(vpb.dot(vpd.cross(vpc)));
        vol_b = orient*(vpa.dot(vpc.cross(vpd)));
        vol_c = orient*(vpa.dot(vpd.cross(vpb)));
        vol_d = orient*(vpa.dot(vpb.cross(vpc)));

        if (pass==0)
        {
          // in the first pass, we require strictly that the point is in the tetrahedron
          if (vol_a>=0 && vol_b>=0 && vol_c>=0 && vol_d>=0)
          {
            return unscaleValue(type, (vol_a*_field[ia] + vol_b*_field[ib] + vol_c*_field[ic] + vol_d*_field[id])/vol);
          }
        }
        else
        {
          // in the second pass, we allow slight extroplation
          double extrap_tol = -1e-2*vol;
          if (vol_a>extrap_tol && vol_b>extrap_tol && vol_c>extrap_tol && vol_d>extrap_tol)
          {
            return unscaleValue(type, (vol_a*_field[ia] + vol_b*_field[ib] + vol_c*_field[ic] + vol_d*_field[id])/vol);
          }
        }
      }
    } // end of loop i
  } // end of loop pass
  if((pt-tetp[0]).size()<1e10) // FIXME: should have a more reasonable threshold here.
  {
    //MESSAGE << "Interpolation: warning: using nearest data point at line " << ia << std::endl; RECORD();
    return unscaleValue(type, _field[ia]);
  }
  else
  {
    MESSAGE << "Interpolation: warning: interpolation failed, assume zero at this point." << std::endl; RECORD();
    return 0.0;
  }
}





