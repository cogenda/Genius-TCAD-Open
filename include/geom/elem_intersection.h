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

#ifndef __elem_intersection_h__
#define __elem_intersection_h__

#include <iostream>

#include "point.h"
#include "elem.h"

namespace ElemIntersection
{

  enum IntersectionState
  {
    /**
     * the ray overlap with edge
     */
    Overlap_Edge,

    /**
     * the ray is on face
     */
    On_Face,

    /**
     * the ray intersect with the element
     */
    Intersect_Body,

    /**
     * the ray missed
     */
    Missed
  };


  enum PointLocation
  {
    on_face,
    on_side,
    on_edge,
    on_vertex
  };


  struct Hit_Point
  {
    /**
     * intersection point
     */
    Point p;

    /**
     * ray parameter
     */
    double t;

    /**
     * indicate where the hit point is
     */
    PointLocation point_location;

    /**
     * which face/edge/vertex the intersection point on
     */
    unsigned int mark;
  };


  struct IntersectionResult
  {
    /**
     * the ray-elem intersection state
     */
    IntersectionState  state;

    /**
     * when IntersectionState is Overlap_Edge/On_Face, record the edge/face index
     */
    unsigned int mark;

    /**
     * record all the intersection points
     */
    std::vector<Hit_Point> hit_points;

    bool is_mark_exist(PointLocation loc, unsigned int mark) const
    {
      for(unsigned int n=0; n<hit_points.size(); ++n)
      {
        if(loc==hit_points[n].point_location && hit_points[n].mark==mark) return true;
      }
      return false;
    }

    void print()
    {
      switch(state)
      {
      case Overlap_Edge       : std::cout<<"Overlap_Edge" <<std::endl; break;
      case On_Face            : std::cout<<"On_Face" <<std::endl; break;
      case Intersect_Body     : std::cout<<"Intersect_Body" <<std::endl; break;
      case Missed             : std::cout<<"Missed" <<std::endl; return;
      }

      std::cout<<"Hit_Point="<<hit_points.size()<<std::endl;
      for(unsigned int n=0; n<hit_points.size(); ++n)
      {
        std::cout<<hit_points[n].p;
        switch(hit_points[n].point_location)
        {
        case on_face   : std::cout<<"on_face   " << hit_points[n].mark<<"  t="<<hit_points[n].t<< std::endl; break;
        case on_side   : std::cout<<"on_side   " << hit_points[n].mark<<"  t="<<hit_points[n].t<< std::endl; break;
        case on_edge   : std::cout<<"on_edge   " << hit_points[n].mark<<"  t="<<hit_points[n].t<< std::endl; break;
        case on_vertex : std::cout<<"on_vertex " << hit_points[n].mark<<"  t="<<hit_points[n].t<< std::endl; break;
        }
      }
    }
  };
}

using namespace ElemIntersection;

#endif

