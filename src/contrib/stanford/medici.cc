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

#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include "medici.h"
#include "material_define.h"


MediciTIF::MediciTIF(const std::string & file)
  :StanfordTIF(file)
{}



bool MediciTIF::read()
{
  std::ifstream ctmp(_file.c_str(), std::ios::in);

  if (!ctmp.good())
  {
    std::cerr<<"Open Medici TIF file error."<<std::endl;
    return false;
  }

  std::string buffer;
  std::string flag;

  while(ctmp >> flag)
  {
    // the version of suprem file
    if (flag == "h")
    {
      // skip file version info
      std::getline(ctmp, buffer);
      continue;
    }

    // node coordinate
    else if (flag == "c")
    {
      Node_t node;
      ctmp >> node.index >> node.x >> node.y >> node.h;
      // correct subscript
      node.index = node.index - 1;
      // save it
      _nodes.push_back(node);
      continue;
    }

    else if (flag == "e")
    {
      Edge_t edge;
      ctmp >> edge.index >> edge.point1 >> edge.point2 >> edge.bcode;
      // correct subscript
      edge.index  = edge.index  -1;
      edge.point1 = edge.point1 -1;
      edge.point2 = edge.point2 -1;
      edge.bc_index = -1024;
      if(edge.point1 > edge.point2) std::swap(edge.point1, edge.point2);
      // save it
      _edges.push_back(edge);
      continue;
    }

    // triangle
    else if (flag == "t")
    {
      Tri_t t;
      ctmp >> t.index >> t.region >> t.c1 >> t.c2 >> t.c3 >> t.t1 >> t.t2 >> t.t3;
      // correct subscript
      t.index = t.index - 1;
      t.region = t.region - 1;
      t.c1 = t.c1 -1;
      t.c2 = t.c2 -1;
      t.c3 = t.c3 -1;
      t.t1 = t.t1 -1;
      t.t2 = t.t2 -1;
      t.t3 = t.t3 -1;
      // save it
      _tris.push_back(t);
      continue;
    }


    //region
    else if (flag == "r")
    {
      Region_t region;
      int material_id;
      ctmp >> region.index >> region.material >> region.name;
      region.index = region.index - 1;
      region.node_num  = 0;
      region.tri_num   = 0;
      region.segment = false;
      //add "r" before region name if it begin with number
      if(isdigit(region.name[0]))
        region.name = std::string("r") + region.name;
      // save it
      _regions.push_back(region);
      continue;
    }
    // region boundary
    else if (flag == "b")
    {
      // last region
      Region_t & region = _regions[_regions.size()-1];
      int boundary_index;
      ctmp >> boundary_index;
      region.boundary.push_back(boundary_index-1);
      continue;
    }

    //interface
    else if (flag == "i")
    {
      Region_t region;
      int material_id;
      ctmp >> region.index >> region.material >> region.name >> region.region;
      region.index = region.index - 1;
      region.node_num  = 0;
      region.tri_num   = 0;
      region.region = region.region-1;
      region.segment = true;
      //add "i" before region name if it begin with number
      if(isdigit(region.name[0]))
        region.name = std::string("i") + region.name;
      // save it
      _regions.push_back(region);
      continue;
    }
    // interface boundary
    else if (flag == "j")
    {
      // last region
      Region_t & region = _regions[_regions.size()-1];
      int boundary_index;
      ctmp >> boundary_index;
      region.boundary.push_back(boundary_index-1);
      continue;
    }

    // solutions
    else if (flag == "s")
    {
      ctmp >> _sol_head.sol_num;
      for(int i = 0; i < _sol_head.sol_num; i++)
      {
        std::string sol_name;
        ctmp >> sol_name;
        _sol_head.sol_name_array.push_back(sol_name);
      }
      continue;
    }

    // solution data
    else if (flag == "n")
    {
      SolData_t solution;
      ctmp >> solution.index >> solution.material;
      solution.index = solution.index - 1;
      solution.region_index=-1;
      //For all data values...
      for(int i = 0; i < _sol_head.sol_num; i++)
      {
        double dval;
        ctmp >> dval;
        solution.data_array.push_back(dval);
      }
      _sol_data.push_back(solution);
      continue;
    }

    else
      std::getline(ctmp, buffer);
  }


  ctmp.close();

  if(_nodes.size()==0 || _tris.size()==0 ) return false;

  //statistic how many triangles in each region
  for(unsigned int n=0; n<_tris.size(); ++n)
    _regions[_tris[n].region].tri_num++;

  _find_solution_region_by_material();

  _acceptor_index = _sol_head.solution_index("Accept");
  _donor_index    = _sol_head.solution_index("Donor");
  _mole_x_index   = _sol_head.solution_index("mole_x");
  _mole_y_index   = _sol_head.solution_index("mole_y");

  return true;

}




