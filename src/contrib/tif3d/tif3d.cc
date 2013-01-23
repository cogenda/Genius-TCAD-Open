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
#include <iostream>
#include <fstream>
#include <sstream>

#include "tif3d.h"


//stuff for the TIF3D read write


TIF3D::TIF3D(const std::string & file) : _file(file)
{

}


int TIF3D::read()
{
  std::ifstream ctmp(_file.c_str(), std::ios::in);

  if (!ctmp.good())
  {
    std::cerr<<"Open TIF3D file error."<<std::endl;
    return 1;
  }


  char flag;
  while(ctmp >> flag)
  {

    if (flag == 'H' || flag == 'h')
    {
      // skip tif file header
      std::string buf;
      std::getline(ctmp, buf);

      if( buf.find("V1.1") == std::string::npos && buf.find("V1.2") == std::string::npos )
      {
        std::cerr<<"TIF3D should have version >= 1.1"<<std::endl;
        return 1;
      }
    }


    // node coordinate
    else if (flag == 'C' || flag == 'c')
    {
      Node_t node;
      ctmp >> node.index >> node.x >> node.y >> node.z;
      // save it
      _nodes.push_back(node);
    }

    // face
    else if (flag == 'F' || flag == 'f')
    {
      std::string buf;
      std::getline(ctmp, buf);
      std::stringstream ss(buf);

      Face_t f;
      ss >> f.index >> f.point1 >> f.point2 >> f.point3 >> f.bc_index ;
      // save it
      _faces.push_back(f);
    }


    // tet
    else if (flag == 'T' || flag == 't')
    {
      std::string buf;
      std::getline(ctmp, buf);
      std::stringstream ss(buf);

      Tet_t t;
      ss >> t.index >> t.region >> t.c1 >> t.c2 >> t.c3 >> t.c4;
      // save it
      _tets.push_back(t);
    }


    //region
    else if (flag == 'R' || flag == 'r')
    {
      std::string buf;
      std::getline(ctmp, buf);
      std::stringstream ss(buf);

      Region_t region;
      ss >> region.index >> region.material >> region.name;
      region.node_num  = 0;
      region.tet_num   = 0;
      // save it
      _regions.push_back(region);
    }

    //face label
    else if (flag == 'I' || flag == 'i')
    {
      int index, bc_index;
      std::string type, name;
      ctmp >> index >> name >> bc_index;
      // save it
      _face_labels.insert(std::make_pair(bc_index, name));
    }

    // solutions
    else if (flag == 'S' || flag == 's')
    {
      ctmp >> _sol_head.sol_num;
      for(int i = 0; i < _sol_head.sol_num; i++)
      {
        std::string sol_name;
        ctmp >> sol_name;
        _sol_head.sol_name_array.push_back(sol_name);
      }
    }

     // solution units
    else if (flag == 'U' || flag == 'u')
    {
      int unit_num;
      ctmp >> unit_num;
      for(int i = 0; i < unit_num; i++)
      {
        std::string sol_unit;
        ctmp >> sol_unit;
        _sol_head.sol_unit_array.push_back(sol_unit);
      }
    }

    // solution data
    else if (flag == 'N' || flag == 'n')
    {
      SolData_t solution;
      ctmp >> solution.index >> solution.region_index;
      //For all data values...
      for(int i = 0; i < _sol_head.sol_num; i++)
      {
        double dval;
        ctmp >> dval;
        solution.data_array.push_back(dval);
      }
      _sol_data.push_back(solution);
    }

    else
    {
      std::string rubbish;
      std::getline(ctmp, rubbish);
    }
  }

  ctmp.close();

  //statistic how many triangles in each region
  for(unsigned int n=0; n<_tets.size(); ++n)
  {
    assert(_tets[n].region < _regions.size());
    _regions[_tets[n].region].tet_num++;
  }

  // check if empty region exist
  for(unsigned int n=0; n<_regions.size(); ++n)
    assert(_regions[n].tet_num > 0);

  return 0;
}

