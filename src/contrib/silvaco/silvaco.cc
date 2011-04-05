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

#include <iostream>
#include <fstream>
#include <sstream>

#include "silvaco.h"


//stuff for the Silvaco read write
SilvacoTIF::SilvacoTIF(const std::string & file) : _file(file)
{
  // init SilImp
  SilImp[   0] = "Vacancy";
  SilImp[   1] = "Interstitial";
  SilImp[   2] = "Arsenic";
  SilImp[   3] = "Phosphorus";
  SilImp[   4] = "Antimony";
  SilImp[   5] = "Boron";
  SilImp[   6] = "Electrons";
  SilImp[   7] = "Holes";
  SilImp[  10] = "Oxidant";
  SilImp[  11] = "Oxidant";
  SilImp[  12] = "Traps";
  SilImp[  14] = "Potential";
  SilImp[  25] = "Doping";
  SilImp[  71] = "Donor";
  SilImp[  72] = "Acceptor";
  SilImp[ 151] = "BoronActive";
  SilImp[ 152] = "PhosphorusActive";
  SilImp[ 153] = "ArsenicActive";
  SilImp[ 154] = "AntimonyActive";

  // init SilMat
  SilMat[   0] = "Gas";
  SilMat[   1] = "Oxide";
  SilMat[   2] = "Nitride";
  SilMat[   3] = "Silicon";
  SilMat[   4] = "Poly";
  SilMat[   5] = "Oxynitride";
  SilMat[   6] = "Aluminum";
  SilMat[   7] = "Photores";
  SilMat[   8] = "GaAs";
  SilMat[  14] = "Titanium";
  SilMat[  28] = "Gas";
  SilMat[  32] = "Ge";
  SilMat[  81] = "TiSi2";
  SilMat[  91] = "Elec";
  SilMat[ 101] = "SiGe";
  SilMat[ 102] = "4H-SiC";
}


//Silvaco storage format
int SilvacoTIF::read()
{
  std::ifstream ctmp(_file.c_str(), std::ios::in);

  if (!ctmp.good())
  {
    std::cerr<<"Open silvaco file error."<<std::endl;
    return 1;
  }

  char flag;

  while(ctmp >> flag)
  {
    switch(flag)
    {
      // node coordinate
    case 'c' :
      {
        Node_t node;
        ctmp >> node.index >> node.x >> node.y >> node.z;
        // correct subscript
        node.index = node.index - 1;
        // save it
        _nodes.push_back(node);
      }
      break;

      // edge
    case 'e' :
      {
        Edge_t edge;
        ctmp >> edge.index >> edge.point1 >> edge.point2 >> edge.bcode;
        // correct subscript
        edge.index  = edge.index  -1;
        edge.point1 = edge.point1 -1;
        edge.point2 = edge.point2 -1;
        if(edge.point1 > edge.point2) std::swap(edge.point1, edge.point2);
        // save it
        _edges.push_back(edge);
      }
      break;

      // triangle
    case 't' :
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
      }
      break;

      //region
    case 'r' :
      {
        Region_t region;
        int material_index;
        ctmp >> region.index >> material_index;
        region.index = region.index - 1;
        region.node_num  = 0;
        region.tri_num   = 0;
        region.electrode = false;
        region.segment   = false;
        if(SilMat.find(material_index)==SilMat.end())
        {
          std::cout<<"Unrecognized material code " << material_index << "in silvaco file." << std::endl;
          return 1;
        }
        region.material  = SilMat[material_index];
        std::stringstream ss;
        ss << region.index;
        ss >> region.name;
        region.name      = "region" + region.name;
        // save it
        _regions.push_back(region);
      }
      break;

      // electrode info?
    case 'x' :
      {
        _electrode_info.clear();
        int n_electrode_infos;
        ctmp >> n_electrode_infos;
        for(int i=0; i<n_electrode_infos; ++i)
        {
          int code;
          ctmp >> code;
          _electrode_info.push_back(code);
        }
      }
      break;

      // electrode
    case 'w' :
      {
        // last region
        Region_t & region = _regions[_regions.size()-1];
        region.electrode = true;

        int region_flag;
        int electrode_index;
        int unknow;
        for(unsigned int i=0; i<_electrode_info.size(); ++i)
        {
          switch(_electrode_info[i])
          {
          case 43 :
            ctmp >> region.name;
	    // remove "" if exist
	    if(region.name.at(0)=='"')
	      region.name = region.name.substr(1, region.name.size()-2);
	    break;
          case 42 :
            ctmp >> region_flag;
            region.segment = (region_flag==2);
            break;
          case 46 :
            ctmp >> electrode_index;
            break;
          case 62 :
          default :
            ctmp >> unknow;
          }
        }
      }
      break;

      // region boundary
    case 'b' :
      {
        // last region
        Region_t & region = _regions[_regions.size()-1];
        int boundary_index;
        ctmp >> boundary_index;
        region.boundary.push_back(boundary_index-1);
      }
      break;

      // solution
    case 's' :
      {
        ctmp >> _sol_head.sol_num;
        for(int i = 0; i < _sol_head.sol_num; i++)
        {
          int loc;
          ctmp >> loc;
          //find the storage index....
          std::string sol_name;
          if( SilImp.find(loc) != SilImp.end() )
            sol_name = SilImp[loc];
          else
            sol_name = "unknow";
          _sol_head.sol_name_array.push_back(sol_name);
        }
      }
      break;


      // solution data
    case 'n' :
      {
        int loc;
        SolData_t solution;
        //find the storage index...
        ctmp >> solution.index >> loc;
        solution.material = SilMat[loc];

        //For all data values...
        for(int i = 0; i < _sol_head.sol_num; i++)
        {
          double dval;
          ctmp >> dval;
          solution.data_array.push_back(dval);
        }
        _sol_data.push_back(solution);
      }
      break;

    default :
      {
        std::string rubbish;
        std::getline(ctmp, rubbish);
      }
    }
  }
  ctmp.close();

  //statistic how many triangles in each region
  for(unsigned int n=0; n<_tris.size(); ++n)
    _regions[_tris[n].region].tri_num++;

  return 0;
}

