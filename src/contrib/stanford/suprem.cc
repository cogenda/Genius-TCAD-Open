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





#include "suprem.h"


SupremTIF::SupremTIF(const std::string & file)
  :StanfordTIF(file)
{
  _init_index_string_map();
}



bool SupremTIF::read()
{
  std::ifstream ctmp(_file.c_str(), std::ios::in);

  if (!ctmp.good())
  {
    std::cerr<<"Open Suprem4GS file error."<<std::endl;
    return false;
  }

  char buffer[1024];
  std::string flag;

  while(ctmp >> flag)
  {

    // the version of suprem file
    if (flag == "v")
    {
      // skip suprem file version info
      ctmp.getline(buffer, 1024);
    }

    // dimension info
    else if (flag == "D")
    {
      int mode, nvrt, nedg;
      ctmp >> mode >> nvrt >> nedg;
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
    }


    // triangle
    else if (flag == "t")
    {
      Tri_t t;
      ctmp >> t.index >> t.region >> t.c1 >> t.c2 >> t.c3 >> t.t1 >> t.t2 >> t.t3 >> t.father >> t.offspr;
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


    //region
    else if (flag == "r")
    {
      Region_t region;
      int material_id;
      ctmp >> region.index >> material_id;
      region.index = region.index - 1;
      region.node_num  = 0;
      region.tri_num   = 0;
      region.material = _material_index_to_string[material_id];
      region.segment = false;
      std::stringstream ss;
      ss << "region_" << region.index;
      region.name = ss.str();
      // save it
      _regions.push_back(region);
    }


    // solutions
    else if (flag == "s")
    {
      ctmp >> _sol_head.sol_num;
      for(int i = 0; i < _sol_head.sol_num; i++)
      {
        int sol_index;
        ctmp >> sol_index;
        std::string sol_name = _doping_index_to_string[sol_index];
        _sol_head.sol_name_array.push_back(sol_name);
      }
    }

    // solution data
    else if (flag == "n")
    {
      SolData_t solution;
      int material_id;
      ctmp >> solution.index >> material_id;
      solution.material = _material_index_to_string[material_id];
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
      ctmp.getline(buffer, BUFSIZ);
  }


  ctmp.close();

  if(_nodes.size()==0 || _tris.size()==0 ) return false;

  //statistic how many triangles in each region
  for(unsigned int n=0; n<_tris.size(); ++n)
    _regions[_tris[n].region].tri_num++;

  _find_solution_region_by_material();

  return true;

}



void SupremTIF::_init_index_string_map()
{
  {
    /*top surface material number - gas contact*/
    _material_string_to_index["Vacuum"  ] = 0;

    _material_string_to_index["SiO2" ] = 1;
    _material_string_to_index["Si3N4"] = 2;       /*nitride*/
    _material_string_to_index["OxNi" ] = 5 ;      /*oxynitride*/
    _material_string_to_index["PhRs" ] = 7;       /*photoresist*/


    /*define numbers for different types of semiconductors*/
    _material_string_to_index["Si"   ] = 3;
    _material_string_to_index["Poly" ] = 4;
    _material_string_to_index["GaAs" ] = 8;

    /*let's not leave metals out*/
    _material_string_to_index["Al"  ] = 6;

    std::map<std::string, int>::const_iterator it = _material_string_to_index.begin();
    for(; it != _material_string_to_index.end(); ++it)
    {
      _material_index_to_string[it->second] = it->first;
    }
  }

  {

    _doping_string_to_index["V"   ] = 0;          /*Vacancies*/
    _doping_string_to_index["I"   ] = 1;          /*Interstitials*/
    _doping_string_to_index["As"  ] = 2;          /*Arsenic*/
    _doping_string_to_index["P"   ] = 3;          /*Phosphorus*/
    _doping_string_to_index["Sb"  ] = 4;          /*antimony*/
    _doping_string_to_index["B"   ] = 5;          /*Boron*/
    _doping_string_to_index["N"   ] = 6;          /*electron concentration*/
    _doping_string_to_index["H"   ] = 7;          /*hole concentration*/
    _doping_string_to_index["XVEL"] = 8;          /*X velocity*/
    _doping_string_to_index["YVEL"] = 9;          /*Y velocity*/
    _doping_string_to_index["O2"  ] = 10;         /*Dry O2*/
    _doping_string_to_index["H2O" ] = 11;         /*Wet O2*/
    _doping_string_to_index["T"   ] = 12;           /*Interstitial Traps*/
    _doping_string_to_index["Au"  ] = 13;         /*Gold*/
    _doping_string_to_index["Psi" ] = 14;         /*Potential*/
    _doping_string_to_index["Sxx" ] = 15;         /* Components of stress - not solution variables*/
    _doping_string_to_index["Syy" ] = 16;         /* but too expensive to recompute for plotting */
    _doping_string_to_index["Sxy" ] = 17;
    _doping_string_to_index["Cs"  ] = 18;         /* Cesium for oxide charges */
    _doping_string_to_index["DELA"] = 19;         /*change in interface area - not solution variable*/
    _doping_string_to_index["Asa" ] = 20;         /*Arsenic active concentration*/
    _doping_string_to_index["Pa"  ] = 21;         /*Phosphorus active concentration*/
    _doping_string_to_index["Sba" ] = 22;         /*antimony active concentration*/
    _doping_string_to_index["Ba"  ] = 23;         /*Boron active concentration*/
    _doping_string_to_index["GRN" ] = 24;         /*Polysilicon grain size*/
    _doping_string_to_index["Ga"  ] = 25;         /*p-type tracer to help model miyake*/
    _doping_string_to_index["iBe" ] = 31;         /*Beryllium impurity*/
    _doping_string_to_index["iBea"] = 32;         /*Beryllium active concentration*/
    _doping_string_to_index["iMg" ] = 33;         /*Magnesium impurity*/
    _doping_string_to_index["iMga"] = 34;         /*Magnesium active concentration*/
    _doping_string_to_index["iSe" ] = 35;         /*Selenium impurity*/
    _doping_string_to_index["iSea"] = 36;         /*Selenium active concentration*/
    _doping_string_to_index["iSi" ] = 37;         /*Silicon impurity*/
    _doping_string_to_index["iSia"] = 38;         /*Silcon active concentration*/
    _doping_string_to_index["iSn" ] = 39;         /*Tin impurity*/
    _doping_string_to_index["iSna"] = 40;         /*Tin active concentration*/
    _doping_string_to_index["iGe" ] = 41;         /*Germanium impurity*/
    _doping_string_to_index["iGea"] = 42;         /*Germanium active concentration*/
    _doping_string_to_index["iZn" ] = 43;         /*Zinc impurity*/
    _doping_string_to_index["iZna"] = 44;         /*Zinc active concentration*/
    _doping_string_to_index["iC"  ] = 45;         /*Carbon impurity*/
    _doping_string_to_index["iCa" ] = 46;         /*Carbon active concentration*/
    _doping_string_to_index["iG"  ] = 47;         /*Generic impurity*/
    _doping_string_to_index["iGa" ] = 48;         /*Generic active concentration*/


    std::map<std::string, int>::const_iterator it = _doping_string_to_index.begin();
    for(; it != _doping_string_to_index.end(); ++it)
    {
      _doping_index_to_string[it->second] = it->first;
    }
  }
}



