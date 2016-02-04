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
#include <iomanip>




#include "suprem.h"

SupremTIF::SupremTIF()
{
  _init_index_string_map();
}


SupremTIF::SupremTIF(const std::string & file)
  :StanfordTIF(file)
{
  _init_index_string_map();
}



bool SupremTIF::read(std::string &err)
{
  std::ifstream ctmp(_file.c_str(), std::ios::in);

  if (!ctmp.good())
  {
    err = "Open Suprem4GS file error.";
    return false;
  }

  std::string buffer;
  std::string flag;

  while(ctmp >> flag)
  {

    // the version of suprem file
    if (flag == "v")
    {
      // skip suprem file version info
      std::getline(ctmp, buffer);
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
      if(t.t1>0) t.t1 = t.t1 -1;
      if(t.t2>0) t.t2 = t.t2 -1;
      if(t.t3>0) t.t3 = t.t3 -1;
      t.c4 = -1;
      t.t4 = -1;
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

    else if (flag == "M" || flag == "I")
    {
      std::string info;
      std::getline (ctmp,info);
      _sol_head.extra_infos.push_back(flag + " " + info);
    }

    else
      std::getline(ctmp, buffer);
  }


  ctmp.close();

  if(_nodes.size()==0 || _tris.size()==0 )
  {
    err = "Empty Suprem file.";
    return false;
  }

  //statistic how many triangles in each region
  for(unsigned int n=0; n<_tris.size(); ++n)
    _regions[_tris[n].region].tri_num++;

  //
  if(!_check_duplicate_name())
  {
    err = "Suprem file with duplicate region/boundary names.";
    return false;
  }

  _find_solution_region_by_material();

  return true;

}


void SupremTIF::export_sup(const std::string & file) const
{
    std::ofstream fout;
    fout.open ( file.c_str(), std::ofstream::trunc );
    // set the float number precision
    fout.precision ( 8 );

    // set output width and format
    fout<< std::scientific << std::right;

    fout << "v SUPREM-IV.GS B.9305" << '\n';
    fout << "D 2 3 3"  << '\n';


    //write point
    for ( unsigned int i=0; i<_nodes.size(); ++i )
    {
      fout<< 'c'
      << std::setw ( 8 ) << _nodes[i].index+1 << "   "
      << std::setw ( 15 ) << _nodes[i].x      << "   "
      << std::setw ( 15 ) << _nodes[i].y      << "   "
      << _nodes[i].h << '\n';
    }

    // write region information
    int region_count = 0;
    std::map<int, int> region_map;
    for ( unsigned int r=0; r<_regions.size(); ++r )
    {
      if(_regions[r].segment == false)
      {
        fout<< 'r'               << "   "
            << 1 + region_count  << "   "
            << _material_string_to_index.find(_regions[r].material)->second << '\n';
        region_map[r] = region_count;    
        region_count++;
      }
    }
 

    //write triangles
    for ( unsigned int i=0; i<_tris.size(); ++i )
    {
      if(_tris[i].c4 < 0)
      {
        fout<< 't'
            << std::setw ( 8 )<< _tris[i].index+1 << "  "
            << std::setw ( 8 )<< region_map[_tris[i].region] +1 << "  "
            << std::setw ( 8 )<< _tris[i].c1+1 << "  "
            << std::setw ( 8 )<< _tris[i].c2+1 << "  "
            << std::setw ( 8 )<< _tris[i].c3+1 << "  "
            << std::setw ( 8 )<< (_tris[i].t1>=0 ? _tris[i].t1+1 : _tris[i].t1) << "  "
            << std::setw ( 8 )<< (_tris[i].t2>=0 ? _tris[i].t2+1 : _tris[i].t2) << "  "
            << std::setw ( 8 )<< (_tris[i].t3>=0 ? _tris[i].t3+1 : _tris[i].t3) << "  "
            << std::setw ( 8 )<< _tris[i].father << "  "
            << std::setw ( 8 )<< _tris[i].offspr << "  "
            << '\n';
      }
      else
      {
        fout<< 'q'
            << std::setw ( 8 )<< _tris[i].index+1 << "  "
            << std::setw ( 8 )<< region_map[_tris[i].region] +1 << "  "
            << std::setw ( 8 )<< _tris[i].c1+1 << "  "
            << std::setw ( 8 )<< _tris[i].c2+1 << "  "
            << std::setw ( 8 )<< _tris[i].c3+1 << "  "
            << std::setw ( 8 )<< _tris[i].c4+1 << "  "
            << std::setw ( 8 )<< (_tris[i].t1>=0 ? _tris[i].t1+1 : _tris[i].t1) << "  "
            << std::setw ( 8 )<< (_tris[i].t2>=0 ? _tris[i].t2+1 : _tris[i].t2) << "  "
            << std::setw ( 8 )<< (_tris[i].t3>=0 ? _tris[i].t3+1 : _tris[i].t3) << "  "
            << std::setw ( 8 )<< (_tris[i].t4>=0 ? _tris[i].t4+1 : _tris[i].t4) << "  "
            << std::setw ( 8 )<< _tris[i].father << "  "
            << std::setw ( 8 )<< _tris[i].offspr << "  "
            << '\n';
      }
    }



    // write profile information
    if(_sol_head.sol_num)
    {
      fout<< "s    " << _sol_head.sol_name_array.size() << "  ";
      for(unsigned int s=0; s<_sol_head.sol_name_array.size(); ++s)
        fout<< _doping_string_to_index.find(_sol_head.sol_name_array[s])->second << "  ";
      fout<<std::endl;
      
      std::map<std::pair<int, std::string>, unsigned int>  region_solution_map; 
      for(unsigned int i=0; i<_sol_data.size(); ++i)
      {
        int index = _sol_data[i].index;
        std::string material = _sol_data[i].material;
        region_solution_map[ std::make_pair(index, material) ] = i;
      }
      
      
      std::map<std::pair<int, std::string>, unsigned int>::const_iterator it=region_solution_map.begin();
      for(; it!=region_solution_map.end(); it++)
      {
        int index = it->first.first;
        std::string material = it->first.second;
        fout<< 'n' << "   "
        << std::setw(8) << index << " "
        << std::setw(15)<< _material_string_to_index.find(material)->second << "   ";

        unsigned int d = it->second;
        for(unsigned int s=0; s<_sol_data[d].data_array.size(); ++s)
          fout<< std::setw(15) << _sol_data[d].data_array[s] << "  ";

        fout<<std::endl;
      }
    }
    
    for(unsigned int n=0; n<_sol_head.extra_infos.size(); ++n)
      fout<< _sol_head.extra_infos[n] << '\n';

    fout.close();
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
    _material_string_to_index["SiC6H" ] = 9;
    _material_string_to_index["SiC4H" ] = 10;
    _material_string_to_index["SiC3C" ] = 11;

    /*let's not leave metals out*/
    _material_string_to_index["Al"  ] = 6;

    std::map<std::string, int>::const_iterator it = _material_string_to_index.begin();
    for(; it != _material_string_to_index.end(); ++it)
    {
      _material_index_to_string[it->second] = it->first;
    }
  }

  {
     _doping_string_to_index["Vacancie"    ] = 0;          /*Vacancies*/
     _doping_string_to_index["Interstitial"] = 1;          /*Interstitials*/
     _doping_string_to_index["Arsenic"     ] = 2;          /*Arsenic*/
     _doping_string_to_index["Phosphorus"  ] = 3;          /*Phosphorus*/
     _doping_string_to_index["Antimony"    ] = 4;          /*antimony*/
     _doping_string_to_index["Boron"       ] = 5;          /*Boron*/
     _doping_string_to_index["BF2"         ] = 6;          /*BF2*/
     _doping_string_to_index["Aluminum"    ] = 7;          /*Aluminum*/
     _doping_string_to_index["Nitrogen"    ] = 8;          /*Nitrogen*/


     _doping_string_to_index["O2"          ] = 10;         /*Dry O2*/
     _doping_string_to_index["H2O"         ] = 11;         /*Wet O2*/
     _doping_string_to_index["Trap"        ] = 12;         /*Interstitial Traps*/
     _doping_string_to_index["Au"          ] = 13;         /*Gold*/

     _doping_string_to_index["Sxx" ] = 15;                 /* Components of stress - not solution variables*/
     _doping_string_to_index["Syy" ] = 16;                 /* but too expensive to recompute for plotting */
     _doping_string_to_index["Sxy" ] = 17;
     _doping_string_to_index["Cs"  ] = 18;                 /* Cesium for oxide charges */
     _doping_string_to_index["DELA"] = 19;                 /* change in interface area - not solution variable*/
     _doping_string_to_index["ArsenicActive"    ] = 20;       /*Arsenic active concentration*/
     _doping_string_to_index["PhosphorusActive" ] = 21;       /*Phosphorus active concentration*/
     _doping_string_to_index["AntimonyActive"   ] = 22;       /*antimony active concentration*/
     _doping_string_to_index["BoronActive"      ] = 23;       /*Boron active concentration*/

     _doping_string_to_index["AluminumActive" ] = 24;         /*Aluminum active concentration*/
     _doping_string_to_index["NitrogenActive" ] = 25;         /*Nitrogen active concentration*/

     _doping_string_to_index["Beryllium"      ] = 31;         /*Beryllium impurity*/
     _doping_string_to_index["BerylliumActive"] = 32;         /*Beryllium active concentration*/
     _doping_string_to_index["Magnesium"      ] = 33;         /*Magnesium impurity*/
     _doping_string_to_index["MagnesiumActive"] = 34;         /*Magnesium active concentration*/
     _doping_string_to_index["Selenium"       ] = 35;         /*Selenium impurity*/
     _doping_string_to_index["SeleniumActive" ] = 36;         /*Selenium active concentration*/
     _doping_string_to_index["Silicon"        ] = 37;         /*Silicon impurity*/
     _doping_string_to_index["SiliconActive"  ] = 38;         /*Silcon active concentration*/
     _doping_string_to_index["Tin"            ] = 39;         /*Tin impurity*/
     _doping_string_to_index["TinActive"      ] = 40;         /*Tin active concentration*/
     _doping_string_to_index["Germanium"      ] = 41;         /*Germanium impurity*/
     _doping_string_to_index["GermaniumActive"] = 42;         /*Germanium active concentration*/
     _doping_string_to_index["Zinc"           ] = 43;         /*Zinc impurity*/
     _doping_string_to_index["ZincActive"     ] = 44;         /*Zinc active concentration*/
     _doping_string_to_index["Carbon"         ] = 45;         /*Carbon impurity*/
     _doping_string_to_index["CarbonActive"   ] = 46;         /*Carbon active concentration*/
     _doping_string_to_index["Generic"        ] = 47;         /*Generic impurity*/
     _doping_string_to_index["GenericActive"  ] = 48;         /*Generic active concentration*/

     _doping_string_to_index["GRN" ] = 50;          /*Polysilicon grain size*/
     _doping_string_to_index["Ga"  ] = 51;          /*p-type tracer to help model miyake*/

     _doping_string_to_index["Psi" ] = 60;          /*Potential*/
     _doping_string_to_index["Elec"] = 61;          /*electron concentration*/
     _doping_string_to_index["Hole"] = 62;          /*hole concentration*/

     _doping_string_to_index["XVEL"] = 70;          /*X velocity*/
     _doping_string_to_index["YVEL"] = 71;          /*Y velocity*/



    std::map<std::string, int>::const_iterator it = _doping_string_to_index.begin();
    for(; it != _doping_string_to_index.end(); ++it)
    {
      _doping_index_to_string[it->second] = it->first;
    }
  }
}



