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

#include "silvaco.h"

SilvacoTIF::SilvacoTIF()
: _dim(2)
{
  _init_index_string_map();
}


SilvacoTIF::SilvacoTIF(const std::string & file) : StanfordTIF(file), _dim(2)
{
  _init_index_string_map();
}


//Silvaco storage format
bool SilvacoTIF::read(std::string &err)
{
  std::ifstream ctmp(_file.c_str(), std::ios::in);

  if (!ctmp.good())
  {
    err = "Open silvaco file error.";
    return false;
  }

  char flag;

  while(ctmp >> flag)
  {
    switch(flag)
    {
    case 'v' :
      {
        ctmp >> _version;
        if( _version != "ATLAS" && _version != "ATHENA" && _version != "DEVEDIT")
        {
          err = "Only support silvaco ATLAS/ATHENA/DEVEDIT file, but get " + _version + " file.";
          return false;
        }
        break;
      }
    case 'k' :
      {
        if( _version == "ATLAS" )
        {
          unsigned int dmesh;
          double dm1, dm2;
          ctmp >> _dim >> dmesh >> dm1 >> dm2;
        }

        if( _version == "ATHENA" )
        {
          _dim = 2;
          std::string rubbish;
          std::getline(ctmp, rubbish);
        }

        if( _version == "DEVEDIT" )
        {
          _dim = 2;
          std::string rubbish;
          std::getline(ctmp, rubbish);
        }

        break;
      }
      // node coordinate
    case 'c' :
      {
        Node_t node;
        ctmp >> node.index >> node.x >> node.y >> node.z;
        node.h = node.z;
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
        t.c4 = -1;
        t.t4 = -1;
        // save it
        _tris.push_back(t);
      }
      break;

    case 'Z' :
      {
        int index;
        double z;
        ctmp >> index >> z;
        _z_slice.push_back(z);
        break;
      }

    case 'P' :
      {
        Prism_t Prism;
        int z1, z2;
        ctmp >> Prism.index >> Prism.tri >> Prism.region >> z1 >> z2;
        Prism.index = Prism.index-1;
        Prism.tri = Prism.tri-1;
        Prism.region = Prism.region-1;
        Prism.z1 = z1-1;
        Prism.z2 = z2-1;

        _prisms.push_back(Prism);
        break;
      }
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
          std::stringstream ss;
          ss << "Unrecognized material code " << material_index << " in silvaco file." << std::endl;
          err = ss.str();
          return false;
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
        int z;
        int unknow;
        for(unsigned int i=0; i<_electrode_info.size(); ++i)
        {
          switch(_electrode_info[i])
          {
          case 42 :
              ctmp >> region_flag;
              region.segment = (region_flag==2);
              break;
          case 43 :
            ctmp >> region.name;
	    // remove "" if exist
	    if(region.name.at(0)=='"')
	      region.name = region.name.substr(1, region.name.size()-2);
	    break;
          case 44 :
            ctmp >> region.z2;
            region.z2 = region.z2 -1;
            break;
          case 45:
            ctmp >> region.z1;
            region.z1 = region.z1 -1;
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
        solution.region_index = -1;
        //For all data values...
        for(int i = 0; i < _sol_head.sol_num; i++)
        {
          if( _sol_head.sol_name_array[i] == "Region" )
          {
            ctmp >> solution.region_index ;
            solution.region_index =  solution.region_index-1;
            solution.data_array.push_back(solution.region_index);
            continue;
          }

          if( _sol_head.sol_name_array[i] == "ZPlaneIndex" )
          {
            ctmp >> solution.zplane ;
            solution.zplane =  solution.zplane-1;
            solution.data_array.push_back(solution.zplane);
            continue;
          }

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

  if(_dim == 2)
  {
    //statistic how many triangles in each region
    for(unsigned int n=0; n<_tris.size(); ++n)
      _regions[_tris[n].region].tri_num++;

    _find_solution_region_by_material();
  }

  //
  if(!_check_duplicate_name())
  {
    err = "Silvaco file with duplicate region/boundary names.";
    return false;
  }

  if(_sol_head.has_solution("Acceptor"))
    _solution_index_map.insert(std::make_pair("acceptor_index", _sol_head.solution_index("Acceptor")));

  if(_sol_head.has_solution("Donor"))
    _solution_index_map.insert(std::make_pair("donor_index", _sol_head.solution_index("Donor")));

  if(_sol_head.has_solution("mole_x"))
    _solution_index_map.insert(std::make_pair("mole_x_index", _sol_head.solution_index("mole_x")));

  if(_sol_head.has_solution("mole_y"))
    _solution_index_map.insert(std::make_pair("mole_y_index", _sol_head.solution_index("mole_y")));

  if(_sol_head.has_solution("Arsenic"))
    _solution_index_map.insert(std::make_pair("As_index", _sol_head.solution_index("Arsenic")));

  if(_sol_head.has_solution("Phosphorus"))
    _solution_index_map.insert(std::make_pair("P_index", _sol_head.solution_index("Phosphorus")));

  if(_sol_head.has_solution("Antimony"))
    _solution_index_map.insert(std::make_pair("Sb_index", _sol_head.solution_index("Antimony")));

  if(_sol_head.has_solution("Boron"))
    _solution_index_map.insert(std::make_pair("B_index", _sol_head.solution_index("Boron")));


  //export_scatter_doping_na();
  //export_scatter_doping_nd();

  return true;
}





void SilvacoTIF::export_scatter_doping_na() const
{
  std::string filename = _file+".na";
  std::ofstream fout( filename.c_str(), std::ofstream::trunc);

  int donor = _sol_head.solution_index("Donor");
  int acceptor = _sol_head.solution_index("Acceptor");
  std::map<int, double> node_doping;

  std::vector<SolData_t>::const_iterator data_it = _sol_data.begin();
  for( ; data_it != _sol_data.end(); ++data_it)
  {
    int index = data_it->index;
    double doping = data_it->data_array[acceptor];
    if(node_doping.find(index)!=node_doping.end())
    {
      double origin_doping = node_doping.find(index)->second;
      if( fabs(doping) > fabs(origin_doping) )
        node_doping.insert(std::make_pair(index, doping));
    }
    else
      node_doping.insert(std::make_pair(index, doping));
  }

  std::vector<Node_t>::const_iterator node_it= _nodes.begin();
  for(; node_it!= _nodes.end(); ++node_it)
  {
    int index = node_it->index;
    assert(node_doping.find(index)!=node_doping.end());
    fout << node_it->x << " " << node_it->y << " " << node_doping.find(index)->second << std::endl;
  }

  fout.close();
}


void SilvacoTIF::export_scatter_doping_nd() const
{
  std::string filename = _file+".nd";
  std::ofstream fout( filename.c_str(), std::ofstream::trunc);

  int donor = _sol_head.solution_index("Donor");
  int acceptor = _sol_head.solution_index("Acceptor");
  std::map<int, double> node_doping;

  std::vector<SolData_t>::const_iterator data_it = _sol_data.begin();
  for( ; data_it != _sol_data.end(); ++data_it)
  {
    int index = data_it->index;
    double doping = data_it->data_array[donor];
    if(node_doping.find(index)!=node_doping.end())
    {
      double origin_doping = node_doping.find(index)->second;
      if( fabs(doping) > fabs(origin_doping) )
        node_doping.insert(std::make_pair(index, doping));
    }
    else
      node_doping.insert(std::make_pair(index, doping));
  }

  std::vector<Node_t>::const_iterator node_it= _nodes.begin();
  for(; node_it!= _nodes.end(); ++node_it)
  {
    int index = node_it->index;
    assert(node_doping.find(index)!=node_doping.end());
    fout << node_it->x << " " << node_it->y << " " << node_doping.find(index)->second << std::endl;
  }

  fout.close();
}


void SilvacoTIF::_init_index_string_map()
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

  SilImp[ 513] = "Region";
  SilImp[ 600] = "ZPlaneIndex";

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
  SilMat[   9] = "Sapphire";
  SilMat[  10] = "Au";
  SilMat[  11] = "Ag";
  SilMat[  12] = "AlSi";
  SilMat[  13] = "W";
  SilMat[  14] = "Ti";
  SilMat[  17] = "Co";
  SilMat[  21] = "Ta";
  SilMat[  24] = "AlGaAs";
  SilMat[  25] = "InGaAs";
  SilMat[  26] = "AlInAs";
  SilMat[  27] = "InP";
  SilMat[  28] = "Vacuum";
  SilMat[  32] = "Ge";
  SilMat[  81] = "TiSi2";
  SilMat[  91] = "Elec";
  SilMat[  99] = "3C-SiC";
  SilMat[ 100] = "Diamond";
  SilMat[ 101] = "SiGe";
  SilMat[ 102] = "6H-SiC";
  SilMat[ 103] = "4H-SiC";
  SilMat[ 104] = "AlP";
  SilMat[ 105] = "AlSb";
  SilMat[ 106] = "GaSb";
  SilMat[ 107] = "GaP";
  SilMat[ 108] = "InSb";
  SilMat[ 109] = "InAs";
  SilMat[ 110] = "ZnS";
  SilMat[ 111] = "ZnSe";
  SilMat[ 112] = "ZnTe";
  SilMat[ 113] = "CdS";
  SilMat[ 114] = "CdSe";
  SilMat[ 115] = "CdTe";
  SilMat[ 116] = "HgS";
  SilMat[ 117] = "HgSe";
  SilMat[ 118] = "HgTe";
  SilMat[ 119] = "PbS";
  SilMat[ 120] = "PbSe";
  SilMat[ 121] = "PbTe";
  SilMat[ 122] = "SnTe";
  SilMat[ 123] = "ScN";
  SilMat[ 124] = "GaN";
  SilMat[ 125] = "AlN";
  SilMat[ 126] = "InN";
  SilMat[ 127] = "BeTe";
  SilMat[ 128] = "InGaP";
  SilMat[ 129] = "GaSbP";
  SilMat[ 130] = "GaSbAs";
  SilMat[ 131] = "InAlAs";
  SilMat[ 132] = "InAsP";
  SilMat[ 133] = "GaAsP";
  SilMat[ 134] = "HgCdTe";
  SilMat[ 135] = "InGaAsP";
  SilMat[ 136] = "AlGaAsP";
  SilMat[ 137] = "AlGaAsSb";
  SilMat[ 140] = "SiN";
  SilMat[ 143] = "Si";
  SilMat[ 144] = "Polymer";
  SilMat[ 145] = "CuInGaSe";
  SilMat[ 146] = "InGaN";
  SilMat[ 147] = "AlGaN";
  SilMat[ 148] = "InAlGaN";
  SilMat[ 149] = "InGaNAs";
  SilMat[ 198] = "ITO";
   
}
 

