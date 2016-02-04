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

#include "medici.h"
#include "material_define.h"


MediciTIF::MediciTIF(const std::string & file)
  :StanfordTIF(file)
{
  _version = "MEDICI";

  // solution name map
  _solution_name_map["Acceptor"]                          = "Accept";

  // map sdevice impurity name to genius impurity name
  _solution_name_map["boron"]                             = "Boron";
  _solution_name_map["BoronConcentration"]                = "Boron";
  _solution_name_map["BoronChemicalConcentration"]        = "Boron";
  _solution_name_map["a_boron"]                           = "BoronActive";
  _solution_name_map["Ba"]                                = "BoronActive";


  _solution_name_map["phosphorus"]                        = "Phosphorus";
  _solution_name_map["PhosphorusConcentration"]           = "Phosphorus";
  _solution_name_map["PhosphorusChemicalConcentration"]   = "Phosphorus";
  _solution_name_map["a_phosphorus"]                      = "PhosphorusActive";
  _solution_name_map["Pa"]                                = "PhosphorusActive";

  
  _solution_name_map["arsenic"]                           = "Arsenic";
  _solution_name_map["ArsenicConcentration"]              = "Arsenic";
  _solution_name_map["ArsenicChemicalConcentration"]      = "Arsenic";
  _solution_name_map["a_arsenic"]                         = "ArsenicActive";
  _solution_name_map["Asa"]                               = "ArsenicActive";

  _solution_name_map["antimony"]                          = "Antimony";
  _solution_name_map["AntimonyConcentration"]             = "Antimony";
  _solution_name_map["AntimonyChemicalConcentration"]     = "Antimony";
  _solution_name_map["a_antimony"]                        = "AntimonyActive";
  _solution_name_map["Sba"]                               = "AntimonyActive";

}



bool MediciTIF::read(std::string &err)
{
  std::ifstream ctmp(_file.c_str(), std::ios::in);

  if (!ctmp.good())
  {
    err = "Open Medici TIF file error.";
    return false;
  }

  std::string buffer;
  std::string flag;

  while(ctmp >> flag)
  {
    // the version of tif file
    // the version of suprem file
    if (flag == "h")
    {
      // skip tif file header
      std::string buf;
      std::getline(ctmp, buf);
      if(buf.find("MEDICI") != std::string::npos) _version = "MEDICI";
      else _version = "TIF";
    }

    else if (flag == "v")
    {
      // skip tif file header
      std::string buf;
      std::getline(ctmp, buf);
      if(buf.find("TMA") != std::string::npos) _version = "TMA";
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
      if(t.t1>0) t.t1 = t.t1 -1;
      if(t.t2>0) t.t2 = t.t2 -1;
      if(t.t3>0) t.t3 = t.t3 -1;
      t.c4 = -1;
      t.t4 = -1;
      // save it
      _tris.push_back(t);
      continue;
    }

    // quadangle
    else if (flag == "q")
    {
      Tri_t t;
      ctmp >> t.index >> t.region >> t.c1 >> t.c2 >> t.c3 >> t.c4 >> t.t1 >> t.t2 >> t.t3 >> t.t4;
      // correct subscript
      t.index = t.index - 1;
      t.region = t.region - 1;
      t.c1 = t.c1 -1;
      t.c2 = t.c2 -1;
      t.c3 = t.c3 -1;
      t.c4 = t.c4 -1;
      if(t.t1>0) t.t1 = t.t1 -1;
      if(t.t2>0) t.t2 = t.t2 -1;
      if(t.t3>0) t.t3 = t.t3 -1;
      if(t.t4>0) t.t4 = t.t4 -1;
      // save it
      _tris.push_back(t);
      continue;
    }


    //region
    else if (flag == "r")
    {
      Region_t region;
      if(_version != "TMA")
      {
        ctmp >> region.index >> region.material >> region.name;
        //add "r" before region name if it begin with number
        if(isdigit(region.name[0]))
          region.name = std::string("r") + region.name;
      }
      else
      {
        ctmp >> region.index >> region.material;
        std::stringstream ss;
        ss << 'r' <<  region.index;
        region.name = ss.str();
      }
      region.index = region.index - 1;
      region.node_num  = 0;
      region.tri_num   = 0;
      region.region = region.index;
      region.segment = false;

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
      if(region.region == 0)
      {
        region.index = region.index - 1;
        region.node_num  = 0;
        region.tri_num   = 0;
        region.region = -1;
        region.segment = true;

        //add "i" before region name if it begin with number
        if(isdigit(region.name[0]))
          region.name = std::string("i") + region.name;
        // save it
        _regions.push_back(region);
      }
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
        _sol_head.sol_name_array.push_back(_solution_name_format(sol_name));
      }
      continue;
    }

    // solution units
    else if (flag == "u")
    {
      int unit_num;
      ctmp >> unit_num;
      for(int i = 0; i < unit_num; i++)
      {
        std::string sol_unit;
        ctmp >> sol_unit;
        _sol_head.sol_unit_array.push_back(sol_unit);
      }
      continue;
    }

    // solution data
    else if (flag == "n")
    {
      SolData_t solution;

      if(_version != "TMA")
      {
        ctmp >> solution.index >> solution.material;
        solution.index = solution.index - 1;
        solution.region_index=-1;
      }
      else
      {
        ctmp >> solution.index >> solution.region_index;
        solution.index = solution.index - 1;
        solution.region_index -= 1;
        if(solution.region_index >= 0)
          solution.material = _regions[solution.region_index].material;
        else
          solution.material = "Vacuum";
      }

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

  if(_nodes.size()==0 || _tris.size()==0 )
  {
    err = "Empty Medici TIF file.";
    return false;
  }

  //statistic how many triangles in each region
  for(unsigned int n=0; n<_tris.size(); ++n)
    _regions[_tris[n].region].tri_num++;

  //
  if(!_check_duplicate_name())
  {
    err = "TIF file with duplicate region/boundary names.";
    return false;
  }

  _find_solution_region_by_material();

  // build Accept and Donor when not exist
  if( _version != "TMA" && !_sol_head.has_solution("Accept") && !_sol_head.has_solution("Donor") )
  {
    if(_sol_head.has_solution("Net")  &&  _sol_head.has_solution("Total") )
    {
      _sol_head.sol_num += 2;
      _sol_head.sol_name_array.push_back("Accept");
      _sol_head.sol_name_array.push_back("Donor");

      unsigned int net_index = _sol_head.solution_index("Net");
      unsigned int tot_index = _sol_head.solution_index("Total");

      for(unsigned int n=0; n<_sol_data.size(); ++n)
      {
        SolData_t & sol_date = _sol_data[n];
        double net = sol_date.data_array[net_index];
        double tot = sol_date.data_array[tot_index];
        sol_date.data_array.push_back(0.5*(tot-net));
        sol_date.data_array.push_back(0.5*(tot+net));
      }
    }
  }

  if(!_sol_head.has_solution("Accept") && !_sol_head.has_solution("Donor"))
  {
    unsigned int boron_index      = _sol_head.solution_index("BoronActive");
    unsigned int aluminum_index   = _sol_head.solution_index("AluminumActive");
    
    unsigned int nitrogen_index   = _sol_head.solution_index("NitrogenActive");
    unsigned int phosphorus_index = _sol_head.solution_index("PhosphorusActive");
    unsigned int arsenic_index    = _sol_head.solution_index("ArsenicActive");
    unsigned int antimony_index   = _sol_head.solution_index("AntimonyActive");

    _sol_head.sol_num += 2;
    _sol_head.sol_name_array.push_back("Accept");
    _sol_head.sol_name_array.push_back("Donor");

    for(unsigned int n=0; n<_sol_data.size(); ++n)
    {
      SolData_t & sol_date = _sol_data[n];
      double nd = 0.0;
      double na = 0.0;
      if(boron_index != static_cast<unsigned int>(-1))
        na += sol_date.data_array[boron_index];
      if(aluminum_index != static_cast<unsigned int>(-1))
        na += sol_date.data_array[aluminum_index];

      if(nitrogen_index != static_cast<unsigned int>(-1))
        nd += sol_date.data_array[nitrogen_index];
      if(phosphorus_index != static_cast<unsigned int>(-1))
        nd += sol_date.data_array[phosphorus_index];
      if(arsenic_index != static_cast<unsigned int>(-1))
        nd += sol_date.data_array[arsenic_index];
      if(antimony_index != static_cast<unsigned int>(-1))
        nd += sol_date.data_array[antimony_index];

      sol_date.data_array.push_back(na);
      sol_date.data_array.push_back(nd);
    }
  }


  if(_sol_head.has_solution("Accept"))
    _solution_index_map.insert(std::make_pair("acceptor_index", _sol_head.solution_index("Accept")));

  if(_sol_head.has_solution("Donor"))
    _solution_index_map.insert(std::make_pair("donor_index", _sol_head.solution_index("Donor")));

  if(_sol_head.has_solution("mole_x"))
    _solution_index_map.insert(std::make_pair("mole_x_index", _sol_head.solution_index("mole_x")));

  if(_sol_head.has_solution("mole_y"))
    _solution_index_map.insert(std::make_pair("mole_y_index", _sol_head.solution_index("mole_y")));

  return true;

}



void MediciTIF::export_tif(const std::string & file) const
{
    std::ofstream fout;
    fout.open ( file.c_str(), std::ofstream::trunc );
    // set the float number precision
    fout.precision ( 8 );

    // set output width and format
    fout<< std::scientific << std::right;

    time_t          _time;
    time ( &_time );

    fout << "h TIF V1.2.1 created by Genius, Copyright (C) by Cogenda Pte. Ltd. Date: " << ctime ( &_time ) ; // ctime with '\n'
    fout << "cd GEN          blnk               blnk          blnk        cart2D    1.00000E+00  0.00000E+00"  << '\n';
    fout << "cg   3.00000E+02" << '\n';

    //write point
    for ( unsigned int i=0; i<_nodes.size(); ++i )
    {
      fout<< 'c'
      << std::setw ( 8 ) << _nodes[i].index+1 << "   "
      << std::setw ( 15 ) << _nodes[i].x      << "   "
      << std::setw ( 15 ) << _nodes[i].y      << "   "
      << _nodes[i].h << '\n';
    }
    //fout<<'\n';

    //write edge
    for ( unsigned int i=0; i<_edges.size(); ++i )
    {
      fout<< 'e'
      << std::setw ( 8 ) << _edges[i].index+1    << "   "
      << std::setw ( 8 ) << _edges[i].point1 + 1 << "   "
      << std::setw ( 8 ) << _edges[i].point2 + 1 << "   "
      << 0 << '\n';
    }
    //fout<<'\n';

    // write region information
    int region_count = 0;
    int segment_count = 0;
    std::map<int, int> region_map;
    for ( unsigned int r=0; r<_regions.size(); ++r )
    {
      if(_regions[r].segment == false)
      {
        fout<< 'r'                  << "   "
            << 1 + region_count  << "   "
            << _regions[r].material << "   "
            << _regions[r].name << '\n';
        region_map[r] = region_count;
        const std::vector<int> & boundary = _regions[r].boundary;
        for(unsigned int b=0; b<boundary.size(); ++b)
          fout<< " b" << "   " << boundary[b]+1 << '\n';

        region_count++;
      }
      else
      {
        fout<< 'i'                  << "   "
            << 1 + segment_count    << "   "
            << _regions[r].material << "   "
            << _regions[r].name     << "   "
            << 0 << '\n';

        const std::vector<int> & boundary = _regions[r].boundary;
        for(unsigned int b=0; b<boundary.size(); ++b)
          fout<< " j" << "   " << boundary[b]+1 << '\n';

        segment_count++;
      }
    }
    //fout<<'\n';


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
            << std::setw ( 8 )<< (_tris[i].t3>=0 ? _tris[i].t3+1 : _tris[i].t3) << '\n';
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
            << std::setw ( 8 )<< (_tris[i].t4>=0 ? _tris[i].t4+1 : _tris[i].t4) << '\n';
      }
    }
    //fout<<'\n';


    // write profile information
    if(_sol_head.sol_num)
    {
      fout<< "s    " << _sol_head.sol_name_array.size() << "  ";
      for(unsigned int s=0; s<_sol_head.sol_name_array.size(); ++s)
        fout<< _sol_head.sol_name_array[s] << "  ";
      fout<<std::endl;
      
      if(!_sol_head.sol_unit_array.empty())
      {
        fout<< "u    " << _sol_head.sol_unit_array.size() << "  ";
        for(unsigned int u=0; u<_sol_head.sol_unit_array.size(); ++u)
          fout<< _sol_head.sol_unit_array[u] << "  ";
        fout<<std::endl;
      }
      
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
        << std::setw(8) << index+1 << " "
        << std::setw(15)<< material << "   ";

        unsigned int d = it->second;
        for(unsigned int s=0; s<_sol_data[d].data_array.size(); ++s)
          fout<< std::setw(15) << _sol_data[d].data_array[s] << "  ";

        fout<<std::endl;
      }
    }

    fout.close();
}


std::string MediciTIF::_solution_name_format( const std::string & name) const
{
  if( _solution_name_map.find(name) != _solution_name_map.end() )
    return _solution_name_map.find(name)->second;

  return name;
}



