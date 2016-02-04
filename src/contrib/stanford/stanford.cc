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
#include <ctime>

#include <map>
#include <set>
#include <fstream>

#include "stanford.h"
#include "parallel.h"

void StanfordTIF::broadcast(unsigned int root)
{

 // broadcast SolHead_t to all processors
 Parallel::broadcast(_sol_head.sol_num, root);
 Parallel::broadcast(_sol_head.sol_name_array, root);

 //broadcast SolData to all processors
 {
   unsigned int n_solution = _sol_data.size();
   Parallel::broadcast(n_solution);

   std::vector<int> solution_index;
   std::vector<int> solution_region;
   std::vector<std::string> solution_material;
   std::vector<double> solution_array;
   if(Genius::processor_id() == 0)
   {
     for(unsigned int n=0; n < n_solution; ++n)
     {
       const SolData_t & data = _sol_data[n];
       solution_index.push_back(data.index);
       solution_region.push_back(data.region_index);
       solution_material.push_back(data.material);
       solution_array.insert(solution_array.end(), data.data_array.begin(), data.data_array.end());
     }
   }
   Parallel::broadcast(solution_index );
   Parallel::broadcast(solution_region );
   Parallel::broadcast(solution_material );
   Parallel::broadcast(solution_array );

   if(Genius::processor_id() != 0)
   {
     int sol = _sol_head.sol_num;
     for(unsigned int n=0; n < n_solution; ++n)
     {
       SolData_t data;
       data.index = solution_index[n];
       data.region_index = solution_region[n];
       data.material = solution_material[n];
       for(int s=0; s<sol; s++)
         data.data_array.push_back(solution_array[n*sol+s]);
       _sol_data.push_back(data);
       std::pair<unsigned int, std::string> key = std::make_pair(data.index, data.material);
       _solution_map.insert(std::make_pair(key, n));
     }
   }
 }

 Parallel::broadcast(_solution_index_map, root);
 
 Parallel::broadcast(_region_material, root);
}



bool StanfordTIF::_check_duplicate_name() const
{
  std::set<std::string> region_set;
  std::set<std::string> boundary_set;
  for(unsigned int n=0; n<_regions.size(); n++)
  {
    if( _regions[n].segment == false )
      region_set.insert(_regions[n].name);
    else
      boundary_set.insert(_regions[n].name);
  }

  return (region_set.size() + boundary_set.size()) == _regions.size();
}



void StanfordTIF::_find_solution_region_by_material()
{
  for(unsigned int n=0; n<_regions.size(); ++n)
  {
    if(!_regions[n].segment)
    {
      _region_material[_regions[n].index] = _regions[n].material;
    }
  }

  for(unsigned int n=0; n<this->sol_data_array().size(); ++n)
  {
    unsigned int node = _sol_data[n].index;
    std::string material = _sol_data[n].material;

    std::pair<unsigned int, std::string> key = std::make_pair(node, material);
    _solution_map.insert(std::make_pair(key, n));
  }

}



double StanfordTIF::solution(unsigned int sol_index, unsigned int data_index) const
{
  assert(data_index < _sol_data.size());
  if( sol_index != static_cast<unsigned int>(-1) )
    return _sol_data[data_index].data_array[sol_index];
  return 0.0;
}


double StanfordTIF::solution(unsigned int sol_index, unsigned int r, unsigned int n) const
{
  if(_region_material.find(r)==_region_material.end()) return 0.0;
  
  std::string material = _region_material.find(r)->second;
  std::pair<unsigned int, std::string> key = std::make_pair(n, material);
  if(_solution_map.find(key) == _solution_map.end())
    return 0.0;

  unsigned int data_index = _solution_map.find(key)->second;
  return _sol_data[data_index].data_array[sol_index];
}



double StanfordTIF::acceptor(unsigned int data_index) const
{
  assert(data_index < _sol_data.size());
  if( _solution_index_map.find("acceptor_index") != _solution_index_map.end() )
  {
    return _sol_data[data_index].data_array[_solution_index_map.find("acceptor_index")->second];
  }

  double a = 0.0;
  if( _solution_index_map.find("B_index") != _solution_index_map.end() )
    a += _sol_data[data_index].data_array[_solution_index_map.find("B_index")->second];

  return a;
}


double StanfordTIF::donor(unsigned int data_index) const
{
  assert(data_index < _sol_data.size());
  if( _solution_index_map.find("donor_index") != _solution_index_map.end() )
    return _sol_data[data_index].data_array[_solution_index_map.find("donor_index")->second];

  double d = 0.0;
  if( _solution_index_map.find("As_index") != _solution_index_map.end() )
    d += _sol_data[data_index].data_array[_solution_index_map.find("As_index")->second];
  if( _solution_index_map.find("P_index") != _solution_index_map.end() )
    d += _sol_data[data_index].data_array[_solution_index_map.find("P_index")->second];
  if( _solution_index_map.find("Sb_index") != _solution_index_map.end() )
    d += _sol_data[data_index].data_array[_solution_index_map.find("Sb_index")->second];

  return d;
}

double StanfordTIF::mole_x(unsigned int data_index) const
{
  assert(data_index < _sol_data.size());
  if( _solution_index_map.find("mole_x_index") != _solution_index_map.end() )
    return _sol_data[data_index].data_array[_solution_index_map.find("mole_x_index")->second];
  return 0.0;
}


double StanfordTIF::mole_y(unsigned int data_index) const
{
  assert(data_index < _sol_data.size());
  if( _solution_index_map.find("mole_y_index") != _solution_index_map.end() )
    return _sol_data[data_index].data_array[_solution_index_map.find("mole_y_index")->second];
  return 0.0;
}


double StanfordTIF::acceptor(unsigned int r, unsigned int n) const
{
  if(_region_material.find(r)==_region_material.end()) return 0.0;
  
  std::string material = _region_material.find(r)->second;
  std::pair<unsigned int, std::string> key = std::make_pair(n, material);
  if(_solution_map.find(key) == _solution_map.end())
    return 0.0;

  unsigned int data_index = _solution_map.find(key)->second;
  return acceptor(data_index);
}


double StanfordTIF::donor(unsigned int r, unsigned int n) const
{
  if(_region_material.find(r)==_region_material.end()) return 0.0;
  
  std::string material = _region_material.find(r)->second;
  std::pair<unsigned int, std::string> key = std::make_pair(n, material);
  if(_solution_map.find(key) == _solution_map.end())
    return 0.0;

  unsigned int data_index = _solution_map.find(key)->second;
  return donor(data_index);
}


double StanfordTIF::mole_x(unsigned int r, unsigned int n) const
{
  if(_region_material.find(r)==_region_material.end()) return 0.0;
  
  std::string material = _region_material.find(r)->second;
  std::pair<unsigned int, std::string> key = std::make_pair(n, material);
  if(_solution_map.find(key) == _solution_map.end())
    return 0.0;

  unsigned int data_index = _solution_map.find(key)->second;
  return mole_x(data_index);
}


double StanfordTIF::mole_y(unsigned int r, unsigned int n) const
{
  if(_region_material.find(r)==_region_material.end()) return 0.0;
  
  std::string material = _region_material.find(r)->second;
  std::pair<unsigned int, std::string> key = std::make_pair(n, material);
  if(_solution_map.find(key) == _solution_map.end())
    return 0.0;

  unsigned int data_index = _solution_map.find(key)->second;
  return mole_y(data_index);
}


StanfordTIF * StanfordTIF::merge( const std::vector<const StanfordTIF *> & meshes)
{
  StanfordTIF * tif = new StanfordTIF;

  std::vector<unsigned int> node_offset(1, 0);
  std::vector<unsigned int> edge_offset(1, 0);
  std::vector<unsigned int> elem_offset(1, 0);
  std::vector<unsigned int> region_offset(1, 0);
  std::vector<unsigned int> data_offset(1, 0);

  for(unsigned int m=0; m<meshes.size(); ++m)
  {
    const StanfordTIF * mesh = meshes[m];
    node_offset.push_back(node_offset[m] + mesh->_nodes.size());
    edge_offset.push_back(edge_offset[m] + mesh->_edges.size());
    elem_offset.push_back(elem_offset[m] + mesh->_tris.size());
    region_offset.push_back(region_offset[m] + mesh->_regions.size());
    data_offset.push_back(data_offset[m] + mesh->_sol_data.size());
  }

  for(unsigned int m=0; m<meshes.size(); ++m)
  {
    const StanfordTIF * mesh = meshes[m];
    // node
    for(unsigned int n=0; n<mesh->_nodes.size(); ++n)
    {
      Node_t node = mesh->_nodes[n];
      node.index += node_offset[m];
      tif->_nodes.push_back(node);
    }

    // edge
    for(unsigned int e=0; e<mesh->_edges.size(); ++e)
    {
      Edge_t edge = mesh->_edges[e];
      edge.index += edge_offset[m];
      edge.point1 += node_offset[m];
      edge.point2 += node_offset[m];
      tif->_edges.push_back(edge);
    }

    // region
    for(unsigned int r=0; r<mesh->_regions.size(); ++r)
    {
      Region_t region = mesh->_regions[r];
      region.index += region_offset[m];
      if(region.region != -1)
        region.region = region.index;
      for(unsigned int b=0; b<region.boundary.size(); b++)
      {
        region.boundary[b] +=  edge_offset[m];
      }
      tif->_regions.push_back(region);
    }

    // elem
    for(unsigned int t=0; t<mesh->_tris.size(); ++t)
    {
      Tri_t tri = mesh->_tris[t];
      tri.index += elem_offset[m];
      tri.region += region_offset[m];
      tri.c1 += node_offset[m];
      tri.c2 += node_offset[m];
      tri.c3 += node_offset[m];
      if(tri.c4 != -1)
        tri.c4 += node_offset[m];
      if(tri.t1 >= 0)
        tri.t1 += elem_offset[m];
      if(tri.t2 >= 0)
        tri.t2 += elem_offset[m];
      if(tri.t3 >= 0)
        tri.t3 += elem_offset[m];
      if(tri.t4 >= 0)
        tri.t4 += elem_offset[m];
      tif->_tris.push_back(tri);
    }
  }

  tif->_sol_head.sol_num = 4;
  tif->_sol_head.sol_name_array.push_back("Accept");
  tif->_sol_head.sol_name_array.push_back("Donor");
  tif->_sol_head.sol_name_array.push_back("mole_x");
  tif->_sol_head.sol_name_array.push_back("mole_y");


  for(unsigned int m=0; m<meshes.size(); ++m)
  {
    const StanfordTIF * mesh = meshes[m];
    for(unsigned int s=0; s<mesh->_sol_data.size(); ++s)
    {
       SolData_t sol = mesh->_sol_data[s];
       sol.index += node_offset[m];
       sol.region_index += region_offset[m];
       sol.data_array.clear();
       sol.data_array.push_back(mesh->acceptor(s));
       sol.data_array.push_back(mesh->donor(s));
       sol.data_array.push_back(mesh->mole_x(s));
       sol.data_array.push_back(mesh->mole_y(s));
       tif->_sol_data.push_back(sol);
    }
  }

  // post process
  tif->_find_solution_region_by_material();
  tif->_solution_index_map.insert(std::make_pair("acceptor_index", tif->_sol_head.solution_index("Accept")));
  tif->_solution_index_map.insert(std::make_pair("donor_index", tif->_sol_head.solution_index("Donor")));
  tif->_solution_index_map.insert(std::make_pair("mole_x_index", tif->_sol_head.solution_index("mole_x")));
  tif->_solution_index_map.insert(std::make_pair("mole_y_index", tif->_sol_head.solution_index("mole_y")));

  return tif;
}





