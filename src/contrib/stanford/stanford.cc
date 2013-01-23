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
    std::vector<double> solution_array;
    if(Genius::processor_id() == 0)
    {
      for(unsigned int n=0; n < n_solution; ++n)
      {
        const SolData_t & data = _sol_data[n];
        solution_index.push_back(data.index);
        solution_region.push_back(data.region_index);
        solution_array.insert(solution_array.end(), data.data_array.begin(), data.data_array.end());
      }
    }
    Parallel::broadcast(solution_index );
    Parallel::broadcast(solution_region );
    Parallel::broadcast(solution_array );

    if(Genius::processor_id() != 0)
    {
      int sol = _sol_head.sol_num;
      for(unsigned int n=0; n < n_solution; ++n)
      {
        SolData_t data;
        data.index = solution_index[n];
        data.region_index = solution_region[n];
        for(int s=0; s<sol; s++)
          data.data_array.push_back(solution_array[n*sol+s]);
        _sol_data.push_back(data);
      }
    }
  }

  Parallel::broadcast(_acceptor_index, root);
  Parallel::broadcast(_donor_index, root);
  Parallel::broadcast(_mole_x_index, root);
  Parallel::broadcast(_mole_y_index, root);
}


void StanfordTIF::_find_solution_region_by_material()
{
  std::map<int, std::set<int> > node_region_map; // <node, <region> >
  for(unsigned int n=0; n<_tris.size(); ++n)
  {
    int region = _tris[n].region;
    node_region_map[_tris[n].c1].insert(region);
    node_region_map[_tris[n].c2].insert(region);
    node_region_map[_tris[n].c3].insert(region);
  }

  std::map<int, std::string> region_material;
  for(unsigned int n=0; n<_regions.size(); ++n)
  {
    if(!_regions[n].segment)
      region_material[_regions[n].index] = _regions[n].material;
  }

  for(unsigned int n=0; n<_sol_data.size(); ++n)
  {
    int node = _sol_data[n].index;
    std::set<int> & regions =  node_region_map.find(node)->second;

    for(std::set<int>::const_iterator it = regions.begin(); it != regions.end(); ++it)
    {
      if( region_material.find(*it)->second == _sol_data[n].material )
      {
        _sol_data[n].region_index = *it;
        regions.erase(*it); // remove candidate region
        break;
      }
    }
  }
}



double StanfordTIF::acceptor(unsigned int data_index) const
{
  assert(data_index < _sol_data.size());
  if( _acceptor_index != static_cast<unsigned int>(-1) )
    return _sol_data[data_index].data_array[_acceptor_index];
  return 0.0;
}


double StanfordTIF::donor(unsigned int data_index) const
{
  assert(data_index < _sol_data.size());
  if( _donor_index != static_cast<unsigned int>(-1) )
    return _sol_data[data_index].data_array[_donor_index];
  return 0.0;
}

double StanfordTIF::mole_x(unsigned int data_index) const
{
  assert(data_index < _sol_data.size());
  if( _mole_x_index != static_cast<unsigned int>(-1) )
    return _sol_data[data_index].data_array[_mole_x_index];
  return 0.0;
}


double StanfordTIF::mole_y(unsigned int data_index) const
{
  assert(data_index < _sol_data.size());
  if( _mole_y_index != static_cast<unsigned int>(-1) )
    return _sol_data[data_index].data_array[_mole_y_index];
  return 0.0;
}




