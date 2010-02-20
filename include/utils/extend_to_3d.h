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

#ifndef __extend_to_3d_h__
#define __extend_to_3d_h__

#include <vector>
#include <map>
#include <string>

#include "point.h"
#include "key.h"

class SimulationSystem;


/**
 * give a 2D mesh (tri/quad), extend it to 3D prismatic mesh.
 * and rebuild the simulation system
 */
class ExtendTo3D
{
public:

  ExtendTo3D(SimulationSystem & system, const Parser::Card & c)
  : _system(system), _card(c)
  {}

  /**
   * nothing to do
   */
  ~ExtendTo3D() {}

  /**
   * convert it here!
   */
  void operator() ();

private:

  SimulationSystem & _system;

  const Parser::Card & _card;

  /**
   * magic number of old mesh
   */
  unsigned int magic_num;

  /**
   * node number of old mesh
   */
  unsigned int n_nodes;

  /**
   * elem number of old mesh
   */
  unsigned int n_elem;

  /**
   * subdomain number of old mesh
   */
  unsigned int n_subs;

  /**
   * point of old mesh
   */
  std::vector<Point>          mesh_points;

  /**
   * subdomain label of old mesh
   */
  std::vector<int>            mesh_conn;

  /**
   * subdomain material of old mesh
   */
  std::vector<std::string>    subdomain_label;

  /**
   * elem info of old mesh
   */
  std::vector<std::string>    subdomain_material;

  /**
   * id of boundary elem in old mesh
   */
  std::vector<unsigned int>       bd_elems;

  /**
   * side of boundary elem in old mesh
   */
  std::vector<unsigned short int> bd_sides;

  /**
   * boundary id of elem-side pair
   */
  std::vector<short int>          bd_ids;

  /**
   * map boundary label to boundary id
   */
  std::map<const std::string, short int> bd_map;

  typedef std::map<const std::string, short int>::iterator Bd_It;


  typedef std::map<unsigned int, double> variable_map;

  /**
   * old system variables
   */
  std::map<std::string, std::vector<variable_map> > variables;

  /**
   * save the old 2d system
   */
  void save_old_system();

  /**
   * Packs the element \p elem at the end of vector \p conn.
   */
  void pack_element (std::vector<int> &conn, const Elem* elem) const;

  /**
   * aux function to sync solutions between different processors
   */
  void sync_solution(const std::string &sol);

  /**
   * map renumbered node id to original node id
   */
  std::map<unsigned int, unsigned int> _node_to_old_node_id_map;

  /**
   * setup new 3d system
   */
  void set_new_system();

  /**
   * setup variables in new system regions
   */
  void set_new_regions();

};

#endif

