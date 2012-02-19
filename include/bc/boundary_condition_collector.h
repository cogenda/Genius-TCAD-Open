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

//  $Id: boundary_condition_collector.h,v 1.22 2008/07/09 05:58:16 gdiso Exp $

#ifndef __bc_collector_h__
#define __bc_collector_h__

#include <ctype.h>
#include <string.h>


#include "parser.h"
#include "physical_unit.h"
#include "boundary_condition.h"
#include "simulation_region.h"
#include "mxml.h"

/**
 * Store and manage all the Boundary Conditions
 */
class BoundaryConditionCollector
{
public:

  /**
   *  constructor
   */
  BoundaryConditionCollector(SimulationSystem & system, Parser::InputParser & decks);

  /**
   *  destructor
   */
  ~BoundaryConditionCollector()
  {
    this->clear();
  }

  /**
   * clear all
   */
  void clear()
  {
    for(size_t nbc=0; nbc<_bcs.size(); nbc++ )
      delete _bcs[nbc];
    _bcs.clear();
    _bc_label_to_bc.clear();
    _bd_id_to_bc_index.clear();
    _bc_index_to_bd_id.clear();
  }

  /**
   * @return the boundary condition number it holds
   */
  unsigned int n_bcs() const
  { return  _bcs.size(); }

  /**
   * @return the electrode boundary condition number it holds
   */
  unsigned int n_electrode_bcs() const
  {
    unsigned int n_electrode=0;
    for(unsigned int n=0; n<_bcs.size(); ++n)
      if(_bcs[n]->is_electrode()) n_electrode++;
    return  n_electrode;
  }

  /**
   * add bc to boundary condition collector
   */
  void add_bc(BoundaryCondition * bc)
  {
    _bcs.push_back(bc);
    _bc_label_to_bc[bc->label()] = bc;
  }

  /**
   * @return the const bc pointer by its index
   */
  const BoundaryCondition * get_bc(unsigned int n) const
    { genius_assert(n<n_bcs()); return  _bcs[n];}

  /**
   * @return the writable bc pointer by its index
   */
  BoundaryCondition * get_bc(unsigned int n)
  { genius_assert(n<n_bcs()); return  _bcs[n];}

  /**
   * @return the const bc pointer by boundary id
   */
  const BoundaryCondition * get_bc_by_bd_id(short int id) const
    {  return  _bcs[ get_bc_index_by_bd_id(id) ];}

  /**
   * @return the writable bc pointer by boundary id
   */
  BoundaryCondition * get_bc_by_bd_id(short int id)
  {  return  _bcs[ get_bc_index_by_bd_id(id) ];}


  /**
   * @return the const bc pointer by its label
   */
  const BoundaryCondition * get_bc(const std::string & label) const
  {
    if(  _bc_label_to_bc.find(label) !=_bc_label_to_bc.end() )
      return  (*_bc_label_to_bc.find(label)).second ;
    return NULL;
  }

  /**
   * @return the writable bc pointer by its label
   */
  BoundaryCondition * get_bc(const std::string & label)
  {
    if(  _bc_label_to_bc.find(label) !=_bc_label_to_bc.end() )
      return  (*_bc_label_to_bc.find(label)).second ;
    return NULL;
  }

  /**
   * set the bc index and boundary id mapping
   */
  void build_bc_index_to_bd_id_map(const std::set<short int> & bd_id)
  {
    std::set<short int>::const_iterator it = bd_id.begin();
    for(unsigned int index=0; it!=bd_id.end(); ++index, ++it )
    {
      _bc_index_to_bd_id[index] = *it;
      _bd_id_to_bc_index[*it]   = index;
    }
  }


  /**
   * @return the writable bc pointer by its case insensitivity label
   */
  BoundaryCondition * get_bc_nocase(const std::string & label);

  /**
   * @return the writable bc pointer by its case insensitivity label
   */
  const BoundaryCondition * get_bc_nocase(const std::string & label) const;


  /**
   * find the bc by its case insensitivity label and assign vapp to it
   */
  void set_vapp_nocase(const std::string & label, double vapp)
  {
    BoundaryCondition * bc = get_bc_nocase(label);
    genius_assert( bc != NULL );
    genius_assert( bc->is_electrode() );
    bc->ext_circuit()->Vapp() = vapp;
  }


  /**
   * @return bc index by boundary id stored in mesh.boundary_info
   */
  unsigned int get_bc_index_by_bd_id(short int bd_id) const;

  /**
   * @return boundary id by bc index
   */
  short int get_bd_id_by_bc_index(unsigned int bc_index) const;

  /**
   * set up the boundary condition here
   */
  int bc_setup();

  /**
   * set up the boundary nodes for PMI
   */
  void pmi_init_bc(const std::string &region_label, const std::string &type);

  /**
   * force all the electrode with zero Vapp and Iapp
   */
  void all_electrode_ground()
  {
    for(unsigned n=0; n<_bcs.size(); n++)
      if( _bcs[n]->is_electrode() )
        _bcs[n]->ext_circuit()->ground();
  }

  /**
   * export all the boundary condition to a file
   */
  void export_boundary_condition(const std::string& filename) const;

  /**
   * get all terminal voltage/current as an dom element
   */
  mxml_node_t* get_dom_terminal_info() const;

private:

  /**
   * the reference to input deck
   */
  Parser::InputParser & _decks;

  /**
   * the reference to corresponding SimulationSystem
   * since bc may contain physical equations, it is important
   * for bc to access region information
   */
  SimulationSystem    & _system;

  /**
   * the writable reference to mesh,
   * we need to reset boundary condition index to boundary_info
   */
  MeshBase & _mesh;

  /**
   * create resistive metal regions instead of electrode region
   */
  bool _resistive_metal_mode;


  /**
   *  the boundary condition array
   */
  std::vector<BoundaryCondition *>   _bcs;

  /**
   *  map boundary label to boundary condition
   */
  std::map<std::string, BoundaryCondition * >  _bc_label_to_bc;

  /**
   * map boundary id to boundary condition index
   */
  std::map<short int, unsigned int> _bd_id_to_bc_index;

  /**
   * map boundary condition index to boundary id
   */
  std::map<unsigned int, short int> _bc_index_to_bd_id;

  /**
   * @return all the interfaces between region r1 and region r2
   */
  std::vector<std::string> _find_interface_by_2_regions(const std::string &r1, const std::string &r2) const;

  /**
   * get bc index from card
   */
  unsigned int get_bc_from_card (const Parser::Card &c, std::string & bc_label);

  // private function, setting each bcs

  int Set_BC_NeumannBoundary(const Parser::Card &c);

  int Set_BC_OhmicContact(const Parser::Card &c);

  int Set_BC_IF_Metal_Ohmic(const Parser::Card &c);

  int Set_BC_SchottkyContact(const Parser::Card &c);

  int Set_BC_IF_Metal_Schottky(const Parser::Card &c);

  int Set_BC_GateContact(const Parser::Card &c);

  int Set_BC_SimpleGateContact(const Parser::Card &c);

  int Set_BC_InsulatorInterface(const Parser::Card &c);

  int Set_BC_HeteroInterface(const Parser::Card &c);

  int Set_BC_Solderpad(const Parser::Card &c);

  int Set_BC_ChargedContact(const Parser::Card &c);


  int Set_Electrode_OhmicContact(const Parser::Card &c);

  int Set_Electrode_SchottkyContact(const Parser::Card &c);

  int Set_Electrode_GateContact(const Parser::Card &c);

  int Set_SetFloatMetal(const Parser::Card &c);

  int Set_InterConnect(const Parser::Card &c);

  int Set_Charge(const Parser::Card &c);

  int Set_BC_AbsorbingBoundary(const Parser::Card &c);
};

#endif //#ifndef __bc_collector_h__


