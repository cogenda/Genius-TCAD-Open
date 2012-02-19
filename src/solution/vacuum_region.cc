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


#include "elem.h"
#include "vacuum_region.h"
#include "fvm_node_data_vacuum.h"
#include "fvm_cell_data_vacuum.h"

using PhysicalUnit::cm;
using PhysicalUnit::m;
using PhysicalUnit::s;
using PhysicalUnit::V;
using PhysicalUnit::C;
using PhysicalUnit::K;
using PhysicalUnit::g;
using PhysicalUnit::A;
using PhysicalUnit::eV;
using PhysicalUnit::eps0;
using PhysicalUnit::mu0;




VacuumSimulationRegion::VacuumSimulationRegion(const std::string &name, const std::string &material, const double T, const double z)
:SimulationRegion(name, material, T, z)
{
  this->set_region_variables();
  mt = new Material::MaterialVacuum(this);
}


void VacuumSimulationRegion::insert_cell (const Elem * e)
{
  // not a local element
  if( !e->on_local() ) return;

  // insert into region element vector
  _region_cell.push_back(e);
  _region_cell_data.push_back(new FVM_Vacuum_CellData(&_cell_data_storage, _region_cell_variables) );
}


void VacuumSimulationRegion::insert_fvm_node(FVM_Node * fn)
{
  // node (or ghost node) belongs to this processor
  // we should build FVM_NodeData structure for it

  if  ( fn->root_node()->on_local() )
    fn->hold_node_data( new FVM_Vacuum_NodeData(&_node_data_storage, _region_point_variables) );

  _region_node[fn->root_node()->id()] = fn;
}


void VacuumSimulationRegion::set_region_variables()
{
  // preset named variales
  _node_data_storage.allocate_scalar_variable( std::vector<bool>(FVM_Vacuum_NodeData::n_scalar(),false));
  _node_data_storage.allocate_complex_variable( std::vector<bool>(FVM_Vacuum_NodeData::n_complex(),false));
  _node_data_storage.allocate_vector_variable( std::vector<bool>(FVM_Vacuum_NodeData::n_vector(),false));
  _node_data_storage.allocate_tensor_variable( std::vector<bool>(FVM_Vacuum_NodeData::n_tensor(),false));

  _cell_data_storage.allocate_scalar_variable( std::vector<bool>(FVM_Vacuum_CellData::n_scalar(),false));
  _cell_data_storage.allocate_complex_variable( std::vector<bool>(FVM_Vacuum_CellData::n_complex(),false));
  _cell_data_storage.allocate_vector_variable( std::vector<bool>(FVM_Vacuum_CellData::n_vector(),false));
  _cell_data_storage.allocate_tensor_variable( std::vector<bool>(FVM_Vacuum_CellData::n_tensor(),false));

  // define _region_point_variables
  _region_point_variables["potential"     ] = SimulationVariable("potential", SCALAR, POINT_CENTER, "V", FVM_Vacuum_NodeData::_psi_, true);
  _region_point_variables["density"       ] = SimulationVariable("density", SCALAR, POINT_CENTER, "g/cm^3", FVM_Vacuum_NodeData::_density_, true);
  _region_point_variables["affinity"      ] = SimulationVariable("affinity", SCALAR, POINT_CENTER, "eV", FVM_Vacuum_NodeData::_affinity_, true);
  _region_point_variables["eps"           ] = SimulationVariable("eps", SCALAR, POINT_CENTER, "C/V/m", FVM_Vacuum_NodeData::_eps_, true);
  _region_point_variables["mu"            ] = SimulationVariable("mu", SCALAR, POINT_CENTER, "s^2*V/C/m", FVM_Vacuum_NodeData::_mu_, true);
  _region_point_variables["potential.last"] = SimulationVariable("potential.last", SCALAR, POINT_CENTER, "V", FVM_Vacuum_NodeData::_psi_last_, true);

  _region_point_variables["efield"        ] = SimulationVariable("efield", VECTOR, POINT_CENTER, "V/cm", FVM_Vacuum_NodeData::_E_, true);

  _region_point_variables["optical_efield"] = SimulationVariable("optical_efield", COMPLEX, POINT_CENTER, "V/cm", FVM_Vacuum_NodeData::_OpE_complex_, false);
  _region_point_variables["optical_hfield"] = SimulationVariable("optical_hfield", COMPLEX, POINT_CENTER, "A/cm", FVM_Vacuum_NodeData::_OpH_complex_, false);

  // define _region_cell_variables
  _region_cell_variables["efield"        ] = SimulationVariable("efield", VECTOR, CELL_CENTER, "V/cm", FVM_Vacuum_CellData::_E_, true);

  // allocate variables
  std::map<std::string, SimulationVariable>::iterator it;
  for( it = _region_point_variables.begin(); it !=  _region_point_variables.end(); ++it)
    this->add_variable(it->second);
  for( it = _region_cell_variables.begin(); it !=  _region_cell_variables.end(); ++it)
    this->add_variable(it->second);
}



void VacuumSimulationRegion::init(PetscScalar T_external)
{

  //init FVM_NodeData
  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;
    FVM_NodeData * node_data = fvm_node->node_data();

    // map current node to material buffer
    mt->mapping(fvm_node->root_node(), node_data, 0.0);

    // set aux data for insulator
    node_data->affinity() = mt->basic->Affinity(T_external);
    node_data->density()  = mt->basic->Density(T_external);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();

    // the initial vacuum potential is equal to electron affinity potential
    node_data->psi() = -node_data->affinity();
  }

}



void VacuumSimulationRegion::reinit_after_import()
{
  //init FVM_NodeData
  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;
    FVM_NodeData * node_data = fvm_node->node_data();

    // map current node to material buffer
    mt->mapping(fvm_node->root_node(), node_data, 0.0);

    // lattice temperature.
    PetscScalar T = 1.0;

    // set aux data for insulator
    node_data->affinity() = mt->basic->Affinity(T);
    node_data->density()  = mt->basic->Density(T);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
  }
}

