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

//  $Id: insulator_region.cc,v 1.6 2008/07/09 05:58:16 gdiso Exp $

#include "elem.h"
#include "insulator_region.h"
#include "fvm_node_data_insulator.h"
#include "fvm_cell_data_insulator.h"


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




InsulatorSimulationRegion::InsulatorSimulationRegion(const std::string &name, const std::string &material, const double T, const double z)
:SimulationRegion(name, material, T, z)
{
  // material should be initializted after region variables
  this->set_region_variables();
  mt = new Material::MaterialInsulator(this);
}


void InsulatorSimulationRegion::insert_cell (const Elem * e)
{
  // not a local element
  if( !e->on_local() ) return;

  // insert into region element vector
  _region_cell.push_back(e);
  _region_cell_data.push_back( new FVM_Insulator_CellData(&_cell_data_storage, _region_cell_variables) );
}


void InsulatorSimulationRegion::insert_fvm_node(FVM_Node * fn)
{
    // node (or ghost node) belongs to this processor
    // we should build FVM_NodeData structure for it

  if  ( fn->root_node()->on_local() )
    fn->hold_node_data( new FVM_Insulator_NodeData(&_node_data_storage, _region_point_variables) );

  _region_node[fn->root_node()->id()] = fn;
}


void InsulatorSimulationRegion::set_region_variables()
{
  // preset named variales
  _node_data_storage.allocate_scalar_variable( std::vector<bool>(FVM_Insulator_NodeData::n_scalar(),false));
  _node_data_storage.allocate_complex_variable( std::vector<bool>(FVM_Insulator_NodeData::n_complex(),false));
  _node_data_storage.allocate_vector_variable( std::vector<bool>(FVM_Insulator_NodeData::n_vector(),false));
  _node_data_storage.allocate_tensor_variable( std::vector<bool>(FVM_Insulator_NodeData::n_tensor(),false));

  _cell_data_storage.allocate_scalar_variable( std::vector<bool>(FVM_Insulator_CellData::n_scalar(),false));
  _cell_data_storage.allocate_complex_variable( std::vector<bool>(FVM_Insulator_CellData::n_complex(),false));
  _cell_data_storage.allocate_vector_variable( std::vector<bool>(FVM_Insulator_CellData::n_vector(),false));
  _cell_data_storage.allocate_tensor_variable( std::vector<bool>(FVM_Insulator_CellData::n_tensor(),false));

  // define _region_point_variables
  _region_point_variables["potential"       ] = SimulationVariable("potential", SCALAR, POINT_CENTER, "V", FVM_Insulator_NodeData::_psi_, true);
  _region_point_variables["temperature"     ] = SimulationVariable("temperature", SCALAR, POINT_CENTER, "K", FVM_Insulator_NodeData::_T_, true);
  _region_point_variables["density"         ] = SimulationVariable("density", SCALAR, POINT_CENTER, "g/cm^3", FVM_Insulator_NodeData::_density_, true);
  _region_point_variables["affinity"        ] = SimulationVariable("affinity", SCALAR, POINT_CENTER, "eV", FVM_Insulator_NodeData::_affinity_, true);
  _region_point_variables["ec"              ] = SimulationVariable("ec", SCALAR, POINT_CENTER, "eV", FVM_Insulator_NodeData::_Ec_, true);
  _region_point_variables["ev"              ] = SimulationVariable("ev", SCALAR, POINT_CENTER, "eV", FVM_Insulator_NodeData::_Ev_, true);
  _region_point_variables["eg"              ] = SimulationVariable("eg", SCALAR, POINT_CENTER, "eV", FVM_Insulator_NodeData::_Eg_, true);
  _region_point_variables["eps"             ] = SimulationVariable("eps", SCALAR, POINT_CENTER, "C/V/m", FVM_Insulator_NodeData::_eps_, true);
  _region_point_variables["mu"              ] = SimulationVariable("mu", SCALAR, POINT_CENTER, "s^2*V/C/m", FVM_Insulator_NodeData::_mu_, true);

  _region_point_variables["optical_energy"  ] = SimulationVariable("optical_energy", SCALAR, POINT_CENTER, "eV", FVM_Insulator_NodeData::_OptE_, true);
  _region_point_variables["particle_energy" ] = SimulationVariable("particle_energy", SCALAR, POINT_CENTER, "eV", FVM_Insulator_NodeData::_PatE_, true);

  _region_point_variables["potential.last"  ] = SimulationVariable("potential.last", SCALAR, POINT_CENTER, "V", FVM_Insulator_NodeData::_psi_last_, true);
  _region_point_variables["potential.old"     ] = SimulationVariable("potential.old", SCALAR, POINT_CENTER, "V", FVM_Insulator_NodeData::_psi_old_, true);
  _region_point_variables["temperature.last"] = SimulationVariable("temperature.last", SCALAR, POINT_CENTER, "K", FVM_Insulator_NodeData::_T_last_, true);

  _region_point_variables["efield"        ] = SimulationVariable("efield", VECTOR, POINT_CENTER, "V/cm", FVM_Insulator_NodeData::_E_, true);

  _region_point_variables["potential.ac"  ] = SimulationVariable("potential.ac", COMPLEX, POINT_CENTER, "V", FVM_Insulator_NodeData::_psi_ac_, false);
  _region_point_variables["temperature.ac"] = SimulationVariable("temperature.ac", COMPLEX, POINT_CENTER, "K", FVM_Insulator_NodeData::_T_ac_, false);
  _region_point_variables["optical_efield"] = SimulationVariable("optical_efield", COMPLEX, POINT_CENTER, "V/cm", FVM_Insulator_NodeData::_OpE_complex_, false);
  _region_point_variables["optical_hfield"] = SimulationVariable("optical_hfield", COMPLEX, POINT_CENTER, "A/cm", FVM_Insulator_NodeData::_OpH_complex_, false);

  // define _region_cell_variables
  _region_cell_variables["efield"        ] = SimulationVariable("efield", VECTOR, CELL_CENTER, "V/cm", FVM_Insulator_CellData::_E_, true);
  _region_cell_variables["elec_current"  ] = SimulationVariable("elec_current", VECTOR, CELL_CENTER, "A/cm", FVM_Insulator_CellData::_Jn_, true);
  _region_cell_variables["hole_current"  ] = SimulationVariable("hole_current", VECTOR, CELL_CENTER, "A/cm", FVM_Insulator_CellData::_Jp_, true);

    // allocate variables
  std::map<std::string, SimulationVariable>::iterator it;
  for( it = _region_point_variables.begin(); it !=  _region_point_variables.end(); ++it)
    this->add_variable(it->second);
  for( it = _region_cell_variables.begin(); it !=  _region_cell_variables.end(); ++it)
    this->add_variable(it->second);
}

void InsulatorSimulationRegion::init(PetscScalar T_external)
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

    // set the initial temperature of lattice to external temperature
    node_data->T()  =  T_external;

    // set aux data for insulator
    node_data->affinity() = mt->basic->Affinity(T_external);
    node_data->Eg()       = mt->band->Eg(T_external);
    node_data->density()  = mt->basic->Density(T_external);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();

    // the initial vacuum potential is equal to electron affinity potential
    node_data->psi() = -node_data->affinity();
  }

}



void InsulatorSimulationRegion::reinit_after_import()
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

    // lattice temperature, have been read from data file!
    PetscScalar T = node_data->T();

    // set aux data for insulator
    node_data->affinity() = mt->basic->Affinity(T);
    node_data->Eg()       = mt->band->Eg(T);
    node_data->density()  = mt->basic->Density(T);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
  }

  //NOTE T and psi are read from data file!
}


void InsulatorSimulationRegion::set_pmi(const std::string &type, const std::string &model_name, std::vector<Parser::Parameter> & pmi_parameters)
{
  get_material_base()->set_pmi(type,model_name,pmi_parameters);

  local_node_iterator it = on_local_nodes_begin();
  for ( ; it!=on_local_nodes_end(); ++it)
  {
    FVM_Node * node = (*it);
    FVM_NodeData * node_data = (*it)->node_data();
    genius_assert(node_data!=NULL);
    get_material_base()->init_node(type, (*it)->root_node(), node_data);
  }

  // update buffered value as if changed by PMI
  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;
    FVM_NodeData * node_data = fvm_node->node_data();

    // map current node to material buffer
    mt->mapping(fvm_node->root_node(), node_data, 0.0);

    // lattice temperature, have been read from data file!
    PetscScalar T = node_data->T();

    // set aux data for insulator
    node_data->affinity() = mt->basic->Affinity(T);
    node_data->Eg()       = mt->band->Eg(T);
    node_data->density()  = mt->basic->Density(T);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
  }
}



