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

// C++ includes
#include <queue>
#include <algorithm>


#include "elem.h"
#include "resistance_region.h"
#include "fvm_node_data_resistance.h"
#include "fvm_cell_data_resistance.h"
#include "boundary_info.h"
#include "boundary_condition.h"

#include "parallel.h"


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


// static variable
double MetalSimulationRegion::_aux_resistance = 1e100;
double MetalSimulationRegion::_aux_capacitance = 0;

void MetalSimulationRegion::set_aux_parasitic_parameter(double res, double cap)
{
  _aux_resistance = res;
  _aux_capacitance = cap;
}


MetalSimulationRegion::MetalSimulationRegion(const std::string &name, const std::string &material, const double T, const unsigned int dim, const double z)
  :SimulationRegion(name, material, T, dim, z)
{
  // material should be initializted after region variables
  this->set_region_variables();
  mt = new Material::MaterialConductor(this);

}


void MetalSimulationRegion::insert_cell (const Elem * e)
{
  // not a local element
  if( !e->on_local() ) return;

  // insert into region element vector
  _region_cell.push_back(e);
  _region_cell_data.push_back(new FVM_Resistance_CellData(&_cell_data_storage, _region_cell_variables));
}


void MetalSimulationRegion::insert_fvm_node(FVM_Node * fn)
{

  // node (or ghost node) belongs to this processor
  // we should build FVM_NodeData structure for it

  if  ( fn->root_node()->on_local() )
    fn->hold_node_data( new FVM_Resistance_NodeData(&_node_data_storage, _region_point_variables) );

  _region_node[fn->root_node()->id()] = fn;
}


void MetalSimulationRegion::set_region_variables()
{
  // preset named variales
  _node_data_storage.allocate_scalar_variable( std::vector<bool>(FVM_Resistance_NodeData::n_scalar(),false));
  _node_data_storage.allocate_complex_variable( std::vector<bool>(FVM_Resistance_NodeData::n_complex(),false));
  _node_data_storage.allocate_vector_variable( std::vector<bool>(FVM_Resistance_NodeData::n_vector(),false));
  _node_data_storage.allocate_tensor_variable( std::vector<bool>(FVM_Resistance_NodeData::n_tensor(),false));

  _cell_data_storage.allocate_scalar_variable( std::vector<bool>(FVM_Resistance_CellData::n_scalar(),false));
  _cell_data_storage.allocate_complex_variable( std::vector<bool>(FVM_Resistance_CellData::n_complex(),false));
  _cell_data_storage.allocate_vector_variable( std::vector<bool>(FVM_Resistance_CellData::n_vector(),false));
  _cell_data_storage.allocate_tensor_variable( std::vector<bool>(FVM_Resistance_CellData::n_tensor(),false));

  // define _region_point_variables
  _region_point_variables["dmin"            ] = SimulationVariable("dmin", SCALAR, POINT_CENTER, "1", FVM_Resistance_NodeData::_dmin_, true);
  _region_point_variables["electron"        ] = SimulationVariable("electron", SCALAR, POINT_CENTER, "cm^-3", FVM_Resistance_NodeData::_n_, true);
  _region_point_variables["potential"       ] = SimulationVariable("potential", SCALAR, POINT_CENTER, "V", FVM_Resistance_NodeData::_psi_, true);
  _region_point_variables["temperature"     ] = SimulationVariable("temperature", SCALAR, POINT_CENTER, "K", FVM_Resistance_NodeData::_T_, true);
  _region_point_variables["density"         ] = SimulationVariable("density", SCALAR, POINT_CENTER, "g/cm^3", FVM_Resistance_NodeData::_density_, true);
  _region_point_variables["affinity"        ] = SimulationVariable("affinity", SCALAR, POINT_CENTER, "eV", FVM_Resistance_NodeData::_affinity_, true);
  _region_point_variables["Ec"              ] = SimulationVariable("Ec", SCALAR, POINT_CENTER, "eV", FVM_Resistance_NodeData::_Ec_, true);
  _region_point_variables["Ev"              ] = SimulationVariable("Ev", SCALAR, POINT_CENTER, "eV", FVM_Resistance_NodeData::_Ev_, true);
  _region_point_variables["eps"             ] = SimulationVariable("eps", SCALAR, POINT_CENTER, "C/V/m", FVM_Resistance_NodeData::_eps_, true);
  _region_point_variables["mu"              ] = SimulationVariable("mu", SCALAR, POINT_CENTER, "s^2*V/C/m", FVM_Resistance_NodeData::_mu_, true);

  _region_point_variables["optical_energy"  ] = SimulationVariable("optical_energy", SCALAR, POINT_CENTER, "eV", FVM_Resistance_NodeData::_OptE_, true);
  _region_point_variables["particle_energy" ] = SimulationVariable("particle_energy", SCALAR, POINT_CENTER, "eV", FVM_Resistance_NodeData::_PatE_, true);

  _region_point_variables["electron.last"   ] = SimulationVariable("electron.last", SCALAR, POINT_CENTER, "cm^-3", FVM_Resistance_NodeData::_n_last_, true);
  _region_point_variables["potential.last"  ] = SimulationVariable("potential.last", SCALAR, POINT_CENTER, "V", FVM_Resistance_NodeData::_psi_last_, true);
  _region_point_variables["potential.old"   ] = SimulationVariable("potential.old", SCALAR, POINT_CENTER, "V", FVM_Resistance_NodeData::_psi_old_, true);
  _region_point_variables["temperature.last"] = SimulationVariable("temperature.last", SCALAR, POINT_CENTER, "K", FVM_Resistance_NodeData::_T_last_, true);

  _region_point_variables["efield"        ] = SimulationVariable("efield", VECTOR, POINT_CENTER, "V/cm", FVM_Resistance_NodeData::_E_, true);

  _region_point_variables["potential.ac"  ] = SimulationVariable("potential.ac", COMPLEX, POINT_CENTER, "V", FVM_Resistance_NodeData::_psi_ac_, false);
  _region_point_variables["temperature.ac"] = SimulationVariable("temperature.ac", COMPLEX, POINT_CENTER, "K", FVM_Resistance_NodeData::_T_ac_, false);
  _region_point_variables["optical_efield"] = SimulationVariable("optical_efield", COMPLEX, POINT_CENTER, "V/cm", FVM_Resistance_NodeData::_OpE_complex_, false);
  _region_point_variables["optical_hfield"] = SimulationVariable("optical_hfield", COMPLEX, POINT_CENTER, "A/cm", FVM_Resistance_NodeData::_OpH_complex_, false);

  // define _region_cell_variables
  _region_cell_variables["efield"        ] = SimulationVariable("efield", VECTOR, CELL_CENTER, "V/cm", FVM_Resistance_CellData::_E_, true);

  // allocate variables
  std::map<std::string, SimulationVariable>::iterator it;
  for( it = _region_point_variables.begin(); it !=  _region_point_variables.end(); ++it)
    this->add_variable(it->second);
  for( it = _region_cell_variables.begin(); it !=  _region_cell_variables.end(); ++it)
    this->add_variable(it->second);

}


void MetalSimulationRegion::init(PetscScalar T_external)
{
  _conductance = mt->basic->Conductance();

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

    // set aux data for metal
    node_data->n()        = mt->basic->IonDensity(T_external);
    node_data->affinity() = mt->basic->Affinity(T_external);
    node_data->psi()      = -mt->basic->Affinity(T_external);
    node_data->density()  = mt->basic->Density(T_external);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();


    // the initial vacuum potential is equal to electron affinity potential
    node_data->psi() = -node_data->affinity();
  }

  find_total_nodes_in_connected_resistance_region();
}


void MetalSimulationRegion::reinit_after_import()
{

  _conductance = mt->basic->Conductance();

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
    PetscScalar T =node_data->T();

    // set aux data for metal
    node_data->n()        = mt->basic->IonDensity(T);
    node_data->affinity() = mt->basic->Affinity(T);
    node_data->density()  = mt->basic->Density(T);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();

  }


  //NOTE T and psi are read from data file!

  find_total_nodes_in_connected_resistance_region();
}


void MetalSimulationRegion::set_pmi(const std::string &type, const std::string &model_name, std::vector<Parser::Parameter> & pmi_parameters)
{
  get_material_base()->set_pmi(type,model_name,pmi_parameters);

  _conductance = mt->basic->Conductance();

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
    PetscScalar T =node_data->T();

    // set aux data for metal
    node_data->affinity() = mt->basic->Affinity(T);
    node_data->density()  = mt->basic->Density(T);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
  }
}



void MetalSimulationRegion::find_total_nodes_in_connected_resistance_region()
{
  std::set<const SimulationRegion *> metal_regions;

  // find metal regions connected to me
  {
    std::queue<const SimulationRegion *> Q;
    std::vector<bool> visit_flag(_subdomain_id_to_region_map.size(), false);

    Q.push( this );
    while(!Q.empty())
    {
      const SimulationRegion * region = Q.front();
      Q.pop();

      visit_flag[region->subdomain_id()] = true;

      if(region->type() == MetalRegion)
        metal_regions.insert(region);

      const std::vector<SimulationRegion *> & region_neighbors = region->region_neighbors();
      for(unsigned int n=0; n<region_neighbors.size(); ++n)
      {
        const SimulationRegion * region_neighbor = region_neighbors[n];
        if( visit_flag[region_neighbor->subdomain_id()] == false && region_neighbor->type() == MetalRegion)
          Q.push(region_neighbor);
      }
    }
  }

  // we should at least find one (this region)
  genius_assert( metal_regions.size() );

  _total_nodes_in_connected_resistance_region = 0;
  std::set<const SimulationRegion *>::const_iterator region_it = metal_regions.begin();
  for(; region_it != metal_regions.end(); ++region_it)
  {
    _total_nodes_in_connected_resistance_region += (*region_it)->n_node();
  }
}


void MetalSimulationRegion::find_low_resistance_solderpad()
{
  std::set<const SimulationRegion *> metal_regions;

  // find metal regions connected to me
  {
    std::queue<const SimulationRegion *> Q;
    std::vector<bool> visit_flag(_subdomain_id_to_region_map.size(), false);

    Q.push( this );
    while(!Q.empty())
    {
      const SimulationRegion * region = Q.front();
      Q.pop();

      visit_flag[region->subdomain_id()] = true;

      if(region->type() == MetalRegion)
        metal_regions.insert(region);

      const std::vector<SimulationRegion *> & region_neighbors = region->region_neighbors();
      for(unsigned int n=0; n<region_neighbors.size(); ++n)
      {
        const SimulationRegion * region_neighbor = region_neighbors[n];
        if( visit_flag[region_neighbor->subdomain_id()] == false && region_neighbor->type() == MetalRegion)
          Q.push(region_neighbor);
      }
    }
  }

  // we should at least find one (this region)
  genius_assert( metal_regions.size() );


  _connect_to_low_resistance_solderpad = false;
  std::set<const SimulationRegion *>::const_iterator region_it = metal_regions.begin();
  for(; region_it != metal_regions.end(); ++region_it)
  {
    const SimulationRegion * region = *region_it;
    std::map<short int, BoundaryCondition *>::const_iterator bc_it = region->region_boundaries().begin();
    for(; bc_it != region->region_boundaries().end(); ++bc_it)
    {
      const BoundaryCondition * bc = bc_it->second;
      if( bc->bc_type() == SolderPad )
      {
        // SolderPad with low resistance
        if(bc->ext_circuit()->is_voltage_driven()==true && bc->ext_circuit()->serial_resistance() < 1e3*V/A)
        {  _connect_to_low_resistance_solderpad = true; return ; }
      }
    }
  }

}




