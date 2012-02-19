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

//  $Id: semiconductor_region.cc,v 1.18 2008/07/09 05:58:16 gdiso Exp $

#include "asinh.hpp" // for asinh

#include "elem.h"
#include "mesh_tools.h"
#include "boundary_info.h"
#include "semiconductor_region.h"
#include "fvm_node_data_semiconductor.h"
#include "fvm_cell_data_semiconductor.h"
#include "parallel.h"
#include "log.h"

using PhysicalUnit::cm;
using PhysicalUnit::nm;
using PhysicalUnit::m;
using PhysicalUnit::s;
using PhysicalUnit::V;
using PhysicalUnit::C;
using PhysicalUnit::K;
using PhysicalUnit::g;
using PhysicalUnit::A;
using PhysicalUnit::eV;
using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::eps0;
using PhysicalUnit::mu0;



SemiconductorSimulationRegion::SemiconductorSimulationRegion(const std::string &name, const std::string &material, const double T, const double z)
    :SimulationRegion(name, material, T, z)
{
  // material should be initializted after region variables
  this->set_region_variables();
  mt = new Material::MaterialSemiconductor(this);
}


void SemiconductorSimulationRegion::clear()
{
  SimulationRegion::clear();

  // clear previous value
  _elem_on_insulator_interface.clear();
  _elem_in_mos_channel.clear();
  _nearest_interface_normal.clear();
  _elem_touch_boundary.clear();
}

void SemiconductorSimulationRegion::insert_cell (const Elem * e)
{
  // not a local element
  if( !e->on_local() ) return;

  // insert into region element vector
  _region_cell.push_back(e);
  _region_cell_data.push_back(new FVM_Semiconductor_CellData(&_cell_data_storage, _region_cell_variables) );
}


void SemiconductorSimulationRegion::insert_fvm_node(FVM_Node * fn)
{
    // node (or ghost node) belongs to this processor
    // we should build FVM_NodeData structure for it
  if  ( fn->root_node()->on_local() )
    fn->hold_node_data( new FVM_Semiconductor_NodeData(&_node_data_storage, _region_point_variables) );

  _region_node[fn->root_node()->id()] = fn;
}


void SemiconductorSimulationRegion::set_region_variables()
{
  // preset named variales
  _node_data_storage.allocate_scalar_variable( std::vector<bool>(FVM_Semiconductor_NodeData::n_scalar(),false));
  _node_data_storage.allocate_complex_variable( std::vector<bool>(FVM_Semiconductor_NodeData::n_complex(),false));
  _node_data_storage.allocate_vector_variable( std::vector<bool>(FVM_Semiconductor_NodeData::n_vector(),false));
  _node_data_storage.allocate_tensor_variable( std::vector<bool>(FVM_Semiconductor_NodeData::n_tensor(),false));

  _cell_data_storage.allocate_scalar_variable( std::vector<bool>(FVM_Semiconductor_CellData::n_scalar(),false));
  _cell_data_storage.allocate_complex_variable( std::vector<bool>(FVM_Semiconductor_CellData::n_complex(),false));
  _cell_data_storage.allocate_vector_variable( std::vector<bool>(FVM_Semiconductor_CellData::n_vector(),false));
  _cell_data_storage.allocate_tensor_variable( std::vector<bool>(FVM_Semiconductor_CellData::n_tensor(),false));

  // define _region_point_variables
  _region_point_variables["electron"           ] = SimulationVariable("electron", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_n_, true);
  _region_point_variables["hole"               ] = SimulationVariable("hole", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_p_, true);
  _region_point_variables["potential"          ] = SimulationVariable("potential", SCALAR, POINT_CENTER, "V", FVM_Semiconductor_NodeData::_psi_, true);
  _region_point_variables["potential_correct"  ] = SimulationVariable("potential_correct", SCALAR, POINT_CENTER, "V", FVM_Semiconductor_NodeData::_dpsi_, true);
  _region_point_variables["temperature"        ] = SimulationVariable("temperature", SCALAR, POINT_CENTER, "K", FVM_Semiconductor_NodeData::_T_, true);
  _region_point_variables["elec_temperature"   ] = SimulationVariable("elec_temperature", SCALAR, POINT_CENTER, "K", FVM_Semiconductor_NodeData::_Tn_, true);
  _region_point_variables["hole_temperature"   ] = SimulationVariable("hole_temperature", SCALAR, POINT_CENTER, "K", FVM_Semiconductor_NodeData::_Tp_, true);

  _region_point_variables["density"            ] = SimulationVariable("density", SCALAR, POINT_CENTER, "g/cm^3", FVM_Semiconductor_NodeData::_density_, true);
  _region_point_variables["affinity"           ] = SimulationVariable("affinity", SCALAR, POINT_CENTER, "eV", FVM_Semiconductor_NodeData::_affinity_, true);
  _region_point_variables["ec"                 ] = SimulationVariable("ec", SCALAR, POINT_CENTER, "eV", FVM_Semiconductor_NodeData::_Ec_, true);
  _region_point_variables["ev"                 ] = SimulationVariable("ev", SCALAR, POINT_CENTER, "eV", FVM_Semiconductor_NodeData::_Ev_, true);
  _region_point_variables["eg"                 ] = SimulationVariable("eg", SCALAR, POINT_CENTER, "eV", FVM_Semiconductor_NodeData::_Eg_, true);
  _region_point_variables["qfn"                ] = SimulationVariable("qfn", SCALAR, POINT_CENTER, "eV", FVM_Semiconductor_NodeData::_qFn_, true);
  _region_point_variables["qfp"                ] = SimulationVariable("qfp", SCALAR, POINT_CENTER, "eV", FVM_Semiconductor_NodeData::_qFp_, true);

  _region_point_variables["nc"                 ] = SimulationVariable("nc", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_Nc_, true);
  _region_point_variables["nv"                 ] = SimulationVariable("nv", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_Nv_, true);

  _region_point_variables["na"                 ] = SimulationVariable("na", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_Na_, true);
  _region_point_variables["nd"                 ] = SimulationVariable("nd", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_Nd_, true);


  _region_point_variables["mole_x"             ] = SimulationVariable("mole_x", SCALAR, POINT_CENTER, "1", FVM_Semiconductor_NodeData::_mole_x_, true);
  _region_point_variables["mole_y"             ] = SimulationVariable("mole_y", SCALAR, POINT_CENTER, "1", FVM_Semiconductor_NodeData::_mole_y_, true);

  _region_point_variables["mun"                ] = SimulationVariable("mun", SCALAR, POINT_CENTER, "cm^2/V/s", FVM_Semiconductor_NodeData::_mun_, true);
  _region_point_variables["mup"                ] = SimulationVariable("mup", SCALAR, POINT_CENTER, "cm^2/V/s", FVM_Semiconductor_NodeData::_mup_, true);

  _region_point_variables["recombination"      ] = SimulationVariable("recombination", SCALAR, POINT_CENTER, "cm^-3/s", FVM_Semiconductor_NodeData::_Recomb_, true);
  _region_point_variables["recombination_dir"  ] = SimulationVariable("recombination_dir", SCALAR, POINT_CENTER, "cm^-3/s", FVM_Semiconductor_NodeData::_Recomb_Dir_, true);
  _region_point_variables["recombination_srh"  ] = SimulationVariable("recombination_srh", SCALAR, POINT_CENTER, "cm^-3/s", FVM_Semiconductor_NodeData::_Recomb_SRH_, true);
  _region_point_variables["recombination_auger"] = SimulationVariable("recombination_auger", SCALAR, POINT_CENTER, "cm^-3/s", FVM_Semiconductor_NodeData::_Recomb_Auger_, true);

  _region_point_variables["impact_ionization"  ] = SimulationVariable("impact_ionization", SCALAR, POINT_CENTER, "cm^-3/s", FVM_Semiconductor_NodeData::_ImpactIonization_, true);

  _region_point_variables["field_generation"   ] = SimulationVariable("field_generation", SCALAR, POINT_CENTER, "cm^-3/s", FVM_Semiconductor_NodeData::_Field_G_, true);
  _region_point_variables["optical_generation" ] = SimulationVariable("optical_generation", SCALAR, POINT_CENTER, "cm^-3/s", FVM_Semiconductor_NodeData::_OptG_, true);
  _region_point_variables["optical_heat"       ] = SimulationVariable("optical_heat", SCALAR, POINT_CENTER, "eV/s", FVM_Semiconductor_NodeData::_OptQ_, true);
  _region_point_variables["optical_energy"     ] = SimulationVariable("optical_energy", SCALAR, POINT_CENTER, "eV", FVM_Semiconductor_NodeData::_OptE_, true);
  _region_point_variables["particle_generation"] = SimulationVariable("particle_generation", SCALAR, POINT_CENTER, "cm^-3/s", FVM_Semiconductor_NodeData::_PatG_, true);
  _region_point_variables["particle_energy"    ] = SimulationVariable("particle_energy", SCALAR, POINT_CENTER, "eV", FVM_Semiconductor_NodeData::_PatE_, true);

  _region_point_variables["elec_injection"     ] = SimulationVariable("elec_injection", SCALAR, POINT_CENTER, "A", FVM_Semiconductor_NodeData::_EIn_, true);
  _region_point_variables["hole_injection"     ] = SimulationVariable("hole_injection", SCALAR, POINT_CENTER, "A", FVM_Semiconductor_NodeData::_HIn_, true);

  _region_point_variables["eps"                ] = SimulationVariable("eps", SCALAR, POINT_CENTER, "C/V/m", FVM_Semiconductor_NodeData::_eps_, true);
  _region_point_variables["mu"                 ] = SimulationVariable("mu", SCALAR, POINT_CENTER, "s^2*V/C/m", FVM_Semiconductor_NodeData::_mu_, true);

  _region_point_variables["electron.last"      ] = SimulationVariable("electron.last", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_n_last_, true);
  _region_point_variables["hole.last"          ] = SimulationVariable("hole.last", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_p_last_, true);
  _region_point_variables["potential.last"     ] = SimulationVariable("potential.last", SCALAR, POINT_CENTER, "V", FVM_Semiconductor_NodeData::_psi_last_, true);
  _region_point_variables["temperature.last"   ] = SimulationVariable("temperature.last", SCALAR, POINT_CENTER, "K", FVM_Semiconductor_NodeData::_T_last_, true);
  _region_point_variables["elec_temperature.last"] = SimulationVariable("elec_temperature.last", SCALAR, POINT_CENTER, "K", FVM_Semiconductor_NodeData::_Tn_last_, true);
  _region_point_variables["hole_temperature.last"] = SimulationVariable("hole_temperature.last", SCALAR, POINT_CENTER, "K", FVM_Semiconductor_NodeData::_Tp_last_, true);
  _region_point_variables["electron.old"      ] = SimulationVariable("electron.old", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_n_old_, true);
  _region_point_variables["hole.old"          ] = SimulationVariable("hole.old", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_p_old_, true);
  _region_point_variables["potential.old"     ] = SimulationVariable("potential.old", SCALAR, POINT_CENTER, "V", FVM_Semiconductor_NodeData::_psi_old_, true);

  _region_point_variables["charge_density"     ] = SimulationVariable("charge_density", SCALAR, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_rho_, true);

  _region_point_variables["efield"             ] = SimulationVariable("efield", VECTOR, POINT_CENTER, "V/cm", FVM_Semiconductor_NodeData::_E_, true);
  _region_point_variables["elec_current"       ] = SimulationVariable("elec_current", VECTOR, POINT_CENTER, "A/cm", FVM_Semiconductor_NodeData::_Jn_, true);
  _region_point_variables["hole_current"       ] = SimulationVariable("hole_current", VECTOR, POINT_CENTER, "A/cm", FVM_Semiconductor_NodeData::_Jp_, true);

  _region_point_variables["electron.ac"        ] = SimulationVariable("electron.ac", COMPLEX, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_n_ac_, false);
  _region_point_variables["hole.ac"            ] = SimulationVariable("hole.ac", COMPLEX, POINT_CENTER, "cm^-3", FVM_Semiconductor_NodeData::_p_ac_, false);
  _region_point_variables["potential.ac"       ] = SimulationVariable("potential.ac", COMPLEX, POINT_CENTER, "V", FVM_Semiconductor_NodeData::_psi_ac_, false);
  _region_point_variables["temperature.ac"     ] = SimulationVariable("temperature.ac", COMPLEX, POINT_CENTER, "K", FVM_Semiconductor_NodeData::_T_ac_, false);
  _region_point_variables["elec_temperature.ac"] = SimulationVariable("elec_temperature.ac", COMPLEX, POINT_CENTER, "K", FVM_Semiconductor_NodeData::_Tn_ac_, false);
  _region_point_variables["hole_temperature.ac"] = SimulationVariable("hole_temperature.ac", COMPLEX, POINT_CENTER, "K", FVM_Semiconductor_NodeData::_Tp_ac_, false);

  _region_point_variables["optical_efield"     ] = SimulationVariable("optical_efield", COMPLEX, POINT_CENTER, "V/cm", FVM_Semiconductor_NodeData::_OpE_complex_, false);
  _region_point_variables["optical_hfield"     ] = SimulationVariable("optical_hfield", COMPLEX, POINT_CENTER, "A/cm", FVM_Semiconductor_NodeData::_OpH_complex_, false);

  // define _region_cell_variables
  _region_cell_variables["efield"        ] = SimulationVariable("efield", VECTOR, CELL_CENTER, "V/cm", FVM_Semiconductor_CellData::_E_, true);
  _region_cell_variables["elec_current"  ] = SimulationVariable("elec_current", VECTOR, CELL_CENTER, "A/cm", FVM_Semiconductor_CellData::_Jn_, true);
  _region_cell_variables["hole_current"  ] = SimulationVariable("hole_current", VECTOR, CELL_CENTER, "A/cm", FVM_Semiconductor_CellData::_Jp_, true);


  // allocate variables
  std::map<std::string, SimulationVariable>::iterator it;
  for( it = _region_point_variables.begin(); it !=  _region_point_variables.end(); ++it)
    this->add_variable(it->second);
  for( it = _region_cell_variables.begin(); it !=  _region_cell_variables.end(); ++it)
    this->add_variable(it->second);

}




void SemiconductorSimulationRegion::init(PetscScalar T_external)
{
  //init FVM_NodeData
  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;
    FVM_NodeData * node_data = fvm_node->node_data();

    // set the initial temperature of lattice, electron and hole to external temperature
    node_data->T()  =  T_external;
    node_data->Tn() =  T_external;
    node_data->Tp() =  T_external;

    // map current node to material buffer
    mt->mapping(fvm_node->root_node(), node_data, 0.0);

    // we can get some parameter only related with temperature
    PetscScalar ni  = mt->band->ni(T_external);
    PetscScalar Nc  = mt->band->Nc(T_external);
    PetscScalar Nv  = mt->band->Nv(T_external);

    // get net doping concentration
    PetscScalar Net_Doping = node_data->Net_doping();

    // compute electron and hole concentration for thermal equilibrium state
    if( Net_Doping > 0 )       //n-type
    {
      node_data->n_last() = node_data->n() = (Net_Doping + std::sqrt(std::pow(Net_Doping,2) + 4*ni*ni))/2;
      node_data->p_last() = node_data->p() = ni*ni/node_data->n();
    }
    else                       //p-type
    {
      node_data->p_last() = node_data->p() = (-Net_Doping + std::sqrt(std::pow(Net_Doping,2) + 4*ni*ni))/2;
      node_data->n_last() = node_data->n() = ni*ni/node_data->p();
    }

    // init more physical parameters
    node_data->affinity() = mt->basic->Affinity(T_external);
    node_data->density()  = mt->basic->Density(T_external);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
    node_data->Eg()       = mt->band->Eg(T_external);
    node_data->Nc()       = mt->band->Nc(T_external);
    node_data->Nv()       = mt->band->Nv(T_external);
    node_data->mun()      = mt->mob->ElecMob(node_data->p(), node_data->n(), T_external, 0, 0, T_external);
    node_data->mup()      = mt->mob->HoleMob(node_data->p(), node_data->n(), T_external, 0, 0, T_external);

    // the electrostatic potential is initializted as vacuum level
    node_data->psi() = kb*T_external/e*boost::math::asinh(Net_Doping/(2*ni))
                       - node_data->Eg()/(2*e)
                       - kb*T_external/(2*e)*log(Nc/Nv)
                       - node_data->affinity()
                       - mt->band->EgNarrowToEc(node_data->p(), node_data->n(), T_external);

    node_data->Ec() = -(e*node_data->psi() + node_data->affinity());
    node_data->Ev() = -(e*node_data->psi() + node_data->affinity()
                        - mt->band->EgNarrowToEv(node_data->p(), node_data->n(), T_external)
                        + node_data->Eg()
                       );
    /*
    // the quantum conduction band is initializted as conduction band
    node_data->Eqc() = -(e*node_data->psi()
                         + node_data->affinity()
                         + mt->band->EgNarrowToEc(node_data->p(), node_data->n(), T_external)
                         + kb*T_external*log(node_data->Nc()));

    // the quantum valence band is initializted as valence band
    node_data->Eqv() = -(e*node_data->psi()
                         + node_data->affinity()
                         - mt->band->EgNarrowToEv(node_data->p(), node_data->n(), T_external)
                         - kb*T_external*log(node_data->Nv())
                         + node_data->Eg());
    */
  }

  // build data structure for insulator interface
  find_elem_on_insulator_interface();
  find_nearest_interface_normal();
  find_elem_touch_boundary();
  check_fvm_cell();
}



void SemiconductorSimulationRegion::reinit_after_import()
{
  //init FVM_NodeData
  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;
    FVM_NodeData * node_data = fvm_node->node_data();

    // lattice temperature, have been read from data file!
    PetscScalar T = node_data->T();

    // map current node to material buffer
    mt->mapping(fvm_node->root_node(), node_data, 0.0);

    // init more physical parameters
    node_data->n_last()   = node_data->n();
    node_data->p_last()   = node_data->p();
    node_data->affinity() = mt->basic->Affinity(T);
    node_data->density()  = mt->basic->Density(T);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
    node_data->Eg()       = mt->band->Eg(T);
    node_data->Nc()       = mt->band->Nc(T);
    node_data->Nv()       = mt->band->Nv(T);
    node_data->mun()      = mt->mob->ElecMob(node_data->p(), node_data->n(), T, 0, 0, T);
    node_data->mup()      = mt->mob->HoleMob(node_data->p(), node_data->n(), T, 0, 0, T);
  }

  // build data structure for insulator interface
  find_elem_on_insulator_interface();
  find_nearest_interface_normal();
  find_elem_touch_boundary();
  check_fvm_cell();
}



void SemiconductorSimulationRegion::set_pmi(const std::string &type, const std::string &model_name, std::vector<Parser::Parameter> & pmi_parameters)
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

    // lattice temperature, have been read from data file!
    PetscScalar T = node_data->T();

    // map current node to material buffer
    mt->mapping(fvm_node->root_node(), node_data, 0.0);

    PetscScalar  d_affinity =  mt->basic->Affinity(T) - node_data->affinity();

    node_data->affinity() = mt->basic->Affinity(T);
    node_data->density()  = mt->basic->Density(T);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
    node_data->Eg()       = mt->band->Eg(T);
    node_data->Nc()       = mt->band->Nc(T);
    node_data->Nv()       = mt->band->Nv(T);
    node_data->mun()      = mt->mob->ElecMob(node_data->p(), node_data->n(), T, 0, 0, T);
    node_data->mup()      = mt->mob->HoleMob(node_data->p(), node_data->n(), T, 0, 0, T);

    node_data->psi() = node_data->psi() - d_affinity;

    node_data->Ec() = -(e*node_data->psi() + node_data->affinity());
    node_data->Ev() = -(e*node_data->psi() + node_data->affinity()
                      - mt->band->EgNarrowToEv(node_data->p(), node_data->n(), T)
                      + node_data->Eg() );
  }
}


void SemiconductorSimulationRegion::find_elem_on_insulator_interface()
{
  const_element_iterator it = elements_begin();
  const_element_iterator it_end = elements_end();
  for(; it!=it_end; ++it)
  {
    const Elem * elem = (*it);
    if(!elem->on_interface()) continue;

    for(unsigned int n=0; n<elem->n_sides(); ++n)
    {
      const Elem * neighbor = elem->neighbor(n);
      if(neighbor && elem->subdomain_id()!=neighbor->subdomain_id())
      {
        genius_assert(neighbor->on_local());
        SimulationRegion * region = _subdomain_id_to_region_map.find(neighbor->subdomain_id())->second;
        if(region->type()==InsulatorRegion)
        {
          std::pair<unsigned int, SimulationRegion *> val = std::make_pair(n, region);
          _elem_on_insulator_interface.insert(std::make_pair(elem->id(), val));
        }
      }
    }
  }
}


void SemiconductorSimulationRegion::find_elem_touch_boundary()
{
  std::set<const Node *> boundary_node;

  const_element_iterator it = elements_begin();
  const_element_iterator it_end = elements_end();
  for(; it!=it_end; ++it)
  {
    const Elem * elem = (*it);
    if(elem->on_boundary() || elem->on_interface())
    {
      for(unsigned int n=0; n<elem->n_nodes(); ++n)
        boundary_node.insert(elem->get_node(n));
    }
  }

  for(it = elements_begin(); it!=it_end; ++it)
  {
    const Elem * elem = (*it);
    for(unsigned int n=0; n<elem->n_nodes(); ++n)
      if( boundary_node.find(elem->get_node(n)) != boundary_node.end())
        _elem_touch_boundary.insert(elem);
  }

}




void SemiconductorSimulationRegion::find_nearest_interface_normal()
{



}


bool SemiconductorSimulationRegion::is_elem_on_insulator_interface(const Elem *elem) const
{ return _elem_on_insulator_interface.find(elem->id())!=_elem_on_insulator_interface.end(); }


void SemiconductorSimulationRegion::elem_on_insulator_interface(const Elem *elem,
    std::vector<unsigned int> & sides, std::vector<SimulationRegion *> &regions) const
{
  typedef _multimap_elem_on_interface_type::const_iterator It;
  std::pair<It, It> pos = _elem_on_insulator_interface.equal_range(elem->id());
  while(pos.first!=pos.second)
  {
    sides.push_back(pos.first->second.first);
    regions.push_back(pos.first->second.second);
    ++pos.first;
  }
}


bool SemiconductorSimulationRegion::highfield_mobility() const
{
  return ( _advanced_model.HighFieldMobility || _advanced_model.ImpactIonization || _advanced_model.BandBandTunneling) ;
}


Real SemiconductorSimulationRegion::truncated_partial_area(const Elem * elem, unsigned int ne) const
{
  //return std::abs(elem->partial_area_with_edge(ne));
#if 1
  Real partial_area =  elem->partial_area_with_edge_truncated(ne);

  std::pair<unsigned int, unsigned int> edge_nodes;
  elem->nodes_on_edge(ne, edge_nodes);
  const FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes.first);   // fvm_node of node1
  const FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes.second);  // fvm_node of node2

  unsigned int dim = elem->dim();
  Real min_area = 0;
  if( dim == 2)
    min_area = 0.1*nm;
  if( dim == 3)
    min_area = 0.01*nm*nm;

  // however, zero truncated partial area on surface will cut the current path, we here set an average "partial area"
  if(fvm_n1->cv_surface_area(fvm_n2->root_node()) < min_area)
    partial_area = elem->partial_area_with_edge_average();
  if(partial_area<min_area && (fvm_n1->boundary_id() != BoundaryInfo::invalid_id && fvm_n2->boundary_id() != BoundaryInfo::invalid_id) )
    partial_area = elem->partial_area_with_edge_average();

  return std::max(min_area, partial_area);
#endif
}



void SemiconductorSimulationRegion::zero_node_current()
{
  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;
    FVM_NodeData * node_data = fvm_node->node_data();

    node_data->Jn() = 0.0;
    node_data->Jp() = 0.0;
  }
}


bool SemiconductorSimulationRegion::check_fvm_cell() const
{
  std::map<const FVM_Node *, std::map<const FVM_Node *, double> > truncated_partial_area_map;

  const_element_iterator elem_it = elements_begin();
  const_element_iterator elem_it_end = elements_end();
  for(; elem_it!=elem_it_end; ++elem_it)
  {
    const Elem * elem = *elem_it;

    for(unsigned int ne=0; ne<elem->n_edges(); ++ne )
    {
      std::pair<unsigned int, unsigned int> edge_nodes;
      elem->nodes_on_edge(ne, edge_nodes);
      const FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes.first);   // fvm_node of node1
      const FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes.second);  // fvm_node of node2

      double truncated_partial_area =  this->truncated_partial_area(elem, ne);;
      truncated_partial_area_map[fvm_n1][fvm_n2] += truncated_partial_area;
      truncated_partial_area_map[fvm_n2][fvm_n1] += truncated_partial_area;
    }
  }

  Real fvm_cell_quality = 1.0;
  std::map<const FVM_Node *, std::map<const FVM_Node *, double> >::const_iterator fvm_it = truncated_partial_area_map.begin();
  for( ; fvm_it != truncated_partial_area_map.end(); ++fvm_it)
  {
    const FVM_Node * fvm_node = fvm_it->first;
    const std::map<const FVM_Node *, double> & truncated_partial_area_with_neighbor = fvm_it->second;

    unsigned int nonzero_partial_area = 0;
    double total_area = 0.0;
    double min_area   = 1e100;
    double max_area   = 0.0;
    std::map<const FVM_Node *, double>::const_iterator it = truncated_partial_area_with_neighbor.begin();
    for(; it != truncated_partial_area_with_neighbor.end(); ++it)
    {
      double area = it->second;
      if(area > 0.0)
      {
        nonzero_partial_area++;
        total_area += area;
        min_area = std::min(min_area, area);
        max_area = std::max(max_area, area);
      }
    }

    if(nonzero_partial_area < 2)
      std::cerr << "Warning: Bad voronoi cell find within region " << this->name() << std::endl;

    fvm_cell_quality = std::min(fvm_cell_quality, min_area/max_area);
  }

  Parallel::min(fvm_cell_quality);

  if(fvm_cell_quality < 1e-4)
  {
    MESSAGE << "Warning: Voronoi cell quality " << fvm_cell_quality << " (the ratio of min/max surface area) of region " << this->name() << " is poor. "
            << "It may cause inaccurate result or even convergence difficults. Try to improve mesh quality if possible!" << std::endl;
    RECORD();
  }
  return true;

}


