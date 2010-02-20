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
#include "simulation_system.h"
#include "semiconductor_region.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::eps0;
using PhysicalUnit::mu0;



SemiconductorSimulationRegion::SemiconductorSimulationRegion(const std::string &name, const std::string &material, SimulationSystem & system)
    :SimulationRegion(name, material, system), mt( new Material::MaterialSemiconductor(material, name) )
{
}


void SemiconductorSimulationRegion::insert_cell (const Elem * e)
{
    // not a local element
  if( !e->on_local() ) return;

    // insert into region element vector
  _region_cell.push_back(e);
  _region_cell_data.push_back(new FVM_Semiconductor_CellData);
}


void SemiconductorSimulationRegion::init(PetscScalar T_external)
{
  //init FVM_NodeData
  node_iterator it= _region_node.begin();

  for ( ; it!=_region_node.end(); ++it)
  {

    // when node_data is not NULL, this node belongs to current process
    // or at least it is a ghost node
    if((*it).second->node_data() == NULL) continue;

    FVM_NodeData * node_data = (*it).second->node_data();

    // set the initial temperature of lattice, electron and hole to external temperature
    node_data->T()  =  T_external;
    node_data->Tn() =  T_external;
    node_data->Tp() =  T_external;

    // map current node to material buffer
    mt->mapping((*it).second->root_node(), node_data, 0.0);

    // we can get some parameter only related with temperature
    PetscScalar ni  = mt->band->ni(T_external);
    PetscScalar Nc  = mt->band->Nc(T_external);
    PetscScalar Nv  = mt->band->Nv(T_external);

    // get net doping concentration
    PetscScalar Net_Doping = node_data->Net_doping();

    // compute electron and hole concentration for thermal equilibrium state
    if( Net_Doping > 0 )       //n-type
    {
      node_data->n() = (Net_Doping + std::sqrt(std::pow(Net_Doping,2) + 4*ni*ni))/2;
      node_data->p() = ni*ni/node_data->n();
    }
    else                       //p-type
    {
      node_data->p() = (-Net_Doping + std::sqrt(std::pow(Net_Doping,2) + 4*ni*ni))/2;
      node_data->n() = ni*ni/node_data->p();
    }

    // init more physical parameters
    node_data->affinity() = mt->basic->Affinity(T_external);
    node_data->density()  = mt->basic->Density(T_external);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
    node_data->Eg()       = mt->band->Eg(T_external);
    node_data->Nc()       = mt->band->Nc(T_external);
    node_data->Nv()       = mt->band->Nv(T_external);


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

  }


  find_elem_on_insulator_interface();
  find_nearest_interface_normal();
}



void SemiconductorSimulationRegion::reinit_after_import()
{
  //init FVM_NodeData
  node_iterator it= _region_node.begin();

  for ( ; it!=_region_node.end(); ++it)
  {

    // when node_data is not NULL, this node belongs to current process
    // or at least it is a ghost node
    if((*it).second->node_data() == NULL) continue;

    FVM_NodeData * node_data = (*it).second->node_data();

    // lattice temperature, have been read from data file!
    PetscScalar T = node_data->T();

    // map current node to material buffer
    mt->mapping((*it).second->root_node(), node_data, 0.0);

    // init more physical parameters
    node_data->affinity() = mt->basic->Affinity(T);
    node_data->density()  = mt->basic->Density(T);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
    node_data->Eg()       = mt->band->Eg(T);
    node_data->Nc()       = mt->band->Nc(T);
    node_data->Nv()       = mt->band->Nv(T);
  }

  find_elem_on_insulator_interface();
  find_nearest_interface_normal();
}


void SemiconductorSimulationRegion::find_elem_on_insulator_interface()
{
  // clear previous value
  _elem_on_insulator_interface.clear();

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
        SimulationRegion * region = _subdomain_id_to_region_map.find(neighbor->subdomain_id())->second;
        if(region->type()==InsulatorRegion)
        {
          std::pair<unsigned int, SimulationRegion *> val = std::make_pair(n, region);
          _elem_on_insulator_interface.insert(std::make_pair(elem, val));
        }
      }
    }
  }
}


void SemiconductorSimulationRegion::find_nearest_interface_normal()
{
  // clear previous value
  _nearest_interface_normal.clear();


}



void SemiconductorSimulationRegion::zero_node_current()
{
  node_iterator it = nodes_begin();
  node_iterator it_end = nodes_end();
  for(; it!=it_end; ++it)
  {
    FVM_Node * fvm_node = (*it).second;

    // NOTE: here, solution for all the local node should be updated!
    if( !fvm_node->root_node()->on_local() ) continue;

    FVM_NodeData * node_data = fvm_node->node_data();

    node_data->Jn() = 0.0;
    node_data->Jp() = 0.0;
  }
}

