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
#include "simulation_system.h"
#include "insulator_region.h"

using PhysicalUnit::eps0;
using PhysicalUnit::mu0;




InsulatorSimulationRegion::InsulatorSimulationRegion(const std::string &name, const std::string &material, SimulationSystem & system)
:SimulationRegion(name, material, system) , mt( new Material::MaterialInsulator(material, name) )
{}


void InsulatorSimulationRegion::insert_cell (const Elem * e)
{
    // not a local element
  if( !e->on_local() ) return;

    // insert into region element vector
  _region_cell.push_back(e);
  _region_cell_data.push_back(new FVM_Insulator_CellData);
}


void InsulatorSimulationRegion::init(PetscScalar T_external)
{

  //init FVM_NodeData
  node_iterator it= _region_node.begin();

  for ( ; it!=_region_node.end(); ++it)
  {

    // when node_data is not NULL, this node belongs to current process
    // or at least it is a ghost node
    if((*it).second->node_data() == NULL) continue;

    FVM_Insulator_NodeData * node_data = dynamic_cast<FVM_Insulator_NodeData *> ((*it).second->node_data());

    // map current node to material buffer
    mt->mapping((*it).second->root_node(), node_data, 0.0);

    // set the initial temperature of lattice to external temperature
    node_data->T()  =  T_external;

    // set aux data for insulator
    node_data->affinity() = mt->basic->Affinity(T_external);
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
  node_iterator it= _region_node.begin();

  for ( ; it!=_region_node.end(); ++it)
  {

    // when node_data is not NULL, this node belongs to current process
    // or at least it is a ghost node
    if((*it).second->node_data() == NULL) continue;

    FVM_Insulator_NodeData * node_data = dynamic_cast<FVM_Insulator_NodeData *> ((*it).second->node_data());

    // map current node to material buffer
    mt->mapping((*it).second->root_node(), node_data, 0.0);

    // lattice temperature, have been read from data file!
    PetscScalar T = node_data->T();

    // set aux data for insulator
    node_data->affinity() = mt->basic->Affinity(T);
    node_data->Eg()       = mt->basic->Eg(T);
    node_data->density()  = mt->basic->Density(T);
    node_data->eps()      = eps0*mt->basic->Permittivity();
    node_data->mu()       = mu0*mt->basic->Permeability();
  }

  //NOTE T and psi are read from data file!


}

