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


// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_neumann.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;





///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void NeumannBC::DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Neumann boundary condition is processed here
  bool surface_recomb = this->flag("surface.recombination");
  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");


  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }


  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const SimulationRegion * region = (*region_node_begin(*node_it)).second.first;
    const FVM_Node * fvm_node = (*region_node_begin(*node_it)).second.second;

    switch ( region->type() )
    {
    case SemiconductorRegion:
      {
        SemiconductorSimulationRegion * sregion = (SemiconductorSimulationRegion *) region;

        PetscScalar n   =  x[fvm_node->local_offset()+1];                         // electron density
        PetscScalar p   =  x[fvm_node->local_offset()+2];                         // hole density
        // process governing equation of T, which should consider heat exchange to environment
        PetscScalar T = x[fvm_node->local_offset()+3];  // lattice temperature

        PetscScalar boundary_area = fvm_node->outside_boundary_surface_area();

        if(surface_recomb)
        {
          // surface recombination
          Material::MaterialSemiconductor *mt =  sregion->material();
          mt->mapping(fvm_node->root_node(), fvm_node->node_data(), SolverSpecify::clock);
          PetscScalar GSurf = - mt->band->R_Surf(p, n, T) * boundary_area; //generation due to SRH

          VecSetValue(f, fvm_node->global_offset()+1, GSurf, ADD_VALUES);
          VecSetValue(f, fvm_node->global_offset()+2, GSurf, ADD_VALUES);
        }

        // add heat flux out of Neumann boundary to lattice temperature equatiuon
        PetscScalar fT = Heat_Transfer*(T_external()-T)*boundary_area;
        VecSetValue(f, fvm_node->global_offset()+3, fT, ADD_VALUES);

        break;
      }
    case MetalRegion     :
    case ElectrodeRegion :
    case InsulatorRegion :
      {
        // process governing equation of T, which should consider heat exchange to environment
        PetscScalar T = x[fvm_node->local_offset()+1];  // lattice temperature
        // add heat flux out of Neumann boundary to lattice temperature equatiuon
        PetscScalar boundary_area  = fvm_node->outside_boundary_surface_area();
        PetscScalar fT = Heat_Transfer*(T_external()-T)*boundary_area;
        VecSetValue(f, fvm_node->global_offset()+1, fT, ADD_VALUES);

        break;
      }
    case VacuumRegion:
      break;
    default: genius_error(); //we should never reach here
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}








/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void NeumannBC::DDM2_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // Neumann boundary condition is processed here
  bool surface_recomb = this->flag("surface.recombination");
  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");


  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const SimulationRegion * region = (*region_node_begin(*node_it)).second.first;
    const FVM_Node * fvm_node = (*region_node_begin(*node_it)).second.second;

    switch ( region->type() )
    {
    case SemiconductorRegion:
      {
        // process interface traps
        SemiconductorSimulationRegion * sregion = (SemiconductorSimulationRegion *) region;

        //the indepedent variable number, 4 for each node
        adtl::AutoDScalar::numdir = 4;
        PetscInt index[4] = {fvm_node->global_offset()+0,
                             fvm_node->global_offset()+1,
                             fvm_node->global_offset()+2,
                             fvm_node->global_offset()+3};

        AutoDScalar n   = x[fvm_node->local_offset()+1];   n.setADValue(1, 1.0);              // electron density
        AutoDScalar p   = x[fvm_node->local_offset()+2];   p.setADValue(2, 1.0);              // hole density
        // process governing equation of T, which should consider heat exchange to environment
        AutoDScalar T   = x[fvm_node->local_offset()+3];   T.setADValue(3, 1.0); // lattice temperature

        PetscScalar boundary_area = fvm_node->outside_boundary_surface_area();

        Material::MaterialSemiconductor *mt =  sregion->material();
        mt->mapping(fvm_node->root_node(), fvm_node->node_data(), SolverSpecify::clock);
        //synchronize with material database
        mt->set_ad_num(adtl::AutoDScalar::numdir);

        if(surface_recomb)
        {
          // surface recombination
          AutoDScalar GSurf = - mt->band->R_Surf(p, n, T) * boundary_area;
          jac->add_row(  index[1],  4,  &index[0],  GSurf.getADValue() );
          jac->add_row(  index[2],  4,  &index[0],  GSurf.getADValue() );
        }

        // add heat flux out of Neumann boundary to lattice temperature equatiuon
        AutoDScalar fT = Heat_Transfer*(T_external()-T)*boundary_area;
        jac->add( index[3],  index[3],  fT.getADValue(3) );

        break;
      }
    case MetalRegion     :
    case ElectrodeRegion :
    case InsulatorRegion :
      {
        //the indepedent variable number, we only need 1 here.
        adtl::AutoDScalar::numdir=1;

        // process governing equation of T, which should consider heat exchange to environment
        AutoDScalar T = x[fvm_node->local_offset()+1];  T.setADValue(0, 1.0); // lattice temperature

        // add heat flux out of Neumann boundary to lattice temperature equatiuon
        PetscScalar boundary_area  = fvm_node->outside_boundary_surface_area();
        AutoDScalar fT = Heat_Transfer*(T_external()-T)*boundary_area;

        jac->add( fvm_node->global_offset()+1,  fvm_node->global_offset()+1,  fT.getADValue(0) );

        break;
      }
    case VacuumRegion:
      break;
    default: genius_error(); //we should never reach here
    }


  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}

