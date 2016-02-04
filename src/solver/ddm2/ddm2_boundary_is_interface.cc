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




#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_is.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;



///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * do pre-process to function for DDM2 solver
 */
void InsulatorSemiconductorInterfaceBC::DDM2_Function_Preprocess(PetscScalar * ,Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // Insulator-Semiconductor interface at Semiconductor side, do nothing
          case SemiconductorRegion:  break;

          // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
          case InsulatorRegion:
          {
            // record the source row and dst row
            src_row.push_back(fvm_nodes[i]->global_offset());
            dst_row.push_back(fvm_nodes[0]->global_offset());
            clear_row.push_back(fvm_nodes[i]->global_offset());

            src_row.push_back(fvm_nodes[i]->global_offset()+1);
            dst_row.push_back(fvm_nodes[0]->global_offset()+3);
            clear_row.push_back(fvm_nodes[i]->global_offset()+1);
            break;
          }
          case VacuumRegion:
          break;

          default: genius_error(); //we should never reach here
      }
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void InsulatorSemiconductorInterfaceBC::DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // Insulator-Semiconductor interface is processed here

  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // buffer for Vec location
  std::vector<PetscInt> iy;
  iy.reserve(2*n_nodes());
  // buffer for Vec value
  std::vector<PetscScalar> y;
  y.reserve(2*n_nodes());

  const PetscScalar qf = this->scalar("qf");

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // Insulator-Semiconductor interface at Semiconductor side, do nothing
          case SemiconductorRegion:
          {
            // semiconductor region should be the first region
            genius_assert(i==0);

            // process interface traps
            SemiconductorSimulationRegion * sregion = (SemiconductorSimulationRegion *) regions[i];

            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();
            genius_assert(node_data);

            PetscScalar n   =  x[fvm_nodes[i]->local_offset()+1];                         // electron density
            PetscScalar p   =  x[fvm_nodes[i]->local_offset()+2];                         // hole density
            PetscScalar T   =  x[fvm_nodes[i]->local_offset()+3];

            // process interface fixed charge density
            PetscScalar boundary_area = fvm_nodes[i]->outside_boundary_surface_area();
            VecSetValue(f, fvm_nodes[i]->global_offset(), (qf+node_data->interface_charge())*boundary_area, ADD_VALUES);

            {
              // surface recombination
              Material::MaterialSemiconductor *mt =  sregion->material();
              mt->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);
              mt->band->nie(p, n, T);
              PetscScalar GSurf = - mt->band->R_Surf(p, n, T) * boundary_area; //generation due to SRH
              PetscScalar HR = -GSurf*(node_data->Eg() + 3*PhysicalUnit::kb*T); // heat transferred to lattice

              VecSetValue(f, fvm_nodes[i]->global_offset()+1, GSurf, ADD_VALUES);
              VecSetValue(f, fvm_nodes[i]->global_offset()+2, GSurf, ADD_VALUES);
              VecSetValue(f, fvm_nodes[i]->global_offset()+3, HR, ADD_VALUES);
            }

            if (sregion->get_advanced_model()->Trap)
            {
              // process interface traps
              sregion->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

              // calculate interface trap occupancy
              PetscScalar ni = sregion->material()->band->nie(p, n, T);
              sregion->material()->trap->Calculate(false,p,n,ni,T);

              // contribution of trapped charge to Poisson's eqn
              PetscScalar TrappedC = sregion->material()->trap->Charge(false) * boundary_area;
              if (TrappedC !=0)
                VecSetValue(f, fvm_nodes[i]->global_offset(), TrappedC, ADD_VALUES);

              // electron/hole capture rates and contribution to continuity equations
              PetscScalar TrapElec = sregion->material()->trap->ElectronTrapRate(false,n,ni,T) * boundary_area;
              PetscScalar TrapHole = sregion->material()->trap->HoleTrapRate    (false,p,ni,T) * boundary_area;
              if (TrapElec != 0)
                VecSetValue(f, fvm_nodes[i]->global_offset()+1, -TrapElec, ADD_VALUES);
              if (TrapHole != 0)
                VecSetValue(f, fvm_nodes[i]->global_offset()+2, -TrapHole, ADD_VALUES);

              PetscScalar EcEi = 0.5*node_data->Eg() - kb*T*log(node_data->Nc()/node_data->Nv());
              PetscScalar EiEv = 0.5*node_data->Eg() + kb*T*log(node_data->Nc()/node_data->Nv());
              PetscScalar H    = sregion->material()->trap->TrapHeat(false,p,n,ni,T,T,T,EcEi,EiEv);
              VecSetValue(f, fvm_nodes[i]->global_offset()+3, H*boundary_area, ADD_VALUES);
            }
            break;
          }

          // Insulator-Semiconductor interface at Insulator side, we get the previous f for insulator region,
          // as well as surface charge density, plus to corresponding function of semiconductor. and then
          // force the potential equal to corresponding point in semiconductor
          case InsulatorRegion:
          {
            genius_assert(fvm_nodes[i]->root_node()->processor_id() == fvm_nodes[0]->root_node()->processor_id() );

            // the governing equation of this fvm node

            iy.push_back(fvm_nodes[i]->global_offset()+0);
            // psi of this node
            PetscScalar V = x[fvm_nodes[i]->local_offset()+0];
            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            //genius_assert( regions[0]->type()==SemiconductorRegion );
            PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+0];
            // the psi of this node is equal to corresponding psi of semiconductor node
            // since psi should be continuous for the interface
            PetscScalar ff1 = V - V_semi;
            y.push_back(ff1);

            iy.push_back(fvm_nodes[i]->global_offset()+1);
            // T of this node
            PetscScalar T = x[fvm_nodes[i]->local_offset()+1];
            // T for corresponding semiconductor region
            PetscScalar T_semi = x[fvm_nodes[0]->local_offset()+3];
            // the T of this node is equal to corresponding T of semiconductor node
            // by assuming no heat resistance between 2 region
            PetscScalar ff2 = T - T_semi;

            y.push_back(ff2);

            break;
          }
          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }
  }


  // set new value to row
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y[0]), ADD_VALUES);

  add_value_flag = ADD_VALUES;


  // gate current
  _gate_current_function(x, f, add_value_flag);

}






/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML2 solver
 */
void InsulatorSemiconductorInterfaceBC::DDM2_Jacobian_Preprocess(PetscScalar *,SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // Insulator-Semiconductor interface at Semiconductor side, do nothing
          case SemiconductorRegion:  break;

          // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
          case InsulatorRegion:
          {
            // record the source row and dst row
            src_row.push_back(fvm_nodes[i]->global_offset());
            dst_row.push_back(fvm_nodes[0]->global_offset());
            clear_row.push_back(fvm_nodes[i]->global_offset());

            src_row.push_back(fvm_nodes[i]->global_offset()+1);
            dst_row.push_back(fvm_nodes[0]->global_offset()+3);
            clear_row.push_back(fvm_nodes[i]->global_offset()+1);
            break;
          }
          case VacuumRegion:
          break;

          default: genius_error(); //we should never reach here
      }
    }

  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void InsulatorSemiconductorInterfaceBC::DDM2_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // after that, set values to source rows
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // Insulator-Semiconductor interface at Semiconductor side, do nothing
          case SemiconductorRegion:
          {
            //the indepedent variable number, 4 for each node
            adtl::AutoDScalar::numdir = 4;

            // process interface traps
            SemiconductorSimulationRegion * sregion = (SemiconductorSimulationRegion *) regions[i];

            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();
            genius_assert(node_data);

            {
              PetscInt index[4] = {fvm_nodes[i]->global_offset()+0, fvm_nodes[i]->global_offset()+1, fvm_nodes[i]->global_offset()+2, fvm_nodes[i]->global_offset()+3};
              AutoDScalar n   =  x[fvm_nodes[i]->local_offset()+1];   n.setADValue(1, 1.0);              // electron density
              AutoDScalar p   =  x[fvm_nodes[i]->local_offset()+2];   p.setADValue(2, 1.0);              // hole density
              AutoDScalar T   =  x[fvm_nodes[i]->local_offset()+3];   T.setADValue(3, 1.0);
              PetscScalar boundary_area = fvm_nodes[i]->outside_boundary_surface_area();

              Material::MaterialSemiconductor *mt =  sregion->material();
              //synchronize with material database
              mt->set_ad_num(adtl::AutoDScalar::numdir);


              mt->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

              AutoDScalar ni = mt->band->nie(p, n, T);

              {
                // surface recombination
                AutoDScalar GSurf = - mt->band->R_Surf(p, n, T) * boundary_area;
                AutoDScalar HR = -GSurf*(node_data->Eg() + 3*PhysicalUnit::kb*T); // heat transferred to lattice

                jac->add_row(  index[1],  4,  &index[0],  GSurf.getADValue() );
                jac->add_row(  index[2],  4,  &index[0],  GSurf.getADValue() );
                jac->add_row(  index[3],  4,  &index[0],  HR.getADValue() );
              }

              if (sregion->get_advanced_model()->Trap)
              {

                mt->trap->Calculate(false,p,n,ni,T);
                AutoDScalar TrappedC = mt->trap->ChargeAD(false) * boundary_area;
                jac->add_row(  index[0],  4,  &index[0],  TrappedC.getADValue() );

                AutoDScalar GElec = - mt->trap->ElectronTrapRate(false,n,ni,T) * boundary_area;
                AutoDScalar GHole = - mt->trap->HoleTrapRate    (false,p,ni,T) * boundary_area;
                jac->add_row(  index[1],  4,  &index[0],  GElec.getADValue() );
                jac->add_row(  index[2],  4,  &index[0],  GHole.getADValue() );

                AutoDScalar EcEi = 0.5*node_data->Eg() - kb*T*log(node_data->Nc()/node_data->Nv());
                AutoDScalar EiEv = 0.5*node_data->Eg() + kb*T*log(node_data->Nc()/node_data->Nv());
                AutoDScalar H = mt->trap->TrapHeat(false,p,n,ni,T,T,T,EcEi,EiEv);
                jac->add_row(  index[3],  4,  &index[0],  (H*boundary_area).getADValue() );

              }
            }
            break;
          }
          // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
          case InsulatorRegion:
          {

            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            //genius_assert( regions[0]->type()==SemiconductorRegion );
            AutoDScalar  V_semi = x[fvm_nodes[0]->local_offset()]; V_semi.setADValue(1,1.0);

            // the psi of this node is equal to corresponding psi of semiconductor node
            AutoDScalar  ff1 = V - V_semi;

            // set Jacobian of governing equation ff
            jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[i]->global_offset(),  ff1.getADValue(0) );
            jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[0]->global_offset(),  ff1.getADValue(1) );


            // T of this node
            AutoDScalar  T = x[fvm_nodes[i]->local_offset()+1]; T.setADValue(0,1.0);

            // T for corresponding semiconductor region
            AutoDScalar  T_semi = x[fvm_nodes[0]->local_offset()+3]; T_semi.setADValue(1,1.0);

            // the T of this node is equal to corresponding T of semiconductor node
            // we assuming T is continuous for the interface
            AutoDScalar ff2 = T - T_semi;

            // set Jacobian of governing equation ff2
            jac->add( fvm_nodes[i]->global_offset()+1,  fvm_nodes[i]->global_offset()+1,  ff2.getADValue(0) );
            jac->add( fvm_nodes[i]->global_offset()+1,  fvm_nodes[0]->global_offset()+3,  ff2.getADValue(1) );

            break;

          }
          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

  // gate current
  _gate_current_jacobian(x, jac, add_value_flag);

}

