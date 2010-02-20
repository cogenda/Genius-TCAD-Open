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
#include <numeric>

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition.h"
#include "parallel.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML1 solver
 */
void SchottkyContactBC::Mix_DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Schottky boundary condition is processed here

  // note, we will use INSERT_VALUES to set values of vec f
  // if the previous operator is not insert_VALUES, we should assembly the vec
  if( (add_value_flag != INSERT_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // the electrode potential in current iteration
  PetscScalar Ve = ext_circuit()->Vapp();

  // data buffer for ADD_VALUES
  std::vector<int> ix;
  std::vector<PetscScalar> y;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();
  std::vector<double> current_buffer;

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
        // Semiconductor Region of course owns Schottky Contact BC
      case SemiconductorRegion:
        {
          // semiconductor region should be the first region
          genius_assert(i==0);

          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

          // mapping this node to material library
          semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

          PetscScalar V = x[fvm_nodes[i]->local_offset()+0];  // psi of this node
          PetscScalar n = x[fvm_nodes[i]->local_offset()+1];  // electron density
          PetscScalar p = x[fvm_nodes[i]->local_offset()+2];  // hole density
          PetscScalar T = T_external();

          //Schotty Barrier Lowerring
          PetscScalar deltaVB = semi_region->material()->band->SchottyBarrierLowerring(node_data->eps(), node_data->E().size());
          //Schottky current
          PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
          PetscScalar Fn = semi_region->material()->band->SchottyJsn(n, T, this->Work_Function() - node_data->affinity() - deltaVB) * S;
          PetscScalar Fp = semi_region->material()->band->SchottyJsp(p, T, this->Work_Function() - node_data->affinity() + deltaVB) * S;

          PetscScalar ff = V + this->Work_Function() - deltaVB - Ve;

          // set governing equation to boundary condition of poisson's equation
          VecSetValue(f, fvm_nodes[i]->global_offset()+0, ff, INSERT_VALUES);
          // add Schottky current to continuity equations
          // we buffer this operator
          ix.push_back(fvm_nodes[i]->global_offset()+1);
          y.push_back( Fn);
          ix.push_back(fvm_nodes[i]->global_offset()+2);
          y.push_back(-Fp);

          // compute the current flow out of schottky electrode

          // schottky thermal emit current
          current_buffer.push_back(-(Fn+Fp)*current_scale);

          // displacement current
          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
          for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
          {
            FVM_Node *nb_node = (*nb_it).second;
            FVM_NodeData * nb_node_data = nb_node->node_data();
            // the psi of neighbor node
            PetscScalar V_nb = x[nb_node->local_offset()+0];
            // distance from nb node to this node
            PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
            // area of out surface of control volume related with neighbor node
            PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node->root_node());
            PetscScalar dEdt;
            if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false) //second order
            {
              PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
              dEdt = ( (2-r)/(1-r)*(V-V_nb)
                       - 1.0/(r*(1-r))*(node_data->psi()-nb_node_data->psi())
                       + (1-r)/r*(node_data->psi_last()-nb_node_data->psi_last()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
            }
            else//first order
            {
              dEdt = ((V-V_nb)-(node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
            }

            current_buffer.push_back( cv_boundary*node_data->eps()*dEdt*current_scale );
          }

          break;
        }
        // conductor region which has an interface with Schottky Contact boundary to semiconductor region
      case ConductorRegion:
        {
          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()];

          // since the region is sorted, we know region[0] is semiconductor region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
          //genius_assert( regions[0]->type()==SemiconductorRegion );
          PetscScalar V_semi = x[fvm_nodes[0]->local_offset()];

          // the psi of this node is equal to corresponding psi of semiconductor node
          PetscScalar ff = V - V_semi;
          // set governing equation to function vector
          VecSetValue(f, fvm_nodes[i]->global_offset(), ff, INSERT_VALUES);

          break;
        }
        // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
        // the boundary for the corner point may be Schottky. (not a nice behavier)
      case InsulatorRegion:
        {
          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()];

          // since the region is sorted, we know region[0] is semiconductor region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
          //genius_assert( regions[0]->type()==SemiconductorRegion );
          PetscScalar V_semi = x[fvm_nodes[0]->local_offset()];

          // the psi of this node is equal to corresponding psi of semiconductor node
          PetscScalar ff = V - V_semi;
          // set governing equation to function vector
          VecSetValue(f, fvm_nodes[i]->global_offset(), ff, INSERT_VALUES);

          break;
        }
      case VacuumRegion:
        break;
      default: genius_error(); //we should never reach here
      }
    }

  }


  // process buffered ADD_VALUES operator
  VecAssemblyBegin(f);
  VecAssemblyEnd(f);

  if(ix.size()) VecSetValues(f, ix.size(), &ix[0], &y[0], ADD_VALUES) ;

  // we gather the electrode current here
  Parallel::allgather(current_buffer);

  // for get the current, we must sum all the terms in current_buffer
  PetscScalar current = std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  // save the IV of current iteration
  ext_circuit()->current() = current;
  ext_circuit()->potential() = Ve;


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for Mixed DDML1 solver
 */
void SchottkyContactBC::Mix_DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

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

      case SemiconductorRegion:
        {
          break;
        }
      case ConductorRegion:
        {
          // insert none zero pattern
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), 0, ADD_VALUES);

          break;
        }
      case InsulatorRegion:
        {
          // insert none zero pattern
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), 0, ADD_VALUES);

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

}






/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML1 solver
 */
void SchottkyContactBC::Mix_DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of Schottky boundary condition is processed here
  // we use AD again. no matter it is overkill here.

  PetscScalar T = T_external();
  // the electrode potential in current iteration
  PetscScalar Ve = ext_circuit()->Vapp();

  // first, we zero all the rows corresponding to poisson's equation of Schottky bc.
  {
    // the indicator which rows should be set to zero
    std::vector<PetscInt> id;
    id.reserve(n_nodes());

    //note! MatZeroRows should be excuted on all the processor
    //no matter whether it owns this row!
    MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);

    BoundaryCondition::const_node_iterator node_it = nodes_begin();
    BoundaryCondition::const_node_iterator end_it = nodes_end();

    for(; node_it!=end_it; ++node_it )
    {

      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      // search all the fvm_node which has *node_it as root node
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

      // should clear all the rows related with this boundary condition
      for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
      {
        const FVM_Node *  fvm_node = (*rnode_it).second.second ;

        PetscInt row = fvm_node->global_offset();
        id.push_back(row);
      }
    }

    // for efficient resion, we separate MatAssemblyBegin and MatAssemblyEnd
    MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);

    MatZeroRows(*jac, id.size(), id.empty() ? NULL : &id[0], 0.0);
  }

  // after the clean, we should do schottky boundary process here

  // loop again
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
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
        // Semiconductor Region of course owns Schottky Contact BC
      case SemiconductorRegion:
        {
          // semiconductor region should be the first region
          genius_assert(i==0);

          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

          // mapping this node to material library
          semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

          //the indepedent variable number, we only need 3 here.
          adtl::AutoDScalar::numdir = 3;
          //synchronize with material database
          semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

          AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];    V.setADValue(0, 1.0);  // psi of this node
          AutoDScalar n = x[fvm_nodes[i]->local_offset()+1];    n.setADValue(1, 1.0);  // electron density
          AutoDScalar p = x[fvm_nodes[i]->local_offset()+2];    p.setADValue(2, 1.0);  // hole density

          //Schotty Barrier Lowerring
          PetscScalar deltaVB = semi_region->material()->band->SchottyBarrierLowerring(node_data->eps(), node_data->E().size());

          //Schottky current
          PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
          AutoDScalar Fn = semi_region->material()->band->SchottyJsn(n, T, this->Work_Function() - node_data->affinity() - deltaVB) * S;
          AutoDScalar Fp = semi_region->material()->band->SchottyJsp(p, T, this->Work_Function() - node_data->affinity() + deltaVB) * S;

          // schottky boundary condition of poisson's equation
          AutoDScalar ff = V + this->Work_Function() - deltaVB - Ve;

          // the insert position
          std::vector<PetscInt> row, col;
          row.push_back(fvm_nodes[i]->global_offset()+0);
          row.push_back(fvm_nodes[i]->global_offset()+1);
          row.push_back(fvm_nodes[i]->global_offset()+2);
          col = row;

          // set Jacobian of governing equation ff
          MatSetValues(*jac, 1, &row[0], col.size(), &col[0], ff.getADValue(), ADD_VALUES);
          // process the Jacobian of Schottky current
          MatSetValues(*jac, 1, &row[1], col.size(), &col[0], Fn.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[2], col.size(), &col[0], (-Fp).getADValue(), ADD_VALUES);

          break;
        }
        // conductor region which has an interface with Schottky Contact boundary to semiconductor region
      case ConductorRegion:
        // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
        // the boundary for the corner point may be Schottky .
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
          AutoDScalar  ff = V - V_semi;

          // set Jacobian of governing equation ff
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), ff.getADValue(1), ADD_VALUES);

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

}



void SchottkyContactBC::Mix_DDM1_Electrode_Load(const Vec lx, const Mat *, double & I, PetscScalar & pdI_pdV, Vec & pdI_pdx, Vec & pdF_pdV)
{
  I = this->ext_circuit()->current();
  pdI_pdV = 0;

  VecZeroEntries(pdI_pdx);
  VecZeroEntries(pdF_pdV);

  PetscScalar * xx;
  VecGetArray(lx, &xx);

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();
  PetscScalar T = T_external();

  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // get the derivative of electrode current to schottky node
    const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, SemiconductorRegion);
    const FVM_NodeData * node_data = fvm_node->node_data();
    const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(get_fvm_node_region(*node_it, SemiconductorRegion));

    std::vector<PetscInt>    index(2);
    index[0] = fvm_node->global_offset()+1;
    index[1] = fvm_node->global_offset()+2;

    //the indepedent variable number, we need 2 here.
    adtl::AutoDScalar::numdir=2;
    //synchronize with material database
    semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

    AutoDScalar n = xx[fvm_node->local_offset()+1];   n.setADValue(0, 1.0);  // electron density of node
    AutoDScalar p = xx[fvm_node->local_offset()+2];   p.setADValue(1, 1.0);  // hole density of node

    //Schottky current
    PetscScalar S  = fvm_node->outside_boundary_surface_area();
    PetscScalar deltaVB = semi_region->material()->band->SchottyBarrierLowerring( node_data->eps(), node_data->E().size() );
    AutoDScalar Fn = semi_region->material()->band->SchottyJsn(n, T, this->Work_Function() - node_data->affinity() - deltaVB) * S;
    AutoDScalar Fp = semi_region->material()->band->SchottyJsp(p, T, this->Work_Function() - node_data->affinity() + deltaVB) * S;

    VecSetValues( pdI_pdx, index.size(), &index[0], (-Fn*current_scale).getADValue(), ADD_VALUES);
    VecSetValues( pdI_pdx, index.size(), &index[0], (-Fp*current_scale).getADValue(), ADD_VALUES);

    /*
     * displacement current
     */

    //the indepedent variable number, we need 2 here.
    adtl::AutoDScalar::numdir=2;

    FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
    for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
    {
      const FVM_Node *nb_node = (*nb_it).second;
      const FVM_NodeData * nb_node_data = nb_node->node_data();

      // potential of node
      AutoDScalar V = xx[fvm_node->local_offset()+0];   V.setADValue(0, 1.0);
      // the psi of neighbor node
      AutoDScalar V_nb = xx[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);

      // distance from nb node to this node
      PetscScalar distance = (*(fvm_node->root_node()) - *(nb_node->root_node())).size();
      // area of out surface of control volume related with neighbor node
      PetscScalar cv_boundary = fvm_node->cv_surface_area(nb_node->root_node());
      AutoDScalar dEdt;
      if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false) //second order
      {
        PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
        dEdt = ( (2-r)/(1-r)*(V-V_nb)
                 - 1.0/(r*(1-r))*(node_data->psi()-nb_node_data->psi())
                 + (1-r)/r*(node_data->psi_last()-nb_node_data->psi_last()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
      }
      else//first order
      {
        dEdt = ((V-V_nb)-(node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
      }

      AutoDScalar current_disp = cv_boundary*node_data->eps()*dEdt*current_scale;

      VecSetValue( pdI_pdx,  fvm_node->global_offset(), current_disp.getADValue(0), ADD_VALUES);
      VecSetValue( pdI_pdx,  nb_node->global_offset(),  current_disp.getADValue(1), ADD_VALUES);
    }

    {
      VecSetValue( pdF_pdV, fvm_node->global_offset(), 1.0, ADD_VALUES);
    }

  }

  VecAssemblyBegin(pdI_pdx);
  VecAssemblyBegin(pdF_pdV);

  VecAssemblyEnd(pdI_pdx);
  VecAssemblyEnd(pdF_pdV);

  VecRestoreArray(lx, &xx);
}
