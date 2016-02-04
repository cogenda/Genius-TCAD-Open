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
#include "boundary_condition_schottky.h"
#include "parallel.h"

using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::kb;
using PhysicalUnit::e;



///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for Mixed DDML1 solver
 */
void SchottkyContactBC::Mix_DDM1_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // only process nodes belong to this processor
    if( (*node_it)->processor_id() != Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with this boundary condition
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;

      PetscInt row = fvm_node->global_offset();
      clear_row.push_back(row);
    }
  }

}

/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML1 solver
 */
void SchottkyContactBC::Mix_DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Schottky boundary condition is processed here

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // data buffer for ADD_VALUES
  std::vector<int> ix;
  std::vector<PetscScalar> y;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;
  std::vector<double> current_buffer;

  const PetscScalar Work_Function = this->scalar("workfunction");

  // the electrode potential in current iteration
  PetscScalar Ve = x[this->local_offset()];

  // search and process all the nodes belongs to this bc
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
            PetscScalar Fn = semi_region->material()->band->SchottyJsn(n, T, Work_Function - node_data->affinity() - deltaVB) * S;
            PetscScalar Fp = semi_region->material()->band->SchottyJsp(p, T, Work_Function - node_data->affinity() + deltaVB) * S;

            PetscScalar ff = V + Work_Function - deltaVB - Ve;

            // set governing equation to boundary condition of poisson's equation
            VecSetValue(f, fvm_nodes[i]->global_offset()+0, ff, ADD_VALUES);
            // add Schottky current to continuity equations
            // we buffer this operator
            ix.push_back(fvm_nodes[i]->global_offset()+1);
            y.push_back( Fn);
            ix.push_back(fvm_nodes[i]->global_offset()+2);
            y.push_back(-Fp);

            // compute the current flow out of schottky electrode

            PetscScalar inject_current = -(Fn+Fp);
            // compute Ic
            if(SolverSpecify::TimeDependent == true)
            {
              //second order
              if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false)
              {
                PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
                PetscScalar Tn = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                                 / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_nodes[i]->volume();
                PetscScalar Tp = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                                 / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_nodes[i]->volume();
                inject_current += -(Tn+Tp);
              }
              else //first order
              {
                PetscScalar Tn = -(n - node_data->n())/SolverSpecify::dt*fvm_nodes[i]->volume();
                PetscScalar Tp = -(p - node_data->p())/SolverSpecify::dt*fvm_nodes[i]->volume();
                inject_current += -(Tn+Tp);
              }
            }
            current_buffer.push_back(inject_current);

            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                const FVM_Node *nb_node = (*nb_it).first;
                const FVM_NodeData * nb_node_data = nb_node->node_data();
                // the psi of neighbor node
                PetscScalar V_nb = x[nb_node->local_offset()+0];
                // distance from nb node to this node
                PetscScalar distance = fvm_nodes[i]->distance(nb_node);
                // area of out surface of control volume related with neighbor node
                PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node);
                PetscScalar dEdt;
                if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false) //second order
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

                current_buffer.push_back( cv_boundary*node_data->eps()*dEdt );
              }
            }
            break;
          }
          // conductor region which has an interface with Schottky Contact boundary to semiconductor region
          case ElectrodeRegion:
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
            VecSetValue(f, fvm_nodes[i]->global_offset(), ff, ADD_VALUES);

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
            VecSetValue(f, fvm_nodes[i]->global_offset(), ff, ADD_VALUES);

            break;
          }
          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }

  }

  if(ix.size()) VecSetValues(f, ix.size(), &ix[0], &y[0], ADD_VALUES) ;

  // for get the current, we must sum all the terms in current_buffer
  PetscScalar current = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );
  this->current() = current*A;

  // Add current to spice node
  VecSetValue(f, this->global_offset(), current, ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}






/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for Mixed DDML1 solver
 */
void SchottkyContactBC::Mix_DDM1_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // only process nodes belong to this processor
    if( (*node_it)->processor_id() != Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with this boundary condition
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;

      PetscInt row = fvm_node->global_offset();
      clear_row.push_back(row);
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML1 solver
 */
void SchottkyContactBC::Mix_DDM1_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{

  // the Jacobian of Schottky boundary condition is processed here


  // we use AD again. no matter it is overkill here.

  PetscInt bc_global_offset = this->global_offset();

  const PetscScalar T = T_external();

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;

  const PetscScalar Work_Function = this->scalar("workfunction");

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
            const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

            // mapping this node to material library
            semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

            //the indepedent variable number, we only need 4 here.
            adtl::AutoDScalar::numdir=4;
            //synchronize with material database
            semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

            AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];    V.setADValue(0, 1.0);  // psi of this node
            AutoDScalar n = x[fvm_nodes[i]->local_offset()+1];    n.setADValue(1, 1.0);  // electron density
            AutoDScalar p = x[fvm_nodes[i]->local_offset()+2];    p.setADValue(2, 1.0);  // hole density

            // the electrode potential in current iteration
            genius_assert( local_offset()!=invalid_uint );
            AutoDScalar Ve = x[this->local_offset()];             Ve.setADValue(3, 1.0);

            //Schotty Barrier Lowerring
            PetscScalar deltaVB = semi_region->material()->band->SchottyBarrierLowerring(node_data->eps(), node_data->E().size());

            //Schottky current
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            AutoDScalar Fn = semi_region->material()->band->SchottyJsn(n, T, Work_Function - node_data->affinity() - deltaVB) * S;
            AutoDScalar Fp = semi_region->material()->band->SchottyJsp(p, T, Work_Function - node_data->affinity() + deltaVB) * S;

            // schottky boundary condition of poisson's equation
            AutoDScalar ff = V + Work_Function - deltaVB - Ve;

            // the insert position
            std::vector<PetscInt> row, col;
            row.push_back(fvm_nodes[i]->global_offset()+0);
            row.push_back(fvm_nodes[i]->global_offset()+1);
            row.push_back(fvm_nodes[i]->global_offset()+2);
            col = row;
            col.push_back(this->global_offset()); // the position of electrode equation

            // set Jacobian of governing equation ff
            jac->add_row(  row[0],  col.size(),  &col[0],  ff.getADValue() );

            // process the Jacobian of Schottky current
            jac->add_row(  row[1],  col.size(),  &col[0],  Fn.getADValue() );
            jac->add_row(  row[2],  col.size(),  &col[0],  (-Fp).getADValue() );


            // process the Jacobian of current flow out of schottky electrode

            // compute the schottky thermal emit current
            AutoDScalar current_emit = -(Fn+Fp);

            if(SolverSpecify::TimeDependent == true)
            {
              //second order
              if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false)
              {
                PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
                AutoDScalar Tn = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                                 / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_nodes[i]->volume();
                AutoDScalar Tp = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                                 / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_nodes[i]->volume();
                current_emit += -(Tn+Tp);
              }
              else //first order
              {
                AutoDScalar Tn = -(n - node_data->n())/SolverSpecify::dt*fvm_nodes[i]->volume();
                AutoDScalar Tp = -(p - node_data->p())/SolverSpecify::dt*fvm_nodes[i]->volume();
                current_emit += -(Tn+Tp);
              }
            }

            jac->add_row(  bc_global_offset,  col.size(),  &(col[0]),  (current_emit*current_scale).getADValue() );

            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                FVM_Node *nb_node = (*nb_it).first;
                FVM_NodeData * nb_node_data = nb_node->node_data();
                // the psi of neighbor node
                AutoDScalar V_nb = x[nb_node->local_offset()+0];  V_nb.setADValue(1, 1.0);
                // distance from nb node to this node
                PetscScalar distance = fvm_nodes[i]->distance(nb_node);
                // area of out surface of control volume related with neighbor node
                PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node);

                AutoDScalar dEdt;
                if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false) //second order
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

                jac->add( bc_global_offset,  fvm_nodes[i]->global_offset()+0,  current_disp.getADValue(0) );
                jac->add( bc_global_offset,  nb_node->global_offset()+0,  current_disp.getADValue(1) );
              }
            }

            break;
          }
          // conductor region which has an interface with Schottky Contact boundary to semiconductor region
          case ElectrodeRegion:
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
            jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[i]->global_offset(),  ff.getADValue(0) );
            jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[0]->global_offset(),  ff.getADValue(1) );

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
 * update electrode IV
 */
void SchottkyContactBC::Mix_DDM1_Update_Solution(PetscScalar *)
{
  Parallel::sum(this->current());
}


