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
#include "petsc_utils.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::Ohm;

/*---------------------------------------------------------------------
 * fill Schottky electrode potential into initial vector
 */
void SchottkyContactBC::DDM1_Fill_Value(Vec x, Vec L)
{
  const PetscScalar current_scale = this->z_width();

  if(Genius::processor_id() == Genius::n_processors() -1)
  {
    VecSetValue(x, this->global_offset(), this->ext_circuit()->potential(), INSERT_VALUES);

    if(this->is_inter_connect_bc())
    {
      VecSetValue(L, this->global_offset(), 1.0, INSERT_VALUES);
    }
    //for stand alone electrode
    else
    {
      const PetscScalar s = ext_circuit()->electrode_scaling(SolverSpecify::dt);
      VecSetValue(L, this->global_offset(), s, INSERT_VALUES);
    }

  }
}



///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * do pre-process to function for DDML1 solver
 */
void SchottkyContactBC::DDM1_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
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
 * build function and its jacobian for DDML1 solver
 */
void SchottkyContactBC::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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
  std::vector<PetscInt> ix;
  std::vector<PetscScalar> y;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();
  std::vector<double> current_buffer;

  const PetscScalar Work_Function = this->scalar("workfunction");

  // the electrode potential in current iteration
  genius_assert( local_offset()!=invalid_uint );
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
      const unsigned int global_offset = fvm_nodes[i]->global_offset();
      const unsigned int local_offset = fvm_nodes[i]->local_offset();

      switch ( regions[i]->type() )
      {
          // Semiconductor Region of course owns Schottky Contact BC
          case SemiconductorRegion:
          {
            // semiconductor region should be the first region
            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();
            const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);
            // mapping this node to material library
            semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

            PetscScalar V = x[local_offset+0];  // psi of this node
            PetscScalar n = x[local_offset+1];  // electron density
            PetscScalar p = x[local_offset+2];  // hole density
            PetscScalar T = T_external();

            //Schotty Barrier Lowerring
            PetscScalar deltaVB = semi_region->material()->band->SchottyBarrierLowerring(node_data->eps(), node_data->E().size());
            //Schottky current
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            PetscScalar Fn = semi_region->material()->band->SchottyJsn(n, T, Work_Function - node_data->affinity() - deltaVB) * S;
            PetscScalar Fp = semi_region->material()->band->SchottyJsp(p, T, Work_Function - node_data->affinity() + deltaVB) * S;

            PetscScalar f_psi = V + Work_Function - deltaVB - Ve;

            // set governing equation to boundary condition of poisson's equation
            ix.push_back(global_offset+0);
            y.push_back(f_psi);
            // add Schottky current to continuity equations
            // we buffer this operator
            ix.push_back(global_offset+1);
            y.push_back( Fn);
            ix.push_back(global_offset+2);
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

            // schottky thermal emit current
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
                PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
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
          // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
          // the boundary for the corner point may be Schottky. (not a nice behavier)
          case InsulatorRegion:
          {
            // psi of this node
            PetscScalar V = x[local_offset];

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            //genius_assert( regions[0]->type()==SemiconductorRegion );
            PetscScalar V_semi = x[fvm_nodes[0]->local_offset()];

            // the psi of this node is equal to corresponding psi of semiconductor node
            PetscScalar f_psi = V - V_semi;
            // set governing equation to function vector
            y.push_back( f_psi );
            ix.push_back(global_offset);

            break;
          }

          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }

  }

  if(ix.size()) VecSetValues(f, ix.size(), &ix[0], &y[0], ADD_VALUES) ;


  // the extra equation of schottky boundary
  // For voltage driven
  //
  //          _____                Ve
  //    -----|_____|----/\/\/\/\-------> to schottky electrode (Ve, I)
  //    | +     R          L       |
  //   Vapp                     C ===
  //    | -                        |
  //    |__________________________|

  //           GND
  //
  // And for current driven
  //                               Ve
  //    -->-----------------------------> to schottky electrode (Ve, I)
  //    |                          |
  //   Iapp                     C ===
  //    |__________________________|
  //           GND
  //
  // Or for inter connect
  //
  //          _____                Ve
  //    -----|_____|-------------------> to schottky electrode (Ve, I)
  //    |       R
  //    |
  // V_inter_connect
  //


  // for get the current, we must sum all the terms in current_buffer
  // NOTE: only statistic current flow belongs to on processor node
  PetscScalar current = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  ext_circuit()->potential() = Ve;
  ext_circuit()->current() = current;

  PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);

  //for inter connect electrode
  if(this->is_inter_connect_bc())
  {
    PetscScalar R = ext_circuit()->inter_connect_resistance();                               // resistance
    PetscScalar f_ext = R*current;
    VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
  }
  // for stand alone electrode
  else
  {
    PetscScalar f_ext = mna_scaling*current;
    VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
  }

  if(Genius::is_last_processor())
  {
    //for inter connect electrode
    if(this->is_inter_connect_bc())
    {
      PetscScalar V_ic = x[this->inter_connect_hub()->local_offset()];  // potential at inter connect node
      PetscScalar f_ext = Ve - V_ic;
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
    // for stand alone electrode
    else
    {
      PetscScalar f_ext = ext_circuit()->mna_function(SolverSpecify::dt);
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}







/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void SchottkyContactBC::DDM1_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
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
    for(; rnode_it!=end_rnode_it; ++rnode_it  )
    {
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;

      PetscInt row = fvm_node->global_offset();
      clear_row.push_back(row);
    }
  }

}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void SchottkyContactBC::DDM1_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // the Jacobian of Schottky boundary condition is processed here
  // we use AD again. no matter it is overkill here.

  PetscInt bc_global_offset = this->global_offset();


  const PetscScalar Work_Function = this->scalar("workfunction");

  // do schottky boundary process here

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();


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

            PetscScalar T = T_external();

            //Schotty Barrier Lowerring
            PetscScalar deltaVB = semi_region->material()->band->SchottyBarrierLowerring(node_data->eps(), node_data->E().size());

            //Schottky current
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            AutoDScalar Fn = semi_region->material()->band->SchottyJsn(n, T, Work_Function - node_data->affinity() - deltaVB) * S;
            AutoDScalar Fp = semi_region->material()->band->SchottyJsp(p, T, Work_Function - node_data->affinity() + deltaVB) * S;

            // schottky boundary condition of poisson's equation
            AutoDScalar f_psi = V + Work_Function - deltaVB - Ve;

            // the insert position
            PetscInt row[3], col[4];
            col[0] = row[0] = fvm_nodes[i]->global_offset()+0;
            col[1] = row[1] = fvm_nodes[i]->global_offset()+1;
            col[2] = row[2] = fvm_nodes[i]->global_offset()+2;
            col[3] = this->global_offset(); // the position of electrode equation

            // set Jacobian of governing equation ff
            jac->add_row(  row[0],  4,  &col[0],  f_psi.getADValue() );

            // process the Jacobian of Schottky current
            jac->add_row(  row[1],  4,  &col[0],  Fn.getADValue() );
            jac->add_row(  row[2],  4,  &col[0],  (-Fp).getADValue() );


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

            //for inter connect electrode
            if(this->is_inter_connect_bc())
            {
              PetscScalar R = ext_circuit()->inter_connect_resistance();
              jac->add_row(  bc_global_offset,  4,  &(col[0]),  (R*current_emit*current_scale).getADValue() );
            }
                // for stand alone electrode
            else
            {
              PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);
              jac->add_row(  bc_global_offset,  4,  &(col[0]),  (mna_scaling*current_emit*current_scale).getADValue() );
            }

            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                FVM_Node *nb_node = (*nb_it).first;
                FVM_NodeData * nb_node_data = nb_node->node_data();
                // the psi of neighbor node
                AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);
                // distance from nb node to this node
                PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
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

                //for inter connect electrode
                if(this->is_inter_connect_bc())
                {
                  PetscScalar R = ext_circuit()->inter_connect_resistance();
                  current_disp = R*current_disp;
                }
                // for stand alone electrode
                else
                {
                  PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);
                  current_disp = mna_scaling*current_disp;
                }

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

  // the extra equation of schottky boundary
  // For voltage driven
  //
  //          _____                Ve
  //    -----|_____|----/\/\/\/\-------> to schottky electrode (Ve, I)
  //    |       R          L       |
  //   Vapp                     C ===
  //    |__________________________|
  //           GND
  //
  // And for current driven
  //                               Ve
  //    --------------------------------> to schottky electrode (Ve, I)
  //    |                          |
  //   Iapp                     C ===
  //    |__________________________|
  //           GND
  //
  // Or for inter connect
  //
  //          _____                Ve
  //    -----|_____|-------------------> to schottky electrode (Ve, I)
  //    |       R
  //    |
  // V_inter_connect
  //
  //


  if(Genius::is_last_processor())
  {
    //for inter connect electrode
    if(this->is_inter_connect_bc())
    {
      // the external electrode equation is:
      // f_ext = Ve - V_ic + R*current;

      // d(f_ext)/d(Ve)
      jac->add( bc_global_offset,  bc_global_offset,  1.0 );
      // d(f_ext)/d(V_ic)
      jac->add( bc_global_offset,  this->inter_connect_hub()->global_offset(),  -1.0 );
    }
    //for stand alone electrode
    else
    {
      ext_circuit()->potential() = x[this->local_offset()];
      jac->add( bc_global_offset,  bc_global_offset,  ext_circuit()->mna_jacobian(SolverSpecify::dt) );
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


void SchottkyContactBC::DDM1_Electrode_Trace(Vec lx, SparseMatrix<PetscScalar> *jac, Vec pdI_pdx, Vec pdF_pdV)
{
  VecZeroEntries(pdI_pdx);
  VecZeroEntries(pdF_pdV);

  PetscScalar * xx;
  VecGetArray(lx, &xx);

  const PetscScalar Work_Function = this->scalar("workfunction");

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
    AutoDScalar Fn = semi_region->material()->band->SchottyJsn(n, T, Work_Function - node_data->affinity() - deltaVB) * S;
    AutoDScalar Fp = semi_region->material()->band->SchottyJsp(p, T, Work_Function - node_data->affinity() + deltaVB) * S;

    VecSetValues( pdI_pdx, index.size(), &index[0], (-Fn*current_scale).getADValue(), ADD_VALUES);
    VecSetValues( pdI_pdx, index.size(), &index[0], (-Fp*current_scale).getADValue(), ADD_VALUES);

    {
      VecSetValue( pdF_pdV, fvm_node->global_offset(), 1.0, ADD_VALUES);
    }

  }

  VecAssemblyBegin(pdI_pdx);
  VecAssemblyBegin(pdF_pdV);

  VecAssemblyEnd(pdI_pdx);
  VecAssemblyEnd(pdF_pdV);

  VecRestoreArray(lx, &xx);

  //delete electrode current equation, omit the effect of external resistance
  PetscInt bc_global_offset = this->global_offset();
  jac->clear_row(bc_global_offset, 1.0);
}


/*---------------------------------------------------------------------
 * update electrode IV
 */
void SchottkyContactBC::DDM1_Update_Solution(PetscScalar *)
{
  Parallel::sum(ext_circuit()->current());
  this->ext_circuit()->update();
}

