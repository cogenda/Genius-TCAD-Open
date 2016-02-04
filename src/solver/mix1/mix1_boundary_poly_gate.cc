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
#include "boundary_condition_poly_gate.h"
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
void PolyGateContactBC::Mix_DDM1_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
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
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;

      switch ( region->type() )
      {
      case SemiconductorRegion:
        {
          PetscInt row = fvm_node->global_offset();
          clear_row.push_back(row+0);
          clear_row.push_back(row+1);
          clear_row.push_back(row+2);
          break;
        }
      case InsulatorRegion:
        {
          PetscInt row = fvm_node->global_offset();
          clear_row.push_back(row);
          break;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML1 solver
 */
void PolyGateContactBC::Mix_DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Gate boundary condition is processed here

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }


  const PetscScalar T = T_external();

  // the electrode current, since the electrode may be partitioned into several processor,
  // we should collect it.
  std::vector<double> current_buffer;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;

  const PetscScalar Work_Function = this->scalar("workfunction");

  // the electrode potential in current iteration
  PetscScalar Ve = x[this->local_offset()];

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
        // semiconductor region.
      case SemiconductorRegion:
        {

          // semiconductor region should be the first region
          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();
          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);


          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()];
          PetscScalar n = x[fvm_nodes[i]->local_offset()+1];
          PetscScalar p = x[fvm_nodes[i]->local_offset()+2];

          // mapping this node to material library
          semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

          PetscScalar nie = semi_region->material()->band->nie(p, n, T);

          // the governing equation
          PetscScalar f_psi = V + Work_Function - Ve;

          PetscScalar  electron_density;
          PetscScalar  hole_density;
          PetscScalar  net_dpoing = node_data->Net_doping();
          if( net_dpoing <0 )                   //p-type
          {
            hole_density = (-net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
            electron_density = nie*nie/hole_density;
          }
          else                                  //n-type
          {
            electron_density = (net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
            hole_density = nie*nie/electron_density;
          }

          // set governing equation to function vector
          VecSetValue(f, fvm_nodes[i]->global_offset()+0, f_psi, ADD_VALUES);
          VecSetValue(f, fvm_nodes[i]->global_offset()+1, n - electron_density, ADD_VALUES);
          VecSetValue(f, fvm_nodes[i]->global_offset()+2, p - hole_density, ADD_VALUES);


          // FIXME: tunneling current.

          break;
        }


        // insulator region
      case InsulatorRegion:
        {
          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()];

          PetscScalar V_0 = x[fvm_nodes[0]->local_offset()];

          // the psi of this node is equal to corresponding psi of semiconductor node
          PetscScalar ff = V - V_0;

          // set governing equation to function vector
          VecSetValue(f, fvm_nodes[i]->global_offset(), ff, ADD_VALUES);

          // MOS gate can have displacement current and tunneling current

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

      case VacuumRegion:
        break;
      default: genius_error(); //we should never reach here
      }
    }

  }


  // for get the current, we must sum all the terms in current_buffer
  PetscScalar current = std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );
  this->current() = current*A;

  // Add current to spice node
  VecSetValue(f, this->global_offset(), current, ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}






/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for  Mixed DDML1 solver
 */
void PolyGateContactBC::Mix_DDM1_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
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
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;

      switch ( region->type() )
      {
      case SemiconductorRegion:
        {
          PetscInt row = fvm_node->global_offset();
          clear_row.push_back(row+0);
          clear_row.push_back(row+1);
          clear_row.push_back(row+2);
          break;
        }
      case InsulatorRegion:
        {
          PetscInt row = fvm_node->global_offset();
          clear_row.push_back(row);
          break;
        }
      }
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML1 solver
 */
void PolyGateContactBC::Mix_DDM1_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // the Jacobian of GateContact boundary condition is processed here
  // we use AD again. no matter it is overkill here.

  const PetscScalar T = T_external();

  PetscInt bc_global_offset = this->global_offset();


  // do gate boundary process here

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
      case SemiconductorRegion:
        {
          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();
          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // psi of this node
          AutoDScalar  V = x[fvm_nodes[i]->local_offset()];   V.setADValue(0, 1.0);
          AutoDScalar  n = x[fvm_nodes[i]->local_offset()+1]; n.setADValue(0, 1.0);
          AutoDScalar  p = x[fvm_nodes[i]->local_offset()+2]; p.setADValue(0, 1.0);

          // the electrode potential in current iteration
          genius_assert( local_offset()!=invalid_uint );
          AutoDScalar Ve = x[this->local_offset()];         Ve.setADValue(1, 1.0);

          // the governing equation of gate electrode
          AutoDScalar f_psi = V + Work_Function - Ve;


          semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

          AutoDScalar nie = semi_region->material()->band->nie(p, n, T);

          AutoDScalar  electron_density;
          AutoDScalar  hole_density;
          PetscScalar  net_dpoing = node_data->Net_doping();
          if( net_dpoing <0 )   //p-type
          {
            hole_density = (-net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
            electron_density = nie*nie/hole_density;
          }
          else                               //n-type
          {
            electron_density = (net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
            hole_density = nie*nie/electron_density;
          }

          AutoDScalar fn =  n - electron_density;  //governing equation for electron density
          AutoDScalar fp =  p - hole_density;      //governing equation for hole density

          // set Jacobian of governing equation f_psi
          jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[i]->global_offset(),  f_psi.getADValue(0) );
          jac->add( fvm_nodes[i]->global_offset(),  this->global_offset(),  f_psi.getADValue(1) );

          jac->add( fvm_nodes[i]->global_offset()+1,  fvm_nodes[i]->global_offset()+1,  fn.getADValue(0) );
          jac->add( fvm_nodes[i]->global_offset()+2,  fvm_nodes[i]->global_offset()+2,  fp.getADValue(0) );

          //FIXME tunneling current should be considered here

          break;
        }
      case InsulatorRegion:
        {

          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // psi of this node
          AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

          AutoDScalar  V_0 = x[fvm_nodes[0]->local_offset()]; V_0.setADValue(1,1.0);

          // the psi of this node is equal to corresponding psi of semiconductor node
          AutoDScalar  ff = V - V_0;

          // set Jacobian of governing equation ff
          jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[i]->global_offset(),  ff.getADValue(0) );
          jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[0]->global_offset(),  ff.getADValue(1) );


          /*
           * process the Jacobian of current flow out of gate electrode
           */

          // compute displacement current

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // displacement current
          if(SolverSpecify::TimeDependent == true)
          {
            FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
            for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
            {
              const FVM_Node *nb_node = (*nb_it).first;
              const FVM_NodeData * nb_node_data = nb_node->node_data();

              // the psi of this node
              AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0, 1.0);
              // the psi of neighbor node
              AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);

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

              jac->add( bc_global_offset,  fvm_nodes[i]->global_offset(),  current_disp.getADValue(0) );
              jac->add( bc_global_offset,  nb_node->global_offset(),  current_disp.getADValue(1) );
            }
          }
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
void PolyGateContactBC::Mix_DDM1_Update_Solution(PetscScalar *)
{
  Parallel::sum(this->current());
}



