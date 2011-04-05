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

//  $Id: ddm1_boundary_gate.cc,v 1.5 2008/07/09 05:58:16 gdiso Exp $

// C++ includes
#include <numeric>

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_gate.h"
#include "parallel.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * fill gate electrode potential into initial vector
 */
void GateContactBC::DDM1_Fill_Value(Vec x, Vec L)
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
      if(ext_circuit()->is_voltage_driven())
      {
        const PetscScalar R = this->ext_circuit()->R();
        VecSetValue(L, this->global_offset(), 1.0/((1.0+R)*current_scale), INSERT_VALUES);
      }

      if(ext_circuit()->is_current_driven())
      {
        VecSetValue(L, this->global_offset(), 1.0/(current_scale), INSERT_VALUES);
      }
    }

  }
}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////



/*---------------------------------------------------------------------
 * do pre-process to function for DDML1 solver
 */
void GateContactBC::DDM1_Function_Preprocess(Vec f, std::vector<PetscInt> &src_row,
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
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;

      PetscInt row = fvm_node->global_offset();
      clear_row.push_back(row);
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void GateContactBC::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Gate boundary condition is processed here

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // the electrode current, since the electrode may be partitioned into several processor,
  // we should collect it.
  std::vector<double> current_buffer;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();

  // the electrode potential in current iteration
  genius_assert( local_offset()!=invalid_uint );
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
          // insulator region.
          case InsulatorRegion:
          {

            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

            // psi of this node
            PetscScalar V = x[fvm_nodes[i]->local_offset()];

            // the governing equation
            PetscScalar ff = V + Work_Function() - Ve;

            // set governing equation to function vector
            VecSetValue(f, fvm_nodes[i]->global_offset(), ff, ADD_VALUES);

            // MOS gate can have displacement current and tunneling current

            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                const FVM_Node *nb_node = (*nb_it).second;
                const FVM_NodeData * nb_node_data = nb_node->node_data();
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

                current_buffer.push_back( cv_boundary*node_data->eps()*dEdt );
              }
            }
            // FIXME: tunneling current.

            break;
          }


          // conductor region
          case ElectrodeRegion:
          {

            // psi of this node
            PetscScalar V = x[fvm_nodes[i]->local_offset()];

            // since the region is sorted, we know region[0] is Insulator region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding Insulator region

            PetscScalar V_in = x[fvm_nodes[0]->local_offset()];

            // the psi of this node is equal to corresponding psi of Insulator node
            PetscScalar ff = V - V_in;

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

  // the extra equation of gate boundary
  // For voltage driven
  //
  //          _____                Ve
  //    -----|_____|----/\/\/\/\-------> to gate electrode (Ve, I)
  //    | +     R          L       |
  //   Vapp                     C ===
  //    | -                        |
  //    |__________________________|

  //           GND
  //
  // And for current driven
  // NOTE: It is dangerous to attach current source to MOS gate!
  //                               Ve
  //    -->-----------------------------> to gate electrode (Ve, I)
  //    |                          |
  //   Iapp                     C ===
  //    |__________________________|
  //           GND

  // Or for inter connect
  //
  //          _____                Ve
  //    -----|_____|-------------------> to gate electrode (Ve, I)
  //    |       R
  //    |
  // V_inter_connect
  //
  //

  // for get the current, we must sum all the terms in current_buffer
  // NOTE: only statistic current flow belongs to on processor node
  PetscScalar current = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  //for inter connect electrode
  if(this->is_inter_connect_bc())
  {
    PetscScalar R = ext_circuit()->R();                               // resistance
    PetscScalar f_ext = R*current;
    VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
  }
  // for stand alone electrode
  else
  {
    if(ext_circuit()->is_voltage_driven())
    {
      PetscScalar R = ext_circuit()->R();             // resistance
      PetscScalar L = ext_circuit()->L();             // inductance
      PetscScalar dt = SolverSpecify::dt;
      PetscScalar f_ext = (L/dt+R)*current;
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }

    if(ext_circuit()->is_current_driven())
    {
      PetscScalar f_ext = current;
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
  }

  if(Genius::is_last_processor())
  {
    //for inter connect electrode
    if(this->is_inter_connect_bc())
    {
      PetscScalar V_ic = x[this->inter_connect_hub()->local_offset()];  // potential at inter connect node
      PetscScalar R = ext_circuit()->R();                               // resistance
      PetscScalar f_ext = Ve - V_ic;
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
    // for stand alone electrode
    else
    {
      if(ext_circuit()->is_voltage_driven())
      {
        PetscScalar Vapp = ext_circuit()->Vapp();       // application voltage
        PetscScalar R = ext_circuit()->R();             // resistance
        PetscScalar C = ext_circuit()->C();             // capacitance
        PetscScalar L = ext_circuit()->L();             // inductance
        PetscScalar I = ext_circuit()->current();       // the previous step current flow into electrode
        PetscScalar Ic = ext_circuit()->cap_current();  // the previous step current flow pass though cap to ground.
        PetscScalar P  = ext_circuit()->potential();    // the previous step potential of the electrode
        PetscScalar dt = SolverSpecify::dt;
        PetscScalar f_ext = (Ve-Vapp) + (L/dt+R)*C/dt*Ve - (L/dt+R)*C/dt*P - L/dt*(I+Ic);
        VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
      }

      if(ext_circuit()->is_current_driven())
      {
        PetscScalar Iapp = ext_circuit()->Iapp();         // application current
        PetscScalar Ic   = ext_circuit()->cap_current();  // the previous step current flow pass though cap to ground.
        PetscScalar f_ext = Ic - Iapp;
        VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
      }
    }
  }


  // save the IV of current iteration
  ext_circuit()->current_itering() = current;
  ext_circuit()->potential_itering() = Ve;

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML1 solver
 */
void GateContactBC::DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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
          case InsulatorRegion:
          {
            // reserve for electrode potential item
            MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, this->global_offset(), 0, ADD_VALUES);

            break;
          }
          case ElectrodeRegion:
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

  // reserve jacobian entries for the circuit equation of gate electrode
  {

    std::vector<PetscInt> bc_node_reserve;
    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // get the derivative of electrode current to ohmic node
      const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, InsulatorRegion);
      if(fvm_node->on_processor())
      {
        bc_node_reserve.push_back(fvm_node->global_offset()+0);

        // get the derivative of electrode current to neighbors of bc node
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
        for(; nb_it != nb_it_end; ++nb_it)
        {
          const FVM_Node *  fvm_nb_node = (*nb_it).second;
          bc_node_reserve.push_back(fvm_nb_node->global_offset()+0);
        }
      }
    }
    Parallel::allgather(bc_node_reserve);

    if(Genius::processor_id() == Genius::n_processors()-1)
    {
      PetscInt bc_global_offset = this->global_offset();

      MatSetValue(*jac, bc_global_offset, bc_global_offset, 0, ADD_VALUES);

      if(this->is_inter_connect_bc())
        MatSetValue(*jac, bc_global_offset, this->inter_connect_hub()->global_offset(), 0, ADD_VALUES);

      if(bc_node_reserve.size())
      {
        std::vector<PetscScalar> bc_node_reserve_zero(bc_node_reserve.size(), 0.0);
        MatSetValues(*jac, 1, &bc_global_offset, bc_node_reserve.size(), &bc_node_reserve[0], &bc_node_reserve_zero[0], ADD_VALUES);
      }
    }

  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void GateContactBC::DDM1_Jacobian_Preprocess(Mat *jac, std::vector<PetscInt> &src_row,
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
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;

      PetscInt row = fvm_node->global_offset();
      clear_row.push_back(row);
    }
  }
}



/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void GateContactBC::DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of GateContact boundary condition is processed here
  // we use AD again. no matter it is overkill here.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  PetscInt bc_global_offset = this->global_offset();
  PetscScalar R = ext_circuit()->R();             // resistance
  PetscScalar C = ext_circuit()->C();             // capacitance
  PetscScalar L = ext_circuit()->L();             // inductance
  PetscScalar dt = SolverSpecify::dt;


  // do gate boundary process here

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
          case InsulatorRegion:
          {

            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0, 1.0);

            // the electrode potential in current iteration
            genius_assert( local_offset()!=invalid_uint );
            AutoDScalar Ve = x[this->local_offset()];         Ve.setADValue(1, 1.0);

            // the governing equation of gate electrode
            AutoDScalar ff = V + Work_Function() - Ve;

            // the insert position
            std::vector<PetscInt> row, col;
            row.push_back(fvm_nodes[i]->global_offset());
            col = row;
            col.push_back(this->global_offset()); // the position of electrode equation

            // set Jacobian of governing equation ff
            MatSetValues(*jac, 1, &row[0], col.size(), &col[0], ff.getADValue(), ADD_VALUES);

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
                const FVM_Node *nb_node = (*nb_it).second;
                const FVM_NodeData * nb_node_data = nb_node->node_data();

                // the psi of this node
                AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0, 1.0);
                // the psi of neighbor node
                AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);

                // distance from nb node to this node
                PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();

                // area of out surface of control volume related with neighbor node
                PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node->root_node());
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

                // consider electrode connect

                //for inter connect electrode
                if(this->is_inter_connect_bc())
                  current_disp = R*current_disp;
                // for stand alone electrode
                else
                {
                  if(ext_circuit()->is_voltage_driven())
                    current_disp = (L/dt+R)*current_disp;
                  else
                    current_disp = current_disp;
                }

                MatSetValue(*jac, bc_global_offset, fvm_nodes[i]->global_offset(), current_disp.getADValue(0), ADD_VALUES);
                MatSetValue(*jac, bc_global_offset, nb_node->global_offset(), current_disp.getADValue(1), ADD_VALUES);
              }
            }

            //FIXME tunneling current should be considered here

            break;
          }
          // conductor region (gate) which has an interface with insulator region
          case ElectrodeRegion:
          {

            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

            // since the region is sorted, we know region[0] is Insulator region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding Insulator region
            AutoDScalar  V_in = x[fvm_nodes[0]->local_offset()]; V_in.setADValue(1,1.0);

            // the psi of this node is equal to corresponding psi of Insulator node
            AutoDScalar  ff = V - V_in;

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


  // the extra equation of gate boundary
  // For voltage driven
  //
  //          _____                Ve
  //    -----|_____|----/\/\/\/\-------> to gate electrode (Ve, I)
  //    |       R          L       |
  //   Vapp                     C ===
  //    |__________________________|
  //           GND
  //
  // And for current driven
  // NOTE: It is dangerous to attach current source to MOS gate!
  //                               Ve
  //    --------------------------------> to gate electrode (Ve, I)
  //    |                          |
  //   Iapp                     C ===
  //    |__________________________|
  //           GND
  //
  // Or for inter connect
  //
  //          _____                Ve
  //    -----|_____|-------------------> to gate electrode (Ve, I)
  //    |       R
  //    |
  // V_inter_connect
  //
  //

  if(Genius::is_last_processor())
  {
    // here we process the external circuit, we do not use AD here
    // NOTE current item such as (L/dt+R)*current and current has already been processed before

    //for inter connect electrode
    if(this->is_inter_connect_bc())
    {
      // the external electrode equation is:
      // f_ext = Ve - V_ic + R*current;

      // d(f_ext)/d(Ve)
      MatSetValue(*jac, bc_global_offset, bc_global_offset, 1.0, ADD_VALUES);
      // d(f_ext)/d(V_ic)
      MatSetValue(*jac, bc_global_offset, this->inter_connect_hub()->global_offset(), -1.0, ADD_VALUES);
    }
    //for stand alone electrode
    else
    {
      if(ext_circuit()->is_voltage_driven())
      {
        // the external electrode equation is:
        // f_ext = (L/dt+R)*current + (Ve-Vapp) + (L/dt+R)*C/dt*Ve - (L/dt+R)*C/dt*P - L/dt*(I+Ic);

        // d(f_ext)/d(Ve)
        MatSetValue(*jac, bc_global_offset, bc_global_offset, 1+(L/dt+R)*C/dt, ADD_VALUES);
      }

      if(ext_circuit()->is_current_driven())
      {
        // the external electrode equation is:
        // f_ext = current + Ic - Iapp;
        // so nothing to do
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * update electrode IV
 */
void GateContactBC::DDM1_Update_Solution(PetscScalar *)
{
  Parallel::sum(ext_circuit()->current_itering());
  this->ext_circuit()->update();
}
