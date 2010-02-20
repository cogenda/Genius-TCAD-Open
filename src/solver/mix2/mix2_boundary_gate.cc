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
#include "petsc_utils.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;





///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////



/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML2 solver
 */
void GateContactBC::Mix_DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Gate boundary condition is processed here

  // note, we will use INSERT_VALUES to set values of vec f
  // if the previous operator is not insert_VALUES, we should assembly the vec
  if( (add_value_flag != INSERT_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // data buffer for add row to row location
  std::vector<PetscInt> src_row;
  std::vector<PetscInt> dst_row;

  // data buffer for INSERT_VALUES
  std::vector<int> iy;
  std::vector<PetscScalar> y;


  // the electrode current, since the electrode may be partitioned into several processor,
  // we should collect it.
  std::vector<double> current_buffer;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();

  // the electrode potential in current iteration
  PetscScalar Ve = ext_circuit()->Vapp();

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

          PetscScalar V = x[fvm_nodes[i]->local_offset()+0];  // psi of this node
          PetscScalar T = x[fvm_nodes[i]->local_offset()+1];  // lattice temperature

          // the governing equation
          PetscScalar ff = V + Work_Function() - Ve;

          // set governing equation to function vector
          iy.push_back(fvm_nodes[i]->global_offset());
          y.push_back( ff );


          // add heat flux out of gate boundary to lattice temperature equatiuon
          // when this bc is external boundary
          if(node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion))
          {
            PetscScalar h = this->Heat_Transfer();
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            PetscScalar fT = h*(T_external()-T)*S;
            VecSetValue(f, fvm_nodes[i]->global_offset()+1, fT, ADD_VALUES);
          }


          // MOS gate can have displacement current and tunneling current

          // displacement current
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

            current_buffer.push_back( cv_boundary*node_data->eps()*dEdt*current_scale );
          }

          // FIXME: tunneling current.

          break;
        }


        // conductor region
      case ConductorRegion:
        {

          PetscScalar V = x[fvm_nodes[i]->local_offset()+0];   // psi of this node
          PetscScalar T = x[fvm_nodes[i]->local_offset()+1];   // lattice temperature of this node

          // since the region is sorted, we know region[0] is Insulator region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding Insulator region
          PetscScalar V_in = x[fvm_nodes[0]->local_offset()+0];
          PetscScalar T_in = x[fvm_nodes[0]->local_offset()+1];

          // let insulator node process the complete governing equation of heat transfer at interface
          // then we can set the node temperature of conductor region equal to node temperature at insulator region
          src_row.push_back(fvm_nodes[i]->global_offset()+1);
          dst_row.push_back(fvm_nodes[0]->global_offset()+1);

          // the psi of this node is equal to corresponding psi of Insulator node
          PetscScalar ff1 = V - V_in;
          // the T of this node is equal to corresponding T of Insulator node
          PetscScalar ff2 = T - T_in;

          // set governing equation to function vector
          iy.push_back(fvm_nodes[i]->global_offset()+0);
          y.push_back( ff1 );

          iy.push_back(fvm_nodes[i]->global_offset()+1);
          y.push_back( ff2 );

          break;
        }
      case VacuumRegion:
        break;
      default: genius_error(); //we should never reach here
      }
    }

  }

  // add src row to dst row, it will assemble vec automatically
  PetscUtils::VecAddRowToRow(f, src_row, dst_row);

  if(iy.size()) VecSetValues(f, iy.size(), &iy[0], &y[0], INSERT_VALUES) ;


  // we first gather the electrode current
  Parallel::allgather(current_buffer);

  // for get the current, we must sum all the terms in current_buffer
  PetscScalar current = std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  // save the IV of current iteration
  ext_circuit()->current() = current;
  ext_circuit()->potential() = Ve;

  // the last operator is INSERT_VALUES
  add_value_flag = INSERT_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for Mixed DDML2 solver
 */
void GateContactBC::Mix_DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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
          FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
          FVM_Node::fvm_ghost_node_iterator gn_it_end = fvm_nodes[i]->ghost_node_end();
          for(; gn_it != gn_it_end; ++gn_it)
          {
            const FVM_Node * ghost_fvm_node = (*gn_it).first;
            // skip NULL neighbor which means the node is on Neumann boundary
            if(ghost_fvm_node==NULL) continue;

            MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, ghost_fvm_node->global_offset()+1, 0, ADD_VALUES);
          }
          break;
        }
      case ConductorRegion:
        {
          // insert none zero pattern
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, fvm_nodes[0]->global_offset()+0, 0, ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[0]->global_offset()+1, 0, ADD_VALUES);
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
 * build function and its jacobian for Mixed DDML2 solver
 */
void GateContactBC::Mix_DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of GateContact boundary condition is processed here
  // we use AD again. no matter it is overkill here.



  // first, we zero all the rows corresponding to GateContact bc.
  {
    // data buffer for add row to rwo location
    std::vector<PetscInt> src_row;
    std::vector<PetscInt> dst_row;

    // the indicator which rows should be set to zero
    std::vector<PetscInt> row_for_clear;


    BoundaryCondition::const_node_iterator node_it = nodes_begin();
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for(; node_it!=end_it; ++node_it )
    {

      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      std::vector<const FVM_Node *> fvm_nodes;

      // search all the fvm_node which has *node_it as root node
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
      for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
      {
        const SimulationRegion * region = (*rnode_it).second.first;
        fvm_nodes.push_back((*rnode_it).second.second);

        switch ( region->type() )
        {
        case InsulatorRegion:
          row_for_clear.push_back(fvm_nodes[i]->global_offset() + 0);
          break;
        case ConductorRegion:
          // if code reaches here, then the gate bc is an interface to conductor region
          src_row.push_back(fvm_nodes[i]->global_offset()+1);
          dst_row.push_back(fvm_nodes[0]->global_offset()+1);

          row_for_clear.push_back(fvm_nodes[i]->global_offset() + 0);
          row_for_clear.push_back(fvm_nodes[i]->global_offset() + 1);
          break;
        case VacuumRegion:
          break;
        default: genius_error();
        }
      }
    }

    //ok, we add source rows to destination rows
    PetscUtils::MatAddRowToRow(*jac, src_row, dst_row);

    MatZeroRows(*jac, row_for_clear.size(), row_for_clear.empty() ? NULL : &row_for_clear[0], 0.0);

  }


  // the electrode potential in current iteration
  PetscScalar Ve = ext_circuit()->Vapp();

  // after that, we should do gate boundary process here

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

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;


          AutoDScalar  V = x[fvm_nodes[i]->local_offset()+0]; V.setADValue(0, 1.0); // psi of this node
          AutoDScalar  T = x[fvm_nodes[i]->local_offset()+1]; T.setADValue(1, 1.0); // psi of this node

          // the governing equation of potential
          AutoDScalar ff = V + Work_Function() - Ve;

          // the insert position
          std::vector<PetscInt> row, col;
          row.push_back(fvm_nodes[i]->global_offset()+0);
          row.push_back(fvm_nodes[i]->global_offset()+1);
          col = row;

          // process the Jacobian of governing equation of potential
          MatSetValues(*jac, 1, &row[0], col.size(), &col[0], ff.getADValue(), ADD_VALUES);

          // process the Jacobian of equation of T
          // if this gate bc is external boundary, set heat flux here
          if(node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion))
          {
            PetscScalar h = this->Heat_Transfer();
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            AutoDScalar fT = h*(T_external()-T)*S;
            MatSetValues(*jac, 1, &row[1], col.size(), &col[0], fT.getADValue(), ADD_VALUES);
          }

          break;
        }
        // conductor region (gate) which has an interface with insulator region
      case ConductorRegion:
        {

          //the indepedent variable number, we need 4 here.
          adtl::AutoDScalar::numdir=4;

          AutoDScalar  V = x[fvm_nodes[i]->local_offset()+0]; V.setADValue(0,1.0); // psi of this node
          AutoDScalar  T = x[fvm_nodes[i]->local_offset()+1]; T.setADValue(1,1.0); // lattice temperature of this node

          // since the region is sorted, we know region[0] is insulator region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding insulator region
          AutoDScalar  V_in = x[fvm_nodes[0]->local_offset()+0]; V_in.setADValue(2,1.0);
          AutoDScalar  T_in = x[fvm_nodes[0]->local_offset()+1]; T_in.setADValue(3,1.0);

          // the psi of this node is equal to corresponding psi of insulator node
          AutoDScalar  ff1 = V - V_in;
          // the T of this node is equal to corresponding T of insulator node
          AutoDScalar  ff2 = T - T_in;

          // set Jacobian of governing equation ff1 and ff2
          std::vector<PetscInt> rows,cols;
          rows.push_back(fvm_nodes[i]->global_offset()+0);
          rows.push_back(fvm_nodes[i]->global_offset()+1);
          cols = rows;
          cols.push_back(fvm_nodes[0]->global_offset()+0);
          cols.push_back(fvm_nodes[0]->global_offset()+1);

          MatSetValues(*jac, 1, &rows[0], cols.size(), &cols[0], ff1.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &rows[1], cols.size(), &cols[0], ff2.getADValue(), ADD_VALUES);

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




void GateContactBC::Mix_DDM2_Electrode_Load(const Vec lx, const Mat *, double & I, PetscScalar & pdI_pdV, Vec & pdI_pdx, Vec & pdF_pdV)
{
  I = this->ext_circuit()->current();
  pdI_pdV = 0;

  VecZeroEntries(pdI_pdx);
  VecZeroEntries(pdF_pdV);

  PetscScalar * xx;
  VecGetArray(lx, &xx);

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();


  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // get the derivative of electrode current to gate node
    const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, InsulatorRegion);
    const FVM_NodeData * node_data = fvm_node->node_data();

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
