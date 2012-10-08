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

#include "asinh.hpp" // for asinh

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_ohmic.h"
#include "spice_ckt.h"
#include "parallel.h"
#include "mathfunc.h"

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
void OhmicContactBC::MixA_DDM1_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  this->_current_buffer.clear();

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const PetscScalar current_scale = this->z_width()/A;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with ohmic boundary condition
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

            // for conduction current
            {
              PetscInt    ix[2] = {fvm_node->global_offset()+1, fvm_node->global_offset()+2};
              // I={In, Ip} the electron and hole current flow into this boundary cell.
              // NOTE: although In has dn/dt and R items, they are zero since n is const and n=n0 holds
              // so does Ip
              PetscScalar I[2];

              VecGetValues(f, 2, ix, I);

              // the current = In - Ip;
              this->_current_buffer.push_back((I[0] - I[1])*current_scale);
            }
            break;
          }
          case ElectrodeRegion:
          case InsulatorRegion:
          {
            PetscInt row = fvm_node->global_offset();
            clear_row.push_back(row);
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
 * build function and its jacobian for Mixed DDML1 solver
 */
void OhmicContactBC::MixA_DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Ohmic boundary condition is processed here.
  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // data buffer for mesh nodes
  std::vector<PetscInt> iy;
  std::vector<PetscScalar> y;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;

  // the electrode potential in current iteration
  PetscScalar Ve = x[ckt->local_offset_x(spice_node_index)];

  // search for all the boundary node
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
          // Semiconductor Region of course owns OhmicContactBC
          case SemiconductorRegion:
          {
            // semiconductor region should be the first region
            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();
            const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

            PetscScalar V = x[local_offset+0];  // psi of this node
            PetscScalar n = x[local_offset+1];  // electron density
            PetscScalar p = x[local_offset+2];  // hole density

            PetscScalar T = T_external();

            // mapping this node to material library
            semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

            PetscScalar nie = semi_region->material()->band->nie(p, n, T);
            PetscScalar Nc  = semi_region->material()->band->Nc(T);
            PetscScalar Nv  = semi_region->material()->band->Nv(T);
            PetscScalar Eg  = semi_region->material()->band->Eg(T);

            //governing equation for Ohmic contact boundary
            if(semi_region->get_advanced_model()->Fermi) //Fermi
            {
              PetscScalar Ec =  -(e*V + node_data->affinity() + semi_region->material()->band->EgNarrowToEc(p, n, T) + kb*T*log(Nc));
              PetscScalar Ev =  -(e*V + node_data->affinity() - semi_region->material()->band->EgNarrowToEv(p, n, T) - kb*T*log(Nv) + Eg);

              // the quasi-fermi potential equals to electrode Vapp
              PetscScalar phin = Ve;
              PetscScalar phip = Ve;

              PetscScalar etan = (-e*phin-Ec)/kb/T;
              PetscScalar etap = (Ev+e*phip)/kb/T;

              y.push_back( Nc*fermi_half(etan) - Nv*fermi_half(etap) -node_data->Net_doping() );
              y.push_back( n - Nc*fermi_half(etan) );
              y.push_back( p - Nv*fermi_half(etap) );

            }
            else     //Boltzmann
            {
              y.push_back  ( V - kb*T/e*boost::math::asinh(node_data->Net_doping()/(2*nie))
                             + Eg/(2*e)
                             + kb*T*log(Nc/Nv)/(2*e)
                             + node_data->affinity()/e
                             - Ve
                           );                        //governing equation for psi

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
              y.push_back( n - electron_density );  //governing equation for electron density
              y.push_back( p - hole_density );      //governing equation for hole density
            }

            // save insert position
            iy.push_back(global_offset+0);
            iy.push_back(global_offset+1);
            iy.push_back(global_offset+2);


            // here we should calculate current flow into this cell

            //for displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                FVM_Node *nb_node = (*nb_it).first;
                FVM_NodeData * nb_node_data = nb_node->node_data();
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
                this->_current_buffer.push_back( cv_boundary*node_data->eps()*dEdt*current_scale );
              }
            }
            break;
          }
          // conductor region which has an interface with OhmicContact boundary to semiconductor region
          case ElectrodeRegion:
          {

            // psi of this node
            PetscScalar V = x[local_offset];

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            //genius_assert( regions[0]->type()==SemiconductorRegion );
            PetscScalar V_semi = x[fvm_nodes[0]->local_offset()];

            // the psi of this node is equal to corresponding psi of semiconductor node
            y.push_back( V - V_semi );

            // save insert position
            iy.push_back(global_offset);

            break;
          }
          // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
          // the boundary for the corner point may be Ohmic. (not a nice behavier)
          case InsulatorRegion:
          {

            // psi of this node
            PetscScalar V = x[local_offset];

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            //genius_assert( regions[0]->type()==SemiconductorRegion );
            PetscScalar V_semi = x[fvm_nodes[0]->local_offset()];

            // the psi of this node is equal to corresponding psi of semiconductor node
            y.push_back( V - V_semi );

            // save insert position
            iy.push_back(global_offset);

            break;
          }
          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }
  }

  // prevent insert zero length vector
  if(iy.size())  VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);

  // for get the current, we must sum all the terms in current_buffer
  PetscScalar current = std::accumulate(this->_current_buffer.begin(), this->_current_buffer.end(), 0.0 );

  // Add current to spice node
  VecSetValue(f, this->global_offset(), current, ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for Mixed DDML1 solver
 */
void OhmicContactBC::MixA_DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{
  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);

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
            // insert none zero pattern
            // none zero pattern includes bd node and their neighbors!

            // bd node, psi = Ve
            MatSetValue(*jac, fvm_nodes[i]->global_offset(), this->global_offset(), 0, ADD_VALUES);

            break;
          }
          case ElectrodeRegion:
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


  // reserve jacobian entries for the circuit equation of ohmic electrode
  {

    std::vector<PetscInt> bc_node_reserve;
    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // get the derivative of electrode current to ohmic node
      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
      // but belong to different regions in logic.
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
      for(; rnode_it!=end_rnode_it; ++rnode_it  )
      {
        const SimulationRegion * region = (*rnode_it).second.first;
        const FVM_Node *  fvm_node = (*rnode_it).second.second;

        if ( region->type() == SemiconductorRegion)
        {
          bc_node_reserve.push_back(fvm_node->global_offset()+0);
          bc_node_reserve.push_back(fvm_node->global_offset()+1);
          bc_node_reserve.push_back(fvm_node->global_offset()+2);

          // get the derivative of electrode current to neighbors of ohmic node

          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
          FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
          for(; nb_it != nb_it_end; ++nb_it)
          {
            const FVM_Node *  fvm_nb_node = (*nb_it).first;
            bc_node_reserve.push_back(fvm_nb_node->global_offset()+0);
            bc_node_reserve.push_back(fvm_nb_node->global_offset()+1);
            bc_node_reserve.push_back(fvm_nb_node->global_offset()+2);
          }
        }
      }
    }
    Parallel::allgather(bc_node_reserve);

    if(Genius::is_last_processor())
    {
      PetscInt bc_global_offset = this->global_offset();

      MatSetValue(*jac, bc_global_offset, this->global_offset(), 0, ADD_VALUES);
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
 * do pre-process to jacobian matrix for  Mixed DDML1 solver
 */
void OhmicContactBC::MixA_DDM1_Jacobian_Preprocess(PetscScalar *, Mat *jac, std::vector<PetscInt> &src_row,
                                                   std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  _buffer_cols.clear();
  _buffer_jacobian_entries.clear();

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;

  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // only process nodes belong to this processor
    if( (*node_it)->processor_id() != Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(; rnode_it!=end_rnode_it; ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node *  fvm_node = (*rnode_it).second.second;

      if ( region->type() == SemiconductorRegion)
      {
        const FVM_NodeData * node_data = fvm_node->node_data();

        std::vector<PetscScalar> A1(3), A2(3), JM(3), JN(3);
        std::vector<PetscInt>    row(3);
        row[0] = fvm_node->global_offset()+0;
        row[1] = fvm_node->global_offset()+1;
        row[2] = fvm_node->global_offset()+2;

        //NOTE MatGetValues only get value from local block!
        MatGetValues(*jac, 1, &row[1], 3, &row[0], &A1[0]);
        MatGetValues(*jac, 1, &row[2], 3, &row[0], &A2[0]);
        JM[0] = (A1[0]-A2[0])*current_scale;
        JM[1] = (A1[1]-A2[1])*current_scale;
        JM[2] = (A1[2]-A2[2])*current_scale;


        // get the derivative of electrode current to neighbors of ohmic node
        // NOTE neighbors and ohmic bc node may on different processor!
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
        for(; nb_it != nb_it_end; ++nb_it)
        {
          const FVM_Node *  fvm_nb_node = (*nb_it).first;

          // distance from nb node to this node
          PetscScalar distance = fvm_node->distance(fvm_nb_node);
          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = fvm_node->cv_surface_area(fvm_nb_node);

          std::vector<PetscInt>    col(3);
          col[0] = fvm_nb_node->global_offset()+0;
          col[1] = fvm_nb_node->global_offset()+1;
          col[2] = fvm_nb_node->global_offset()+2;

          MatGetValues(*jac, 1, &row[1], 3, &col[0], &A1[0]);
          MatGetValues(*jac, 1, &row[2], 3, &col[0], &A2[0]);

          JN[0] = (A1[0]-A2[0])*current_scale;
          JN[1] = (A1[1]-A2[1])*current_scale;
          JN[2] = (A1[2]-A2[2])*current_scale;

          _buffer_cols.push_back(col);
          _buffer_jacobian_entries.push_back(JN);
        }

        _buffer_cols.push_back(row);
        _buffer_jacobian_entries.push_back(JM);
      }

    }
  }

  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {

    // only process nodes belong to this processor
    if( (*node_it)->processor_id() != Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with ohmic boundary condition
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
          case ElectrodeRegion:
          case InsulatorRegion:
          {
            PetscInt row = fvm_node->global_offset();
            clear_row.push_back(row);
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
 * build function and its jacobian for Mixed DDML1 solver
 */
void OhmicContactBC::MixA_DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  // the Jacobian of Ohmic boundary condition is processed here
  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);

  PetscInt bc_global_offset = this->global_offset();

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;

  // d(current)/d(independent variables of bd node and its neighbors)
  for(unsigned int n=0; n<_buffer_cols.size(); ++n)
  {
    MatSetValues(*jac, 1, &bc_global_offset, _buffer_cols[n].size(), &(_buffer_cols[n])[0], &(_buffer_jacobian_entries[n])[0], ADD_VALUES);
  }

  // after that, we set Jacobian entries for ohmic boundary condition
  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

      switch ( regions[i]->type() )
      {
          // Semiconductor Region of course owns OhmicContactBC
          case SemiconductorRegion:
          {
            const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

            //the indepedent variable number, we only need 4 here.
            adtl::AutoDScalar::numdir=4;
            //synchronize with material database
            semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

            AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];    V.setADValue(0, 1.0);  // psi of this node
            AutoDScalar n = x[fvm_nodes[i]->local_offset()+1];    n.setADValue(1, 1.0);  // electron density
            AutoDScalar p = x[fvm_nodes[i]->local_offset()+2];    p.setADValue(2, 1.0);  // hole density

            // the electrode potential in current iteration
            AutoDScalar Ve = x[ckt->local_offset_x(spice_node_index)];             Ve.setADValue(3, 1.0);

            PetscScalar T = T_external();

            semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

            AutoDScalar nie = semi_region->material()->band->nie(p, n, T);
            PetscScalar Nc  = semi_region->material()->band->Nc(T);
            PetscScalar Nv  = semi_region->material()->band->Nv(T);
            PetscScalar Eg  = semi_region->material()->band->Eg(T);

            //governing equation for Ohmic contact boundary
            AutoDScalar ff1,ff2,ff3;
            if(semi_region->get_advanced_model()->Fermi) //Fermi
            {
              AutoDScalar Ec =  -(e*V + node_data->affinity() + semi_region->material()->band->EgNarrowToEc(p, n, T) + kb*T*log(Nc));
              AutoDScalar Ev =  -(e*V + node_data->affinity() - semi_region->material()->band->EgNarrowToEv(p, n, T) - kb*T*log(Nv) + Eg);

              // the quasi-fermi potential equals to electrode Vapp
              AutoDScalar phin = Ve;
              AutoDScalar phip = Ve;

              AutoDScalar etan = (-e*phin-Ec)/kb/T;
              AutoDScalar etap = (Ev+e*phip)/kb/T;

              ff1 =  Nc*fermi_half(etan) - Nv*fermi_half(etap) -node_data->Net_doping();
              ff2 =  n - Nc*fermi_half(etan);
              ff3 =  p - Nv*fermi_half(etap);
            }
            else //Boltzmann
            {

              ff1 = V - kb*T/e*adtl::asinh(node_data->Net_doping()/(2*nie))
                    + Eg/(2*e)
                    + kb*T*log(Nc/Nv)/(2*e)
                    + node_data->affinity()/e
                    - Ve;

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

              ff2 =  n - electron_density;  //governing equation for electron density
              ff3 =  p - hole_density;      //governing equation for hole density
            }
            // the insert position
            PetscInt row[3], col[4];
            col[0] = row[0] = fvm_nodes[i]->global_offset()+0;
            col[1] = row[1] = fvm_nodes[i]->global_offset()+1;
            col[2] = row[2] = fvm_nodes[i]->global_offset()+2;
            col[3] = this->global_offset(); // the position of electrode equation

            // set Jacobian of governing equations
            MatSetValues(*jac, 1, &row[0], 4, &col[0], ff1.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[1], 4, &col[0], ff2.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[2], 4, &col[0], ff3.getADValue(), ADD_VALUES);

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

                MatSetValue(*jac, bc_global_offset, fvm_nodes[i]->global_offset()+0, current_disp.getADValue(0), ADD_VALUES);
                MatSetValue(*jac, bc_global_offset, nb_node->global_offset()+0, current_disp.getADValue(1), ADD_VALUES);
              }
            }

            break;
          }
          // conductor region which has an interface with OhmicContact boundary to semiconductor region
          case ElectrodeRegion:
          {
            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            AutoDScalar  V_semi = x[fvm_nodes[0]->local_offset()]; V_semi.setADValue(1,1.0);

            // the psi of this node is equal to corresponding psi of semiconductor node
            AutoDScalar  ff = V - V_semi;

            // set Jacobian of governing equation ff
            MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), ff.getADValue(1), ADD_VALUES);

            break;
          }
          // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
          // the boundary for the corner point may be Ohmic.
          case InsulatorRegion:
          {
            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
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




