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
#include "boundary_condition.h"
#include "parallel.h"
#include "mathfunc.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;



///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML2 solver
 */
void OhmicContactBC::Mix_DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Ohmic boundary condition is processed here.

  // since we will use INSERT_VALUES operat, check the vec state.
  if( (add_value_flag != INSERT_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }
  // we should do two things here. one is performance ohmic bc to corresponding mesh nodes.
  // another problem is the ohmic bc has an extra external circuit equation. we must deal with it.


  // note, we will use VecGetValues to get some values from vec f
  // as a result, we should assembly the vec first.
  VecAssemblyBegin(f);
  VecAssemblyEnd(f);

  // data buffer for add row to row location
  std::vector<PetscInt> src_row;
  std::vector<PetscInt> dst_row;

  // data buffer for mesh nodes
  std::vector<PetscInt> iy;
  std::vector<PetscScalar> y;

  // the electrode current, since the electrode may be partitioned into several processor,
  // we should collect it.
  std::vector<double> current_buffer;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();

  // the electrode potential in current iteration
  PetscScalar Ve = ext_circuit()->Vapp();

  // search for all the boundary node
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

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
          // semiconductor region should be the first region
          genius_assert(i==0);

          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          PetscScalar V = x[fvm_nodes[i]->local_offset()+0];  // psi of this node
          PetscScalar n = x[fvm_nodes[i]->local_offset()+1];  // electron density
          PetscScalar p = x[fvm_nodes[i]->local_offset()+2];  // hole density
          PetscScalar T = x[fvm_nodes[i]->local_offset()+3];  // lattice temperature

          // mapping this node to material library
          semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

          PetscScalar nie = semi_region->material()->band->nie(p, n, T);
          PetscScalar Nc  = semi_region->material()->band->Nc(T);
          PetscScalar Nv  = semi_region->material()->band->Nv(T);
          PetscScalar Eg  = semi_region->material()->band->Eg(T);

          //governing equation of psi/n/p for Ohmic contact boundary
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
                           + node_data->affinity()
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
          iy.push_back(fvm_nodes[i]->global_offset()+0);
          iy.push_back(fvm_nodes[i]->global_offset()+1);
          iy.push_back(fvm_nodes[i]->global_offset()+2);



          /*
           * process governing equation of T, which should consider heat exchange to entironment
           */

          // if this ohmic bc is external boundary, set heat flux here
          if(node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion))
          {
            PetscScalar fT;
            PetscInt    iT = fvm_nodes[i]->global_offset()+3;
            // get value of lattice temperature equatiuon
            VecGetValues(f, 1, &iT, &fT);
            // add heat flux out of ohmic boundary to lattice temperature equatiuon
            PetscScalar h = this->Heat_Transfer();
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            y.push_back( fT + h*(T_external()-T)*S );
            iy.push_back(iT);
          }


          // here we should calculate current flow into this cell

          // for conduction current
          {
            PetscInt    ix[2] = {fvm_nodes[i]->global_offset()+1, fvm_nodes[i]->global_offset()+2};
            // I={In, Ip} the electron and hole current flow into this boundary cell.
            PetscScalar I[2];

            VecGetValues(f, 2, ix, I);

            // the current = In - Ip;
            current_buffer.push_back((I[0] - I[1])*current_scale);
          }



          //for displacement current
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
        // conductor region which has an interface with OhmicContact boundary to semiconductor region
      case ConductorRegion:
        // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
        // the boundary for the corner point may be Ohmic. (not a nice behavier)
      case InsulatorRegion:
        {

          PetscScalar V = x[fvm_nodes[i]->local_offset()+0]; // psi of this node
          PetscScalar T = x[fvm_nodes[i]->local_offset()+1]; // lattice temperature of this node

          // since the region is sorted, we know region[0] is semiconductor region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
          genius_assert( regions[0]->type()==SemiconductorRegion );
          PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+0];
          PetscScalar T_semi = x[fvm_nodes[0]->local_offset()+3];

          // let semiconductor node process the complete governing equation of heat transfer at interface
          // then we can set the node temperature of insulator/conductor region equal to node temperature at semiconductor region
          src_row.push_back(fvm_nodes[i]->global_offset()+1);
          dst_row.push_back(fvm_nodes[0]->global_offset()+3);

          // the psi of this node is equal to corresponding psi of semiconductor node
          y.push_back( V - V_semi );
          // the T of this node is equal to corresponding T of semiconductor node
          y.push_back( T - T_semi );

          // save insert position
          iy.push_back(fvm_nodes[i]->global_offset()+0);
          iy.push_back(fvm_nodes[i]->global_offset()+1);

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

  // prevent insert zero length vector
  if(iy.size())  VecSetValues(f, iy.size(), &iy[0], &y[0], INSERT_VALUES);

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
void OhmicContactBC::Mix_DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
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
          genius_assert(i==0);

          FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
          FVM_Node::fvm_ghost_node_iterator gn_it_end = fvm_nodes[i]->ghost_node_end();
          for(; gn_it != gn_it_end; ++gn_it)
          {
            const FVM_Node * ghost_fvm_node = (*gn_it).first;
            // skip NULL neighbor which means the node is on Neumann boundary
            if(ghost_fvm_node==NULL) continue;

            MatSetValue(*jac, fvm_nodes[i]->global_offset()+3, ghost_fvm_node->global_offset()+1, 0, ADD_VALUES);
          }

          break;
        }
      case ConductorRegion:
        {
          // insert none zero pattern
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, fvm_nodes[0]->global_offset()+0, 0, ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[0]->global_offset()+3, 0, ADD_VALUES);

          break;
        }
      case InsulatorRegion:
        {
          // insert none zero pattern
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, fvm_nodes[0]->global_offset()+0, 0, ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[0]->global_offset()+3, 0, ADD_VALUES);

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
void OhmicContactBC::Mix_DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of Ohmic boundary condition is processed here

  // first, we zero all the rows corresponding to ohmic bc.
  {

    // data buffer for add row to rwo location
    std::vector<PetscInt> src_row;
    std::vector<PetscInt> dst_row;

    // the indicator which rows should be set to zero
    std::vector<PetscInt> row_for_clear;

    BoundaryCondition::const_node_iterator node_it;
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      std::vector<const FVM_Node *> fvm_nodes;

      // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
      // but in different region.
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

      // should clear all the rows related with ohmic boundary condition
      for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
      {
        const SimulationRegion * region = (*rnode_it).second.first;
        fvm_nodes.push_back((*rnode_it).second.second);

        switch ( region->type() )
        {
        case SemiconductorRegion:
          {
            row_for_clear.push_back(fvm_nodes[i]->global_offset()+0);
            row_for_clear.push_back(fvm_nodes[i]->global_offset()+1);
            row_for_clear.push_back(fvm_nodes[i]->global_offset()+2);
            break;
          }
        case ConductorRegion:
        case InsulatorRegion:
          {
            // if code reaches here, then the ohmic bc is an interface to conductor region

            src_row.push_back(fvm_nodes[i]->global_offset()+1);
            dst_row.push_back(fvm_nodes[0]->global_offset()+3);

            row_for_clear.push_back(fvm_nodes[i]->global_offset()+0);
            row_for_clear.push_back(fvm_nodes[i]->global_offset()+1);
            break;
          }
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

  // after that, we set Jacobian entries for ohmic boundary condition
  BoundaryCondition::const_node_iterator node_it;
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

      const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

      switch ( regions[i]->type() )
      {
        // Semiconductor Region of course owns OhmicContactBC
      case SemiconductorRegion:
        {
          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);
          // semiconductor region should be the first region
          genius_assert(i==0);

          //the indepedent variable number, we only need 4 here.
          adtl::AutoDScalar::numdir=4;
          //synchronize with material database
          semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

          AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];    V.setADValue(0, 1.0);  // psi of this node
          AutoDScalar n = x[fvm_nodes[i]->local_offset()+1];    n.setADValue(1, 1.0);  // electron density
          AutoDScalar p = x[fvm_nodes[i]->local_offset()+2];    p.setADValue(2, 1.0);  // hole density
          AutoDScalar T = x[fvm_nodes[i]->local_offset()+3];    T.setADValue(3, 1.0);  // lattice temperature

          // the insert position
          std::vector<PetscInt> row, col;
          row.push_back(fvm_nodes[i]->global_offset()+0);
          row.push_back(fvm_nodes[i]->global_offset()+1);
          row.push_back(fvm_nodes[i]->global_offset()+2);
          row.push_back(fvm_nodes[i]->global_offset()+3);
          col = row;

          semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);


          AutoDScalar nie = semi_region->material()->band->nie(p, n, T);
          AutoDScalar Nc  = semi_region->material()->band->Nc(T);
          AutoDScalar Nv  = semi_region->material()->band->Nv(T);
          AutoDScalar Eg  = semi_region->material()->band->Eg(T);



          //governing equation of pis/n/p for Ohmic contact boundary
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
                  + node_data->affinity()
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

          // set Jacobian of governing equations
          MatSetValues(*jac, 1, &row[0], col.size(), &col[0], ff1.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[1], col.size(), &col[0], ff2.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[2], col.size(), &col[0], ff3.getADValue(), ADD_VALUES);

          //governing equation of T
          // if this ohmic bc is external boundary, set heat flux here
          if(node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion))
          {
            PetscScalar h = this->Heat_Transfer();
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            AutoDScalar fT = h*(T_external()-T)*S;
            MatSetValues(*jac, 1, &row[3], col.size(), &col[0], fT.getADValue(),  ADD_VALUES);
          }

          break;
        }
        // conductor region which has an interface with OhmicContact boundary to semiconductor region
      case ConductorRegion:
        // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
        // the boundary for the corner point may be Ohmic.
      case InsulatorRegion:
        {
          //the indepedent variable number, we need 4 here.
          adtl::AutoDScalar::numdir=4;

          AutoDScalar  V = x[fvm_nodes[i]->local_offset()+0]; V.setADValue(0,1.0); // psi of this node
          AutoDScalar  T = x[fvm_nodes[i]->local_offset()+1]; T.setADValue(1,1.0); // lattice temperature of this node

          // since the region is sorted, we know region[0] is semiconductor region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
          AutoDScalar  V_semi = x[fvm_nodes[0]->local_offset()+0]; V_semi.setADValue(2,1.0);
          AutoDScalar  T_semi = x[fvm_nodes[0]->local_offset()+3]; T_semi.setADValue(3,1.0);

          // the psi of this node is equal to corresponding psi of semiconductor node
          AutoDScalar  ff1 = V - V_semi;
          // the T of this node is equal to corresponding T of semiconductor node
          AutoDScalar  ff2 = T - T_semi;

          // set Jacobian of governing equation ff1 and ff2
          std::vector<PetscInt> rows,cols;
          rows.push_back(fvm_nodes[i]->global_offset()+0);
          rows.push_back(fvm_nodes[i]->global_offset()+1);
          cols = rows;
          cols.push_back(fvm_nodes[0]->global_offset()+0);
          cols.push_back(fvm_nodes[0]->global_offset()+3);

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


  // the last operator is INSERT_VALUES
  add_value_flag = INSERT_VALUES;

}




void OhmicContactBC::Mix_DDM2_Electrode_Load(const Vec , const Mat * jac, double & I, PetscScalar & pdI_pdV, Vec & pdI_pdx, Vec & pdF_pdV)
{
  I = this->ext_circuit()->current();
  pdI_pdV = 0;

  VecZeroEntries(pdI_pdx);
  VecZeroEntries(pdF_pdV);


  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();


  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // get the derivative of electrode current to ohmic node
    const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, SemiconductorRegion);
    const FVM_NodeData * node_data = fvm_node->node_data();

    // displacement current
    PetscScalar disp_t;
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false) //second order
    {
      PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
      disp_t =  (2-r)/(1-r)/(SolverSpecify::dt_last+SolverSpecify::dt);
    }
    else//first order
    {
      disp_t = 1.0/SolverSpecify::dt;
    }

    std::vector<PetscScalar> A1(4), A2(4), JM(4), JN(4);
    std::vector<PetscInt>    row(4);
    row[0] = fvm_node->global_offset()+0;
    row[1] = fvm_node->global_offset()+1;
    row[2] = fvm_node->global_offset()+2;
    row[3] = fvm_node->global_offset()+3;

    MatGetValues(*jac, 1, &row[1], 4, &row[0], &A1[0]);
    MatGetValues(*jac, 1, &row[2], 4, &row[0], &A2[0]);

    JM[0] = (A1[0]-A2[0])*current_scale;
    JM[1] = (A1[1]-A2[1])*current_scale;
    JM[2] = (A1[2]-A2[2])*current_scale;
    JM[3] = (A1[3]-A2[3])*current_scale;



    // get the derivative of electrode current to neighbors of ohmic node

    FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
    FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();

    for(; nb_it != nb_it_end; ++nb_it)
    {
      const FVM_Node *  fvm_nb_node = (*nb_it).second;
      // distance from nb node to this node
      PetscScalar distance = (*(fvm_node->root_node()) - *(fvm_nb_node->root_node())).size();
      // area of out surface of control volume related with neighbor node
      PetscScalar cv_boundary = fvm_node->cv_surface_area(fvm_nb_node->root_node());

      std::vector<PetscInt>    col(4);
      col[0] = fvm_nb_node->global_offset()+0;
      col[1] = fvm_nb_node->global_offset()+1;
      col[2] = fvm_nb_node->global_offset()+2;
      col[3] = fvm_nb_node->global_offset()+3;

      MatGetValues(*jac, 1, &row[1], 4, &col[0], &A1[0]);
      MatGetValues(*jac, 1, &row[2], 4, &col[0], &A2[0]);

      JN[0] = (A1[0]-A2[0])*current_scale;
      JN[1] = (A1[1]-A2[1])*current_scale;
      JN[2] = (A1[2]-A2[2])*current_scale;
      JN[3] = (A1[3]-A2[3])*current_scale;

      // displacement current
      JM[0]+= cv_boundary*node_data->eps()/distance*disp_t*current_scale;
      JN[0]-= cv_boundary*node_data->eps()/distance*disp_t*current_scale;

      VecSetValues( pdI_pdx, col.size(), &col[0], &JN[0], ADD_VALUES);
    }

    VecSetValues( pdI_pdx, row.size(), &row[0], &JM[0], ADD_VALUES);

    {
      VecSetValue( pdF_pdV, row[0], 1.0, ADD_VALUES);
    }

  }

  VecAssemblyBegin(pdI_pdx);
  VecAssemblyBegin(pdF_pdV);

  VecAssemblyEnd(pdI_pdx);
  VecAssemblyEnd(pdF_pdV);
}

