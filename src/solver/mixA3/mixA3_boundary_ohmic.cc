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
#include "spice_ckt.h"
#include "parallel.h"
#include "mathfunc.h"
#include "petsc_utils.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::V;
using PhysicalUnit::A;



///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * build function and its jacobian for EBM solver
 */
void OhmicContactBC::MixA_EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Ohmic boundary condition is processed here.
  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);

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
  PetscScalar current_scale = this->z_width()/A;


  // the electrode potential in current iteration
  PetscScalar Ve = x[ckt->local_offset(spice_node_index)];

  // search for all the boundary node
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // region and fvm node buffer
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

      const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

      switch ( regions[i]->type() )
      {
        // Semiconductor Region of course owns OhmicContactBC
      case SemiconductorRegion:
        {
          // semiconductor region should be the first region
          genius_assert(i==0);

          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          // find the node variable offset
          unsigned int node_psi_offset = semi_region->ebm_variable_offset(POTENTIAL);
          unsigned int node_n_offset   = semi_region->ebm_variable_offset(ELECTRON);
          unsigned int node_p_offset   = semi_region->ebm_variable_offset(HOLE);
          unsigned int node_Tl_offset  = semi_region->ebm_variable_offset(TEMPERATURE);
          unsigned int node_Tn_offset  = semi_region->ebm_variable_offset(E_TEMP);
          unsigned int node_Tp_offset  = semi_region->ebm_variable_offset(H_TEMP);

          PetscScalar V = x[fvm_nodes[i]->local_offset() + node_psi_offset];  // psi of this node
          PetscScalar n = x[fvm_nodes[i]->local_offset() + node_n_offset];  // electron density
          PetscScalar p = x[fvm_nodes[i]->local_offset() + node_p_offset];  // hole density
          PetscScalar T = T_external();

          // lattice temperature
          if(semi_region->get_advanced_model()->enable_Tl())
            T = x[fvm_nodes[i]->local_offset() + node_Tl_offset];

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
          iy.push_back(fvm_nodes[i]->global_offset() + node_psi_offset);
          iy.push_back(fvm_nodes[i]->global_offset() + node_n_offset);
          iy.push_back(fvm_nodes[i]->global_offset() + node_p_offset);


          // process governing equation of T if required,
          // we should consider heat exchange to entironment if this ohmic bc is external boundary
          if( semi_region->get_advanced_model()->enable_Tl() && (node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion)) )
          {
            PetscScalar fT;
            PetscInt    iT = fvm_nodes[i]->global_offset()+node_Tl_offset;
            // get value of lattice temperature equatiuon
            VecGetValues(f, 1, &iT, &fT);
            // add heat flux out of ohmic boundary to lattice temperature equatiuon
            PetscScalar h = this->Heat_Transfer();
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            y.push_back( fT + h*(T_external()-T)*S );
            iy.push_back(iT);
          }

          // electron temperature if required
          if(semi_region->get_advanced_model()->enable_Tn())
          {
            PetscScalar Tn = x[fvm_nodes[i]->local_offset() + node_Tn_offset]/n;
            y.push_back(n*(Tn - T));
            iy.push_back(fvm_nodes[i]->global_offset() + node_Tn_offset);
          }

          // hole temperature if required
          if(semi_region->get_advanced_model()->enable_Tp())
          {
            PetscScalar Tp = x[fvm_nodes[i]->local_offset() + node_Tp_offset]/p;
            y.push_back(p*(Tp - T));
            iy.push_back(fvm_nodes[i]->global_offset() + node_Tp_offset);
          }

          // here we should calculate current flow into this cell

          // for conduction current
          {
            PetscInt    ix[2] = {fvm_nodes[i]->global_offset()+node_n_offset, fvm_nodes[i]->global_offset()+node_p_offset};
            // I={In, Ip} the electron and hole current flow into this boundary cell.
            PetscScalar I[2];

            VecGetValues(f, 2, ix, I);

            // the current = In - Ip;
            current_buffer.push_back((I[0] - I[1])*current_scale);
          }


          //for displacement current, only first order in time.
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

          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          unsigned int semiregion_node_psi_offset = regions[0]->ebm_variable_offset(POTENTIAL);
          unsigned int semiregion_node_Tl_offset  = regions[0]->ebm_variable_offset(TEMPERATURE);

          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()+node_psi_offset];

          // lattice temperature
          PetscScalar T = T_external();
          if(regions[i]->get_advanced_model()->enable_Tl())
            T = x[fvm_nodes[i]->local_offset() + node_Tl_offset];

          // since the region is sorted, we know region[0] is semiconductor region
          genius_assert( regions[0]->type()==SemiconductorRegion );

          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
          PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+semiregion_node_psi_offset];
          // the psi of this node is equal to corresponding psi of semiconductor node
          y.push_back( V - V_semi );
          iy.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            // let semiconductor node process the complete governing equation of heat transfer at interface
            // then we can set the node temperature of insulator/conductor region equal to node temperature at semiconductor region
            src_row.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);
            dst_row.push_back(fvm_nodes[0]->global_offset()+semiregion_node_Tl_offset);

            // the T of this node is equal to corresponding T of semiconductor node
            PetscScalar T_semi = x[fvm_nodes[0]->local_offset()+semiregion_node_Tl_offset];
            y.push_back( T - T_semi );
            iy.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);
          }

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

  //Add current to spice node
  VecAssemblyBegin(f);
  VecAssemblyEnd(f);
  if(Genius::is_last_processor())
    VecSetValue(f, this->global_offset(), current, ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for EBM3 solver
 */
void OhmicContactBC::MixA_EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.
  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);

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
          // insert none zero pattern
          // none zero pattern includes bd node and their neighbors!
          unsigned int n_node_var  = regions[i]->ebm_n_variables();
          unsigned int global_offset   = fvm_nodes[i]->global_offset();
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          for(unsigned int nv=0; nv<n_node_var; ++nv)
            MatSetValue(*jac, fvm_nodes[i]->global_offset()+nv, ckt->global_offset(spice_node_index), 0, ADD_VALUES);

          // reserve for heat transport equation
          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
            FVM_Node::fvm_ghost_node_iterator gn_it_end = fvm_nodes[i]->ghost_node_end();
            for(; gn_it != gn_it_end; ++gn_it)
            {
              const FVM_Node * ghost_fvm_node = (*gn_it).first;
              // skip NULL neighbor which means the node is on Neumann boundary
              if(ghost_fvm_node==NULL) continue;

              const SimulationRegion * ghost_region = this->system().region((*gn_it).second.first);
              genius_assert(ghost_region!=NULL);
              unsigned int ghostregion_node_Tl_offset  = ghost_region->ebm_variable_offset(TEMPERATURE);

              MatSetValue(*jac, global_offset+node_Tl_offset, ghost_fvm_node->global_offset()+ghostregion_node_Tl_offset, 0,ADD_VALUES);

              FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
              for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
                MatSetValue(*jac, global_offset+node_Tl_offset, (*gnb_it).second->global_offset()+ghostregion_node_Tl_offset, 0, ADD_VALUES);
            }
          }

          break;
        }
      case ConductorRegion:
      case InsulatorRegion:
        {
          unsigned int global_offset   = fvm_nodes[i]->global_offset();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          unsigned int semiregion_node_psi_offset = regions[0]->ebm_variable_offset(POTENTIAL);
          unsigned int semiregion_node_Tl_offset  = regions[0]->ebm_variable_offset(TEMPERATURE);

          // insert none zero pattern
          MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[0]->global_offset()+semiregion_node_psi_offset, 0, ADD_VALUES);

          if(regions[i]->get_advanced_model()->enable_Tl())
            MatSetValue(*jac, global_offset+node_Tl_offset,  fvm_nodes[0]->global_offset()+semiregion_node_Tl_offset, 0, ADD_VALUES);

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
      // reserve for displacement current
      const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, SemiconductorRegion);
      const SimulationRegion * region = get_fvm_node_region(*node_it, SemiconductorRegion);

      if(fvm_node->on_processor())
      {
        for(unsigned int nv=0; nv<region->ebm_n_variables(); ++nv)
          bc_node_reserve.push_back(fvm_node->global_offset()+nv);

        FVM_Node::fvm_neighbor_node_iterator nb_it     =  fvm_node->neighbor_node_begin();
        FVM_Node::fvm_neighbor_node_iterator nb_it_end =  fvm_node->neighbor_node_end();
        for(; nb_it!=nb_it_end; ++nb_it)
        {
          const FVM_Node *  fvm_nb_node = (*nb_it).second;
          for(unsigned int nv=0; nv<region->ebm_n_variables(); ++nv)
            bc_node_reserve.push_back(fvm_nb_node->global_offset()+nv);
        }
      }
    }
    Parallel::allgather(bc_node_reserve);

    if(Genius::processor_id() == Genius::n_processors()-1)
    {
      PetscInt bc_global_offset = this->global_offset();

      MatSetValue(*jac, bc_global_offset, ckt->global_offset(spice_node_index), 0, ADD_VALUES);

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
 * build function and its jacobian for EBM3 solver
 */
void OhmicContactBC::MixA_EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of Ohmic boundary condition is processed here
  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);
  PetscInt bc_global_offset = this->global_offset();

  // first, we should set Jacobian entries related with obmic electrode current
  // which is required by external electrode circuit
  MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);
  {
    std::vector< std::vector<PetscInt> >    cols;
    std::vector< std::vector<PetscScalar> > current_jacobian;

    // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
    PetscScalar current_scale = this->z_width()/A;

    BoundaryCondition::const_node_iterator node_it;
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      // get the derivative of electrode current to ohmic node
      const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, SemiconductorRegion);
      const FVM_NodeData * node_data = fvm_node->node_data();

      const SimulationRegion * region = get_fvm_node_region(*node_it, SemiconductorRegion);
      unsigned int n_node_var  = region->ebm_n_variables();
      unsigned int node_n_offset = region->ebm_variable_offset(ELECTRON);
      unsigned int node_p_offset = region->ebm_variable_offset(HOLE);

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

      std::vector<PetscScalar> A1(n_node_var), A2(n_node_var), JM(n_node_var), JN(n_node_var);
      std::vector<PetscInt>    row(n_node_var);
      for(unsigned int nv=0; nv<n_node_var; ++nv)
        row[nv] = fvm_node->global_offset()+nv;

      MatGetValues(*jac, 1, &row[node_n_offset], n_node_var, &row[0], &A1[0]);
      MatGetValues(*jac, 1, &row[node_p_offset], n_node_var, &row[0], &A2[0]);


      for(unsigned int nv=0; nv<n_node_var; ++nv)
        JM[nv] = (A1[nv]-A2[nv])*current_scale;


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

        std::vector<PetscInt>    col(n_node_var);
        for(unsigned int i=0; i<n_node_var; ++i)
          col[i] = fvm_nb_node->global_offset()+i;

        MatGetValues(*jac, 1, &row[node_n_offset], n_node_var, &col[0], &A1[0]);
        MatGetValues(*jac, 1, &row[node_p_offset], n_node_var, &col[0], &A2[0]);

        for(unsigned int nv=0; nv<n_node_var; ++nv)
          JN[nv] = (A1[nv]-A2[nv])*current_scale;
        // displacement current
        JM[0]+= cv_boundary*node_data->eps()/distance*disp_t*current_scale;
        JN[0]-= cv_boundary*node_data->eps()/distance*disp_t*current_scale;

        cols.push_back(col);
        current_jacobian.push_back(JN);
      }

      cols.push_back(row);
      current_jacobian.push_back(JM);
    }

    // d(current)/d(independent variables of bd node and its neighbors)
    for(unsigned int n=0; n<cols.size(); ++n)
      MatSetValues(*jac, 1, &bc_global_offset, cols[n].size(), &(cols[n])[0], &(current_jacobian[n])[0], ADD_VALUES);
  }


  // then, we can zero all the rows corresponding to Ohmic bc since they have been filled previously.
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

      std::vector<const SimulationRegion *> regions;
      std::vector<const FVM_Node *> fvm_nodes;

      // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
      // but in different region.
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
      for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
      {
        // should clear all the rows related with ohmic boundary condition
        regions.push_back( (*rnode_it).second.first );
        fvm_nodes.push_back( (*rnode_it).second.second );

        switch ( regions[i]->type() )
        {
        case SemiconductorRegion:
          {
            PetscInt row = fvm_nodes[i]->global_offset();

            row_for_clear.push_back(row+regions[i]->ebm_variable_offset(POTENTIAL));
            row_for_clear.push_back(row+regions[i]->ebm_variable_offset(ELECTRON));
            row_for_clear.push_back(row+regions[i]->ebm_variable_offset(HOLE));

            if(regions[i]->get_advanced_model()->enable_Tn())
              row_for_clear.push_back(row+regions[i]->ebm_variable_offset(E_TEMP));
            if(regions[i]->get_advanced_model()->enable_Tp())
              row_for_clear.push_back(row+regions[i]->ebm_variable_offset(H_TEMP));
            break;
          }
        case ConductorRegion:
        case InsulatorRegion:
          {
            PetscInt row = fvm_nodes[i]->global_offset();

            row_for_clear.push_back(row+regions[i]->ebm_variable_offset(POTENTIAL));

            if(regions[i]->get_advanced_model()->enable_Tl())
            {
              // if code reaches here, then the ohmic bc is an interface to conductor region
              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));

              row_for_clear.push_back(row+regions[i]->ebm_variable_offset(TEMPERATURE));
            }

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
          // semiconductor region should be the first region
          genius_assert(i==0);

          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          // find the node variable offset
          unsigned int n_node_var      = semi_region->ebm_n_variables();
          unsigned int node_psi_offset = semi_region->ebm_variable_offset(POTENTIAL);
          unsigned int node_n_offset   = semi_region->ebm_variable_offset(ELECTRON);
          unsigned int node_p_offset   = semi_region->ebm_variable_offset(HOLE);
          unsigned int node_Tl_offset  = semi_region->ebm_variable_offset(TEMPERATURE);
          unsigned int node_Tn_offset  = semi_region->ebm_variable_offset(E_TEMP);
          unsigned int node_Tp_offset  = semi_region->ebm_variable_offset(H_TEMP);

          //the indepedent variable number, we only need n_node_var plus electrode potential here.
          adtl::AutoDScalar::numdir = n_node_var + 1;
          //synchronize with material database
          semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

          AutoDScalar V = x[fvm_nodes[i]->local_offset()+node_psi_offset];  V.setADValue(node_psi_offset, 1.0);  // psi of this node
          AutoDScalar n = x[fvm_nodes[i]->local_offset()+node_n_offset];    n.setADValue(node_n_offset, 1.0);  // electron density
          AutoDScalar p = x[fvm_nodes[i]->local_offset()+node_p_offset];    p.setADValue(node_p_offset, 1.0);  // hole density

          AutoDScalar T  =  T_external();
          AutoDScalar Tn =  T_external();
          AutoDScalar Tp =  T_external();

          // lattice temperature if required
          if(semi_region->get_advanced_model()->enable_Tl())
          {
            T =  x[fvm_nodes[i]->local_offset() + node_Tl_offset];
            T.setADValue(node_Tl_offset, 1.0);
          }

          // electron temperature if required
          if(semi_region->get_advanced_model()->enable_Tn())
          {
            AutoDScalar nTn = x[fvm_nodes[i]->local_offset() + node_Tn_offset];
            nTn.setADValue(node_Tn_offset, 1.0);
            Tn = nTn/n;
          }

          // hole temperature if required
          if(semi_region->get_advanced_model()->enable_Tp())
          {
            AutoDScalar pTp = x[fvm_nodes[i]->local_offset() + node_Tp_offset];
            pTp.setADValue(node_Tp_offset, 1.0);
            Tp = pTp/p;
          }

          // the electrode potential in current iteration
          AutoDScalar Ve = x[ckt->local_offset(spice_node_index)];             Ve.setADValue(n_node_var, 1.0);

          // the insert position
          std::vector<PetscInt> row, col;
          for(unsigned int nv=0; nv<n_node_var; ++nv)
            row.push_back(fvm_nodes[i]->global_offset()+nv);
          col = row;
          col.push_back(ckt->global_offset(spice_node_index)); // the position of electrode equation

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

          // set Jacobian of governing equations of psi/n/p
          MatSetValues(*jac, 1, &row[node_psi_offset], col.size(), &col[0], ff1.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[node_n_offset],   col.size(), &col[0], ff2.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[node_p_offset],   col.size(), &col[0], ff3.getADValue(), ADD_VALUES);


          //governing equation of T, if this ohmic bc is external boundary, set heat flux here
          if(regions[i]->get_advanced_model()->enable_Tl() && (node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion)) )
          {
            PetscScalar h = this->Heat_Transfer();
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            AutoDScalar fT = h*(T_external()-T)*S;
            MatSetValues(*jac, 1, &row[node_Tl_offset], col.size(), &col[0], fT.getADValue(),  ADD_VALUES);
          }

          // electron temperature if required
          if(regions[i]->get_advanced_model()->enable_Tn())
          {
            AutoDScalar fTn = (n*(Tn - T));
            MatSetValues(*jac, 1, &row[node_Tn_offset], col.size(), &col[0], fTn.getADValue(),  ADD_VALUES);
          }

          // hole temperature if required
          if(regions[i]->get_advanced_model()->enable_Tp())
          {
            AutoDScalar fTp = (p*(Tp - T));
            MatSetValues(*jac, 1, &row[node_Tp_offset], col.size(), &col[0], fTp.getADValue(),  ADD_VALUES);
          }

          break;
        }
        // conductor region which has an interface with OhmicContact boundary to semiconductor region
      case ConductorRegion:
        // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
        // the boundary for the corner point may be Ohmic.
      case InsulatorRegion:
        {
          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          unsigned int semiregion_node_psi_offset = regions[0]->ebm_variable_offset(POTENTIAL);
          unsigned int semiregion_node_Tl_offset  = regions[0]->ebm_variable_offset(TEMPERATURE);

          {
            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()+node_psi_offset];  V.setADValue(0,1.0);
            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            AutoDScalar  V_semi = x[fvm_nodes[0]->local_offset()+semiregion_node_psi_offset]; V_semi.setADValue(1,1.0);
            // the psi of this node is equal to corresponding psi of semiconductor node
            AutoDScalar  ff1 = V - V_semi;
            PetscInt row = fvm_nodes[i]->global_offset()+node_psi_offset;
            PetscInt col[2] = {fvm_nodes[i]->global_offset()+node_psi_offset, fvm_nodes[0]->global_offset()+semiregion_node_psi_offset};
            MatSetValues(*jac, 1, &row, 2, &col[0], ff1.getADValue(), ADD_VALUES);
          }

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            AutoDScalar  T = x[fvm_nodes[i]->local_offset()+node_Tl_offset]; T.setADValue(0,1.0); // lattice temperature of this node
            AutoDScalar  T_semi = x[fvm_nodes[0]->local_offset()+semiregion_node_Tl_offset]; T_semi.setADValue(1,1.0);
            // the T of this node is equal to corresponding T of semiconductor node
            AutoDScalar  ff2 = T - T_semi;
            PetscInt row = fvm_nodes[i]->global_offset()+node_Tl_offset;
            PetscInt col[2] = {fvm_nodes[i]->global_offset()+node_Tl_offset, fvm_nodes[0]->global_offset()+semiregion_node_Tl_offset};
            MatSetValues(*jac, 1, &row, 2, &col[0], ff2.getADValue(), ADD_VALUES);
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






