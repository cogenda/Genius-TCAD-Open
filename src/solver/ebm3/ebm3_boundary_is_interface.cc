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
#include "insulator_region.h"
#include "boundary_condition_is.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for EBM3 solver
 */
void InsulatorSemiconductorInterfaceBC::EBM3_Function_Preprocess(PetscScalar *,Vec f, std::vector<PetscInt> &src_row,
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

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
            dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
            clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
          }
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
 * build function and its jacobian for EBM3 solver
 */
void InsulatorSemiconductorInterfaceBC::EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Insulator-Semiconductor interface is processed here

  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // buffer for Vec location
  std::vector<PetscInt> iy;

  // buffer for Vec value
  std::vector<PetscScalar> y_new;

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
        // Insulator-Semiconductor interface at Semiconductor side, do nothing
      case SemiconductorRegion:
        {
          // semiconductor region should be the first region
          genius_assert(i==0);

          const SemiconductorSimulationRegion * sregion = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();
          genius_assert(node_data);

          unsigned int node_psi_offset = sregion->ebm_variable_offset(POTENTIAL);
          unsigned int node_n_offset   = sregion->ebm_variable_offset(ELECTRON);
          unsigned int node_p_offset   = sregion->ebm_variable_offset(HOLE);
          unsigned int node_Tl_offset  = sregion->ebm_variable_offset(TEMPERATURE);
          unsigned int node_Tn_offset  = sregion->ebm_variable_offset(E_TEMP);
          unsigned int node_Tp_offset  = sregion->ebm_variable_offset(H_TEMP);

          PetscScalar n   =  x[fvm_nodes[i]->local_offset()+node_n_offset];       // electron density
          PetscScalar p   =  x[fvm_nodes[i]->local_offset()+node_p_offset];       // hole density
          PetscScalar T   = T_external();
          PetscScalar Tn  = T_external();
          PetscScalar Tp  = T_external();

          // lattice temperature if required
          if(sregion->get_advanced_model()->enable_Tl())
            T =  x[fvm_nodes[i]->local_offset() + node_Tl_offset];

          // electron temperature if required
          if(sregion->get_advanced_model()->enable_Tn())
            Tn = x[fvm_nodes[i]->local_offset() + node_Tn_offset]/n;

          // hole temperature if required
          if(sregion->get_advanced_model()->enable_Tp())
            Tp = x[fvm_nodes[i]->local_offset() + node_Tp_offset]/p;

          PetscScalar Eg = sregion->material()->band->Eg(T);

          // process interface fixed charge density
          PetscScalar boundary_area = fvm_nodes[i]->outside_boundary_surface_area();
          VecSetValue(f, fvm_nodes[i]->global_offset()+node_psi_offset, qf*boundary_area, ADD_VALUES);

          {
            // surface recombination
            Material::MaterialSemiconductor *mt =  sregion->material();
            mt->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);
            mt->band->nie(p, n, T);
            PetscScalar GSurf = - mt->band->R_Surf(p, n, T) * boundary_area; //generation due to SRH

            VecSetValue(f, fvm_nodes[i]->global_offset()+node_n_offset, GSurf, ADD_VALUES);
            VecSetValue(f, fvm_nodes[i]->global_offset()+node_p_offset, GSurf, ADD_VALUES);

            if (sregion->get_advanced_model()->enable_Tn())
            {
              PetscScalar Hn = 1.5*kb*Tn*GSurf;
              VecSetValue(f, fvm_nodes[i]->global_offset()+node_Tn_offset, Hn, ADD_VALUES);
            }
            if (sregion->get_advanced_model()->enable_Tp())
            {
              PetscScalar Hp = 1.5*kb*Tn*GSurf;
              VecSetValue(f, fvm_nodes[i]->global_offset()+node_Tp_offset, Hp, ADD_VALUES);
            }
            if (sregion->get_advanced_model()->enable_Tl())
            {
              PetscScalar Eg = mt->band->Eg(T);
              PetscScalar H  = - GSurf*(Eg+1.5*kb*Tn+1.5*kb*Tp);
              VecSetValue(f, fvm_nodes[i]->global_offset()+node_Tl_offset, H, ADD_VALUES);
            }
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

            VecSetValue(f, fvm_nodes[i]->global_offset()+node_n_offset, -TrapElec, ADD_VALUES);
            VecSetValue(f, fvm_nodes[i]->global_offset()+node_p_offset, -TrapHole, ADD_VALUES);

            if(sregion->get_advanced_model()->enable_Tn())
              VecSetValue(f, fvm_nodes[i]->global_offset()+node_Tn_offset, -1.5*kb*Tn*TrapHole, ADD_VALUES);
            if(sregion->get_advanced_model()->enable_Tp())
              VecSetValue(f, fvm_nodes[i]->global_offset()+node_Tp_offset, -1.5*kb*Tp*TrapHole, ADD_VALUES);
            if(sregion->get_advanced_model()->enable_Tl())
            {
              PetscScalar EcEi = 0.5*Eg - kb*T*log(node_data->Nc()/node_data->Nv());
              PetscScalar EiEv = 0.5*Eg + kb*T*log(node_data->Nc()/node_data->Nv());
              PetscScalar H    = sregion->material()->trap->TrapHeat(false,p,n,ni,Tp,Tn,T,EcEi,EiEv);
              VecSetValue(f, fvm_nodes[i]->global_offset()+node_Tl_offset, H*boundary_area, ADD_VALUES);
            }
          }

          break;
        }

        // Insulator-Semiconductor interface at Insulator side, we get the previous f for insulator region,
        // as well as surface charge density, plus to corresponding function of semiconductor. and then
        // force the potential equal to corresponding point in semiconductor
      case InsulatorRegion:
        {

          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          // the governing equation of this fvm node

          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()+0];

          // since the region is sorted, we know region[0] is semiconductor region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
          //genius_assert( regions[0]->type()==SemiconductorRegion );
          PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+0];

          // the psi of this node is equal to corresponding psi of semiconductor node
          // since psi should be continuous for the interface
          PetscScalar ff1 = V - V_semi;

          // find the position ff will be add to
          iy.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);
          y_new.push_back(ff1);


          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            // T of this node
            PetscScalar T = x[fvm_nodes[i]->local_offset()+node_Tl_offset];

            // T for corresponding semiconductor region
            PetscScalar T_semi = x[fvm_nodes[0]->local_offset()+regions[0]->ebm_variable_offset(TEMPERATURE)];

            // the T of this node is equal to corresponding T of semiconductor node
            // by assuming no heat resistance between 2 region
            PetscScalar ff2 = T - T_semi;

            iy.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);
            y_new.push_back(ff2);
          }

          break;
        }
      case VacuumRegion:
        break;
      default: genius_error(); //we should never reach here
      }
    }

  }

  // insert new value
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y_new[0]), ADD_VALUES);

  add_value_flag = ADD_VALUES;

  // gate current
  _gate_current_function(x, f, add_value_flag);
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for EBM3 solver
 */
void InsulatorSemiconductorInterfaceBC::EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  // search for all the node with this boundary type
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

    std::vector<SimulationRegion *> regions;
    std::vector<FVM_Node *> fvm_nodes;

    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
        // Insulator-Semiconductor interface at Semiconductor side, we should reserve entrance for later add operator
      case SemiconductorRegion:
        {
          // semiconductor region should be the first region
          genius_assert(i==0);

          unsigned int global_offset   = fvm_nodes[i]->global_offset();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
          for( ; gn_it!=fvm_nodes[i]->ghost_node_end(); ++gn_it )
          {
            // ghost node is semiconductor region node
            const FVM_Node * ghost_fvm_node = (*gn_it).first;

            const SimulationRegion * ghost_region = this->system().region((*gn_it).second.first);
            unsigned int ghostregion_node_psi_offset = ghost_region->ebm_variable_offset(POTENTIAL);
            unsigned int ghostregion_node_Tl_offset  = ghost_region->ebm_variable_offset(TEMPERATURE);

            MatSetValue(*jac, global_offset+node_psi_offset, ghost_fvm_node->global_offset()+ghostregion_node_psi_offset, 0, ADD_VALUES);
            if(regions[i]->get_advanced_model()->enable_Tl())
              MatSetValue(*jac, global_offset+node_Tl_offset, ghost_fvm_node->global_offset()+ghostregion_node_Tl_offset, 0, ADD_VALUES);

            FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
            for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
            {
              MatSetValue(*jac, global_offset+node_psi_offset, (*gnb_it).first->global_offset()+ghostregion_node_psi_offset, 0, ADD_VALUES);

              if(regions[i]->get_advanced_model()->enable_Tl())
                MatSetValue(*jac, global_offset+node_Tl_offset, (*gnb_it).first->global_offset()+ghostregion_node_Tl_offset, 0, ADD_VALUES);
            }

          }

          break;
        }

        // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
      case InsulatorRegion:
        {

          unsigned int global_offset   = fvm_nodes[i]->global_offset();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          unsigned int semiconductor_node_psi_offset = regions[0]->ebm_variable_offset(POTENTIAL);
          unsigned int semiconductor_node_Tl_offset  = regions[0]->ebm_variable_offset(TEMPERATURE);

          // reserve for later operator
          MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[0]->global_offset()+semiconductor_node_psi_offset, 0, ADD_VALUES);

          if(regions[i]->get_advanced_model()->enable_Tl())
            MatSetValue(*jac, global_offset+node_Tl_offset, fvm_nodes[0]->global_offset()+semiconductor_node_Tl_offset, 0, ADD_VALUES);

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
  _gate_current_jacobian_reserve(jac, add_value_flag);

}




/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for EBM3 solver
 */
void InsulatorSemiconductorInterfaceBC::EBM3_Jacobian_Preprocess(PetscScalar * ,Mat *jac, std::vector<PetscInt> &src_row,
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

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
            dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
            clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
          }
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
 * build function and its jacobian for EBM3 solver
 */
void InsulatorSemiconductorInterfaceBC::EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  // after that, set values to source rows
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
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
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
        // Insulator-Semiconductor interface at Semiconductor side, do nothing
      case SemiconductorRegion:
        {
          // process interface traps
          const SemiconductorSimulationRegion * sregion = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          unsigned int n_node_var      = sregion->ebm_n_variables();
          unsigned int node_psi_offset = sregion->ebm_variable_offset(POTENTIAL);
          unsigned int node_n_offset   = sregion->ebm_variable_offset(ELECTRON);
          unsigned int node_p_offset   = sregion->ebm_variable_offset(HOLE);
          unsigned int node_Tl_offset  = sregion->ebm_variable_offset(TEMPERATURE);
          unsigned int node_Tn_offset  = sregion->ebm_variable_offset(E_TEMP);
          unsigned int node_Tp_offset  = sregion->ebm_variable_offset(H_TEMP);


          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();
          genius_assert(node_data);


          {
            //the indepedent variable number
            adtl::AutoDScalar::numdir = n_node_var;

            Material::MaterialSemiconductor *mt =  sregion->material();

            //synchronize with material database
            mt->set_ad_num(adtl::AutoDScalar::numdir);
            mt->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

            std::vector<PetscInt> index;
            for(unsigned int nv=0; nv<n_node_var; ++nv)  index.push_back( fvm_nodes[i]->global_offset()+nv );


            AutoDScalar n = x[fvm_nodes[i]->local_offset() + node_n_offset];
            n.setADValue(node_n_offset, 1.0);     // electron density

            AutoDScalar p = x[fvm_nodes[i]->local_offset() + node_p_offset];
            p.setADValue(node_p_offset, 1.0);     // hole density

            AutoDScalar T  =  T_external();
            AutoDScalar Tn =  T_external();
            AutoDScalar Tp =  T_external();

            // lattice temperature if required
            if(sregion->get_advanced_model()->enable_Tl())
            {
              T =  x[fvm_nodes[i]->local_offset() + node_Tl_offset];
              T.setADValue(node_Tl_offset, 1.0);
            }

            // electron temperature if required
            if(sregion->get_advanced_model()->enable_Tn())
            {
              AutoDScalar nTn = x[fvm_nodes[i]->local_offset() + node_Tn_offset];
              nTn.setADValue(node_Tn_offset, 1.0);
              Tn = nTn/n;
            }

            // hole temperature if required
            if(sregion->get_advanced_model()->enable_Tp())
            {
              AutoDScalar pTp = x[fvm_nodes[i]->local_offset() + node_Tp_offset];
              pTp.setADValue(node_Tp_offset, 1.0);
              Tp = pTp/p;
            }
            AutoDScalar Eg = mt->band->Eg(T);
            AutoDScalar ni = sregion->material()->band->nie(p, n, T);

            PetscScalar boundary_area = fvm_nodes[i]->outside_boundary_surface_area();

            // surface recombination
            AutoDScalar GSurf = - mt->band->R_Surf(p, n, T) * boundary_area; //generation due to SRH

            MatSetValues(*jac, 1, &index[node_n_offset], n_node_var, &index[0], GSurf.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &index[node_p_offset], n_node_var, &index[0], GSurf.getADValue(), ADD_VALUES);

            if (sregion->get_advanced_model()->enable_Tn())
            {
              AutoDScalar Hn = 1.5*kb*Tn*GSurf;
              MatSetValues(*jac, 1, &index[node_Tn_offset], n_node_var, &index[0], Hn.getADValue(), ADD_VALUES);
            }
            if (sregion->get_advanced_model()->enable_Tp())
            {
              AutoDScalar Hp = 1.5*kb*Tn*GSurf;
              MatSetValues(*jac, 1, &index[node_Tp_offset], n_node_var, &index[0], Hp.getADValue(), ADD_VALUES);
            }
            if (sregion->get_advanced_model()->enable_Tl())
            {
              AutoDScalar H  = - GSurf*(Eg+1.5*kb*Tn+1.5*kb*Tp);
              MatSetValues(*jac, 1, &index[node_Tl_offset], n_node_var, &index[0], H.getADValue(), ADD_VALUES);
            }

            if (sregion->get_advanced_model()->Trap)
            {
              mt->trap->Calculate(false,p,n,ni,T);
              AutoDScalar TrappedC = mt->trap->ChargeAD(false) * boundary_area;
              MatSetValues(*jac, 1, &index[node_psi_offset], n_node_var, &index[0], TrappedC.getADValue(), ADD_VALUES);

              AutoDScalar TrapElec = mt->trap->ElectronTrapRate(false,n,ni,T);
              AutoDScalar TrapHole = mt->trap->HoleTrapRate    (false,p,ni,T);
              MatSetValues(*jac, 1, &index[node_n_offset], n_node_var, &index[0], (-TrapElec*boundary_area).getADValue(), ADD_VALUES);
              MatSetValues(*jac, 1, &index[node_p_offset], n_node_var, &index[0], (-TrapElec*boundary_area).getADValue(), ADD_VALUES);

              if(sregion->get_advanced_model()->enable_Tn())
              {
                AutoDScalar Hn = - 1.5 * kb*Tn * TrapElec;
                MatSetValues(*jac, 1, &index[node_Tn_offset], n_node_var, &index[0], (Hn*boundary_area).getADValue(),   ADD_VALUES);
              }

              if(sregion->get_advanced_model()->enable_Tp())
              {
                AutoDScalar Hp = - 1.5 * kb*Tp * TrapHole;
                MatSetValues(*jac, 1, &index[node_Tp_offset], n_node_var, &index[0], (Hp*boundary_area).getADValue(),   ADD_VALUES);
              }

              if(sregion->get_advanced_model()->enable_Tl())
              {
                AutoDScalar EcEi = 0.5*Eg - kb*T*log(node_data->Nc()/node_data->Nv());
                AutoDScalar EiEv = 0.5*Eg + kb*T*log(node_data->Nc()/node_data->Nv());
                AutoDScalar H = mt->trap->TrapHeat(false,p,n,ni,Tp,Tn,T,EcEi,EiEv);
                MatSetValues(*jac, 1, &index[node_Tl_offset], n_node_var, &index[0], (H*boundary_area).getADValue(),   ADD_VALUES);
              }
            }
          }
          break;
        }
        // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
      case InsulatorRegion:
        {
          unsigned int global_offset   = fvm_nodes[i]->global_offset();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          unsigned int semiconductor_node_psi_offset = regions[0]->ebm_variable_offset(POTENTIAL);
          unsigned int semiconductor_node_Tl_offset  = regions[0]->ebm_variable_offset(TEMPERATURE);

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // psi of this node
          AutoDScalar  V = x[fvm_nodes[i]->local_offset()+node_psi_offset]; V.setADValue(0,1.0);

          // since the region is sorted, we know region[0] is semiconductor region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
          //genius_assert( regions[0]->type()==SemiconductorRegion );
          AutoDScalar  V_semi = x[fvm_nodes[0]->local_offset()+semiconductor_node_psi_offset]; V_semi.setADValue(1,1.0);

          // the psi of this node is equal to corresponding psi of semiconductor node
          AutoDScalar  ff1 = V - V_semi;

          // set Jacobian of governing equation ff
          MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[i]->global_offset()+node_psi_offset, ff1.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[0]->global_offset()+semiconductor_node_psi_offset, ff1.getADValue(1), ADD_VALUES);

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            // T of this node
            AutoDScalar  T = x[fvm_nodes[i]->local_offset()+node_Tl_offset]; T.setADValue(0,1.0);

            // T for corresponding semiconductor region
            AutoDScalar  T_semi = x[fvm_nodes[0]->local_offset()+semiconductor_node_Tl_offset]; T_semi.setADValue(1,1.0);

            // the T of this node is equal to corresponding T of semiconductor node
            // we assuming T is continuous for the interface
            AutoDScalar ff2 = T - T_semi;

            // set Jacobian of governing equation ff2
            MatSetValue(*jac, global_offset+node_Tl_offset, fvm_nodes[i]->global_offset()+node_Tl_offset, ff2.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, global_offset+node_Tl_offset, fvm_nodes[0]->global_offset()+semiconductor_node_Tl_offset, ff2.getADValue(1), ADD_VALUES);
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

  // gate current
  _gate_current_jacobian(x, jac, add_value_flag);
}
