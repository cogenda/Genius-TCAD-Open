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

// fill AC matrix of node on Insulator - Semiconductor Interface
// the AC matrix should always keep synchronous with Jacobian matrix of EBM

void InsulatorSemiconductorInterfaceBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag )
{

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

          const SemiconductorSimulationRegion * sregion = dynamic_cast<const SemiconductorSimulationRegion *>(regions[0]);

          /**
           * load Jacobian entry of this node from J and fill into AC matrix A
           */
          sregion->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], A, b, J, omega, add_value_flag);


          /**
           * process interface traps
           */
          const FVM_NodeData * node_data = fvm_nodes[0]->node_data();
          sregion->material()->mapping(fvm_nodes[0]->root_node(), node_data, SolverSpecify::clock);

          unsigned int n_node_var      = sregion->ebm_n_variables();
          unsigned int node_psi_offset = sregion->ebm_variable_offset(POTENTIAL);
          unsigned int node_n_offset   = sregion->ebm_variable_offset(ELECTRON);
          unsigned int node_p_offset   = sregion->ebm_variable_offset(HOLE);
          unsigned int node_Tl_offset  = sregion->ebm_variable_offset(TEMPERATURE);
          unsigned int node_Tn_offset  = sregion->ebm_variable_offset(E_TEMP);
          unsigned int node_Tp_offset  = sregion->ebm_variable_offset(H_TEMP);

          if (sregion->get_advanced_model()->Trap)
          {
            //the indepedent variable number
            adtl::AutoDScalar::numdir = n_node_var;

            //synchronize with material database
            sregion->material()->set_ad_num(adtl::AutoDScalar::numdir);

            std::vector<PetscInt> index_re, index_im;
            for(unsigned int nv=0; nv<n_node_var; ++nv)
            {
              index_re.push_back( fvm_nodes[i]->global_offset()+nv );
              index_im.push_back( fvm_nodes[i]->global_offset()+n_node_var+nv );
            }

            AutoDScalar n = node_data->n();
            n.setADValue(node_n_offset, 1.0);     // electron density

            AutoDScalar p = node_data->p();
            p.setADValue(node_p_offset, 1.0);     // hole density

            AutoDScalar T  =  T_external();
            AutoDScalar Tn =  T_external();
            AutoDScalar Tp =  T_external();

            // lattice temperature if required
            if(sregion->get_advanced_model()->enable_Tl())
            {
              T =  node_data->T();
              T.setADValue(node_Tl_offset, 1.0);
            }

            // electron temperature if required
            if(sregion->get_advanced_model()->enable_Tn())
            {
              AutoDScalar nTn = node_data->n()*node_data->Tn();
              nTn.setADValue(node_Tn_offset, 1.0);
              Tn = nTn/n;
            }

            // hole temperature if required
            if(sregion->get_advanced_model()->enable_Tp())
            {
              AutoDScalar pTp = node_data->p()*node_data->Tp();
              pTp.setADValue(node_Tp_offset, 1.0);
              Tp = pTp/p;
            }


            AutoDScalar Eg = sregion->material()->band->Eg(T);
            PetscScalar boundary_area = fvm_nodes[i]->outside_boundary_surface_area();
            AutoDScalar ni = sregion->material()->band->nie(p, n, T);
            sregion->material()->trap->Calculate(false,p,n,ni,T);
            AutoDScalar TrappedC = sregion->material()->trap->ChargeAD(false) * boundary_area;
            MatSetValues(A, 1, &index_re[node_psi_offset], n_node_var, &index_re[0], TrappedC.getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &index_im[node_psi_offset], n_node_var, &index_im[0], TrappedC.getADValue(), ADD_VALUES);


            AutoDScalar TrapElec = sregion->material()->trap->ElectronTrapRate(false,n,ni,T);
            AutoDScalar TrapHole = sregion->material()->trap->HoleTrapRate    (false,p,ni,T);
            MatSetValues(A, 1, &index_re[node_n_offset], n_node_var, &index_re[0], (-TrapElec*boundary_area).getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &index_im[node_n_offset], n_node_var, &index_im[0], (-TrapElec*boundary_area).getADValue(), ADD_VALUES);

            MatSetValues(A, 1, &index_re[node_p_offset], n_node_var, &index_re[0], (-TrapElec*boundary_area).getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &index_im[node_p_offset], n_node_var, &index_im[0], (-TrapElec*boundary_area).getADValue(), ADD_VALUES);


            if(sregion->get_advanced_model()->enable_Tn())
            {
              AutoDScalar Hn = - 1.5 * kb*Tn * TrapElec;
              MatSetValues(A, 1, &index_re[node_Tn_offset], n_node_var, &index_re[0], (Hn*boundary_area).getADValue(),   ADD_VALUES);
              MatSetValues(A, 1, &index_im[node_Tn_offset], n_node_var, &index_im[0], (Hn*boundary_area).getADValue(),   ADD_VALUES);
            }

            if(sregion->get_advanced_model()->enable_Tp())
            {
              AutoDScalar Hp = - 1.5 * kb*Tp * TrapHole;
              MatSetValues(A, 1, &index_re[node_Tp_offset], n_node_var, &index_re[0], (Hp*boundary_area).getADValue(),   ADD_VALUES);
              MatSetValues(A, 1, &index_im[node_Tp_offset], n_node_var, &index_im[0], (Hp*boundary_area).getADValue(),   ADD_VALUES);
            }

            if(sregion->get_advanced_model()->enable_Tl())
            {
              AutoDScalar EcEi = 0.5*Eg - kb*T*log(node_data->Nc()/node_data->Nv());
              AutoDScalar EiEv = 0.5*Eg + kb*T*log(node_data->Nc()/node_data->Nv());
              AutoDScalar H = sregion->material()->trap->TrapHeat(true,p,n,ni,Tp,Tn,T,EcEi,EiEv);
              MatSetValues(A, 1, &index_re[node_Tl_offset], n_node_var, &index_re[0], (H*boundary_area).getADValue(),   ADD_VALUES);
              MatSetValues(A, 1, &index_im[node_Tl_offset], n_node_var, &index_im[0], (H*boundary_area).getADValue(),   ADD_VALUES);
            }

          }
          break;
        }
        // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
      case InsulatorRegion:
        {


          /**
           * load Jacobian entry from J and fill into AC matrix A.
           * Please note, we add the entries to location of semiconductor node
           */
          regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], A, b, J, omega, add_value_flag, regions[0], fvm_nodes[0]);


          /**
           *  let psi and T of node in this region equal to psi(T) of node in the semiconductor region
           */
          unsigned int n_variables     = regions[i]->ebm_n_variables();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          {
            // psi of this node
            AutoDScalar  V = fvm_nodes[i]->node_data()->psi(); V.setADValue(0,1.0);

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            //genius_assert( regions[0]->type()==SemiconductorRegion );
            AutoDScalar  V_semi = fvm_nodes[0]->node_data()->psi(); V_semi.setADValue(1,1.0);

            // the psi of this node is equal to corresponding psi of insulator node in the other region
            AutoDScalar  ff1 = V - V_semi;

            // set Jacobian of governing equation ff1
            PetscInt real_row = fvm_nodes[i]->global_offset() + node_psi_offset;
            PetscInt imag_row = fvm_nodes[i]->global_offset() + n_variables + node_psi_offset;
            PetscInt real_col[2]={fvm_nodes[i]->global_offset() + node_psi_offset,
                                  fvm_nodes[0]->global_offset() + regions[0]->ebm_variable_offset(POTENTIAL)};
            PetscInt imag_col[2]={fvm_nodes[i]->global_offset() + n_variables + node_psi_offset,
                                  fvm_nodes[0]->global_offset() + regions[0]->ebm_n_variables() + regions[0]->ebm_variable_offset(POTENTIAL)};

            MatSetValues(A, 1, &real_row, 2, real_col, ff1.getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &imag_row, 2, imag_col, ff1.getADValue(), ADD_VALUES);
          }

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            // T of this node
            AutoDScalar  T = fvm_nodes[i]->node_data()->T(); T.setADValue(0,1.0);

            // T of insulator node in the other region
            AutoDScalar  T_semi = fvm_nodes[0]->node_data()->T(); T_semi.setADValue(1,1.0);

            // the T of this node is equal to corresponding T of insulator node in the other region
            AutoDScalar ff2 = T - T_semi;

            // set Jacobian of governing equation ff1
            PetscInt real_row = fvm_nodes[i]->global_offset() + node_Tl_offset;
            PetscInt imag_row = fvm_nodes[i]->global_offset() + n_variables + node_Tl_offset;
            PetscInt real_col[2]={fvm_nodes[i]->global_offset() + node_Tl_offset,
                                  fvm_nodes[0]->global_offset() + regions[0]->ebm_variable_offset(TEMPERATURE)};
            PetscInt imag_col[2]={fvm_nodes[i]->global_offset() + n_variables + node_Tl_offset,
                                  fvm_nodes[0]->global_offset() + regions[0]->ebm_n_variables() + regions[0]->ebm_variable_offset(TEMPERATURE)};

            MatSetValues(A, 1, &real_row, 2, real_col, ff2.getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &imag_row, 2, imag_col, ff2.getADValue(), ADD_VALUES);
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
