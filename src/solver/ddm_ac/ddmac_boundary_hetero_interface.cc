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
#include "boundary_condition_hetero.h"
#include "petsc_utils.h"
#include "mathfunc.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;

void HeteroInterfaceBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const PetscScalar omega, InsertMode & add_value_flag )
{
  //the indepedent variable number, we need max 12 here.
  adtl::AutoDScalar::numdir=12;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // the variable for first region
    AutoDScalar V0;
    AutoDScalar n0;
    AutoDScalar p0;
    AutoDScalar T0;
    AutoDScalar Tn0;
    AutoDScalar Tp0;
    Material::MaterialSemiconductor * mt0;
    AutoDScalar Ec0,Ev0;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      if(!fvm_node->is_valid()) continue;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0)
      {
        /*
         * load Jacobian entry of this node from J and fill into AC matrix A
         */
        regions[0]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], POTENTIAL, A, b, J, omega, add_value_flag);
        regions[0]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], ELECTRON, A, b, J, omega, add_value_flag);
        regions[0]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], HOLE, A, b, J, omega, add_value_flag);

        if(regions[i]->get_advanced_model()->enable_Tl())
          regions[0]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], TEMPERATURE, A, b, J, omega, add_value_flag);

        // calculate variables for first region
        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

        const FVM_NodeData * n0_data = fvm_nodes[i]->node_data();

        unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
        unsigned int node_n_offset   = regions[i]->ebm_variable_offset(ELECTRON);
        unsigned int node_p_offset   = regions[i]->ebm_variable_offset(HOLE);
        unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);
        unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset(E_TEMP);
        unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset(H_TEMP);

        V0 = fvm_nodes[i]->node_data()->psi();  V0.setADValue(node_psi_offset,1.0);
        n0 = fvm_nodes[i]->node_data()->n();    n0.setADValue(node_n_offset,1.0);  // electron density
        p0 = fvm_nodes[i]->node_data()->p();    p0.setADValue(node_p_offset,1.0);  // hole density

        if(regions[i]->get_advanced_model()->enable_Tl())
        {
          T0 = fvm_nodes[i]->node_data()->T();
          T0.setADValue(node_Tl_offset,1.0);
        }
        else
          T0 = T_external();

        if(regions[i]->get_advanced_model()->enable_Tn())
        {
          Tn0 = fvm_nodes[i]->node_data()->Tn();
          Tn0.setADValue(node_Tn_offset, 1.0);
        }
        else
          Tn0 = T0;

        if(regions[i]->get_advanced_model()->enable_Tp())
        {
          Tp0 = fvm_nodes[i]->node_data()->Tp();
          Tp0.setADValue(node_Tp_offset, 1.0);
        }
        else
          Tp0 = T0;

        mt0 = semi_region->material();
        mt0->mapping(fvm_nodes[i]->root_node(), n0_data, SolverSpecify::clock);

        Ec0 =  -(e*V0 + n0_data->affinity() + mt0->band->EgNarrowToEc(p0, n0, T0) + kb*T0*log(n0_data->Nc()));
        Ev0 =  -(e*V0 + n0_data->affinity() - mt0->band->EgNarrowToEv(p0, n0, T0) - kb*T0*log(n0_data->Nv()) + mt0->band->Eg(T0));
        if(semi_region->get_advanced_model()->Fermi)
        {
          Ec0 = Ec0 - kb*T0*log(gamma_f(fabs(n0)/n0_data->Nc()));
          Ev0 = Ev0 + kb*T0*log(gamma_f(fabs(p0)/n0_data->Nv()));
        }

      }

      // other region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);
              const FVM_NodeData * n_data = fvm_nodes[i]->node_data();

              unsigned int n_node_var_0    = regions[0]->ebm_n_variables();
              unsigned int n_node_var      = regions[i]->ebm_n_variables();
              unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
              unsigned int node_n_offset   = regions[i]->ebm_variable_offset(ELECTRON);
              unsigned int node_p_offset   = regions[i]->ebm_variable_offset(HOLE);
              unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);
              unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset(E_TEMP);
              unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset(H_TEMP);


              AutoDScalar V = fvm_nodes[i]->node_data()->psi();  V.setADValue(n_node_var_0+node_psi_offset, 1.0); // psi of this node
              AutoDScalar n = fvm_nodes[i]->node_data()->n();    n.setADValue(n_node_var_0+node_n_offset, 1.0);   // electron density
              AutoDScalar p = fvm_nodes[i]->node_data()->p();    p.setADValue(n_node_var_0+node_p_offset, 1.0);   // hole density

              AutoDScalar T  =  T_external();
              AutoDScalar Tn =  T_external();
              AutoDScalar Tp =  T_external();

              // lattice temperature if required
              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                T =  fvm_nodes[i]->node_data()->T();
                T.setADValue(n_node_var_0+node_Tl_offset, 1.0);
              }

              // electron temperature if required
              if(regions[i]->get_advanced_model()->enable_Tn())
              {
                Tn = fvm_nodes[i]->node_data()->Tn();
                Tn.setADValue(n_node_var_0+node_Tn_offset, 1.0);
              }

              // hole temperature if required
              if(regions[i]->get_advanced_model()->enable_Tp())
              {
                Tp = fvm_nodes[i]->node_data()->Tp();
                Tp.setADValue(n_node_var_0+node_Tp_offset, 1.0);
              }

              // mapping this node to material library
              Material::MaterialSemiconductor *mt = semi_region->material();
              mt->mapping(fvm_nodes[i]->root_node(), n_data, SolverSpecify::clock);
              AutoDScalar Ec =  -(e*V + n_data->affinity() + mt->band->EgNarrowToEc(p, n, T) + kb*T*log(n_data->Nc()));
              AutoDScalar Ev =  -(e*V + n_data->affinity() - mt->band->EgNarrowToEv(p, n, T) - kb*T*log(n_data->Nv()) + mt->band->Eg(T));
              if(semi_region->get_advanced_model()->Fermi)
              {
                Ec = Ec - kb*T*log(gamma_f(fabs(n)/n_data->Nc()));
                Ev = Ev + kb*T*log(gamma_f(fabs(p)/n_data->Nv()));
              }

              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_nodes[i]->outside_boundary_surface_area();

              {

                regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], POTENTIAL, A, b, J, omega, add_value_flag, regions[0], fvm_nodes[0]);
                regions[i]->DDMAC_Force_equal(fvm_nodes[i], POTENTIAL, A, add_value_flag, regions[0], fvm_nodes[0]);
              }

              {

                regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], ELECTRON, A, b, J, omega, add_value_flag);
                regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], HOLE, A, b, J, omega, add_value_flag);
              }


              std::vector<int> rows_re, rows_im, cols_re, cols_im;
              for(unsigned int nv=0; nv<n_node_var_0; ++nv)
              {
                rows_re.push_back(fvm_nodes[0]->global_offset()+nv);
                rows_im.push_back(fvm_nodes[0]->global_offset()+n_node_var_0+nv);
              }
              for(unsigned int nv=0; nv<n_node_var; ++nv)
              {
                rows_re.push_back(fvm_nodes[i]->global_offset()+nv);
                rows_im.push_back(fvm_nodes[i]->global_offset()+n_node_var+nv);
              }
              cols_re = rows_re;
              cols_im = rows_im;

              // thermal emit current
              AutoDScalar Jn=0, Jp=0;
              AutoDScalar Sn=0, Sp=0;

              if(Ec0 > Ec)
              {
                // Jn is enter region 0
                AutoDScalar pm = mt0->band->EffecElecMass(T)/mt->band->EffecElecMass(T);
                Jn =  2*e*(mt0->band->ThermalVn(T)*n0 - pm*mt->band->ThermalVn(T)*n*exp(-(Ec0-Ec)/(kb*T)))*cv_boundary;
                Sn = -2*e*(mt0->band->ThermalVn(T)*n0*2.5*kb*Tn0 - pm*mt->band->ThermalVn(T)*n*exp(-(Ec0-Ec)/(kb*T))*2.5*kb*Tn)*cv_boundary;
              }
              else
              {
                // Jn is left region 0
                AutoDScalar pm = mt->band->EffecElecMass(T)/mt0->band->EffecElecMass(T);
                Jn = -2*e*(mt->band->ThermalVn(T)*n - pm*mt0->band->ThermalVn(T)*n0*exp(-(Ec-Ec0)/(kb*T)))*cv_boundary;
                Sn =  2*e*(mt->band->ThermalVn(T)*n*2.5*kb*Tn - pm*mt0->band->ThermalVn(T)*n0*exp(-(Ec-Ec0)/(kb*T))*2.5*kb*Tn0)*cv_boundary;
              }

              if(Ev0 < Ev)
              {
                // Jp is enter region 0
                AutoDScalar pm = mt0->band->EffecHoleMass(T)/mt->band->EffecHoleMass(T);
                Jp = -2*e*(mt0->band->ThermalVp(T)*p0 - pm*mt->band->ThermalVp(T)*p*exp((Ev0-Ev)/(kb*T)))*cv_boundary;
                Sp =  2*e*(mt0->band->ThermalVp(T)*p0*2.5*kb*Tp0 - pm*mt->band->ThermalVp(T)*p*exp((Ev0-Ev)/(kb*T))*2.5*kb*Tp)*cv_boundary;
              }
              else
              {
                AutoDScalar pm = mt->band->EffecHoleMass(T)/mt0->band->EffecHoleMass(T);
                Jp =  2*e*(mt->band->ThermalVp(T)*p - pm*mt0->band->ThermalVp(T)*p0*exp((Ev-Ev0)/(kb*T)))*cv_boundary;
                Sp = -2*e*(mt->band->ThermalVp(T)*p*2.5*kb*Tp - pm*mt0->band->ThermalVp(T)*p0*exp((Ev-Ev0)/(kb*T))*2.5*kb*Tp0)*cv_boundary;
              }


              // current flow
              MatSetValues(A, 1, &rows_re[n_node_var_0+node_n_offset], cols_re.size(), &cols_re[0], (Jn).getADValue(), ADD_VALUES);
              MatSetValues(A, 1, &rows_re[node_n_offset], cols_re.size(), &cols_re[0], (-Jn).getADValue(), ADD_VALUES);
              MatSetValues(A, 1, &rows_im[n_node_var_0+node_n_offset], cols_im.size(), &cols_im[0], (Jn).getADValue(), ADD_VALUES);
              MatSetValues(A, 1, &rows_im[node_n_offset], cols_im.size(), &cols_im[0], (-Jn).getADValue(), ADD_VALUES);

              MatSetValues(A, 1, &rows_re[n_node_var_0+node_p_offset], cols_re.size(), &cols_re[0], (-Jp).getADValue(), ADD_VALUES);
              MatSetValues(A, 1, &rows_re[node_p_offset], cols_re.size(), &cols_re[0], (Jp).getADValue(), ADD_VALUES);
              MatSetValues(A, 1, &rows_im[n_node_var_0+node_p_offset], cols_im.size(), &cols_im[0], (-Jp).getADValue(), ADD_VALUES);
              MatSetValues(A, 1, &rows_im[node_p_offset], cols_im.size(), &cols_im[0], (Jp).getADValue(), ADD_VALUES);



              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], TEMPERATURE, A, b, J, omega, add_value_flag, regions[0], fvm_nodes[0]);
                regions[i]->DDMAC_Force_equal(fvm_nodes[i], TEMPERATURE, A, add_value_flag, regions[0], fvm_nodes[0]);
              }

              // FIXME Sn, Sp

            }


            case InsulatorRegion:
            {
              //load Jacobian entry of this node from J and fill into AC matrix A with the location of fvm_node[0]
              regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], A, b, J, omega, add_value_flag, regions[0], fvm_nodes[0]);

              // let indepedent variable of node in this region equal to indepedent variable of node in the first semiconductor region
              regions[i]->DDMAC_Force_equal(fvm_nodes[i], POTENTIAL, A, add_value_flag, regions[0], fvm_nodes[0]);
              if(regions[i]->get_advanced_model()->enable_Tl())
                regions[i]->DDMAC_Force_equal(fvm_nodes[i], TEMPERATURE, A, add_value_flag, regions[0], fvm_nodes[0]);
              break;
            }
            default: genius_error();
        }
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}

