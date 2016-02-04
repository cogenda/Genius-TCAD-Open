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
#include "solver_specify.h"
#include "boundary_condition_hetero.h"
#include "petsc_utils.h"
#include "mathfunc.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * do pre-process to function for DDML1 solver
 */
void HeteroInterfaceBC::DDM1_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
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
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0) continue;

      // other  region
      else
      {
        switch( region->type() )
        {
          case SemiconductorRegion :
          {
              // record the source row and dst row
            src_row.push_back(fvm_nodes[i]->global_offset()+0);
            dst_row.push_back(fvm_nodes[0]->global_offset()+0);
            clear_row.push_back(fvm_nodes[i]->global_offset()+0);
            break;
          }
          case InsulatorRegion:
          {
              // record the source row and dst row
            src_row.push_back(fvm_nodes[i]->global_offset()+0);
            dst_row.push_back(fvm_nodes[0]->global_offset()+0);
            clear_row.push_back(fvm_nodes[i]->global_offset()+0);
            break;
          }
          default: genius_error();
        }
      }
    }

  }

}


//FIXME how to compute current pass through hetero junction?
// Medici says there are thermal emit and tunneling current
// However, dessis only considers thermal emit current.
// what shoud I do?

/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void HeteroInterfaceBC::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  const PetscScalar T = T_external();
  const PetscScalar Vt  = kb*T/e;
  const PetscScalar qf = this->scalar("qf");

  // buffer for Vec value
  std::vector<PetscInt> iy;
  std::vector<PetscScalar> y_new;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // the variable for first region
    const FVM_Node * fvm_node0;
    PetscScalar V0;
    PetscScalar n0;
    PetscScalar p0;
    Material::MaterialSemiconductor * mt0;
    PetscScalar Ec0,Ev0;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      // the first semiconductor region
      if(i==0)
      {
        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
        const FVM_NodeData * n0_data = fvm_node->node_data();

        fvm_node0 = fvm_node;
        V0 = x[fvm_node->local_offset()+0];
        n0 = x[fvm_node->local_offset()+1];  // electron density
        p0 = x[fvm_node->local_offset()+2];  // hole density

        mt0 = semi_region->material();
        mt0->mapping(fvm_node->root_node(), n0_data, SolverSpecify::clock);

        Ec0 =  -(e*V0 + n0_data->affinity() + mt0->band->EgNarrowToEc(p0, n0, T) + kb*T*log(n0_data->Nc()));
        Ev0 =  -(e*V0 + n0_data->affinity() - mt0->band->EgNarrowToEv(p0, n0, T) - kb*T*log(n0_data->Nv()) + mt0->band->Eg(T));
        if(semi_region->get_advanced_model()->Fermi)
        {
          Ec0 = Ec0 - e*Vt*log(gamma_f(fabs(n0)/n0_data->Nc()));
          Ev0 = Ev0 + e*Vt*log(gamma_f(fabs(p0)/n0_data->Nv()));
        }

        PetscScalar boundary_area = std::abs(fvm_node->outside_boundary_surface_area());
        VecSetValue(f, fvm_node->global_offset(), qf*boundary_area, ADD_VALUES);
      }

      // other semiconductor region
      else
      {
        // the ghost node should have the same processor_id with me
        genius_assert(fvm_node->root_node()->processor_id() == fvm_node0->root_node()->processor_id() );

        switch( region->type() )
        {
            case SemiconductorRegion :
            {

              const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
              const FVM_NodeData * n_data = fvm_node->node_data();

              PetscScalar V = x[fvm_node->local_offset()+0];  // psi of this node
              PetscScalar n = x[fvm_node->local_offset()+1];  // electron density
              PetscScalar p = x[fvm_node->local_offset()+2];  // hole density

              // mapping this node to material library
              Material::MaterialSemiconductor *mt = semi_region->material();
              mt->mapping(fvm_node->root_node(), n_data, SolverSpecify::clock);
              PetscScalar Ec =  -(e*V + n_data->affinity() + mt->band->EgNarrowToEc(p, n, T) + kb*T*log(n_data->Nc()));
              PetscScalar Ev =  -(e*V + n_data->affinity() - mt->band->EgNarrowToEv(p, n, T) - kb*T*log(n_data->Nv()) + mt->band->Eg(T));
              if(semi_region->get_advanced_model()->Fermi)
              {
                Ec = Ec - e*Vt*log(gamma_f(fabs(n)/n_data->Nc()));
                Ev = Ev + e*Vt*log(gamma_f(fabs(p)/n_data->Nv()));
              }


              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              PetscScalar ff1 = V - V0;
              iy.push_back(fvm_node->global_offset()+0);
              y_new.push_back(ff1);

              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_node->outside_boundary_surface_area();

              // thermal emit current
              if(Ec0 > Ec)
              {
                // electrons are leaving region 0
                PetscScalar pm = mt0->band->EffecElecMass(T)/mt->band->EffecElecMass(T);
                PetscScalar Jn = 2.*(mt0->band->ThermalVn(T)*n0 - pm*mt->band->ThermalVn(T)*(n*exp(-(Ec0-Ec)/(kb*T))))*cv_boundary;
                VecSetValue(f, fvm_node->global_offset()+1,   Jn, ADD_VALUES);
                VecSetValue(f, fvm_node0->global_offset()+1, -Jn, ADD_VALUES);
              }
              else
              {
                // electrons are leaving this region
                PetscScalar pm = mt->band->EffecElecMass(T)/mt0->band->EffecElecMass(T);
                PetscScalar Jn = 2.*(mt->band->ThermalVn(T)*n - pm*mt0->band->ThermalVn(T)*(n0*exp(-(Ec-Ec0)/(kb*T))))*cv_boundary;
                VecSetValue(f, fvm_node->global_offset()+1,  -Jn, ADD_VALUES);
                VecSetValue(f, fvm_node0->global_offset()+1,  Jn, ADD_VALUES);
              }

              if(Ev0 < Ev)
              {
                // holes are leaving region 0
                PetscScalar pm = mt0->band->EffecHoleMass(T)/mt->band->EffecHoleMass(T);
                PetscScalar Jp = 2.*(mt0->band->ThermalVp(T)*p0 - pm*mt->band->ThermalVp(T)*(p*exp((Ev0-Ev)/(kb*T))))*cv_boundary;
                VecSetValue(f, fvm_node->global_offset()+2,   Jp, ADD_VALUES);
                VecSetValue(f, fvm_node0->global_offset()+2, -Jp, ADD_VALUES);
              }
              else
              {
                // holes are leaving this region
                PetscScalar pm = mt->band->EffecHoleMass(T)/mt0->band->EffecHoleMass(T);
                PetscScalar Jp = 2.*(mt->band->ThermalVp(T)*p - pm*mt0->band->ThermalVp(T)*(p0*exp((Ev-Ev0)/(kb*T))))*cv_boundary;
                VecSetValue(f, fvm_node->global_offset()+2,  -Jp, ADD_VALUES);
                VecSetValue(f, fvm_node0->global_offset()+2,  Jp, ADD_VALUES);
              }
              break;
            }


            case InsulatorRegion:
            {

              // the governing equation of this fvm node
              PetscScalar V = x[fvm_node->local_offset()+0];  // psi of this node
              PetscScalar V_semi = x[fvm_node0->local_offset()+0];

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              PetscScalar ff1 = V - V_semi;

              iy.push_back(fvm_node->global_offset()+0);
              y_new.push_back(ff1);
              break;
            }
            default: genius_error();
        }
      }
    }
  }


  // set new value
  if( iy.size() )
    VecSetValues(f, iy.size(), &iy[0], &y_new[0], ADD_VALUES);

  add_value_flag = ADD_VALUES;


}







/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void HeteroInterfaceBC::DDM1_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
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
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0) continue;

      // other  region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+0);
              break;
            }
            case InsulatorRegion:
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+0);
              break;
            }
            default: genius_error();
        }
      }
    }

  }

}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void HeteroInterfaceBC::DDM1_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{

  const PetscScalar T = T_external();
  const PetscScalar Vt  = kb*T/e;
  //the indepedent variable number, we need 6 here.
  adtl::AutoDScalar::numdir=6;

  // after that, set values to source rows
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // the variable for first region
    const FVM_Node * fvm_node0;
    AutoDScalar V0;
    AutoDScalar n0;
    AutoDScalar p0;
    Material::MaterialSemiconductor * mt0;
    AutoDScalar Ec0,Ev0;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    for(unsigned int i=0; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      // the first semiconductor region
      if(i==0)
      {
        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
        const FVM_NodeData * n0_data = fvm_node->node_data();

        mt0 = semi_region->material();
        mt0->set_ad_num(adtl::AutoDScalar::numdir);
        mt0->mapping(fvm_node->root_node(), n0_data, SolverSpecify::clock);

        fvm_node0 = fvm_node;
        V0 = x[fvm_node->local_offset()+0];  V0.setADValue(0,1.0);
        n0 = x[fvm_node->local_offset()+1];  n0.setADValue(1,1.0);  // electron density
        p0 = x[fvm_node->local_offset()+2];  p0.setADValue(2,1.0);  // hole density



        Ec0 =  -(e*V0 + n0_data->affinity() + mt0->band->EgNarrowToEc(p0, n0, T) + kb*T*log(n0_data->Nc()));
        Ev0 =  -(e*V0 + n0_data->affinity() - mt0->band->EgNarrowToEv(p0, n0, T) - kb*T*log(n0_data->Nv()) + mt0->band->Eg(T));
        if(semi_region->get_advanced_model()->Fermi)
        {
          Ec0 = Ec0 - e*Vt*log(gamma_f(fabs(n0)/n0_data->Nc()));
          Ev0 = Ev0 + e*Vt*log(gamma_f(fabs(p0)/n0_data->Nv()));
        }
      }

      // other semiconductor region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              std::vector<int> rows, cols;
              rows.push_back(fvm_node0->global_offset()+0);
              rows.push_back(fvm_node0->global_offset()+1);
              rows.push_back(fvm_node0->global_offset()+2);
              rows.push_back(fvm_node->global_offset()+0);
              rows.push_back(fvm_node->global_offset()+1);
              rows.push_back(fvm_node->global_offset()+2);
              cols = rows;

              const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
              const FVM_NodeData * n_data = fvm_node->node_data();

              AutoDScalar V = x[fvm_node->local_offset()+0];  V.setADValue(3,1.0); // psi of this node
              AutoDScalar n = x[fvm_node->local_offset()+1];  n.setADValue(4,1.0); // electron density
              AutoDScalar p = x[fvm_node->local_offset()+2];  p.setADValue(5,1.0); // hole density

              // mapping this node to material library
              Material::MaterialSemiconductor *mt = semi_region->material();
              mt->set_ad_num(adtl::AutoDScalar::numdir);
              mt->mapping(fvm_node->root_node(), n_data, SolverSpecify::clock);
              AutoDScalar Ec =  -(e*V + n_data->affinity() + mt->band->EgNarrowToEc(p, n, T) + kb*T*log(n_data->Nc()));
              AutoDScalar Ev =  -(e*V + n_data->affinity() - mt->band->EgNarrowToEv(p, n, T) - kb*T*log(n_data->Nv()) + mt->band->Eg(T));
              if(semi_region->get_advanced_model()->Fermi)
              {
                Ec = Ec - e*Vt*log(gamma_f(fabs(n)/n_data->Nc()));
                Ev = Ev + e*Vt*log(gamma_f(fabs(p)/n_data->Nv()));
              }

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              AutoDScalar ff1 = V - V0;
              jac->add_row(  rows[3],  cols.size(),  &cols[0],  ff1.getADValue() );

              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_node->outside_boundary_surface_area();

              // thermal emit current
              if(Ec0 > Ec)
              {
                // electrons are leaving region 0
                PetscScalar pm = mt0->band->EffecElecMass(T)/mt->band->EffecElecMass(T);
                AutoDScalar Jn = 2.*(mt0->band->ThermalVn(T)*n0 - pm*mt->band->ThermalVn(T)*(n*exp(-(Ec0-Ec)/(kb*T))))*cv_boundary;
                jac->add_row(  rows[4],  cols.size(),  &cols[0],  Jn.getADValue() );
                jac->add_row(  rows[1],  cols.size(),  &cols[0],  (-Jn).getADValue() );
              }
              else
              {
                // electrons are leaving this region
                PetscScalar pm = mt->band->EffecElecMass(T)/mt0->band->EffecElecMass(T);
                AutoDScalar Jn = 2.*(mt->band->ThermalVn(T)*n - pm*mt0->band->ThermalVn(T)*(n0*exp(-(Ec-Ec0)/(kb*T))))*cv_boundary;
                jac->add_row(  rows[4],  cols.size(),  &cols[0],  (-Jn).getADValue() );
                jac->add_row(  rows[1],  cols.size(),  &cols[0],  Jn.getADValue() );
              }

              if(Ev0 < Ev)
              {
                // holes are leaving region 0
                PetscScalar pm = mt0->band->EffecHoleMass(T)/mt->band->EffecHoleMass(T);
                AutoDScalar Jp = 2.*(mt0->band->ThermalVp(T)*p0 - pm*mt->band->ThermalVp(T)*(p*exp((Ev0-Ev)/(kb*T))))*cv_boundary;
                jac->add_row(  rows[5],  cols.size(),  &cols[0],  Jp.getADValue() );
                jac->add_row(  rows[2],  cols.size(),  &cols[0],  (-Jp).getADValue() );
              }
              else
              {
                // holes are leaving this region
                PetscScalar pm = mt->band->EffecHoleMass(T)/mt0->band->EffecHoleMass(T);
                AutoDScalar Jp = 2.*(mt->band->ThermalVp(T)*p - pm*mt0->band->ThermalVp(T)*(p0*exp((Ev-Ev0)/(kb*T))))*cv_boundary;
                jac->add_row(  rows[5],  cols.size(),  &cols[0],  (-Jp).getADValue() );
                jac->add_row(  rows[2],  cols.size(),  &cols[0],  Jp.getADValue() );
              }
              break;
            }
            case InsulatorRegion :
            {
              // the governing equation of this fvm node
              AutoDScalar V = x[fvm_node->local_offset()+0];  V.setADValue(3,1.0); // psi of this node

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              AutoDScalar ff1 = V - V0;
              jac->add( fvm_node->global_offset(),  fvm_node0->global_offset(),  ff1.getADValue(0) );
              jac->add( fvm_node->global_offset(),  fvm_node->global_offset(),  ff1.getADValue(3) );
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




