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
#include "boundary_condition.h"
#include "petsc_utils.h"
#include "mathfunc.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

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

  PetscScalar T = T_external();
  PetscScalar Vt  = kb*T/e;

  // buffer for Vec location
  std::vector<PetscInt> src_row;
  std::vector<PetscInt> dst_row;

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


    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

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

      // the variable for first region
      PetscScalar V0;
      PetscScalar n0;
      PetscScalar p0;
      Material::MaterialSemiconductor * mt0;
      PetscScalar Ec0,Ev0;

      // the first semiconductor region
      if(i==0)
      {
        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

        const FVM_NodeData * n0_data = fvm_nodes[i]->node_data();

        V0 = x[fvm_nodes[i]->local_offset()+0];
        n0 = x[fvm_nodes[i]->local_offset()+1];  // electron density
        p0 = x[fvm_nodes[i]->local_offset()+2];  // hole density

        mt0 = semi_region->material();
        mt0->mapping(fvm_nodes[i]->root_node(), n0_data, SolverSpecify::clock);

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
        // the ghost node should have the same processor_id with me
        genius_assert(fvm_nodes[i]->root_node()->processor_id() == fvm_nodes[0]->root_node()->processor_id() );

        // record the source row and dst row
        src_row.push_back(fvm_nodes[i]->global_offset()+0);
        dst_row.push_back(fvm_nodes[0]->global_offset()+0);


        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);
        const FVM_NodeData * n_data = fvm_nodes[i]->node_data();

        PetscScalar V = x[fvm_nodes[i]->local_offset()+0];  // psi of this node
        PetscScalar n = x[fvm_nodes[i]->local_offset()+1];  // electron density
        PetscScalar p = x[fvm_nodes[i]->local_offset()+2];  // hole density

        // mapping this node to material library
        Material::MaterialSemiconductor *mt = semi_region->material();
        mt->mapping(fvm_nodes[i]->root_node(), n_data, SolverSpecify::clock);
        PetscScalar Ec =  -(e*V + n_data->affinity() + mt->band->EgNarrowToEc(p, n, T) + kb*T*log(n_data->Nc()));
        PetscScalar Ev =  -(e*V + n_data->affinity() - mt->band->EgNarrowToEv(p, n, T) - kb*T*log(n_data->Nv()) + mt->band->Eg(T));
        if(semi_region->get_advanced_model()->Fermi)
        {
          Ec = Ec - e*Vt*log(gamma_f(fabs(n)/n_data->Nc()));
          Ev = Ev + e*Vt*log(gamma_f(fabs(p)/n_data->Nv()));
        }


        // the solution value of this node is equal to corresponding node value in the first semiconductor region
        PetscScalar ff1 = V - V0;
        iy.push_back(fvm_nodes[i]->global_offset()+0);
        y_new.push_back(ff1);


        // thermal emit current
        PetscScalar Jn=0,Jp=0;

        if(Ec0 > Ec)
        {
          PetscScalar pm = mt0->band->EffecElecMass(T)/mt->band->EffecElecMass(T);
          Jn = -2*e*(mt0->band->ThermalVn(T)*n0 - pm*mt->band->ThermalVn(T)*n*exp(-(Ec0-Ec)/(kb*T)));
        }
        else
        {
          PetscScalar pm = mt->band->EffecElecMass(T)/mt0->band->EffecElecMass(T);
          Jn = -2*e*(mt->band->ThermalVn(T)*n - pm*mt0->band->ThermalVn(T)*n0*exp(-(Ec-Ec0)/(kb*T)));
        }

        if(Ev0 > Ev)
        {
          PetscScalar pm = mt0->band->EffecHoleMass(T)/mt->band->EffecHoleMass(T);
          Jp = 2*e*(mt0->band->ThermalVp(T)*p0 - pm*mt->band->ThermalVp(T)*p*exp(-(Ev0-Ev)/(kb*T)));
        }
        else
        {
          PetscScalar pm = mt->band->EffecHoleMass(T)/mt0->band->EffecHoleMass(T);
          Jp = 2*e*(mt->band->ThermalVp(T)*p - pm*mt0->band->ThermalVp(T)*p0*exp(-(Ev-Ev0)/(kb*T)));
        }

        // area of out surface of control volume related with neighbor node
        PetscScalar cv_boundary = fvm_nodes[i]->outside_boundary_surface_area();

        // current flow
        VecSetValue(f, fvm_nodes[i]->global_offset()+1,  Jn*cv_boundary, ADD_VALUES);
        VecSetValue(f, fvm_nodes[0]->global_offset()+1, -Jn*cv_boundary, ADD_VALUES);

        VecSetValue(f, fvm_nodes[i]->global_offset()+2, -Jp*cv_boundary, ADD_VALUES);
        VecSetValue(f, fvm_nodes[0]->global_offset()+2,  Jp*cv_boundary, ADD_VALUES);

      }
    }

  }


  // add src row to dst row, it will assemble vec automatically
  PetscUtils::VecAddRowToRow(f, src_row, dst_row);

  // insert new value
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y_new[0]), INSERT_VALUES);

  add_value_flag = INSERT_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML1 solver
 */
void HeteroInterfaceBC::DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

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

      {
        // reserve items for all the ghost nodes
        std::vector<int> rows, cols;
        rows.push_back(fvm_nodes[i]->global_offset()+0);
        rows.push_back(fvm_nodes[i]->global_offset()+1);
        rows.push_back(fvm_nodes[i]->global_offset()+2);

        FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
        for(; gn_it != fvm_nodes[i]->ghost_node_end(); ++gn_it)
        {
          const FVM_Node * ghost_fvm_node = (*gn_it).first;
          genius_assert(ghost_fvm_node!=NULL);

          cols.push_back(ghost_fvm_node->global_offset()+0);
          cols.push_back(ghost_fvm_node->global_offset()+1);
          cols.push_back(ghost_fvm_node->global_offset()+2);

          FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
          for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
          {
            cols.push_back((*gnb_it).second->global_offset()+0);
            cols.push_back((*gnb_it).second->global_offset()+1);
            cols.push_back((*gnb_it).second->global_offset()+2);
          }
        }

        std::vector<PetscScalar> value(rows.size()*cols.size(),0);

        MatSetValues(*jac, rows.size(), &rows[0], cols.size(), &cols[0], &value[0], ADD_VALUES);

      }


    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}






/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void HeteroInterfaceBC::DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // here we do several things:
  // add some row to other, clear some row, insert some value to row
  // I wonder if there are some more efficient way to do these.

  {

    // buffer for mat rows which should be added to other row
    std::vector<PetscInt> src_row;
    std::vector<PetscInt> dst_row;

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
        const FVM_Node * fvm_node = (*rnode_it).second.second;
        if(!fvm_node->is_valid()) continue;

        fvm_nodes.push_back( fvm_node );

        // the first semiconductor region
        if(i==0) continue;

        // other semiconductor region
        else
        {
          src_row.push_back(fvm_nodes[i]->global_offset()+0);
          dst_row.push_back(fvm_nodes[0]->global_offset()+0);
        }
      }

    }

    //ok, we add source rows to destination rows
    PetscUtils::MatAddRowToRow(*jac, src_row, dst_row);

    // clear src_row
    MatZeroRows(*jac, src_row.size(), src_row.empty() ? NULL : &src_row[0], 0.0);

  }

  PetscScalar T = T_external();
  PetscScalar Vt  = kb*T/e;

  // after that, set values to source rows
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    //the indepedent variable number, we need 6 here.
    adtl::AutoDScalar::numdir=6;

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

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

      // the variable for first region
      AutoDScalar V0;
      AutoDScalar n0;
      AutoDScalar p0;
      Material::MaterialSemiconductor * mt0;
      AutoDScalar Ec0,Ev0;

      // the first semiconductor region
      if(i==0)
      {
        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

        const FVM_NodeData * n0_data = fvm_nodes[i]->node_data();

        V0 = x[fvm_nodes[i]->local_offset()+0];  V0.setADValue(0,1.0);
        n0 = x[fvm_nodes[i]->local_offset()+1];  n0.setADValue(1,1.0);  // electron density
        p0 = x[fvm_nodes[i]->local_offset()+2];  p0.setADValue(2,1.0);  // hole density

        mt0 = semi_region->material();
        mt0->mapping(fvm_nodes[i]->root_node(), n0_data, SolverSpecify::clock);

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
        // the ghost node should have the same processor_id with me
        genius_assert(fvm_nodes[i]->root_node()->processor_id() == fvm_nodes[0]->root_node()->processor_id() );

        std::vector<int> rows, cols;
        rows.push_back(fvm_nodes[0]->global_offset()+0);
        rows.push_back(fvm_nodes[0]->global_offset()+1);
        rows.push_back(fvm_nodes[0]->global_offset()+2);
        rows.push_back(fvm_nodes[i]->global_offset()+0);
        rows.push_back(fvm_nodes[i]->global_offset()+1);
        rows.push_back(fvm_nodes[i]->global_offset()+2);
        cols = rows;

        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);
        const FVM_NodeData * n_data = fvm_nodes[i]->node_data();

        AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];  V.setADValue(3,1.0); // psi of this node
        AutoDScalar n = x[fvm_nodes[i]->local_offset()+1];  n.setADValue(4,1.0); // electron density
        AutoDScalar p = x[fvm_nodes[i]->local_offset()+2];  p.setADValue(5,1.0); // hole density

        // mapping this node to material library
        Material::MaterialSemiconductor *mt = semi_region->material();
        mt->mapping(fvm_nodes[i]->root_node(), n_data, SolverSpecify::clock);
        AutoDScalar Ec =  -(e*V + n_data->affinity() + mt->band->EgNarrowToEc(p, n, T) + kb*T*log(n_data->Nc()));
        AutoDScalar Ev =  -(e*V + n_data->affinity() - mt->band->EgNarrowToEv(p, n, T) - kb*T*log(n_data->Nv()) + mt->band->Eg(T));
        if(semi_region->get_advanced_model()->Fermi)
        {
          Ec = Ec - e*Vt*log(gamma_f(fabs(n)/n_data->Nc()));
          Ev = Ev + e*Vt*log(gamma_f(fabs(p)/n_data->Nv()));
        }

        // the solution value of this node is equal to corresponding node value in the first semiconductor region
        AutoDScalar ff1 = V - V0;
        MatSetValues(*jac, 1, &rows[3], cols.size(), &cols[0], ff1.getADValue(), ADD_VALUES);

        // thermal emit current
        AutoDScalar Jn=0,Jp=0;

        if(Ec0 > Ec)
        {
          PetscScalar pm = mt0->band->EffecElecMass(T)/mt->band->EffecElecMass(T);
          Jn = -2*e*(mt0->band->ThermalVn(T)*n0 - pm*mt->band->ThermalVn(T)*n*exp(-(Ec0-Ec)/(kb*T)));
        }
        else
        {
          PetscScalar pm = mt->band->EffecElecMass(T)/mt0->band->EffecElecMass(T);
          Jn = -2*e*(mt->band->ThermalVn(T)*n - pm*mt0->band->ThermalVn(T)*n0*exp(-(Ec-Ec0)/(kb*T)));
        }

        if(Ev0 > Ev)
        {
          PetscScalar pm = mt0->band->EffecHoleMass(T)/mt->band->EffecHoleMass(T);
          Jp = 2*e*(mt0->band->ThermalVp(T)*p0 - pm*mt->band->ThermalVp(T)*p*exp(-(Ev0-Ev)/(kb*T)));
        }
        else
        {
          PetscScalar pm = mt->band->EffecHoleMass(T)/mt0->band->EffecHoleMass(T);
          Jp = 2*e*(mt->band->ThermalVp(T)*p - pm*mt0->band->ThermalVp(T)*p0*exp(-(Ev-Ev0)/(kb*T)));
        }

        // area of out surface of control volume related with neighbor node
        PetscScalar cv_boundary = fvm_nodes[i]->outside_boundary_surface_area();

        // current flow
        AutoDScalar Fn = Jn*cv_boundary;
        MatSetValues(*jac, 1, &rows[4], cols.size(), &cols[0], Fn.getADValue(), ADD_VALUES);
        MatSetValues(*jac, 1, &rows[1], cols.size(), &cols[0], (-Fn).getADValue(), ADD_VALUES);

        AutoDScalar Fp = -Jp*cv_boundary;
        MatSetValues(*jac, 1, &rows[5], cols.size(), &cols[0], Fp.getADValue(), ADD_VALUES);
        MatSetValues(*jac, 1, &rows[2], cols.size(), &cols[0], (-Fp).getADValue(), ADD_VALUES);

      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
