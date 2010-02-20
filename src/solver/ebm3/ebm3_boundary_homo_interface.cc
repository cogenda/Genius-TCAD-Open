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
#include "boundary_condition.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void HomoInterfaceBC::EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // buffer for Vec location
  std::vector<PetscInt> src_row;
  std::vector<PetscInt> dst_row;

  // buffer for Vec value
  std::vector<PetscInt>    iy;
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

      // the first semiconductor region
      if(i==0)
      {
        // do nothing.
        // however, we will add fvm integral of other regions to it.
      }

      // other semiconductor region
      else
      {
        unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
        unsigned int node_n_offset   = regions[i]->ebm_variable_offset(ELECTRON);
        unsigned int node_p_offset   = regions[i]->ebm_variable_offset(HOLE);
        unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);
        unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset(E_TEMP);
        unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset(H_TEMP);


        // record the source row and dst row
        src_row.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);
        dst_row.push_back(fvm_nodes[0]->global_offset()+node_psi_offset);

        src_row.push_back(fvm_nodes[i]->global_offset()+node_n_offset);
        dst_row.push_back(fvm_nodes[0]->global_offset()+node_n_offset);

        src_row.push_back(fvm_nodes[i]->global_offset()+node_p_offset);
        dst_row.push_back(fvm_nodes[0]->global_offset()+node_p_offset);

        // when have the lattice temperature equation, we add dst to src
        if(regions[i]->get_advanced_model()->enable_Tl())
        {
          // lattice temperature equation should be global
          genius_assert(regions[0]->get_advanced_model()->enable_Tl());
          src_row.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);
          dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
        }

        // when both 2 regions have the electron/hole temperature equation, we add dst to src
        if(regions[0]->get_advanced_model()->enable_Tn() && regions[i]->get_advanced_model()->enable_Tn())
        {
          src_row.push_back(fvm_nodes[i]->global_offset()+node_Tn_offset);
          dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP));
        }

        if(regions[0]->get_advanced_model()->enable_Tp() && regions[i]->get_advanced_model()->enable_Tp())
        {
          src_row.push_back(fvm_nodes[i]->global_offset()+node_Tp_offset);
          dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP));
        }


        // the governing equation of this fvm node

        // the solution value of this node is equal to corresponding node value in the first semiconductor region

        // psi of this node
        PetscScalar V = x[fvm_nodes[i]->local_offset()+node_psi_offset];
        PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+node_psi_offset];
        PetscScalar ff1 = V - V_semi;
        iy.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);
        y_new.push_back(ff1);

        // electron density
        PetscScalar n = x[fvm_nodes[i]->local_offset()+node_n_offset];
        PetscScalar n_semi = x[fvm_nodes[0]->local_offset()+node_n_offset];
        PetscScalar ff2 = n - n_semi;
        iy.push_back(fvm_nodes[i]->global_offset()+node_n_offset);
        y_new.push_back(ff2);

        // hole density
        PetscScalar p = x[fvm_nodes[i]->local_offset()+node_p_offset];
        PetscScalar p_semi = x[fvm_nodes[0]->local_offset()+node_p_offset];  // hole density
        PetscScalar ff3 = p - p_semi;
        iy.push_back(fvm_nodes[i]->global_offset()+node_p_offset);
        y_new.push_back(ff3);

        // variables for first semiconductor region
        PetscScalar T_semi  = T_external();
        PetscScalar Tn_semi = T_external();
        PetscScalar Tp_semi = T_external();

        if(regions[0]->get_advanced_model()->enable_Tl())
          T_semi =  x[fvm_nodes[0]->local_offset() + regions[0]->ebm_variable_offset(TEMPERATURE)];

        if(regions[0]->get_advanced_model()->enable_Tn())
          Tn_semi = x[fvm_nodes[0]->local_offset() + regions[0]->ebm_variable_offset(E_TEMP)]/n_semi;

        if(regions[0]->get_advanced_model()->enable_Tp())
          Tp_semi = x[fvm_nodes[0]->local_offset() + regions[0]->ebm_variable_offset(H_TEMP)]/p_semi;


        // extra equations for this region
        PetscScalar T  = T_external();
        PetscScalar Tn = T_external();
        PetscScalar Tp = T_external();

        if(regions[i]->get_advanced_model()->enable_Tl())
        {
          T =  x[fvm_nodes[i]->local_offset() + node_Tl_offset];
          PetscScalar ff4 = T - T_semi;
          iy.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);
          y_new.push_back(ff4);
        }

        if(regions[i]->get_advanced_model()->enable_Tn())
        {
          Tn = x[fvm_nodes[i]->local_offset() + node_Tn_offset]/n;
          PetscScalar ff5 = n*(Tn - Tn_semi);
          iy.push_back(fvm_nodes[i]->global_offset()+node_Tn_offset);
          y_new.push_back(ff5);
        }

        if(regions[i]->get_advanced_model()->enable_Tp())
        {
          Tp = x[fvm_nodes[i]->local_offset() + node_Tp_offset]/p;
          PetscScalar ff6 = p*(Tp - Tp_semi);
          iy.push_back(fvm_nodes[i]->global_offset()+node_Tp_offset);
          y_new.push_back(ff6);
        }

      }
    }

  }


  // add src row to dst row, it will assemble vec automatically
  PetscUtils::VecAddRowToRow(f, src_row, dst_row);

  // insert new value to src row
  if( src_row.size() )
    VecSetValues(f, src_row.size(), &(iy[0]), &(y_new[0]), INSERT_VALUES);

  add_value_flag = INSERT_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for EBM3 solver
 */
void HomoInterfaceBC::EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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

      // the first semiconductor region
      if(i==0)
      {
        // reserve items for all the ghost nodes
        std::vector<int> rows, cols;
        for(unsigned int nv=0; nv<regions[0]->ebm_n_variables(); ++nv)
          rows.push_back(fvm_nodes[0]->global_offset()+nv);

        FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
        for(; gn_it != fvm_nodes[i]->ghost_node_end(); ++gn_it)
        {
          const FVM_Node * ghost_fvm_node = (*gn_it).first;
          genius_assert(ghost_fvm_node!=NULL);

          for(unsigned int nv=0; nv<regions[0]->ebm_n_variables(); ++nv)
            cols.push_back(ghost_fvm_node->global_offset()+nv);

          FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
          for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
          {
            for(unsigned int nv=0; nv<regions[0]->ebm_n_variables(); ++nv)
              cols.push_back((*gnb_it).second->global_offset()+nv);
          }
        }

        std::vector<PetscScalar> value(rows.size()*cols.size(),0);

        MatSetValues(*jac, rows.size(), &rows[0], cols.size(), &cols[0], &value[0], ADD_VALUES);

      }

      // other semiconductor region
      else
      {
        // reserve for later operator
        for(unsigned int nv=0; nv<regions[i]->ebm_n_variables(); ++nv)
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+nv, fvm_nodes[0]->global_offset()+nv, 0, ADD_VALUES);
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}






/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void HomoInterfaceBC::EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
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

        // the first semiconductor region
        if(i==0) continue;

        // other semiconductor region
        else
        {
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_n_offset   = regions[i]->ebm_variable_offset(ELECTRON);
          unsigned int node_p_offset   = regions[i]->ebm_variable_offset(HOLE);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);
          unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset(E_TEMP);
          unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset(H_TEMP);


          // record the source row and dst row

          src_row.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);
          dst_row.push_back(fvm_nodes[0]->global_offset()+node_psi_offset);

          src_row.push_back(fvm_nodes[i]->global_offset()+node_n_offset);
          dst_row.push_back(fvm_nodes[0]->global_offset()+node_n_offset);

          src_row.push_back(fvm_nodes[i]->global_offset()+node_p_offset);
          dst_row.push_back(fvm_nodes[0]->global_offset()+node_p_offset);

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            src_row.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);
            dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
          }

          if(regions[0]->get_advanced_model()->enable_Tn() && regions[i]->get_advanced_model()->enable_Tn())
          {
            src_row.push_back(fvm_nodes[i]->global_offset()+node_Tn_offset);
            dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP));
          }

          if(regions[0]->get_advanced_model()->enable_Tp() && regions[i]->get_advanced_model()->enable_Tp())
          {
            src_row.push_back(fvm_nodes[i]->global_offset()+node_Tp_offset);
            dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP));
          }
        }
      }

    }

    //ok, we add source rows to destination rows
    PetscUtils::MatAddRowToRow(*jac, src_row, dst_row);

    // clear src rows
    MatZeroRows(*jac, src_row.size(), src_row.empty() ? NULL : &src_row[0], 0.0);

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
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      if(!fvm_node->is_valid()) continue;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0) continue;

      // other semiconductor region
      else
      {
        unsigned int global_offset   = fvm_nodes[i]->global_offset();
        unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
        unsigned int node_n_offset   = regions[i]->ebm_variable_offset(ELECTRON);
        unsigned int node_p_offset   = regions[i]->ebm_variable_offset(HOLE);
        unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);
        unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset(E_TEMP);
        unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset(H_TEMP);

        //the indepedent variable number, we need 4 here.
        adtl::AutoDScalar::numdir=4;

        // psi of this node
        AutoDScalar V = x[fvm_nodes[i]->local_offset()+node_psi_offset];       V.setADValue(0,1.0);
        AutoDScalar V_semi = x[fvm_nodes[0]->local_offset()+node_psi_offset];  V_semi.setADValue(1,1.0);
        AutoDScalar ff1 = V - V_semi;
        MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[i]->global_offset()+node_psi_offset, ff1.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[0]->global_offset()+node_psi_offset, ff1.getADValue(1), ADD_VALUES);

        // electron density
        AutoDScalar n = x[fvm_nodes[i]->local_offset()+node_n_offset];         n.setADValue(0,1.0);
        AutoDScalar n_semi = x[fvm_nodes[0]->local_offset()+node_n_offset];    n_semi.setADValue(1,1.0);
        AutoDScalar ff2 = n - n_semi;
        MatSetValue(*jac, global_offset+node_n_offset, fvm_nodes[i]->global_offset()+node_n_offset, ff2.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, global_offset+node_n_offset, fvm_nodes[0]->global_offset()+node_n_offset, ff2.getADValue(1), ADD_VALUES);


        // hole density
        AutoDScalar p = x[fvm_nodes[i]->local_offset()+node_p_offset];         p.setADValue(0,1.0);
        AutoDScalar p_semi = x[fvm_nodes[0]->local_offset()+node_p_offset];    p_semi.setADValue(1,1.0);  // hole density
        AutoDScalar ff3 = p - p_semi;
        MatSetValue(*jac, global_offset+node_p_offset, fvm_nodes[i]->global_offset()+node_p_offset, ff3.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, global_offset+node_p_offset, fvm_nodes[0]->global_offset()+node_p_offset, ff3.getADValue(1), ADD_VALUES);


        // variables for first semiconductor region
        AutoDScalar T_semi   = T_external();
        AutoDScalar Tn_semi  = T_external();
        AutoDScalar Tp_semi  = T_external();

        // extra equations for this region
        AutoDScalar T  = T_external();
        AutoDScalar Tn = T_external();
        AutoDScalar Tp = T_external();

        if(regions[i]->get_advanced_model()->enable_Tl())
        {
          T =  x[fvm_nodes[i]->local_offset() + node_Tl_offset];  T.setADValue(0,1.0);
          genius_assert(regions[0]->get_advanced_model()->enable_Tl());
          AutoDScalar T_semi = x[fvm_nodes[0]->local_offset() + regions[0]->ebm_variable_offset(TEMPERATURE)]; T_semi.setADValue(1,1.0);
          AutoDScalar ff4 = T - T_semi;
          MatSetValue(*jac, global_offset+node_Tl_offset, fvm_nodes[i]->global_offset()+node_Tl_offset, ff4.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, global_offset+node_Tl_offset, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE), ff4.getADValue(1), ADD_VALUES);
        }

        if(regions[i]->get_advanced_model()->enable_Tn())
        {
          AutoDScalar nTn = x[fvm_nodes[i]->local_offset() + node_Tn_offset]; nTn.setADValue(2,1.0);
          Tn = nTn/n;
          if(regions[0]->get_advanced_model()->enable_Tn())
          {
            AutoDScalar nTn_semi = x[fvm_nodes[0]->local_offset() + regions[0]->ebm_variable_offset(E_TEMP)]; nTn_semi.setADValue(3,1.0);
            Tn_semi = nTn_semi/n_semi;
          }
          AutoDScalar ff5 = n*(Tn - Tn_semi);

          MatSetValue(*jac, global_offset+node_Tn_offset, fvm_nodes[i]->global_offset()+node_n_offset,  ff5.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, global_offset+node_Tn_offset, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(ELECTRON), ff5.getADValue(1), ADD_VALUES);
          MatSetValue(*jac, global_offset+node_Tn_offset, fvm_nodes[i]->global_offset()+node_Tn_offset, ff5.getADValue(2), ADD_VALUES);
          MatSetValue(*jac, global_offset+node_Tn_offset, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP), ff5.getADValue(3), ADD_VALUES);

        }

        if(regions[i]->get_advanced_model()->enable_Tp())
        {
          AutoDScalar pTp = x[fvm_nodes[i]->local_offset() + node_Tp_offset]; pTp.setADValue(2,1.0);
          Tp = pTp/p;
          if(regions[0]->get_advanced_model()->enable_Tp())
          {
            AutoDScalar pTp_semi = x[fvm_nodes[0]->local_offset() + regions[0]->ebm_variable_offset(H_TEMP)]; pTp_semi.setADValue(3,1.0);
            Tp_semi = pTp_semi/p_semi;
          }
          AutoDScalar ff6 = p*(Tp - Tp_semi);

          MatSetValue(*jac, global_offset+node_Tp_offset, fvm_nodes[i]->global_offset()+node_p_offset,  ff6.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, global_offset+node_Tp_offset, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(HOLE), ff6.getADValue(1), ADD_VALUES);
          MatSetValue(*jac, global_offset+node_Tp_offset, fvm_nodes[i]->global_offset()+node_Tp_offset, ff6.getADValue(2), ADD_VALUES);
          MatSetValue(*jac, global_offset+node_Tp_offset, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP), ff6.getADValue(3), ADD_VALUES);

        }

      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
