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
#include "boundary_condition_homo.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for EBM3 solver
 */
void HomoInterfaceBC::EBM3_Function_Preprocess(PetscScalar *x, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  _ddm_current_interface(x, f);

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

      // other  region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(POTENTIAL));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(POTENTIAL));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(POTENTIAL));

              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(ELECTRON));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(ELECTRON));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(ELECTRON));

              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(HOLE));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(HOLE));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(HOLE));

              // when have the lattice temperature equation, we add dst to src
              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                // lattice temperature equation should be global
                genius_assert(regions[0]->get_advanced_model()->enable_Tl());
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }

              // when both 2 regions have the electron/hole temperature equation, we add dst to src
              if(regions[i]->get_advanced_model()->enable_Tn())
              {
                if(regions[0]->get_advanced_model()->enable_Tn())
                {
                  src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP));
                  dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP));
                }
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP));
              }

              if(regions[i]->get_advanced_model()->enable_Tp())
              {
                if(regions[0]->get_advanced_model()->enable_Tp())
                {
                  src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP));
                  dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP));
                }
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP));
              }

              break;
            }

            case InsulatorRegion:
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+0);

              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }
              break;
            }
            default: genius_error();
        }
      }
    }

  }

}



/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void HomoInterfaceBC::EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

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

      // other region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
              unsigned int node_n_offset   = regions[i]->ebm_variable_offset(ELECTRON);
              unsigned int node_p_offset   = regions[i]->ebm_variable_offset(HOLE);
              unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);
              unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset(E_TEMP);
              unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset(H_TEMP);


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
              break;
            }

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
            default: genius_error();
        }
      }
    }

  }

  // insert new value to src row
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y_new[0]), ADD_VALUES);

  add_value_flag = ADD_VALUES;
}






/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for EBM3 solver
 */
void HomoInterfaceBC::EBM3_Jacobian_Preprocess(PetscScalar * ,SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
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

      // other  region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(POTENTIAL));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(POTENTIAL));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(POTENTIAL));

              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(ELECTRON));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(ELECTRON));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(ELECTRON));

              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(HOLE));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(HOLE));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(HOLE));

              // when have the lattice temperature equation, we add dst to src
              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                // lattice temperature equation should be global
                genius_assert(regions[0]->get_advanced_model()->enable_Tl());
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }

              // when both 2 regions have the electron/hole temperature equation, we add dst to src
              if(regions[i]->get_advanced_model()->enable_Tn())
              {
                if(regions[0]->get_advanced_model()->enable_Tn())
                {
                  src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP));
                  dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP));
                }
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP));
              }

              if(regions[i]->get_advanced_model()->enable_Tp())
              {
                if(regions[0]->get_advanced_model()->enable_Tp())
                {
                  src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP));
                  dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP));
                }
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP));
              }

              break;
            }

            case InsulatorRegion:
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+0);

              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }
              break;
            }
            default: genius_error();
        }
      }
    }

  }
}




/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void HomoInterfaceBC::EBM3_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{

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

      // other region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
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
              jac->add( global_offset+node_psi_offset,  fvm_nodes[i]->global_offset()+node_psi_offset,  ff1.getADValue(0) );
              jac->add( global_offset+node_psi_offset,  fvm_nodes[0]->global_offset()+node_psi_offset,  ff1.getADValue(1) );

              // electron density
              AutoDScalar n = x[fvm_nodes[i]->local_offset()+node_n_offset];         n.setADValue(0,1.0);
              AutoDScalar n_semi = x[fvm_nodes[0]->local_offset()+node_n_offset];    n_semi.setADValue(1,1.0);
              AutoDScalar ff2 = n - n_semi;
              jac->add( global_offset+node_n_offset,  fvm_nodes[i]->global_offset()+node_n_offset,  ff2.getADValue(0) );
              jac->add( global_offset+node_n_offset,  fvm_nodes[0]->global_offset()+node_n_offset,  ff2.getADValue(1) );


              // hole density
              AutoDScalar p = x[fvm_nodes[i]->local_offset()+node_p_offset];         p.setADValue(0,1.0);
              AutoDScalar p_semi = x[fvm_nodes[0]->local_offset()+node_p_offset];    p_semi.setADValue(1,1.0);  // hole density
              AutoDScalar ff3 = p - p_semi;
              jac->add( global_offset+node_p_offset,  fvm_nodes[i]->global_offset()+node_p_offset,  ff3.getADValue(0) );
              jac->add( global_offset+node_p_offset,  fvm_nodes[0]->global_offset()+node_p_offset,  ff3.getADValue(1) );


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
                jac->add( global_offset+node_Tl_offset,  fvm_nodes[i]->global_offset()+node_Tl_offset,  ff4.getADValue(0) );
                jac->add( global_offset+node_Tl_offset,  fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE),  ff4.getADValue(1) );
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

                jac->add( global_offset+node_Tn_offset,  fvm_nodes[i]->global_offset()+node_n_offset,   ff5.getADValue(0) );
                jac->add( global_offset+node_Tn_offset,  fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(ELECTRON),  ff5.getADValue(1) );
                jac->add( global_offset+node_Tn_offset,  fvm_nodes[i]->global_offset()+node_Tn_offset,  ff5.getADValue(2) );
                jac->add( global_offset+node_Tn_offset,  fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP),  ff5.getADValue(3) );

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

                jac->add( global_offset+node_Tp_offset,  fvm_nodes[i]->global_offset()+node_p_offset,   ff6.getADValue(0) );
                jac->add( global_offset+node_Tp_offset,  fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(HOLE),  ff6.getADValue(1) );
                jac->add( global_offset+node_Tp_offset,  fvm_nodes[i]->global_offset()+node_Tp_offset,  ff6.getADValue(2) );
                jac->add( global_offset+node_Tp_offset,  fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP),  ff6.getADValue(3) );

              }

              break;
            }

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
              jac->add( global_offset+node_psi_offset,  fvm_nodes[i]->global_offset()+node_psi_offset,  ff1.getADValue(0) );
              jac->add( global_offset+node_psi_offset,  fvm_nodes[0]->global_offset()+semiconductor_node_psi_offset,  ff1.getADValue(1) );

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
                jac->add( global_offset+node_Tl_offset,  fvm_nodes[i]->global_offset()+node_Tl_offset,  ff2.getADValue(0) );
                jac->add( global_offset+node_Tl_offset,  fvm_nodes[0]->global_offset()+semiconductor_node_Tl_offset,  ff2.getADValue(1) );
              }

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

