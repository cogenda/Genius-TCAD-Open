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
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_homo.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * do pre-process to function for DDML2 solver
 */
void HomoInterfaceBC::DDM2_Function_Preprocess(PetscScalar * ,Vec f, std::vector<PetscInt> &src_row,
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

      // other regions
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              src_row.push_back(fvm_nodes[i]->global_offset()+1);
              src_row.push_back(fvm_nodes[i]->global_offset()+2);
              src_row.push_back(fvm_nodes[i]->global_offset()+3);

              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+1);
              dst_row.push_back(fvm_nodes[0]->global_offset()+2);
              dst_row.push_back(fvm_nodes[0]->global_offset()+3);

              clear_row.push_back(fvm_nodes[i]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+1);
              clear_row.push_back(fvm_nodes[i]->global_offset()+2);
              clear_row.push_back(fvm_nodes[i]->global_offset()+3);
              break;
            }
            case InsulatorRegion:
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              src_row.push_back(fvm_nodes[i]->global_offset()+1);
              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+3);
              clear_row.push_back(fvm_nodes[i]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+1);
              break;
            }
            default: genius_error();
        }
      }

    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void HomoInterfaceBC::DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

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

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
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
      // other regions
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // the governing equation of this fvm node

              PetscScalar V = x[fvm_nodes[i]->local_offset()+0];  // psi of this node
              PetscScalar n = x[fvm_nodes[i]->local_offset()+1];  // electron density
              PetscScalar p = x[fvm_nodes[i]->local_offset()+2];  // hole density
              PetscScalar T = x[fvm_nodes[i]->local_offset()+3];  // lattice temperature

              PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+0];
              PetscScalar n_semi = x[fvm_nodes[0]->local_offset()+1];  // electron density
              PetscScalar p_semi = x[fvm_nodes[0]->local_offset()+2];  // hole density
              PetscScalar T_semi = x[fvm_nodes[0]->local_offset()+3];  // lattice temperature

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              PetscScalar ff1 = V - V_semi;
              PetscScalar ff2 = n - n_semi;
              PetscScalar ff3 = p - p_semi;
              PetscScalar ff4 = T - T_semi;

              iy.push_back(fvm_nodes[i]->global_offset()+0);
              iy.push_back(fvm_nodes[i]->global_offset()+1);
              iy.push_back(fvm_nodes[i]->global_offset()+2);
              iy.push_back(fvm_nodes[i]->global_offset()+3);
              y_new.push_back(ff1);
              y_new.push_back(ff2);
              y_new.push_back(ff3);
              y_new.push_back(ff4);

              genius_assert(iy.size()==y_new.size());
              break;
            }

            case InsulatorRegion:
            {
              // the governing equation of this fvm node
              PetscScalar V = x[fvm_nodes[i]->local_offset()+0];  // psi of this node
              PetscScalar T = x[fvm_nodes[i]->local_offset()+1];  // lattice temperature

              PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+0];
              PetscScalar T_semi = x[fvm_nodes[0]->local_offset()+3];  // lattice temperature

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              PetscScalar ff1 = V - V_semi;
              PetscScalar ff2 = T - T_semi;

              iy.push_back(fvm_nodes[i]->global_offset()+0);
              y_new.push_back(ff1);
              iy.push_back(fvm_nodes[i]->global_offset()+1);
              y_new.push_back(ff2);

              genius_assert(iy.size()==y_new.size());
              break;
            }
            default: genius_error();
        }

      }
    }
  }

  // set new value to row
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y_new[0]), ADD_VALUES);

  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML2 solver
 */
void HomoInterfaceBC::DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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
      if(i==0)
      {
        genius_assert( region->type() == SemiconductorRegion );
        // do nothing.
      }

      // other semiconductor region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // reserve items for all the ghost nodes
              std::vector<int> rows, cols;
              rows.push_back(fvm_nodes[0]->global_offset()+0);
              rows.push_back(fvm_nodes[0]->global_offset()+1);
              rows.push_back(fvm_nodes[0]->global_offset()+2);
              rows.push_back(fvm_nodes[0]->global_offset()+3);

              cols.push_back(fvm_nodes[i]->global_offset()+0);
              cols.push_back(fvm_nodes[i]->global_offset()+1);
              cols.push_back(fvm_nodes[i]->global_offset()+2);
              cols.push_back(fvm_nodes[i]->global_offset()+3);

              FVM_Node::fvm_neighbor_node_iterator  nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                cols.push_back((*nb_it).first->global_offset()+0);
                cols.push_back((*nb_it).first->global_offset()+1);
                cols.push_back((*nb_it).first->global_offset()+2);
                cols.push_back((*nb_it).first->global_offset()+3);
              }

              std::vector<PetscScalar> value(rows.size()*cols.size(),0);

              MatSetValues(*jac, rows.size(), &rows[0], cols.size(), &cols[0], &value[0], ADD_VALUES);

              // reserve for later operator
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, fvm_nodes[0]->global_offset()+0, 0, ADD_VALUES);
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[0]->global_offset()+1, 0, ADD_VALUES);
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+2, fvm_nodes[0]->global_offset()+2, 0, ADD_VALUES);
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+3, fvm_nodes[0]->global_offset()+3, 0, ADD_VALUES);
              break;
            }
            case InsulatorRegion:
            {
              // reserve items for all the ghost nodes
              std::vector<int> rows, cols;
              rows.push_back(fvm_nodes[0]->global_offset()+0);
              rows.push_back(fvm_nodes[0]->global_offset()+3);
              cols.push_back(fvm_nodes[i]->global_offset()+0);
              cols.push_back(fvm_nodes[i]->global_offset()+1);

              FVM_Node::fvm_neighbor_node_iterator  nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                cols.push_back((*nb_it).first->global_offset()+0);
                cols.push_back((*nb_it).first->global_offset()+1);
              }

              std::vector<PetscScalar> value(rows.size()*cols.size(),0);

              MatSetValues(*jac, rows.size(), &rows[0], cols.size(), &cols[0], &value[0], ADD_VALUES);

              MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, fvm_nodes[0]->global_offset()+0, 0, ADD_VALUES);
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[0]->global_offset()+3, 0, ADD_VALUES);
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



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML2 solver
 */
void HomoInterfaceBC::DDM2_Jacobian_Preprocess(PetscScalar *,Mat *jac, std::vector<PetscInt> &src_row,
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

      // other semiconductor region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              src_row.push_back(fvm_nodes[i]->global_offset()+1);
              src_row.push_back(fvm_nodes[i]->global_offset()+2);
              src_row.push_back(fvm_nodes[i]->global_offset()+3);

              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+1);
              dst_row.push_back(fvm_nodes[0]->global_offset()+2);
              dst_row.push_back(fvm_nodes[0]->global_offset()+3);

              clear_row.push_back(fvm_nodes[i]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+1);
              clear_row.push_back(fvm_nodes[i]->global_offset()+2);
              clear_row.push_back(fvm_nodes[i]->global_offset()+3);
              break;
            }
            case InsulatorRegion:
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              src_row.push_back(fvm_nodes[i]->global_offset()+1);
              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+3);
              clear_row.push_back(fvm_nodes[i]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+1);
              break;
            }
            default: genius_error();
        }
      }

    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void HomoInterfaceBC::DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of HomoInterface boundary condition is processed here

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  // after that, set values to source rows
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
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              //the indepedent variable number, we need 8 here.
              adtl::AutoDScalar::numdir=8;

              std::vector<int> rows, cols;

              // the governing equation of this fvm node
              AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];  V.setADValue(0,1.0);  // psi of this node
              AutoDScalar n = x[fvm_nodes[i]->local_offset()+1];  n.setADValue(1,1.0);  // electron density
              AutoDScalar p = x[fvm_nodes[i]->local_offset()+2];  p.setADValue(2,1.0);  // hole density
              AutoDScalar T = x[fvm_nodes[i]->local_offset()+3];  T.setADValue(3,1.0);  // lattice temperature

              AutoDScalar V_semi = x[fvm_nodes[0]->local_offset()+0]; V_semi.setADValue(4,1.0);
              AutoDScalar n_semi = x[fvm_nodes[0]->local_offset()+1]; n_semi.setADValue(5,1.0);  // electron density
              AutoDScalar p_semi = x[fvm_nodes[0]->local_offset()+2]; p_semi.setADValue(6,1.0);  // hole density
              AutoDScalar T_semi = x[fvm_nodes[0]->local_offset()+3]; T_semi.setADValue(7,1.0);  // lattice temperature

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              AutoDScalar ff1 = V - V_semi;
              AutoDScalar ff2 = n - n_semi;
              AutoDScalar ff3 = p - p_semi;
              AutoDScalar ff4 = T - T_semi;

              rows.push_back(fvm_nodes[i]->global_offset()+0);
              rows.push_back(fvm_nodes[i]->global_offset()+1);
              rows.push_back(fvm_nodes[i]->global_offset()+2);
              rows.push_back(fvm_nodes[i]->global_offset()+3);

              cols = rows;
              cols.push_back(fvm_nodes[0]->global_offset()+0);
              cols.push_back(fvm_nodes[0]->global_offset()+1);
              cols.push_back(fvm_nodes[0]->global_offset()+2);
              cols.push_back(fvm_nodes[0]->global_offset()+3);

              // set Jacobian of governing equations
              MatSetValues(*jac, 1, &rows[0], cols.size(), &cols[0], ff1.getADValue(), ADD_VALUES);
              MatSetValues(*jac, 1, &rows[1], cols.size(), &cols[0], ff2.getADValue(), ADD_VALUES);
              MatSetValues(*jac, 1, &rows[2], cols.size(), &cols[0], ff3.getADValue(), ADD_VALUES);
              MatSetValues(*jac, 1, &rows[3], cols.size(), &cols[0], ff4.getADValue(), ADD_VALUES);
              break;
            }
            case InsulatorRegion :
            {
              //the indepedent variable number, we need 2 here.
              adtl::AutoDScalar::numdir=2;

              // the governing equation of this fvm node
              AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];  V.setADValue(0,1.0); // psi of this node
              AutoDScalar T = x[fvm_nodes[i]->local_offset()+1];  T.setADValue(0,1.0);// lattice temperature
              AutoDScalar V_semi = x[fvm_nodes[0]->local_offset()+0]; V_semi.setADValue(1,1.0);
              AutoDScalar T_semi = x[fvm_nodes[0]->local_offset()+3]; T_semi.setADValue(1,1.0); // lattice temperature
              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              AutoDScalar ff1 = V - V_semi;
              AutoDScalar ff2 = T - T_semi;

              PetscInt row_psi = fvm_nodes[i]->global_offset()+0;
              PetscInt cols_psi[2] = { fvm_nodes[i]->global_offset()+0,  fvm_nodes[0]->global_offset()+0};
              MatSetValues(*jac, 1, &row_psi, 2, cols_psi, ff1.getADValue(), ADD_VALUES);

              PetscInt row_t = fvm_nodes[i]->global_offset()+1;
              PetscInt cols_t[2] = { fvm_nodes[i]->global_offset()+1,  fvm_nodes[0]->global_offset()+3};
              MatSetValues(*jac, 1, &row_t, 2, cols_t, ff2.getADValue(), ADD_VALUES);

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
