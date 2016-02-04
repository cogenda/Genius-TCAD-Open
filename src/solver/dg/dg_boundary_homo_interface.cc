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
 * do pre-process to function for Density Gradient solver
 */
void HomoInterfaceBC::DG_Function_Preprocess(PetscScalar *x, Vec f, std::vector<PetscInt> &src_row,
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

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const FVM_Node *> fvm_nodes;
    std::vector<const SimulationRegion *> regions;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      fvm_nodes.push_back(fvm_node);
      regions.push_back(region);

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

            dst_row.push_back(fvm_nodes[0]->global_offset()+0);
            dst_row.push_back(fvm_nodes[0]->global_offset()+1);
            dst_row.push_back(fvm_nodes[0]->global_offset()+2);

            clear_row.push_back(fvm_nodes[i]->global_offset()+0);
            clear_row.push_back(fvm_nodes[i]->global_offset()+1);
            clear_row.push_back(fvm_nodes[i]->global_offset()+2);

            if(regions[0]->get_advanced_model()->QNEnabled && regions[i]->get_advanced_model()->QNEnabled)
            {
              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQC));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->dg_variable_offset(EQC));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQC));
            }

            if(regions[0]->get_advanced_model()->QPEnabled && regions[i]->get_advanced_model()->QPEnabled)
            {
              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQV));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->dg_variable_offset(EQV));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQV));
            }

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
 * build function and its jacobian for Density Gradient solver
 */
void HomoInterfaceBC::DG_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0)
      {
        genius_assert( regions[0]->type() == SemiconductorRegion );
        // do nothing.
        // however, we will add fvm integral of other regions to it.
      }

      // other semiconductor region
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

              PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+0];
              PetscScalar n_semi = x[fvm_nodes[0]->local_offset()+1];  // electron density
              PetscScalar p_semi = x[fvm_nodes[0]->local_offset()+2];  // hole density

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              PetscScalar ff1 = V - V_semi;
              PetscScalar ff2 = n - n_semi;
              PetscScalar ff3 = p - p_semi;

              iy.push_back(fvm_nodes[i]->global_offset()+0);
              iy.push_back(fvm_nodes[i]->global_offset()+1);
              iy.push_back(fvm_nodes[i]->global_offset()+2);
              y_new.push_back(ff1);
              y_new.push_back(ff2);
              y_new.push_back(ff3);

              if(regions[0]->get_advanced_model()->QNEnabled && regions[i]->get_advanced_model()->QNEnabled)
              {
                PetscScalar Eqc = x[fvm_nodes[i]->local_offset() + regions[i]->dg_variable_offset(EQC)];
                PetscScalar Eqc_semi = x[fvm_nodes[0]->local_offset() + regions[0]->dg_variable_offset(EQC)];
                PetscScalar ff = (Eqc - Eqc_semi);
                iy.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQC));
                y_new.push_back(ff);
              }

              if(regions[0]->get_advanced_model()->QPEnabled && regions[i]->get_advanced_model()->QPEnabled)
              {
                PetscScalar Eqv = x[fvm_nodes[i]->local_offset() + regions[i]->dg_variable_offset(EQV)];
                PetscScalar Eqv_semi = x[fvm_nodes[0]->local_offset() + regions[0]->dg_variable_offset(EQV)];
                PetscScalar ff = (Eqv - Eqv_semi);
                iy.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQV));
                y_new.push_back(ff);
              }

              break;
            }
            case InsulatorRegion:
            {
              // the governing equation of this fvm node
              PetscScalar V = x[fvm_nodes[i]->local_offset()+0];  // psi of this node
              PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+0];

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              PetscScalar ff1 = V - V_semi;

              iy.push_back(fvm_nodes[i]->global_offset()+0);
              y_new.push_back(ff1);
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
 * do pre-process to jacobian matrix for Density Gradient solver
 */
void HomoInterfaceBC::DG_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
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
    std::vector<const SimulationRegion *> regions;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      fvm_nodes.push_back(fvm_node);
      regions.push_back(region);

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

              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+1);
              dst_row.push_back(fvm_nodes[0]->global_offset()+2);

              clear_row.push_back(fvm_nodes[i]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+1);
              clear_row.push_back(fvm_nodes[i]->global_offset()+2);

              if(regions[0]->get_advanced_model()->QNEnabled && regions[i]->get_advanced_model()->QNEnabled)
              {
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQC));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->dg_variable_offset(EQC));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQC));
              }

              if(regions[0]->get_advanced_model()->QPEnabled && regions[i]->get_advanced_model()->QPEnabled)
              {
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQV));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->dg_variable_offset(EQV));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->dg_variable_offset(EQV));
              }

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
 * build function and its jacobian for Density Gradient solver
 */
void HomoInterfaceBC::DG_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // the Jacobian of HomoInterface boundary condition is processed here

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const FVM_Node *> fvm_nodes;
    std::vector<const SimulationRegion *> regions;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      fvm_nodes.push_back(fvm_node);
      regions.push_back(region);

      // the first semiconductor region
      if(i==0) continue;

      // other semiconductor region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              //the indepedent variable number, we need 2 here.
              adtl::AutoDScalar::numdir=2;

              // the governing equation of this fvm node
              AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];  V.setADValue(0,1.0);  // psi of this node
              AutoDScalar n = x[fvm_nodes[i]->local_offset()+1];  n.setADValue(0,1.0);  // electron density
              AutoDScalar p = x[fvm_nodes[i]->local_offset()+2];  p.setADValue(0,1.0);  // hole density

              AutoDScalar V_semi = x[fvm_nodes[0]->local_offset()+0]; V_semi.setADValue(1,1.0);
              AutoDScalar n_semi = x[fvm_nodes[0]->local_offset()+1]; n_semi.setADValue(1,1.0);  // electron density
              AutoDScalar p_semi = x[fvm_nodes[0]->local_offset()+2]; p_semi.setADValue(1,1.0);  // hole density

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              AutoDScalar ff1 = V - V_semi;
              jac->add( fvm_nodes[i]->global_offset()+0,  fvm_nodes[i]->global_offset()+0,  ff1.getADValue(0) );
              jac->add( fvm_nodes[i]->global_offset()+0,  fvm_nodes[0]->global_offset()+0,  ff1.getADValue(1) );

              AutoDScalar ff2 = n - n_semi;
              jac->add( fvm_nodes[i]->global_offset()+1,  fvm_nodes[i]->global_offset()+1,  ff2.getADValue(0) );
              jac->add( fvm_nodes[i]->global_offset()+1,  fvm_nodes[0]->global_offset()+1,  ff2.getADValue(1) );

              AutoDScalar ff3 = p - p_semi;
              jac->add( fvm_nodes[i]->global_offset()+2,  fvm_nodes[i]->global_offset()+2,  ff3.getADValue(0) );
              jac->add( fvm_nodes[i]->global_offset()+2,  fvm_nodes[0]->global_offset()+2,  ff3.getADValue(1) );

              if(regions[0]->get_advanced_model()->QNEnabled && regions[i]->get_advanced_model()->QNEnabled)
              {
                unsigned int node_eqc0_offset = regions[0]->dg_variable_offset(EQC);
                unsigned int node_eqc_offset = regions[i]->dg_variable_offset(EQC);
                AutoDScalar Eqc = x[fvm_nodes[i]->local_offset() + node_eqc_offset];      Eqc.setADValue(0,1.0);
                AutoDScalar Eqc_semi = x[fvm_nodes[0]->local_offset() + node_eqc0_offset]; Eqc_semi.setADValue(1,1.0);
                AutoDScalar ff = (Eqc - Eqc_semi);
                jac->add( fvm_nodes[i]->global_offset()+node_eqc_offset,  fvm_nodes[i]->global_offset()+node_eqc_offset,  ff.getADValue(0) );
                jac->add( fvm_nodes[i]->global_offset()+node_eqc_offset,  fvm_nodes[0]->global_offset()+node_eqc0_offset,  ff.getADValue(1) );
              }

              if(regions[0]->get_advanced_model()->QPEnabled && regions[i]->get_advanced_model()->QPEnabled)
              {
                unsigned int node_eqv0_offset = regions[0]->dg_variable_offset(EQV);
                unsigned int node_eqv_offset = regions[i]->dg_variable_offset(EQV);
                AutoDScalar Eqv = x[fvm_nodes[i]->local_offset() + node_eqv_offset];      Eqv.setADValue(0,1.0);
                AutoDScalar Eqv_semi = x[fvm_nodes[0]->local_offset() + node_eqv0_offset]; Eqv_semi.setADValue(1,1.0);
                AutoDScalar ff = (Eqv - Eqv_semi);
                jac->add( fvm_nodes[i]->global_offset()+node_eqv_offset,  fvm_nodes[i]->global_offset()+node_eqv_offset,  ff.getADValue(0) );
                jac->add( fvm_nodes[i]->global_offset()+node_eqv_offset,  fvm_nodes[0]->global_offset()+node_eqv0_offset,  ff.getADValue(1) );
              }

              break;
            }

            case InsulatorRegion :
            {
              //the indepedent variable number, we need 2 here.
              adtl::AutoDScalar::numdir=2;

              // the governing equation of this fvm node
              AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];  V.setADValue(0,1.0);  // psi of this node
              AutoDScalar V_semi = x[fvm_nodes[0]->local_offset()+0]; V_semi.setADValue(1,1.0);

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              AutoDScalar ff1 = V - V_semi;

              PetscInt row = fvm_nodes[i]->global_offset()+0;
              PetscInt cols[2] = { fvm_nodes[i]->global_offset()+0,  fvm_nodes[0]->global_offset()+0};
              jac->add_row(  row,  2,  cols,  ff1.getADValue() );

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

