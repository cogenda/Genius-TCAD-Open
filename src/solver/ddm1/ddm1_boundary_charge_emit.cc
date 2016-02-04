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


// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "boundary_condition_charge_emit.h"
#include "parallel.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::A;
using PhysicalUnit::Ohm;



///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * do pre-process to function for DDML1 solver
 */
void ChargeEmitBC::DDM1_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  this->_current_buffer.clear();
  this->_electron_current_buffer.clear();
  this->_hole_current_buffer.clear();

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with ohmic boundary condition
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;
      switch ( region->type() )
      {
      case SemiconductorRegion:
        {
          PetscInt row = fvm_node->global_offset();
          clear_row.push_back(row+1);

          // for conduction current
          {
            PetscInt    ix[2] = {fvm_node->global_offset()+1, fvm_node->global_offset()+2};
            // I={In, Ip} the electron and hole current flow into this boundary cell.
            // NOTE: although In has dn/dt and R items, they are zero since n is const and n=n0 holds
            // so does Ip
            PetscScalar I[2];

            VecGetValues(f, 2, ix, I);

            // the current = In - Ip;
            this->_current_buffer.push_back((I[0] - I[1]));
            this->_electron_current_buffer.push_back(I[0]);
            this->_hole_current_buffer.push_back(-I[1]);
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
 * build function and its jacobian for DDML1 solver
 */
void ChargeEmitBC::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Ohmic boundary condition is processed here.

  // we should do two things here. one is performance ohmic bc to corresponding mesh nodes.
  // another problem is the ohmic bc has an extra external circuit equation.
  // we must compute the total current flow in/out the ohmic bc.

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  const PetscScalar T = T_external();

  // data buffer for mesh nodes
  std::vector<PetscInt> iy;
  std::vector<PetscScalar> y;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const PetscScalar current_scale = this->z_width();

  // search and process all the boundary nodes
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

      const unsigned int global_offset = fvm_nodes[i]->global_offset();
      const unsigned int local_offset = fvm_nodes[i]->local_offset();

      switch ( regions[i]->type() )
      {
        // Semiconductor Region of course owns OhmicContactBC
      case SemiconductorRegion:
        {
          // semiconductor region should be the first region
          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();
          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          PetscScalar V = x[local_offset+0];  // psi of this node
          PetscScalar n = x[local_offset+1];  // electron density
          PetscScalar p = x[local_offset+2];  // hole density

          // mapping this node to material library
          semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

          PetscScalar ni  = semi_region->material()->band->ni(T);
          PetscScalar nie = semi_region->material()->band->nie(p, n, T);
          PetscScalar Nc  = semi_region->material()->band->Nc(T);
          PetscScalar Nv  = semi_region->material()->band->Nv(T);
          PetscScalar Eg  = semi_region->material()->band->Eg(T);
          PetscScalar dEg = semi_region->material()->band->EgNarrow(p, n, T);

          //Boltzmann
          {

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
          }

          // save insert position
          iy.push_back(global_offset+1);


          break;
        }

      case VacuumRegion:   break;
      default: genius_error(); //we should never reach here
      }
    }

  }

  // prevent insert zero length vector
  if(iy.size())  VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);


  const std::vector<PetscScalar> & current_buffer = this->_current_buffer;
  this->current()  = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void ChargeEmitBC::DDM1_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {

    // only process nodes belong to this processor
    if( (*node_it)->processor_id() != Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with ohmic boundary condition
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;
      switch ( region->type() )
      {
      case SemiconductorRegion:
        {
          PetscInt row = fvm_node->global_offset();
          clear_row.push_back(row+1);
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
 * build function and its jacobian for DDML1 solver
 */
void ChargeEmitBC::DDM1_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  const PetscScalar T = T_external();

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
          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

          //the indepedent variable number, we only need 4 here.
          adtl::AutoDScalar::numdir=4;
          //synchronize with material database
          semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

          AutoDScalar V = x[fvm_nodes[i]->local_offset()+0];    V.setADValue(0, 1.0);  // psi of this node
          AutoDScalar n = x[fvm_nodes[i]->local_offset()+1];    n.setADValue(1, 1.0);  // electron density
          AutoDScalar p = x[fvm_nodes[i]->local_offset()+2];    p.setADValue(2, 1.0);  // hole density

          semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

          PetscScalar ni  = semi_region->material()->band->ni(T);
          AutoDScalar nie = semi_region->material()->band->nie(p, n, T);
          PetscScalar Nc  = semi_region->material()->band->Nc(T);
          PetscScalar Nv  = semi_region->material()->band->Nv(T);
          PetscScalar Eg  = semi_region->material()->band->Eg(T);
          AutoDScalar dEg = semi_region->material()->band->EgNarrow(p, n, T);

          //governing equation for Ohmic contact boundary
          AutoDScalar ff2;
          //Boltzmann
          {
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
          }
          // the insert position
          PetscInt row[3], col[3];
          col[0] = row[0] = fvm_nodes[i]->global_offset()+0;
          col[1] = row[1] = fvm_nodes[i]->global_offset()+1;
          col[2] = row[2] = fvm_nodes[i]->global_offset()+2;

          // set Jacobian of governing equations
          jac->add_row(  row[1],  3,  &col[0],  ff2.getADValue() );

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


/*---------------------------------------------------------------------
 * update electrode IV
 */
void ChargeEmitBC::DDM1_Update_Solution(PetscScalar *)
{
  Parallel::sum(this->current());
}


