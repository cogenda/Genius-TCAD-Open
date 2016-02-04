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

#include "simulation_system.h"
#include "semiconductor_region.h"
#include "boundary_condition_homo.h"
#include "petsc_utils.h"
#include "parallel.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;

void HomoInterfaceBC::_ddm_current_interface(PetscScalar *x, Vec f)
{
  std::vector<PetscScalar>   electron_current_buffer;
  std::vector<PetscScalar>   hole_current_buffer;
  std::vector<PetscScalar>   displacement_current_buffer;

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
      const FVM_NodeData * node_data = fvm_node->node_data();

      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0) continue;

      if(region->type() != SemiconductorRegion) continue;

      // for conduction current
      {
        PetscInt    ix[2] = {fvm_nodes[i]->global_offset()+1, fvm_nodes[i]->global_offset()+2};
        // I={In, Ip} the electron and hole current flow into this boundary cell.
        PetscScalar I[2];

        VecGetValues(f, 2, ix, I);

        // the current = In - Ip;
        electron_current_buffer.push_back(I[0]);
        hole_current_buffer.push_back(-I[1]);
      }

      // displacement current
      if(SolverSpecify::TimeDependent == true)
      {
        PetscScalar I_displacement = 0.0;
        // the psi of this node
        PetscScalar V = x[fvm_nodes[i]->local_offset()+0];
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
        for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
        {
          const FVM_Node *nb_node = (*nb_it).first;
          const FVM_NodeData * nb_node_data = nb_node->node_data();
          // the psi of neighbor node
          PetscScalar V_nb = x[nb_node->local_offset()+0];
          // distance from nb node to this node
          PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node);
          PetscScalar dEdt;
          if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false) //second order
          {
            PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
            dEdt = ( (2-r)/(1-r)*(V-V_nb)
                     - 1.0/(r*(1-r))*(node_data->psi()-nb_node_data->psi())
                     + (1-r)/r*(node_data->psi_last()-nb_node_data->psi_last()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
          }
          else//first order
          {
            dEdt = ((V-V_nb)-(node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
          }
          I_displacement += cv_boundary*fvm_nodes[i]->node_data()->eps() *dEdt;
        }
        displacement_current_buffer.push_back( I_displacement );
      }
    }
  }

  const PetscScalar current_scale = this->z_width();
  PetscScalar electron_current  = current_scale*std::accumulate(electron_current_buffer.begin(), electron_current_buffer.end(), 0.0 );
  PetscScalar hole_current  = current_scale*std::accumulate(hole_current_buffer.begin(), hole_current_buffer.end(), 0.0 );
  PetscScalar displacement_current  = current_scale*std::accumulate(displacement_current_buffer.begin(), displacement_current_buffer.end(), 0.0 );

  Parallel::sum(electron_current);
  Parallel::sum(hole_current);
  Parallel::sum(displacement_current);

  this->scalar("electron_current") = electron_current;
  this->scalar("hole_current") = hole_current;
  this->scalar("displacement_current") = displacement_current;

}

