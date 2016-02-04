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

//  $Id: ddm1_boundary_ohmic.cc,v 1.19 2008/07/09 07:45:27 gdiso Exp $

// C++ includes
#include <numeric>

#include "asinh.hpp" // for asinh

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_ohmic.h"
#include "parallel.h"
#include "mathfunc.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::A;
using PhysicalUnit::Ohm;

/*---------------------------------------------------------------------
 * fill ohmic electrode potential into initial vector
 */
void OhmicContactBC::DDM1_Fill_Value(Vec x, Vec L)
{
  const double current_scale = this->z_width();

  if(Genius::processor_id() == Genius::n_processors() -1)
  {
    VecSetValue(x, this->global_offset(), this->ext_circuit()->potential(), INSERT_VALUES);

    if(this->is_inter_connect_bc())
    {
      VecSetValue(L, this->global_offset(), 1.0, INSERT_VALUES);
    }
    //for stand alone electrode
    else
    {
      const PetscScalar s = ext_circuit()->electrode_scaling(SolverSpecify::dt);
      VecSetValue(L, this->global_offset(), s, INSERT_VALUES);
    }

  }

}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * do pre-process to function for DDML1 solver
 */
void OhmicContactBC::DDM1_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
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
            clear_row.push_back(row+0);
            clear_row.push_back(row+1);
            clear_row.push_back(row+2);

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
          case ElectrodeRegion:
          case MetalRegion:
          case InsulatorRegion:
          {
            PetscInt row = fvm_node->global_offset();
            clear_row.push_back(row);
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
void OhmicContactBC::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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
  std::vector<PetscScalar>   displacement_current_buffer;

  // the electrode potential in current iteration
  genius_assert( local_offset()!=invalid_uint );
  const PetscScalar Ve = x[this->local_offset()];

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
            PetscScalar dEc = semi_region->material()->band->EgNarrowToEc(p, n, T);
            PetscScalar dEv = semi_region->material()->band->EgNarrowToEv(p, n, T);
            //governing equation for Ohmic contact boundary
            if(semi_region->get_advanced_model()->Fermi) //Fermi
            {
              PetscScalar Ec =  -(e*V + node_data->affinity() + dEc);
              PetscScalar Ev =  -(e*V + node_data->affinity() + Eg - dEv);

              // the quasi-fermi potential equals to electrode Vapp
              PetscScalar phin = Ve;
              PetscScalar phip = Ve;

              PetscScalar etan = (-e*phin-Ec)/kb/T;
              PetscScalar etap = (Ev+e*phip)/kb/T;

              y.push_back( Nc*fermi_half(etan) - Nv*fermi_half(etap) -node_data->Net_doping() );
              y.push_back( n - Nc*fermi_half(etan) );
              y.push_back( p - Nv*fermi_half(etap) );
            }
            else     //Boltzmann
            {
              y.push_back  ( V - kb*T/e*asinh(node_data->Net_doping()/(2*nie))
                             + Eg/(2*e)
                             + kb*T*log(Nc/Nv)/(2*e)
                             + node_data->affinity()/e
                             - Ve
                           );                        //governing equation for psi


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
              y.push_back( p - hole_density );      //governing equation for hole density
            }

            // save insert position
            iy.push_back(global_offset+0);
            iy.push_back(global_offset+1);
            iy.push_back(global_offset+2);


            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              PetscScalar I_displacement = 0.0;
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
                I_displacement += cv_boundary*node_data->eps() *dEdt;
              }
              displacement_current_buffer.push_back( I_displacement );
            }
            break;
          }
          // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
          // the boundary for the corner point may be Ohmic. (not a nice behavier)
          case InsulatorRegion:
          {
            // psi of this node
            PetscScalar V = x[local_offset];

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            PetscScalar V_ref = x[fvm_nodes[0]->local_offset()];

            // the psi of this node is equal to corresponding psi of semiconductor node
            y.push_back( V - V_ref );

            // save insert position
            iy.push_back(global_offset);

            break;
          }
          // conductor region which has an interface with OhmicContact boundary to semiconductor region
          case MetalRegion:
          case ElectrodeRegion:
          {

            // psi of this node
            PetscScalar V = x[local_offset];

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            PetscScalar V_ref = x[fvm_nodes[0]->local_offset()];

            // the psi of this node is equal to corresponding psi of semiconductor node
            y.push_back( V - V_ref );

            // save insert position
            iy.push_back(global_offset);

            break;
          }

          case VacuumRegion:   break;
          default: genius_error(); //we should never reach here
      }
    }

  }

  // prevent insert zero length vector
  if(iy.size())  VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);


  // the extra equation of ohmic boundary
  // For voltage driven
  //
  //          _____                Ve
  //    -----|_____|----/\/\/\/\-------> to ohmic electrode (Ve, I)
  //    | +     R          L       |
  //   Vapp                     C ===
  //    | -                        |
  //    |__________________________|
  //           GND
  //
  // And for current driven
  //                               Ve
  //    -->-----------------------------> to ohmic electrode (Ve, I)
  //    |                          |
  //   Iapp                     C ===
  //    |__________________________|
  //           GND
  //
  // Or for inter connect
  //
  //          _____                Ve
  //    -----|_____|-------------------> to ohmic electrode (Ve, I)
  //    |       R
  //    |
  // V_inter_connect
  //
  //


  // for get the current, we must sum all the terms in current_buffer
  // NOTE: only statistic current flow belongs to on processor node

  const std::vector<PetscScalar> & current_buffer = this->_current_buffer;
  PetscScalar current_conductance  = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );
  PetscScalar current_displacement = current_scale*std::accumulate(displacement_current_buffer.begin(), displacement_current_buffer.end(), 0.0 );
  PetscScalar current = current_conductance + current_displacement;

  ext_circuit()->potential() = Ve;
  ext_circuit()->current() = current;
  ext_circuit()->current_displacement() =  current_displacement;
  ext_circuit()->current_electron() =  current_scale*std::accumulate(this->_electron_current_buffer.begin(), this->_electron_current_buffer.end(), 0.0 );
  ext_circuit()->current_hole() =  current_scale*std::accumulate(this->_hole_current_buffer.begin(), this->_hole_current_buffer.end(), 0.0 );

  PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);

  //for inter connect electrode
  if(this->is_inter_connect_bc())
  {
    PetscScalar R = ext_circuit()->inter_connect_resistance();                               // resistance
    PetscScalar f_ext = R*current;
    VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
  }
  // for stand alone electrode
  else
  {
    PetscScalar f_ext = mna_scaling*current;
    VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
  }

  if(Genius::is_last_processor())
  {
    //for inter connect electrode
    if(this->is_inter_connect_bc())
    {
      PetscScalar V_ic = x[this->inter_connect_hub()->local_offset()];  // potential at inter connect node
      PetscScalar f_ext = Ve - V_ic;
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
    // for stand alone electrode
    else
    {
      PetscScalar f_ext = ext_circuit()->mna_function(SolverSpecify::dt);
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
  }



  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}







/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void OhmicContactBC::DDM1_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
                                              std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  _buffer_cols.clear();
  _buffer_jacobian_entries.clear();

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();

  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // only process nodes belong to this processor
    if( (*node_it)->processor_id() != Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(; rnode_it!=end_rnode_it; ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node *  fvm_node = (*rnode_it).second.second;

      if ( region->type() == SemiconductorRegion)
      {
        const FVM_NodeData * node_data = fvm_node->node_data();

        PetscScalar A1[3], A2[3];
        std::vector<PetscScalar> JM(3), JN(3);
        std::vector<PetscInt>    row(3);
        row[0] = fvm_node->global_offset()+0;
        row[1] = fvm_node->global_offset()+1;
        row[2] = fvm_node->global_offset()+2;

        //NOTE only get value from local block!
        jac->get_row(  row[1],  3,  &row[0],  &A1[0]);
        jac->get_row(  row[2],  3,  &row[0],  &A2[0]);

        JM[0] = (A1[0]-A2[0]);
        JM[1] = (A1[1]-A2[1]);
        JM[2] = (A1[2]-A2[2]);

        // get the derivative of electrode current to neighbors of ohmic node
        // NOTE neighbors and ohmic bc node may on different processor!
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
        for(; nb_it != nb_it_end; ++nb_it)
        {
          const FVM_Node *  fvm_nb_node = (*nb_it).first;

          std::vector<PetscInt>    col(3);
          col[0] = fvm_nb_node->global_offset()+0;
          col[1] = fvm_nb_node->global_offset()+1;
          col[2] = fvm_nb_node->global_offset()+2;

          jac->get_row(  row[1],  3,  &col[0],  &A1[0]);
          jac->get_row(  row[2],  3,  &col[0],  &A2[0]);

          JN[0] = (A1[0]-A2[0]);
          JN[1] = (A1[1]-A2[1]);
          JN[2] = (A1[2]-A2[2]);

          _buffer_cols.push_back(col);
          _buffer_jacobian_entries.push_back(JN);
        }

        _buffer_cols.push_back(row);
        _buffer_jacobian_entries.push_back(JM);
      }
    }
  }

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
            clear_row.push_back(row+0);
            clear_row.push_back(row+1);
            clear_row.push_back(row+2);
            break;
          }
          case MetalRegion:
          case ElectrodeRegion:
          case InsulatorRegion:
          {
            PetscInt row = fvm_node->global_offset();
            clear_row.push_back(row);
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
void OhmicContactBC::DDM1_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // the Jacobian of Ohmic boundary condition is processed here

  PetscInt bc_global_offset = this->global_offset();

  const PetscScalar T = T_external();

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const PetscScalar current_scale = this->z_width();

  // d(current)/d(independent variables of bd node and its neighbors)
  {
    PetscScalar bc_current_scale = 1.0;
    if(this->is_inter_connect_bc())
    {
      PetscScalar R = ext_circuit()->inter_connect_resistance();
      bc_current_scale *= R*current_scale;
    }
    // for stand alone electrode
    else
    {
      PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);
      bc_current_scale *= mna_scaling*current_scale;
    }
    for(unsigned int n=0; n<_buffer_cols.size(); ++n)
    {
      std::vector<PetscScalar> bc_current_jacobian = _buffer_jacobian_entries[n];
      for(unsigned int k=0; k<bc_current_jacobian.size(); ++k)
        bc_current_jacobian[k] *= bc_current_scale;
      jac->add_row(  bc_global_offset,  _buffer_cols[n].size(),  &(_buffer_cols[n])[0],  &(bc_current_jacobian)[0] );
    }
  }


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

            // the electrode potential in current iteration
            AutoDScalar Ve = x[this->local_offset()];             Ve.setADValue(3, 1.0);

            semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

            PetscScalar ni  = semi_region->material()->band->ni(T);
            AutoDScalar nie = semi_region->material()->band->nie(p, n, T);
            PetscScalar Nc  = semi_region->material()->band->Nc(T);
            PetscScalar Nv  = semi_region->material()->band->Nv(T);
            PetscScalar Eg  = semi_region->material()->band->Eg(T);
            AutoDScalar dEg = semi_region->material()->band->EgNarrow(p, n, T);
            AutoDScalar dEc = semi_region->material()->band->EgNarrowToEc(p, n, T);
            AutoDScalar dEv = semi_region->material()->band->EgNarrowToEv(p, n, T);

            //governing equation for Ohmic contact boundary
            AutoDScalar ff1,ff2,ff3;
            if(semi_region->get_advanced_model()->Fermi) //Fermi
            {
              AutoDScalar Ec =  -(e*V + node_data->affinity() + dEc);
              AutoDScalar Ev =  -(e*V + node_data->affinity() + Eg - dEv);

              // the quasi-fermi potential equals to electrode Vapp
              AutoDScalar phin = Ve;
              AutoDScalar phip = Ve;

              AutoDScalar etan = (-e*phin-Ec)/kb/T;
              AutoDScalar etap = (Ev+e*phip)/kb/T;

              ff1 =  Nc*fermi_half(etan) - Nv*fermi_half(etap) -node_data->Net_doping();
              ff2 =  n - Nc*fermi_half(etan);
              ff3 =  p - Nv*fermi_half(etap);
            }
            else //Boltzmann
            {

              ff1 = V - kb*T/e*adtl::asinh(node_data->Net_doping()/(2*nie))
                    + Eg/(2*e)
                    + kb*T*log(Nc/Nv)/(2*e)
                    + node_data->affinity()/e
                    - Ve;

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
              ff3 =  p - hole_density;      //governing equation for hole density
            }
            // the insert position
            PetscInt row[3], col[4];
            col[0] = row[0] = fvm_nodes[i]->global_offset()+0;
            col[1] = row[1] = fvm_nodes[i]->global_offset()+1;
            col[2] = row[2] = fvm_nodes[i]->global_offset()+2;
            col[3] = this->global_offset(); // the position of electrode equation

            // set Jacobian of governing equations
            jac->add_row(  row[0],  4,  &col[0],  ff1.getADValue() );
            jac->add_row(  row[1],  4,  &col[0],  ff2.getADValue() );
            jac->add_row(  row[2],  4,  &col[0],  ff3.getADValue() );

            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                FVM_Node *nb_node = (*nb_it).first;
                FVM_NodeData * nb_node_data = nb_node->node_data();
                // the psi of neighbor node
                AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);
                // distance from nb node to this node
                PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
                // area of out surface of control volume related with neighbor node
                PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node);

                AutoDScalar dEdt;
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

                AutoDScalar current_disp = cv_boundary*node_data->eps()*dEdt*current_scale;

                //for inter connect electrode
                if(this->is_inter_connect_bc())
                {
                  PetscScalar R = ext_circuit()->inter_connect_resistance();
                  current_disp = R*current_disp;
                }
                // for stand alone electrode
                else
                {
                  PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);
                  current_disp = mna_scaling*current_disp;
                }

                jac->add( bc_global_offset,  fvm_nodes[i]->global_offset()+0,  current_disp.getADValue(0) );
                jac->add( bc_global_offset,  nb_node->global_offset()+0,  current_disp.getADValue(1) );
              }
            }

            break;
          }
          // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
          // the boundary for the corner point may be Ohmic.
          case InsulatorRegion:
          {
            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

            AutoDScalar  V_ref = x[fvm_nodes[0]->local_offset()]; V_ref.setADValue(1,1.0);

            // the psi of this node is equal to corresponding psi of semiconductor node
            AutoDScalar  ff = V - V_ref;

            // set Jacobian of governing equation ff
            jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[i]->global_offset(),  ff.getADValue(0) );
            jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[0]->global_offset(),  ff.getADValue(1) );

            break;
          }
          // conductor region which has an interface with OhmicContact boundary to semiconductor region
          case MetalRegion:
          case ElectrodeRegion:
          {
            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

            AutoDScalar  V_ref = x[fvm_nodes[0]->local_offset()]; V_ref.setADValue(1,1.0);

            // the psi of this node is equal to corresponding psi of semiconductor node
            AutoDScalar  ff = V - V_ref;

            // set Jacobian of governing equation ff
            jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[i]->global_offset(),  ff.getADValue(0) );
            jac->add( fvm_nodes[i]->global_offset(),  fvm_nodes[0]->global_offset(),  ff.getADValue(1) );

            break;
          }

          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }

  }


  // the extra equation of ohmic boundary
  // For voltage driven
  //
  //          _____                Ve
  //    -----|_____|----/\/\/\/\-------> to ohmic electrode (Ve, I)
  //    | +     R          L       |
  //   Vapp                     C ===
  //    | -                        |
  //    |__________________________|
  //           GND
  //
  // And for current driven
  //                               Ve
  //    -->-----------------------------> to ohmic electrode (Ve, I)
  //    |                          |
  //   Iapp                     C ===
  //    |__________________________|
  //           GND
  //
  // Or for inter connect
  //
  //          _____                Ve
  //    -----|_____|-------------------> to ohmic electrode (Ve, I)
  //    |       R
  //    |
  // V_inter_connect
  //
  //

  if(Genius::is_last_processor())
  {
    //for inter connect electrode
    if(this->is_inter_connect_bc())
    {
      // the external electrode equation is:
      // f_ext = Ve - V_ic + R*current;

      // d(f_ext)/d(Ve)
      jac->add( bc_global_offset,  bc_global_offset,  1.0 );
      // d(f_ext)/d(V_ic)
      jac->add( bc_global_offset,  this->inter_connect_hub()->global_offset(),  -1.0 );
    }
    //for stand alone electrode
    else
    {
      ext_circuit()->potential() = x[this->local_offset()];
      jac->add( bc_global_offset,  bc_global_offset,  ext_circuit()->mna_jacobian(SolverSpecify::dt) );
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


void OhmicContactBC::DDM1_Electrode_Trace(Vec, SparseMatrix<PetscScalar> *jac, Vec pdI_pdx, Vec pdF_pdV)
{
  VecZeroEntries(pdI_pdx);
  VecZeroEntries(pdF_pdV);

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();
  for(unsigned int n=0; n<_buffer_cols.size(); ++n)
  {
    std::vector<PetscScalar> bc_current_jacobian = _buffer_jacobian_entries[n];
    for(unsigned int k=0; k<bc_current_jacobian.size(); ++k)
      bc_current_jacobian[k] *= current_scale;
    VecSetValues(pdI_pdx, _buffer_cols[n].size(), &(_buffer_cols[n])[0], &(bc_current_jacobian)[0], ADD_VALUES);
  }


  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // get the derivative of electrode current to ohmic node
    const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, SemiconductorRegion);
    const FVM_NodeData * node_data = fvm_node->node_data();

    VecSetValue( pdF_pdV, fvm_node->global_offset(), 1.0, ADD_VALUES);
  }

  VecAssemblyBegin(pdI_pdx);
  VecAssemblyBegin(pdF_pdV);

  VecAssemblyEnd(pdI_pdx);
  VecAssemblyEnd(pdF_pdV);


  //delete electrode current equation, omit the effect of external resistance
  PetscInt bc_global_offset = this->global_offset();
  jac->clear_row(bc_global_offset, 1.0);
}


/*---------------------------------------------------------------------
 * update electrode IV
 */
void OhmicContactBC::DDM1_Update_Solution(PetscScalar *)
{
  Parallel::sum(ext_circuit()->current());
  this->ext_circuit()->update();

  // also statistic displacement current, electron current and hole current
  Parallel::sum(ext_circuit()->current_displacement());
  Parallel::sum(ext_circuit()->current_electron());
  Parallel::sum(ext_circuit()->current_hole());
}


