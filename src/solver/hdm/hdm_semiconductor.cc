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

#include "node.h"
#include "elem.h"
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "solver_specify.h"

#include "hdm_flux.h"
#include <TNT/jama_lu.h>

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::cm;


void SemiconductorSimulationRegion::HDM_Fill_Value(Vec x, Vec vol)
{

  PetscInt n_local_dofs;
  VecGetLocalSize(x, &n_local_dofs);

  std::vector<int> index;
  std::vector<PetscScalar> x_init;
  std::vector<PetscScalar> vol_init;

  index.reserve(n_local_dofs);
  x_init.reserve(n_local_dofs);
  vol_init.reserve(n_local_dofs);

  const PetscScalar mn = mt->band->EffecElecMass(T_external());
  const PetscScalar mp = mt->band->EffecHoleMass(T_external());
  const PetscScalar kT = kb*T_external();

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * node = *node_it;

    const FVM_NodeData * node_data = node->node_data();

    PetscScalar n = node_data->n();
    PetscScalar p = node_data->p();
    int seq=0;

    // the first variable, electron density n
    index.push_back(node->global_offset()+seq++);
    x_init.push_back(n);

    // the second vector, mnv
    index.push_back(node->global_offset()+seq++);
    index.push_back(node->global_offset()+seq++);
    index.push_back(node->global_offset()+seq++);
    x_init.push_back(mn*n*0.0);
    x_init.push_back(mn*n*0.0);
    x_init.push_back(mn*n*0.0);

    // energy
    //index.push_back(node->global_offset()+seq++);
    //x_init.push_back(1.5*kT*n + 0.5*mn*n*node_data->Vn()*node_data->Vn() );


    // the first variable, hole density p
    index.push_back(node->global_offset()+seq++);
    x_init.push_back(p);

    // the second vector, mpv
    index.push_back(node->global_offset()+seq++);
    index.push_back(node->global_offset()+seq++);
    index.push_back(node->global_offset()+seq++);
    x_init.push_back(mp*p*0.0);
    x_init.push_back(mp*p*0.0);
    x_init.push_back(mp*p*0.0);

    // energy
    //index.push_back(node->global_offset()+seq++);
    //x_init.push_back(1.5*kT*p + 0.5*mp*p*node_data->Vp()*node_data->Vp() );


    //volume, stupid code
    vol_init.push_back(1.0/node->volume());
    vol_init.push_back(1.0/node->volume());
    vol_init.push_back(1.0/node->volume());
    vol_init.push_back(1.0/node->volume());
    vol_init.push_back(1.0/node->volume());
    vol_init.push_back(1.0/node->volume());
    vol_init.push_back(1.0/node->volume());
    vol_init.push_back(1.0/node->volume());
  }

  if( index.size() )
  {
    VecSetValues(x, index.size(), &index[0], &x_init[0], INSERT_VALUES) ;
    VecSetValues(vol, index.size(), &index[0], &vol_init[0], INSERT_VALUES) ;
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}


void SemiconductorSimulationRegion::HDM_Flux(const PetscScalar * x, Vec flux, Vec t)
{
  const PetscScalar kT = kb*T_external();
  const PetscScalar mn = mt->band->EffecElecMass(T_external());
  const PetscScalar mp = mt->band->EffecHoleMass(T_external());

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * node = *node_it;

    Real dt=1e38;
    PetscInt loc[8];
    for(int i=0; i<8; ++i)
      loc[i] = node->global_offset() + i;

    HDMVector Un1(&x[node->local_offset()]);
    HDMVector Up1(&x[node->local_offset()+4]);


    FVM_Node::fvm_neighbor_node_iterator  neighbor_begin = node->neighbor_node_begin();
    FVM_Node::fvm_neighbor_node_iterator  neighbor_end = node->neighbor_node_end();
    for(; neighbor_begin!=neighbor_end; ++neighbor_begin )
    {
      const FVM_Node * neigbor_node = (*neighbor_begin).first;
      Real d = node->distance(neigbor_node);
      Real S = node->cv_surface_area(neigbor_node);
      VectorValue<double> dir = (*(neigbor_node->root_node()) - *(node->root_node())).unit();
      PetscScalar local_dt1, local_dt2;

      HDMVector Un2(&x[neigbor_node->local_offset()]);
      HDMVector Up2(&x[neigbor_node->local_offset()+4]);

      HDMVector fn = AUSM_if_flux(Un1, Un2, mn, kb, d, S, dir, local_dt1);
      HDMVector fp = AUSM_if_flux(Up1, Up2, mp, kb, d, S, dir, local_dt2);

      VecSetValues(flux, 4, &loc[0], &(fn[0]), ADD_VALUES);
      VecSetValues(flux, 4, &loc[4], &(fp[0]), ADD_VALUES);

      dt = std::min(dt, std::min(local_dt1, local_dt2));
    }
    std::vector<PetscScalar> dt_vector(8, dt);
    VecSetValues(t, 8, &loc[0], &dt_vector[0], INSERT_VALUES);
  }


#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}


void SemiconductorSimulationRegion::HDM_Source(const PetscScalar * lx, const PetscScalar * lt, Vec x)
{
  const PetscScalar T = T_external();
  const PetscScalar mn = mt->band->EffecElecMass(T);
  const PetscScalar mp = mt->band->EffecHoleMass(T);
  const PetscScalar damping_density = 1e17*std::pow(cm, -3);

  TNT::Array2D<Real> I(4, 4, 0.0);
  I[0][0]=I[1][1]=I[2][2]=I[3][3]=1.0;

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * node = *node_it;

    const FVM_NodeData * node_data = node->node_data();

    PetscInt loc[8];
    for(int i=0; i<8; ++i)
      loc[i] = node->global_offset() + i;

    PetscScalar dt        = lt[node->local_offset()];

    PetscScalar n         = lx[node->local_offset()];
    PetscScalar p         = lx[node->local_offset()+4];
    PetscScalar R         = mt->band->Recomb(p, n, T);
    PetscScalar mun       = mt->mob->ElecMob(p, n, T, 0.0, 0.0, T);
    PetscScalar mup       = mt->mob->HoleMob(p, n, T, 0.0, 0.0, T);

    PetscScalar carrier = node_data->n() + node_data->p();
    PetscScalar damping = carrier > damping_density ? damping_density/carrier : 1.0;
    PetscScalar taon = mn*mun/e;
    PetscScalar taop = mp*mup/e;

    {
      TNT::Array2D<Real> An (4, 4, 0.0);
      TNT::Array1D<Real> bn (4);
      TNT::Array1D<Real> rn (4);

      An[0][0] = -R/n;
      An[1][0] = -node_data->E()(0);
      An[1][1] = -1/taon;
      An[2][0] = -node_data->E()(1);
      An[2][2] = -1/taon;
      An[3][0] = -node_data->E()(2);
      An[3][3] = -1/taon;

      bn[0] = lx[node->local_offset()+0];
      bn[1] = lx[node->local_offset()+1];
      bn[2] = lx[node->local_offset()+2];
      bn[3] = lx[node->local_offset()+3];

      rn = dt*(An*bn);
      An = I - dt*An;

      JAMA::LU<Real> solver(An);
      TNT::Array1D<Real> dx = solver.solve(rn);
      //dx = 0.1*dx;
      VecSetValues(x, 4, &loc[0], &dx[0], ADD_VALUES);
    }

    {
      TNT::Array2D<Real> Ap (4, 4, 0.0);
      TNT::Array1D<Real> bp (4);
      TNT::Array1D<Real> rp (4);

      Ap[0][0] = -R/p;
      Ap[1][0] = node_data->E()(0);
      Ap[1][1] = -1/taop;
      Ap[2][0] = node_data->E()(1);
      Ap[2][2] = -1/taop;
      Ap[3][0] = node_data->E()(2);
      Ap[3][3] = -1/taop;

      bp[0] = lx[node->local_offset()+4];
      bp[1] = lx[node->local_offset()+5];
      bp[2] = lx[node->local_offset()+6];
      bp[3] = lx[node->local_offset()+7];

      rp = dt*(Ap*bp);
      Ap = I - dt*Ap;

      JAMA::LU<Real> solver(Ap);
      TNT::Array1D<Real> dx = solver.solve(rp);
      //dx = 0.1*dx;
      VecSetValues(x, 4, &loc[4], &dx[0], ADD_VALUES);
    }
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}



#if 0
void SemiconductorSimulationRegion::HDM_Source(const PetscScalar * lx, const PetscScalar * lt, Vec x)
{
  const PetscScalar T = T_external();
  const PetscScalar Vt = kb*T/e;
  const PetscScalar mn = mt->band->EffecElecMass(T);
  const PetscScalar mp = mt->band->EffecHoleMass(T);

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * node = *node_it;
    const FVM_NodeData * node_data = node->node_data();

    PetscInt loc[8];
    for(int i=0; i<8; ++i)
      loc[i] = node->global_offset() + i;

    PetscScalar V         = node_data->psi();
    PetscScalar dt        = lt[node->local_offset()];
    PetscScalar n         = lx[node->local_offset()];
    PetscScalar p         = lx[node->local_offset()+4];
    VectorValue<PetscScalar> vn(&lx[node->local_offset()+1]);
    VectorValue<PetscScalar> vp(&lx[node->local_offset()+5]);
    PetscScalar R         = mt->band->Recomb(p, n, T);
    PetscScalar mun       = mt->mob->ElecMob(p, n, T, 0.0, 0.0, T);
    PetscScalar mup       = mt->mob->HoleMob(p, n, T, 0.0, 0.0, T);

    PetscScalar taon = mn*mun/e;
    PetscScalar taop = mp*mup/e;

    PetscScalar Sn = n - R*dt;
    PetscScalar Sp = p - R*dt;

    VectorValue<PetscScalar> vne, vpe;
    FVM_Node::fvm_neighbor_node_iterator  neighbor_begin = node->neighbor_node_begin();
    FVM_Node::fvm_neighbor_node_iterator  neighbor_end = node->neighbor_node_end();
    for(; neighbor_begin!=neighbor_end; ++neighbor_begin )
    {
      const FVM_Node * neigbor_node = (*neighbor_begin).second;
      const FVM_NodeData * neigbor_node_data = neigbor_node->node_data();
      Real S = node->cv_surface_area(neigbor_node);
      VectorValue<double> dir = (*(neigbor_node->root_node()) - *(node->root_node())).unit();
      PetscScalar VV         = neigbor_node_data->psi();
      PetscScalar nn         = lx[neigbor_node->local_offset()];
      PetscScalar pp         = lx[neigbor_node->local_offset()+4];

      PetscScalar nmid = n;//nmid_dd(Vt, V, VV, n, nn);
      PetscScalar pmid = p;//pmid_dd(Vt, V, VV, p, pp);

      vne +=  dt*nmid*0.5*(V+VV)*S*dir;
      vpe += -dt*pmid*0.5*(V+VV)*S*dir;
    }

    vne /= node->volume();
    vpe /= node->volume();

    VectorValue<PetscScalar> Svn = (vn+vne)/(1+dt/taon);
    VectorValue<PetscScalar> Svp = (vp+vpe)/(1+dt/taop);

    VecSetValues(x, 1, &loc[0], &Sn, INSERT_VALUES);
    VecSetValues(x, 3, &loc[1], &Svn[0], INSERT_VALUES);
    VecSetValues(x, 1, &loc[4], &Sp, INSERT_VALUES);
    VecSetValues(x, 3, &loc[5], &Svp[0], INSERT_VALUES);
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}

#endif


void SemiconductorSimulationRegion::HDM_Update_Solution(const PetscScalar * x)
{
  const PetscScalar mn = mt->band->EffecElecMass(T_external());
  const PetscScalar mp = mt->band->EffecHoleMass(T_external());

  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * node = *node_it;
    FVM_NodeData * node_data = node->node_data();

    node_data->n_last()   = node_data->n();
    node_data->n()        = x[node->local_offset()];
    node_data->Jn()       = VectorValue<Real>(&x[node->local_offset()+1])/mn;

    node_data->p_last()   = node_data->p();
    node_data->p()        = x[node->local_offset()+4];
    node_data->Jp()       = VectorValue<Real>(&x[node->local_offset()+5])/mp;

    node_data->rho()      = node_data->Net_doping() - node_data->n() + node_data->p();
  }
}


