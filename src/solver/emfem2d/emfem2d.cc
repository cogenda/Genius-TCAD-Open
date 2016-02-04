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

#include <fstream>
#include <sstream>

#include "mesh_base.h"
#include "emfem2d/emfem2d.h"
#include "petsc_type.h"
#include "fe_type.h"
#include "fe_base.h"
#include "quadrature_gauss.h"
#include "dense_vector.h"
#include "dense_matrix.h"
#include "petsc_utils.h"
#include "spline.h"
#include "parallel.h"

using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::J;
using PhysicalUnit::s;
using PhysicalUnit::um;
using PhysicalUnit::cm;
using PhysicalUnit::eps0;
using PhysicalUnit::mu0;
using PhysicalUnit::h;

/*------------------------------------------------------------------
 * create the poisson solver contex
 */
int EMFEM2DSolver::create_solver()
{
#if 0
  MESSAGE<< '\n' << "EM FEM 2D Solver init..." << std::endl;
  RECORD();

  // set ac variables for each region
  set_variables();

  //parse the command card here
  setup_solver_parameters();

  // must set linear matrix/vector here!
  setup_linear_data();

  // rtol   = 1e-10*n_global_dofs  - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol = 1e-20*n_global_dofs  - the absolute convergence tolerance (absolute size of the residual norm)
  KSPSetTolerances(ksp, 1e-10*n_global_dofs, 1e-20*n_global_dofs, PETSC_DEFAULT, n_global_dofs/10);

  // user can do further adjusment from command line
  KSPSetFromOptions (ksp);

#endif

  return 0;
}


/*------------------------------------------------------------------
 * prepare solution and aux variables used by this solver
 */
int EMFEM2DSolver::set_variables()
{
#if 0
  for ( unsigned int n=0; n<_system.n_regions(); n++ )
  {
    SimulationRegion * region = _system.region ( n );
    region->add_variable("optical_efield", POINT_CENTER);
    region->add_variable("optical_hfield", POINT_CENTER);
  }
#endif
  return 0;
}


/*------------------------------------------------------------------
 *
 */
int EMFEM2DSolver::solve()
{
#if 0
  START_LOG("EM FEM 2D Linear Solver", "solve");

  // for each wave length, solve 2d fem probelm
  for(unsigned int n=0; n<_optical_sources.size(); ++n)
  {
    if(_optical_sources[n].TE_weight>0)
    {
      solve_TE_scatter_problem(_optical_sources[n].wave_length,
                               _optical_sources[n].power*_optical_sources[n].TE_weight,
                               _optical_sources[n].phi_TE);

      //update solution
      save_TE_solution(_optical_sources[n].wave_length,
                       _optical_sources[n].power*_optical_sources[n].TE_weight,
                       _optical_sources[n].phi_TE,
                       _optical_sources[n].eta,
                       _optical_sources[n].eta_auto);
    }

    if(_optical_sources[n].TM_weight>0)
    {
      solve_TM_scatter_problem(_optical_sources[n].wave_length,
                               _optical_sources[n].power*_optical_sources[n].TM_weight,
                               _optical_sources[n].phi_TM);

      //update solution
      save_TM_solution(_optical_sources[n].wave_length,
                       _optical_sources[n].power*_optical_sources[n].TM_weight,
                       _optical_sources[n].phi_TM,
                       _optical_sources[n].eta,
                       _optical_sources[n].eta_auto);
    }

  }

  STOP_LOG("EM FEM 2D Linear Solver", "solve");
#endif
  return 0;
}



int EMFEM2DSolver::destroy_solver()
{
  // clear linear contex
  clear_linear_data();
  return 0;
}



void EMFEM2DSolver::solve_TM_scatter_problem(double lambda, double power, double phase0)
{
#if 0
  MESSAGE<<"Solve TM Mode. WaveLength = "<<lambda/um<< " um, Power = " <<power/(J/s/cm/cm)<<" W/(cm^2)." << std::endl; RECORD();

  build_TM_matrix_rhs(lambda, power, phase0);

#if PETSC_VERSION_GE(3,5,0)
  KSPSetOperators(ksp,A,A);
#else
  KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
#endif
  KSPSolve(ksp,b,x);

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  PetscInt   its;
  KSPGetIterationNumber(ksp, &its);

  PetscReal  rnorm;
  KSPGetResidualNorm(ksp, &rnorm);

  MESSAGE<<"------> residual norm = "<<rnorm<<" its = "<<its<<" with "<<KSPConvergedReasons[reason]<<"\n\n";
  RECORD();
#endif
}



void EMFEM2DSolver::solve_TE_scatter_problem(double lambda, double power, double phase0)
{
#if 0
  MESSAGE<<"Solve TE Mode. WaveLength = "<<lambda/um<< " um, Power = " <<power/(J/s/cm/cm)<<" W/(cm^2)." << std::endl; RECORD();

  build_TE_matrix_rhs(lambda, power, phase0);

#if PETSC_VERSION_GE(3,5,0)
  KSPSetOperators(ksp,A,A);
#else
  KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
#endif
  KSPSolve(ksp,b,x);

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  PetscInt   its;
  KSPGetIterationNumber(ksp, &its);

  PetscReal  rnorm;
  KSPGetResidualNorm(ksp, &rnorm);

  MESSAGE<<"------> residual norm = "<<rnorm<<" its = "<<its<<" with "<<KSPConvergedReasons[reason]<<"\n\n";
  RECORD();
#endif
}


void EMFEM2DSolver::build_TE_matrix_rhs(double lambda, double power, double phase0)
{
#if 0
  //wave vector
  double k = 2*M_PI/lambda;

  // wave magnitude, compute from power
  double H = sqrt(power*sqrt(eps0/mu0));

  Complex j(0,1);

  const MeshBase& mesh = _system.mesh();
  const unsigned int dim = mesh.mesh_dimension();
  genius_assert(dim==2);

  // the integral order used in gauss intergral. default value= 2*fe_order+1
  Order int_order=SECOND;
  FEType fe_type;

  VecZeroEntries(x);
  VecZeroEntries(b);
  MatZeroEntries(A);

  // process scatter field
  {
    AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

    // Gauss quadrature rule for numerical integration.
    QGauss qrule (dim, int_order);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<Real> >& dphidx       = fe->get_dphidx();
    const std::vector<std::vector<Real> >& dphidy       = fe->get_dphidy();
    const std::vector<std::vector<Real> >& dphidz       = fe->get_dphidz();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.
    DenseMatrix<Complex> Ke;
    DenseVector<Complex> Fe;

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
    for(unsigned int n=0; n<_system.n_regions(); ++n)
    {
      SimulationRegion * region = _system.region(n);

      Complex r = region->get_optical_refraction(lambda);
      Complex eps(r.real()*r.real()-r.imag()*r.imag(), -2*r.real()*r.imag()) ;
      double mu=1.0;

      //for all the elements in this region
      // note, they are all local element, thus must be processed
      SimulationRegion::element_iterator it = region->elements_begin();
      SimulationRegion::element_iterator it_end = region->elements_end();
      for(; it!=it_end; ++it)
      {
        const Elem* elem  = *it;
        genius_assert(elem->active());
        if(elem->processor_id()!=Genius::processor_id()) continue;

        std::vector<PetscInt> dof_indices;
        this->build_dof_indices(elem, dof_indices);

        fe->reinit (elem);

        // one complex variable per noe
        Ke.resize (elem->n_nodes(), elem->n_nodes());
        Fe.resize (elem->n_nodes());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          for (unsigned int m=0; m<phi.size(); m++)
          {
            for (unsigned int n=0; n<phi.size(); n++)
            {
              // \nabla \cdot frac{1}{eps} \nabla H^_{sc}
              Ke(m,n) += - JxW[qp]/eps*(  dphidx[m][qp]*dphidx[n][qp]
                                          + dphidy[m][qp]*dphidy[n][qp]
                                          + dphidz[m][qp]*dphidz[n][qp]);
              // k^2*eps*H^_{sc}
              Ke(m,n) += JxW[qp]*k*k*mu*phi[m][qp]*phi[n][qp];
            }
            // source item
            const Node * node = elem->get_node(m);
            double phase = phase0 - k*(node->x()*cos(_incidence_angle)+node->y()*sin(_incidence_angle));
            Complex H_inc = H*std::exp(j*phase);
            Complex F_inc = -k*k*(1.0/eps-mu)*H_inc;
            Fe(m) += JxW[qp]*F_inc*phi[m][qp];
          }
        }
        PetscUtils::VecAdd(b, Fe, dof_indices);
        PetscUtils::MatAdd(A, Ke, dof_indices);
      }
    }
  }


  // process external absobing boundary
  {
    // Declare a special finite element object for boundary integration.
    AutoPtr<FEBase> fe_face (FEBase::build(dim-1, fe_type));

    // Boundary integration requires one quadraure rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    QGauss qface(dim-1, int_order);

    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule (&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe_face->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe_face->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<Real> >& dphidxi      = fe_face->get_dphidxi();

    // Define data structures to contain the element matrix
    DenseMatrix<Complex> Ke;

    for(unsigned e=0; e<absorb_edge_chain.size(); ++e)
    {
      const Elem * boundary_elem = absorb_edge_chain[e].first;
      unsigned int f = absorb_edge_chain[e].second;

      if(boundary_elem->processor_id()!=Genius::processor_id()) continue;

      Point norm = boundary_elem->outside_unit_normal(f);
      Point unit = (boundary_elem->point(1)-boundary_elem->point(0)).unit();

      AutoPtr<Elem> boundary_face=boundary_elem->build_side(f);

      std::vector<PetscInt> dof_indices;
      this->build_dof_indices(boundary_face.get(), dof_indices);

      fe_face->reinit (boundary_face.get());

      // one complex variable per node
      Ke.resize (boundary_face->n_nodes(), boundary_face->n_nodes());

      std::vector<double> curvatures = curvature_at_edge(e, boundary_face.get());

      for (unsigned int qp=0; qp<qface.n_points(); qp++)
      {
        for (unsigned int m=0; m<phi.size(); m++)
        {

          for (unsigned int n=0; n<phi.size(); n++)
          {
            if(_abc_type == FirstOrder)
              Ke(m,n) += JxW[qp]*(j*k+0.5*curvatures[n])*phi[m][qp]*phi[n][qp];

            if(_abc_type == SecondOrder)
            {
              Complex r1 = j*k + 0.5*curvatures[n] - j*curvatures[n]*curvatures[n]/(8.0*(j*curvatures[n])-k);
              Complex r2 = -j/(2.0*(j*curvatures[n]-k));
              Ke(m,n) += JxW[qp]*(r1*phi[m][qp]*phi[n][qp]);
              Ke(m,n) += (r2*dphidxi[m][qp]*dphidxi[n][qp])/boundary_face->volume();
            }
          }
        }
      }

      PetscUtils::MatAdd(A, Ke, dof_indices);
    }
  }


  // assemble matrix and vec
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);

  //MatView(A, PETSC_VIEWER_DRAW_WORLD);
  //getchar();
  // All done!
#endif
}




void EMFEM2DSolver::build_TM_matrix_rhs(double lambda, double power, double phase0)
{
#if 0
  //wave vector
  double k = 2*M_PI/lambda;

  // wave magnitude, compute from power
  double E = sqrt(power*sqrt(mu0/eps0));

  Complex j(0,1);

  const MeshBase& mesh = _system.mesh();
  const unsigned int dim = mesh.mesh_dimension();
  genius_assert(dim==2);

  // the integral order used in gauss intergral. default value= 2*fe_order+1
  Order int_order=SECOND;
  FEType fe_type;

  VecZeroEntries(x);
  VecZeroEntries(b);
  MatZeroEntries(A);

  // process scatter field
  {
    AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

    // Gauss quadrature rule for numerical integration.
    QGauss qrule (dim, int_order);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<Real> >& dphidx       = fe->get_dphidx();
    const std::vector<std::vector<Real> >& dphidy       = fe->get_dphidy();
    const std::vector<std::vector<Real> >& dphidz       = fe->get_dphidz();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.
    DenseMatrix<Complex> Ke;
    DenseVector<Complex> Fe;

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
    for(unsigned int n=0; n<_system.n_regions(); ++n)
    {
      SimulationRegion * region = _system.region(n);

      Complex r = region->get_optical_refraction(lambda);
      Complex eps(r.real()*r.real()-r.imag()*r.imag(), -2*r.real()*r.imag()) ;
      double mu=1.0;

      //for all the elements in this region
      // note, they are all local element, thus must be processed
      SimulationRegion::element_iterator it = region->elements_begin();
      SimulationRegion::element_iterator it_end = region->elements_end();
      for(; it!=it_end; ++it)
      {
        const Elem* elem  = *it;
        genius_assert(elem->active());
        if(elem->processor_id()!=Genius::processor_id()) continue;

        std::vector<PetscInt> dof_indices;
        this->build_dof_indices(elem, dof_indices);

        fe->reinit (elem);

        // one complex variable per noe
        Ke.resize (elem->n_nodes(), elem->n_nodes());
        Fe.resize (elem->n_nodes());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          for (unsigned int m=0; m<phi.size(); m++)
          {
            for (unsigned int n=0; n<phi.size(); n++)
            {
              // \nabla^2 E^_{sc}
              Ke(m,n) += - JxW[qp]/mu*( dphidx[m][qp]*dphidx[n][qp]
                                        + dphidy[m][qp]*dphidy[n][qp]
                                        + dphidz[m][qp]*dphidz[n][qp]);
              // k^2*eps*E^_{sc}
              Ke(m,n) += JxW[qp]*k*k*eps*phi[m][qp]*phi[n][qp];
            }
            // source item, F_inc = \nabla^2 E^_{inc} + k^2*eps*E^_{inc}
            const Node * node = elem->get_node(m);
            double phase = phase0 - k*(node->x()*cos(_incidence_angle)+node->y()*sin(_incidence_angle));
            Complex E_inc = E*std::exp(j*phase);
            Complex F_inc = -k*k*(1.0/mu-eps)*E_inc;
            Fe(m) += JxW[qp]*F_inc*phi[m][qp];
          }
        }
        PetscUtils::VecAdd(b, Fe, dof_indices);
        PetscUtils::MatAdd(A, Ke, dof_indices);
      }
    }
  }


  // process external absobing boundary
  {
    // Declare a special finite element object for boundary integration.
    AutoPtr<FEBase> fe_face (FEBase::build(dim-1, fe_type));

    // Boundary integration requires one quadraure rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    QGauss qface(dim-1, int_order);

    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule (&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe_face->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe_face->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<Real> >& dphidxi      = fe_face->get_dphidxi();

    // Define data structures to contain the element matrix
    DenseMatrix<Complex> Ke;

    for(unsigned e=0; e<absorb_edge_chain.size(); ++e)
    {
      const Elem * boundary_elem = absorb_edge_chain[e].first;
      unsigned int f = absorb_edge_chain[e].second;

      if(boundary_elem->processor_id()!=Genius::processor_id()) continue;

      Point norm = boundary_elem->outside_unit_normal(f);
      Point unit = (boundary_elem->point(1)-boundary_elem->point(0)).unit();

      AutoPtr<Elem> boundary_face=boundary_elem->build_side(f);

      std::vector<PetscInt> dof_indices;
      this->build_dof_indices(boundary_face.get(), dof_indices);

      fe_face->reinit (boundary_face.get());

      // one complex variable per node
      Ke.resize (boundary_face->n_nodes(), boundary_face->n_nodes());

      std::vector<double> curvatures = curvature_at_edge(e, boundary_face.get());

      for (unsigned int qp=0; qp<qface.n_points(); qp++)
      {
        for (unsigned int m=0; m<phi.size(); m++)
        {

          for (unsigned int n=0; n<phi.size(); n++)
          {
            // (n*k0+curvature/2)*E^_{sc}
            if(_abc_type == FirstOrder)
              Ke(m,n) += JxW[qp]*(j*k+0.5*curvatures[n])*phi[m][qp]*phi[n][qp];

            if(_abc_type == SecondOrder)
            {
              Complex r1 = j*k + 0.5*curvatures[n] - j*curvatures[n]*curvatures[n]/(8.0*(j*curvatures[n])-k);
              Complex r2 = -j/(2.0*(j*curvatures[n]-k));
              Ke(m,n) += JxW[qp]*(r1*phi[m][qp]*phi[n][qp]);
              Ke(m,n) += (r2*dphidxi[m][qp]*dphidxi[n][qp])/boundary_face->volume();
            }
          }
        }
      }

      PetscUtils::MatAdd(A, Ke, dof_indices);
    }
  }


  // assemble matrix and vec
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);

  //MatView(A, PETSC_VIEWER_DRAW_WORLD);
  //getchar();
  // All done!
#endif
}




void EMFEM2DSolver::save_TE_solution(double lambda, double power, double phase0, double eta, bool eta_auto, bool append)
{
#if 0
  //wave vector
  double k = 2*M_PI/lambda;
  double c = 1.0/sqrt(eps0*mu0);
  double omega = 2*M_PI*c/lambda;

  // wave magnitude, compute from power
  double H = sqrt(power*sqrt(eps0/mu0));

  Complex j(0,1);

  // the solution is in vec x
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  VecGetArray(lx, &lxx);

  // calculate Exy from Hz, E=1/(j*omega*eps*eps0)*(\nabla cross Hz)
  std::map<const Node *, Complex> node_Exy_map;
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    Complex r = region->get_optical_refraction(lambda);
    Complex eps(r.real()*r.real()-r.imag()*r.imag(), -2*r.real()*r.imag());

    SimulationRegion::element_iterator elem_it = region->elements_begin();
    SimulationRegion::element_iterator elem_it_end = region->elements_end();
    for(; elem_it!=elem_it_end; ++elem_it)
    {
      const Elem* elem  = *elem_it;
      std::vector<Complex> H_element;
      for(unsigned int i=0; i<elem->n_nodes(); ++i)
      {
        const Node * node = elem->get_node(i);
        double phase = phase0 - k*(node->x()*cos(_incidence_angle)+node->y()*sin(_incidence_angle));
        Complex H_inc = H*std::exp(Complex(0,phase));
        unsigned int local_offset = node->local_dof_id();
        Complex H_field = H_inc - Complex(lxx[local_offset], lxx[local_offset+1]);
        H_element.push_back(H_field);
      }
      VectorValue<Complex> gradH = elem->gradient(H_element);
      Complex Ex = 1.0/(j*eps*eps0*omega)*gradH[1]; //dH/dy
      Complex Ey = 1.0/(j*eps*eps0*omega)*gradH[0]; //dH/dx
      Complex E  = std::sqrt(Ex*Ex+Ey*Ey);

      for(unsigned int i=0; i<elem->n_nodes(); ++i)
      {
        const Node * node = elem->get_node(i);
        Complex E_weight = E*elem->partial_volume_truncated(i)/elem->volume();
        if(node_Exy_map.find(node)!=node_Exy_map.end())
          node_Exy_map[node] += E_weight;
        else
          node_Exy_map[node]  = E_weight;
      }
    }
  }

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    Complex r = region->get_optical_refraction(lambda);
    Complex eps(r.real()*r.real()-r.imag()*r.imag(), -2*r.real()*r.imag()) ;

    // calculate optical gen
    if(eta_auto)
    {
      double Eg=region->get_optical_Eg();
      double E_phono = h*c/lambda;
      eta = floor(E_phono/Eg)*Eg/E_phono;
    }

    SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
    SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = *node_it;

      const Node * node = fvm_node->root_node();
      double phase = phase0 - k*(node->x()*cos(_incidence_angle)+node->y()*sin(_incidence_angle));
      Complex H_inc = H*std::exp(Complex(0,phase));

      unsigned int local_offset = node->local_dof_id();

      FVM_NodeData * fvm_node_data = fvm_node->node_data();
      Complex H_field = H_inc - Complex(lxx[local_offset], lxx[local_offset+1]);
      Complex E_field = node_Exy_map[node];
      if(append)
      {
        fvm_node_data->OptH_complex() += H_field;
        fvm_node_data->OptG()         += eta*M_PI/h*eps0*(-eps.imag())*std::abs(E_field)*std::abs(E_field);
      }
      else
      {
        fvm_node_data->OptH_complex()  = H_field;
        fvm_node_data->OptG()          = eta*M_PI/h*eps0*(-eps.imag())*std::abs(E_field)*std::abs(E_field);
      }
    }
  }

  VecRestoreArray(lx, &lxx);
  
#endif
}



void EMFEM2DSolver::save_TM_solution(double lambda, double power, double phase0, double eta, bool eta_auto, bool append)
{
#if 0
  //wave vector
  double k = 2*M_PI/lambda;
  double c = 1.0/sqrt(eps0*mu0);

  // wave magnitude, compute from power
  double E = sqrt(power*sqrt(mu0/eps0));

  // the solution is in vec x
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  VecGetArray(lx, &lxx);

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    Complex r = region->get_optical_refraction(lambda);
    Complex eps(r.real()*r.real()-r.imag()*r.imag(), -2*r.real()*r.imag()) ;

    // calculate optical gen
    if(eta_auto)
    {
      double Eg=region->get_optical_Eg();
      double E_phono = h*c/lambda;
      eta = floor(E_phono/Eg)*Eg/E_phono;
    }

    SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
    SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = *node_it;

      const Node * node = fvm_node->root_node();
      double phase = phase0 - k*(node->x()*cos(_incidence_angle)+node->y()*sin(_incidence_angle));
      Complex E_inc = E*std::exp(Complex(0,phase));

      unsigned int local_offset = node->local_dof_id();

      FVM_NodeData * fvm_node_data = fvm_node->node_data();
      Complex E_field = E_inc - Complex(lxx[local_offset], lxx[local_offset+1]);

      if(append)
      {
        fvm_node_data->OptE_complex() += E_field;
        fvm_node_data->OptG()         += eta*M_PI/h*eps0*(-eps.imag())*std::abs(E_field)*std::abs(E_field);
      }
      else
      {
        fvm_node_data->OptE_complex()  = E_field;
        fvm_node_data->OptG()          = eta*M_PI/h*eps0*(-eps.imag())*std::abs(E_field)*std::abs(E_field);
      }
    }
  }

  VecRestoreArray(lx, &lxx);
#endif
}




#define ORDER 6
std::vector<double> EMFEM2DSolver::curvature_at_edge(unsigned int i, const Elem * edge)
{
#if 0
  // find the ORDER neighbor points on absorbing boundary
  const Node * p[ORDER];

  p[0] = absorb_node_chain[loop_count(absorb_node_chain.size(), i, -2)].first;
  p[1] = absorb_node_chain[loop_count(absorb_node_chain.size(), i, -1)].first;
  p[2] = absorb_node_chain[i].first;
  p[3] = absorb_node_chain[i].second;
  p[4] = absorb_node_chain[loop_count(absorb_node_chain.size(), i,  1)].second;
  p[5] = absorb_node_chain[loop_count(absorb_node_chain.size(), i,  2)].second;

  // when we know the shape is circle, things are easy
  if(_abc_shape == Circle ) return curvature_of_circle(*p[1], *p[2], *p[3]);

  //or we have to do a interpolation with spline
  bool inverse=false;
  if(edge->get_node(0)==p[3] && edge->get_node(1)==p[2])
    inverse=!inverse;

  //--------------------------------
  // test if we can use y=f(x)
  //--------------------------------

  // we need to have an increasement x

  if(p[ORDER-1]->x() < p[0]->x())
  {
    std::swap(p[0], p[5]);
    std::swap(p[1], p[4]);
    std::swap(p[2], p[3]);
    inverse = !inverse;
  }

  // test if monotone in x
  bool bad_x=false;
  for(unsigned int n=0; n<ORDER-1; ++n)
  {
    if(p[n+1]->x()<=p[n]->x()) bad_x = true;
  }

  //--------------------------------
  // or we have to use x=f(y)
  //--------------------------------
  bool bad_y=false;
  if(bad_x)
  {
    if(p[ORDER-1]->y() < p[0]->y())
    {
      std::swap(p[0], p[5]);
      std::swap(p[1], p[4]);
      std::swap(p[2], p[3]);
      inverse = !inverse;
    }

    for(unsigned int n=0; n<ORDER-1; ++n)
    {
      if(p[n+1]->y()<=p[n]->y()) bad_y = true;
    }
  }

  // non of them can be used? too had
  if(bad_x && bad_y)
  {
    MESSAGE<<"ERROR at " <<_card.get_fileline()<< " EMFEM2D: Absorbing boundary is not smooth enough for curvature operation." << std::endl; RECORD();
    genius_error();
  }


  double t[ORDER];
  double y[ORDER];
  double *ypp;
  double t1var, t2var;
  if(!bad_x)
  {
    for(unsigned int n=0; n<ORDER; ++n)
    {
      t[n] = p[n]->x();
      y[n] = p[n]->y();
    }

    t1var = p[2]->x();
    t2var = p[3]->x();
    if(inverse) std::swap(t1var, t2var);
  }
  else
  {
    genius_assert(bad_y==false);
    for(unsigned int n=0; n<ORDER; ++n)
    {
      t[n] = p[n]->y();
      y[n] = p[n]->x();
    }

    t1var = p[2]->y();
    t2var = p[3]->y();
    if(inverse) std::swap(t1var, t2var);
  }

  ypp = SPLINE::spline_cubic_set ( ORDER, t, y, 2,  0.0, 2, 0.0 );

  double t1varp, t1varpp;
  SPLINE::spline_cubic_val ( ORDER, t, t1var, y, ypp, &t1varp, &t1varpp );
  double curvature1 = std::abs(t1varpp)/std::pow(1+t1varp*t1varp, 1.5);

  double t2varp, t2varpp;
  SPLINE::spline_cubic_val ( ORDER, t, t2var, y, ypp, &t2varp, &t2varpp );
  double curvature2 = std::abs(t2varpp)/std::pow(1+t2varp*t2varp, 1.5);

  delete [] ypp;

  std::vector<double> curvatures(2);
  curvatures[0]=curvature1;
  curvatures[1]=curvature2;

  return curvatures;
#endif
}



std::vector<double> EMFEM2DSolver::curvature_of_circle(const Point &p0, const Point &p1, const Point &p2)
{
#if 0
  // Vector pointing from A to C
  Point AC ( p2 - p0 );

  // Vector pointing from A to B
  Point AB ( p1 - p0 );

  // Vector pointing from B to C
  Point BC ( p2 - p1 );

  // Vector normal to plan ABC
  Point n = AB.cross(AC);

  // Vector in plan ABC and normal to AB
  Point a = n.cross(AB);

  // the angle between AB and AC
  double angle_a = AB.angle(AC);

  // the radius of circumcircle
  double R = BC.size() / (2*sin(angle_a));

  std::vector<double> curvatures(2);
  curvatures[0]=1.0/R;
  curvatures[1]=1.0/R;

  return curvatures;
#endif
}


void EMFEM2DSolver::build_absorb_chain()
{
#if 0
  const MeshBase & mesh = system.mesh();
  const BoundaryConditionCollector  * bcs = _system.get_bcs();

  std::vector< std::pair<const Elem *, unsigned int> > edge_chain;

  // collect all the absorb boundary
  for(unsigned b=0; b<bcs->n_bcs(); ++b)
  {
    const BoundaryCondition * bc = bcs->get_bc(b);
    if(bc->bc_type()!=AbsorbingBoundary) continue;

    for(unsigned n=0; n<bc->n_elems(); ++n)
      edge_chain.push_back(bc->get_elem(n));
  }
  if(!edge_chain.size()) return;

  // reorder
  std::vector<bool> visit_flag(edge_chain.size(), false);

  absorb_edge_chain.push_back(edge_chain[0]);
  visit_flag[0]=true;
  AutoPtr<Elem> first_elem = edge_chain[0].first->build_side(edge_chain[0].second);
  const Node * first_point = first_elem->get_node(0);
  const Node * current_point = first_elem->get_node(1);
  absorb_node_chain.push_back(std::make_pair(first_point, current_point));

  for(unsigned int m=1; m<edge_chain.size(); ++m) //loop edge_chain.size() -1 times
    for(unsigned int n=1; n<edge_chain.size(); ++n)
    {
      if(!visit_flag[n]) //not visited
      {
        const Elem * boundary_elem = edge_chain[n].first;
        unsigned int f = edge_chain[n].second;
        AutoPtr<Elem> boundary_face = boundary_elem->build_side(f);
        if(boundary_face->get_node(0) == current_point)
        {
          visit_flag[n] = true;
          current_point = boundary_face->get_node(1);
          absorb_edge_chain.push_back(edge_chain[n]);
          absorb_node_chain.push_back(std::make_pair(boundary_face->get_node(0), boundary_face->get_node(1)));
          break;
        }
        if(boundary_face->get_node(1) == current_point)
        {
          visit_flag[n] = true;
          current_point = boundary_face->get_node(0);
          absorb_edge_chain.push_back(edge_chain[n]);
          absorb_node_chain.push_back(std::make_pair(boundary_face->get_node(1), boundary_face->get_node(0)));
          break;
        }
      }
    }

  genius_assert(absorb_edge_chain.size()==edge_chain.size());
  genius_assert(absorb_node_chain[0].first==absorb_node_chain[absorb_node_chain.size()-1].second);
#endif
}




void EMFEM2DSolver::setup_solver_parameters()
{
#if 0
  // optical wave is defined by command line
  if(_card.is_parameter_exist("lambda")||_card.is_parameter_exist("wavelength"))
  {
    OpticalSource source;
    source.wave_length = _card.get_real("lambda", 0.532, "wavelength")*um;// wave length of light source
    source.power       = _card.get_real("intensity", 0.0)*J/s/cm/cm;      // incident power
    source.TE_weight   = _card.get_real("wte", 0.5);
    source.TM_weight   = _card.get_real("wtm", 0.5);
    source.phi_TM      = _card.get_real("phase.tm", 0.0)*M_PI/180;
    source.phi_TE      = _card.get_real("phase.te", 0.0)*M_PI/180;
    source.eta_auto    = !_card.is_parameter_exist("quan.eff");
    source.eta         = _card.get_real("quan.eff", 1.0);

    if(std::abs(source.TE_weight+source.TM_weight-1.0)>1e-3)
    {
      MESSAGE<<"ERROR at " <<_card.get_fileline()<< " EMFEM2D: Total weight of TE and TM is greater than 1.0." << std::endl; RECORD();
      genius_error();
    }

    _optical_sources.push_back(source);
  }



  // read optical wave from spectrum file
  if(_card.is_parameter_exist("spectrumfile"))
  {
    parse_spectrum_file(_card.get_string("spectrumfile", ""));
  }

  _abc_type = SecondOrder;
  if(_card.is_parameter_exist("abc.type"))
  {
    if( _card.is_enum_value("abc.type","firstorder"))
      _abc_type = FirstOrder;
    if( _card.is_enum_value("abc.type","secondorder"))
      _abc_type = SecondOrder;
    if( _card.is_enum_value("abc.type","pml"))
      _abc_type = PML;
  }

  _abc_shape = UNKNOWN_SHAPE;
  if(_card.is_parameter_exist("abc.shape"))
  {
    _card.is_enum_value("abc.shape", "circle");
    _abc_shape = Circle;
  }

  // the angle of incident wave
  _incidence_angle = _card.get_real("angle", 90.0)*M_PI/180;

  // set linear solver type
  SolverSpecify::LS = SolverSpecify::linear_solver_type(_card.get_string("ls", "bcgs"));

  // set preconditioner type
  SolverSpecify::PC = SolverSpecify::preconditioner_type(_card.get_string("pc", "asm"));

  // build the absorbing boundary
  build_absorb_chain();
#endif
}



void EMFEM2DSolver::parse_spectrum_file(const std::string & filename)
{
#if 0
  // only processor 0 read the spectrum file
  std::vector<OpticalSource> _opt_srcs;
  if(Genius::processor_id() == 0)
  {
    std::ifstream in(filename.c_str());
    genius_assert(in.good());

    std::string line;
    while(!in.eof())
    {
      std::getline(in, line);
      // drop blank/table/cr chars
      {
        std::string filt_elems(" \t\r\n");
        std::string::size_type pos = 0;
        while (( pos = line.find_first_of( filt_elems, pos )) != std::string::npos )
        {
          line.erase(pos, 1);
        }
      }
      // skip empty lines
      if(line.size()==0) continue;
      // skip the line begin with '#'
      if(line.find('#')==0) continue;

      std::stringstream ss;
      ss << line;

      OpticalSource source;
      // read wave length and power information from line
      ss >> source.wave_length;
      ss >> source.power;
      // test if it contains extra parameter for quantum efficiency
      ss >> source.eta;
      // if no quantum efficiency parameter find, set eta_auto flag to true
      source.eta_auto = ss.fail();

      _opt_srcs.push_back(source);
    }
    in.close();
  }

  // broadcast to all the processors
  unsigned int n_sources = _opt_srcs.size();
  if (Genius::n_processors() > 1)
  {
    Parallel::broadcast(n_sources);
    if (Genius::processor_id() != 0)
      _opt_srcs.resize(n_sources);
    for(unsigned int n=0; n<n_sources; ++n)
    {
      std::vector<double> package(4);
      if (Genius::processor_id() == 0)
      {
        package[0] = _opt_srcs[n].wave_length;
        package[1] = _opt_srcs[n].power;
        package[2] = double(_opt_srcs[n].eta_auto);
        package[3] = _opt_srcs[n].eta;
      }

      Parallel::broadcast(package);

      if (Genius::processor_id() != 0)
      {
        _opt_srcs[n].wave_length = package[0];
        _opt_srcs[n].power       = package[1];
        _opt_srcs[n].eta_auto    = package[2]>0.5;
        _opt_srcs[n].eta         = package[3];
      }
    }
  }

  //post process
  genius_assert(n_sources>=2);

  //this three parameters are still defined by command line
  double TE_weight   = _card.get_real("wte", 0.5);
  double TM_weight   = _card.get_real("wtm", 0.5);
  double phi_TM    = _card.get_real("phase.tm", 0.0)*M_PI/180;
  double phi_TE    = _card.get_real("phase.te", 0.0)*M_PI/180;

  if(std::abs(TE_weight+TM_weight-1.0)>1e-3)
  {
    MESSAGE<<"ERROR at " <<_card.get_fileline()<< " EMFEM2D: Total weight of TE and TM is greater than 1.0." << std::endl; RECORD();
    genius_error();
  }

  double total_power = 0.0;
  for(unsigned int n=0; n<n_sources; ++n)
  {
    double lambda_distance;
    if(n==0)
      lambda_distance = 0.5*(_opt_srcs[1].wave_length -_opt_srcs[0].wave_length);
    else if(n==n_sources-1)
      lambda_distance = 0.5*(_opt_srcs[n_sources-1].wave_length -_opt_srcs[n_sources-2].wave_length);
    else
      lambda_distance = 0.5*(_opt_srcs[n+1].wave_length -_opt_srcs[n-1].wave_length);

    OpticalSource source;
    source.wave_length = _opt_srcs[n].wave_length*um;
    source.power = (_opt_srcs[n].power*lambda_distance)*J/s/(cm*cm);
    source.TE_weight   = TE_weight;
    source.TM_weight   = TM_weight;
    source.phi_TM      = phi_TM;
    source.phi_TE      = phi_TE;
    source.eta_auto    = _opt_srcs[n].eta_auto;
    source.eta         = _opt_srcs[n].eta;
    _optical_sources.push_back(source);

    total_power += _opt_srcs[n].power*lambda_distance;
  }

  MESSAGE<<"Read "<<n_sources<<" spectrums from file "<<filename<<".\n"
  <<"Total power intensity is approximate "<<std::fixed<<total_power<<" W/(cm^2)"
  <<std::endl;

  RECORD();
#endif

}
