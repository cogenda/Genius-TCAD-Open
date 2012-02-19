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

//  $Id: poisson.cc,v 1.36 2008/07/09 07:53:36 gdiso Exp $


#include "elem.h"
#include "mesh_base.h"
#include "stress_solver/stress_solver.h"
#include "solver_specify.h"
#include "fe_type.h"
#include "fe_base.h"
#include "quadrature_gauss.h"
#include "tensor_value.h"


/*------------------------------------------------------------------
 * create the poisson solver contex
 */
int StressSolver::create_solver()
{
  MESSAGE<< '\n' << "Stress Solver init..." << std::endl;
  RECORD();


  // adjust default snes/ksp/pc/mat/vec settings

  // use BCGS as default solver
  KSPSetType(ksp, KSPBCGS);

  // Set user-specified  solver and preconditioner types
  set_linear_solver_type    ( SolverSpecify::LS );
  set_preconditioner_type   ( SolverSpecify::PC );

  // must set linear matrix/vector here!
  setup_linear_data();

  //SNESLineSearchSet(snes,SNESLineSearchNo, 0);

  // rtol   = 1e-10*n_global_dofs  - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol = 1e-20*n_global_dofs  - the absolute convergence tolerance (absolute size of the residual norm)
  KSPSetTolerances(ksp, 1e-10*n_global_dofs, 1e-20*n_global_dofs, PETSC_DEFAULT, n_global_dofs/10);

  // user can do further adjusment from command line
  KSPSetFromOptions (ksp);

  return 0;

}





/*------------------------------------------------------------------
 *
 */
int StressSolver::solve()
{
  START_LOG("StressSolver_Linear()", "StressSolver");

  build_matrix(A, A);

  build_rhs(b);

  KSPSolve(ksp,b,x);

  STOP_LOG("StressSolver_Linear()", "StressSolver");

  return 0;
}

int StressSolver::destroy_solver()
{

  // clear linear contex
  clear_linear_data();

  return 0;
}



void StressSolver::build_matrix(Mat A, Mat )
{
  // the integral order used in gauss intergral. default value= 2*fe_order+1
  Order int_order=FIFTH;

  const MeshBase& mesh = _system.mesh();

  const unsigned int dim = mesh.mesh_dimension();

  FEType fe_type;

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  // Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, int_order);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for boundary integration.
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));

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
  const std::vector<Real>& JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  //const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<std::vector<Real> >& dphidx       = fe->get_dphidx();
  const std::vector<std::vector<Real> >& dphidy       = fe->get_dphidy();
  const std::vector<std::vector<Real> >& dphidz       = fe->get_dphidz();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".  These datatypes are templated on
  //  Number, which allows the same code to work for real
  // or complex numbers.
  //DenseMatrix<Number> Ke;
  //DenseVector<Number> Fe;


  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  //std::vector<unsigned int> dof_indices;

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element iterators are a nice way to iterate through
  // all the elements, or all the elements that have some property.
  // There are many types of element iterators, but here we will
  // use the most basic type, the  const_elem_iterator.  The iterator
  //  el will iterate from the first to the last element.  The
  // iterator  end_el tells us when to stop.  It is smart to make
  // this one  const so that we don't accidentally mess it up!
  //   const_elem_iterator           el (mesh.elements_begin());
  //   const const_elem_iterator end_el (mesh.elements_end());

  //for(unsigned int n=0; n<_system.n_regions(); ++n)
  //{
  //	SimulationRegion * region = _system.region(n);
  //}

  MeshBase::const_element_iterator       el     = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();

  //this should be function of node. for convinent, we take parameter as const in every element
  // you should re_value it for your own use.
  //VectorValue<double> E;      //yang's modul
  //TensorValue<double> G;   //shear modul
  //TensorValue<double> nu;  //possion ratio

  //the stiff matrix. a general stiff matrix named C of 3D is a 6x6 symmetry matrx
  // so it has 21 independent parameter. in some special case the independent parameter number will be reduce.
  // note that in 2D case, plane stress and plane strain are two deferent case.
  double C[6][6]={0.};

  //std::cout<<"k="<<k<<std::endl;
  // Loop over the elements.  Note that  ++el is preferred to
  // el++ since the latter requires an unnecessary temporary
  // object.
  for ( ; el != end_el ; ++el)
  {
    //material initialization for every element, now only a sample is taken
    C[0][0]=C[1][1]=C[2][2]=1.e9;
    C[3][3]=C[4][4]=C[5][5]=1.e8;

    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem* elem = *el;

    //matrix K and F are matrix of Ax=b in every element Kx=F, its size are changed with node in element,
    // so initialize it in every element.
    unsigned int n_node=elem->n_nodes();
    unsigned int dim=3;   // 3 is 3D, 2D and 1D are regard as reduced 3D problem.
    unsigned int dim_stress=6; //6 is the dim of stress. stress is 3x3 symmetry matrix and has 6 independent element.

    unsigned int n_node2=n_node*dim;

    double* K  =new double[n_node2*n_node2];   //K=B'CB
    double* F  =new double[n_node2];
    double* CB =new double[dim_stress*n_node2]; //a matrix=C*B

    for(unsigned int ii=0;ii<n_node2*n_node2;++ii)
    {
      K[ii]=0.;
    }

    for(unsigned int ii=0;ii<n_node2;++ii)
    {
      F[ii]=0.;
    }


    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    //	dof_map.dof_indices (elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit (elem);

    double B[6][81]={0.}; //81=27*3, and 27 is the max n_node in an element.


    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).

    // The  DenseMatrix::resize() and the  DenseVector::resize()
    // members will automatically zero out the matrix  and vector.
    //	Ke.resize (dof_indices.size(),
    //			dof_indices.size());

    //	Fe.resize (dof_indices.size());

    // Now loop over the quadrature points.  This handles
    // the numeric integration.

    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {

      // Now we will build the element matrix.  This involves
      // a double loop to integrate the test funcions (i) against
      // the trial functions (j).
      for (unsigned int i=0; i<phi.size(); i++)
      {
        for (unsigned int j=0; j<phi.size(); j++)
        {
          //	Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
        }
      }

      for(unsigned int i=0;i<phi.size();i++)
      {
        B[0][i*3+0]=dphidx[i][qp];  //vd0[0][i];

        B[1][i*3+1]=dphidy[i][qp];  //vd0[1][i];

        B[2][i*3+2]=dphidz[i][qp];  //vd0[2][i];

        B[3][i*3+0]=dphidy[i][qp];  //vd0[1][i];
        B[3][i*3+1]=dphidx[i][qp];  //vd0[0][i];

        B[4][i*3+1]=dphidz[i][qp];  //vd0[2][i];
        B[4][i*3+2]=dphidy[i][qp];  //vd0[1][i];

        B[5][i*3+0]=dphidz[i][qp];  //vd0[2][i];
        B[5][i*3+2]=dphidx[i][qp];  //vd0[0][i];
      }

      for(unsigned int ii=0;ii<dim_stress;ii++)
      {
        for(unsigned int jj=0;jj<n_node2;jj++)
        {
          CB[ii*n_node2+jj]=0.0;
          for(unsigned int kk=0;kk<dim_stress;kk++)
          {
            CB[ii*n_node2+jj]+=C[ii][kk]*B[kk][jj];
          }
        }
      }

      for(unsigned int ii=0;ii<n_node2;ii++)
      {
        for(unsigned int jj=0;jj<n_node2;jj++)
        {
          for(unsigned int kk=0;kk<dim_stress;kk++)
          {
            K[ii*n_node2+jj]+=B[kk][ii]*CB[kk*n_node2+jj]*JxW[qp];
          }
          //	cout<<element_k[ii]<<endl;
        }
      }

      double f[3]={1.e3,0.,0.};
      for (unsigned int ii=0; ii<n_node; ii++)
      {
        for(unsigned int jj=0; jj<dim; jj++)
        {
          F[ii*dim+jj] += JxW[qp]*f[jj]*phi[ii][qp];
        }
      }
    }

    // We have now reached the end of the RHS summation,
    // and the end of quadrature point loop, so
    // the interior element integration has
    // been completed.  However, we have not yet addressed
    // boundary conditions.  For this example we will only
    // consider simple Dirichlet boundary conditions.
    //
    // There are several ways Dirichlet boundary conditions
    // can be imposed.  A simple approach, which works for
    // interpolary bases like the standard Lagrange polynomials,
    // is to assign function values to the
    // degrees of freedom living on the domain boundary. This
    // works well for interpolary bases, but is more difficult
    // when non-interpolary (e.g Legendre or Hierarchic) bases
    // are used.
    //
    // Dirichlet boundary conditions can also be imposed with a
    // "penalty" method.  In this case essentially the L2 projection
    // of the boundary values are added to the matrix. The
    // projection is multiplied by some large factor so that, in
    // floating point arithmetic, the existing (smaller) entries
    // in the matrix and right-hand-side are effectively ignored.
    //
    // This amounts to adding a term of the form (in latex notation)
    //
    // \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i
    //
    // where
    //
    // \frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1
    /*
    {

    	// The following loop is over the sides of the element.
    	// If the element has no neighbor on a side then that
    	// side MUST live on a boundary of the domain.
    	for (unsigned int side=0; side<elem->n_sides(); side++)
    	{
    		if (elem->neighbor(side) == NULL)
    		{
    			// The value of the shape functions at the quadrature
    			// points.
    			const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();

    			// The Jacobian * Quadrature Weight at the quadrature
    			// points on the face.
    			const std::vector<Real>& JxW_face = fe_face->get_JxW();

    			// The XYZ locations (in physical space) of the
    			// quadrature points on the face.  This is where
    			// we will interpolate the boundary value function.
    			const std::vector<Point >& qface_point = fe_face->get_xyz();

    			// Compute the shape function values on the element
    			// face.
    			fe_face->reinit(elem, side);

    			// Loop over the face quadrature points for integration.
    			for (unsigned int qp=0; qp<qface.n_points(); qp++)
    			{

    				// The location on the boundary of the current
    				// face quadrature point.
    				const Real xf = qface_point[qp](0);
    				const Real yf = qface_point[qp](1);

    				// The penalty value.  \frac{1}{\epsilon}
    				// in the discussion above.
    				const Real penalty = 1.e10;

    				// The boundary value.
    				const Real value = 20.;//exact_solution(xf, yf);

    				// Matrix contribution of the L2 projection.
    				for (unsigned int i=0; i<phi_face.size(); i++)
    				for (unsigned int j=0; j<phi_face.size(); j++)
    				Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

    				// Right-hand-side contribution of the L2
    				// projection.
    				for (unsigned int i=0; i<phi_face.size(); i++)
    				Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
    			}
    		}
    	}
    }/*/

    // We have now finished the quadrature point loop,
    // and have therefore applied all the boundary conditions.
    //
    // The element matrix and right-hand-side are now built
    // for this element.  Add them to the global matrix and
    // right-hand-side vector.  The  SparseMatrix::add_matrix()
    // and  NumericVector::add_vector() members do this for us.
    //	system.matrix->add_matrix (Ke, dof_indices);
    //	system.rhs->add_vector    (Fe, dof_indices);

    delete []K;
    delete []F;
    delete []CB;

  }

  // All done!
}

void StressSolver::build_rhs(Vec b)
{
  //std::cout<<"here b"<<std::endl;
}

