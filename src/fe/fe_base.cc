// $Id: fe_base.C 2674 2008-02-17 20:41:09Z benkirk $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local includes
#include "fe.h"
#include "perf_log.h"
// For projection code:
//#include "boundary_info.h"
//#include "mesh_base.h"
//#include "dense_matrix.h"
//#include "dense_vector.h"
//#include "dof_map.h"
//#include "elem.h"
//#include "fe_interface.h"
//#include "numeric_vector.h"
//#include "quadrature.h"
//#include "quadrature_gauss.h"
//#include "threads.h"



// ------------------------------------------------------------
// FEBase class members
AutoPtr<FEBase> FEBase::build (const unsigned int dim,
                               const FEType& fet)
{
  // The stupid AutoPtr<FEBase> ap(); return ap;
  // construct is required to satisfy IBM's xlC

  switch (dim)
  {
    // 1D
  case 1:
    {
      switch (fet.family)
      {
      case LAGRANGE:
        {
          AutoPtr<FEBase> ap(new FE<1,LAGRANGE>(fet));
          return ap;
        }

      default:
        std::cout << "ERROR: Bad FEType.family= " << fet.family << std::endl;
        genius_error();
      }
    }


    // 2D
  case 2:
    {
      switch (fet.family)
      {
      case LAGRANGE:
        {
          AutoPtr<FEBase> ap(new FE<2,LAGRANGE>(fet));
          return ap;
        }

      default:
        std::cout << "ERROR: Bad FEType.family= " << fet.family << std::endl;
        genius_error();
      }
    }


    // 3D
  case 3:
    {
      switch (fet.family)
      {
      case LAGRANGE:
        {
          AutoPtr<FEBase> ap(new FE<3,LAGRANGE>(fet));
          return ap;
        }

      default:
        std::cout << "ERROR: Bad FEType.family= " << fet.family << std::endl;
        genius_error();
      }
    }

  default:
    genius_error();
  }

  genius_error();
  AutoPtr<FEBase> ap(NULL);
  return ap;
}








void FEBase::compute_shape_functions (const Elem*)
{
  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Start logging the shape function computation
  START_LOG("compute_shape_functions()", "FE");

  calculations_started = true;

  calculate_phi = calculate_dphi = true;

  // Compute the value of the derivative shape function i at quadrature point p
  switch (dim)
  {

  case 1:
    {
      if (calculate_dphi)
        for (unsigned int i=0; i<dphi.size(); i++)
          for (unsigned int p=0; p<dphi[i].size(); p++)
          {
            // dphi/dx    = (dphi/dxi)*(dxi/dx)
            dphi[i][p](0) = dphidx[i][p] = dphidxi[i][p]*dxidx_map[p];

            dphi[i][p](1) = dphidy[i][p] = 0.;
            dphi[i][p](2) = dphidz[i][p] = 0.;
          }


      // All done
      break;
    }

  case 2:
    {
      if (calculate_dphi)
        for (unsigned int i=0; i<dphi.size(); i++)
          for (unsigned int p=0; p<dphi[i].size(); p++)
          {
            // dphi/dx    = (dphi/dxi)*(dxi/dx) + (dphi/deta)*(deta/dx)
            dphi[i][p](0) =
              dphidx[i][p] = (dphidxi[i][p]*dxidx_map[p] +
                              dphideta[i][p]*detadx_map[p]);

            // dphi/dy    = (dphi/dxi)*(dxi/dy) + (dphi/deta)*(deta/dy)
            dphi[i][p](1) =
              dphidy[i][p] = (dphidxi[i][p]*dxidy_map[p] +
                              dphideta[i][p]*detady_map[p]);

            // dphi/dz    = (dphi/dxi)*(dxi/dz) + (dphi/deta)*(deta/dz)
#if DIM == 3
            dphi[i][p](2) = // can only assign to the Z component if DIM==3
#endif
              dphidz[i][p] = (dphidxi[i][p]*dxidz_map[p] +
                              dphideta[i][p]*detadz_map[p]);
          }

      // All done
      break;
    }

  case 3:
    {
      if (calculate_dphi)
        for (unsigned int i=0; i<dphi.size(); i++)
          for (unsigned int p=0; p<dphi[i].size(); p++)
          {
            // dphi/dx    = (dphi/dxi)*(dxi/dx) + (dphi/deta)*(deta/dx) + (dphi/dzeta)*(dzeta/dx);
            dphi[i][p](0) =
              dphidx[i][p] = (dphidxi[i][p]*dxidx_map[p] +
                              dphideta[i][p]*detadx_map[p] +
                              dphidzeta[i][p]*dzetadx_map[p]);

            // dphi/dy    = (dphi/dxi)*(dxi/dy) + (dphi/deta)*(deta/dy) + (dphi/dzeta)*(dzeta/dy);
            dphi[i][p](1) =
              dphidy[i][p] = (dphidxi[i][p]*dxidy_map[p] +
                              dphideta[i][p]*detady_map[p] +
                              dphidzeta[i][p]*dzetady_map[p]);

            // dphi/dz    = (dphi/dxi)*(dxi/dz) + (dphi/deta)*(deta/dz) + (dphi/dzeta)*(dzeta/dz);
            dphi[i][p](2) =
              dphidz[i][p] = (dphidxi[i][p]*dxidz_map[p] +
                              dphideta[i][p]*detadz_map[p] +
                              dphidzeta[i][p]*dzetadz_map[p]);
          }


      // All done
      break;
    }

  default:
    {
      genius_error();
    }
  }

  // Stop logging the shape function computation
  STOP_LOG("compute_shape_functions()", "FE");
}




/*
void FEBase::print_JxW(std::ostream& os) const
{
  for (unsigned int i=0; i<JxW.size(); ++i)
    os << JxW[i] << std::endl;
}




void FEBase::print_phi(std::ostream& os) const
{
  for (unsigned int i=0; i<phi.size(); ++i)
    for (unsigned int j=0; j<phi[i].size(); ++j)
      os << " phi[" << i << "][" << j << "]=" << phi[i][j] << std::endl;
}




void FEBase::print_dphi(std::ostream& os) const
{
  for (unsigned int i=0; i<dphi.size(); ++i)
    for (unsigned int j=0; j<dphi[i].size(); ++j)
      os << " dphi[" << i << "][" << j << "]=" << dphi[i][j];
}





void FEBase::print_xyz(std::ostream& os) const
{
  for (unsigned int i=0; i<xyz.size(); ++i)
    os << xyz[i];
}




void FEBase::print_info(std::ostream& os) const
{
  os << "Shape functions at the Gauss pts." << std::endl;
  this->print_phi(os);

  os << "Shape function gradients at the Gauss pts." << std::endl;
  this->print_dphi(os);

  os << "XYZ locations of the Gauss pts." << std::endl;
  this->print_xyz(os);

  os << "Values of JxW at the Gauss pts." << std::endl;
  this->print_JxW(os);
}




std::ostream& operator << (std::ostream& os, const FEBase& fe)
{
  fe.print_info(os);
  return os;
}*/
