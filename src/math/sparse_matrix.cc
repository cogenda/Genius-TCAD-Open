// $Id: sparse_matrix.C 2866 2008-06-16 01:11:52Z benkirk $

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



// C++ includes
#include <numeric>


// Local Includes
#include "sparse_matrix.h"
#include "petsc_matrix.h"
#include "symbolic_matrix.h"
#include "parallel.h"

//------------------------------------------------------------------
// SparseMatrix Methods

template <typename T>
SparseMatrix<T>::SparseMatrix (const unsigned int m,   const unsigned int n,
                               const unsigned int m_l, const unsigned int n_l)
  :_m_global(m), _n_global(n), _m_local(m_l), _n_local(n_l), _is_initialized(false)
{
  std::vector<int> block_size;
  block_size.push_back(_m_local);
  Parallel::allgather(block_size);

  // the offset of local block at global dof array
  _global_offset = std::accumulate(block_size.begin(), block_size.begin()+Genius::processor_id(), 0 );
}



template <typename T>
SparseMatrix<T>::~SparseMatrix ()
{}



// Full specialization for Real datatypes
template <typename T>
AutoPtr<SparseMatrix<T> >
SparseMatrix<T>::build(const std::string & solver_package,
                       const unsigned int m,   const unsigned int n,
                       const unsigned int m_l, const unsigned int n_l,
                       const std::vector<int> &n_nz,
                       const std::vector<int> &n_oz)
{
  // Build the appropriate vector
  if (solver_package == "petsc")
  {
#ifdef HAVE_PETSC
    AutoPtr<SparseMatrix<T> > ap(new PetscMatrix<T>(m, n, m_l, n_l, n_nz, n_oz));
    return ap;
#endif
  }

  if (solver_package == "symbolic")
  {
    AutoPtr<SparseMatrix<T> > ap(new SymbolicMatrix<T>(m, n, m_l, n_l));
    return ap;
  }

  AutoPtr<SparseMatrix<T> > ap(NULL);
  return ap;
}


template <typename T>
void SparseMatrix<T>::print(std::ostream& os) const
{
  genius_assert (this->initialized());

  for (unsigned int i=0; i<this->m(); i++)
  {
    for (unsigned int j=0; j<this->n(); j++)
      os << std::setw(8) << (*this)(i,j) << " ";
    os << std::endl;
  }
}



// Full specialization for Complex datatypes
template <>
void SparseMatrix<Complex>::print(std::ostream& os) const
{
  // std::complex<>::operator<<() is defined, but use this form

  std::cout << "Real part:" << std::endl;
  for (unsigned int i=0; i<this->m(); i++)
  {
    for (unsigned int j=0; j<this->n(); j++)
      os << std::setw(8) << (*this)(i,j).real() << " ";
    os << std::endl;
  }

  os << std::endl << "Imaginary part:" << std::endl;
  for (unsigned int i=0; i<this->m(); i++)
  {
    for (unsigned int j=0; j<this->n(); j++)
      os << std::setw(8) << (*this)(i,j).imag() << " ";
    os << std::endl;
  }
}


//------------------------------------------------------------------
// Explicit instantiations
template class SparseMatrix<Real>;
