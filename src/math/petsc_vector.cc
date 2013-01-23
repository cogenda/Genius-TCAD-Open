// $Id: petsc_vector.C 2789 2008-04-13 02:24:40Z roystgnr $

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

// Local Includes
#include "petsc_vector.h"
#include "petsc_matrix.h"

#ifdef HAVE_PETSC

#include "parallel.h"
#include "dense_vector.h"
#include "petsc_macro.h"



//-----------------------------------------------------------------------
// PetscVector members

// void PetscVector<T>::init (const NumericVector<T>& v, const bool fast)
// {
//   libmesh_error();

//   init (v.local_size(), v.size(), fast);

//   vec = dynamic_cast<const PetscVector<T>&>(v).vec;
// }

template <typename T>
T PetscVector<T>::sum () const
{
  genius_assert(this->closed());

  int ierr=0;
  PetscScalar value=0.;

  ierr = VecSum (_vec, &value);

  return static_cast<T>(value);
}


template <typename T>
Real PetscVector<T>::l1_norm () const
{
  genius_assert(this->closed());

  int ierr=0;
  PetscReal value=0.;

  ierr = VecNorm (_vec, NORM_1, &value);

  return static_cast<Real>(value);
}



template <typename T>
Real PetscVector<T>::l2_norm () const
{
  genius_assert(this->closed());

  int ierr=0;
  PetscReal value=0.;

  ierr = VecNorm (_vec, NORM_2, &value);

  return static_cast<Real>(value);
}




template <typename T>
Real PetscVector<T>::linfty_norm () const
{
  genius_assert(this->closed());

  int ierr=0;
  PetscReal value=0.;

  ierr = VecNorm (_vec, NORM_INFINITY, &value);

  return static_cast<Real>(value);
}




template <typename T>
NumericVector<T>&
PetscVector<T>::operator += (const NumericVector<T>& v)
{
  genius_assert(this->closed());

  this->add(1., v);

  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator -= (const NumericVector<T>& v)
{
  genius_assert(this->closed());

  this->add(-1., v);

  return *this;
}



template <typename T>
void PetscVector<T>::set (const unsigned int i, const T value)
{
  genius_assert(i<size());

  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, INSERT_VALUES);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add (const unsigned int i, const T value)
{
  genius_assert(i<size());

  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, ADD_VALUES);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add_vector (const std::vector<T>& v,
                                 const std::vector<unsigned int>& dof_indices)
{
  genius_assert (v.size() == dof_indices.size());

  for (unsigned int i=0; i<v.size(); i++)
    this->add (dof_indices[i], v[i]);
}



template <typename T>
void PetscVector<T>::add_vector (const NumericVector<T>& V,
                                 const std::vector<unsigned int>& dof_indices)
{
  genius_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::add_vector (const NumericVector<T>& V_in,
                                 const SparseMatrix<T>& A_in)
{
  const PetscVector<T>* V = dynamic_cast<const PetscVector<T>*>(&V_in);
  const PetscMatrix<T>* A = dynamic_cast<const PetscMatrix<T>*>(&A_in);

  genius_assert (V != NULL);
  genius_assert (A != NULL);

  int ierr=0;

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.
  ierr = MatMultAdd(const_cast<PetscMatrix<T>*>(A)->mat(), V->_vec, _vec, _vec);
}



template <typename T>
void PetscVector<T>::add_vector (const DenseVector<T>& V,
                                 const std::vector<unsigned int>& dof_indices)
{
  genius_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::add (const T v_in)
{
  int ierr=0;
  PetscScalar* values;
  const PetscScalar v = static_cast<PetscScalar>(v_in);
  const int n   = static_cast<int>(this->local_size());
  const int fli = static_cast<int>(this->first_local_index());

  for (int i=0; i<n; i++)
  {
    ierr = VecGetArray (_vec, &values);

    int ig = fli + i;

    PetscScalar value = (values[ig] + v);

    ierr = VecRestoreArray (_vec, &values);

    ierr = VecSetValues (_vec, 1, &ig, &value, INSERT_VALUES);
  }

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add (const NumericVector<T>& v)
{
  this->add (1., v);
}



template <typename T>
void PetscVector<T>::add (const T a_in, const NumericVector<T>& v_in)
{
  int ierr = 0;
  PetscScalar a = static_cast<PetscScalar>(a_in);

  const PetscVector<T>* v = dynamic_cast<const PetscVector<T>*>(&v_in);

  genius_assert (v != NULL);
  genius_assert(this->size() == v->size());

  ierr = VecAXPY(_vec, a, v->_vec);

}



template <typename T>
void PetscVector<T>::insert (const std::vector<T>& v,
                             const std::vector<unsigned int>& dof_indices)
{
  genius_assert (v.size() == dof_indices.size());

  for (unsigned int i=0; i<v.size(); i++)
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void PetscVector<T>::insert (const NumericVector<T>& V,
                             const std::vector<unsigned int>& dof_indices)
{
  genius_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::insert (const DenseVector<T>& V,
                             const std::vector<unsigned int>& dof_indices)
{
  genius_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::scale (const T factor_in)
{
  int ierr = 0;
  PetscScalar factor = static_cast<PetscScalar>(factor_in);

  ierr = VecScale(_vec, factor);
}




template <typename T>
T PetscVector<T>::dot (const NumericVector<T>& V) const
{
  // Error flag
  int ierr = 0;

  // Return value
  PetscScalar value=0.;

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T>* v = dynamic_cast<const PetscVector<T>*>(&V);
  genius_assert (v != NULL);

  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecDot(this->_vec, v->_vec, &value);

  return static_cast<T>(value);
}




template <typename T>
NumericVector<T>&
PetscVector<T>::operator = (const T s_in)
{
  genius_assert(this->closed());

  int ierr = 0;
  PetscScalar s = static_cast<PetscScalar>(s_in);

  if (this->size() != 0)
  {
    ierr = VecSet(_vec, s);
  }

  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator = (const NumericVector<T>& v_in)
{
  const PetscVector<T>* v = dynamic_cast<const PetscVector<T>*>(&v_in);

  genius_assert (v != NULL);

  *this = *v;

  return *this;
}



template <typename T>
PetscVector<T>&
PetscVector<T>::operator = (const PetscVector<T>& v)
{
  if (v.initialized())
  {
    this->init (v.size(), v.local_size());
    this->_is_closed      = v._is_closed;
    this->_is_initialized = v._is_initialized;

    if (v.size() != 0)
    {
      int ierr = 0;

      ierr = VecCopy (v._vec, this->_vec);
    }
  }

  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator = (const std::vector<T>& v)
{
  const unsigned int nl   = this->local_size();
  const unsigned int ioff = this->first_local_index();
  int ierr=0;
  PetscScalar* values;

  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (this->size() == v.size())
  {
    ierr = VecGetArray (_vec, &values);

    for (unsigned int i=0; i<nl; i++)
      values[i] =  static_cast<PetscScalar>(v[i+ioff]);

    ierr = VecRestoreArray (_vec, &values);
  }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else
  {
    genius_assert (this->local_size() == v.size());

    ierr = VecGetArray (_vec, &values);

    for (unsigned int i=0; i<nl; i++)
      values[i] = static_cast<PetscScalar>(v[i]);

    ierr = VecRestoreArray (_vec, &values);
  }

  return *this;
}




template <typename T>
void PetscVector<T>::localize (std::vector<T>& v_local) const
{
  // This function must be run on all processors at once
  parallel_only();

  int ierr=0;
  const int n = this->size();
  const int nl = this->local_size();
  PetscScalar *values;

  v_local.clear();
  v_local.resize(n, 0.);

  ierr = VecGetArray (_vec, &values);

  unsigned int ioff = first_local_index();

  for (int i=0; i<nl; i++)
    v_local[i+ioff] = static_cast<T>(values[i]);

  ierr = VecRestoreArray (_vec, &values);

  Parallel::sum(v_local);
}



// Full specialization for Real datatypes
template <>
void PetscVector<Real>::localize_to_one (std::vector<Real>& v_local,
    const unsigned int pid) const
{
  int ierr=0;
  const int n  = size();
  const int nl = local_size();
  PetscScalar *values;


  v_local.resize(n);


  // only one processor
  if (n == nl)
  {
    ierr = VecGetArray (_vec, &values);

    for (int i=0; i<n; i++)
      v_local[i] = static_cast<Real>(values[i]);

    ierr = VecRestoreArray (_vec, &values);
  }

  // otherwise multiple processors
  else
  {
    unsigned int ioff = this->first_local_index();
    std::vector<Real> local_values (n, 0.);

    {
      ierr = VecGetArray (_vec, &values);

      for (int i=0; i<nl; i++)
        local_values[i+ioff] = static_cast<Real>(values[i]);

      ierr = VecRestoreArray (_vec, &values);
    }

    // FIXME: MPI_DOUBLE will fail if REAL=longdouble
    MPI_Reduce (&local_values[0], &v_local[0], n, MPI_DOUBLE, MPI_SUM,
                pid, PETSC_COMM_WORLD);
  }
}




// Full specialization for Complex datatypes
#ifdef USE_COMPLEX_NUMBERS

template <>
void PetscVector<Complex>::localize_to_one (std::vector<Complex>& v_local,
    const unsigned int pid) const
{
  int ierr=0;
  const int n  = size();
  const int nl = local_size();
  PetscScalar *values;


  v_local.resize(n);


  for (int i=0; i<n; i++)
    v_local[i] = 0.;

  // only one processor
  if (n == nl)
  {
    ierr = VecGetArray (_vec, &values);

    for (int i=0; i<n; i++)
      v_local[i] = static_cast<Complex>(values[i]);

    ierr = VecRestoreArray (_vec, &values);
  }

  // otherwise multiple processors
  else
  {
    unsigned int ioff = this->first_local_index();

    /* in here the local values are stored, acting as send buffer for MPI
     * initialize to zero, since we collect using MPI_SUM
     */
    std::vector<Real> real_local_values(n, 0.);
    std::vector<Real> imag_local_values(n, 0.);

    {
      ierr = VecGetArray (_vec, &values);

      // provide my local share to the real and imag buffers
      for (int i=0; i<nl; i++)
      {
        real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
        imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
      }

      ierr = VecRestoreArray (_vec, &values);
    }

    /* have buffers of the real and imaginary part of v_local.
     * Once MPI_Reduce() collected all the real and imaginary
     * parts in these std::vector<double>, the values can be 
     * copied to v_local
     */
    std::vector<Real> real_v_local(n);
    std::vector<Real> imag_v_local(n);

    // collect entries from other proc's in real_v_local, imag_v_local
    MPI_Reduce (&real_local_values[0], &real_v_local[0], n,
                MPI_DOUBLE, MPI_SUM,
                pid, libMesh::COMM_WORLD);

    MPI_Reduce (&imag_local_values[0], &imag_v_local[0], n,
                MPI_DOUBLE, MPI_SUM,
                pid, libMesh::COMM_WORLD);

    // copy real_v_local and imag_v_local to v_local
    for (int i=0; i<n; i++)
      v_local[i] = Complex(real_v_local[i], imag_v_local[i]);
  }
}

#endif



template <typename T>
void PetscVector<T>::print_matlab (const std::string name) const
{
  genius_assert (this->initialized());
  genius_assert (this->closed());

  int ierr=0;
  PetscViewer petsc_viewer;


  ierr = PetscViewerCreate (PETSC_COMM_WORLD,
                            &petsc_viewer);

  /**
   * Create an ASCII file containing the matrix
   * if a filename was provided.  
   */
  if (name != "NULL")
  {
    ierr = PetscViewerASCIIOpen( PETSC_COMM_WORLD,
                                 name.c_str(),
                                 &petsc_viewer);

    ierr = PetscViewerSetFormat (petsc_viewer,
                                 PETSC_VIEWER_ASCII_MATLAB);

    ierr = VecView (_vec, petsc_viewer);
  }

  /**
   * Otherwise the matrix will be dumped to the screen.
   */
  else
  {
    ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
                                 PETSC_VIEWER_ASCII_MATLAB);

    ierr = VecView (_vec, PETSC_VIEWER_STDOUT_WORLD);
  }


  /**
   * Destroy the viewer.
   */
  ierr = PetscViewerDestroy (PetscDestroyObject(petsc_viewer));
}




//------------------------------------------------------------------
// Explicit instantiations
template class PetscVector<Real>;



#endif // #ifdef HAVE_PETSC
