// $Id: petsc_matrix.C 2789 2008-04-13 02:24:40Z roystgnr $

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
#include "config.h"

#include <map>

#ifdef HAVE_PETSC

// Local includes
#include "petsc_matrix.h"
#include "dense_matrix.h"



//-----------------------------------------------------------------------
// PetscMatrix members

//-----------------------------------------------------------------------
// PetscMatrix inline members
template <typename T>
PetscMatrix<T>::PetscMatrix(const unsigned int m,   const unsigned int n,
                            const unsigned int m_l, const unsigned int n_l,
                            const std::vector<int> &n_nz, const std::vector<int> &n_oz)
  : SparseMatrix<T>(m,n,m_l,n_l), _add_value_flag(NOT_SET_VALUES), _destroy_mat_on_exit(true)
{
  if ((m==0) || (n==0))
    return;

  int ierr     = 0;
  int m_global = static_cast<int>(m);
  int n_global = static_cast<int>(n);
  int m_local  = static_cast<int>(m_l);
  int n_local  = static_cast<int>(n_l);

  ierr = MatCreate(PETSC_COMM_WORLD,&_mat); genius_assert(!ierr);
  ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global); genius_assert(!ierr);

  // create a sequential matrix on one processor
  if ((m_l == m) && (n_l == n))
  {
    ierr = MatSetType(_mat, MATSEQAIJ); genius_assert(!ierr);
    // alloc memory for sequence matrix here
    ierr = MatSeqAIJSetPreallocation(_mat, 0, &n_nz[0]); genius_assert(!ierr);
  }
  else
  {
    ierr = MatSetType(_mat, MATMPIAIJ); genius_assert(!ierr);
    // alloc memory for parallel matrix here
    ierr = MatMPIAIJSetPreallocation(_mat, 0, &n_nz[0], 0, &n_oz[0]); genius_assert(!ierr);
  }

  // indicates when PetscUtils::MatZeroRows() is called the zeroed entries are kept in the nonzero structure
#if PETSC_VERSION_GE(3,1,0)
  ierr = MatSetOption(_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); genius_assert(!ierr);
#endif

#if PETSC_VERSION_EQ(3,0,0)
  ierr = MatSetOption(_mat, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE); genius_assert(!ierr);
#endif

  ierr = MatSetFromOptions(_mat); genius_assert(!ierr);

  SparseMatrix<T>::_is_initialized = true;
}


template <typename T>
PetscMatrix<T>::~PetscMatrix()
{
  this->clear();
}


template <typename T>
void PetscMatrix<T>::set (const unsigned int i, const unsigned int j, const T value)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==INSERT_VALUES || _add_value_flag==NOT_SET_VALUES);

  int ierr=0, i_val=i, j_val=j;
  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val, &petsc_value, INSERT_VALUES); genius_assert(!ierr);
}



template <typename T>
void PetscMatrix<T>::add (const unsigned int i, const unsigned int j, const T value)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==ADD_VALUES || _add_value_flag==NOT_SET_VALUES);

  int ierr=0, i_val=i, j_val=j;
  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val, &petsc_value, ADD_VALUES); genius_assert(!ierr);
}


template <typename T>
void PetscMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
                                const std::vector<unsigned int>& rows,
                                const std::vector<unsigned int>& cols)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==ADD_VALUES || _add_value_flag==NOT_SET_VALUES);

  const unsigned int m = dm.m();
  const unsigned int n = dm.n();

  genius_assert (rows.size() == m);
  genius_assert (cols.size() == n);

  int ierr=0;

  ierr = MatSetValues(_mat,
                      m, (int*) &rows[0],
                      n, (int*) &cols[0],
                      (PetscScalar*) &dm.get_values()[0],
                      ADD_VALUES);
  genius_assert(!ierr);
}



template <typename T>
void PetscMatrix<T>::add_matrix (const T* dm,
                   unsigned int m, unsigned int * rows,
                   unsigned int n, unsigned int * cols)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==ADD_VALUES || _add_value_flag==NOT_SET_VALUES);

  int ierr=0;

  ierr = MatSetValues(_mat,
                      m, (int*) rows,
                      n, (int*) cols,
                      (PetscScalar*) dm,
                      ADD_VALUES);
  genius_assert(!ierr);
}




template <typename T>
bool PetscMatrix<T>::closed() const
{
  genius_assert (this->initialized());

  int ierr=0;
  PetscBool assembled;

  ierr = MatAssembled(_mat, &assembled);

  return (assembled == PETSC_TRUE);
}




template <typename T>
void PetscMatrix<T>::init ()
{

}


template <typename T>
void PetscMatrix<T>::close (bool )
{
  // BSK - 1/19/2004
  // strictly this check should be OK, but it seems to
  // fail on matrix-free matrices.  Do they falsely
  // state they are assembled?  Check with the developers...
  //   if (this->closed())
  //     return;

  int ierr=0;

  ierr = MatAssemblyBegin (_mat, MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd   (_mat, MAT_FINAL_ASSEMBLY);
}



template <typename T>
void PetscMatrix<T>::zero ()
{
  genius_assert (this->initialized());

  int ierr=0;

  ierr = MatZeroEntries(_mat);
}



template <typename T>
void PetscMatrix<T>::clear ()
{
  int ierr=0;

  if ((this->initialized()) && (this->_destroy_mat_on_exit))
  {
    ierr = MatDestroy (PetscDestroyObject(_mat));

    this->_is_initialized = false;
  }
}



template <typename T>
Real PetscMatrix<T>::l1_norm () const
{
  genius_assert (this->initialized());

  int ierr=0;
  PetscReal petsc_value;
  Real value;

  genius_assert (this->closed());

  ierr = MatNorm(_mat, NORM_1, &petsc_value);

  value = static_cast<Real>(petsc_value);

  return value;
}



template <typename T>
Real PetscMatrix<T>::linfty_norm () const
{
  genius_assert (this->initialized());

  int ierr=0;
  PetscReal petsc_value;
  Real value;

  genius_assert (this->closed());

  ierr = MatNorm(_mat, NORM_INFINITY, &petsc_value);

  value = static_cast<Real>(petsc_value);

  return value;
}


template <typename T>
T PetscMatrix<T>::operator () (const unsigned int i, const unsigned int j) const
{
  genius_assert (this->initialized());

  const PetscScalar *petsc_row;
  const PetscInt    *petsc_cols;


  T value=0.;

  int ierr=0;
  int ncols=0;
  int i_val=static_cast<int>(i);
  int j_val=static_cast<int>(j);


  ierr = MatGetRow(_mat, i_val, &ncols, &petsc_cols, &petsc_row);

  // Perform a binary search to find the contiguous index in
  // petsc_cols (resp. petsc_row) corresponding to global index j_val
  std::pair<const int*, const int*> p =
    std::equal_range (&petsc_cols[0], &petsc_cols[0] + ncols, j_val);

  // Found an entry for j_val
  if (p.first != p.second)
  {
    // The entry in the contiguous row corresponding
    // to the j_val column of interest
    const int j = std::distance (const_cast<int*>(&petsc_cols[0]),
                                 const_cast<int*>(p.first));

    genius_assert (j < ncols);
    genius_assert (petsc_cols[j] == j_val);

    value = static_cast<T> (petsc_row[j]);

    ierr  = MatRestoreRow(_mat, i_val, &ncols, &petsc_cols, &petsc_row);

    return value;
  }

  // Otherwise the entry is not in the sparse matrix,
  // i.e. it is 0.
  return 0.;
}







template <typename T>
void PetscMatrix<T>::add_row_to_row(const std::vector<unsigned int> &src_rows,
                                    const std::vector<unsigned int> &dst_rows)
{

  // test if the matrix is assembled
  // note: the test is not work properly! if it is a bug...

  //PetscBool assembled;
  //MatAssembled(mat, &assembled);

  //if( !assembled )
  {
    MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
  }

  if( src_rows.size() )
  {
    genius_assert(src_rows.size() == dst_rows.size());

    std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > > row_add_to_buffer;

    // read the row
    for(unsigned int nrow=0; nrow<src_rows.size(); nrow++)
    {
      PetscInt ncols;
      const PetscInt * row_cols_pointer;
      const PetscScalar * row_vals_pointer;

      MatGetRow(_mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);

      // save the values
      std::vector<PetscInt>    row_cols;
      std::vector<PetscScalar> row_vals;
      for(PetscInt i=0;i <ncols; i++)
      {
        row_cols.push_back(row_cols_pointer[i]);
        row_vals.push_back(row_vals_pointer[i]);
      }

      row_add_to_buffer.insert(std::pair< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >
                               (dst_rows[nrow], std::pair<std::vector<PetscInt>, std::vector<PetscScalar> >(row_cols, row_vals)));

      // restore pointers
      MatRestoreRow(_mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);
    }

    //ok, we add rows to destination
    std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >::iterator row_add_it =  row_add_to_buffer.begin();
    for(; row_add_it !=  row_add_to_buffer.end(); ++row_add_it)
    {
      MatSetValues(_mat, 1, &((*row_add_it).first) ,
                   (*row_add_it).second.first.size(), &((*row_add_it).second.first[0]) ,
                   &((*row_add_it).second.second[0]), ADD_VALUES);

    }

  }

  MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
}


template <typename T>
void PetscMatrix<T>::clear_row(const std::vector<unsigned int> &rows, const T diag)
{
#if PETSC_VERSION_GE(3,2,0)
  MatZeroRows(_mat, rows.size(), (int*)&rows[0], diag, PETSC_NULL, PETSC_NULL);
#else
  MatZeroRows(_mat, rows.size(), (int*)&rows[0], diag);
#endif
}



template <typename T>
void PetscMatrix<T>::print_personal(std::ostream& os) const
{
  genius_assert (this->initialized());

  int ierr=0;

  ierr = MatView(_mat, PETSC_VIEWER_STDOUT_SELF);
}


template <typename T>
void PetscMatrix<T>::print_matlab (const std::string name) const
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

    ierr = MatView (_mat, petsc_viewer);
  }

  /**
   * Otherwise the matrix will be dumped to the screen.
   */
  else
  {
    ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
                                 PETSC_VIEWER_ASCII_MATLAB);

    ierr = MatView (_mat, PETSC_VIEWER_STDOUT_WORLD);
  }


  /**
   * Destroy the viewer.
   */
  ierr = PetscViewerDestroy (PetscDestroyObject(petsc_viewer));
}




//------------------------------------------------------------------
// Explicit instantiations
template class PetscMatrix<Real>;


#endif // #ifdef HAVE_PETSC
