// $Id: petsc_matrix.h 2789 2008-04-13 02:24:40Z roystgnr $

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



#ifndef __petsc_matrix_h__
#define __petsc_matrix_h__

#include "genius_common.h"

#ifdef HAVE_PETSC

// C++ includes
#include <algorithm>


// Local includes
#include "genius_petsc.h"
#include "sparse_matrix.h"
#include "petsc_macro.h"

// Forward Declarations
template <typename T> class DenseMatrix;


/**
 * Petsc include files.
 */
#ifndef USE_COMPLEX_NUMBERS
extern "C"
{
# include <petscmat.h>
}
#else
# include <petscmat.h>
#endif



/**
 * Petsc matrix. Provides a nice interface to the
 * Petsc C-based data structures for parallel,
 * sparse matrices.
 *
 * @author Benjamin S. Kirk, 2002
 */

template <typename T>
class PetscMatrix : public SparseMatrix<T>
{
public:
  /**
   * Constructor; initializes the matrix to
   * be empty, without any structure, i.e.
   * the matrix is not usable at all. This
   * constructor is therefore only useful
   * for matrices which are members of a
   * class. All other matrices should be
   * created at a point in the data flow
   * where all necessary information is
   * available.
   *
   * You have to initialize
   * the matrix before usage with
   * \p init(...).
   */
  PetscMatrix (const unsigned int m,   const unsigned int n,
               const unsigned int m_l, const unsigned int n_l,
               const std::vector<int> &n_nz=std::vector<int>(),
               const std::vector<int> &n_oz=std::vector<int>());


  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  ~PetscMatrix ();

  /**
   * set the matrix sparsity pattern. When your \p SparseMatrix<T>
   * implementation does not need this data simply do
   * not overload this method.
   */
  void set_sparsity_pattern (const std::vector<int> &n_nz, const std::vector<int> &n_oz);

  /**
   * Initialize a Petsc matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.
   */
  void init ();


  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. 
   */
  void clear ();

  /**
   * Set all entries to 0. This method retains 
   * sparsity structure.
   */
  void zero ();

  /**
   * Call the Petsc assemble routines.
   * sends necessary messages to other
   * processors
   */
  void close (bool final) ;


  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  void set (const unsigned int i,
            const unsigned int j,
            const T value);

  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  void add (const unsigned int i,
              const unsigned int j,
              const T value);

  /**
   * Add the full matrix to the
   * Petsc matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */

  void add_matrix (const DenseMatrix<T> &dm,
                   const std::vector<unsigned int> &rows,
                   const std::vector<unsigned int> &cols);

  /**
   * Add the full matrix to the
   * Sparse matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
  void add_matrix (const T* dm,
                   unsigned int m, unsigned int * rows,
                   unsigned int n, unsigned int * cols);

  /**
   * Add the source row entries to the destination row
   */
  void add_row_to_row(const std::vector<unsigned int> &src_rows,
                      const std::vector<unsigned int> &dst_rows);

  /**
   * clear the given row and fill diag with given value
   */
  void clear_row(const std::vector<unsigned int> &rows, const T diag=T(0.0) );


  /**
   * Return the value of the entry
   * \p (i,j).  This may be an
   * expensive operation and you
   * should always take care where
   * to call this function.  In
   * order to avoid abuse, this
   * function throws an exception
   * if the required element does
   * not exist in the matrix.
   *
   * In case you want a function
   * that returns zero instead (for
   * entries that are not in the
   * sparsity pattern of the
   * matrix), use the \p el
   * function.
   */
  T operator () (const unsigned int i,
                 const unsigned int j) const;

  /**
   * Return the l1-norm of the matrix, that is
   * \f$|M|_1=max_{all columns j}\sum_{all 
   * rows i} |M_ij|\f$,
   * (max. sum of columns).
   * This is the
   * natural matrix norm that is compatible
   * to the l1-norm for vectors, i.e.
   * \f$|Mv|_1\leq |M|_1 |v|_1\f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  Real l1_norm () const;

  /**
   * Return the linfty-norm of the
   * matrix, that is
   * \f$|M|_infty=max_{all rows i}\sum_{all 
   * columns j} |M_ij|\f$,
   * (max. sum of rows).
   * This is the
   * natural matrix norm that is compatible
   * to the linfty-norm of vectors, i.e.
   * \f$|Mv|_infty \leq |M|_infty |v|_infty\f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  Real linfty_norm () const;

  /**
   * see if Petsc matrix has been closed
   * and fully assembled yet
   */
  bool closed() const;

  /**
   * Print the contents of the matrix to the screen
   * with the PETSc viewer.  This function only allows
   * printing to standard out, this is because we have
   * limited ourselves to one PETSc implementation for
   * writing.
   */
  void print_personal(std::ostream& os=std::cout) const;

  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  void print_matlab(const std::string name="NULL") const;


  /**
   * Returns the raw PETSc matrix context pointer.  Note this is generally
   * not required in user-level code. Just don't do anything crazy like
   * calling MatDestroy()!
   */
  Mat mat () { return _mat; }


private:


  /**
   * Petsc matrix datatype to store values
   */
  Mat _mat;


  /**
   * Petsc add or insert flag
   */
  InsertMode _add_value_flag;

  /**
   * This boolean value should only be set to false
   * for the constructor which takes a PETSc Mat object. 
   */
  bool _destroy_mat_on_exit;
};






#endif // #ifdef HAVE_PETSC
#endif // #ifdef __petsc_matrix_h__
