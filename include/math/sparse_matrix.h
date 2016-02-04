// $Id: sparse_matrix.h 2789 2008-04-13 02:24:40Z roystgnr $

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



#ifndef __sparse_matrix_h__
#define __sparse_matrix_h__


// C++ includes
#include <iomanip>
#include <vector>

// Local includes
#include "genius_common.h"
#include "auto_ptr.h"

// forward declarations
template <typename T> class SparseMatrix;
template <typename T> inline std::ostream& operator << (std::ostream& os, const SparseMatrix<T>& m);



/**
 * Generic sparse matrix. This class contains
 * pure virtual members that must be overloaded
 * in derived classes.  Using a derived class
 * allows for uniform access to sparse matrices
 * from various different solver packages in
 * different formats.
 *
 * @author Benjamin S. Kirk, 2003
 */

template <typename T>
class SparseMatrix
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
  SparseMatrix (const unsigned int m,   const unsigned int n,
                const unsigned int m_l, const unsigned int n_l);

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  virtual ~SparseMatrix ();

  /**
   * Builds a \p SparseMatrix<T> using the linear solver package specified by
   * \p solver_package
   */
  static AutoPtr<SparseMatrix<T> >  build(const std::string & solver_package,
                                          const unsigned int m,   const unsigned int n,
                                          const unsigned int m_l, const unsigned int n_l);

  /**
   * @returns true if the matrix has been initialized,
   * false otherwise.
   */
  virtual bool initialized() const { return _is_initialized; }


  /**
   * get the matrix sparsity pattern. When your \p SparseMatrix<T>
   * implementation does not need this data simply do
   * not overload this method.
   */
  virtual void get_sparsity_pattern (std::vector<int> &n_nz, std::vector<int> &n_oz) {}
  
  
  /**
   * get the matrix nonzero pattern. When your \p SparseMatrix<T>
   * implementation does not need this data simply do
   * not overload this method.
   */
  virtual void get_nonzero_pattern (std::vector< std::vector<int> > & nz) {}

  
  /**
   * Initialize a Sparse matrix
   */
  virtual void init() = 0;

  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. 
   */
  virtual void clear () = 0;

  /**
   * Set all entries to 0.
   */
  virtual void zero () = 0;

  /**
   * Call the Sparse assemble routines.
   * sends necessary messages to other
   * processors
   */
  virtual void close (bool final) = 0;

  /**
   * @returns \p m, the row-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */
  unsigned int m () const { return _m_global; }

  /**
   * @returns \p n, the column-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */
  unsigned int n () const { return _n_global; }

  /**
   * return row_start, the index of the first
   * matrix row stored on this processor
   */
  unsigned int row_start () const { return _global_offset; }

  /**
   * @return row_stop, the index of the last
   * matrix row (+1) stored on this processor
   */
  unsigned int row_stop () const { return _global_offset + _m_local; }

  /**
   * @return true iff the given row in global index on this processor
   */
  bool row_on_processor(unsigned int row) const
  { return row >= _global_offset && row < _global_offset + _m_local; }


  /**
   * @return true iff the given col in global index on this processor
   */
  bool col_on_processor(unsigned int col) const
  { return col >= _global_offset && col < _global_offset + _n_local; }

  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  virtual void set (  const unsigned int i,
                      const unsigned int j,
                      const T value) = 0;

  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  virtual void add (const unsigned int i,
                    const unsigned int j,
                    const T value) = 0;

                    
  /**
   * Add a row to the Sparse matrix.  
   */
  virtual void add_row (unsigned int row, 
                        const std::vector<unsigned int> &cols, 
                        const T* dm) = 0;

  /**
   * Add a row to the Sparse matrix.  
   */
  virtual void add_row (unsigned int row, 
                        unsigned int n, const unsigned int * cols, 
                        const T* dm) = 0;

  /**
   * Add a row to the Sparse matrix.  
   */
  virtual void add_row (unsigned int row, 
                        int n, const int * cols, 
                        const T* dm) = 0;

  /**
   * Add the full matrix to the
   * Sparse matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
  virtual void add_matrix (const std::vector<unsigned int> &rows,
                           const std::vector<unsigned int> &cols,
                           const T* dm) = 0;

  /**
   * Add the full matrix to the
   * Sparse matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
  virtual void add_matrix (unsigned int m, unsigned int * rows,
                           unsigned int n, unsigned int * cols,
                           const T* dm) = 0;


  /**
   * Add the source row entries to the destination row
   */
  virtual void add_row_to_row(const std::vector<int> &src_rows,
                              const std::vector<int> &dst_rows) = 0;

                              
  /**
   * clear the given row and fill diag with given value
   */
  virtual void clear_row(int row, const T diag=T(0.0) ) = 0;
  
  /**
   * clear the given row and fill diag with given value
   */
  virtual void clear_row(const std::vector<int> &rows, const T diag=T(0.0) ) = 0;

  /**
   * get a local row from Sparse matrix.  
   */                      
  virtual void get_row (unsigned int row, int n, const int * cols, T* dm) = 0;

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
  virtual T operator () (const unsigned int i,
                         const unsigned int j) const = 0;


  /**
   * see if Sparse matrix has been closed
   * and fully assembled yet
   */
  virtual bool closed() const = 0;

  /**
   * Print the contents of the matrix to the screen
   * in a uniform style, regardless of matrix/solver
   * package being used.
   */
  void print(std::ostream& os=std::cout) const;

  /**
   * Same as the print method above, but allows you
   * to print to a stream in the standard syntax.
   */
  template <typename U>
  friend std::ostream& operator << (std::ostream& os, const SparseMatrix<U>& m);

  /**
   * Print the contents of the matrix to the screen
   * in a package-personalized style, if available.
   */
  virtual void print_personal(std::ostream& os=std::cout) const = 0;

  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  virtual void print_matlab(const std::string name="NULL") const
  {
    std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
    std::cerr << "ERROR writing MATLAB file " << name << std::endl;
    genius_error();
  }


protected:

  /**
   * global row size
   */
  unsigned int _m_global;

  /**
   * global column size
   */
  unsigned int _n_global;

  /**
   * row size of local block
   */
  unsigned int _m_local;

  /**
   * column size of local block
   */
  unsigned int _n_local;

  /**
   * the offset of local block of dofs at global dofs
   */
  unsigned int _global_offset;

  /**
   * Flag indicating whether or not the matrix
   * has been initialized.
   */
  bool _is_initialized;
};



//-----------------------------------------------------------------------
// SparseMatrix inline members



// For SGI MIPSpro this implementation must occur after
// the partial specialization of the print() member.
template <typename T>
inline
std::ostream& operator << (std::ostream& os, const SparseMatrix<T>& m)
{
  m.print(os);
  return os;
}



#endif // #ifndef __sparse_matrix_h__
