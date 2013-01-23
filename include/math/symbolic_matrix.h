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



#ifndef __symbolic_matrix_h__
#define __symbolic_matrix_h__

#include "genius_common.h"
#include "sparse_matrix.h"

// C++ includes
#include <algorithm>

#if defined(HAVE_TR1_UNORDERED_SET)
#include <tr1/unordered_set>
#elif defined(HAVE_TR1_UNORDERED_SET_WITH_STD_HEADER) || defined(HAVE_UNORDERED_SET)
#include <unordered_set>
#else
#include <set>
#endif
#include <functional>

// Forward Declarations
template <typename T> class DenseMatrix;



/**
 * Symbolic matrix. Create fill-in pattern of Sparse Matrix
 */
template <typename T>
class SymbolicMatrix : public SparseMatrix<T>
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
  SymbolicMatrix (const unsigned int m,   const unsigned int n,
                  const unsigned int m_l, const unsigned int n_l);


  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  ~SymbolicMatrix ();

  /**
   * Initialize a Petsc matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.
   */
  void init ();


  /**
   * get the matrix sparsity pattern. When your \p SparseMatrix<T>
   * implementation does not need this data simply do
   * not overload this method.
   */
  void get_sparsity_pattern (std::vector<int> &n_nz, std::vector<int> &n_oz);


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
  void zero () {}

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
   * matrix.  This is useful
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
  void clear_row(const std::vector<unsigned int> &, const T =T(0.0) ) {}

  /**
   * Return the value of the entry \p (i,j).
   */
  T operator () (const unsigned int i, const unsigned int j) const { return 0.0; }

  /**
   * Return the l1-norm of the matrix
   */
  Real l1_norm () const { return 0.0; }

  /**
   * Return the linfty-norm of the matrix
   */
  Real linfty_norm () const { return 0.0; }

  /**
   * see if  matrix has been closed
   * and fully assembled yet
   */
  bool closed() const { return _closed; }

  /**
   * Print the contents of the matrix to the screen
   */
  void print_personal(std::ostream& os=std::cout) const {}

private:

  struct aij_type
  {
    /**
     * the index of column entries of each local row
     */
    std::vector<unsigned int> i;

    /**
     * the max number of column entries of each local row
     */
    std::vector<unsigned int> imax;

    /**
     * the occupied number of column entries of each local row
     */
    std::vector<unsigned int> iocp;

    /**
     * index of all the local column entries
     */
    std::vector<unsigned int> j;
  };



  typedef std::pair<unsigned int, unsigned int> key_type;

  struct key_hash
  {
    size_t operator()(const key_type & key) const
    {
      return ((key.first << 19) | (key.second << 7));
    }
  };

#if defined(HAVE_UNORDERED_SET)
  typedef std::unordered_set<key_type, key_hash> set_type;
#elif defined(HAVE_TR1_UNORDERED_SET) || defined(HAVE_TR1_UNORDERED_SET_WITH_STD_HEADER)
  typedef std::tr1::unordered_set<key_type, key_hash> set_type;
#else
  typedef std::set<key_type>  set_type;
#endif

private:

  /**
   * approx matrix band width, default is 30
   */
  unsigned int _nz;


  /**
   * local entries
   */
  aij_type _local_entries; 


  /**
   * extra local entries
   */
  set_type _extra_local_entries;


  /**
   * off processor entries
   */
  set_type _off_processor_entries;


  /**
   * nonlocal entries
   */
  set_type _nonlocal_entries;



  /**
   * set a entry of sparse matrix
   */
  void _set_entry(unsigned int i, unsigned int j);


  /**
   * set entry of aij struct
   */
  void _set_aij_entry(aij_type & aij_entries, unsigned int i, unsigned int j);


  /**
   * hash value of matrix pattern
   */
  size_t _hash_entry;

  /**
   * build the hash value of all the entries 
   */
  size_t _build_hash_entry() const;

  /**
   * flag for same matrix pattern
   */
  bool _same_pattern;

  /**
   * matrix is ready
   */
  bool _closed;

};




//-----------------------------------------------------------------------
// PetscMatrix inline members
template <typename T>
inline SymbolicMatrix<T>::SymbolicMatrix(const unsigned int m,   const unsigned int n,
                                         const unsigned int m_l, const unsigned int n_l)
  : SparseMatrix<T>(m,n,m_l,n_l), _nz(30), _hash_entry(0), _same_pattern(false), _closed(false)
{}



template <typename T>
inline SymbolicMatrix<T>::~SymbolicMatrix()
{
  this->clear();
}





#endif // #ifdef __symbolic_matrix_h__
