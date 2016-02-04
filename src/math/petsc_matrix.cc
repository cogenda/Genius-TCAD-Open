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


#ifdef HAVE_PETSC

// Local includes
#include "petsc_matrix.h"
#include "parallel.h"



//-----------------------------------------------------------------------
// PetscMatrix members

//-----------------------------------------------------------------------
// PetscMatrix inline members
template <typename T>
PetscMatrix<T>::PetscMatrix(const unsigned int m,   const unsigned int n,
                            const unsigned int m_l, const unsigned int n_l)
  : SparseMatrix<T>(m,n,m_l,n_l), 
    _mat_buf_mode(true), 
    _add_value_flag(NOT_SET_VALUES), 
    _closed(false), 
    _destroy_mat_on_exit(false)
{
  if ((m==0) || (n==0))
    return;
  
   _mat_local.resize(m_l);
  
  int ierr = 0;
  ierr = MatCreate(PETSC_COMM_WORLD,&_mat); genius_assert(!ierr);
  ierr = MatSetSizes(_mat, m_l, n_l, m, n); genius_assert(!ierr);

  // create a sequential matrix on one processor
  if (Genius::n_processors()==1)
  {
    ierr = MatSetType(_mat, MATSEQAIJ); genius_assert(!ierr);
  }
  else
  {
    ierr = MatSetType(_mat, MATMPIAIJ); genius_assert(!ierr);
  }
  
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
  
  if(_mat_buf_mode)
  {
    if( SparseMatrix<T>::row_on_processor(i) )
      _mat_local[i-SparseMatrix<T>::_global_offset][j] = value;
    else 
      _mat_nonlocal[std::make_pair(i,j)] = value;
  }
  else
  {
    int ierr=0, i_val=i, j_val=j;
    PetscScalar petsc_value = static_cast<PetscScalar>(value);
    ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val, &petsc_value, INSERT_VALUES); genius_assert(!ierr);
  }
  
  _closed = false;
  _add_value_flag=INSERT_VALUES;
}



template <typename T>
void PetscMatrix<T>::add (const unsigned int i, const unsigned int j, const T value)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==ADD_VALUES || _add_value_flag==NOT_SET_VALUES);

  if(_mat_buf_mode)
  {
    if( SparseMatrix<T>::row_on_processor(i) )
      _mat_local[i-SparseMatrix<T>::_global_offset][j] += value;
    else 
      _mat_nonlocal[std::make_pair(i,j)] += value;
  }
  else
  {
    int ierr=0, i_val=i, j_val=j;
    PetscScalar petsc_value = static_cast<PetscScalar>(value);
    ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val, &petsc_value, ADD_VALUES); genius_assert(!ierr);
  }
  _closed = false;
  _add_value_flag=ADD_VALUES;
}



template <typename T>
void PetscMatrix<T>::add_row (unsigned int row, const std::vector<unsigned int> &cols, const T* dm)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==ADD_VALUES || _add_value_flag==NOT_SET_VALUES);

  if(_mat_buf_mode)
  {
    if( SparseMatrix<T>::row_on_processor(row) )
    {
      for(unsigned int j=0; j<cols.size(); j++)
         _mat_local[row-SparseMatrix<T>::_global_offset][cols[j]] += dm[j];
    }
    else
    {
      for(unsigned int j=0; j<cols.size(); j++)
        _mat_nonlocal[std::make_pair(row,cols[j])] += dm[j];
    }
  }
  else
  {
    int ierr=0;

    ierr = MatSetValues(_mat, 1, (int*) &row, cols.size(), (int*) &cols[0], (PetscScalar*) dm, ADD_VALUES);
    genius_assert(!ierr);
  }
  
  _closed = false;
  _add_value_flag=ADD_VALUES;
}

template <typename T>
void PetscMatrix<T>::add_row (unsigned int row, unsigned int n, const unsigned int * cols, const T* dm)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==ADD_VALUES || _add_value_flag==NOT_SET_VALUES);

  if(_mat_buf_mode)
  {
    if( SparseMatrix<T>::row_on_processor(row) )
    {
      for(unsigned int j=0; j<n; j++)
         _mat_local[row-SparseMatrix<T>::_global_offset][cols[j]] += dm[j];
    }
    else
    {
      for(unsigned int j=0; j<n; j++)
        _mat_nonlocal[std::make_pair(row,cols[j])] += dm[j];
    }
  }
  else
  {
    int ierr=0;

    ierr = MatSetValues(_mat, 1, (int*) &row, n, (int*) &cols[0], (PetscScalar*) dm, ADD_VALUES);
    genius_assert(!ierr);
  }
  
  _closed = false;
  _add_value_flag=ADD_VALUES;
}


template <typename T>
void PetscMatrix<T>::add_row (unsigned int row, int n, const int * cols, const T* dm)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==ADD_VALUES || _add_value_flag==NOT_SET_VALUES);

  if(_mat_buf_mode)
  {
    if( SparseMatrix<T>::row_on_processor(row) )
    {
      for(int j=0; j<n; j++)
         _mat_local[row-SparseMatrix<T>::_global_offset][cols[j]] += dm[j];
    }
    else
    {
      for(int j=0; j<n; j++)
        _mat_nonlocal[std::make_pair(row,cols[j])] += dm[j];
    }
  }
  else
  {
    int ierr=0;

    ierr = MatSetValues(_mat, 1, (int*) &row, n, (int*) &cols[0], (PetscScalar*) dm, ADD_VALUES);
    genius_assert(!ierr);
  }
  
  _closed = false;
  _add_value_flag=ADD_VALUES;
}


template <typename T>
void PetscMatrix<T>::add_matrix(const std::vector<unsigned int>& rows,
                                const std::vector<unsigned int>& cols,
                                const T* dm)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==ADD_VALUES || _add_value_flag==NOT_SET_VALUES);

  const unsigned int m = rows.size();
  const unsigned int n = cols.size();

  if(_mat_buf_mode)
  {
    for(unsigned int i=0; i<m; i++)
    {
      if( SparseMatrix<T>::row_on_processor(rows[i]) )
      {
        for(unsigned int j=0; j<n; j++)
          _mat_local[rows[i]-SparseMatrix<T>::_global_offset][cols[j]] += dm[i*n+j];
      }
      else
        for(unsigned int j=0; j<n; j++)
          _mat_nonlocal[std::make_pair(rows[i],cols[j])] += dm[i*n+j];
    }
  }
  else
  {
    int ierr=0;

    ierr = MatSetValues(_mat, m, (int*) &rows[0], n, (int*) &cols[0], (PetscScalar*) dm, ADD_VALUES);
    genius_assert(!ierr);
  }
  
  _closed = false;
  _add_value_flag=ADD_VALUES;
}




template <typename T>
void PetscMatrix<T>::add_matrix (unsigned int m, unsigned int * rows,
                                 unsigned int n, unsigned int * cols,
                                 const T* dm)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==ADD_VALUES || _add_value_flag==NOT_SET_VALUES);

  if(_mat_buf_mode)
  {
    for(unsigned int i=0; i<m; i++)
    {
      if( SparseMatrix<T>::row_on_processor(rows[i]) )
      {
        for(unsigned int j=0; j<n; j++)
          _mat_local[rows[i]-SparseMatrix<T>::_global_offset][cols[j]] += dm[i*n+j];
      }
      else
        for(unsigned int j=0; j<n; j++)
          _mat_nonlocal[std::make_pair(rows[i],cols[j])] += dm[i*n+j];
    }
  }
  else
  {
    int ierr=0;

    ierr = MatSetValues(_mat, m, (int*) &rows[0], n, (int*) &cols[0], (PetscScalar*) dm, ADD_VALUES);
    genius_assert(!ierr);
  }
  
  _closed = false;
  _add_value_flag=ADD_VALUES;
}




template <typename T>
bool PetscMatrix<T>::closed() const
{
  genius_assert (this->initialized());

  if(_mat_buf_mode)
  {
    return _closed;
  }
  else 
  {
    int ierr=0;
    PetscBool assembled;

    ierr = MatAssembled(_mat, &assembled);

    return (assembled == PETSC_TRUE);
  }
  
  return _closed;
}




template <typename T>
void PetscMatrix<T>::init ()
{

}


template <typename T>
void PetscMatrix<T>::close (bool final)
{
  if(_mat_buf_mode)
  {
    unsigned int nonlocal_entries = _mat_nonlocal.size();
    Parallel::sum(nonlocal_entries);
    if(nonlocal_entries)
    {
      std::vector<unsigned int> rows;
      std::vector<unsigned int> cols;
      std::vector<T> values;
      for(typename std::map< std::pair<unsigned int, unsigned int>, T >::const_iterator it =_mat_nonlocal.begin(); 
        it !=_mat_nonlocal.end(); ++it)
      {
        rows.push_back(it->first.first);
        cols.push_back(it->first.second);
        values.push_back(it->second);
      }
    
      _mat_nonlocal.clear();
    
      Parallel::allgather(rows);
      Parallel::allgather(cols);
      Parallel::allgather(values);
    
      for(unsigned int n=0; n<rows.size(); ++n)
      {
        unsigned int row = rows[n];
        if( !SparseMatrix<T>::row_on_processor(row) ) continue;

        unsigned int col = cols[n];
        T value = values[n];
        _mat_local[row-SparseMatrix<T>::_global_offset][col] += value;
      }
    }
    
    _closed = true;
    
    if(final) flush_buf();
  }
  else
  {
    int ierr=0;
    ierr = MatAssemblyBegin (_mat, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd   (_mat, MAT_FINAL_ASSEMBLY);
  }
  
  
}



template <typename T>
void PetscMatrix<T>::zero ()
{
  if(_mat_buf_mode)
  {
    for(size_t n=0; n<_mat_local.size(); ++n)
    {
      std::map<unsigned int, T> & cols = _mat_local[n];
      for(typename std::map<unsigned int, T>::iterator it=cols.begin(); it!=cols.end(); it++)
        it->second = 0.0;
    }
    
    for(typename std::map< std::pair<unsigned int, unsigned int>, T >::iterator it= _mat_nonlocal.begin();it!=_mat_nonlocal.end(); it++)
      it->second = 0.0;  
  }
  else
  {
    genius_assert (this->initialized());

    int ierr=0;

    ierr = MatZeroEntries(_mat);
  }
}



template <typename T>
void PetscMatrix<T>::clear ()
{
  _mat_local.clear();
  _mat_nonlocal.clear();
  
  int ierr=0;

  if ((this->initialized()) && (this->_destroy_mat_on_exit))
  {
    ierr = MatDestroy (PetscDestroyObject(_mat));

    this->_is_initialized = false;
  }

}


template <typename T>
void PetscMatrix<T>::get_row (unsigned int row, int n, const int * cols, T* dm)
{
  if(_mat_buf_mode)
  {
    const std::map<unsigned int, T> & buf = _mat_local[row-SparseMatrix<T>::_global_offset];
    for(int i=0; i<n; i++)
    {
      dm[i] = buf.find(cols[i])->second;
    }
  }
  else
  {
    MatGetValues(_mat, 1, (int*)&row, n, (int*)cols, (PetscScalar*)dm);
  }
}


template <typename T>
T PetscMatrix<T>::operator () (const unsigned int i, const unsigned int j) const
{
  genius_assert (this->initialized());
  
  if(_mat_buf_mode)
  {
    const std::map<unsigned int, T> & buf = _mat_local[i-SparseMatrix<T>::_global_offset];
  
    typename std::map<unsigned int, T>::const_iterator ent =  buf.find(j);
  
    if(ent != buf.end() ) return ent->second;

    // Otherwise the entry is not in the sparse matrix,
    // i.e. it is 0.
    return 0.0;
  }
  
  // else 

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
  std::pair<const int*, const int*> p = std::equal_range (&petsc_cols[0], &petsc_cols[0] + ncols, j_val);

  // Found an entry for j_val
  if (p.first != p.second)
  {  
    // The entry in the contiguous row corresponding
    // to the j_val column of interest
    const int j = std::distance (const_cast<int*>(&petsc_cols[0]), const_cast<int*>(p.first));

    genius_assert (j < ncols);
    genius_assert (petsc_cols[j] == j_val);

    value = static_cast<T> (petsc_row[j]);

    ierr  = MatRestoreRow(_mat, i_val, &ncols, &petsc_cols, &petsc_row);

    return value;
  }

  // Otherwise the entry is not in the sparse matrix,
  // i.e. it is 0.
  return 0.0;
  
}







template <typename T>
void PetscMatrix<T>::add_row_to_row(const std::vector<int> &src_rows,
                                    const std::vector<int> &dst_rows)
{
  if(_mat_buf_mode)
  {
    genius_assert(_closed);

    for(unsigned int n=0; n<src_rows.size(); n++)
    {
      unsigned int src_row = static_cast<unsigned int>(src_rows[n]);
      unsigned int dst_row = static_cast<unsigned int>(dst_rows[n]);

      genius_assert(SparseMatrix<T>::row_on_processor(src_row));
    
      unsigned int local_src_row = src_row - SparseMatrix<T>::_global_offset;
      const std::map<unsigned int, T> & cols = _mat_local[local_src_row];
      for(typename std::map<unsigned int, T>::const_iterator it=cols.begin(); it!=cols.end(); it++)
      {
        unsigned int col = it->first;
        T value = it->second;
        add(dst_row, col, value);
      }
    }
  
    // sync _mat_nonlocal entries
    close(false);
    
    return;
  }
  
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
void PetscMatrix<T>::clear_row(int row, const T diag )
{
  if(_mat_buf_mode)
  {
    unsigned int local_row = row-SparseMatrix<T>::_global_offset;
    std::map<unsigned int, T> & cols = _mat_local[local_row];
    
    for(typename std::map<unsigned int, T>::iterator it=cols.begin(); it!=cols.end(); it++)
      it->second = 0.0;
    
    cols[row] = diag;
    return;
  }
    
  
#if PETSC_VERSION_GE(3,2,0)
  MatZeroRows(_mat, 1, &row, diag, PETSC_NULL, PETSC_NULL);
#else
  MatZeroRows(_mat, 1, &row, diag);
#endif  
  
}


template <typename T>
void PetscMatrix<T>::clear_row(const std::vector<int> &rows, const T diag)
{
  if(_mat_buf_mode)
  {
    for(unsigned int n=0; n<rows.size(); n++)
    {
      unsigned int local_row = rows[n]-SparseMatrix<T>::_global_offset;
      std::map<unsigned int, T> & cols = _mat_local[local_row];
    
      for(typename std::map<unsigned int, T>::iterator it=cols.begin(); it!=cols.end(); it++)
        it->second = 0.0;
    
      cols[rows[n]] = diag;
    }
    return;
  }
    
  
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



template <typename T>
void PetscMatrix<T>::flush_buf()
{
  genius_assert(_closed);
  genius_assert(_mat_buf_mode);
    
  std::vector<int> n_nz(SparseMatrix<T>::_m_local, 0);
  std::vector<int> n_oz(SparseMatrix<T>::_m_local, 0);
  
  for(size_t n=0; n<_mat_local.size(); ++n)
  {
    const std::map<unsigned int, T> & cols = _mat_local[n];
    int nz = 0;
    int noz = 0;
    for(typename std::map<unsigned int, T>::const_iterator it=cols.begin(); it!=cols.end(); it++)
    {
      unsigned int col = it->first;
      if( SparseMatrix<T>::col_on_processor(col) ) nz++;
      else noz++;
    }
    n_nz[n] = nz;
    n_oz[n] = noz;
  }
  
  int ierr     = 0;

  // create a sequential matrix on one processor
  if (Genius::n_processors()==1)
  {
    // alloc memory for sequence matrix here
    ierr = MatSeqAIJSetPreallocation(_mat, 0, &n_nz[0]); genius_assert(!ierr);
  }
  else
  {
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
  
  // we have to set this flag since preallocation may not be exact
  ierr = MatSetOption(_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); genius_assert(!ierr);
  
  // extra flag
  ierr = MatSetFromOptions(_mat); genius_assert(!ierr);
  
  // set value
  for(size_t n=0; n<_mat_local.size(); ++n)
  {
    unsigned int row = n+SparseMatrix<T>::_global_offset; 
    
    const std::map<unsigned int, T> & col_map = _mat_local[n];
    std::vector<unsigned int> cols;
    std::vector<T> col_values; 
    for(typename std::map<unsigned int, T>::const_iterator it=col_map.begin(); it!=col_map.end(); it++)
    {
      cols.push_back(it->first);
      col_values.push_back(it->second);
    }
    
    ierr = MatSetValues(_mat, 1, (int*) &row, cols.size(), (int*) &cols[0], &col_values[0], ADD_VALUES);
    genius_assert(!ierr);
  }
  
  ierr = MatAssemblyBegin (_mat, MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd   (_mat, MAT_FINAL_ASSEMBLY);
  genius_assert(!ierr);
  
  
  _mat_local.clear();
  _mat_buf_mode = false;
}

//------------------------------------------------------------------
// Explicit instantiations
template class PetscMatrix<PetscScalar>;


#endif // #ifdef HAVE_PETSC
