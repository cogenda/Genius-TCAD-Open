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


// Local includes
#include "csr_matrix.h"
#include "parallel.h"



//-----------------------------------------------------------------------
// CSRMatrix members

//-----------------------------------------------------------------------
// CSRMatrix inline members
template <typename T>
CSRMatrix<T>::CSRMatrix(const unsigned int m,   const unsigned int n,
                        const unsigned int m_l, const unsigned int n_l)
  : SparseMatrix<T>(m,n,m_l,n_l), _add_value_flag(true), _closed(false), _destroy_mat_on_exit(true)
{
  if ((m==0) || (n==0))
    return;

  _mat.resize(m_l);

  SparseMatrix<T>::_is_initialized = true;
}


template <typename T>
CSRMatrix<T>::~CSRMatrix()
{
  this->clear();
}




template <typename T>
void CSRMatrix<T>::set (const unsigned int i, const unsigned int j, const T value)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==false );
  if( SparseMatrix<T>::row_on_processor(i) )
    _mat[i-SparseMatrix<T>::_global_offset][j] = value;
  else 
  {
    _mat_nonlocal[std::make_pair(i,j)] = value;
    _closed = false;
  }
}



template <typename T>
void CSRMatrix<T>::add (const unsigned int i, const unsigned int j, const T value)
{
  genius_assert (this->initialized());
  genius_assert (_add_value_flag==true );
  if( SparseMatrix<T>::row_on_processor(i) )
    _mat[i-SparseMatrix<T>::_global_offset][j] += value;
  else 
  {
    _mat_nonlocal[std::make_pair(i,j)] += value;
    _closed = false;
  }
}

template <typename T>
void CSRMatrix<T>::add_row (unsigned int row, const std::vector<unsigned int> &cols, const T* dm)
{
   unsigned int n = cols.size();
   
   for(unsigned int j=0; j<n; j++)
     add(row, cols[j], dm[j]);
}


template <typename T>
void CSRMatrix<T>::add_row (unsigned int row, unsigned int n, const unsigned int * cols, const T* dm)
{
  for(int j=0; j<n; j++)
     add(row, cols[j], dm[j]);
}


template <typename T>
void CSRMatrix<T>::add_row (unsigned int row,  int n, const int * cols, const T* dm)
{
  for(int j=0; j<n; j++)
     add(row, cols[j], dm[j]);
}


template <typename T>
void CSRMatrix<T>::add_matrix(const std::vector<unsigned int>& rows,
                              const std::vector<unsigned int>& cols,
                              const T* dm)
{
   unsigned int m = rows.size();
   unsigned int n = cols.size();
   
   for(unsigned int i=0; i<m; i++)
    for(unsigned int j=0; j<n; j++)
      add(rows[i], cols[j], dm[i*n+j]);
}



template <typename T>
void CSRMatrix<T>::add_matrix (unsigned int m, unsigned int * rows,
                               unsigned int n, unsigned int * cols,
                               const T* dm)
{
  for(unsigned int i=0; i<m; i++)
    for(unsigned int j=0; j<n; j++)
      add(rows[i], cols[j], dm[i*n+j]);
}




template <typename T>
bool CSRMatrix<T>::closed() const
{
  genius_assert (this->initialized());
  return _closed;
}




template <typename T>
void CSRMatrix<T>::init ()
{

}


template <typename T>
void CSRMatrix<T>::close (bool )
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
      _mat[row-SparseMatrix<T>::_global_offset][col] += value;
    }
  }
  
  _closed = true;
}



template <typename T>
void CSRMatrix<T>::zero ()
{
  genius_assert(_closed);
  
  for(size_t n=0; n<_mat.size(); ++n)
  {
    std::map<unsigned int, T> & cols = _mat[n];
    for(typename std::map<unsigned int, T>::iterator it=cols.begin(); it!=cols.end(); it++)
      it->second = 0.0;
  }
}



template <typename T>
void CSRMatrix<T>::clear ()
{
  _mat.clear();
  _mat_nonlocal.clear();
}


template <typename T>
void CSRMatrix<T>::get_sparsity_pattern (std::vector<int> &n_nz, std::vector<int> &n_oz)
{
  genius_assert(_closed);
  
  n_nz.resize(SparseMatrix<T>::_m_local);
  n_oz.resize(SparseMatrix<T>::_m_local);

  std::fill(n_nz.begin(), n_nz.end(), 0);
  std::fill(n_oz.begin(), n_oz.end(), 0);
  
  for(size_t n=0; n<_mat.size(); ++n)
  {
    const std::map<unsigned int, T> & cols = _mat[n];
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
}

  
template <typename T>
void CSRMatrix<T>::get_nonzero_pattern (std::vector< std::vector<int> > & nz)
{
  genius_assert(_closed);
    
  nz.resize(SparseMatrix<T>::_m_local);
  
  for(size_t n=0; n<_mat.size(); ++n)
  {
    const std::map<unsigned int, T> & cols = _mat[n];
    for(typename std::map<unsigned int, T>::const_iterator it=cols.begin(); it!=cols.end(); it++)
    {
      unsigned int col = it->first;
      nz[n].push_back(col);
    }
  }
}

template <typename T>
void CSRMatrix<T>::get_row (unsigned int row, int n, const int * cols, T* dm)
{
  const std::map<unsigned int, T> & buf = _mat[row-SparseMatrix<T>::_global_offset];
  for(int i=0; i<n; i++)
  {
    dm[i] = buf.find(cols[i])->second;
  }
}


template <typename T>
T CSRMatrix<T>::operator () (const unsigned int i, const unsigned int j) const
{
  genius_assert (this->initialized());
  genius_assert ( i-SparseMatrix<T>::_global_offset  < SparseMatrix<T>::_m_local);
  
  const std::map<unsigned int, T> & buf = _mat[i-SparseMatrix<T>::_global_offset];
  
  typename std::map<unsigned int, T>::const_iterator ent =  buf.find(j);
  
  if(ent != buf.end() ) return ent->second;

  // Otherwise the entry is not in the sparse matrix,
  // i.e. it is 0.
  return 0.0;
}







template <typename T>
void CSRMatrix<T>::add_row_to_row(const std::vector<int> &src_rows,
                                  const std::vector<int> &dst_rows)
{
  genius_assert(_closed);

  for(unsigned int n=0; n<src_rows.size(); n++)
  {
    unsigned int src_row = static_cast<unsigned int>(src_rows[n]);
    unsigned int dst_row = static_cast<unsigned int>(dst_rows[n]);

    genius_assert(SparseMatrix<T>::row_on_processor(src_row));
    
    unsigned  local_src_row = src_row - SparseMatrix<T>::_global_offset;
    const std::map<unsigned int, T> & cols = _mat[local_src_row];
    for(typename std::map<unsigned int, T>::const_iterator it=cols.begin(); it!=cols.end(); it++)
    {
      unsigned int col = it->first;
      T value = it->second;
      add(dst_row, col, value);
    }
  }
  
  close(false);
}


template <typename T>
void CSRMatrix<T>::clear_row(int row, const T diag )
{
  unsigned int local_row = row-SparseMatrix<T>::_global_offset;
  std::map<unsigned int, T> & cols = _mat[local_row];
    
  for(typename std::map<unsigned int, T>::iterator it=cols.begin(); it!=cols.end(); it++)
    it->second = 0.0;
    
  cols[row] = diag;
}


template <typename T>
void CSRMatrix<T>::clear_row(const std::vector<int> &rows, const T diag)
{
  for(unsigned int n=0; n<rows.size(); n++)
  {
    unsigned int local_row = rows[n]-SparseMatrix<T>::_global_offset;
    std::map<unsigned int, T> & cols = _mat[local_row];
    
    for(typename std::map<unsigned int, T>::iterator it=cols.begin(); it!=cols.end(); it++)
      it->second = 0.0;
    
    cols[rows[n]] = diag;
  }
}



template <typename T>
void CSRMatrix<T>::print_personal(std::ostream& os) const
{

}


template <typename T>
void CSRMatrix<T>::print_matlab (const std::string name) const
{

}




//------------------------------------------------------------------
// Explicit instantiations
template class CSRMatrix<PetscScalar>;
