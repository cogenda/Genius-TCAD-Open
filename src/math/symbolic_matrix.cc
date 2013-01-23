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


// C++ includes
#include "symbolic_matrix.h"
#include "parallel.h"


//-----------------------------------------------------------------------
// PetscMatrix members
template <typename T>
void SymbolicMatrix<T>::init ()
{
  _local_entries.i.resize(SparseMatrix<T>::_m_local+1, 0);
  _local_entries.imax.resize(SparseMatrix<T>::_m_local, _nz);
  _local_entries.iocp.resize(SparseMatrix<T>::_m_local, 0);
  _local_entries.j.resize(_nz*SparseMatrix<T>::_m_local, 0);

  for(unsigned int n=1; n<=SparseMatrix<T>::_m_local; n++)
    _local_entries.i[n] = _local_entries.i[n-1] + _local_entries.imax[n-1];
}


template <typename T>
void SymbolicMatrix<T>::close (bool final)
{
  // flush nonlocal entries
  unsigned int nonlocal_entries = _nonlocal_entries.size();
  Parallel::sum(nonlocal_entries);
  if(nonlocal_entries)
  {
    std::vector<unsigned int> rows;
    std::vector<unsigned int> cols;
    for(typename set_type::const_iterator it =_nonlocal_entries.begin(); it !=_nonlocal_entries.end(); ++it)
    {
      rows.push_back((*it).first);
      cols.push_back((*it).second);
    }
    Parallel::allgather(rows);
    Parallel::allgather(cols);

    _nonlocal_entries.clear();

    for(unsigned int n=0; n<rows.size(); ++n)
    {
      unsigned int row = rows[n];
      if( !SparseMatrix<T>::row_on_processor(row) ) continue;

      unsigned int col = cols[n];
      _set_entry(row, col);
    }
  }
  genius_assert(_nonlocal_entries.empty());


  // flush _extra_local_entries
  if(_extra_local_entries.empty())
  {
    aij_type aij;
    aij.imax = _local_entries.iocp;
    for(typename set_type::const_iterator it =_extra_local_entries.begin(); it !=_extra_local_entries.end(); ++it)
    {
      unsigned int row = (*it).first;
      int local_row = row - SparseMatrix<T>::_global_offset;
      aij.imax[local_row]++;
    }

    aij.i.resize(SparseMatrix<T>::_m_local+1, 0);
    for(unsigned int n=1; n<=SparseMatrix<T>::_m_local; n++)
      aij.i[n] = aij.i[n-1] + aij.imax[n-1];

    aij.iocp.resize(SparseMatrix<T>::_m_local, 0);
    aij.j.resize(aij.i.back(), 0);

    for(unsigned int i = 0; i < SparseMatrix<T>::_m_local; i++)
    {
      unsigned int row = i + SparseMatrix<T>::_global_offset;
      for(unsigned int j=0; j<_local_entries.iocp[i]; j++)
      {
        unsigned int col = _local_entries.j[_local_entries.i[i]+j];
        _set_aij_entry(aij, row, col);
      }
    }

    for(typename set_type::const_iterator it =_extra_local_entries.begin(); it !=_extra_local_entries.end(); ++it)
    {
      unsigned int row = (*it).first;
      unsigned int col = (*it).second;
      _set_aij_entry(aij, row, col);
    }

    _local_entries = aij;
    _extra_local_entries.clear();
  }

  if(final)
  {
    size_t new_hash_entry = _build_hash_entry();
    _same_pattern = (new_hash_entry == _hash_entry);
    _hash_entry = new_hash_entry;
  }

  _closed = true;
}





template <typename T>
void SymbolicMatrix<T>::clear ()
{

}




template <typename T>
void SymbolicMatrix<T>::_set_entry (const unsigned int i, const unsigned int j)
{
  // not a local entry, should be processed in close() function
  if( !SparseMatrix<T>::row_on_processor(i) )
  {
    _nonlocal_entries.insert( std::make_pair(i, j) );
    _closed = false;
    return;
  }

  // off processor entry
  if( !SparseMatrix<T>::col_on_processor(j) )
  {
    _off_processor_entries.insert( std::make_pair(i, j) );
  }
  // on processor entry
  else
  {
    int local_row = i - SparseMatrix<T>::_global_offset;
    int local_col = j - SparseMatrix<T>::_global_offset;

    unsigned int col_begin = _local_entries.i[local_row];
    unsigned int col_end   = col_begin + _local_entries.iocp[local_row];
    unsigned int col=0;
    while (col_end-col_begin > 5)
    {
      unsigned int t = (col_begin+col_end)/2;
      if (_local_entries.j[t] > j)
        col_end = t;
      else
        col_begin  = t;
    }

    for (col=col_begin; col<col_end; col++)
    {
      if (_local_entries.j[col] > j) break;
      if (_local_entries.j[col] == j) return; // in aij array
    }

    // in hash table
    if( _extra_local_entries.find(std::make_pair(i, j)) != _extra_local_entries.end() ) return;

    // aij array is full
    if(_local_entries.iocp[local_row] == _local_entries.imax[local_row])
    {
      _extra_local_entries.insert( std::make_pair(i, j) );
      _closed = false;
    }
    else
    {
      // shift up all the later entries in this row
      for (int ii=_local_entries.iocp[local_row]-1; ii>=col; ii--)
      {
        _local_entries.j[ii+1] = _local_entries.j[ii];
      }
      _local_entries.j[col] = j;
      _local_entries.iocp[local_row]++;
    }
  }

}


template <typename T>
void SymbolicMatrix<T>::_set_aij_entry(aij_type & aij_entries, unsigned int i, unsigned int j)
{

  int local_row = i - SparseMatrix<T>::_global_offset;
  int local_col = j - SparseMatrix<T>::_global_offset;

  unsigned int col_begin = aij_entries.i[local_row];
  unsigned int col_end   = col_begin + aij_entries.iocp[local_row];
  unsigned int col=0;
  while (col_end-col_begin > 5)
  {
    unsigned int t = (col_begin+col_end)/2;
    if (aij_entries.j[t] > j)
      col_end = t;
    else
      col_begin  = t;
  }

  for (col=col_begin; col<col_end; col++)
  {
    if (aij_entries.j[col] > j) break;
    if (aij_entries.j[col] == j) return; // in aij array
  }

  // aij array should not full
  genius_assert(aij_entries.iocp[local_row] < aij_entries.imax[local_row]);

  // shift up all the later entries in this row
  for (int ii=aij_entries.iocp[local_row]-1; ii>=col; ii--)
  {
    aij_entries.j[ii+1] = aij_entries.j[ii];
  }

  aij_entries.j[col] = j;
  aij_entries.iocp[local_row]++;
}




template <typename T>
void SymbolicMatrix<T>::set (const unsigned int i, const unsigned int j, const T)
{
  _set_entry(i,j);
}



template <typename T>
void SymbolicMatrix<T>::add (const unsigned int i, const unsigned int j, const T)
{
  _set_entry(i,j);
}



template <typename T>
void SymbolicMatrix<T>::add_matrix(const DenseMatrix<T>& ,
                                   const std::vector<unsigned int>& rows,
                                   const std::vector<unsigned int>& cols)
{
  for(unsigned int i=0; i<rows.size(); i++)
    for(unsigned int j=0; j<cols.size(); j++)
      _set_entry(rows[i], cols[j]);
}


template <typename T>
void  SymbolicMatrix<T>::add_matrix (const T* ,
                                     unsigned int m, unsigned int * rows,
                                     unsigned int n, unsigned int * cols)
{
  for(unsigned int i=0; i<m; i++)
    for(unsigned int j=0; j<n; j++)
      _set_entry(rows[i], cols[j]);
}


template <typename T>
void SymbolicMatrix<T>::add_row_to_row(const std::vector<unsigned int> &src_rows,
                                       const std::vector<unsigned int> &dst_rows)
{
  genius_assert(_closed);

  for(unsigned int n=0; n<src_rows.size(); n++)
  {
    unsigned int src_row = src_rows[n];
    unsigned int dst_row = dst_rows[n];

    genius_assert(SparseMatrix<T>::row_on_processor(src_row));
    int local_src_row = src_row - SparseMatrix<T>::_global_offset;
    for(unsigned int j=0; j<_local_entries.iocp[local_src_row]; j++)
    {
      unsigned int src_col = _local_entries.j[_local_entries.i[local_src_row]+j];
      _set_entry(dst_row, src_col);
    }
  }
}





template <typename T>
void SymbolicMatrix<T>::get_sparsity_pattern (std::vector<int> &n_nz, std::vector<int> &n_oz)
{
  genius_assert(_closed);

  n_nz.resize(SparseMatrix<T>::_m_local);
  n_oz.resize(SparseMatrix<T>::_m_local);

  std::fill(n_nz.begin(), n_nz.end(), 0);
  std::fill(n_oz.begin(), n_oz.end(), 0);

  // local aij
  for(unsigned int i = 0; i < SparseMatrix<T>::_m_local; i++)
  {
    n_nz[i] += _local_entries.iocp[i];
  }

  // off processor
  for(typename set_type::const_iterator it =_off_processor_entries.begin(); it !=_off_processor_entries.end(); ++it)
  {
    unsigned int row = (*it).first;
    int local_row = row - SparseMatrix<T>::_global_offset;
    n_oz[local_row]++;
  }
}



template <typename T>
size_t SymbolicMatrix<T>::_build_hash_entry() const
{
  genius_assert(_closed);

  unsigned int hash_value = 0;

  // local aij
  for(unsigned int i = 0; i < SparseMatrix<T>::_m_local; i++)
  {
    unsigned int row = i + SparseMatrix<T>::_global_offset;
    for(unsigned int j=0; j<_local_entries.iocp[i]; j++)
    {
      unsigned int col = _local_entries.j[_local_entries.i[i]+j];
      hash_value += key_hash()(std::make_pair(row, col));
    }
  }

  // off processor
  for(typename set_type::const_iterator it =_off_processor_entries.begin(); it !=_off_processor_entries.end(); ++it)
  {
    hash_value += key_hash()(*it);
  }

  Parallel::sum(hash_value);

  return static_cast<size_t>(hash_value);
}





//------------------------------------------------------------------
// Explicit instantiations
template class SymbolicMatrix<Real>;


