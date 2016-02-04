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
// SymbolicMatrix members
template <typename T>
void SymbolicMatrix<T>::init ()
{
  // on processor 
  _on_processor_entries.i.resize(SparseMatrix<T>::_m_local+1, 0);
  _on_processor_entries.imax.resize(SparseMatrix<T>::_m_local, _nz);
  _on_processor_entries.iocp.resize(SparseMatrix<T>::_m_local, 0);
  _on_processor_entries.j.resize(_nz*SparseMatrix<T>::_m_local, 0);

  for(unsigned int n=1; n<=SparseMatrix<T>::_m_local; n++)
    _on_processor_entries.i[n] = _on_processor_entries.i[n-1] + _on_processor_entries.imax[n-1];
  
  // off processor 
  _off_processor_entries.i.resize(SparseMatrix<T>::_m_local+1, 0);
  _off_processor_entries.imax.resize(SparseMatrix<T>::_m_local, _noz);
  _off_processor_entries.iocp.resize(SparseMatrix<T>::_m_local, 0);
  _off_processor_entries.j.resize(_noz*SparseMatrix<T>::_m_local, 0);

  for(unsigned int n=1; n<=SparseMatrix<T>::_m_local; n++)
    _off_processor_entries.i[n] = _off_processor_entries.i[n-1] + _off_processor_entries.imax[n-1];
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


  // flush _extra_on_processor_entries not empty or _on_processor_entries not full
  if((!_extra_on_processor_entries.empty()) || (!_on_processor_entries.is_full()))
  {
    aij_type aij;
    aij.imax = _on_processor_entries.iocp;
    for(typename set_type::const_iterator it =_extra_on_processor_entries.begin(); it !=_extra_on_processor_entries.end(); ++it)
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
      for(unsigned int j=0; j<_on_processor_entries.iocp[i]; j++)
      {
        unsigned int col = _on_processor_entries.j[_on_processor_entries.i[i]+j];
        _set_aij_entry(aij, row, col);
      }
    }

    for(typename set_type::const_iterator it =_extra_on_processor_entries.begin(); it !=_extra_on_processor_entries.end(); ++it)
    {
      unsigned int row = (*it).first;
      unsigned int col = (*it).second;
      _set_aij_entry(aij, row, col);
    }

    _on_processor_entries = aij;
    _extra_on_processor_entries.clear();
  }
  
  // flush _extra_off_processor_entries not empty or _off_processor_entries not full
  if((!_extra_off_processor_entries.empty()) || (!_off_processor_entries.is_full()))
  {
    aij_type aij;
    aij.imax = _off_processor_entries.iocp;
    for(typename set_type::const_iterator it =_extra_off_processor_entries.begin(); it !=_extra_off_processor_entries.end(); ++it)
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
      for(unsigned int j=0; j<_off_processor_entries.iocp[i]; j++)
      {
        unsigned int col = _off_processor_entries.j[_off_processor_entries.i[i]+j];
        _set_aij_entry(aij, row, col);
      }
    }

    for(typename set_type::const_iterator it =_extra_off_processor_entries.begin(); it !=_extra_off_processor_entries.end(); ++it)
    {
      unsigned int row = (*it).first;
      unsigned int col = (*it).second;
      _set_aij_entry(aij, row, col);
    }

    _off_processor_entries = aij;
    _extra_off_processor_entries.clear();
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

  // on processor entry
  if( SparseMatrix<T>::col_on_processor(j) )
  {
    _set_local_entry(_on_processor_entries, _extra_on_processor_entries, i, j);
  }
  // off processor entry
  else
  {
     _set_local_entry(_off_processor_entries, _extra_off_processor_entries, i, j);
  }
}


template <typename T>
void SymbolicMatrix<T>::_set_local_entry(aij_type & aij_entries, set_type & set_entries ,unsigned int i, unsigned int j)
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

  // in hash table
  if( set_entries.find(std::make_pair(i, j)) != set_entries.end() ) return;

  // aij array is full
  if(aij_entries.iocp[local_row] == aij_entries.imax[local_row])
  {
    set_entries.insert( std::make_pair(i, j) );
    _closed = false;
  }
  else
  {
    // shift up all the later entries in this row
    for (int ii=aij_entries.iocp[local_row]-1; ii>=col; ii--)
    {
      aij_entries.j[ii+1] = aij_entries.j[ii];
    }
    aij_entries.j[col] = j;
    aij_entries.iocp[local_row]++;
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
void SymbolicMatrix<T>::add_row (unsigned int row, const std::vector<unsigned int> &cols, const T* dm)
{
  for(unsigned int j=0; j<cols.size(); j++)
    _set_entry(row, cols[j]);
}

template <typename T>
void SymbolicMatrix<T>::add_row (unsigned int row, unsigned int n, const unsigned int * cols, const T* dm)
{
  for(unsigned int j=0; j<n; j++)
    _set_entry(row, cols[j]);
}


template <typename T>
void SymbolicMatrix<T>::add_row (unsigned int row, int n, const int * cols, const T* dm)
{
  for(int j=0; j<n; j++)
    _set_entry(row, cols[j]);
}


template <typename T>
void SymbolicMatrix<T>::get_row (unsigned int row, int n, const int * cols, T* dm)
{
  // do nothing
}


template <typename T>
void SymbolicMatrix<T>::add_matrix(const std::vector<unsigned int>& rows,
                                   const std::vector<unsigned int>& cols,
                                   const T*)
{
  for(unsigned int i=0; i<rows.size(); i++)
    for(unsigned int j=0; j<cols.size(); j++)
      _set_entry(rows[i], cols[j]);
}


template <typename T>
void  SymbolicMatrix<T>::add_matrix (unsigned int m, unsigned int * rows,
                                     unsigned int n, unsigned int * cols,
                                     const T*)
{
  for(unsigned int i=0; i<m; i++)
    for(unsigned int j=0; j<n; j++)
      _set_entry(rows[i], cols[j]);
}


template <typename T>
void SymbolicMatrix<T>::add_row_to_row(const std::vector<int> &src_rows,
                                       const std::vector<int> &dst_rows)
{
  genius_assert(_closed);

  for(unsigned int n=0; n<src_rows.size(); n++)
  {
    unsigned int src_row = static_cast<unsigned int>(src_rows[n]);
    unsigned int dst_row = static_cast<unsigned int>(dst_rows[n]);

    genius_assert(SparseMatrix<T>::row_on_processor(src_row));
    
    // on processor
    int local_src_row = src_row - SparseMatrix<T>::_global_offset;
    for(unsigned int j=0; j<_on_processor_entries.iocp[local_src_row]; j++)
    {
      unsigned int src_col = _on_processor_entries.j[_on_processor_entries.i[local_src_row]+j];
      _set_entry(dst_row, src_col);
    }
    
    // off processor
    for(unsigned int j=0; j<_off_processor_entries.iocp[local_src_row]; j++)
    {
      unsigned int src_col = _off_processor_entries.j[_off_processor_entries.i[local_src_row]+j];
      _set_entry(dst_row, src_col);
    }
  }
  
  close(false);
}





template <typename T>
void SymbolicMatrix<T>::get_sparsity_pattern (std::vector<int> &n_nz, std::vector<int> &n_oz)
{
  genius_assert(_closed);

  n_nz.resize(SparseMatrix<T>::_m_local);
  n_oz.resize(SparseMatrix<T>::_m_local);

  std::fill(n_nz.begin(), n_nz.end(), 0);
  std::fill(n_oz.begin(), n_oz.end(), 0);

  for(unsigned int i = 0; i < SparseMatrix<T>::_m_local; i++)
  {
    n_nz[i] += _on_processor_entries.iocp[i];
    n_oz[i] += _off_processor_entries.iocp[i];
  }
}



template <typename T>
void SymbolicMatrix<T>::get_nonzero_pattern (std::vector< std::vector<int> > & nz)
{
  genius_assert(_closed);
  

  nz.resize(SparseMatrix<T>::_m_local);
  
  // on processor
  for(unsigned int i = 0; i < SparseMatrix<T>::_m_local; i++)
  {
    unsigned int col_offset = _on_processor_entries.i[i];
    unsigned int col_num    = _on_processor_entries.iocp[i];
    for(unsigned int j=0; j<col_num; j++)
      nz[i].push_back(_on_processor_entries.j[col_offset+j]);
  }
  
  // off processor
  for(unsigned int i = 0; i < SparseMatrix<T>::_m_local; i++)
  {
    unsigned int col_offset = _off_processor_entries.i[i];
    unsigned int col_num    = _off_processor_entries.iocp[i];
    for(unsigned int j=0; j<col_num; j++)
      nz[i].push_back(_off_processor_entries.j[col_offset+j]);
  }
}



template <typename T>
size_t SymbolicMatrix<T>::_build_hash_entry() const
{
  genius_assert(_closed);

  unsigned int hash_value = 0;

 
  for(unsigned int i = 0; i < SparseMatrix<T>::_m_local; i++)
  {
     // local aij
    unsigned int row = i + SparseMatrix<T>::_global_offset;
    for(unsigned int j=0; j<_on_processor_entries.iocp[i]; j++)
    {
      unsigned int col = _on_processor_entries.j[_on_processor_entries.i[i]+j];
      hash_value += key_hash()(std::make_pair(row, col));
    }
    
    // off processor
    for(unsigned int j=0; j<_off_processor_entries.iocp[i]; j++)
    {
      unsigned int col = _off_processor_entries.j[_off_processor_entries.i[i]+j];
      hash_value += key_hash()(std::make_pair(row, col));
    }
  }

  Parallel::sum(hash_value);

  return static_cast<size_t>(hash_value);
}





//------------------------------------------------------------------
// Explicit instantiations
template class SymbolicMatrix<PetscScalar>;


