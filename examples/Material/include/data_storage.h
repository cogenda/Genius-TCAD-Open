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

#ifndef __data_storage_h__
#define __data_storage_h__

#include <vector>
#include "enum_data_type.h"
#include "vector_value.h"
#include "tensor_value.h"

class DataStorage
{
public:
  DataStorage(unsigned int size=0):_size(size), _reserve_size(0) {}

  ~DataStorage() {}

  /**
   * reserve memory for data array
   */
  void reserve(unsigned int reserve_size)
  {
    _reserve_size = reserve_size;

    for(unsigned int n=0; n<_scalar_fill.size(); ++n)
      if( _scalar_fill[n] )
        _scalar_block[n].reserve(reserve_size);

    for(unsigned int n=0; n<_complex_fill.size(); ++n)
      if( _complex_fill[n] )
        _complex_block[n].reserve(reserve_size);

    for(unsigned int n=0; n<_vector_fill.size(); ++n)
      if( _vector_fill[n] )
        _vector_block[n].reserve(reserve_size);

    for(unsigned int n=0; n<_tensor_fill.size(); ++n)
      if( _tensor_fill[n] )
        _tensor_block[n].reserve(reserve_size);
  }

  /**
   * clear all the data
   */
  void clear()
  {
    _size = 0;

    std::fill(_scalar_fill.begin(), _scalar_fill.end(), false);
    for(unsigned int n=0; n<_scalar_fill.size(); ++n)
      _scalar_block[n].clear();

    std::fill(_complex_fill.begin(), _complex_fill.end(), false);
    for(unsigned int n=0; n<_complex_fill.size(); ++n)
      _complex_block[n].clear();

    std::fill(_vector_fill.begin(), _vector_fill.end(), false);
    for(unsigned int n=0; n<_vector_fill.size(); ++n)
      _vector_block[n].clear();

    std::fill(_tensor_fill.begin(), _tensor_fill.end(), false);
    for(unsigned int n=0; n<_tensor_fill.size(); ++n)
      _tensor_block[n].clear();

  }

  /**
   * increase the size of data array
   */
  unsigned int increase(unsigned int size=1)
  {
    _size += size;

    for(unsigned int n=0; n<_scalar_fill.size(); ++n)
      if( _scalar_fill[n] )
        _scalar_block[n].resize(_size);

    for(unsigned int n=0; n<_complex_fill.size(); ++n)
      if( _complex_fill[n] )
        _complex_block[n].resize(_size);

    for(unsigned int n=0; n<_vector_fill.size(); ++n)
      if( _vector_fill[n] )
        _vector_block[n].resize(_size);

    for(unsigned int n=0; n<_tensor_fill.size(); ++n)
      if( _tensor_fill[n] )
        _tensor_block[n].resize(_size);

    return _size;
  }

  /**
   * @return the current size of data array
   */
  unsigned int size() const
  { return _size; }

  /**
   * allocate variables
   * the total variable number equals to the size of bool vector
   * when the bool value in the given vector is true, a block data will be allocated
   */
  void allocate_scalar_variable(const std::vector<bool> & flags)
  {
    _scalar_fill = flags;
    _scalar_block.resize(flags.size());
    for(unsigned int n=0; n<flags.size(); ++n)
      if(flags[n]) _scalar_block[n].resize(_size, 0.0);
  }

  void allocate_complex_variable(const std::vector<bool> & flags)
  {
    _complex_fill = flags;
    _complex_block.resize(flags.size());
    for(unsigned int n=0; n<flags.size(); ++n)
      if(flags[n]) _complex_block[n].resize(_size);
  }

  void allocate_vector_variable(const std::vector<bool> & flags)
  {
    _vector_fill = flags;
    _vector_block.resize(flags.size());
    for(unsigned int n=0; n<flags.size(); ++n)
      if(flags[n]) _vector_block[n].resize(_size);
  }

  void allocate_tensor_variable(const std::vector<bool> & flags)
  {
    _tensor_fill = flags;
    _tensor_block.resize(flags.size());
    for(unsigned int n=0; n<flags.size(); ++n)
      if(flags[n]) _tensor_block[n].resize(_size);
  }

  /**
   * add a data block for ith variable.
   * when i equals to static_cast<unsigned int>(-1), append to the variable list
   * when flag is true, memory is allocated
   * @return actual variable index
   */
  unsigned int add_scalar_variable(unsigned int v, bool flag)
  {
    if( v == static_cast<unsigned int>(-1) )
      v =  _scalar_fill.size();

    if( v + 1 > _scalar_fill.size())
    {
      _scalar_fill.resize(v+1, false);
      _scalar_block.resize(_scalar_fill.size());
    }

    _scalar_fill[v] = flag;
    if(_scalar_fill[v])
    {
      _scalar_block[v].reserve(_reserve_size);
      _scalar_block[v].resize(_size);
    }

    return v;
  }


  unsigned int add_complex_variable(unsigned int v, bool flag)
  {
    if( v == static_cast<unsigned int>(-1) )
      v =  _complex_fill.size();

    if( v + 1 > _complex_fill.size())
    {
      _complex_fill.resize(v+1, false);
      _complex_block.resize(_complex_fill.size());
    }

    _complex_fill[v] = flag;
    if(_complex_fill[v])
    {
      _complex_block[v].reserve(_reserve_size);
      _complex_block[v].resize(_size);
    }

    return v;
  }

  unsigned int add_vector_variable(unsigned int v, bool flag)
  {
    if( v == static_cast<unsigned int>(-1) )
      v =  _vector_fill.size();

    if( v + 1 > _vector_fill.size())
    {
      _vector_fill.resize(v+1, false);
      _vector_block.resize(_vector_fill.size());
    }

    _vector_fill[v] = flag;
    if(_vector_fill[v])
    {
      _vector_block[v].reserve(_reserve_size);
      _vector_block[v].resize(_size);
    }

    return v;
  }

  unsigned int add_tensor_variable(unsigned int v, bool flag)
  {
    if( v == static_cast<unsigned int>(-1) )
      v =  _tensor_fill.size();

    if( v + 1 > _tensor_fill.size())
    {
      _tensor_fill.resize(v+1, false);
      _tensor_block.resize(_tensor_fill.size());
    }

    _tensor_fill[v] = flag;
    if(_tensor_fill[v])
    {
      _tensor_block[v].reserve(_reserve_size);
      _tensor_block[v].resize(_size);
    }

    return v;
  }

  /**
   * @return the defined scaler variable number
   */
  unsigned int n_scalar() const
    { return _scalar_fill.size(); }

  /**
   * @return the defined complex variable number
   */
  unsigned int n_complex() const
    { return _complex_fill.size(); }

  /**
   * @return the defined vector variable number
   */
  unsigned int n_vector() const
    { return _vector_fill.size(); }

  /**
   * @return the defined tensor variable number
   */
  unsigned int n_tensor() const
    { return _tensor_fill.size(); }

  /**
   * data access function
   */
  const PetscScalar & scalar(const unsigned int v, const unsigned int offset) const
    { return _scalar_block[v][offset]; }

  /**
   * data access function
   */
  PetscScalar & scalar(const unsigned int v, const unsigned int offset)
  { return _scalar_block[v][offset]; }

  /**
   * data access function
   */
  const std::complex<PetscScalar> & complex(const unsigned int v, const unsigned int offset) const
    { return _complex_block[v][offset]; }

  /**
   * data access function
   */
  std::complex<PetscScalar> & complex(const unsigned int v, const unsigned int offset)
  { return _complex_block[v][offset]; }

  /**
   * data access function
   */
  const VectorValue<PetscScalar> & vector(const unsigned int v, const unsigned int offset) const
    { return _vector_block[v][offset]; }

  /**
   * data access function
   */
  VectorValue<PetscScalar> & vector(const unsigned int v, const unsigned int offset)
  { return _vector_block[v][offset]; }

  /**
   * data access function
   */
  const TensorValue<PetscScalar> & tensor(const unsigned int v, const unsigned int offset) const
    { return _tensor_block[v][offset]; }

  /**
   * data access function
   */
  TensorValue<PetscScalar> & tensor(const unsigned int v, const unsigned int offset)
  { return _tensor_block[v][offset]; }


  /**
   * universal data access function via template
   */
  template <typename T>
  T & data(const unsigned int , const unsigned int );

  /**
   * universal data access function via template
   */
  template <typename T>
  const T & data(const unsigned int , const unsigned int ) const;

  /**
   * approx memory usage
   */
  size_t memory_size() const
  {
    size_t counter = sizeof(*this);

    for(unsigned int n=0; n<_scalar_block.size(); ++n)
      counter += _scalar_block[n].capacity()*sizeof(PetscScalar);

    for(unsigned int n=0; n<_complex_block.size(); ++n)
      counter += _complex_block[n].capacity()*sizeof(std::complex<PetscScalar>);

    for(unsigned int n=0; n<_vector_block.size(); ++n)
      counter += _vector_block[n].capacity()*sizeof(VectorValue<PetscScalar>);

    for(unsigned int n=0; n<_tensor_block.size(); ++n)
      counter += _tensor_block[n].capacity()*sizeof(TensorValue<PetscScalar>);

    return counter;
  }

private:

  /**
   * the size of data array
   */
  unsigned int _size;

  /**
   * the reserved size of data array
   */
  unsigned int _reserve_size;

  /**
   * indicator of scalar value
   */
  std::vector<bool> _scalar_fill;
  std::vector< std::vector<PetscScalar> > _scalar_block;

  /**
   * indicator of complex value
   */
  std::vector<bool> _complex_fill;
  std::vector< std::vector<std::complex<PetscScalar> > > _complex_block;

  /**
   * indicator of vector value
   */
  std::vector<bool> _vector_fill;
  std::vector< std::vector<VectorValue<PetscScalar> > > _vector_block;

  /**
   * indicator of tensor value
   */
  std::vector<bool> _tensor_fill;
  std::vector< std::vector<TensorValue<PetscScalar> > > _tensor_block;

};


//explicit instantiation
template <typename T>
inline T & DataStorage::data(const unsigned int , const unsigned int ) { }

template<>
inline PetscScalar & DataStorage::data< PetscScalar >(const unsigned int v, const unsigned int offset)  {  return _scalar_block[v][offset]; }

template<>
inline std::complex<PetscScalar> & DataStorage::data< std::complex<PetscScalar> >(const unsigned int v, const unsigned int offset) { return _complex_block[v][offset]; }

template<>
inline VectorValue<PetscScalar> & DataStorage::data< VectorValue<PetscScalar> >(const unsigned int v, const unsigned int offset) { return _vector_block[v][offset]; }

template<>
inline TensorValue<PetscScalar> & DataStorage::data< TensorValue<PetscScalar> >(const unsigned int v, const unsigned int offset) { return _tensor_block[v][offset]; }



template <typename T>
inline const T & DataStorage::data(const unsigned int , const unsigned int ) const { }

template<>
inline const PetscScalar & DataStorage::data< PetscScalar >(const unsigned int v, const unsigned int offset) const { return _scalar_block[v][offset]; }

template<>
inline const std::complex<PetscScalar> & DataStorage::data< std::complex<PetscScalar> >(const unsigned int v, const unsigned int offset) const { return _complex_block[v][offset]; }

template<>
inline const VectorValue<PetscScalar> & DataStorage::data< VectorValue<PetscScalar> >(const unsigned int v, const unsigned int offset) const { return _vector_block[v][offset]; }

template<>
inline const TensorValue<PetscScalar> & DataStorage::data< TensorValue<PetscScalar> >(const unsigned int v, const unsigned int offset) const { return _tensor_block[v][offset]; }


#endif

