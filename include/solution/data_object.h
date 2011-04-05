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


#ifndef __data_object_h__
#define __data_object_h__

#include <complex>
#include <vector>
#include <map>

#include "data_storage.h"
#include "variable_define.h"

/**
 * The \p DataObject defines an abstract base class for data access interface
 * it is light weighted since data is storaged in class \p DataStorage
 */
class DataObject
{

protected:

  /**
   * Constructor. Protected so that you can't instantiate one of these.
   */
  DataObject(DataStorage * data_storage, const std::map<std::string, SimulationVariable> & variables)
      :_data_storage(data_storage), _variables(variables)
  {
    _offset = data_storage->increase(1) - 1;
  }

public:

  /**
   * Copy-constructor.
   */
  //DataObject (const DataObject&);


  /**
   * destruction
   */
  virtual ~DataObject()
  {}

  /**
   * @return the pointer to data storage class
   */
  DataStorage * data_storage() const
    { return _data_storage; }

  /**
   * @return the offset of this data object in data storage
   */
  unsigned int offset() const
    { return _offset; }


  /**
   * universal data access function by variable index via template
   */
  template <typename T>
  T & data(const unsigned int);


  /**
   * universal data access function by variable index via template
   */
  template <typename T>
  const T & data(const unsigned int) const;


  /**
   * @return true when variable exist
   */
  bool has_variable(const std::string &v, DataType type) const
  {
    if( _variables.find(v) == _variables.end() ) return false;
    return (_variables.find(v)->second.variable_data_type == type);
  }


  /**
   * @return variable index when variable exist. else return static_cast<unsigned int>(-1)
   */
  unsigned int variable_index(const std::string &v) const
  {
    if( _variables.find(v) == _variables.end() ) return static_cast<unsigned int>(-1);
    return _variables.find(v)->second.variable_index;
  }

  /**
   * universal data access function by variable name via template
   */
  template <typename T>
  T & data(const std::string &);

  /**
   * universal data access function by variable name via template
   */
  template <typename T>
  const T & data(const std::string &) const;


protected:

  /**
   * pointer to data storage class
   */
  DataStorage * _data_storage;

  /**
   * const reference to variables
   */
  const std::map<std::string, SimulationVariable>  & _variables;

  /**
   * the offset in data array
   */
  unsigned int _offset;

};


template <typename T>
inline T & data(const unsigned int ) {  }

template<>
inline Real & DataObject::data< Real >(const unsigned int v)  {  return _data_storage->scalar(v, _offset); }

template<>
inline std::complex<Real> & DataObject::data< std::complex<Real> >(const unsigned int v) { return _data_storage->complex(v, _offset); }

template<>
inline VectorValue<Real> & DataObject::data< VectorValue<Real> >(const unsigned int v) { return _data_storage->vector(v, _offset); }

template<>
inline TensorValue<Real> & DataObject::data< TensorValue<Real> >(const unsigned int v) { return _data_storage->tensor(v, _offset); }


template <typename T>
inline const T & DataObject::data(const unsigned int ) const { }

template<>
inline const Real & DataObject::data< Real >(const unsigned int v) const { return _data_storage->scalar(v, _offset); }

template<>
inline const std::complex<Real> & DataObject::data< std::complex<Real> >(const unsigned int v) const { return _data_storage->complex(v, _offset); }

template<>
inline const VectorValue<Real> & DataObject::data< VectorValue<Real> >(const unsigned int v) const { return _data_storage->vector(v, _offset); }

template<>
inline const TensorValue<Real> & DataObject::data< TensorValue<Real> >(const unsigned int v) const { return _data_storage->tensor(v, _offset); }


template <typename T>
inline T & DataObject::data(const std::string &) { }

template<>
inline Real & DataObject::data< Real >(const std::string & v) { return _data_storage->scalar(variable_index(v), _offset); }

template<>
inline std::complex<Real> & DataObject::data< std::complex<Real> >(const std::string & v) { return _data_storage->complex(variable_index(v), _offset); }

template<>
inline VectorValue<Real> & DataObject::data< VectorValue<Real> >(const std::string & v) { return _data_storage->vector(variable_index(v), _offset); }

template<>
inline TensorValue<Real> & DataObject::data< TensorValue<Real> >(const std::string & v) { return _data_storage->tensor(variable_index(v), _offset); }



template <typename T>
inline const T & DataObject::data(const std::string & ) const { }

template<>
inline const Real & DataObject::data< Real >(const std::string & v) const { return _data_storage->scalar(variable_index(v), _offset); }

template<>
inline const std::complex<Real> & DataObject::data< std::complex<Real> >(const std::string & v) const { return _data_storage->complex(variable_index(v), _offset); }

template<>
inline const VectorValue<Real> & DataObject::data< VectorValue<Real> >(const std::string & v) const { return _data_storage->vector(variable_index(v), _offset); }

template<>
inline const TensorValue<Real> & DataObject::data< TensorValue<Real> >(const std::string & v) const { return _data_storage->tensor(variable_index(v), _offset); }

#endif
