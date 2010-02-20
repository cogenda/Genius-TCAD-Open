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

#include "vector_value.h"
#include "tensor_value.h"

/**
 * The \p DataObject defines an abstract base class for objects that
 * have data associated with them.
 */
class DataObject
{

protected:

  /**
   * Constructor. Protected so that you can't instantiate one of these.
   */
  DataObject()
  : _scalar_value(0),
    _aux_scalar_value(0),
    _complex_value(0),
    _vecctor_value(0),
    _tensor_value(0)
  {}

public:

  /**
   * Copy-constructor.
   */
  //DataObject (const DataObject&);


  /**
   * destruction
   */
  virtual ~DataObject()
  {
    delete [] _scalar_value;
    delete [] _aux_scalar_value;
    delete [] _complex_value;
    delete [] _vecctor_value;
    delete [] _tensor_value;

    _scalar_value    = 0;
    _aux_scalar_value  = 0;
    _complex_value = 0;
    _vecctor_value = 0;
    _tensor_value  = 0;
  }


  /**
   * @return the solution variable number
   */
  virtual size_t n_scalar() const=0;

  /**
   * @return the scalar aux variable number
   */
  virtual size_t n_aux_scalar() const=0;

  /**
   * @return the complex variable number
   */
  virtual size_t n_complex() const=0;

  /**
   * @return the vector variable number
   */
  virtual size_t n_vector() const=0;

  /**
   * @return the tensor variable number
   */
  virtual size_t n_tensor() const=0;

  /**
   * create a user-defined value at the node
   */
  void CreateUserScalarValue(const std::string &name, PetscScalar var=0.0)
  {
    if(_user_defined_scalar_value_map.find(name) != _user_defined_scalar_value_map.end()) return;
    _user_defined_scalar_value.push_back(var);
    unsigned int n_ud_var = _user_defined_scalar_value.size();
    _user_defined_scalar_value_map.insert(std::pair<std::string, PetscScalar *>(name, &_user_defined_scalar_value[n_ud_var-1]));
  }

  /**
   * @return the value of a user-defined value at the node
   */
  PetscScalar UserScalarValue(const std::string &name) const
  {
    std::map<std::string, PetscScalar *>::const_iterator it;

    it = _user_defined_scalar_value_map.find(name);
    genius_assert(it != _user_defined_scalar_value_map.end());

    return *(it->second);
  }

  /**
   * @return the reference to a user-defined value at the node
   */
  PetscScalar & UserScalarValue(const std::string &name)
  {
    std::map<std::string,PetscScalar *>::iterator it;

    it = _user_defined_scalar_value_map.find(name);
    genius_assert(it != _user_defined_scalar_value_map.end());

    return *(it->second);
  }

  /**
   * @return the name list of user defined values
   */
  void GetUserScalarList(std::vector<std::string> & name) const
  {
    std::map<std::string, PetscScalar *>::const_iterator it;
    for(it=_user_defined_scalar_value_map.begin(); it!=_user_defined_scalar_value_map.end(); it++)
    {
      name.push_back(it->first);
    }
  }

protected:

  // change vector<PetscScalar> xxx to PetscScalar  * xxx
  // hoping to save memory and speed up data access.

  /**
   * the field value for independent variable
   */
  PetscScalar  * _scalar_value;

  /**
   * the scalar auxiliary value used in solution
   */
  PetscScalar * _aux_scalar_value;

  /**
   * the complex value used in solution
   */
  std::complex<PetscScalar> * _complex_value;

  /**
   * the vector value used in solution
   */
  VectorValue<PetscScalar>  * _vecctor_value;

  /**
   * the tensor value used in solution
   */
  TensorValue<PetscScalar>  * _tensor_value;


  /**
   * user-defined values (PMIs stores data here)
   */
  std::vector<PetscScalar> _user_defined_scalar_value;

  /**
   * a map link value name to value
   */
  std::map<std::string, PetscScalar *> _user_defined_scalar_value_map;

};


#endif
