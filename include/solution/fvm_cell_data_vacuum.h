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



#ifndef __fvm_cell_data_vacuum_h__
#define __fvm_cell_data_vacuum_h__

#include "petsc.h"
#include "fvm_cell_data.h"




/**
 *  FVM cell data for vacuum region
 */
class FVM_Vacuum_CellData : public FVM_CellData
{

public:
  /**
   * the vector auxiliary variable for vacuum region
   */
  enum VacuumAuxVecData
  {
    /**
     * electrical field
     */
    _E_=0,

    /**
     * electron current field
     */
    _Jn_,

    /**
     * electron current field
     */
    _Jp_
  };



public:
  /**
   * constructor
   */
  FVM_Vacuum_CellData()
  {
    _vecctor_value = new VectorValue<PetscScalar>[n_vector()];
    for(unsigned int i=0; i<n_vector(); i++) _vecctor_value[i]=VectorValue<PetscScalar>(0.0, 0.0, 0.0);
  }

  /**
   * destructor
   */
  virtual ~FVM_Vacuum_CellData()  {}

public:

  /**
   * @return the solution variable number
   */
  virtual size_t n_scalar() const
  { return 0; }


  /**
   * @return the scalar aux variable number
   */
  virtual size_t n_aux_scalar() const
  { return 0; }

  /**
   * @return the complex variable number
   */
  virtual size_t n_complex() const
  { return 0; }

  /**
   * @return the vector variable number
   */
  virtual size_t n_vector() const
  { return static_cast<unsigned int>(_E_) +1 ; }

  /**
   * @return the tensor variable number
   */
  virtual size_t n_tensor() const
  { return 0; }

  /**
   * @return the data type
   */
  virtual CellDataType type() const
    { return FVM_CellData::ConductorData; }


public:

  /**
   * @return the electrical field
   */
  virtual VectorValue<PetscScalar> E()       const
  { return _vecctor_value[_E_];}


  /**
   * @return the writable reference to electrical field
   */
  virtual VectorValue<PetscScalar> & E()
  { return _vecctor_value[_E_];}

};




#endif
