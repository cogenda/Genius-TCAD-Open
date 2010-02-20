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



#ifndef __fvm_cell_data_semiconductor_h__
#define __fvm_cell_data_semiconductor_h__

#include "petsc.h"
#include "fvm_cell_data.h"


/**
 *  FVM cell data for semiconductor region
 */
class FVM_Semiconductor_CellData : public FVM_CellData
{

public:

  /**
   * the auxiliary variable for semiconductor region
   */
  enum   SemiconductorAuxData
  {
    /**
     * electron mobility
     */
    _mun_,

    /**
     * hole mobility
     */
    _mup_
  };

  /**
   * the vector auxiliary variable for semiconductor region
   */
  enum SemiconductorAuxVecData
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
  FVM_Semiconductor_CellData()
  {
    _scalar_value = new PetscScalar[n_scalar()];
    for(unsigned int i=0; i<n_scalar(); i++) _scalar_value[i]=0.0;

    _vecctor_value = new VectorValue<PetscScalar>[n_vector()];
    for(unsigned int i=0; i<n_vector(); i++) _vecctor_value[i]=VectorValue<PetscScalar>(0.0, 0.0, 0.0);
  }

  /**
   * destructor
   */
  virtual ~FVM_Semiconductor_CellData()  {}

public:

  /**
   * @return the solution variable number
   */
  virtual size_t n_scalar() const
  {  return static_cast<unsigned int>(_mup_) +1 ; /* return last enum+1*/ }


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
  { return static_cast<unsigned int>(_Jp_) +1 ; }

  /**
   * @return the tensor variable number
   */
  virtual size_t n_tensor() const
  { return 0; }

  /**
   * @return the data type
   */
  virtual CellDataType type() const
    { return FVM_CellData::SemiconductorData; }


public:
  /**
   * @return average electron mobility
   */
  virtual PetscScalar         mun()          const
  { return _scalar_value[_mun_]; }

  /**
   * @return the writable reference to average electron mobility
   */
  virtual PetscScalar &       mun()
  { return _scalar_value[_mun_]; }

  /**
   * @return average hole mobility
   */
  virtual PetscScalar         mup()          const
  { return _scalar_value[_mup_]; }

  /**
   * @return the writable reference to average hole mobility
   */
  virtual PetscScalar &       mup()
  { return _scalar_value[_mup_]; }


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


  /**
   * @return the electron current
   */
  virtual VectorValue<PetscScalar> Jn()       const
  { return _vecctor_value[_Jn_];}


  /**
   * @return the writable reference to electron current
   */
  virtual VectorValue<PetscScalar> & Jn()
  { return _vecctor_value[_Jn_];}


  /**
   * @return the hole current
   */
  virtual VectorValue<PetscScalar> Jp()       const
  { return _vecctor_value[_Jp_];}


  /**
   * @return the writable reference to hole current
   */
  virtual VectorValue<PetscScalar> & Jp()
  { return _vecctor_value[_Jp_];}

};




#endif
