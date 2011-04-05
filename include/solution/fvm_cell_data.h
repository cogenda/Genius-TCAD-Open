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



#ifndef __fvm_cell_data_h__
#define __fvm_cell_data_h__

#include <vector>
#include <map>

#include "data_object.h"
#include "enum_solution.h"

/**
 * the class for store cell data used in FVM solution
 */
class FVM_CellData : public DataObject
{

public:
  /**
   * construction
   */
  FVM_CellData(DataStorage * data_storage, const std::map<std::string, SimulationVariable> & variables)
  : DataObject(data_storage, variables)
  {}

  /**
   * destruction, no pointer here
   */
  virtual ~FVM_CellData() {}


  enum CellDataType {SemiconductorData, InsulatorData, ConductorData, ResistanceData, VacuumData, InvalidData };

  /**
   * @return the data type
   */
  virtual CellDataType type() const=0;

  /**
   * @return the electrical field
   */
  virtual VectorValue<PetscScalar> E()       const
  { return VectorValue<PetscScalar>(0,0,0);}


  /**
   * @return the writable reference to electrical field
   */
  virtual VectorValue<PetscScalar> & E()
  { return _vector_dummy_;}


  /**
   * @return the electron current
   */
  virtual VectorValue<PetscScalar> Jn()       const
  { return VectorValue<PetscScalar>(0,0,0);}


  /**
   * @return the writable reference to electron current
   */
  virtual VectorValue<PetscScalar> & Jn()
  { return _vector_dummy_;}


  /**
   * @return the hole current
   */
  virtual VectorValue<PetscScalar> Jp()       const
  { return VectorValue<PetscScalar>(0,0,0);}


  /**
   * @return the writable reference to hole current
   */
  virtual VectorValue<PetscScalar> & Jp()
  { return _vector_dummy_;}


protected:

  /**
   * dummy scalar parameter to avoid compile problem
   */
  static PetscScalar _scalar_dummy_;

  /**
   * dummy vector parameter to avoid compile problem
   */
  static VectorValue<PetscScalar> _vector_dummy_;
};


#endif
