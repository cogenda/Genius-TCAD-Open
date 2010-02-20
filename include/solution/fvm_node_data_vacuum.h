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



#ifndef __fvm_node_data_vacuum_h__
#define __fvm_node_data_vacuum_h__

#include "petsc.h"
#include "fvm_node_data.h"




/**
 *  FVM nodal data for vacuum region
 */
class FVM_Vacuum_NodeData : public FVM_NodeData
{

public:

  /**
   * the independent variable for vacuum region
   */
  enum   VacuumData
  {
    /**
     * electrostatic potential
     */
    _psi_
  };


  /**
   * the auxiliary variable for vacuum region
   */
  enum   VacuumAuxData
  {
    /**
     * the density of the material
     */
    _density_=0,

    /**
     * electron affinity
     */
    _affinity_,

    /**
     * the dielectric permittivity
     */
    _eps_,

    /**
     * the megnetic permeability
     */
    _mu_,

    /**
     * electrostatic potential at previous time step
     */
    _psi_last_
  };


  /**
   * the vector auxiliary variable for vacuum region
   */
  enum VacuumAuxVecData
  {
    /**
     * electrical field of incident optical wave
     */
    //_OpE_,

    /**
     * magnetic field of incident optical wave
     */
    //_OpH_,

    /**
     * electrical field
     */
    _E_
  };


  /**
   * the complex auxiliary variable for Vacuum region
   */
  enum VacuumAuxComplexData
  {
    /**
     * electrical field of incident optical wave
     */
    _OpE_complex_=0,

    /**
     * magnetic field of incident optical wave
     */
    _OpH_complex_
  };

public:
  /**
   * constructor
   */
  FVM_Vacuum_NodeData()
  {
    _scalar_value = new PetscScalar[n_scalar()];
    for(unsigned int i=0; i<n_scalar(); i++) _scalar_value[i]=0.0;

    _aux_scalar_value = new PetscScalar[n_aux_scalar()];
    for(unsigned int i=0; i<n_aux_scalar(); i++) _aux_scalar_value[i]=0.0;

    _complex_value = new std::complex<PetscScalar>[n_complex()];
    for(unsigned int i=0; i<n_scalar(); i++) _complex_value[i]=std::complex<PetscScalar>(0.0, 0.0);

    _vecctor_value = new VectorValue<PetscScalar>[n_vector()];
    for(unsigned int i=0; i<n_vector(); i++) _vecctor_value[i]=VectorValue<PetscScalar>(0.0, 0.0, 0.0);
  }

  /**
   * destructor
   */
  virtual ~FVM_Vacuum_NodeData()  {}

public:

  /**
   * @return the solution variable number
   */
  virtual size_t n_scalar() const
  { return static_cast<unsigned int>(_psi_) +1 ; /* return last enum+1*/ }


  /**
   * @return the scalar aux variable number
   */
  virtual size_t n_aux_scalar() const
  { return static_cast<unsigned int>(_psi_last_) +1 ; /* return last enum+1*/ }

  /**
   * @return the complex variable number
   */
  virtual size_t n_complex() const
  { return static_cast<unsigned int>(_OpH_complex_) +1 ; /* return last enum+1*/ }

  /**
   * @return the vector variable number
   */
  virtual size_t n_vector() const
  { return static_cast<unsigned int>(_E_) +1 ; /* return last enum+1*/ }

  /**
   * @return the tensor variable number
   */
  virtual size_t n_tensor() const
  { return 0; }

  /**
   * @return the data type
   */
  virtual NodeDataType type() const
    { return FVM_NodeData::VacuumData; }


public:

  /**
   * @return data by enum name
   */
  virtual PetscScalar  get_variable(SolutionVariable variable) const
  {
    switch(variable)
    {
    case POTENTIAL   :  return  psi();                            /* potential */
    case E_FIELD     :  return  _vecctor_value[_E_].size();       /* electric field */
    case QFN         :  return  psi();                            /* electron quasi-Fermi level */
    case QFP         :  return  psi();                            /* hole quasi-Fermi level */
    default          :  return  0.0;
    }
  }

  /**
   * set variable by enum name
   */
  virtual void set_variable(SolutionVariable variable, PetscScalar value)
  {
    switch(variable)
    {
    case POTENTIAL   :  psi() = value;                             /* potential */
    default          :  return;
    }
  }

  /**
   * @return true when this variable valid
   */
  virtual bool is_variable_valid(SolutionVariable variable)  const
  {
    switch(variable)
    {
    case POTENTIAL   :  return true;
    default          :  return false;
    }
  }

  //--------------------------------------------------------------------
  // data access function
  //--------------------------------------------------------------------

  /**
   * @return the statistic potential
   */
  virtual PetscScalar         psi()        const
    { return _scalar_value[_psi_]; }

  /**
   * @return the statistic potential
   */
  virtual PetscScalar &       psi()
  { return _scalar_value[_psi_]; }


  /**
   * @return the statistic potential
   */
  virtual std::complex<PetscScalar>         psi_ac()          const
    { return _complex_value[_psi_]; }

  /**
   * @return the writable reference to statistic potential
   */
  virtual std::complex<PetscScalar> &       psi_ac()
  { return _complex_value[_psi_]; }

  /**
   * @return the statistic potential at previous time step
   */
  virtual PetscScalar         psi_last()          const
    { return _aux_scalar_value[_psi_last_]; }

  /**
   * @return the writable reference to statistic potential at previous time step
   */
  virtual PetscScalar &       psi_last()
  { return _aux_scalar_value[_psi_last_]; }


  /**
   * @return the complex E file. only used by EM FEM solver
   */
  virtual std::complex<PetscScalar>         OptE_complex()          const
  { return _complex_value[_OpE_complex_]; }

  /**
   * @return the writable reference to complex E file. only used by EM FEM solver
   */
  virtual std::complex<PetscScalar> &       OptE_complex()
  { return _complex_value[_OpE_complex_]; }

  /**
   * @return the complex H file. only used by EM FEM solver
   */
  virtual std::complex<PetscScalar>         OptH_complex()          const
  { return _complex_value[_OpH_complex_]; }

  /**
   * @return the writable reference to complex H file. only used by EM FEM solver
   */
  virtual std::complex<PetscScalar> &       OptH_complex()
  { return _complex_value[_OpH_complex_]; }


  /**
   * @return the electron affinity
   */
  virtual PetscScalar         affinity()          const
    { return _aux_scalar_value[_affinity_]; }

  /**
   * @return the writable reference to the electron affinity
   */
  virtual PetscScalar &       affinity()
  { return _aux_scalar_value[_affinity_]; }



  /**
   * @return the mass density of the material
   */
  virtual PetscScalar         density()          const
    { return _aux_scalar_value[_density_]; }

  /**
   * @return the writable reference to the mass density of the material
   */
  virtual PetscScalar &       density()
  { return _aux_scalar_value[_density_]; }



  /**
   * @return the dielectric permittivity
   */
  virtual PetscScalar         eps()          const
    { return _aux_scalar_value[_eps_]; }

  /**
   * @return the writable reference to the dielectric permittivity
   */
  virtual PetscScalar &       eps()
  { return _aux_scalar_value[_eps_]; }



  /**
   * @return the megnetic permeability
   */
  virtual PetscScalar         mu()          const
    { return _aux_scalar_value[_mu_]; }

  /**
   * @return the writable reference to the megnetic permeability
   */
  virtual PetscScalar &       mu()
  { return _aux_scalar_value[_mu_]; }


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
