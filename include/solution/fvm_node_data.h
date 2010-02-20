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

//  $Id: fvm_node_data.h,v 1.19 2008/07/09 05:58:16 gdiso Exp $

#ifndef __fvm_node_data_h__
#define __fvm_node_data_h__

#include <vector>
#include <map>

#include "data_object.h"
#include "enum_solution.h"

/**
 * the class for store nodal data used in FVM solution
 */
class FVM_NodeData : public DataObject
{

public:
  /**
   * construction
   */
  FVM_NodeData() : DataObject()
  {}

  /**
   * destruction, no pointer here
   */
  virtual ~FVM_NodeData() {}


  enum NodeDataType {SemiconductorData, InsulatorData, ConductorData, VacuumData, InvalidData };

  /**
   * @return the data type
   */
  virtual NodeDataType type() const=0;

  /**
   * @return true when this variable valid
   */
  virtual bool is_variable_valid(SolutionVariable variable)  const=0;

  /**
   * @return data by enum name
   */
  virtual PetscScalar  get_variable(SolutionVariable variable) const=0;

  /**
   * set variable by enum name
   */
  virtual void set_variable(SolutionVariable variable, PetscScalar value)=0;

  //--------------------------------------------------------------------
  // these virtual functions should be overloaded by its derived classes
  //--------------------------------------------------------------------

  /**
   * @return the statistic potential
   */
  virtual PetscScalar         psi()          const
    { return 0; }

  /**
   * @return the writable reference to statistic potential
   */
  virtual PetscScalar &       psi()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the lattice temperature
   */
  virtual PetscScalar         T()          const
    { return 0; }

  /**
   * @return the writable reference to lattice temperature
   */
  virtual PetscScalar &       T()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the electron density
   */
  virtual PetscScalar         n()          const
    { return 0; }

  /**
   * @return the writable reference to electron density
   */
  virtual PetscScalar &       n()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the hole density
   */
  virtual PetscScalar         p()          const
    { return 0; }

  /**
   * @return the writable reference to hole density
   */
  virtual PetscScalar &       p()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the electron temperature
   */
  virtual PetscScalar         Tn()          const
  { return 0; }

  /**
   * @return the writable reference to electron temperature
   */
  virtual PetscScalar &       Tn()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the hole temperature
   */
  virtual PetscScalar         Tp()          const
  { return 0; }

  /**
   * @return the writable reference to hole temperature
   */
  virtual PetscScalar &       Tp()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the statistic potential
   */
  virtual std::complex<PetscScalar>         psi_ac()          const
    { return std::complex<PetscScalar>(0.0, 0.0); }

  /**
   * @return the writable reference to statistic potential
   */
  virtual std::complex<PetscScalar> &       psi_ac()
  { return FVM_NodeData::_complex_dummy_; }


  /**
   * @return the lattice temperature
   */
  virtual std::complex<PetscScalar>         T_ac()          const
    { return std::complex<PetscScalar>(0.0, 0.0); }

  /**
   * @return the writable reference to lattice temperature
   */
  virtual std::complex<PetscScalar> &       T_ac()
  { return FVM_NodeData::_complex_dummy_; }



  /**
   * @return the electron density
   */
  virtual std::complex<PetscScalar>         n_ac()          const
    { return std::complex<PetscScalar>(0.0, 0.0); }

  /**
   * @return the writable reference to electron density
   */
  virtual std::complex<PetscScalar> &       n_ac()
  { return FVM_NodeData::_complex_dummy_; }



  /**
   * @return the hole density
   */
  virtual std::complex<PetscScalar>         p_ac()          const
    { return std::complex<PetscScalar>(0.0, 0.0); }

  /**
   * @return the writable reference to hole density
   */
  virtual std::complex<PetscScalar> &       p_ac()
  { return FVM_NodeData::_complex_dummy_; }


  /**
   * @return the electron temperature
   */
  virtual std::complex<PetscScalar>         Tn_ac()          const
  { return std::complex<PetscScalar>(0.0, 0.0); }

  /**
   * @return the writable reference to electron temperature
   */
  virtual std::complex<PetscScalar> &       Tn_ac()
  { return FVM_NodeData::_complex_dummy_; }

  /**
   * @return the hole temperature
   */
  virtual std::complex<PetscScalar>         Tp_ac()          const
  { return std::complex<PetscScalar>(0.0, 0.0); }

  /**
   * @return the writable reference to hole temperature
   */
  virtual std::complex<PetscScalar> &       Tp_ac()
  { return FVM_NodeData::_complex_dummy_; }


  /**
   * @return the complex E file. only used by EM FEM solver
   */
  virtual std::complex<PetscScalar>         OptE_complex()          const
  { return std::complex<PetscScalar>(0.0, 0.0); }

  /**
   * @return the writable reference to complex E file. only used by EM FEM solver
   */
  virtual std::complex<PetscScalar> &       OptE_complex()
  { return FVM_NodeData::_complex_dummy_; }

  /**
   * @return the complex H file. only used by EM FEM solver
   */
  virtual std::complex<PetscScalar>         OptH_complex()          const
  { return std::complex<PetscScalar>(0.0, 0.0); }

  /**
   * @return the writable reference to complex H file. only used by EM FEM solver
   */
  virtual std::complex<PetscScalar> &       OptH_complex()
  { return FVM_NodeData::_complex_dummy_; }


  /**
   * @return the quantum conduction band
   */
  virtual PetscScalar         Eqc()          const
    { return 0; }

  /**
   * @return the writable reference to quantum conduction band
   */
  virtual PetscScalar &       Eqc()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the quantum valence band
   */
  virtual PetscScalar         Eqv()          const
    { return 0; }

  /**
   * @return the writable reference to quantum valence band
   */
  virtual PetscScalar &       Eqv()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the statistic potential at previous time step
   */
  virtual PetscScalar         psi_last()          const
    { return 0; }

  /**
   * @return the writable reference to statistic potential at previous time step
   */
  virtual PetscScalar &       psi_last()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the lattice temperature at previous time step
   */
  virtual PetscScalar         T_last()          const
    { return 0; }

  /**
   * @return the writable reference to lattice temperature at previous time step
   */
  virtual PetscScalar &       T_last()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the electron density at previous time step
   */
  virtual PetscScalar         n_last()          const
    { return 0; }

  /**
   * @return the writable reference to electron density at previous time step
   */
  virtual PetscScalar &       n_last()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the hole density at previous time step
   */
  virtual PetscScalar         p_last()          const
    { return 0; }

  /**
   * @return the writable reference to hole density at previous time step
   */
  virtual PetscScalar &       p_last()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the electron temperature at previous time step
   */
  virtual PetscScalar         Tn_last()          const
  { return 0; }

  /**
   * @return the writable reference to electron temperature at previous time step
   */
  virtual PetscScalar &       Tn_last()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the hole temperature at previous time step
   */
  virtual PetscScalar         Tp_last()          const
  { return 0; }

  /**
   * @return the writable reference to hole temperature at previous time step
   */
  virtual PetscScalar &       Tp_last()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the quantum conduction band at previous time step
   */
  virtual PetscScalar         Eqc_last()          const
    { return 0; }

  /**
   * @return the writable reference to quantum conduction band at previous time step
   */
  virtual PetscScalar &       Eqc_last()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the quantum valence band at previous time step
   */
  virtual PetscScalar         Eqv_last()          const
    { return 0; }

  /**
   * @return the writable reference to quantum valence band at previous time step
   */
  virtual PetscScalar &       Eqv_last()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the electron affinity
   */
  virtual PetscScalar         affinity()          const
    { return 0; }

  /**
   * @return the writable reference to the electron affinity
   */
  virtual PetscScalar &       affinity()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the conduction band
   */
  virtual PetscScalar         Ec()          const
  { return 0; }

  /**
   * @return the writable reference to the conduction band
   */
  virtual PetscScalar &       Ec()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the valance band
   */
  virtual PetscScalar         Ev()          const
  { return 0; }

  /**
   * @return the writable reference to the valance band
   */
  virtual PetscScalar &       Ev()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the mass density of the material
   */
  virtual PetscScalar         density()          const
    { return 0; }

  /**
   * @return the writable reference to the mass density of the material
   */
  virtual PetscScalar &       density()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the dielectric permittivity
   */
  virtual PetscScalar         eps()          const
    { return 0; }

  /**
   * @return the writable reference to the dielectric permittivity
   */
  virtual PetscScalar &       eps()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the megnetic permeability
   */
  virtual PetscScalar         mu()          const
    { return 0; }

  /**
   * @return the writable reference to the megnetic permeability
   */
  virtual PetscScalar &       mu()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the effective density of states in the conduction band
   */
  virtual PetscScalar         Nc()          const
    { return 0; }

  /**
   * @return the writable reference to the effective density of states in the conduction band
   */
  virtual PetscScalar &       Nc()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the effective density of states in the valence band
   */
  virtual PetscScalar         Nv()          const
    { return 0; }

  /**
   * @return the writable reference to the effective density of states in the valence band
   */
  virtual PetscScalar &       Nv()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the bandgap
   */
  virtual PetscScalar         Eg()          const
    { return 0; }

  /**
   * @return the writable reference to the bandgap
   */
  virtual PetscScalar &       Eg()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the mole fraction for single compound material
   */
  virtual PetscScalar         mole_x()          const
    { return 0; }

  /**
   * @return the writable reference to the mole fraction for single compound material
   */
  virtual PetscScalar &       mole_x()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
  * @return the mole fraction for dual compound material
  */
  virtual PetscScalar         mole_y()          const
    { return 0; }

  /**
   * @return the writable reference to the mole fraction for dual compound material
   */
  virtual PetscScalar &       mole_y()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the electron mobility
   */
  virtual PetscScalar         mun()          const
    { return 0; }

  /**
   * @return the writable reference to electron mobility
   */
  virtual PetscScalar &       mun()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the hole mobility
   */
  virtual PetscScalar         mup()          const
    { return 0; }

  /**
   * @return the writable reference to hole mobility
   */
  virtual PetscScalar &       mup()
  { return FVM_NodeData::_scalar_dummy_; }




  /**
   * @return the general doping concentration of acceptor
   */
  virtual PetscScalar         Na()          const
    { return 0; }

  /**
   * @return the writable reference to the general doping concentration of acceptor
   */
  virtual PetscScalar &       Na()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the general doping concentration of donor
   */
  virtual PetscScalar         Nd()          const
    { return 0; }

  /**
   * @return the writable reference to the general doping concentration of donor
   */
  virtual PetscScalar &       Nd()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the concentration of donor atom phosphorus
   */
  virtual PetscScalar         P()          const
  { return 0; }

  /**
   * @return the writable reference to concentration of donor atom phosphorus
   */
  virtual PetscScalar &       P()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the concentration of donor atom arsenic
   */
  virtual PetscScalar         As()          const
  { return 0; }

  /**
   * @return the writable reference to the concentration of donor atom arsenic
   */
  virtual PetscScalar &       As()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the carrier generation ratio due to OptG and PatG
   */
  virtual PetscScalar         Field_G()          const
  { return 0; }

  /**
   * @return the writable carrier generation ratio due to OptG and PatG
   */
  virtual PetscScalar &       Field_G()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the optical generation ratio
   */
  virtual PetscScalar         OptG()          const
  { return 0; }

  /**
   * @return the writable optical generation ratio
   */
  virtual PetscScalar &       OptG()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the heat generation ratio due to optical incident
   */
  virtual PetscScalar         OptQ()          const
  { return 0; }

  /**
   * @return the writable heat generation ratio due to optical incident
   */
  virtual PetscScalar &       OptQ()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the particle generation ratio
   */
  virtual PetscScalar         PatG()          const
  { return 0; }

  /**
   * @return the writable particle generation ratio
   */
  virtual PetscScalar &       PatG()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the concentration of donor atom antimony
   */
  virtual PetscScalar         Sb()          const
  { return 0; }

  /**
   * @return the writable reference to concentration of donor atom antimony
   */
  virtual PetscScalar &       Sb()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the concentration of acceptor atom boron
   */
  virtual PetscScalar         B()          const
  { return 0; }

  /**
   * @return the writable reference to the concentration of acceptor atom boron
   */
  virtual PetscScalar &       B()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the total acceptor concentration
   */
  virtual PetscScalar         Total_Na()     const
    { return 0; }

  /**
   * @return the total donor concentration
   */
  virtual PetscScalar         Total_Nd()     const
    {return 0;}

  /**
   * @return net concentration
   */
  virtual PetscScalar         Net_doping()   const
    {return 0;}

  /**
   * @return the total donor concentration
   */
  virtual PetscScalar         Total_doping() const
    {return 0;}

  /**
   * @return net charge concentration
   */
  virtual PetscScalar         Net_charge()   const
    {return 0;}

  /**
   * @return the intrinsic carrier concentration.
   * @note will not consider bandgap narrowing
   */
  virtual PetscScalar         ni()           const
    { return 0; }

  /**
   * @return the quasi-fermi potential of electron
   */
  virtual PetscScalar         qFn()           const
  { return 0; }

  /**
   * @return the quasi-fermi potential of hole
   */
  virtual PetscScalar         qFp()           const
  { return 0; }

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
   * dummy scalar parameter to avoid compile problem
   */
  static std::complex<PetscScalar> _complex_dummy_;

  /**
   * dummy vector parameter to avoid compile problem
   */
  static VectorValue<PetscScalar> _vector_dummy_;
};


#endif
