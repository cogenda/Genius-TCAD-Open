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
  FVM_NodeData(DataStorage * data_storage, const std::map<std::string, SimulationVariable> & variables)
  : DataObject(data_storage, variables)
  {}

  /**
   * destruction, no pointer here
   */
  virtual ~FVM_NodeData() {}


  enum NodeDataType {SemiconductorData, InsulatorData, ConductorData, ResistanceData, VacuumData, InvalidData };

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
  virtual PetscScalar  get_variable_real(SolutionVariable variable) const=0;

  /**
   * set variable by enum name
   */
  virtual void set_variable_real(SolutionVariable variable, PetscScalar value)=0;


  /**
   * @return the minimal distange to surface
   */
  virtual PetscScalar         dmin()          const
    { return 0; }

  /**
   * @return the writable reference to the minimal distange to surface
   */
  virtual PetscScalar &       dmin()
  { return FVM_NodeData::_scalar_dummy_; }

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
   * @return the correction of potential
   */
  virtual PetscScalar         dpsi()          const
  { return 0; }

  /**
   * @return the writable reference to correction of potential
   */
  virtual PetscScalar &       dpsi()
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
   * @return the H+ densiy
   */
  virtual PetscScalar         HIon()          const
    { return 0; }

  /**
   * @return the writable reference to H+ densiy
   */
  virtual PetscScalar &       HIon()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the A type trap density
   */
  virtual PetscScalar         trap_a()          const
    { return 0; }

  /**
   * @return the writable reference to A type trap density
   */
  virtual PetscScalar &       trap_a()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the B type trap density
   */
  virtual PetscScalar         trap_b()          const
    { return 0; }


  /**
   * @return the writable reference to B type trap density
   */
  virtual PetscScalar &       trap_b()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the neutral state B type trap density
   */
  virtual PetscScalar         trap_bn()          const
    { return 0; }


  /**
   * @return the writable reference to neutral state B type trap density
   */
  virtual PetscScalar &       trap_bn()
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
   * @return the old statistic potential
   */
  virtual PetscScalar         psi_old()          const
  { return 0; }

  /**
   * @return the writable reference to old statistic potential
   */
  virtual PetscScalar &       psi_old()
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
   * @return the old electron density
   */
  virtual PetscScalar         n_old()          const
  { return 0; }

  /**
   * @return the writable reference to old electron density
   */
  virtual PetscScalar &       n_old()
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
   * @return the old hole density
   */
  virtual PetscScalar         p_old()          const
  { return 0; }

  /**
   * @return the writable reference to old hole density
   */
  virtual PetscScalar &       p_old()
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
   * @return the change of conduction band due to strain
   */
  virtual PetscScalar         dEcStrain()          const
  { return 0; }

  /**
   * @return the writable reference to the change of conduction band due to strain
   */
  virtual PetscScalar &       dEcStrain()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the change of valance band due to strain
   */
  virtual PetscScalar         dEvStrain()          const
  { return 0; }

  /**
   * @return the writable reference to the change of valance band due to strain
   */
  virtual PetscScalar &       dEvStrain()
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
   * @return the optical energy
   */
  virtual PetscScalar         OptE()          const
  { return 0; }

  /**
   * @return the writable optical energy
   */
  virtual PetscScalar &       OptE()
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
   * @return the particle energy
   */
  virtual PetscScalar         PatE()          const
  { return 0; }

  /**
   * @return the writable particle generation ratio
   */
  virtual PetscScalar &       PatE()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the dose rate, Gy/s
   */
  virtual PetscScalar         DoseRate()          const
  { return 0; }

  /**
   * @return the writable reference to dose rate
   */
  virtual PetscScalar &       DoseRate()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the electron injected in to the FVM cell.
   * NOTE: it is the flux flow into the FVM cell or total electron density generated in the cell
   */
  virtual PetscScalar         EIn()          const
  { return 0; }

  /**
   * @return the writable reference of electron injected in to the FVM cell.
   */
  virtual PetscScalar &       EIn()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the hole inject injected in to the FVM cell.
   */
  virtual PetscScalar         HIn()          const
  { return 0; }

  /**
   * @return the writable reference of hole injected in to the FVM cell.
   */
  virtual PetscScalar &       HIn()
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
   * @return the quasi-fermi potential of electron
   */
  virtual PetscScalar  &      qFn()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the quasi-fermi potential of hole
   */
  virtual PetscScalar         qFp()           const
  { return 0; }

  /**
   * @return the quasi-fermi potential of hole
   */
  virtual PetscScalar &        qFp()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the charge density
   * this variable is used for data exchange
   * between hdm solver and poisson solver
   */
  virtual PetscScalar         rho()          const
  { return 0; }

  /**
   * @return the writable reference to charge density
   */
  virtual PetscScalar &       rho()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the interface charge density
   */
  virtual PetscScalar         interface_charge()          const
  { return 0; }

  /**
   * @return the writable reference to interface charge density
   */
  virtual PetscScalar &       interface_charge()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the recombnation rate
   */
  virtual PetscScalar         Recomb()          const
  { return 0; }

  /**
   * @return the writable reference to recombnation rate
   */
  virtual PetscScalar &       Recomb()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the direct(optical) recombnation rate
   */
  virtual PetscScalar         Recomb_Dir()          const
  { return 0; }

  /**
   * @return the writable reference to direct(optical) recombnation rate
   */
  virtual PetscScalar &       Recomb_Dir()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the SRH recombnation rate
   */
  virtual PetscScalar         Recomb_SRH()          const
  { return 0; }

  /**
   * @return the writable reference to SRH recombnation rate
   */
  virtual PetscScalar &       Recomb_SRH()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the impact ionization
   */
  virtual PetscScalar         ImpactIonization()          const
  { return 0; }

  /**
   * @return the writable reference to impact ionization
   */
  virtual PetscScalar &       ImpactIonization()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the Auger recombnation rate
   */
  virtual PetscScalar         Recomb_Auger()          const
  { return 0; }

  /**
   * @return the writable reference to Auger recombnation rate
   */
  virtual PetscScalar &       Recomb_Auger()
  { return FVM_NodeData::_scalar_dummy_; }


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


  /**
   * @return the stress tensor
   */
  virtual TensorValue<PetscScalar> stress()       const
  { return TensorValue<PetscScalar>();}


  /**
   * @return the writable reference to stress tensor
   */
  virtual TensorValue<PetscScalar> & stress()
  { return _tensor_dummy_;}


  /**
   * @return the strain tensor
   */
  virtual TensorValue<PetscScalar> strain()       const
  { return TensorValue<PetscScalar>();}


  /**
   * @return the writable reference to strain tensor
   */
  virtual TensorValue<PetscScalar> & strain()
  { return _tensor_dummy_;}

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

  /**
   * dummy tensor parameter to avoid compile problem
   */
  static TensorValue<PetscScalar> _tensor_dummy_;
};


#endif
