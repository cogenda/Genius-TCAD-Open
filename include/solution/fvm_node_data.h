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
  virtual Real  get_variable_real(SolutionVariable variable) const=0;

  /**
   * set variable by enum name
   */
  virtual void set_variable_real(SolutionVariable variable, Real value)=0;

  //--------------------------------------------------------------------
  // these virtual functions should be overloaded by its derived classes
  //--------------------------------------------------------------------

  /**
   * @return the statistic potential
   */
  virtual Real         psi()          const
    { return 0; }

  /**
   * @return the writable reference to statistic potential
   */
  virtual Real &       psi()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the correction of potential
   */
  virtual Real         dpsi()          const
  { return 0; }

  /**
   * @return the writable reference to correction of potential
   */
  virtual Real &       dpsi()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the lattice temperature
   */
  virtual Real         T()          const
    { return 0; }

  /**
   * @return the writable reference to lattice temperature
   */
  virtual Real &       T()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the electron density
   */
  virtual Real         n()          const
    { return 0; }

  /**
   * @return the writable reference to electron density
   */
  virtual Real &       n()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the hole density
   */
  virtual Real         p()          const
    { return 0; }

  /**
   * @return the writable reference to hole density
   */
  virtual Real &       p()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the electron temperature
   */
  virtual Real         Tn()          const
  { return 0; }

  /**
   * @return the writable reference to electron temperature
   */
  virtual Real &       Tn()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the hole temperature
   */
  virtual Real         Tp()          const
  { return 0; }

  /**
   * @return the writable reference to hole temperature
   */
  virtual Real &       Tp()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the statistic potential
   */
  virtual std::complex<Real>         psi_ac()          const
    { return std::complex<Real>(0.0, 0.0); }

  /**
   * @return the writable reference to statistic potential
   */
  virtual std::complex<Real> &       psi_ac()
  { return FVM_NodeData::_complex_dummy_; }


  /**
   * @return the lattice temperature
   */
  virtual std::complex<Real>         T_ac()          const
    { return std::complex<Real>(0.0, 0.0); }

  /**
   * @return the writable reference to lattice temperature
   */
  virtual std::complex<Real> &       T_ac()
  { return FVM_NodeData::_complex_dummy_; }



  /**
   * @return the electron density
   */
  virtual std::complex<Real>         n_ac()          const
    { return std::complex<Real>(0.0, 0.0); }

  /**
   * @return the writable reference to electron density
   */
  virtual std::complex<Real> &       n_ac()
  { return FVM_NodeData::_complex_dummy_; }



  /**
   * @return the hole density
   */
  virtual std::complex<Real>         p_ac()          const
    { return std::complex<Real>(0.0, 0.0); }

  /**
   * @return the writable reference to hole density
   */
  virtual std::complex<Real> &       p_ac()
  { return FVM_NodeData::_complex_dummy_; }


  /**
   * @return the electron temperature
   */
  virtual std::complex<Real>         Tn_ac()          const
  { return std::complex<Real>(0.0, 0.0); }

  /**
   * @return the writable reference to electron temperature
   */
  virtual std::complex<Real> &       Tn_ac()
  { return FVM_NodeData::_complex_dummy_; }

  /**
   * @return the hole temperature
   */
  virtual std::complex<Real>         Tp_ac()          const
  { return std::complex<Real>(0.0, 0.0); }

  /**
   * @return the writable reference to hole temperature
   */
  virtual std::complex<Real> &       Tp_ac()
  { return FVM_NodeData::_complex_dummy_; }


  /**
   * @return the complex E file. only used by EM FEM solver
   */
  virtual std::complex<Real>         OptE_complex()          const
  { return std::complex<Real>(0.0, 0.0); }

  /**
   * @return the writable reference to complex E file. only used by EM FEM solver
   */
  virtual std::complex<Real> &       OptE_complex()
  { return FVM_NodeData::_complex_dummy_; }

  /**
   * @return the complex H file. only used by EM FEM solver
   */
  virtual std::complex<Real>         OptH_complex()          const
  { return std::complex<Real>(0.0, 0.0); }

  /**
   * @return the writable reference to complex H file. only used by EM FEM solver
   */
  virtual std::complex<Real> &       OptH_complex()
  { return FVM_NodeData::_complex_dummy_; }


  /**
   * @return the statistic potential at previous time step
   */
  virtual Real         psi_last()          const
    { return 0; }

  /**
   * @return the writable reference to statistic potential at previous time step
   */
  virtual Real &       psi_last()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the old statistic potential
   */
  virtual Real         psi_old()          const
  { return 0; }

  /**
   * @return the writable reference to old statistic potential
   */
  virtual Real &       psi_old()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the lattice temperature at previous time step
   */
  virtual Real         T_last()          const
    { return 0; }

  /**
   * @return the writable reference to lattice temperature at previous time step
   */
  virtual Real &       T_last()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the electron density at previous time step
   */
  virtual Real         n_last()          const
    { return 0; }

  /**
   * @return the writable reference to electron density at previous time step
   */
  virtual Real &       n_last()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the old electron density
   */
  virtual Real         n_old()          const
  { return 0; }

  /**
   * @return the writable reference to old electron density
   */
  virtual Real &       n_old()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the hole density at previous time step
   */
  virtual Real         p_last()          const
    { return 0; }

  /**
   * @return the writable reference to hole density at previous time step
   */
  virtual Real &       p_last()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the old hole density
   */
  virtual Real         p_old()          const
  { return 0; }

  /**
   * @return the writable reference to old hole density
   */
  virtual Real &       p_old()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the electron temperature at previous time step
   */
  virtual Real         Tn_last()          const
  { return 0; }

  /**
   * @return the writable reference to electron temperature at previous time step
   */
  virtual Real &       Tn_last()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the hole temperature at previous time step
   */
  virtual Real         Tp_last()          const
  { return 0; }

  /**
   * @return the writable reference to hole temperature at previous time step
   */
  virtual Real &       Tp_last()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the electron affinity
   */
  virtual Real         affinity()          const
    { return 0; }

  /**
   * @return the writable reference to the electron affinity
   */
  virtual Real &       affinity()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the conduction band
   */
  virtual Real         Ec()          const
  { return 0; }

  /**
   * @return the writable reference to the conduction band
   */
  virtual Real &       Ec()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the valance band
   */
  virtual Real         Ev()          const
  { return 0; }

  /**
   * @return the writable reference to the valance band
   */
  virtual Real &       Ev()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the mass density of the material
   */
  virtual Real         density()          const
    { return 0; }

  /**
   * @return the writable reference to the mass density of the material
   */
  virtual Real &       density()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the dielectric permittivity
   */
  virtual Real         eps()          const
    { return 0; }

  /**
   * @return the writable reference to the dielectric permittivity
   */
  virtual Real &       eps()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the megnetic permeability
   */
  virtual Real         mu()          const
    { return 0; }

  /**
   * @return the writable reference to the megnetic permeability
   */
  virtual Real &       mu()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the effective density of states in the conduction band
   */
  virtual Real         Nc()          const
    { return 0; }

  /**
   * @return the writable reference to the effective density of states in the conduction band
   */
  virtual Real &       Nc()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the effective density of states in the valence band
   */
  virtual Real         Nv()          const
    { return 0; }

  /**
   * @return the writable reference to the effective density of states in the valence band
   */
  virtual Real &       Nv()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the bandgap
   */
  virtual Real         Eg()          const
    { return 0; }

  /**
   * @return the writable reference to the bandgap
   */
  virtual Real &       Eg()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
   * @return the mole fraction for single compound material
   */
  virtual Real         mole_x()          const
    { return 0; }

  /**
   * @return the writable reference to the mole fraction for single compound material
   */
  virtual Real &       mole_x()
  { return FVM_NodeData::_scalar_dummy_; }



  /**
  * @return the mole fraction for dual compound material
  */
  virtual Real         mole_y()          const
    { return 0; }

  /**
   * @return the writable reference to the mole fraction for dual compound material
   */
  virtual Real &       mole_y()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the electron mobility
   */
  virtual Real         mun()          const
    { return 0; }

  /**
   * @return the writable reference to electron mobility
   */
  virtual Real &       mun()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the hole mobility
   */
  virtual Real         mup()          const
    { return 0; }

  /**
   * @return the writable reference to hole mobility
   */
  virtual Real &       mup()
  { return FVM_NodeData::_scalar_dummy_; }




  /**
   * @return the general doping concentration of acceptor
   */
  virtual Real         Na()          const
    { return 0; }

  /**
   * @return the writable reference to the general doping concentration of acceptor
   */
  virtual Real &       Na()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the general doping concentration of donor
   */
  virtual Real         Nd()          const
    { return 0; }

  /**
   * @return the writable reference to the general doping concentration of donor
   */
  virtual Real &       Nd()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the carrier generation ratio due to OptG and PatG
   */
  virtual Real         Field_G()          const
  { return 0; }

  /**
   * @return the writable carrier generation ratio due to OptG and PatG
   */
  virtual Real &       Field_G()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the optical generation ratio
   */
  virtual Real         OptG()          const
  { return 0; }

  /**
   * @return the writable optical generation ratio
   */
  virtual Real &       OptG()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the heat generation ratio due to optical incident
   */
  virtual Real         OptQ()          const
  { return 0; }

  /**
   * @return the writable heat generation ratio due to optical incident
   */
  virtual Real &       OptQ()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the optical energy
   */
  virtual Real         OptE()          const
  { return 0; }

  /**
   * @return the writable optical energy
   */
  virtual Real &       OptE()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the particle generation ratio
   */
  virtual Real         PatG()          const
  { return 0; }

  /**
   * @return the writable particle generation ratio
   */
  virtual Real &       PatG()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the particle energy
   */
  virtual Real         PatE()          const
  { return 0; }

  /**
   * @return the writable particle generation ratio
   */
  virtual Real &       PatE()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the electron injected in to the FVM cell.
   * NOTE: it is the flux flow into the FVM cell or total electron density generated in the cell
   */
  virtual Real         EIn()          const
  { return 0; }

  /**
   * @return the writable reference of electron injected in to the FVM cell.
   */
  virtual Real &       EIn()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the hole inject injected in to the FVM cell.
   */
  virtual Real         HIn()          const
  { return 0; }

  /**
   * @return the writable reference of hole injected in to the FVM cell.
   */
  virtual Real &       HIn()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the total acceptor concentration
   */
  virtual Real         Total_Na()     const
    { return 0; }

  /**
   * @return the total donor concentration
   */
  virtual Real         Total_Nd()     const
    {return 0;}

  /**
   * @return net concentration
   */
  virtual Real         Net_doping()   const
    {return 0;}

  /**
   * @return the total donor concentration
   */
  virtual Real         Total_doping() const
    {return 0;}

  /**
   * @return net charge concentration
   */
  virtual Real         Net_charge()   const
    {return 0;}

  /**
   * @return the intrinsic carrier concentration.
   * @note will not consider bandgap narrowing
   */
  virtual Real         ni()           const
    { return 0; }

  /**
   * @return the quasi-fermi potential of electron
   */
  virtual Real         qFn()           const
  { return 0; }

  /**
   * @return the quasi-fermi potential of electron
   */
  virtual Real  &      qFn()
  { return FVM_NodeData::_scalar_dummy_; }

  /**
   * @return the quasi-fermi potential of hole
   */
  virtual Real         qFp()           const
  { return 0; }

  /**
   * @return the quasi-fermi potential of hole
   */
  virtual Real &        qFp()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the charge density
   * this variable is used for data exchange
   * between hdm solver and poisson solver
   */
  virtual Real         rho()          const
  { return 0; }

  /**
   * @return the writable reference to charge density
   */
  virtual Real &       rho()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the recombnation rate
   */
  virtual Real         Recomb()          const
  { return 0; }

  /**
   * @return the writable reference to recombnation rate
   */
  virtual Real &       Recomb()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the direct(optical) recombnation rate
   */
  virtual Real         Recomb_Dir()          const
  { return 0; }

  /**
   * @return the writable reference to direct(optical) recombnation rate
   */
  virtual Real &       Recomb_Dir()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the SRH recombnation rate
   */
  virtual Real         Recomb_SRH()          const
  { return 0; }

  /**
   * @return the writable reference to SRH recombnation rate
   */
  virtual Real &       Recomb_SRH()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the impact ionization
   */
  virtual Real         ImpactIonization()          const
  { return 0; }

  /**
   * @return the writable reference to impact ionization
   */
  virtual Real &       ImpactIonization()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the Auger recombnation rate
   */
  virtual Real         Recomb_Auger()          const
  { return 0; }

  /**
   * @return the writable reference to Auger recombnation rate
   */
  virtual Real &       Recomb_Auger()
  { return FVM_NodeData::_scalar_dummy_; }


  /**
   * @return the electrical field
   */
  virtual VectorValue<Real> E()       const
  { return VectorValue<Real>(0,0,0);}


  /**
   * @return the writable reference to electrical field
   */
  virtual VectorValue<Real> & E()
  { return _vector_dummy_;}

  /**
   * @return the electron current
   */
  virtual VectorValue<Real> Jn()       const
  { return VectorValue<Real>(0,0,0);}


  /**
   * @return the writable reference to electron current
   */
  virtual VectorValue<Real> & Jn()
  { return _vector_dummy_;}


  /**
   * @return the hole current
   */
  virtual VectorValue<Real> Jp()       const
  { return VectorValue<Real>(0,0,0);}


  /**
   * @return the writable reference to hole current
   */
  virtual VectorValue<Real> & Jp()
  { return _vector_dummy_;}


protected:

  /**
   * dummy scalar parameter to avoid compile problem
   */
  static Real _scalar_dummy_;

  /**
   * dummy scalar parameter to avoid compile problem
   */
  static std::complex<Real> _complex_dummy_;

  /**
   * dummy vector parameter to avoid compile problem
   */
  static VectorValue<Real> _vector_dummy_;
};


#endif
