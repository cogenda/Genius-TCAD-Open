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



#ifndef __fvm_node_data_semiconductor_h__
#define __fvm_node_data_semiconductor_h__

#include "petsc.h"
#include "fvm_node_data.h"


/**
 *  FVM nodal data for semiconductor region
 */
class FVM_Semiconductor_NodeData : public FVM_NodeData
{

public:

  /**
   * the independent variable for semiconductor region
   */
  enum   SemiconductorData
  {
    /**
     * electron density
     */
    _n_=0,

    /**
     * hole density
     */
    _p_,

    /**
     * electrostatic potential
     */
    _psi_,

    /**
     * lattice temperature
     */
    _T_,

    /**
     * electron temperature, used in energy balance model simulation
     */
    _Tn_,

    /**
     * hole temperature, used in energy balance model simulation
     */
    _Tp_,

    /**
     * quantum conduction band
     */
    _Eqc_,

    /**
     * quantum valence band
     */
    _Eqv_
  };


  /**
   * the auxiliary variable for semiconductor region
   */
  enum   SemiconductorAuxData
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
     * conduction band
     */
    _Ec_,

    /**
     * valence band
     */
    _Ev_,

    /**
     * band gap
     */
    _Eg_,

    /**
     * quasi Fermi potential of electron
     */
    _qFn_,

    /**
     * quasi Fermi potential of hole
     */
    _qFp_,

    /**
     * intrinsic Fermi potential
     */
    //_phi_intrinsic_,

    /**
     * effective density of states in the conduction band
     */
    _Nc_,

    /**
     * effective density of states in the valence band
     */
    _Nv_,

    /**
     * the dielectric permittivity
     */
    _eps_,

    /**
     * the megnetic permeability
     */
    _mu_,

    /**
     * general doping concentration of acceptor
     */
    _Na_,

    /**
     * general doping concentration of donor
     */
    _Nd_,

    /**
     * concentration of donor atom phosphorus
     */
    _P_,

    /**
     * concentration of donor atom arsenic
     */
    _As_,

    /**
     * concentration of donor atom antimony
     */
    _Sb_,

    /**
     * concentration of acceptor atom boron
     */
    _B_,

    /**
     * mole fraction for single compound material
     */
    _mole_x_,

    /**
     * mole fraction for dual compound material
     */
    _mole_y_,

    /**
     * electron mobility
     */
    //_mun_,

    /**
     * hole mobility
     */
    //_mup_,

    /**
     * the _OptG_*time +  _PatG_*time
     */
    _Field_G_,

    /**
     * carrier generation due to incident wave
     */
    _OptG_,

    /**
     * heat generation due to incident wave
     */
    _OptQ_,

    /**
     * carrier generation due to high energy particle
     */
    _PatG_,

    /**
     * electron density at previous time step
     */
    _n_last_,

    /**
     * hole density at previous time step
     */
    _p_last_,

    /**
     * electrostatic potential at previous time step
     */
    _psi_last_,

    /**
     * lattice temperature at previous time step
     */
    _T_last_,

    /**
     * electron temperature, at previous time step
     */
    _Tn_last_,

    /**
     * hole temperature, at previous time step
     */
    _Tp_last_,

    /**
     * quantum conduction band, at previous time step
     */
    _Eqc_last_,

    /**
     * quantum valence band, at previous time step
     */
    _Eqv_last_
  };


  /**
   * the vector auxiliary variable for semiconductor region
   */
  enum SemiconductorAuxVecData
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
    _E_,

    /**
     * electron current
     */
    _Jn_,

    /**
     * hole current
     */
    _Jp_
  };


  /**
   * the complex auxiliary variable for semiconductor region
   */
  enum SemiconductorAuxComplexData
  {
    /**
     * electron density
     */
    _n_ac_=0,

    /**
     * hole density
     */
    _p_ac_,

    /**
     * electrostatic potential
     */
    _psi_ac_,

    /**
     * lattice temperature
     */
    _T_ac_,

    /**
     * electron temperature, used in energy balance model simulation
     */
    _Tn_ac_,

    /**
     * hole temperature, used in energy balance model simulation
     */
    _Tp_ac_,

    /**
     * electrical field of incident optical wave
     */
    _OpE_complex_,

    /**
     * magnetic field of incident optical wave
     */
    _OpH_complex_
  };




public:

  /**
   * constructor
   */
  FVM_Semiconductor_NodeData()
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
virtual ~FVM_Semiconductor_NodeData()  { }

public:
  /**
   * @return the solution variable number
   */
  virtual size_t n_scalar() const
    { return static_cast<unsigned int>(_Eqv_) +1 ; /* return last enum+1*/ }

  /**
   * @return the scalar aux variable number
   */
  virtual size_t n_aux_scalar() const
    { return static_cast<unsigned int>(_Eqv_last_) +1 ; /* return last enum+1*/ }

  /**
   * @return the complex variable number
   */
  virtual size_t n_complex() const
    { return static_cast<unsigned int>(_OpH_complex_) +1 ; /* return last enum+1*/ }

  /**
   * @return the vector variable number
   */
  virtual size_t n_vector() const
    { return static_cast<unsigned int>(_Jp_) +1 ; /* return last enum+1*/ }

  /**
   * @return the tensor variable number
   */
  virtual size_t n_tensor() const
    { return 0; }

  /**
   * @return the data type
   */
  virtual NodeDataType type() const
    { return FVM_NodeData::SemiconductorData; }


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
    case ELECTRON    :  return  n();                              /* electron concentration */
    case HOLE        :  return  p();                              /* hole concentration */
    case TEMPERATURE :  return  T();                              /* lattice temperature */
    case E_TEMP      :  return  Tn();                             /* electron temperature */
    case H_TEMP      :  return  Tp();                             /* hole temperature */
    case DOPING      :  return  Net_doping();                     /* net doping */
    case DOPING_Na   :  return  Total_Na();                       /* acceptor */
    case DOPING_Nd   :  return  Total_Nd();                       /* donor */
    case MOLE_X      :  return  mole_x();
    case MOLE_Y      :  return  mole_y();
    case MIN_CARRIER :  return  Net_doping() > 0 ? n() : p();     /* minority carrier concentration */
    case NET_CARRIER :  return  p() - n();                        /* net carrier concentration */
    case NET_CHARGE  :  return  Net_doping() + p() - n();         /* net charge */
    case OPTICAL_GEN :  return  OptG();                           /* charge genetated by optical ray */
    case OPTICAL_HEAT:  return  OptQ();                           /* heat genetated by optical ray */
    case PARTICLE_GEN:  return  PatG();                           /* charge genetated by particle ray */
    case QFN         :  return  qFn();                            /* electron quasi-Fermi level */
    case QFP         :  return  qFp();                            /* hole quasi-Fermi level */
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
    case ELECTRON    :  n() = value;                               /* electron concentration */
    case HOLE        :  p() = value;                               /* hole concentration */
    case TEMPERATURE :  T() = value;                               /* lattice temperature */
    case E_TEMP      :  Tn() = value;                              /* electron temperature */
    case H_TEMP      :  Tp() = value;                              /* hole temperature */
    case DOPING_Na   :  Na() = value;                              /* acceptor */
    case DOPING_Nd   :  Nd() = value;                              /* donor */
    case OPTICAL_GEN :  OptG() = value;                            /* charge genetated by optical ray */
    case OPTICAL_HEAT:  OptQ() = value;                            /* heat genetated by optical ray */
    case PARTICLE_GEN:  PatG() = value;                            /* charge genetated by particle ray */
    case MOLE_X      :  mole_x() = value;
    case MOLE_Y      :  mole_y() = value;
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
    case POTENTIAL   :
    case ELECTRON    :
    case HOLE        :
    case TEMPERATURE :
    case E_TEMP      :
    case H_TEMP      :
    case DOPING_Na   :
    case DOPING_Nd   :
    case OPTICAL_GEN :
    case OPTICAL_HEAT:
    case PARTICLE_GEN:
    case MOLE_X      :
    case MOLE_Y      :  return true;
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
   * @return the lattice temperature
   */
  virtual PetscScalar         T()          const
    { return _scalar_value[_T_]; }

  /**
   * @return the statistic potential
   */
  virtual PetscScalar &       T()
  { return _scalar_value[_T_]; }



  /**
   * @return the electron density
   */
  virtual PetscScalar         n()          const
    { return _scalar_value[_n_]; }

  /**
   * @return the writable reference to electron density
   */
  virtual PetscScalar &       n()
  { return _scalar_value[_n_]; }



  /**
   * @return the hole density
   */
  virtual PetscScalar         p()          const
    { return _scalar_value[_p_]; }


  /**
   * @return the writable reference to hole density
   */
  virtual PetscScalar &       p()
  { return _scalar_value[_p_]; }


  /**
   * @return the electron temperature
   */
  virtual PetscScalar         Tn()          const
    { return _scalar_value[_Tn_]; }


  /**
   * @return the writable reference to electron temperature
   */
  virtual PetscScalar &       Tn()
  { return _scalar_value[_Tn_]; }


  /**
   * @return the hole temperature
   */
  virtual PetscScalar         Tp()          const
    { return _scalar_value[_Tp_]; }


  /**
   * @return the writable reference to hole temperature
   */
  virtual PetscScalar &       Tp()
  { return _scalar_value[_Tp_]; }


  /**
   * @return the statistic potential
   */
  virtual std::complex<PetscScalar>         psi_ac()          const
    { return _complex_value[_psi_ac_]; }

  /**
   * @return the writable reference to statistic potential
   */
  virtual std::complex<PetscScalar> &       psi_ac()
  { return _complex_value[_psi_ac_]; }


  /**
   * @return the lattice temperature
   */
  virtual std::complex<PetscScalar>         T_ac()          const
    { return  _complex_value[_T_ac_]; }

  /**
   * @return the writable reference to lattice temperature
   */
  virtual std::complex<PetscScalar> &       T_ac()
  { return _complex_value[_T_ac_]; }



  /**
   * @return the electron density
   */
  virtual std::complex<PetscScalar>         n_ac()          const
    { return _complex_value[_n_ac_]; }

  /**
   * @return the writable reference to electron density
   */
  virtual std::complex<PetscScalar> &       n_ac()
  { return _complex_value[_n_ac_]; }



  /**
   * @return the hole density
   */
  virtual std::complex<PetscScalar>         p_ac()          const
    { return _complex_value[_p_ac_]; }

  /**
   * @return the writable reference to hole density
   */
  virtual std::complex<PetscScalar> &       p_ac()
  { return _complex_value[_p_ac_]; }


  /**
   * @return the electron temperature
   */
  virtual std::complex<PetscScalar>         Tn_ac()          const
    { return _complex_value[_Tn_ac_]; }

  /**
   * @return the writable reference to electron temperature
   */
  virtual std::complex<PetscScalar> &       Tn_ac()
  { return _complex_value[_Tn_ac_]; }

  /**
   * @return the hole temperature
   */
  virtual std::complex<PetscScalar>         Tp_ac()          const
    { return _complex_value[_Tp_ac_]; }

  /**
   * @return the writable reference to hole temperature
   */
  virtual std::complex<PetscScalar> &       Tp_ac()
  { return _complex_value[_Tp_ac_]; }


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
   * @return the quantum conduction band
   */
  virtual PetscScalar         Eqc()          const
    { return _scalar_value[_Eqc_]; }

  /**
   * @return the writable reference to quantum conduction band
   */
  virtual PetscScalar &       Eqc()
  { return _scalar_value[_Eqc_]; }



  /**
   * @return the quantum valence band
   */
  virtual PetscScalar         Eqv()          const
    { return _scalar_value[_Eqv_]; }

  /**
   * @return the writable reference to quantum valence band
   */
  virtual PetscScalar &       Eqv()
  { return _scalar_value[_Eqv_]; }



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
   * @return the lattice temperature at previous time step
   */
  virtual PetscScalar         T_last()          const
    { return _aux_scalar_value[_T_last_]; }

  /**
   * @return the writable reference to lattice temperature at previous time step
   */
  virtual PetscScalar &       T_last()
  { return _aux_scalar_value[_T_last_]; }



  /**
   * @return the electron density at previous time step
   */
  virtual PetscScalar         n_last()          const
    { return _aux_scalar_value[_n_last_]; }

  /**
   * @return the writable reference to electron density at previous time step
   */
  virtual PetscScalar &       n_last()
  { return _aux_scalar_value[_n_last_]; }



  /**
   * @return the hole density at previous time step
   */
  virtual PetscScalar         p_last()          const
    { return _aux_scalar_value[_p_last_]; }

  /**
   * @return the writable reference to hole density at previous time step
   */
  virtual PetscScalar &       p_last()
  { return _aux_scalar_value[_p_last_]; }


  /**
   * @return the electron temperature at previous time step
   */
  virtual PetscScalar         Tn_last()          const
    { return _aux_scalar_value[_Tn_last_]; }

  /**
   * @return the writable reference to electron temperature at previous time step
   */
  virtual PetscScalar &       Tn_last()
  { return _aux_scalar_value[_Tn_last_]; }

  /**
   * @return the hole temperature at previous time step
   */
  virtual PetscScalar         Tp_last()          const
    { return _aux_scalar_value[_Tp_last_]; }

  /**
   * @return the writable reference to hole temperature at previous time step
   */
  virtual PetscScalar &       Tp_last()
  { return _aux_scalar_value[_Tp_last_]; }


  /**
   * @return the quantum conduction band at previous time step
   */
  virtual PetscScalar         Eqc_last()          const
    { return _aux_scalar_value[_Eqc_last_]; }

  /**
   * @return the writable reference to quantum conduction band at previous time step
   */
  virtual PetscScalar &       Eqc_last()
  { return _aux_scalar_value[_Eqc_last_]; }



  /**
   * @return the quantum valence band at previous time step
   */
  virtual PetscScalar         Eqv_last()          const
    { return _aux_scalar_value[_Eqv_last_]; }

  /**
   * @return the writable reference to quantum valence band at previous time step
   */
  virtual PetscScalar &       Eqv_last()
  { return _aux_scalar_value[_Eqv_last_]; }




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
   * @return the conduction band
   */
  virtual PetscScalar         Ec()          const
  { return _aux_scalar_value[_Ec_]; }

  /**
   * @return the writable reference to the conduction band
   */
  virtual PetscScalar &       Ec()
  { return _aux_scalar_value[_Ec_]; }


  /**
   * @return the valance band
   */
  virtual PetscScalar         Ev()          const
  { return _aux_scalar_value[_Ev_]; }

  /**
   * @return the writable reference to the valance band
   */
  virtual PetscScalar &       Ev()
  { return _aux_scalar_value[_Ev_]; }


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
   * @return the effective density of states in the conduction band
   */
  virtual PetscScalar         Nc()          const
    { return _aux_scalar_value[_Nc_]; }

  /**
   * @return the writable reference to the effective density of states in the conduction band
   */
  virtual PetscScalar &       Nc()
  { return _aux_scalar_value[_Nc_]; }


  /**
   * @return the effective density of states in the valence band
   */
  virtual PetscScalar         Nv()          const
    { return _aux_scalar_value[_Nv_]; }

  /**
   * @return the writable reference to the effective density of states in the valence band
   */
  virtual PetscScalar &       Nv()
  { return _aux_scalar_value[_Nv_]; }

  /**
   * @return the bandgap
   */
  virtual PetscScalar         Eg()          const
    { return _aux_scalar_value[_Eg_]; }

  /**
   * @return the writable reference to the bandgap
   */
  virtual PetscScalar &       Eg()
  { return _aux_scalar_value[_Eg_]; }



  /**
   * @return the mole fraction for single compound material
   */
  virtual PetscScalar         mole_x()          const
    { return _aux_scalar_value[_mole_x_]; }

  /**
   * @return the writable reference to the mole fraction for single compound material
   */
  virtual PetscScalar &       mole_x()
  { return _aux_scalar_value[_mole_x_]; }



  /**
  * @return the mole fraction for dual compound material
  */
  virtual PetscScalar         mole_y()          const
    { return _aux_scalar_value[_mole_y_]; }

  /**
   * @return the writable reference to the mole fraction for dual compound material
   */
  virtual PetscScalar &       mole_y()
  { return _aux_scalar_value[_mole_y_]; }


  /**
   * @return the general doping concentration of acceptor
   */
  virtual PetscScalar         Na()          const
    { return _aux_scalar_value[_Na_]; }

  /**
   * @return the writable reference to the general doping concentration of acceptor
   */
  virtual PetscScalar &       Na()
  { return _aux_scalar_value[_Na_]; }

  /**
   * @return the general doping concentration of donor
   */
  virtual PetscScalar         Nd()          const
    { return _aux_scalar_value[_Nd_]; }

  /**
   * @return the writable reference to the general doping concentration of donor
   */
  virtual PetscScalar &       Nd()
  { return _aux_scalar_value[_Nd_]; }


  /**
   * @return the concentration of donor atom phosphorus
   */
  virtual PetscScalar         P()          const
    { return _aux_scalar_value[_P_]; }

  /**
   * @return the writable reference to concentration of donor atom phosphorus
   */
  virtual PetscScalar &       P()
  { return _aux_scalar_value[_P_]; }

  /**
   * @return the concentration of donor atom arsenic
   */
  virtual PetscScalar         As()          const
    { return _aux_scalar_value[_As_]; }

  /**
   * @return the writable reference to the concentration of donor atom arsenic
   */
  virtual PetscScalar &       As()
  { return _aux_scalar_value[_As_]; }


  /**
   * @return the concentration of donor atom antimony
   */
  virtual PetscScalar         Sb()          const
    { return _aux_scalar_value[_Sb_]; }

  /**
   * @return the writable reference to concentration of donor atom antimony
   */
  virtual PetscScalar &       Sb()
  { return _aux_scalar_value[_Sb_]; }

  /**
   * @return the concentration of acceptor atom boron
   */
  virtual PetscScalar         B()          const
    { return _aux_scalar_value[_B_]; }

  /**
   * @return the writable reference to the concentration of acceptor atom boron
   */
  virtual PetscScalar &       B()
  { return _aux_scalar_value[_B_]; }


  /**
   * @return the total acceptor concentration
   */
  virtual PetscScalar         Total_Na()     const
    { return _aux_scalar_value[_Na_] + _aux_scalar_value[_B_]; }

  /**
   * @return the total donor concentration
   */
  virtual PetscScalar         Total_Nd()     const
    {return _aux_scalar_value[_Nd_] + _aux_scalar_value[_P_] + _aux_scalar_value[_As_] + _aux_scalar_value[_Sb_];}

  /**
   * @return net concentration
   */
  virtual PetscScalar         Net_doping()   const
    {return Total_Nd()-Total_Na();}

  /**
   * @return the total donor concentration
   */
  virtual PetscScalar         Total_doping() const
    {return Total_Nd()+Total_Na();}

  /**
   * @return net charge concentration
   */
  virtual PetscScalar         Net_charge()   const
    {return (Total_Nd()-n()) + (p()-Total_Na());}

  /**
   * @return intrinsic carrier concentration.
   * @note will not consider bandgap narrowing
   */
  virtual PetscScalar         ni()           const;


  /**
   * @return the quasi-fermi potential of electron
   */
  virtual PetscScalar         qFn()           const;

  /**
   * @return the quasi-fermi potential of hole
   */
  virtual PetscScalar         qFp()           const;



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



  /**
   * @return the carrier generation ratio due to OptG and PatG
   */
  virtual PetscScalar         Field_G()          const
  { return _aux_scalar_value[_Field_G_]; }

  /**
   * @return the writable carrier generation ratio due to OptG and PatG
   */
  virtual PetscScalar &       Field_G()
  { return _aux_scalar_value[_Field_G_]; }

  /**
   * @return the optical generation ratio
   */
  virtual  PetscScalar OptG()       const
  { return _aux_scalar_value[_OptG_];}

  /**
   * @return the writable optical generation ratio
   */
  virtual  PetscScalar & OptG()
  { return _aux_scalar_value[_OptG_];}

  /**
   * @return the heat generation ratio due to optical incident
   */
  virtual PetscScalar         OptQ()          const
  { return _aux_scalar_value[_OptQ_]; }

  /**
   * @return the writable heat generation ratio due to optical incident
   */
  virtual PetscScalar &       OptQ()
  { return _aux_scalar_value[_OptQ_]; }

  /**
   * @return the particle generation ratio
   */
  virtual PetscScalar         PatG()          const
  { return _aux_scalar_value[_PatG_]; }

  /**
   * @return the writable particle generation ratio
   */
  virtual PetscScalar &       PatG()
  { return _aux_scalar_value[_PatG_]; }

};


#endif
