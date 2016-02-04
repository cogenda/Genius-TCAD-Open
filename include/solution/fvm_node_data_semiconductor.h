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
       * minimal distance
       */
      _dmin_ = 0,

      /**
       * electron density
       */
      _n_,

      /**
       * hole density
       */
      _p_,

      /**
       * electrostatic potential
       */
      _psi_,

      /**
       * correction of potential
       */
      _dpsi_,

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
       * quantum conduction band, used in density gradient model simulation
       */
      _Eqc_,

      /**
       * quantum valence band, used in density gradient model simulation
       */
      _Eqv_,

      /**
       * the density of the material
       */
      _density_,

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
       * conduction band change
       */
      _dEcStrain_,

      /**
       * valence band change
       */
      _dEvStrain_,

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
      _mun_,

      /**
       * hole mobility
       */
      _mup_,

      /**
       * total recombination
       */
      _Recomb_,

      /**
       * direct(optical) recombination
       */
      _Recomb_Dir_,

      /**
       * SRH recombination
       */
      _Recomb_SRH_,

      /**
       * Auger recombination
       */
      _Recomb_Auger_,

      /**
       * impact ionization
       */
      _ImpactIonization_,

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
       * energy deposite of incident wave
       */
      _OptE_,

      /**
       * carrier generation due to high energy particle
       */
      _PatG_,

      /**
       * energy deposite of high energy particle
       */
      _PatE_,

      /**
       * electron inject ratio
       */
      _EIn_,

      /**
       * hole inject ratio
       */
      _HIn_,

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
       * old electron density
       */
      _n_old_,

      /**
       * old hole density
       */
      _p_old_,

      /**
       * old electrostatic potential
       */
      _psi_old_,

      /**
       * charge density
       */
      _rho_,

      /**
       * interface charge density
       */
      _Nit_,

      /**
       * last enum number
       */
      ScalarDataCount
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
      _Jp_,

      /**
       * last enum number
       */
      VectorDataCount
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
      _OpH_complex_,

      /**
       * last enum number
       */
      ComplexDataCount
    };

    /**
     * the tensor auxiliary variable for semiconductor region
     */
    enum SemiconductorAuxTensorData
    {
      /**
       * stress
       */
      _stress_=0,

      /**
       * strain
       */
      _strain_,

      /**
       * last enum number
       */
      TensorDataCount
    };


  public:

    /**
     * constructor
     */
    FVM_Semiconductor_NodeData ( DataStorage * data_storage , const std::map<std::string, SimulationVariable> & variables )
    :FVM_NodeData ( data_storage , variables)
    {}

    /**
     * destructor
     */
    virtual ~FVM_Semiconductor_NodeData()  { }

  public:
    /**
     * @return the solution variable number
     */
    static size_t n_scalar()
    { return static_cast<unsigned int> ( ScalarDataCount ) ; /* return last enum*/ }

    /**
     * @return the complex variable number
     */
    static size_t n_complex()
    { return static_cast<unsigned int> ( ComplexDataCount ) ; /* return last enum*/ }

    /**
     * @return the vector variable number
     */
    static size_t n_vector()
    { return static_cast<unsigned int> ( VectorDataCount ) ; /* return last enum*/ }

    /**
     * @return the tensor variable number
     */
    static size_t n_tensor()
    { return static_cast<unsigned int> ( TensorDataCount ) ; /* return last enum*/ }

    /**
     * @return the data type
     */
    virtual NodeDataType type() const
    { return FVM_NodeData::SemiconductorData; }


  public:

    /**
     * @return data by enum name
     */
    virtual PetscScalar  get_variable_real ( SolutionVariable variable ) const
    {
      switch ( variable )
      {
        case POTENTIAL      :  return  psi();                            /* potential */
        case ELECTRON       :  return  n();                              /* electron concentration */
        case HOLE           :  return  p();                              /* hole concentration */
        case TEMPERATURE    :  return  T();                              /* lattice temperature */
        case E_TEMP         :  return  Tn();                             /* electron temperature */
        case H_TEMP         :  return  Tp();                             /* hole temperature */
        case EQC            :  return  Eqc();                            /* quantum conduction band */
        case EQV            :  return  Eqv();                            /* quantum valence band */

        case DOPING         :  return  Net_doping();                     /* net doping */
        case DOPING_Na      :  return  Total_Na();                       /* acceptor */
        case DOPING_Nd      :  return  Total_Nd();                       /* donor */

        case SPECIES_Na     :  return Na();                              /* generic species Na */
        case SPECIES_Nd     :  return Nd();                              /* generic species Na */


        case MIN_CARRIER    :  return  Net_doping() > 0 ? n() : p();     /* minority carrier concentration */
        case NET_CARRIER    :  return  p() - n();                        /* net carrier concentration */
        case NET_CHARGE     :  return  Net_doping() + p() - n();         /* net charge */

        case RECOMBINATION  :  return  Recomb();                         /* total recombination */
        case RECOMB_DIR     :  return  Recomb_Dir();                     /* direct recombination */
        case RECOMB_SHR     :  return  Recomb_SRH();                     /* SHR recombination */
        case RECOMB_AUGER   :  return  Recomb_Auger();                   /* AUGER recombination */

        case MOLE_X         :  return  mole_x();
        case MOLE_Y         :  return  mole_y();

        case OPTICAL_GEN    :  return  OptG();                           /* charge genetated by optical ray */
        case OPTICAL_HEAT   :  return  OptQ();                           /* heat genetated by optical ray */
        case PARTICLE_GEN   :  return  PatG();                           /* charge genetated by particle ray */

        case QFN            :  return  qFn();                            /* electron quasi-Fermi level */
        case QFP            :  return  qFp();                            /* hole quasi-Fermi level */

        default             :  return  0.0;
      }
    }

    /**
     * set variable by enum name
     */
    virtual void set_variable_real ( SolutionVariable variable, PetscScalar value )
    {
      switch ( variable )
      {
        case POTENTIAL     :  psi() = value;                             /* potential */
        case ELECTRON      :  n() = value;                               /* electron concentration */
        case HOLE          :  p() = value;                               /* hole concentration */
        case TEMPERATURE   :  T() = value;                               /* lattice temperature */
        case E_TEMP        :  Tn() = value;                              /* electron temperature */
        case H_TEMP        :  Tp() = value;                              /* hole temperature */
        case DOPING_Na     :  Na() = value;                              /* acceptor */
        case DOPING_Nd     :  Nd() = value;                              /* donor */
        case OPTICAL_GEN   :  OptG() = value;                            /* charge genetated by optical ray */
        case OPTICAL_HEAT  :  OptQ() = value;                            /* heat genetated by optical ray */
        case PARTICLE_GEN  :  PatG() = value;                            /* charge genetated by particle ray */
        case MOLE_X        :  mole_x() = value;
        case MOLE_Y        :  mole_y() = value;
        default            :  return;
      }
    }

    /**
     * @return true when this variable valid
     */
    virtual bool is_variable_valid ( SolutionVariable variable )  const
    {
      switch ( variable )
      {
        case POTENTIAL     :
        case ELECTRON      :
        case HOLE          :
        case TEMPERATURE   :
        case E_TEMP        :
        case H_TEMP        :
        case DOPING_Na     :
        case DOPING_Nd     :
        case OPTICAL_GEN   :
        case OPTICAL_HEAT  :
        case PARTICLE_GEN  :
        case MOLE_X        :
        case MOLE_Y        :  return true;
        default            :  return false;
      }
    }

    //--------------------------------------------------------------------
    // data access function
    //--------------------------------------------------------------------


    /**
     * @return the minimal distange to surface
     */
    virtual PetscScalar         dmin()          const
    { return _data_storage->scalar ( _dmin_, _offset ); }

    /**
     * @return the writable reference to the minimal distange to surface
     */
    virtual PetscScalar &       dmin()
    { return _data_storage->scalar ( _dmin_, _offset ); }


    /**
     * @return the statistic potential
     */
    virtual PetscScalar         psi()        const
    { return _data_storage->scalar ( _psi_, _offset ); }

    /**
     * @return the statistic potential
     */
    virtual PetscScalar &       psi()
    { return _data_storage->scalar ( _psi_, _offset ); }




    /**
     * @return the correction of potential
     */
    virtual PetscScalar         dpsi()          const
    { return _data_storage->scalar ( _dpsi_, _offset ); }

    /**
     * @return the writable reference to correction of potential
     */
    virtual PetscScalar &       dpsi()
    { return _data_storage->scalar ( _dpsi_, _offset ); }


    /**
     * @return the lattice temperature
     */
    virtual PetscScalar         T()          const
    { return _data_storage->scalar ( _T_, _offset ); }

    /**
     * @return the statistic potential
     */
    virtual PetscScalar &       T()
    { return _data_storage->scalar ( _T_, _offset ); }



    /**
     * @return the electron density
     */
    virtual PetscScalar         n()          const
    { return _data_storage->scalar ( _n_, _offset ); }

    /**
     * @return the writable reference to electron density
     */
    virtual PetscScalar &       n()
    { return _data_storage->scalar ( _n_, _offset ); }



    /**
     * @return the hole density
     */
    virtual PetscScalar         p()          const
    { return _data_storage->scalar ( _p_, _offset ); }


    /**
     * @return the writable reference to hole density
     */
    virtual PetscScalar &       p()
    { return _data_storage->scalar ( _p_, _offset ); }


    /**
     * @return the electron temperature
     */
    virtual PetscScalar         Tn()          const
    { return _data_storage->scalar ( _Tn_, _offset ); }


    /**
     * @return the writable reference to electron temperature
     */
    virtual PetscScalar &       Tn()
    { return _data_storage->scalar ( _Tn_, _offset ); }


    /**
     * @return the hole temperature
     */
    virtual PetscScalar         Tp()          const
    { return _data_storage->scalar ( _Tp_, _offset ); }


    /**
     * @return the writable reference to hole temperature
     */
    virtual PetscScalar &       Tp()
    { return _data_storage->scalar ( _Tp_, _offset ); }


    /**
     * @return the quantum conduction band
     */
    virtual PetscScalar         Eqc()          const
    { return _data_storage->scalar ( _Eqc_, _offset ); }

    /**
     * @return the writable reference to quantum conduction band
     */
    virtual PetscScalar &       Eqc()
    { return _data_storage->scalar ( _Eqc_, _offset ); }

    /**
     * @return the quantum valence band
     */
    virtual PetscScalar         Eqv()          const
    { return _data_storage->scalar ( _Eqv_, _offset ); }

    /**
     * @return the writable reference to quantum valence band
     */
    virtual PetscScalar &       Eqv()
    { return _data_storage->scalar ( _Eqv_, _offset ); }


    /**
     * @return the statistic potential
     */
    virtual std::complex<PetscScalar>         psi_ac()          const
    { return _data_storage->complex ( _psi_ac_, _offset ); }

    /**
     * @return the writable reference to statistic potential
     */
    virtual std::complex<PetscScalar> &       psi_ac()
    { return _data_storage->complex ( _psi_ac_, _offset ); }


    /**
     * @return the lattice temperature
     */
    virtual std::complex<PetscScalar>         T_ac()          const
    { return  _data_storage->complex ( _T_ac_, _offset ); }

    /**
     * @return the writable reference to lattice temperature
     */
    virtual std::complex<PetscScalar> &       T_ac()
    { return _data_storage->complex ( _T_ac_, _offset ); }



    /**
     * @return the electron density
     */
    virtual std::complex<PetscScalar>         n_ac()          const
    { return _data_storage->complex ( _n_ac_, _offset ); }

    /**
     * @return the writable reference to electron density
     */
    virtual std::complex<PetscScalar> &       n_ac()
    { return _data_storage->complex ( _n_ac_, _offset ); }



    /**
     * @return the hole density
     */
    virtual std::complex<PetscScalar>         p_ac()          const
    { return _data_storage->complex ( _p_ac_, _offset ); }

    /**
     * @return the writable reference to hole density
     */
    virtual std::complex<PetscScalar> &       p_ac()
    { return _data_storage->complex ( _p_ac_, _offset ); }


    /**
     * @return the electron temperature
     */
    virtual std::complex<PetscScalar>         Tn_ac()          const
    { return _data_storage->complex ( _Tn_ac_, _offset ); }

    /**
     * @return the writable reference to electron temperature
     */
    virtual std::complex<PetscScalar> &       Tn_ac()
    { return _data_storage->complex ( _Tn_ac_, _offset ); }

    /**
     * @return the hole temperature
     */
    virtual std::complex<PetscScalar>         Tp_ac()          const
    { return _data_storage->complex ( _Tp_ac_, _offset ); }

    /**
     * @return the writable reference to hole temperature
     */
    virtual std::complex<PetscScalar> &       Tp_ac()
    { return _data_storage->complex ( _Tp_ac_, _offset ); }


    /**
     * @return the complex E file. only used by EM FEM solver
     */
    virtual std::complex<PetscScalar>         OptE_complex()          const
    { return _data_storage->complex ( _OpE_complex_, _offset ); }

    /**
     * @return the writable reference to complex E file. only used by EM FEM solver
     */
    virtual std::complex<PetscScalar> &       OptE_complex()
    { return _data_storage->complex ( _OpE_complex_, _offset ); }

    /**
     * @return the complex H file. only used by EM FEM solver
     */
    virtual std::complex<PetscScalar>         OptH_complex()          const
    { return _data_storage->complex ( _OpH_complex_, _offset ); }

    /**
     * @return the writable reference to complex H file. only used by EM FEM solver
     */
    virtual std::complex<PetscScalar> &       OptH_complex()
    { return _data_storage->complex ( _OpH_complex_, _offset ); }


    /**
     * @return the statistic potential at previous time step
     */
    virtual PetscScalar         psi_last()          const
    { return _data_storage->scalar ( _psi_last_, _offset ); }

    /**
     * @return the writable reference to statistic potential at previous time step
     */
    virtual PetscScalar &       psi_last()
    { return _data_storage->scalar ( _psi_last_, _offset ); }


    /**
     * @return the old statistic potential
     */
    virtual PetscScalar         psi_old()          const
    { return _data_storage->scalar ( _psi_old_, _offset ); }

    /**
     * @return the writable reference to old statistic potential
     */
    virtual PetscScalar &       psi_old()
    { return _data_storage->scalar ( _psi_old_, _offset ); }



    /**
     * @return the lattice temperature at previous time step
     */
    virtual PetscScalar         T_last()          const
    { return _data_storage->scalar ( _T_last_, _offset ); }

    /**
     * @return the writable reference to lattice temperature at previous time step
     */
    virtual PetscScalar &       T_last()
    { return _data_storage->scalar ( _T_last_, _offset ); }



    /**
     * @return the electron density at previous time step
     */
    virtual PetscScalar         n_last()          const
    { return _data_storage->scalar ( _n_last_, _offset ); }

    /**
     * @return the writable reference to electron density at previous time step
     */
    virtual PetscScalar &       n_last()
    { return _data_storage->scalar ( _n_last_, _offset ); }

    /**
     * @return the old electron density
     */
    virtual PetscScalar         n_old()          const
    { return _data_storage->scalar ( _n_old_, _offset ); }

    /**
     * @return the writable reference to old electron density
     */
    virtual PetscScalar &       n_old()
    { return _data_storage->scalar ( _n_old_, _offset ); }


    /**
     * @return the hole density at previous time step
     */
    virtual PetscScalar         p_last()          const
    { return _data_storage->scalar ( _p_last_, _offset ); }

    /**
     * @return the writable reference to hole density at previous time step
     */
    virtual PetscScalar &       p_last()
    { return _data_storage->scalar ( _p_last_, _offset ); }

    /**
     * @return the old hole density
     */
    virtual PetscScalar         p_old()          const
    { return _data_storage->scalar ( _p_old_, _offset ); }

    /**
     * @return the writable reference to old hole density
     */
    virtual PetscScalar &       p_old()
    { return _data_storage->scalar ( _p_old_, _offset ); }


    /**
     * @return the electron temperature at previous time step
     */
    virtual PetscScalar         Tn_last()          const
    { return _data_storage->scalar ( _Tn_last_, _offset ); }

    /**
     * @return the writable reference to electron temperature at previous time step
     */
    virtual PetscScalar &       Tn_last()
    { return _data_storage->scalar ( _Tn_last_, _offset ); }

    /**
     * @return the hole temperature at previous time step
     */
    virtual PetscScalar         Tp_last()          const
    { return _data_storage->scalar ( _Tp_last_, _offset ); }

    /**
     * @return the writable reference to hole temperature at previous time step
     */
    virtual PetscScalar &       Tp_last()
    { return _data_storage->scalar ( _Tp_last_, _offset ); }




    /**
     * @return the mass density of the material
     */
    virtual PetscScalar         density()          const
    { return _data_storage->scalar ( _density_, _offset ); }

    /**
     * @return the writable reference to the mass density of the material
     */
    virtual PetscScalar &       density()
    { return _data_storage->scalar ( _density_, _offset ); }




    /**
     * @return the electron affinity
     */
    virtual PetscScalar         affinity()          const
    { return _data_storage->scalar ( _affinity_, _offset ); }

    /**
     * @return the writable reference to the electron affinity
     */
    virtual PetscScalar &       affinity()
    { return _data_storage->scalar ( _affinity_, _offset ); }


    /**
     * @return the conduction band
     */
    virtual PetscScalar         Ec()          const
    { return _data_storage->scalar ( _Ec_, _offset ); }

    /**
     * @return the writable reference to the conduction band
     */
    virtual PetscScalar &       Ec()
    { return _data_storage->scalar ( _Ec_, _offset ); }


    /**
     * @return the valance band
     */
    virtual PetscScalar         Ev()          const
    { return _data_storage->scalar ( _Ev_, _offset ); }

    /**
     * @return the writable reference to the valance band
     */
    virtual PetscScalar &       Ev()
    { return _data_storage->scalar ( _Ev_, _offset ); }


    /**
     * @return the change of conduction band due to strain
     */
    virtual PetscScalar         dEcStrain()          const
    { return _data_storage->scalar ( _dEcStrain_, _offset ); }

    /**
     * @return the writable reference to the change of conduction band due to strain
     */
    virtual PetscScalar &       dEcStrain()
    { return _data_storage->scalar ( _dEcStrain_, _offset ); }

    /**
     * @return the change of valance band due to strain
     */
    virtual PetscScalar         dEvStrain()          const
    { return _data_storage->scalar ( _dEvStrain_, _offset ); }

    /**
     * @return the writable reference to the change of valance band due to strain
     */
    virtual PetscScalar &       dEvStrain()
    { return _data_storage->scalar ( _dEvStrain_, _offset ); }

    /**
     * @return the dielectric permittivity
     */
    virtual PetscScalar         eps()          const
    { return _data_storage->scalar ( _eps_, _offset ); }

    /**
     * @return the writable reference to the dielectric permittivity
     */
    virtual PetscScalar &       eps()
    { return _data_storage->scalar ( _eps_, _offset ); }


    /**
     * @return the megnetic permeability
     */
    virtual PetscScalar         mu()          const
    { return _data_storage->scalar ( _mu_, _offset ); }

    /**
     * @return the writable reference to the megnetic permeability
     */
    virtual PetscScalar &       mu()
    { return _data_storage->scalar ( _mu_, _offset ); }


    /**
     * @return the effective density of states in the conduction band
     */
    virtual PetscScalar         Nc()          const
    { return _data_storage->scalar ( _Nc_, _offset ); }

    /**
     * @return the writable reference to the effective density of states in the conduction band
     */
    virtual PetscScalar &       Nc()
    { return _data_storage->scalar ( _Nc_, _offset ); }


    /**
     * @return the effective density of states in the valence band
     */
    virtual PetscScalar         Nv()          const
    { return _data_storage->scalar ( _Nv_, _offset ); }

    /**
     * @return the writable reference to the effective density of states in the valence band
     */
    virtual PetscScalar &       Nv()
    { return _data_storage->scalar ( _Nv_, _offset ); }

    /**
     * @return the bandgap
     */
    virtual PetscScalar         Eg()          const
    { return _data_storage->scalar ( _Eg_, _offset ); }

    /**
     * @return the writable reference to the bandgap
     */
    virtual PetscScalar &       Eg()
    { return _data_storage->scalar ( _Eg_, _offset ); }



    /**
     * @return the mole fraction for single compound material
     */
    virtual PetscScalar         mole_x()          const
    { return _data_storage->scalar ( _mole_x_, _offset ); }

    /**
     * @return the writable reference to the mole fraction for single compound material
     */
    virtual PetscScalar &       mole_x()
    { return _data_storage->scalar ( _mole_x_, _offset ); }



    /**
    * @return the mole fraction for dual compound material
    */
    virtual PetscScalar         mole_y()          const
    { return _data_storage->scalar ( _mole_y_, _offset ); }

    /**
     * @return the writable reference to the mole fraction for dual compound material
     */
    virtual PetscScalar &       mole_y()
    { return _data_storage->scalar ( _mole_y_, _offset ); }

    /**
     * @return the electron mobility
     */
    virtual PetscScalar         mun()          const
    { return _data_storage->scalar ( _mun_, _offset ); }

    /**
     * @return the writable reference to electron mobility
     */
    virtual PetscScalar &       mun()
    { return _data_storage->scalar ( _mun_, _offset ); }

    /**
     * @return the hole mobility
     */
    virtual PetscScalar         mup()          const
    { return _data_storage->scalar ( _mup_, _offset ); }

    /**
     * @return the writable reference to hole mobility
     */
    virtual PetscScalar &       mup()
    { return _data_storage->scalar ( _mup_, _offset ); }

    /**
     * @return the general doping concentration of acceptor
     */
    virtual PetscScalar         Na()          const
    { return _data_storage->scalar ( _Na_, _offset ); }

    /**
     * @return the writable reference to the general doping concentration of acceptor
     */
    virtual PetscScalar &       Na()
    { return _data_storage->scalar ( _Na_, _offset ); }

    /**
     * @return the general doping concentration of donor
     */
    virtual PetscScalar         Nd()          const
    { return _data_storage->scalar ( _Nd_, _offset ); }

    /**
     * @return the writable reference to the general doping concentration of donor
     */
    virtual PetscScalar &       Nd()
    { return _data_storage->scalar ( _Nd_, _offset ); }





    /**
     * @return the total acceptor concentration
     */
    virtual PetscScalar         Total_Na()     const
    {
      return _data_storage->scalar ( _Na_, _offset );
    }

    /**
     * @return the total donor concentration
     */
    virtual PetscScalar         Total_Nd()     const
    {
      return _data_storage->scalar ( _Nd_, _offset );
    }

    /**
     * @return net concentration
     */
    virtual PetscScalar         Net_doping()   const
    {return Total_Nd()-Total_Na();}

    /**
     * @return the total donor concentration
     */
    virtual PetscScalar         Total_doping() const
    {return Total_Nd() +Total_Na();}

    /**
     * @return net charge concentration
     */
    virtual PetscScalar         Net_charge()   const
    {return ( Total_Nd()-n() ) + ( p()-Total_Na() );}

    /**
     * @return intrinsic carrier concentration.
     * @note will not consider bandgap narrowing
     */
    virtual PetscScalar         ni()           const;


    /**
     * @return the quasi-fermi potential of electron
     */
    virtual PetscScalar         qFn()           const
    { return _data_storage->scalar ( _qFn_, _offset ); }

    /**
     * @return the quasi-fermi potential of electron
     */
    virtual PetscScalar &       qFn()
    { return _data_storage->scalar ( _qFn_, _offset ); }


    /**
     * @return the quasi-fermi potential of hole
     */
    virtual PetscScalar         qFp()           const
    { return _data_storage->scalar ( _qFp_, _offset ); }

    /**
     * @return the quasi-fermi potential of hole
     */
    virtual PetscScalar &       qFp()
    { return _data_storage->scalar ( _qFp_, _offset ); }

    /**
     * @return the charge density
     * this variable is used for data exchange
     * between hdm solver and poisson solver
     */
    virtual PetscScalar         rho()          const
    { return _data_storage->scalar ( _rho_, _offset ); }


    /**
     * @return the writable reference to charge density
     */
    virtual PetscScalar &       rho()
    { return _data_storage->scalar ( _rho_, _offset ); }


    /**
     * @return the interface charge density
     */
    virtual PetscScalar         interface_charge()          const
    { return _data_storage->scalar ( _Nit_, _offset ); }


    /**
     * @return the writable reference to interface charge density
     */
    virtual PetscScalar &       interface_charge()
    { return _data_storage->scalar ( _Nit_, _offset ); }



    /**
     * @return the recombnation rate
     */
    virtual PetscScalar         Recomb()          const
    { return _data_storage->scalar ( _Recomb_, _offset ); }

    /**
     * @return the writable reference to recombnation rate
     */
    virtual PetscScalar &       Recomb()
    { return _data_storage->scalar ( _Recomb_, _offset ); }


    /**
     * @return the direct(optical) recombnation rate
     */
    virtual PetscScalar         Recomb_Dir()          const
    { return _data_storage->scalar ( _Recomb_Dir_, _offset ); }

    /**
     * @return the writable reference to direct(optical) recombnation rate
     */
    virtual PetscScalar &       Recomb_Dir()
    { return _data_storage->scalar ( _Recomb_Dir_, _offset ); }


    /**
     * @return the SRH recombnation rate
     */
    virtual PetscScalar         Recomb_SRH()          const
    { return _data_storage->scalar ( _Recomb_SRH_, _offset ); }

    /**
     * @return the writable reference to SRH recombnation rate
     */
    virtual PetscScalar &       Recomb_SRH()
    { return _data_storage->scalar ( _Recomb_SRH_, _offset ); }


    /**
     * @return the Auger recombnation rate
     */
    virtual PetscScalar         Recomb_Auger()          const
    { return _data_storage->scalar ( _Recomb_Auger_, _offset ); }

    /**
     * @return the writable reference to Auger recombnation rate
     */
    virtual PetscScalar &       Recomb_Auger()
    { return _data_storage->scalar ( _Recomb_Auger_, _offset ); }

    /**
     * @return the impact ionization
     */
    virtual PetscScalar         ImpactIonization()          const
    { return _data_storage->scalar ( _ImpactIonization_, _offset ); }

    /**
     * @return the writable reference to impact ionization
     */
    virtual PetscScalar &       ImpactIonization()
    { return _data_storage->scalar ( _ImpactIonization_, _offset ); }


    /**
     * @return the electrical field
     */
    virtual VectorValue<PetscScalar> E()       const
    { return _data_storage->vector ( _E_, _offset );}


    /**
     * @return the writable reference to electrical field
     */
    virtual VectorValue<PetscScalar> & E()
    { return _data_storage->vector ( _E_, _offset );}


    /**
     * @return the electron current
     */
    virtual VectorValue<PetscScalar> Jn()       const
    { return _data_storage->vector ( _Jn_, _offset );}


    /**
     * @return the writable reference to electron current
     */
    virtual VectorValue<PetscScalar> & Jn()
    { return _data_storage->vector ( _Jn_, _offset );}


    /**
     * @return the hole current
     */
    virtual VectorValue<PetscScalar> Jp()       const
    { return _data_storage->vector ( _Jp_, _offset );}


    /**
     * @return the writable reference to hole current
     */
    virtual VectorValue<PetscScalar> & Jp()
    { return _data_storage->vector ( _Jp_, _offset );}


    /**
     * @return the carrier generation ratio due to OptG and PatG
     */
    virtual PetscScalar         Field_G()          const
    { return _data_storage->scalar ( _Field_G_, _offset ); }

    /**
     * @return the writable carrier generation ratio due to OptG and PatG
     */
    virtual PetscScalar &       Field_G()
    { return _data_storage->scalar ( _Field_G_, _offset ); }

    /**
     * @return the optical generation ratio
     */
    virtual  PetscScalar OptG()       const
    { return _data_storage->scalar ( _OptG_, _offset );}

    /**
     * @return the writable optical generation ratio
     */
    virtual  PetscScalar & OptG()
    { return _data_storage->scalar ( _OptG_, _offset );}

    /**
     * @return the heat generation ratio due to optical incident
     */
    virtual PetscScalar         OptQ()          const
    { return _data_storage->scalar ( _OptQ_, _offset ); }

    /**
     * @return the writable heat generation ratio due to optical incident
     */
    virtual PetscScalar &       OptQ()
    { return _data_storage->scalar ( _OptQ_, _offset ); }


    /**
     * @return the optical energy
     */
    virtual PetscScalar         OptE()          const
    { return _data_storage->scalar ( _OptE_, _offset ); }

    /**
     * @return the writable optical energy
     */
    virtual PetscScalar &       OptE()
    { return _data_storage->scalar ( _OptE_, _offset ); }

    /**
     * @return the particle generation ratio
     */
    virtual PetscScalar         PatG()          const
    { return _data_storage->scalar ( _PatG_, _offset ); }

    /**
     * @return the writable particle generation ratio
     */
    virtual PetscScalar &       PatG()
    { return _data_storage->scalar ( _PatG_, _offset ); }


    /**
     * @return the energy deposite of particle
     */
    virtual PetscScalar         PatE()          const
    { return _data_storage->scalar ( _PatE_, _offset ); }

    /**
     * @return the writable energy deposite of particle
     */
    virtual PetscScalar &       PatE()
    { return _data_storage->scalar ( _PatE_, _offset ); }



    /**
     * @return the electron injected in to the FVM cell.
     * NOTE: it is the flux flow into the FVM cell or total electron density generated in the cell
     */
    virtual PetscScalar         EIn()          const
    { return _data_storage->scalar ( _EIn_, _offset ); }

    /**
     * @return the writable reference of electron injected in to the FVM cell.
     */
    virtual PetscScalar &       EIn()
    { return _data_storage->scalar ( _EIn_, _offset ); }


    /**
     * @return the hole inject injected in to the FVM cell.
     */
    virtual PetscScalar         HIn()          const
    { return _data_storage->scalar ( _HIn_, _offset ); }

    /**
     * @return the writable reference of hole injected in to the FVM cell.
     */
    virtual PetscScalar &       HIn()
    { return _data_storage->scalar ( _HIn_, _offset ); }



    /**
     * @return the stress tensor
     */
    virtual TensorValue<PetscScalar> stress()       const
    { return _data_storage->tensor ( _stress_, _offset );}


    /**
     * @return the writable reference to stress tensor
     */
    virtual TensorValue<PetscScalar> & stress()
    { return _data_storage->tensor ( _stress_, _offset );}


    /**
     * @return the strain tensor
     */
    virtual TensorValue<PetscScalar> strain()       const
    { return _data_storage->tensor ( _strain_, _offset );}


    /**
     * @return the writable reference to strain tensor
     */
    virtual TensorValue<PetscScalar> & strain()
    { return _data_storage->tensor ( _strain_, _offset );}

};


#endif
