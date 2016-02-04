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



#ifndef __fvm_node_data_insulator_h__
#define __fvm_node_data_insulator_h__


#include "fvm_node_data.h"




/**
 *  FVM nodal data for insulator region
 */
class FVM_Insulator_NodeData : public FVM_NodeData
{

  public:

    /**
     * the independent variable for insulator region
     */
    enum   InsulatorData
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
       * H+ density
       */
      _HIon_,

      /**
       * A type trap density
       */
      _trap_a_,

      /**
       * B type trap density
       */
      _trap_b_,
      
      /**
       * neutral state B type trap density
       */
      _trap_bn_,

      /**
       * electrostatic potential
       */
      _psi_,

      /**
       * lattice temperature
       */
      _T_,

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
       * band gap
       */
      _Eg_,

      /**
       * the dielectric permittivity
       */
      _eps_,

      /**
       * the megnetic permeability
       */
      _mu_,

      /**
       * the _OptG_*time +  _PatG_*time
       */
      _Field_G_,

      /**
       * energy deposite of high energy particle
       */
      _PatE_,

      /**
       * dose rate
       */
      _DoseRate_,

      /**
       * electrostatic potential at previous time step
       */
      _psi_last_,

      /**
       * electron density at previous time step
       */
      _n_last_,

      /**
       * hole density at previous time step
       */
      _p_last_,

      /**
       * old electrostatic potential
       */
      _psi_old_,

      /**
       * lattice temperature at previous time step
       */
      _T_last_,

      /**
       * last enum number
       */
      ScalarDataCount
    };


    /**
     * the vector auxiliary variable for insulator region
     */
    enum InsulatorAuxVecData
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
       * last enum number
       */
      VectorDataCount
    };

    /**
     * the complex auxiliary variable for insulator region
     */
    enum InsulatorAuxComplexData
    {

      /**
       * electrostatic potential
       */
      _psi_ac_=0,

      /**
       * lattice temperature
       */
      _T_ac_,

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


  public:
    /**
     * constructor
     */
    FVM_Insulator_NodeData ( DataStorage * data_storage , const std::map<std::string, SimulationVariable> & variables )
    :FVM_NodeData ( data_storage , variables)
    {}

    /**
     * destructor
     */
    virtual ~FVM_Insulator_NodeData()  {}

  public:
    /**
     * @return the solution variable number
     */
    static size_t n_scalar()
    { return static_cast<unsigned int> ( ScalarDataCount ) ; /* return last enum */ }

    /**
     * @return the complex variable number
     */
    static size_t n_complex()
    { return static_cast<unsigned int> ( ComplexDataCount ) ; /* return last enum */ }

    /**
     * @return the vector variable number
     */
    static size_t n_vector()
    { return static_cast<unsigned int> ( VectorDataCount ) ; /* return last enum */ }

    /**
     * @return the tensor variable number
     */
    static size_t n_tensor()
    { return 0; }

    /**
     * @return the data type
     */
    virtual NodeDataType type() const
    { return FVM_NodeData::InsulatorData; }


  public:

    /**
     * @return data by enum name
     */
    virtual PetscScalar  get_variable_real ( SolutionVariable variable ) const
    {
      switch ( variable )
      {
        case POTENTIAL   :  return  psi();                            /* potential */
        case ELECTRON    :  return  n();                              /* electron concentration */
        case HOLE        :  return  p();                              /* hole concentration */
        case TEMPERATURE :  return  T();                              /* lattice temperature */
        case E_TEMP      :  return  T();                              /* electron temperature */
        case H_TEMP      :  return  T();                              /* hole temperature */
        case QFN         :  return  psi();                            /* electron quasi-Fermi level */
        case QFP         :  return  psi();                            /* hole quasi-Fermi level */
        default          :  return  0.0;

      }
    }

    /**
     * set variable by enum name
     */
    virtual void set_variable_real ( SolutionVariable variable, PetscScalar value )
    {
      switch ( variable )
      {
        case POTENTIAL   :  psi() = value;                             /* potential */
        case ELECTRON    :  n() = value;                               /* electron concentration */
        case HOLE        :  p() = value;                               /* hole concentration */
        case TEMPERATURE :  T() = value;                               /* lattice temperature */
        default          :  return;
      }
    }

    /**
     * @return true when this variable valid
     */
    virtual bool is_variable_valid ( SolutionVariable variable )  const
    {
      switch ( variable )
      {
        case POTENTIAL   :
        case TEMPERATURE :  return true;
        default          :  return false;
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
     * @return the H+ densiy
     */
    virtual PetscScalar         HIon()          const
    { return _data_storage->scalar ( _HIon_, _offset ); }

    /**
     * @return the writable reference to H+ densiy
     */
    virtual PetscScalar &       HIon()
    { return _data_storage->scalar ( _HIon_, _offset ); }
    
    /**
     * @return the A type trap density
     */
    virtual PetscScalar         trap_a()          const
    { return _data_storage->scalar ( _trap_a_, _offset ); }

    /**
     * @return the writable reference to A type trap density
     */
    virtual PetscScalar &       trap_a()
    { return _data_storage->scalar ( _trap_a_, _offset ); }



    /**
     * @return the B type trap density
     */
    virtual PetscScalar         trap_b()          const
    { return _data_storage->scalar ( _trap_b_, _offset ); }


    /**
     * @return the writable reference to B type trap density
     */
    virtual PetscScalar &       trap_b()
    { return _data_storage->scalar ( _trap_b_, _offset ); }

    
    /**
     * @return the neutral state B type trap density
     */
    virtual PetscScalar         trap_bn()          const
    { return _data_storage->scalar ( _trap_bn_, _offset ); }


    /**
     * @return the writable reference to neutral state B type trap density
     */
    virtual PetscScalar &       trap_bn()
    { return _data_storage->scalar ( _trap_bn_, _offset ); }
  

    /**
     * @return the lattice temperature
     */
    virtual PetscScalar         T()          const
    { return _data_storage->scalar ( _T_, _offset ); }

    /**
     * @return the lattice temperature
     */
    virtual PetscScalar &       T()
    { return _data_storage->scalar ( _T_, _offset ); }


    /**
     * @return the electron temperature, the same as lattice temperature
     */
    virtual PetscScalar         Tn()          const
    { return _data_storage->scalar ( _T_, _offset ); }


    /**
     * @return the writable reference to electron temperature, the same as lattice temperature
     */
    virtual PetscScalar &       Tn()
    { return _data_storage->scalar ( _T_, _offset ); }


    /**
     * @return the hole temperature, the same as lattice temperature
     */
    virtual PetscScalar         Tp()          const
    { return _data_storage->scalar ( _T_, _offset ); }


    /**
     * @return the writable reference to hole temperature, the same as lattice temperature
     */
    virtual PetscScalar &       Tp()
    { return _data_storage->scalar ( _T_, _offset ); }



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
     * @return the particle energy
     */
    virtual PetscScalar         PatE()          const
    { return _data_storage->scalar ( _PatE_, _offset ); }

    /**
     * @return the writable particle energy
     */
    virtual PetscScalar &       PatE()
    { return _data_storage->scalar ( _PatE_, _offset ); }


   /**
    * @return the dose rate
    */
    virtual PetscScalar         DoseRate()          const
    { return _data_storage->scalar ( _DoseRate_, _offset ); }

   /**
    * @return the writable reference to dose rate
    */
    virtual PetscScalar &       DoseRate()
    { return _data_storage->scalar ( _DoseRate_, _offset ); }


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
     * @return the conduction band
     */
    virtual PetscScalar         Ec()          const;

    /**
     * @return the writable reference to the conduction band
     */
    virtual PetscScalar &       Ec()
    { return _data_storage->scalar ( _Ec_, _offset ); }

    /**
     * @return the valance band
     */
    virtual PetscScalar         Ev()          const;

    /**
     * @return the writable reference to the valance band
     */
    virtual PetscScalar &       Ev()
    { return _data_storage->scalar ( _Ev_, _offset ); }
    
    
    /**
     * @return the quantum conduction band
     */
    virtual PetscScalar         Eqc()          const
    { return Ec(); }
   
    
    /**
     * @return the quantum valence band
     */
    virtual PetscScalar         Eqv()          const
    { return Ev(); }


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
    { return _data_storage->vector ( _E_, _offset );}


    /**
     * @return the writable reference to electrical field
     */
    virtual VectorValue<PetscScalar> & E()
    { return _data_storage->vector ( _E_, _offset );}


};




#endif
