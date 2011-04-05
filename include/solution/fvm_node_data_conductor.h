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


#ifndef __fvm_node_data_conductor_h__
#define __fvm_node_data_conductor_h__


#include "fvm_node_data.h"




/**
 *  FVM nodal data for conductor region
 */
class FVM_Conductor_NodeData : public FVM_NodeData
{

  public:

    /**
     * the independent variable for conductor region
     */
    enum   ConductorData
    {
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
      _psi_last_,

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
     * the vector auxiliary variable for conductor region
     */
    enum ConductorAuxVecData
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
     * the complex auxiliary variable for conductor region
     */
    enum ConductorAuxComplexData
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
    FVM_Conductor_NodeData ( DataStorage * data_storage, const std::map<std::string, SimulationVariable> & variables )
    :FVM_NodeData ( data_storage, variables )
    {}

    /**
     * destructor
     */
    virtual ~FVM_Conductor_NodeData() {}

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
    { return FVM_NodeData::ConductorData; }

  public:

    /**
     * @return data by enum name
     */
    virtual Real  get_variable_real ( SolutionVariable variable ) const
    {
      switch ( variable )
      {
        case POTENTIAL   :  return  psi();                            /* potential */
        case ELECTRON    :  return  0.0;                              /* electron concentration */
        case HOLE        :  return  0.0;                              /* hole concentration */
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
    virtual void set_variable_real ( SolutionVariable variable, Real value )
    {
      switch ( variable )
      {
        case POTENTIAL   :  psi() = value;                             /* potential */
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
     * @return the statistic potential
     */
    virtual Real         psi()        const
    { return _data_storage->scalar ( _psi_, _offset ); }

    /**
     * @return the statistic potential
     */
    virtual Real &       psi()
    { return _data_storage->scalar ( _psi_, _offset ); }



    /**
     * @return the lattice temperature
     */
    virtual Real         T()          const
    { return _data_storage->scalar ( _T_, _offset ); }

    /**
     * @return the writable reference to lattice temperature
     */
    virtual Real &       T()
    { return _data_storage->scalar ( _T_, _offset ); }


    /**
     * @return the electron temperature, the same as lattice temperature
     */
    virtual Real         Tn()          const
    { return _data_storage->scalar ( _T_, _offset ); }


    /**
     * @return the writable reference to electron temperature, the same as lattice temperature
     */
    virtual Real &       Tn()
    { return _data_storage->scalar ( _T_, _offset ); }


    /**
     * @return the hole temperature, the same as lattice temperature
     */
    virtual Real         Tp()          const
    { return _data_storage->scalar ( _T_, _offset ); }


    /**
     * @return the writable reference to hole temperature, the same as lattice temperature
     */
    virtual Real &       Tp()
    { return _data_storage->scalar ( _T_, _offset ); }


    /**
     * @return the statistic potential
     */
    virtual std::complex<Real>         psi_ac()          const
    { return _data_storage->complex ( _psi_ac_, _offset ); }

    /**
     * @return the writable reference to statistic potential
     */
    virtual std::complex<Real> &       psi_ac()
    { return _data_storage->complex ( _psi_ac_, _offset ); }

    /**
     * @return the lattice temperature
     */
    virtual std::complex<Real>         T_ac()          const
    { return  _data_storage->complex ( _T_ac_, _offset ); }

    /**
     * @return the writable reference to lattice temperature
     */
    virtual std::complex<Real> &       T_ac()
    { return _data_storage->complex ( _T_ac_, _offset ); }


    /**
     * @return the complex E file. only used by EM FEM solver
     */
    virtual std::complex<Real>         OptE_complex()          const
    { return _data_storage->complex ( _OpE_complex_, _offset ); }

    /**
     * @return the writable reference to complex E file. only used by EM FEM solver
     */
    virtual std::complex<Real> &       OptE_complex()
    { return _data_storage->complex ( _OpE_complex_, _offset ); }

    /**
     * @return the complex H file. only used by EM FEM solver
     */
    virtual std::complex<Real>         OptH_complex()          const
    { return _data_storage->complex ( _OpH_complex_, _offset ); }

    /**
     * @return the writable reference to complex H file. only used by EM FEM solver
     */
    virtual std::complex<Real> &       OptH_complex()
    { return _data_storage->complex ( _OpH_complex_, _offset ); }


    /**
     * @return the statistic potential at previous time step
     */
    virtual Real         psi_last()          const
    { return _data_storage->scalar ( _psi_last_, _offset ); }

    /**
     * @return the writable reference to statistic potential at previous time step
     */
    virtual Real &       psi_last()
    { return _data_storage->scalar ( _psi_last_, _offset ); }



    /**
     * @return the lattice temperature at previous time step
     */
    virtual Real         T_last()          const
    { return _data_storage->scalar ( _T_last_, _offset ); }

    /**
     * @return the writable reference to lattice temperature at previous time step
     */
    virtual Real &       T_last()
    { return _data_storage->scalar ( _T_last_, _offset ); }




    /**
     * @return the electron affinity
     */
    virtual Real         affinity()          const
    { return _data_storage->scalar ( _affinity_, _offset ); }

    /**
     * @return the writable reference to the electron affinity
     */
    virtual Real &       affinity()
    { return _data_storage->scalar ( _affinity_, _offset ); }


    /**
     * @return the conduction band
     */
    virtual Real         Ec()          const;

    /**
     * @return the writable reference to the conduction band
     */
    virtual Real &       Ec()
    { return _data_storage->scalar ( _Ec_, _offset ); }

    /**
     * @return the valance band
     */
    virtual Real         Ev()          const;

    /**
     * @return the writable reference to the valance band
     */
    virtual Real &       Ev()
    { return _data_storage->scalar ( _Ev_, _offset ); }

    /**
     * @return the mass density of the material
     */
    virtual Real         density()          const
    { return _data_storage->scalar ( _density_, _offset ); }

    /**
     * @return the writable reference to the mass density of the material
     */
    virtual Real &       density()
    { return _data_storage->scalar ( _density_, _offset ); }



    /**
     * @return the dielectric permittivity
     */
    virtual Real         eps()          const
    { return _data_storage->scalar ( _eps_, _offset ); }

    /**
     * @return the writable reference to the dielectric permittivity
     */
    virtual Real &       eps()
    { return _data_storage->scalar ( _eps_, _offset ); }



    /**
     * @return the megnetic permeability
     */
    virtual Real         mu()          const
    { return _data_storage->scalar ( _mu_, _offset ); }

    /**
     * @return the writable reference to the megnetic permeability
     */
    virtual Real &       mu()
    { return _data_storage->scalar ( _mu_, _offset ); }


    /**
     * @return the quasi-fermi potential of electron
     */
    virtual Real         qFn()           const;

    /**
     * @return the quasi-fermi potential of hole
     */
    virtual Real         qFp()           const;


    /**
     * @return the electrical field
     */
    virtual VectorValue<Real> E()       const
    { return _data_storage->vector ( _E_, _offset );}


    /**
     * @return the writable reference to electrical field
     */
    virtual VectorValue<Real> & E()
    { return _data_storage->vector ( _E_, _offset );}


};





#endif

