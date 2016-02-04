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

#ifndef __GENIUS_PMI_interface_h__
#define __GENIUS_PMI_interface_h__

#include <string>
#include <complex>
#include <vector>
#include <map>


#include "parser_parameter.h"   // for parameter calibrating from user input file
#include "adolc.h" // for automatic differentiation
#include "atom.h"
#include "variable_define.h"
#include "std_ext.h"
#include "vector_value.h"
#include "tensor_value.h"


using namespace adtl;

//predefine
class Point;
class FVM_NodeData;


//enable calibrate
#define __CALIBRATE__

// This is the head file of physical model interface (PMI), which contains base class of:
//   physical model interface of semiconductor   PMIS
//   physical model interface of insulator       PMII
//   physical model interface of conductor       PMIC
//   physical model interface of vacuum          PMIV
//   physical model interface of PML             PMIP
// It links the main solver and material. The solver load required parameters from
// re-implemented virtual functions.

/**
 * PMI_Environment, this structure will be passed to PMI class when initializing.
 * It contains interface information for linking main genius code to each PMI class
 */
struct PMI_Environment
{
  /**
   * the location of current Point, use "pointer to pointer" method
   * here pp_point point to p_point in the material class which point to
   * current Point as a buffer.
   */
  const Point         **    pp_point;

  /**
   * the location of current node_data, use "pointer to pointer" method
   * here pp_node_data point to pnode_data in the material class which point to
   * current FVM_NodeData as a buffer.
   */
  const FVM_NodeData **    pp_node_data;

  /**
   * the pointer to current time
   */
  const PetscScalar  *     p_clock;

  /**
   * const pointer to region variables
   */
  const std::map<std::string, SimulationVariable>  ** pp_variables;

  /**
   *  the basic length unit
   */
  double   m;

  /**
   * the basic time unit
   */
  double   s;

  /**
   * potential unit
   */
  double   V;

  /**
   * the charge unit
   */
  double   C;

  /**
   * the temperature unit
   */
  double   K;

  /**
   * constructor
   */
  PMI_Environment(const Point** point, const FVM_NodeData **node_data, const PetscScalar *time,
                  const std::map<std::string, SimulationVariable> ** variables,
                  double _m_, double _s_, double _V_, double _C_, double _K_)
  : pp_point(point), pp_node_data(node_data), p_clock(time), pp_variables(variables), m(_m_), s(_s_), V(_V_), C(_C_), K(_K_)
  {}

  /**
   * constructor
   */
  PMI_Environment(double _m_, double _s_, double _V_, double _C_, double _K_)
  : pp_point(0), pp_node_data(0), p_clock(0), pp_variables(0), m(_m_), s(_s_), V(_V_), C(_C_), K(_K_)
  {}

};

/**
 * the parameter structure
 */
struct  PARA
{
  /**
   * the name of the parameter
   */
  std::string   name;

  /**
   * an introduction of the parameter
   */
  std::string   brief_intro;

  /**
   * type of the parameter
   */
  enum ParaType {Real, String} type;

  /**
   * the physical unit of this parameter in string
   */
  std::string   unit_in_string;

  /**
   * the (scaled) physical unit. it will be multiplied to original value
   */
  PetscScalar   unit_in_real;

  /**
   * the pointer to parameter value
   */
  void * value;

  /**
   * empty constructor
   */
  PARA()
  {}

  /**
   * constructor
   */
  PARA(const std::string & _name, const std::string & _intro,
       const std::string & _unit_s, PetscScalar _unit_r, PetscScalar * _value)
      :name(_name),  brief_intro(_intro), type(Real), unit_in_string(_unit_s), unit_in_real(_unit_r), value(_value){}

  PARA(const std::string & _name, const std::string & _intro, std::string * _value)
      :name(_name),  brief_intro(_intro), type(String), unit_in_string(""), unit_in_real(0.0), value(_value){}

};


/**
 * structure for complex material refraction to a specific length of optical wave
 */
struct RefractionItem
{
  /**
   * empty constructor
   */
  RefractionItem() {}

  /**
   * constructor
   */
  RefractionItem(PetscScalar lambda, PetscScalar re, PetscScalar im)
  :wavelength(lambda), RefractionIndexRe(re), RefractionIndexIm(im) {}

  /**
   * wave length
   */
  PetscScalar  wavelength;

  /**
   * real part of refraction index
   */
  PetscScalar  RefractionIndexRe;

  /**
   * image part of refraction index
   */
  PetscScalar  RefractionIndexIm;
};


class RefractionSplineInterp
{
public:

  RefractionSplineInterp() {}

  void add_nk_sorted(PetscScalar wavelength, PetscScalar n, PetscScalar k)
  {
    _w.push_back(wavelength);
    _n.push_back(n);
    _k.push_back(k);
  }

  std::complex<PetscScalar> nk(PetscScalar wavelength) const
  {
    PetscScalar n = std::max(0.0, _eval(wavelength, _w, _n, _npp));
    PetscScalar k = std::max(0.0, _eval(wavelength, _w, _k, _kpp));
    return std::complex<PetscScalar> (n, k);
  }

  void build()
  {
    _build(_w, _n, _npp);
    _build(_w, _k, _kpp);
  }

private:

  void _build(const std::vector<PetscScalar> &t, const std::vector<PetscScalar> &y, std::vector<PetscScalar> &ypp);
  PetscScalar _eval(PetscScalar tval, const std::vector<PetscScalar> &t,  const std::vector<PetscScalar> &y, const std::vector<PetscScalar> &ypp) const;

  // data for interpolation
  std::vector<PetscScalar> _w;
  std::vector<PetscScalar> _n;
  std::vector<PetscScalar> _k;

  // precomputed value
  std::vector<PetscScalar> _npp;
  std::vector<PetscScalar> _kpp;

};


/**
 * PMI_Server, the base class of PMI
 */
class PMI_Server
{
protected:
  //--------------------------------------------------------------------
  // dimension system for material database, it should keep the same as main program

  PetscScalar cm, s, V, C, K;                  // basic unit

  PetscScalar m, um, nm, J, W, kg, g, eV, ps, A, mA, Ohm;    // derived unit

  PetscScalar kb, e, me, eps0, mu0, h, hbar;     // physical constant

  PetscScalar pi;                                // math constant

protected:

  /**
   * const pointer to region variables
   */
  const std::map<std::string, SimulationVariable>  ** pp_variables;

  /**
   * the location of current point, use "pointer to pointer" method
   * here pp_point point to p_point in the material class which point to
   * current Point as a buffer.
   */
  const Point            **pp_point;

  /**
   * the location of current node_data, use "pointer to pointer" method
   * here pp_node_data point to pnode_data in the material class which point to
   * current FVM_NodeData as a buffer.
   */
  const FVM_NodeData    **pp_node_data;

  /**
   * the pointer to current time
   */
  const PetscScalar     *p_clock;

protected:
  /**
   * this map links variable \p name to its \p address
   * user can change the default value by providing variable \p name
   */
  std::map<std::string, PARA >  parameter_map;

  typedef std::pair<const std::string, PARA >  para_item;

  /**
   * an optional string provides information of the PMI library
   */
  std::string PMI_Info;

  std::string _param_string;

  std::string _calibrate_error_info;

public:
  /**
   * aux function return node coordinate.
   */
  void   ReadCoordinate (PetscScalar& x, PetscScalar& y, PetscScalar& z) const;

  /**
   * aux function return current time.
   */
  PetscScalar ReadTime () const;

  /**
   * check iff given variable eixst
   */
  bool HasVariable(const std::string &, DataType t=SCALAR) const;

  /**
   * @return index of given variable
   */
  unsigned int VariableIndex(const std::string &) const;

  /**
   * aux function return scalar value of given variable.
   */
  PetscScalar ReadRealVariable (const unsigned int) const;

  /**
   * aux function return scalar value of given variable.
   */
  PetscScalar ReadRealVariable (const std::string &) const;

  /**
   * initialize node_data and node-specific PMI data
   */
  virtual void init_node() {}

  /**
   * initialize node_data and node-specific PMI data for boundary node
   * @param bc_label  node with which label should be initialized
   */
  virtual void init_bc_node(const std::string & ) {}

  /**
   * set numeric parameter value by its name.
   * @return 0 for success, 1 for variable no find.
   */
  int calibrate_real_parameter(const std::string & var_name, PetscScalar var_value);

  /**
   * set string parameter value by its name.
   * @return 0 for success, 1 for variable no find.
   */
  int calibrate_string_parameter(const std::string & var_name, const std::string &var_value);

  /**
   * set numeric and string parameters value by its name.
   */
  virtual int calibrate(std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * post calibrate process, default do nothing
   */
  virtual void post_calibrate_process() {}

  /**
   * an interface for main code to access the parameter information in the material database
   */
  std::map<std::string, PARA > & get_parameter_info();

  /**
   * @return the string information of calibrate error
   */
  std::string & calibrate_error_info() {return _calibrate_error_info;}

  /**
   * an interface for main code to access the library information in the material database
   */
  const std::string & get_PMI_info();

  /**
   * an interface for main code to get a string representation of the
   * current parameter values in the material database
   */
  virtual const std::string& get_parameter_string(const int verbosity=0);

public:
  /**
   * constructor, link PMI object to material class
   * also set the physical constants
   */
  PMI_Server(const PMI_Environment &env);

  /**
   * destructor, seems nothing to do
   */
  virtual ~PMI_Server(){}
}
;

/*****************************************************************************
 *               Physical Model Interface for Semiconductor
 ****************************************************************************/

typedef PMI_Environment PMIS_Environment;

/**
 * PMIS_Server, the derived class of PMI_Server for semiconductor
 */
class PMIS_Server : public PMI_Server
{
private:
  // debug
  PetscScalar _Na, _Nd, _mole_x, _mole_y, _dmin;
  TensorValue<PetscScalar> _strain;
public:

  /**
   * constructor
   */
  PMIS_Server(const PMIS_Environment &env)
  : PMI_Server(env), _Na(0.0), _Nd(0.0), _mole_x(0.0), _mole_y(0.0), _dmin(0.0) {}

  /**
   * destructor
   */
  virtual ~PMIS_Server(){}

  /**
   * set the fake doping environment, debug only
   */
  void SetFakeDopingEnvironment(PetscScalar Na, PetscScalar Nd);

  /**
   * set the fake mole environment, debug only
   */
  void SetFakeMoleEnvironment(PetscScalar mole_x, PetscScalar mole_y=0.0);

  /**
   * set the fake dmin environment, debug only
   */
  void SetFakeDminEnvironment(PetscScalar dmin);


  /**
   * set the fake strain environment, debug only
   */
  void SetFakeStrainEnvironment(const TensorValue<PetscScalar> &T);

  /**
   * aux function return first mole function of current node.
   */
  PetscScalar ReadxMoleFraction () const;

  /**
   * aux function return first mole function with upper and lower bind of current node.
   */
  PetscScalar ReadxMoleFraction (const PetscScalar mole_xmin, const PetscScalar mole_xmax) const;

  /**
   * aux function return second mole function of current node.
   */
  PetscScalar ReadyMoleFraction () const;

  /**
   * aux function return second mole function with upper and lower bind of current node.
   */
  PetscScalar ReadyMoleFraction (const PetscScalar mole_ymin, const PetscScalar mole_ymax) const;

  /**
   * aux function return total Acceptor concentration of current node
   */
  PetscScalar ReadDopingNa () const;

  /**
   * aux function return total Donor concentration of current node
   */
  PetscScalar ReadDopingNd () const;

  /**
   * aux function return minimal distance to surface
   */
  PetscScalar ReadDmin () const;

  /**
   * aux function return strain tensor
   */
  TensorValue<PetscScalar> ReadStrain() const;
};


/**
 * The PMIS interface for basic physical parameters of
 * semiconductor material. User should implement each pure virtual functions.
 */
class PMIS_BasicParameter : public PMIS_Server
{
public:

  /**
   * constructor
   */
  PMIS_BasicParameter(const PMIS_Environment &env):PMIS_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIS_BasicParameter() {}

  /**
   * @return the mass density [g cm^-3] of material
   */
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;

  /**
   * @return the \p relative \p permittivity of material
   */
  virtual PetscScalar Permittivity  ()                      const  { return 1.0; }

  /**
   * @return the \p relative \p permeability of material
   */
  virtual PetscScalar Permeability  ()                      const  { return 1.0; }

  /**
   * @return strain tensor by stress tensor
   */
  virtual TensorValue<PetscScalar> Strain(const TensorValue<PetscScalar> & stress) const  { return TensorValue<PetscScalar>(); }

  /**
   * @return the affinity energy [eV] of material
   */
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;

  /**
   * get the atom fraction of this material.
   */
  virtual void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const = 0;
};



/**
 * PMIS_BandStructure. The PMIS interface for band structure of
 * semiconductor material. User should implement each pure virtual functions.
 */
class PMIS_BandStructure : public PMIS_Server
{
public:
  /**
   * constructor
   */
  PMIS_BandStructure(const PMIS_Environment &env):PMIS_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIS_BandStructure() {}

  /**
   * @return band gap of semiconductor
   */
  virtual PetscScalar Eg             (const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of band gap to lattice temperature by Automatic Differentiation
   */
  virtual AutoDScalar Eg             (const AutoDScalar &Tl) =0;

  /**
   * @return band gap narrowing due to heavy doping
   */
  virtual PetscScalar EgNarrow       (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of band gap narrowing to lattice temperature by Automatic Differentiation
   */
  virtual AutoDScalar EgNarrow       (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;

  /**
   * @return conduction band shift due to band gap narrowing
   */
  virtual PetscScalar EgNarrowToEc   (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return valence band shift due to band gap narrowing
   */
  virtual PetscScalar EgNarrowToEv   (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of conduction band shift due to band gap narrowing
   * to lattice temperature by Automatic Differentiation
   */
  virtual AutoDScalar EgNarrowToEc   (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;

  /**
   * @return partial derivatives of valence band shift due to band gap narrowing
   * to lattice temperature by Automatic Differentiation
   */
  virtual AutoDScalar EgNarrowToEv   (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;

  /**
   * @return conduction band shift due to strain
   */
  virtual PetscScalar dEcStrain   () { return 0.0; }

  /**
   * @return valence band shift due to strain
   */
  virtual PetscScalar dEvStrain   () { return 0.0; }

  /**
   * @return effective electron mass
   */
  virtual PetscScalar EffecElecMass  (const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of effective electron mass
   */
  virtual AutoDScalar EffecElecMass  (const AutoDScalar &Tl) =0;

  /**
   * @return effective hole mass
   */
  virtual PetscScalar EffecHoleMass  (const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of effective hole mass
   */
  virtual AutoDScalar EffecHoleMass  (const AutoDScalar &Tl) =0;

  /**
   * @return effective density of states in the conduction band
   */
  virtual PetscScalar Nc             (const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of effective density of states in the conduction band
   * to lattice temperature by Automatic Differentiation
   */
  virtual AutoDScalar Nc             (const AutoDScalar &Tl) =0;

  /**
   * @return effective density of states in the valence band
   */
  virtual PetscScalar Nv             (const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of effective density of states in the valence band
   * to lattice temperature by Automatic Differentiation
   */
  virtual AutoDScalar Nv             (const AutoDScalar &Tl) =0;

  /**
   * @return intrinsic carrier concentration
   */
  virtual PetscScalar ni            ( const PetscScalar &Tl) =0;

  /**
   * @return effective intrinsic carrier concentration
   */
  virtual PetscScalar nie            (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of effective intrinsic carrier concentration
   * to lattice temperature by Automatic Differentiation
   */
  virtual AutoDScalar nie            (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;

  /**
   * @return particle energy to elec-hole pare generation rate
   * @ref Meier, Dirk. CVD diamond sensors for particle detection and tracking. Diss. CERN, 1999.
   */
  virtual PetscScalar ParticleQuantumEffect(const PetscScalar &Tl) {return 1.76*eV + 1.84*Eg(Tl);}

  /**
   * @return the ion type by given species, the return value is defined as P-type < 0 and N-type >0
   * each semiconductor material can derive this function
   */
  virtual int IonType( const std::string & )                                                   { return 0; }

  /**
   * @return concentration of Na with incomplete ionization
   */
  virtual PetscScalar Na_II          (const PetscScalar &p, const PetscScalar &Tl, bool fermi) { return ReadDopingNa (); }

  /**
   * @return concentration of Na with incomplete ionization
   */
  virtual AutoDScalar Na_II          (const AutoDScalar &p, const AutoDScalar &Tl, bool fermi) { return ReadDopingNa (); }


  /**
   * @return concentration of Nd with incomplete ionization
   */
  virtual PetscScalar Nd_II          (const PetscScalar &n, const PetscScalar &Tl, bool fermi) { return ReadDopingNd (); }

  /**
   * @return concentration of Nd with incomplete ionization
   */
  virtual AutoDScalar Nd_II          (const AutoDScalar &n, const AutoDScalar &Tl, bool fermi) { return ReadDopingNd (); }

  /**
   * @return direct Recombination rate
   */
  virtual PetscScalar CDIR           (const PetscScalar &Tl) =0;

  /**
   * @return electron lift time in SHR Recombination
   */
  virtual PetscScalar TAUN           (const PetscScalar &Tl) =0;

  /**
   * @return hole lift time in SHR Recombination
   */
  virtual PetscScalar TAUP           (const PetscScalar &Tl) =0;


  /**
   * @return electron Auger Recombination rate
   */
  virtual PetscScalar AUGERN           (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return hole Auger Recombination rate
   */
  virtual PetscScalar AUGERP           (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;


  /**
   * @return electron fit parameter of Density Gradient solver
   */
  virtual PetscScalar Gamman         () {return 1.0;}

  /**
   * @return hole fit parameter of Density Gradient solver
   */
  virtual PetscScalar Gammap         () {return 1.0;}


  /**
   * @return direct recombination rate
   */
  virtual PetscScalar R_Direct     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of direct recombination rate
   */
  virtual AutoDScalar R_Direct     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;

  /**
   * @return Auger recombination rate
   */
  virtual PetscScalar R_Auger      (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;
  virtual PetscScalar R_Auger_N    (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;
  virtual PetscScalar R_Auger_P    (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of Auger recombination rate
   */
  virtual AutoDScalar R_Auger      (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;
  virtual AutoDScalar R_Auger_N    (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;
  virtual AutoDScalar R_Auger_P    (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;

  /**
   * @return SHR recombination rate
   */
  virtual PetscScalar R_SHR        (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of SHR recombination rate
   */
  virtual AutoDScalar R_SHR        (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;

  /**
   * @return SHR recombination rate at surface
   */
  virtual PetscScalar R_Surf       (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of SHR recombination rate at surface
   */
  virtual AutoDScalar R_Surf       (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;

  /**
   * @return bulk recombination rate
   */
  virtual PetscScalar Recomb       (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of bulk recombination rate
   */
  virtual AutoDScalar Recomb       (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;



  /**
   * @return electron energy relaxation time for EBM simulation
   */
  virtual PetscScalar ElecEnergyRelaxTime(const PetscScalar &Tn, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of electron energy relaxation time for EBM simulation
   */
  virtual AutoDScalar ElecEnergyRelaxTime(const AutoDScalar &Tn, const AutoDScalar &Tl) =0;

  /**
   * @return hole energy relaxation time for EBM simulation
   */
  virtual PetscScalar HoleEnergyRelaxTime(const PetscScalar &Tp, const PetscScalar &Tl) =0;

  /**
   * @return partial derivatives of hole energy relaxation time for EBM simulation
   */
  virtual AutoDScalar HoleEnergyRelaxTime(const AutoDScalar &Tp, const AutoDScalar &Tl) =0;



  /**
   * @return electron current at Schottky contact
   */
  virtual PetscScalar SchottyJsn (PetscScalar n,PetscScalar Tl,PetscScalar Vb)=0;

  /**
   * @return partial derivatives of electron current at Schottky contact
   */
  virtual AutoDScalar SchottyJsn (AutoDScalar n,AutoDScalar Tl,AutoDScalar Vb)=0;

  /**
   * @return hole current at Schottky contact
   */
  virtual PetscScalar SchottyJsp (PetscScalar p,PetscScalar Tl,PetscScalar Vb)=0;

  /**
   * @return partial derivatives of hole current at Schottky contact
   */
  virtual AutoDScalar SchottyJsp (AutoDScalar p,AutoDScalar Tl,AutoDScalar Vb)=0;

  /**
   * @return electron Richardson constant
   */
  virtual PetscScalar ARichN()=0;

  /**
   * @return hole Richardson constant
   */
  virtual PetscScalar ARichP()=0;

  /**
   * @return Schottky barrier lowerring due to local electrical field
   */
  virtual PetscScalar SchottyBarrierLowerring (PetscScalar eps, PetscScalar E)=0;



  /**
   * @return electron thermal emit velocity at hetero-junction
   * @note   D.Schroeder, Modelling of Interface Carrier Transport for Device Simulation, Springer, 1994.
   */
  virtual PetscScalar ThermalVn (PetscScalar Tl)=0;

  /**
   * @return partial derivatives of electron thermal emit velocity at hetero-junction
   */
  virtual AutoDScalar ThermalVn (AutoDScalar Tl)=0;

  /**
   * @return hole thermal emit velocity at hetero-junction
   * @note   D.Schroeder, Modelling of Interface Carrier Transport for Device Simulation, Springer, 1994.
   */
  virtual PetscScalar ThermalVp (PetscScalar Tl)=0;

  /**
   * @return partial derivatives of hole thermal emit velocity at hetero-junction
   */
  virtual AutoDScalar ThermalVp (AutoDScalar Tl)=0;




  /**
   * Hot Carrier Injection: probability that an electron will not be scattered in the semiconductor before reaching the interface
   * @param dis distance from the point to the interface.
   */
  virtual PetscScalar HCI_Probability_Semiconductor_n(const PetscScalar &dis) { return 0.0; }

  /**
   * Hot Carrier Injection: probability that a hole will not be scattered in the semiconductor before reaching the interface
   * @param dis distance from the point to the interface.
   */
  virtual PetscScalar HCI_Probability_Semiconductor_p(const PetscScalar &dis) { return 0.0; }

  /**
   * Hot Carrier Injection: Fiegna integral over electron energy distribution
   * @param phin semiconductor-insulator potential barrier to electron
   * @param Eeff electric field in the direction of electron current flow
   */
  virtual PetscScalar HCI_Integral_Fiegna_n(const PetscScalar &phin, const PetscScalar &Eeff) { return 0.0; }

  /**
   * Hot Carrier Injection: Fiegna integral over hole energy distribution
   * @param phip semiconductor-insulator potential barrier to hole
   * @param Eeff electric field in the direction of hole current flow
   */
  virtual PetscScalar HCI_Integral_Fiegna_p(const PetscScalar &phip, const PetscScalar &Eeff) { return 0.0; }

  /**
   * Hot Carrier Injection: Classical integral over electron energy distribution
   * @param phin semiconductor-insulator potential barrier to electron
   * @param Eeff electric field in the direction of electron current flow
   */
  virtual PetscScalar HCI_Integral_Classical_n(const PetscScalar &phin, const PetscScalar &Eeff) { return 0.0; }

  /**
   * Hot Carrier Injection: Classical integral over hole energy distribution
   * @param phip semiconductor-insulator potential barrier to hole
   * @param Eeff electric field in the direction of hole current flow
   */
  virtual PetscScalar HCI_Integral_Classical_p(const PetscScalar &phip, const PetscScalar &Eeff) { return 0.0; }


  /**
   * @return band to band tunneling rate
   */
  virtual PetscScalar BB_Tunneling(const PetscScalar &Tl, const PetscScalar &E) =0;

  /**
   * @return partial derivatives of band to band tunneling rate
   */
  virtual AutoDScalar BB_Tunneling(const AutoDScalar &Tl, const AutoDScalar &E) =0;


};



/**
 * The PMIS interface for semiconductor mobility.
 * User should implement each pure virtual functions.
 */
class PMIS_Mobility : public PMIS_Server
{
public:
  /**
   * constructor
   */
  PMIS_Mobility(const PMIS_Environment &env):PMIS_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIS_Mobility() {}

  /**
   * A factor used in determining the effective electric field at interfaces
   * used in the field-dependent mobility models for electrons.
   */
  virtual PetscScalar ZETAN() { return 0.0; }

  /**
   * A factor used in determining the effective electric field at interfaces
   * used in the field-dependent mobility models for electrons.
   */
  virtual PetscScalar ETAN()  { return 0.0; }

  /**
   * A factor used in determining the effective electric field at interfaces
   * used in field-dependent mobility models for holes.
   */
  virtual PetscScalar ZETAP() { return 0.0; }

  /**
   * A factor used in determining the effective electric field at interfaces
   * used in field-dependent mobility models for holes.
   */
  virtual PetscScalar ETAP()  { return 0.0; }

  /**
   * Hall mobility factor  for electrons
   */
  virtual PetscScalar RH_ELEC()  { return 1.0; }

  /**
   * Hall mobility factor  for holes
   */
  virtual PetscScalar RH_HOLE()  { return 1.0; }

  /**
   * @return the electron mobility
   */
  virtual PetscScalar ElecMob (const PetscScalar &p,  const PetscScalar &n,  const PetscScalar &Tl,
                               const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const=0;

  /**
   * @return the hole mobility
   */
  virtual PetscScalar HoleMob (const PetscScalar &p,  const PetscScalar &n,  const PetscScalar &Tl,
                               const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const=0;

  /**
   * @return the partial derivatives of electron mobility by Automatic Differentiation
   */
  virtual AutoDScalar ElecMob (const AutoDScalar &p,  const AutoDScalar &n,  const AutoDScalar &Tl,
                               const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const=0;

  /**
   * @return the partial derivatives of hole mobility by Automatic Differentiation
   */
  virtual AutoDScalar HoleMob (const AutoDScalar &p,  const AutoDScalar &n,  const AutoDScalar &Tl,
                               const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const=0;

};



/**
 * The PMIS interface for semiconductor impact ionization parameter.
 * User should implement each pure virtual functions.
 */
class PMIS_Avalanche : public PMIS_Server
{
public:
  /**
   * constructor
   */
  PMIS_Avalanche(const PMIS_Environment &env):PMIS_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIS_Avalanche() {}

  /**
   * @return the electron generation rate for DDM simulation
   */
  virtual PetscScalar ElecGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const=0;

  /**
   * @return the hole generation rate for DDM simulation
   */
  virtual PetscScalar HoleGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const=0;

  //Automatic Differentiation version for DDM

  /**
   * @return the partial derivatives of electron generation rate for DDM simulation by Automatic Differentiation
   */
  virtual AutoDScalar ElecGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const=0;

  /**
   * @return the partial derivatives of hole generation rate for DDM simulation by Automatic Differentiation
   */
  virtual AutoDScalar HoleGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const=0;

  /**
   * @return the electron generation rate for EBM simulation
   */
  virtual PetscScalar ElecGenRateEBM (const PetscScalar &Tn,const PetscScalar &Tl,const PetscScalar &Eg) const=0;

  /**
   * @return the hole generation rate for EBM simulation
   */
  virtual PetscScalar HoleGenRateEBM (const PetscScalar &Tp,const PetscScalar &Tl,const PetscScalar &Eg) const=0;

  //Automatic Differentiation version for EBM

  /**
   * @return the partial derivatives of electron generation rate for EBM simulation by Automatic Differentiation
   */
  virtual AutoDScalar ElecGenRateEBM (const AutoDScalar &Tn,const AutoDScalar &Tl,const AutoDScalar &Eg) const=0;

  /**
   * @return the partial derivatives of hole generation rate for EBM simulation by Automatic Differentiation
   */
  virtual AutoDScalar HoleGenRateEBM (const AutoDScalar &Tp,const AutoDScalar &Tl,const AutoDScalar &Eg) const=0;


};

class PMIS_Trap : public PMIS_Server
{
public:
  enum TrapChargeType {Neutral, Acceptor, Donor};
  enum TrapType {Bulk, Interface};

  class TrapLocation
  {
    public:
      double point[3];
      TrapType type;
      TrapLocation(double x, double y, double z, TrapType t) : type(t)
      {
        point[0] = x;
        point[1] = y;
        point[2] = z;
      }
  };

  class TrapLocationComp
  {
    public:
      bool operator() (const TrapLocation &lhs, const TrapLocation &rhs) const
      {
        if (lhs.point[0] < rhs.point[0])
          return true;
        if (lhs.point[0] > rhs.point[0])
          return false;

        if (lhs.point[1] < rhs.point[1])
          return true;
        if (lhs.point[1] > rhs.point[1])
          return false;

        if (lhs.point[2] < rhs.point[2])
          return true;
        if (lhs.point[2] > rhs.point[2])
          return false;
        return lhs.type<rhs.type;
      }
  };

public:
  /**
   * constructor
   */
  PMIS_Trap(const PMIS_Environment &env):PMIS_Server(env) { }
  /**
   * destructor
   */
  virtual ~PMIS_Trap() {}

  /**
   * returns the electric charge density due to trapped charge at this node
   * one should call Calculate() to calculate the electron occupancy before calling this function
   */
  virtual PetscScalar Charge(const bool flag_bulk) = 0;

  /**
   * returns the partial derivatives of electric charge density
   * w.r.t. the local V,n,p
   * one should call CalculateAD() to calculate the electron occupancy before calling this function
   */
  virtual AutoDScalar ChargeAD(const bool flag_bulk) = 0;

  /**
   * Calculates the electron trapping rate
   * one should call Calculate() to calculate the electron occupancy before calling this function
   */
  virtual PetscScalar ElectronTrapRate(const bool flag_bulk, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tl) = 0;

  /**
   * Calculates the hole trapping rate
   * one should call Calculate() to calculate the electron occupancy before calling this function
   */
  virtual PetscScalar HoleTrapRate(const bool flag_bulk, const PetscScalar &p, const PetscScalar &ni, const PetscScalar &Tl) = 0;

  /**
   * Calculates the partial derivatives of electron trapping rate
   * one should call CalculateAD() to calculate the electron occupancy before calling this function
   */
  virtual AutoDScalar ElectronTrapRate(const bool flag_bulk, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tl) = 0;

  /**
   * Calculates the partial derivatives of hole trapping rate
   * one should call CalculateAD() to calculate the electron occupancy before calling this function
   */
  virtual AutoDScalar HoleTrapRate(const bool flag_bulk, const AutoDScalar &p, const AutoDScalar &ni, const AutoDScalar &Tl) = 0;

  /**
   * Calculate the heat transferred to lattice in trapping
   * one should call Calculate() to calculate the electron occupancy before calling this function
   * @EcEi   Energy difference between Conduction band and Intrinsic level
   */
  virtual PetscScalar TrapHeat(const bool flag_bulk, const PetscScalar &p, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tp, const PetscScalar &Tn, const PetscScalar &Tl, const PetscScalar &EcEi, const PetscScalar &EiEv) = 0;

  /**
   * Calculate the partial derivatives of heat transferred to lattice in trapping
   * one should call CalculateAD() to calculate the electron occupancy before calling this function
   * @EcEi   Energy difference between Conduction band and Intrinsic level
   */
  virtual AutoDScalar TrapHeat(const bool flag_bulk, const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tp, const AutoDScalar &Tn, const AutoDScalar &Tl, const AutoDScalar &EcEi, const AutoDScalar &EiEv) = 0;

  /**
   * Calculate trap occupancy.
   */
  virtual void Calculate(const bool flag_bulk, const PetscScalar &p, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tl) = 0;

  /**
   * partial derivatives of trap occupancy w.r.t. local V,n,p,
   */
  virtual void Calculate(const bool flag_bulk, const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tl) = 0;

  /**
   * process the parameters of the PMI command
   */
  virtual void Update(const bool flag_bulk, const PetscScalar &p, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tl) = 0;

};


/**
 * The PMIS interface for semiconductor thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMIS_Thermal : public PMIS_Server
{
public:
  /**
   * constructor
   */
  PMIS_Thermal(const PMIS_Environment &env):PMIS_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIS_Thermal() {}

  /**
   * @return the heat capacity [J/(K*cm^3)] of the material
   */
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;

  /**
   * @return the partial derivatives of heat capacity by AD
   */
  virtual AutoDScalar HeatCapacity  (const AutoDScalar &Tl) const=0;

  /**
   * @return the heat conduction [W/cm/K] of the material
   */
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;

  /**
   * @return the partial derivatives of heat conduction by AD
   */
  virtual AutoDScalar HeatConduction(const AutoDScalar &Tl) const=0;

};



/**
 * The PMIS interface for semiconductor optical refraction index
 * User should implement each pure virtual functions.
 */
class PMIS_Optical : public PMIS_Server
{
protected:

  std::string _refraction_data_file;

  std::vector<RefractionItem> _wave_table;

  /**
   * when refraction_data_file is not empty, read from it
   */
  virtual void post_calibrate_process();

public:

  /**
   * constructor
   */
  PMIS_Optical(const PMIS_Environment &env):PMIS_Server(env)
  {
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("refraction", PARA("refraction", "The refraction data file", &_refraction_data_file)) );
#endif
  }

  /**
   * destructor
   */
  virtual ~PMIS_Optical() {}


  /**
   * an interface for main code to get a string representation of the
   * current parameter values in the material database
   */
  virtual const std::string& get_parameter_string(const int verbosity=0);

  /**
   * @return the complex refraction index of material
   */
  virtual std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Tl, PetscScalar Eg=0) const=0;

  /**
   * @return the free carrier absorption coefficient
   */
  virtual PetscScalar FreeCarrierAbsorption(PetscScalar lamda, PetscScalar n, PetscScalar p, PetscScalar Tl) const { return 0.0; }
};


/*****************************************************************************
 *               Physical Model Interface for Insulator
 ****************************************************************************/

typedef PMI_Environment PMII_Environment;

/**
 * PMII_Server, the derived class of PMI_Server for insulator
 * however, no extra server function
 */
class PMII_Server : public PMI_Server
{
public:
  /**
   * constructor
   */
  PMII_Server(const PMII_Environment &env): PMI_Server(env)
  {}

  /**
   * destructor
   */
  virtual ~PMII_Server() {}
};


/**
 * The PMII interface for basic physical parameters of
 * insulator material. User should implement each pure virtual functions.
 */
class PMII_BasicParameter : public PMII_Server
{
public:

  /**
   * constructor
   */
  PMII_BasicParameter(const PMII_Environment &env):PMII_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMII_BasicParameter() {}

  /**
   * @return the mass density [g cm^-3] of material
   */
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;

  /**
   * @return the \p relative \p permittivity of material
   */
  virtual PetscScalar Permittivity  ()                      const=0;

  /**
   * @return the \p relative \p permeability of material
   */
  virtual PetscScalar Permeability  ()                      const=0;

  /**
   * @return the affinity energy [eV] of material
   */
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;

  /**
   * @return the electrical conductance of material
   */
  virtual PetscScalar Conductance   (const PetscScalar &Tl, const PetscScalar &E)  const=0;

  /**
   * @return the electrical conductance of material
   */
  virtual AutoDScalar Conductance   (const AutoDScalar &Tl, const AutoDScalar &E)  const=0;

  /**
   * @return the radiative generation rate #/cm^-3/rad, also consider recombination yield under E field
   */
  virtual PetscScalar RadGenRate   (const PetscScalar &E) const=0;

  /**
   * @return the radiation induced conductance of material, DRate in the unit of Gy/s
   */
  virtual PetscScalar RadConductance   (const PetscScalar &DRate)  const=0;

  /**
   * @return the electron mobility of material
   */
  virtual PetscScalar ElecMobility   (const PetscScalar &Tl) const=0;

  /**
   * @return the hole mobility of material
   */
  virtual PetscScalar HoleMobility   (const PetscScalar &Tl) const=0;

  /**
   * @return the H+ mobility of material
   */
  virtual PetscScalar HIonMobility   (const PetscScalar &Tl) const { return 0.0; }

  /**
   * get the atom fraction of this material.
   */
  virtual void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const = 0;
};



/**
 * PMII_BandStructure. The PMII interface for band structure of
 * insulator material. User should implement each pure virtual functions.
 */
class PMII_BandStructure : public PMII_Server
{
  public:

  /**
   * constructor
   */
  PMII_BandStructure(const PMII_Environment &env):PMII_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMII_BandStructure() {}

  /**
   * @return band gap of semiconductor
   */
  virtual PetscScalar Eg             (const PetscScalar &Tl) const = 0;

  /**
   * @return particle energy to elec-hole pare generation rate
   */
  virtual PetscScalar ParticleQuantumEffect(const PetscScalar &Tl) const {return 3.0*Eg(Tl);}

  /**
   * @return effective electron mass
   */
  virtual PetscScalar EffecElecMass  (const PetscScalar &Tl) const=0;

  /**
   * @return effective hole mass
   */
  virtual PetscScalar EffecHoleMass  (const PetscScalar &Tl) const=0;


  /**
   * @return trapA density
   */
  virtual PetscScalar TrapADensity() const { return 0.0; }

  /**
   * @return electron trapA capture cross section
   */
  virtual PetscScalar TrapACaptureElecCS(const PetscScalar &Tl) const { return 0.0; }

  /**
   * @return hole trapA capture cross section
   */
  virtual PetscScalar TrapACaptureHoleCS(const PetscScalar &Tl) const { return 0.0; }


  /**
   * @return trapB density
   */
  virtual PetscScalar TrapBDensity() const { return 0.0; }

  /**
   * @return electron trapB capture cross section
   */
  virtual PetscScalar TrapBCaptureElecCS(const PetscScalar &Tl) const { return 0.0; }

  /**
   * @return hole trapB capture cross section
   */
  virtual PetscScalar TrapBCaptureHoleCS(const PetscScalar &Tl) const { return 0.0; }

  /**
   * @return H+ trapB release rate
   */
  virtual PetscScalar TrapBReleaseHIonRate(const PetscScalar &Tl) const { return 0.0; }


  /**
   * @return interface state density
   */
  virtual PetscScalar InterfaceStateDensity() const { return 0.0; }

  /**
   * @return H+ interface reaction cross section
   */
  virtual PetscScalar HIonInterfaceTrapCS(const PetscScalar &Tl) const { return 0.0; }

  /**
   * @return electron fit parameter of Density Gradient solver
   */
  virtual PetscScalar Gamman         () const {return 1.0;}

  /**
   * @return hole fit parameter of Density Gradient solver
   */
  virtual PetscScalar Gammap         () const {return 1.0;}


  /**
   * @return electron Richardson constant
   */
  virtual PetscScalar ARichardson() const = 0;


  /**
   * @return electron Richardson constant
   */
  virtual PetscScalar ElecInject(const PetscScalar &W, const PetscScalar &E, const PetscScalar &Tl) const {return 0.0;}

  /**
   * @return electron Richardson constant
   */
  virtual AutoDScalar ElecInject(const AutoDScalar &W, const AutoDScalar &E, const AutoDScalar &Tl) const {return 0.0;}


  /**
   * Hot Carrier Injection: effective semiconductor-insulator interface barrier to electron
   */
  virtual PetscScalar HCI_Barrier_n(const PetscScalar &affinity_semi, const PetscScalar & Eg_semi, const PetscScalar &affinity_ins,
                                    const PetscScalar &t_ins, const PetscScalar &E_ins) const = 0;

  /**
   * Hot Carrier Injection: effective semiconductor-insulator interface barrier to hole
   */
  virtual PetscScalar HCI_Barrier_p(const PetscScalar &affinity_semi, const PetscScalar & Eg_semi, const PetscScalar &affinity_ins,
                                    const PetscScalar &t_ins, const PetscScalar &E_ins) const = 0;

  /**
   * Hot Carrier Injection: probability that an electron will not be scattered in the insulator
   */
  virtual PetscScalar HCI_Probability_Insulator_n(const PetscScalar &t_ins, const PetscScalar &E_ins) const = 0;

  /**
   * Hot Carrier Injection: probability that a hole will not be scattered in the insulator
   */
  virtual PetscScalar HCI_Probability_Insulator_p(const PetscScalar &t_ins, const PetscScalar &E_ins) const = 0;

  /**
   * Fowler-Nordheim tunneling
   */
  virtual PetscScalar J_FN_Tunneling(const PetscScalar &E_ins, const PetscScalar &alpha) const = 0;

  /**
   * Conduction band electron tunneling
   */
  virtual PetscScalar J_CBET_Tunneling(const PetscScalar &m, const PetscScalar &Tl,
                                       const PetscScalar &Efn1, const PetscScalar &Efn2,
                                       const PetscScalar &Ec1,  const PetscScalar &Ec2,
                                       const PetscScalar &B1,   const PetscScalar &B2,
                                       const PetscScalar &t) const = 0;

  virtual AutoDScalar J_CBET_Tunneling(const PetscScalar &m, const AutoDScalar &Tl,
                                       const AutoDScalar &Efn1, const AutoDScalar &Efn2,
                                       const AutoDScalar &Ec1,  const AutoDScalar &Ec2,
                                       const AutoDScalar &B1,   const AutoDScalar &B2,
                                       const PetscScalar &t) const = 0;


  /**
   * Valence band hole tunneling
   */
  virtual PetscScalar J_VBHT_Tunneling(const PetscScalar &m, const PetscScalar &Tl,
                                       const PetscScalar &Efn1, const PetscScalar &Efn2,
                                       const PetscScalar &Ec1,  const PetscScalar &Ec2,
                                       const PetscScalar &B1,   const PetscScalar &B2,
                                       const PetscScalar &t) const = 0;

  virtual AutoDScalar J_VBHT_Tunneling(const PetscScalar &m, const AutoDScalar &Tl,
                                       const AutoDScalar &Efn1, const AutoDScalar &Efn2,
                                       const AutoDScalar &Ec1,  const AutoDScalar &Ec2,
                                       const AutoDScalar &B1,   const AutoDScalar &B2,
                                       const PetscScalar &t) const = 0;

  /**
   * Valence band electron tunneling
   */
  virtual PetscScalar J_VBET_Tunneling(const PetscScalar &m, const PetscScalar &Tl,
                                       const PetscScalar &Efn1, const PetscScalar &Efn2,
                                       const PetscScalar &Ec1,  const PetscScalar &Ec2,
                                       const PetscScalar &Ev1,  const PetscScalar &Ev2,
                                       const PetscScalar &B1,   const PetscScalar &B2,
                                       const PetscScalar &t) const = 0;

  virtual AutoDScalar J_VBET_Tunneling(const PetscScalar &m, const AutoDScalar &Tl,
                                       const AutoDScalar &Efn1, const AutoDScalar &Efn2,
                                       const AutoDScalar &Ec1,  const AutoDScalar &Ec2,
                                       const AutoDScalar &Ev1,  const AutoDScalar &Ev2,
                                       const AutoDScalar &B1,   const AutoDScalar &B2,
                                       const PetscScalar &t) const = 0;

};


/**
 * The PMII interface for insulator thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMII_Thermal : public PMII_Server
{
public:
  /**
   * constructor
   */
  PMII_Thermal(const PMII_Environment &env):PMII_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMII_Thermal() {}

  /**
   * @return the heat capacity [J/(K*cm^3)] of the material
   */
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;

  /**
   * @return the heat conduction [W/cm/K] of the material
   */
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;

};


/**
 * The PMII interface for insulator optical refraction index
 * User should implement each pure virtual functions.
 */
class PMII_Optical : public PMII_Server
{
protected:

  std::string _refraction_data_file;

  std::vector<RefractionItem> _wave_table;

  /**
   * when refraction_data_file is not empty, read from it
   */
  virtual void post_calibrate_process();

public:
  /**
   * constructor
   */
  PMII_Optical(const PMII_Environment &env):PMII_Server(env)
  {
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("refraction", PARA("refraction", "The refraction data file", &_refraction_data_file)) );
#endif
  }

  /**
   * destructor
   */
  virtual ~PMII_Optical() {}

  /**
   * @return the complex refraction index of material
   */
  virtual std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Tl, PetscScalar Eg=0) const=0;
};


/*****************************************************************************
 *               Physical Model Interface for Conductor
 ****************************************************************************/

typedef PMI_Environment PMIC_Environment;


/**
 * PMIC_Server, the derived class of PMI_Server for conductor
 * however, no extra server function
 */
class PMIC_Server : public PMI_Server
{
public:
  /**
   * constructor
   */
  PMIC_Server(const PMIC_Environment &env): PMI_Server(env)
  {}

  /**
   * destructor
   */
  virtual ~PMIC_Server(){}
};


/**
 * The PMIC interface for basic physical parameters of
 * conductor material. User should implement each pure virtual functions.
 */
class PMIC_BasicParameter : public PMIC_Server
{
public:
  /**
   * constructor
   */
  PMIC_BasicParameter(const PMIC_Environment &env):PMIC_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIC_BasicParameter() {}

  /**
   * @return the mass density [g cm^-3] of material
   */
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;

  /**
   * @return the ion density [cm^-3] of material
   */
  virtual PetscScalar IonDensity    (const PetscScalar &Tl) const=0;

  /**
   * @return the \p relative \p permittivity of material
   */
  virtual PetscScalar Permittivity  ()                      const=0;

  /**
   * @return the \p relative \p permeability of material
   */
  virtual PetscScalar Permeability  ()                      const=0;

  /**
   * @return the affinity energy [eV] of material
   */
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;

  /**
   * @return the electrical conductance of material
   */
  virtual PetscScalar Conductance   ()                      const=0;
   
  /**
   * @return the current density under given E and Tl
   */ 
  virtual PetscScalar CurrentDensity(const PetscScalar &E, const PetscScalar &Tl) const=0; 
     
  /**
   * @return the current density under given E and Tl
   */ 
  virtual AutoDScalar CurrentDensity(const AutoDScalar &E, const AutoDScalar &Tl) const=0; 
  
  /**
   * @return the thermal emit velocity of material
   */
  virtual PetscScalar ThermalVn   (const PetscScalar &Tl)   const=0;

  /**
   * get the atom fraction of this material.
   */
  virtual void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const = 0;

};


/**
 * The PMIC interface for conductor thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMIC_Thermal : public PMIC_Server
{
public:
  /**
   * constructor
   */
  PMIC_Thermal(const PMIC_Environment &env):PMIC_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIC_Thermal() {}

  /**
   * @return the heat capacity [J/(K*cm^3)] of the material
   */
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;

  /**
   * @return the heat conduction [W/cm/K] of the material
   */
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;

};


/**
 * The PMIC interface for conductor optical refraction index
 * User should implement each pure virtual functions.
 */
class PMIC_Optical : public PMIC_Server
{
protected:

  std::string _refraction_data_file;

  std::vector<RefractionItem> _wave_table;

  /**
   * when refraction_data_file is not empty, read from it
   */
  virtual void post_calibrate_process();

public:
  /**
   * constructor
   */
  PMIC_Optical(const PMIC_Environment &env):PMIC_Server(env)
  {
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("refraction", PARA("refraction", "The refraction data file", &_refraction_data_file)) );
#endif
  }

  /**
   * destructor
   */
  virtual ~PMIC_Optical() {}

  /**
   * @return the complex refraction index of material
   */
  virtual std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Tl, PetscScalar Eg=0) const=0;
};






/*****************************************************************************
 * Physical Model Interface for Vacuum, only used in EM FEM solver.
 ****************************************************************************/

typedef PMI_Environment PMIV_Environment;

/**
 * PMIC_Server, the derived class of PMI_Server for Vacuum
 * however, no extra server function
 */
class PMIV_Server : public PMI_Server
{
public:
  /**
   * constructor
   */
  PMIV_Server(const PMIV_Environment &env): PMI_Server(env)
  {}

  /**
   * destructor
   */
  virtual ~PMIV_Server(){}
};


/**
 * The PMIV interface for basic physical parameters of
 * vacuum. User should implement each pure virtual functions.
 */
class PMIV_BasicParameter : public PMIV_Server
{
public:
  /**
   * constructor
   */
  PMIV_BasicParameter(const PMIV_Environment &env):PMIV_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIV_BasicParameter() {}

  /**
   * @return the mass density [g cm^-3] of material
   */
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;

  /**
   * @return the \p relative \p permittivity of material
   */
  virtual PetscScalar Permittivity  ()                      const=0;

  /**
   * @return the \p relative \p permeability of material
   */
  virtual PetscScalar Permeability  ()                      const=0;

  /**
   * @return the affinity energy [eV] of material
   */
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;

  /**
   * get the atom fraction of this material.
   */
  virtual void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const = 0;

};


/**
 * The PMIV interface for vacuum thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMIV_Thermal : public PMIV_Server
{
public:
  /**
   * constructor
   */
  PMIV_Thermal(const PMIV_Environment &env):PMIV_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIV_Thermal() {}

  /**
   * @return the heat capacity [J/(K*cm^3)] of the material
   */
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;

  /**
   * @return the heat conduction [W/cm/K] of the material
   */
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;

};



/**
 * The PMIV interface for conductor optical refraction index
 * User should implement each pure virtual functions.
 */
class PMIV_Optical : public PMIV_Server
{
protected:

  std::string _refraction_data_file;

  std::vector<RefractionItem> _wave_table;

public:
  /**
   * constructor
   */
  PMIV_Optical(const PMIC_Environment &env):PMIV_Server(env)
  {
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("refraction", PARA("refraction", "The refraction data file", &_refraction_data_file)) );
#endif
  }

  /**
   * destructor
   */
  virtual ~PMIV_Optical() {}

  /**
   * @return the complex refraction index of material
   */
  virtual std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Tl, PetscScalar Eg=0) const=0;
};



/*****************************************************************************
 * Physical Model Interface for PML Boundary, only used in EM FEM solver.
 ****************************************************************************/

typedef PMI_Environment PMIP_Environment;

/**
 * PMIP_Server, the derived class of PMI_Server for PML material
 * however, no extra server function
 */
class PMIP_Server : public PMI_Server
{
public:
  /**
   * constructor
   */
  PMIP_Server(const PMIP_Environment &env) : PMI_Server(env)
  {}

  /**
   * destructor
   */
  virtual ~PMIP_Server(){}
};


/**
 * The PMIP interface for basic physical parameters of
 * PML. User should implement each pure virtual functions.
 */
class PMIP_BasicParameter : public PMIP_Server
{
public:
  /**
   * constructor
   */
  PMIP_BasicParameter(const PMIP_Environment &env):PMIP_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIP_BasicParameter() {}

  /**
   * @return the mass density [g cm^-3] of material
   */
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;

  /**
   * @return the \p relative \p permittivity of material
   */
  virtual PetscScalar Permittivity  ()                      const=0;

  /**
   * @return the \p relative \p permeability of material
   */
  virtual PetscScalar Permeability  ()                      const=0;

  /**
   * @return the affinity energy [eV] of material
   */
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;

  /**
   * get the atom fraction of this material.
   */
  virtual void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const = 0;

};


/**
 * The PMIP interface for PML thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMIP_Thermal : public PMIP_Server
{
public:
  /**
   * constructor
   */
  PMIP_Thermal(const PMIP_Environment &env):PMIP_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIP_Thermal() {}

  /**
   * @return the heat capacity [J/(K*cm^3)] of the material
   */
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;

  /**
   * @return the heat conduction [W/cm/K] of the material
   */
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;

};

#endif

