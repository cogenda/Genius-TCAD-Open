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

#include "point.h"
#include "key.h"   // for parameter calibrating from user input file
#include "adolc.h" // for automatic differentiation
using namespace adtl;

//predefine
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
                  double _m_, double _s_, double _V_, double _C_, double _K_)
      :pp_point(point), pp_node_data(node_data), p_clock(time), m(_m_), s(_s_), V(_V_), C(_C_), K(_K_)
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
  PetscScalar  wavelength;
  PetscScalar  RefractionIndexRe;
  PetscScalar  RefractionIndexIm;
};


/**
 * PMI_Server, the base class of PMI
 */
class PMI_Server
{
protected:
  //--------------------------------------------------------------------
  // dimension system for material database, it should keep the same as main program

  PetscScalar cm,s,V,C,K;                  // basic unit

  PetscScalar m,um,J,W,kg,g,eV,ps,A,mA;    // derived unit

  PetscScalar kb,e,me,eps0,mu0,h,hbar;     // hpysical constant

protected:
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
  std::map<const std::string, PARA >  parameter_map;

  typedef std::pair<const std::string, PARA >  para_item;

  /**
   * an optional string provides information of the PMI library
   */
  std::string PMI_Info;

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
   * aux function to return user defined scalar value
   */
  PetscScalar ReadUserScalarValue (std::string &name) const;

  /**
   * initialize node_data and node-specific PMI data
   */
  virtual void init_node() {};

  /**
   * initialize node_data and node-specific PMI data for boundary node
   * @param bc_label  node with which label should be initialized
   */
  virtual void init_bc_node(const std::string & ) {};

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
  virtual int calibrate(const std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * an interface for main code to access the parameter information in the material database
   */
  std::map<const std::string, PARA > & get_parameter_info();

  /**
   * an interface for main code to access the library information in the material database
   */
  const std::string & get_PMI_info();

  /**
   * an interface for main code to get a string representation of the
   * current parameter values in the material database
   */
  virtual std::string get_parameter_string(const int verbosity=0);

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
public:

  /**
   * constructor
   */
  PMIS_Server(const PMIS_Environment &env) : PMI_Server(env) {}

  /**
   * destructor
   */
  virtual ~PMIS_Server(){}

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
  virtual ~PMIS_BasicParameter() {};

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
  virtual void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction) const = 0;
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
  virtual ~PMIS_BandStructure() {};

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
   * @return electron lift time in SHR Recombination
   */
  virtual PetscScalar TAUN           (const PetscScalar &Tl) =0;

  /**
   * @return hole lift time in SHR Recombination
   */
  virtual PetscScalar TAUP           (const PetscScalar &Tl) =0;

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
   * @return partial derivatives of electron current at Schottky contact to electron density
   * may be replaced by AD later
   */
  virtual PetscScalar pdSchottyJsn_pdn(PetscScalar n,PetscScalar Tl,PetscScalar Vb)=0;

  /**
   * @return partial derivatives of hole current at Schottky contact to hole density
   * may be replaced by AD later
   */
  virtual PetscScalar pdSchottyJsp_pdp(PetscScalar p,PetscScalar Tl,PetscScalar Vb)=0;

  /**
   * @return partial derivatives of electron current at Schottky contact to lattice temperature
   * may be replaced by AD later
   */
  virtual PetscScalar pdSchottyJsn_pdTl(PetscScalar n,PetscScalar Tl,PetscScalar Vb)=0;

  /**
   * @return partial derivatives of hole current at Schottky contact to lattice temperature
   * may be replaced by AD later
   */
  virtual PetscScalar pdSchottyJsp_pdTl(PetscScalar p,PetscScalar Tl,PetscScalar Vb)=0;



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
   * @return partial derivatives of electron thermal emit velocity at hetero-junction to lattice temperature
   * may be replaced by AD later
   */
  virtual PetscScalar pdThermalVn_pdTl (PetscScalar Tl)=0;

  /**
   * @return partial derivatives of hole thermal emit velocity at hetero-junction to lattice temperature
   * may be replaced by AD later
   */
  virtual PetscScalar pdThermalVp_pdTl (PetscScalar Tl)=0;



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
  virtual ~PMIS_Mobility() {};

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
  virtual ~PMIS_Avalanche() {};

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
      Point point;
      TrapType type;
      TrapLocation(const Point &p, TrapType t) : point(p), type(t) {};
  };

  class TrapLocationComp
  {
    public:
      bool operator() (const TrapLocation &lhs, const TrapLocation &rhs) const
      {
        if (lhs.point(0) < rhs.point(0))
          return true;
        if (lhs.point(0) > rhs.point(0))
          return false;

        if (lhs.point(1) < rhs.point(1))
          return true;
        if (lhs.point(1) > rhs.point(1))
          return false;

        if (lhs.point(2) < rhs.point(2))
          return true;
        if (lhs.point(2) > rhs.point(2))
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
  virtual ~PMIS_Trap() {};

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
  virtual ~PMIS_Thermal() {};

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
public:

  /**
   * constructor
   */
  PMIS_Optical(const PMIS_Environment &env):PMIS_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIS_Optical() {}

  /**
   * @return the complex refraction index of material
   */
  virtual std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0, PetscScalar Tl=1) const=0;
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
  virtual ~PMII_Server(){};
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
  virtual ~PMII_BasicParameter() {};

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
   * @return band gap of insulator
   */
  virtual PetscScalar Eg             (const PetscScalar &Tl) const=0;

  /**
   * get the atom fraction of this material.
   */
  virtual void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction) const = 0;
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
  virtual ~PMII_Thermal() {};

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
public:
  /**
   * constructor
   */
  PMII_Optical(const PMII_Environment &env):PMII_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMII_Optical() {}

  /**
   * @return the complex refraction index of material
   */
  virtual std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0, PetscScalar Tl=1) const=0;
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
  virtual ~PMIC_Server(){};
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
  virtual ~PMIC_BasicParameter() {};

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
  virtual PetscScalar Conductance   ()                      const=0;

  /**
   * get the atom fraction of this material.
   */
  virtual void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction) const = 0;

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
  virtual ~PMIC_Thermal() {};

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
public:
  /**
   * constructor
   */
  PMIC_Optical(const PMIC_Environment &env):PMIC_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIC_Optical() {}

  /**
   * @return the complex refraction index of material
   */
  virtual std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0, PetscScalar Tl=1) const=0;
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
  virtual ~PMIV_Server(){};
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
  virtual ~PMIV_BasicParameter() {};

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
  virtual void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction) const = 0;

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
  virtual ~PMIV_Thermal() {};

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
public:
  /**
   * constructor
   */
  PMIV_Optical(const PMIC_Environment &env):PMIV_Server(env) { }

  /**
   * destructor
   */
  virtual ~PMIV_Optical() {}

  /**
   * @return the complex refraction index of material
   */
  virtual std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0, PetscScalar Tl=1) const=0;
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
  virtual ~PMIP_Server(){};
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
  virtual ~PMIP_BasicParameter() {};

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
  virtual void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction) const = 0;

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
  virtual ~PMIP_Thermal() {};

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

