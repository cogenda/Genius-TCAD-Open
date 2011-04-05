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

//  $Id: material.h,v 1.10 2008/07/09 05:58:16 gdiso Exp $

#ifndef __material_h__
#define __material_h__

#include <string>
#include <map>

#include "genius_common.h"

#include "material_define.h"
#include "physical_unit.h"
#include "PMI.h"

#ifdef CYGWIN
  class HINSTANCE__; // Forward or never
  typedef HINSTANCE__* HINSTANCE;
#endif

class SimulationRegion;

//using namespace Material;

namespace Material {

  enum PMI_Type {Basic=0, Band, Mobility, Impact, Thermal, Optical, Trap, Invalid_PMI};

  /**
   * convert string to enum
   */
  extern PMI_Type PMI_Type_string_to_enum(const std::string &);


/**
 * the base class of material database interface,
 * it holds the pointer to dynamic loadable library
 * which contains material information
 */
class MaterialBase
{
public:
  /**
   * constructor
   */
  MaterialBase(const SimulationRegion * reg);

  /**
   * destructor, free the dll file pointer
   */
  virtual ~MaterialBase();

  /**
   * mapping Point, its Data and current time to internal image.
   * the PMI has second order pointer to these internal image.
   * so PMI can read information
   */
  void mapping(const Point* point, const FVM_NodeData* node_data, PetscScalar time)
  {
    p_point = point;
    p_node_data = node_data;
    clock = time;
  }

  /**
   * @return PMI_Environment
   */
  PMI_Environment build_PMI_Environment() ;


  /**
   * load the material dll file
   */
  void load_material( const std::string & material );


  /**
   * function pointer to set_ad_number, set the independent variable
   * number of automatically differentiation
   */
  void *              (*set_ad_num)(const unsigned int);

  /**
   * virtual function for set different model, calibrate parameters to PMI
   */
  virtual void set_pmi(const std::string &type, const std::string &model_name,
                       std::vector<Parser::Parameter> & pmi_parameters) = 0;

  /**
   * virtual function for init nodal data
   */
  virtual void init_node(const std::string &type, const Point* point, FVM_NodeData* node_data) = 0;

  /**
   * virtual function for init nodal data with special boundary condition
   */
  virtual void init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data) = 0;

  /**
   * currently activated PMIS models
   */
  std::map<PMI_Type, std::string> active_models;

  /**
   * get an information string of the PMI models
   */
  virtual std::string get_pmi_info(const std::string& type, const int verbosity = 0) = 0;

protected:

  /**
   * const pointer to region
   */
  const SimulationRegion * region;

  /**
   * the material name
   */
  const std::string          material;

  /**
   * pointer to current point, which is updated by mapping function
   */
  const Point                *p_point;

  /**
   * pointer to data of current node, which is updated by mapping function
   */
  const FVM_NodeData         *p_node_data;

  /**
   * current time, which is updated by mapping function
   */
  PetscScalar                clock;

  /**
   * region point based variables
   */
  const std::map<std::string, SimulationVariable>  * point_variables;

  /**
   * region cell based variables
   */
  const std::map<std::string, SimulationVariable>  * cell_variables;

  /**
   * pointer to dynamic loaded library file
   */
#ifdef CYGWIN
  HINSTANCE                  dll_file;
#else
  void                      *dll_file;
#endif

};


/**
 * the derived class for semiconductor material interface
 */
class MaterialSemiconductor: public MaterialBase
{

public:
  /**
   * function pointer to basic semiconductor parameter
   */
  PMIS_BasicParameter *basic;

  /**
   * function pointer to semiconductor band structure parameter
   */
  PMIS_BandStructure  *band;

  /**
   * function pointer to semiconductor mobility parameter
   */
  PMIS_Mobility       *mob;

  /**
   * function pointer to semiconductor impact ionization parameter
   */
  PMIS_Avalanche      *gen;

  /**
   * function pointer to semiconductor thermal parameter
   */
  PMIS_Thermal        *thermal;

  /**
   * function pointer to semiconductor optical parameter
   */
  PMIS_Optical        *optical;

  /**
   * function pointer to semiconductor trapping parameter
   */
  PMIS_Trap           *trap;

public:
  /**
   * constructor, open the material file and point each pointer to corresponding functions
   */
  MaterialSemiconductor(const SimulationRegion * reg);

  /**
   * destructor, unload each pointer
   */
  ~MaterialSemiconductor();

  /**
   * set different model, calibrate parameters in PMI
   */
  void set_pmi(const std::string &type, const std::string &model_name,
               std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * init nodal data for semiconductor region
   */
  void init_node(const std::string &type, const Point* point, FVM_NodeData* node_data);

  /**
   * init nodal data with special boundary condition for semiconductor region
   */
  void init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data);

  /**
   * get an information string of the PMI models
   */
  std::string get_pmi_info(const std::string& type, const int verbosity = 0) ;

};



/**
 * the derived class for insulator material interface
 */
class MaterialInsulator: public MaterialBase
{

public:
  /**
   * function pointer to basic insulator parameter
   */
  PMII_BasicParameter *basic;

  /**
   * function pointer to insulator band structure parameter
   */
  PMII_BandStructure  *band;

  /**
   * function pointer to insulator thermal parameter
   */
  PMII_Thermal        *thermal;

  /**
   * function pointer to insulator optical parameter
   */
  PMII_Optical        *optical;

public:
  /**
   * constructor, open the material file and point each pointer to corresponding functions
   */
  MaterialInsulator(const SimulationRegion * reg);

  /**
   * destructor, unload each pointer
   */
  ~MaterialInsulator();

  /**
   * set different model, calibrate parameters in PMI
   */
  void set_pmi(const std::string &type, const std::string &model_name,
               std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * init nodal data for insulator region
   */
  void init_node(const std::string &type, const Point* point, FVM_NodeData* node_data);

  /**
   * init nodal data with special boundary condition for insulator region
   */
  void init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data);

  /**
   * get an information string of the PMI models
   */
  std::string get_pmi_info(const std::string& type, const int verbosity = 0) ;

};



/**
 * the derived class for conductor material interface
 */
class MaterialConductor: public MaterialBase
{

public:
  /**
   * function pointer to basic conductor parameter
   */
  PMIC_BasicParameter *basic;

  /**
   * function pointer to conductor thermal parameter
   */
  PMIC_Thermal        *thermal;

  /**
   * function pointer to conductor optical parameter
   */
  PMIC_Optical        *optical;

public:
  /**
   * constructor, open the material file and point each pointer to corresponding functions
   */
  MaterialConductor(const SimulationRegion * reg);

  /**
   * destructor, unload each pointer
   */
  ~MaterialConductor();

  /**
   * set different model, calibrate parameters in PMI
   */
  void set_pmi(const std::string &type, const std::string &model_name,
               std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * init nodal data for conductor region
   */
  void init_node(const std::string &type, const Point* point, FVM_NodeData* node_data);

  /**
   * init nodal data with special boundary condition for conductor region
   */
  void init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data);

  /**
   * get an information string of the PMI models
   */
  std::string get_pmi_info(const std::string& type, const int verbosity = 0) ;
};



/**
 * the derived class for vacuum, only used in optical simulation
 */
class MaterialVacuum: public MaterialBase
{

public:
  /**
   * function pointer to basic parameters of vacuum
   */
  PMIV_BasicParameter *basic;

  /**
   * function pointer to thermal parameters of vacuum
   */
  PMIV_Thermal        *thermal;

  /**
   * function pointer to optical parameters of vacuum
   */
  PMIV_Optical        *optical;

public:
  /**
   * constructor, open the material file and point each pointer to corresponding functions
   */
  MaterialVacuum(const SimulationRegion * reg);

  /**
   * destructor, unload each pointer
   */
  ~MaterialVacuum();

  /**
   * set different model, calibrate parameters in PMI
   */
  void set_pmi(const std::string &type, const std::string &model_name,
               std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * init nodal data for conductor region
   */
  void init_node(const std::string &type, const Point* point, FVM_NodeData* node_data);

  /**
   * init nodal data with special boundary condition for conductor region
   */
  void init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data);

  /**
   * get an information string of the PMI models
   */
  std::string get_pmi_info(const std::string& type, const int verbosity = 0) ;
};


/**
 * the derived class for PML, only used in optical simulation
 */
class MaterialPML: public MaterialBase
{

public:
  /**
   * function pointer to basic parameters of PML
   */
  PMIP_BasicParameter *basic;

  /**
   * function pointer to thermal parameters of PML
   */
  PMIP_Thermal        *thermal;

public:
  /**
   * constructor, open the material file and point each pointer to corresponding functions
   */
  MaterialPML(const SimulationRegion * reg);

  /**
   * destructor, unload each pointer
   */
  ~MaterialPML();

  /**
   * set different model, calibrate parameters in PMI
   */
  void set_pmi(const std::string &type, const std::string &model_name,
               std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * init nodal data for conductor region
   */
  void init_node(const std::string &type, const Point* point, FVM_NodeData* node_data);

  /**
   * init nodal data with special boundary condition for conductor region
   */
  void init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data);

  /**
   * get an information string of the PMI models
   */
  std::string get_pmi_info(const std::string& type, const int verbosity = 0) ;

};

} // namespace Material

#endif
