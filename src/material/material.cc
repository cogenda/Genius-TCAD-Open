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

//  $Id: material.cc,v 1.12 2008/07/09 05:58:16 gdiso Exp $

#include <cstdlib>
#include <cmath>
#include <iomanip>

#include "genius_env.h"
#include "genius_common.h"
#include "log.h"
#include "material_define.h"
#include "material.h"


#ifdef CYGWIN
  #define LDFUN GetProcAddress
#else
  #define LDFUN dlsym
#endif


namespace Material
{


  MaterialSemiconductor::MaterialSemiconductor(const std::string & mat, const std::string &reg)
      : MaterialBase(mat, reg)
  {

    PMIS_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);

    std::string _material = FormatMaterialString(material);
    std::string filename =  Genius::genius_dir() + "/lib/lib" + _material + ".so";
    std::string model_fun_name;

#ifdef CYGWIN
    dll_file = LoadLibrary(filename.c_str());
#else
    dll_file = dlopen(filename.c_str(), RTLD_LAZY);
#endif

    if(dll_file==NULL)
    { MESSAGE<<"Open material file lib"<< _material <<".so error." << '\n'; RECORD(); genius_error();}

    PMIS_BasicParameter*(*wbasic)    (const PMIS_Environment& env);
    PMIS_BandStructure* (*wband)     (const PMIS_Environment& env);
    PMIS_Mobility*      (*wmob)      (const PMIS_Environment& env);
    PMIS_Avalanche*     (*wgen)      (const PMIS_Environment& env);
    PMIS_Thermal*       (*wthermal)  (const PMIS_Environment& env);
    PMIS_Optical*       (*woptical)  (const PMIS_Environment& env);
    PMIS_Trap*          (*wtrap)     (const PMIS_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMIS AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    // init basic parameters for the material
    model_fun_name = "PMIS_" + _material + "_BasicParameter_Default";
    wbasic = (PMIS_BasicParameter* (*) (const PMIS_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMIS "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;

    // init band structure model
    model_fun_name = "PMIS_" + _material + "_BandStructure_Default";
    wband =  (PMIS_BandStructure* (*) (const PMIS_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wband) { MESSAGE<<"Open PMIS "<< material <<" BandStructure function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Band] = model_fun_name;


    // init mobility model
    model_fun_name = "PMIS_" + _material + "_Mob_Default";
    wmob  =  (PMIS_Mobility* (*) (const PMIS_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wmob) { MESSAGE<<"Open PMIS "<< material <<" Mobility function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Mobility] = model_fun_name;


    // init Avalanche generation model
    model_fun_name = "PMIS_" + _material + "_Avalanche_Default";
    wgen  =  (PMIS_Avalanche* (*) (const PMIS_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wgen) { MESSAGE<<"Open PMIS "<< material <<" Avalanche function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Impact] = model_fun_name;


    // init Thermal model for lattice temperature equation
    model_fun_name = "PMIS_" + _material + "_Thermal_Default";
    wthermal  = (PMIS_Thermal* (*) (const PMIS_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wthermal) { MESSAGE<<"Open PMIS "<< material <<" Thermal function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Thermal] = model_fun_name;


    // init optical data
    model_fun_name = "PMIS_" + _material + "_Optical_Default";
    woptical  = (PMIS_Optical* (*) (const PMIS_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!woptical) { MESSAGE<<"Open PMIS "<< material <<" Optical function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Optical] = model_fun_name;

    // init trap data
    model_fun_name = "PMIS_" + _material + "_Trap_Default";
    wtrap = (PMIS_Trap* (*) (const PMIS_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wtrap ) { MESSAGE<<"Open PMIS "<< material <<" Trap function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Trap] = model_fun_name;

    basic = wbasic(env);
    band  = wband(env);
    mob   = wmob(env);
    gen   = wgen(env);
    thermal  = wthermal(env);
    optical  = woptical(env);
    trap  = wtrap(env);

  }


  MaterialSemiconductor::~MaterialSemiconductor()
  {
    delete basic;
    delete band;
    delete mob;
    delete gen;
    delete thermal;
    delete optical;
    delete trap;
  }

  void MaterialSemiconductor::init_node(const std::string &type, const Point* point, FVM_NodeData* node_data)
  {
    mapping(point, node_data, clock);
    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic:
      basic->init_node();
      break;
    case Band:
      band->init_node();
      break;
    case Mobility:
      mob->init_node();
      break;
    case Impact:
      gen->init_node();
      break;
    case Thermal:
      thermal->init_node();
      break;
    case Optical:
      optical->init_node();
      break;
    case Trap:
      trap->init_node();
      break;
    default: genius_error();
    }
  }

  void MaterialSemiconductor::init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data)
  {
    this->mapping(point,node_data,clock);

    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
    case Band:
    case Mobility:
    case Impact:
    case Thermal:
    case Optical:
      break;
    case Trap:
      trap->init_bc_node(bc_label);
      break;
    default: genius_error();
    }

  }

  std::string MaterialSemiconductor::get_pmi_info(const std::string& type, const int verbosity)
  {
    std::string _material = FormatMaterialString(material);
    std::stringstream output;

    PMI_Server* pmi;
    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
      pmi = basic; break;
    case Band:
      pmi = band; break;
    case Mobility:
      pmi = mob; break;
    case Impact:
      pmi = gen; break;
    case Thermal:
      pmi = thermal; break;
    case Optical:
      pmi = optical; break;
    case Trap:
      pmi = trap; break;
    default: genius_error();
    }
    output << pmi->get_PMI_info() << std::endl;
    output << pmi->get_parameter_string(verbosity) ;

    return output.str();
  }

  void MaterialSemiconductor::set_pmi(const std::string &type, const std::string &model_name,
                                      const std::vector<Parser::Parameter> & pmi_parameters)
  {

    std::string _material = FormatMaterialString(material);

    PMIS_BasicParameter*(*wbasic)    (const PMIS_Environment& env);
    PMIS_BandStructure* (*wband)     (const PMIS_Environment& env);
    PMIS_Mobility*      (*wmob)      (const PMIS_Environment& env);
    PMIS_Avalanche*     (*wgen)      (const PMIS_Environment& env);
    PMIS_Thermal*       (*wthermal)  (const PMIS_Environment& env);
    PMIS_Optical*       (*woptical)  (const PMIS_Environment& env);
    PMIS_Trap*          (*wtrap)     (const PMIS_Environment& env);

    PMIS_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);

    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMIS_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMIS_BasicParameter* (*) (const PMIS_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wbasic) { MESSAGE<<"Open PMIS "<< material <<" BasicParameter function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete basic;
          // get a new one and do calibrate
          basic = wbasic(env);
          active_models[Basic] = model_fun_name;
        }
        if(basic->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Basic Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Band      :
      {
        std::string model_fun_name = "PMIS_" + _material + "_BandStructure_" + model_name;
        if (active_models[Band] != model_fun_name)
        {
          wband = (PMIS_BandStructure* (*) (const PMIS_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wband) { MESSAGE<<"Open PMIS "<< material <<" BandStructure function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete band;
          // get a new one and do calibrate
          band = wband(env);
          active_models[Band] = model_fun_name;
        }
        if(band->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Band Structure calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Mobility  :
      {
        std::string model_fun_name = "PMIS_" + _material + "_Mob_" + model_name;
        if (active_models[Mobility] != model_fun_name)
        {
          wmob = (PMIS_Mobility* (*) (const PMIS_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wmob) { MESSAGE<<"Open PMIS "<< material <<" Mobility function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete mob;
          // get a new one and do calibrate
          mob = wmob(env);
          active_models[Mobility] = model_fun_name;
        }
        if(mob->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Mobility calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Impact    :
      {
        std::string model_fun_name = "PMIS_" + _material + "_Avalanche_" + model_name;
        if (active_models[Impact] != model_fun_name)
        {
          wgen = (PMIS_Avalanche* (*) (const PMIS_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wgen) { MESSAGE<<"Open PMIS "<< material <<" Avalanche function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete gen;
          // get a new one and do calibrate
          gen = wgen(env);
          active_models[Impact] = model_fun_name;
        }
        if(gen->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Impact Ionization calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Thermal   :
      {
        std::string model_fun_name = "PMIS_" + _material + "_Thermal_" + model_name;
        if (active_models[Thermal] != model_fun_name)
        {
          wthermal = (PMIS_Thermal* (*) (const PMIS_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wthermal) { MESSAGE<<"Open PMIS "<< material <<" Thermal function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete thermal;
          // get a new one and do calibrate
          thermal = wthermal(env);
          active_models[Thermal] = model_fun_name;
        }
        if(thermal->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Thermal Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Optical   :
      {
        std::string model_fun_name = "PMIS_" + _material + "_Optical_" + model_name;
        if (active_models[Optical] != model_fun_name)
        {
          woptical = (PMIS_Optical* (*) (const PMIS_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!woptical) { MESSAGE<<"Open PMIS "<< material <<" Optical function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete optical;
          // get a new one and do calibrate
          optical = woptical(env);
          active_models[Optical] = model_fun_name;
        }
        if(optical->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Optical Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Trap      :
      {
        std::string model_fun_name = "PMIS_" + _material + "_Trap_" + model_name;
        if (active_models[Trap] != model_fun_name)
        {
          wtrap = (PMIS_Trap* (*) (const PMIS_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wtrap) { MESSAGE<<"Open PMIS "<< material <<" Trap function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete trap;
          // get a new one and do calibrate
          trap = wtrap(env);
          active_models[Trap] = model_fun_name;
        }
        if(trap->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Trap Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    default: genius_error();
    }
  }

  //-----------------------------------------------------------------------------------------------------------



  MaterialInsulator::MaterialInsulator(const std::string & mat, const std::string &reg)
      : MaterialBase(mat, reg)
  {


    PMII_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);

    std::string _material = FormatMaterialString(material);
    std::string filename =  Genius::genius_dir() + "/lib/lib" + _material + ".so";
    std::string model_fun_name;

#ifdef CYGWIN
    dll_file = LoadLibrary(filename.c_str());
#else
    dll_file = dlopen(filename.c_str(), RTLD_LAZY);
#endif

    if(dll_file==NULL)
    { MESSAGE<<"Open material file lib"<< _material <<".so error."<<'\n'; RECORD(); genius_error();}

    PMII_BasicParameter*(*wbasic)    (const PMII_Environment& env);
    PMII_Thermal*       (*wthermal)  (const PMII_Environment& env);
    PMII_Optical*       (*woptical)  (const PMII_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMII AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    // init basic parameters for the material
    model_fun_name = "PMII_" + _material + "_BasicParameter_Default";
    wbasic = (PMII_BasicParameter* (*) (const PMII_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMII "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;


    // init Thermal model for lattice temperature equation
    model_fun_name = "PMII_" + _material + "_Thermal_Default";
    wthermal  = (PMII_Thermal* (*) (const PMII_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wthermal) { MESSAGE<<"Open PMII "<< material <<" Thermal function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Thermal] = model_fun_name;

    // init optical data
    model_fun_name = "PMII_" + _material + "_Optical_Default";
    woptical  = (PMII_Optical* (*) (const PMII_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!woptical) { MESSAGE<<"Open PMII "<< material <<" Optical function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Optical] = model_fun_name;

    basic = wbasic(env);
    thermal  = wthermal(env);
    optical  = woptical(env);
  }


  MaterialInsulator::~MaterialInsulator()
  {
    delete basic;
    delete thermal;
    delete optical;
  }

  void MaterialInsulator::init_node(const std::string &type, const Point* point, FVM_NodeData* node_data)
  {
    mapping(point, node_data, clock);
    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic:
      basic->init_node();
      break;
    case Thermal:
      thermal->init_node();
      break;
    case Optical:
      optical->init_node();
      break;
    default: genius_error();
    }
  }

  void MaterialInsulator::init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data)
  {
    genius_assert(bc_label.length()); //prevent compiler warning
    this->mapping(point, node_data, clock);

    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
    case Thermal:
    case Optical:
      break;
    default: genius_error();
    }
  }

  std::string MaterialInsulator::get_pmi_info(const std::string& type, const int verbosity)
  {
    std::string _material = FormatMaterialString(material);
    std::stringstream output;

    PMI_Server* pmi;
    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
      pmi = basic; break;
    case Thermal:
      pmi = thermal; break;
    case Optical:
      pmi = optical; break;
    default: genius_error();
    }
    output << pmi->get_PMI_info() << std::endl;
    output << pmi->get_parameter_string(verbosity) ;

    return output.str();
  }

  void MaterialInsulator::set_pmi(const std::string &type, const std::string &model_name,
                                  const std::vector<Parser::Parameter> & pmi_parameters)
  {
    std::string _material = FormatMaterialString(material);

    PMII_BasicParameter*(*wbasic)    (const PMII_Environment& env);
    PMII_Thermal*       (*wthermal)  (const PMII_Environment& env);
    PMII_Optical*       (*woptical)  (const PMII_Environment& env);

    PMII_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);
    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMII_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMII_BasicParameter* (*) (const PMII_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wbasic) { MESSAGE<<"Open PMII "<< material <<" BasicParameter function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete basic;
          // get a new one and do calibrate
          basic = wbasic(env);
          active_models[Basic] = model_fun_name;
        }
        if(basic->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Basic Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Thermal   :
      {
        std::string model_fun_name = "PMII_" + _material + "_Thermal_" + model_name;
        if (active_models[Thermal] != model_fun_name)
        {
          wthermal = (PMII_Thermal* (*) (const PMII_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wthermal) { MESSAGE<<"Open PMII "<< material <<" Thermal function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete thermal;
          // get a new one and do calibrate
          thermal = wthermal(env);
          active_models[Thermal] = model_fun_name;
        }
        if(thermal->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Thermal Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Optical   :
      {
        std::string model_fun_name = "PMII_" + _material + "_Optical_" + model_name;
        if (active_models[Optical] != model_fun_name)
        {
          woptical = (PMII_Optical* (*) (const PMII_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!woptical) { MESSAGE<<"Open PMII "<< material <<" Optical function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete optical;
          // get a new one and do calibrate
          optical = woptical(env);
          active_models[Optical] = model_fun_name;
        }
        if(optical->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Optical Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    default: genius_error();
    }
  }

  //-----------------------------------------------------------------------------------------------------------



  MaterialConductor::MaterialConductor(const std::string & mat, const std::string &reg)
      : MaterialBase(mat, reg)
  {

    PMIC_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);

    std::string _material = FormatMaterialString(material);
    std::string filename =  Genius::genius_dir() + "/lib/lib" + _material + ".so";
    std::string model_fun_name;

#ifdef CYGWIN
    dll_file = LoadLibrary(filename.c_str());
#else
    dll_file = dlopen(filename.c_str(), RTLD_LAZY);
#endif

    if(dll_file==NULL)
    { MESSAGE<<"Open material file lib"<< _material <<".so error."<<'\n'; RECORD(); genius_error();}

    PMIC_BasicParameter*(*wbasic)    (const PMIC_Environment& env);
    PMIC_Thermal*       (*wthermal)  (const PMIC_Environment& env);
    PMIC_Optical*       (*woptical)  (const PMIC_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMIC AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    // init basic parameters for the material
    model_fun_name = "PMIC_" + _material + "_BasicParameter_Default";
    wbasic = (PMIC_BasicParameter* (*) (const PMIC_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMIC "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;

    // init Thermal model for lattice temperature equation
    model_fun_name = "PMIC_" + _material + "_Thermal_Default";
    wthermal  = (PMIC_Thermal* (*) (const PMIC_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wthermal) { MESSAGE<<"Open PMIC "<< material <<" Thermal function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Thermal] = model_fun_name;

    // init optical data
    model_fun_name = "PMIC_" + _material + "_Optical_Default";
    woptical  = (PMIC_Optical* (*) (const PMIC_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!woptical) { MESSAGE<<"Open PMIC "<< material <<" Optical function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Optical] = model_fun_name;

    basic = wbasic(env);
    thermal  = wthermal(env);
    optical  = woptical(env);
  }


  MaterialConductor::~MaterialConductor()
  {
    delete basic;
    delete thermal;
    delete optical;
  }

  void MaterialConductor::init_node(const std::string &type, const Point* point, FVM_NodeData* node_data)
  {
    mapping(point, node_data, clock);
    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic:
      basic->init_node();
      break;
    case Thermal:
      thermal->init_node();
      break;
    case Optical:
      optical->init_node();
      break;
    default: genius_error();
    }
  }

  void MaterialConductor::init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data)
  {
    genius_assert(bc_label.length()); //prevent compiler warning

    this->mapping(point,node_data,clock);

    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
    case Thermal:
    case Optical:
      break;
    default: genius_error();
    }
  }

  std::string MaterialConductor::get_pmi_info(const std::string& type, const int verbosity)
  {
    std::string _material = FormatMaterialString(material);
    std::stringstream output;

    PMI_Server* pmi;
    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
      pmi = basic; break;
    case Thermal:
      pmi = thermal; break;
    case Optical:
      pmi = optical; break;
    default: genius_error();
    }
    output << pmi->get_PMI_info() << std::endl;
    output << pmi->get_parameter_string(verbosity) ;

    return output.str();
  }


  void MaterialConductor::set_pmi(const std::string &type, const std::string &model_name,
                                  const std::vector<Parser::Parameter> & pmi_parameters)
  {

    std::string _material = FormatMaterialString(material);

    PMIC_BasicParameter*(*wbasic)    (const PMIC_Environment& env);
    PMIC_Thermal*       (*wthermal)  (const PMIC_Environment& env);
    PMIC_Optical*       (*woptical)  (const PMIC_Environment& env);

    PMIC_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);

    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMIC_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMIC_BasicParameter* (*) (const PMIC_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wbasic) { MESSAGE<<"Open PMIC "<< material <<" BasicParameter function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete basic;
          // get a new one and do calibrate
          basic = wbasic(env);
          active_models[Basic] = model_fun_name;
        }
        if(basic->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Basic Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Thermal   :
      {
        std::string model_fun_name = "PMIC_" + _material + "_Thermal_" + model_name;
        if (active_models[Thermal] != model_fun_name)
        {
          wthermal = (PMIC_Thermal* (*) (const PMIC_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wthermal) { MESSAGE<<"Open PMIC "<< material <<" Thermal function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete thermal;
          // get a new one and do calibrate
          thermal = wthermal(env);
          active_models[Thermal] = model_fun_name;
        }
        if(thermal->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Thermal Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Optical   :
      {
        std::string model_fun_name = "PMIC_" + _material + "_Optical_" + model_name;
        if (active_models[Optical] != model_fun_name)
        {
          woptical = (PMIC_Optical* (*) (const PMIC_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!woptical) { MESSAGE<<"Open PMIC "<< material <<" Optical function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete optical;
          // get a new one and do calibrate
          optical = woptical(env);
          active_models[Optical] = model_fun_name;
        }
        if(optical->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Optical Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    default: genius_error();
    }

  }



  //-----------------------------------------------------------------------------------------------------------


  MaterialVacuum::MaterialVacuum(const std::string & mat, const std::string &reg)
      : MaterialBase(mat, reg)
  {

    PMIV_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);

    std::string _material = FormatMaterialString(material);
    std::string filename =  Genius::genius_dir() + "/lib/lib" + _material + ".so";
    std::string model_fun_name;

#ifdef CYGWIN
    dll_file = LoadLibrary(filename.c_str());
#else
    dll_file = dlopen(filename.c_str(), RTLD_LAZY);
#endif

    if(dll_file==NULL)
    { MESSAGE<<"Open material file lib"<< _material <<".so error."<<'\n'; RECORD(); genius_error();}

    PMIV_BasicParameter*(*wbasic)    (const PMIV_Environment& env);
    PMIV_Thermal*       (*wthermal)  (const PMIV_Environment& env);
    PMIV_Optical*       (*woptical)  (const PMIV_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMIV AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    // init basic parameters for the material
    model_fun_name = "PMIV_" + _material + "_BasicParameter_Default";
    wbasic = (PMIV_BasicParameter* (*) (const PMIV_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMIV "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;

    // init Thermal model for lattice temperature equation
    model_fun_name = "PMIV_" + _material + "_Thermal_Default";
    wthermal  = (PMIV_Thermal* (*) (const PMIV_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wthermal) { MESSAGE<<"Open PMIV "<< material <<" Thermal function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Thermal] = model_fun_name;

    // init optical data
    model_fun_name = "PMIV_" + _material + "_Optical_Default";
    woptical  = (PMIV_Optical* (*) (const PMIV_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!woptical) { MESSAGE<<"Open PMIV "<< material <<" Optical function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Optical] = model_fun_name;

    basic = wbasic(env);
    thermal  = wthermal(env);
    optical  = woptical(env);
  }


  MaterialVacuum::~MaterialVacuum()
  {
    delete basic;
    delete thermal;
    delete optical;
  }

  void MaterialVacuum::init_node(const std::string &type, const Point* point, FVM_NodeData* node_data)
  {
    mapping(point, node_data, clock);
    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic:
      basic->init_node();
      break;
    case Thermal:
      thermal->init_node();
      break;
    case Optical:
      optical->init_node();
      break;
    default: genius_error();
    }
  }

  void MaterialVacuum::init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data)
  {
    genius_assert(bc_label.length()); //prevent compiler warning

    this->mapping(point,node_data,clock);

    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
    case Thermal:
    case Optical:
      break;
    default: genius_error();
    }
  }

  std::string MaterialVacuum::get_pmi_info(const std::string& type, const int verbosity)
  {
    std::string _material = FormatMaterialString(material);
    std::stringstream output;

    PMI_Server* pmi;
    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
      pmi = basic; break;
    case Thermal:
      pmi = thermal; break;
    case Optical:
      pmi = optical; break;
    default: genius_error();
    }
    output << pmi->get_PMI_info() << std::endl;
    output << pmi->get_parameter_string(verbosity) ;

    return output.str();
  }

  void MaterialVacuum::set_pmi(const std::string &type, const std::string &model_name,
                               const std::vector<Parser::Parameter> & pmi_parameters)
  {

    std::string _material = FormatMaterialString(material);

    PMIV_BasicParameter*(*wbasic)    (const PMIV_Environment& env);
    PMIV_Thermal*       (*wthermal)  (const PMIV_Environment& env);
    PMIV_Optical*       (*woptical)  (const PMIV_Environment& env);

    PMIV_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);

    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMIV_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMIV_BasicParameter* (*) (const PMIV_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wbasic) { MESSAGE<<"Open PMIV "<< material <<" BasicParameter function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete basic;
          // get a new one and do calibrate
          basic = wbasic(env);
          active_models[Basic] = model_fun_name;
        }
        if(basic->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Basic Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Thermal   :
      {
        std::string model_fun_name = "PMIV_" + _material + "_Thermal_" + model_name;
        if (active_models[Thermal] != model_fun_name)
        {
          wthermal = (PMIV_Thermal* (*) (const PMIV_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wthermal) { MESSAGE<<"Open PMIV "<< material <<" Thermal function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete thermal;
          // get a new one and do calibrate
          thermal = wthermal(env);
          active_models[Thermal] = model_fun_name;
        }
        if(thermal->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Thermal Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Optical   :
      {
        std::string model_fun_name = "PMIV_" + _material + "_Optical_" + model_name;
        if (active_models[Optical] != model_fun_name)
        {
          woptical = (PMIV_Optical* (*) (const PMIV_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!woptical) { MESSAGE<<"Open PMIV "<< material <<" Optical function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete optical;
          // get a new one and do calibrate
          optical = woptical(env);
          active_models[Optical] = model_fun_name;
        }
        if(optical->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Optical Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    default: genius_error();
    }

  }



  //-----------------------------------------------------------------------------------------------------------


  MaterialPML::MaterialPML(const std::string & mat, const std::string &reg)
      : MaterialBase(mat, reg)
  {

    PMIP_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);

    std::string _material = FormatMaterialString(material);
    std::string filename =  Genius::genius_dir() + "/lib/lib" + _material + ".so";
    std::string model_fun_name;

#ifdef CYGWIN
    dll_file = LoadLibrary(filename.c_str());
#else
    dll_file = dlopen(filename.c_str(), RTLD_LAZY);
#endif

    if(dll_file==NULL)
    { MESSAGE<<"Open material file lib"<< _material <<".so error."<<'\n'; RECORD(); genius_error();}

    PMIP_BasicParameter*(*wbasic)    (const PMIP_Environment& env);
    PMIP_Thermal*       (*wthermal)  (const PMIP_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMIP AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    // init basic parameters for the material
    model_fun_name = "PMIP_" + _material + "_BasicParameter_Default";
    wbasic = (PMIP_BasicParameter* (*) (const PMIP_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMIP "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;

    // init Thermal model for lattice temperature equation
    model_fun_name = "PMIP_" + _material + "_Thermal_Default";
    wthermal  = (PMIP_Thermal* (*) (const PMIP_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wthermal) { MESSAGE<<"Open PMIP "<< material <<" Thermal function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Thermal] = model_fun_name;

    basic = wbasic(env);
    thermal  = wthermal(env);
  }


  MaterialPML::~MaterialPML()
  {
    delete basic;
    delete thermal;
  }

  void MaterialPML::init_node(const std::string &type, const Point* point, FVM_NodeData* node_data)
  {
    mapping(point, node_data, clock);
    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic:
      basic->init_node();
      break;
    case Thermal:
      thermal->init_node();
      break;
    default: genius_error();
    }
  }

  void MaterialPML::init_bc_node(const std::string &type, const std::string & bc_label, const Point* point, FVM_NodeData* node_data)
  {
    genius_assert(bc_label.length()); //prevent compiler warning

    this->mapping(point,node_data,clock);

    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
    case Thermal:
      break;
    default: genius_error();
    }
  }

  std::string MaterialPML::get_pmi_info(const std::string& type, const int verbosity)
  {
    std::string _material = FormatMaterialString(material);
    std::stringstream output;

    PMI_Server* pmi;
    switch(PMI_Type_string_to_enum(type))
    {
    case Basic:
      pmi = basic; break;
    case Thermal:
      pmi = thermal; break;
    default: genius_error();
    }
    output << pmi->get_PMI_info() << std::endl;
    output << pmi->get_parameter_string(verbosity) ;

    return output.str();
  }

  void MaterialPML::set_pmi(const std::string &type, const std::string &model_name,
                               const std::vector<Parser::Parameter> & pmi_parameters)
  {

    std::string _material = FormatMaterialString(material);

    PMIP_BasicParameter*(*wbasic)    (const PMIP_Environment& env);
    PMIP_Thermal*       (*wthermal)  (const PMIP_Environment& env);

    PMIP_Environment env(&p_point, &p_node_data, &clock, PhysicalUnit::m, PhysicalUnit::s,
                         PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);

    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMIP_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMIP_BasicParameter* (*) (const PMIP_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wbasic) { MESSAGE<<"Open PMIP "<< material <<" BasicParameter function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete basic;
          // get a new one and do calibrate
          basic = wbasic(env);
          active_models[Basic] = model_fun_name;
        }
        if(basic->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Basic Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Thermal   :
      {
        std::string model_fun_name = "PMIP_" + _material + "_Thermal_" + model_name;
        if (active_models[Thermal] != model_fun_name)
        {
          wthermal = (PMIP_Thermal* (*) (const PMIP_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wthermal) { MESSAGE<<"Open PMIP "<< material <<" Thermal function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete thermal;
          // get a new one and do calibrate
          thermal = wthermal(env);
          active_models[Thermal] = model_fun_name;
        }
        if(thermal->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" Thermal Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    default: genius_error();
    }

  }

  //-------------------------------------------------------------------------------------------------------
  std::map<const std::string, PMI_Type> PMI_name_to_PMI_type;

  static void init_PMI_name_to_PMI_type()
  {
    if( PMI_name_to_PMI_type.empty() )
    {
      PMI_name_to_PMI_type["basic"    ] = Basic;
      PMI_name_to_PMI_type["band"     ] = Band;
      PMI_name_to_PMI_type["mobility" ] = Mobility;
      PMI_name_to_PMI_type["impact"   ] = Impact;
      PMI_name_to_PMI_type["thermal"  ] = Thermal;
      PMI_name_to_PMI_type["optical"  ] = Optical;
      PMI_name_to_PMI_type["trap"     ] = Trap;
    }
  }

  PMI_Type PMI_Type_string_to_enum(const std::string & PMI_name)
  {
    init_PMI_name_to_PMI_type();

    if( PMI_name_to_PMI_type.find(PMI_name)!=PMI_name_to_PMI_type.end() )
      return PMI_name_to_PMI_type[PMI_name];
    else
      return Invalid_PMI;
  }

}
