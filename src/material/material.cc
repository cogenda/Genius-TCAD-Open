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

#ifdef CYGWIN
  #include <Windows.h>
  #undef max
  #undef min
  #define LDFUN GetProcAddress
#else
  #include <dlfcn.h>
  #define LDFUN dlsym
#endif


#include "material.h"
#include "simulation_region.h"




namespace Material
{

  MaterialBase::MaterialBase(const SimulationRegion * reg)
  : set_ad_num(0),  region(reg) , material(reg->material()), p_point(0), p_node_data(0), dll_file(0)
  {
    point_variables = &(region->region_point_variables());
    cell_variables = &(region->region_cell_variables());
  }

  MaterialBase::~MaterialBase()
  {
#ifdef CYGWIN
    FreeLibrary(dll_file);
#else
    if ( dll_file )
      dlclose( dll_file );
#endif
  }


  PMI_Environment MaterialBase::build_PMI_Environment()
  {
     PMI_Environment env(  &p_point, &p_node_data, &clock, &point_variables,
                            PhysicalUnit::m, PhysicalUnit::s, PhysicalUnit::V, PhysicalUnit::C, PhysicalUnit::K);
     return env;
  }


  void MaterialBase::load_material( const std::string & _material )
  {
#ifdef CYGWIN
    std::string filename =  Genius::genius_dir() + "\\lib\\lib" + _material + ".dll";
#else
    std::string filename =  Genius::genius_dir() + "/lib/lib" + _material + ".so";
#endif

#ifdef CYGWIN
    dll_file = LoadLibrary(filename.c_str());
    if(dll_file==NULL)
    {
      MESSAGE<<"Open material file lib"<< _material <<".dll error." << '\n'; RECORD();
      MESSAGE<<"Error code: " << GetLastError() << '\n'; RECORD();
      genius_error();
    }
#else
    dll_file = dlopen(filename.c_str(), RTLD_LAZY);
    if(dll_file==NULL)
    {
      MESSAGE<<"Open material file lib"<< _material <<".so error." << '\n'; RECORD();
      MESSAGE<<"Error code: " << dlerror() << '\n'; RECORD();
      genius_error();
    }
#endif
  }


  MaterialSemiconductor::MaterialSemiconductor(const SimulationRegion * reg)
      : MaterialBase(reg)
  {
    // open material data base
    std::string _material = FormatMaterialString(material);
    load_material(_material);

    PMI_Environment env = build_PMI_Environment();

    PMIS_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMIS_BandStructure* (*wband)     (const PMI_Environment& env);
    PMIS_Mobility*      (*wmob)      (const PMI_Environment& env);
    PMIS_Avalanche*     (*wgen)      (const PMI_Environment& env);
    PMIS_Thermal*       (*wthermal)  (const PMI_Environment& env);
    PMIS_Optical*       (*woptical)  (const PMI_Environment& env);
    PMIS_Trap*          (*wtrap)     (const PMI_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMIS AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    std::string model_fun_name;

    // init basic parameters for the material
    model_fun_name = "PMIS_" + _material + "_BasicParameter_Default";
    wbasic = (PMIS_BasicParameter* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMIS "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;

    // init band structure model
    model_fun_name = "PMIS_" + _material + "_BandStructure_Default";
    wband =  (PMIS_BandStructure* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wband) { MESSAGE<<"Open PMIS "<< material <<" BandStructure function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Band] = model_fun_name;


    // init mobility model
    model_fun_name = "PMIS_" + _material + "_Mob_Default";
    wmob  =  (PMIS_Mobility* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wmob) { MESSAGE<<"Open PMIS "<< material <<" Mobility function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Mobility] = model_fun_name;


    // init Avalanche generation model
    model_fun_name = "PMIS_" + _material + "_Avalanche_Default";
    wgen  =  (PMIS_Avalanche* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wgen) { MESSAGE<<"Open PMIS "<< material <<" Avalanche function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Impact] = model_fun_name;


    // init Thermal model for lattice temperature equation
    model_fun_name = "PMIS_" + _material + "_Thermal_Default";
    wthermal  = (PMIS_Thermal* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wthermal) { MESSAGE<<"Open PMIS "<< material <<" Thermal function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Thermal] = model_fun_name;


    // init optical data
    model_fun_name = "PMIS_" + _material + "_Optical_Default";
    woptical  = (PMIS_Optical* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!woptical) { MESSAGE<<"Open PMIS "<< material <<" Optical function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Optical] = model_fun_name;

    // init trap data
    model_fun_name = "PMIS_" + _material + "_Trap_Default";
    wtrap = (PMIS_Trap* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
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
                                      std::vector<Parser::Parameter> & pmi_parameters)
  {

    std::string _material = FormatMaterialString(material);

    PMIS_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMIS_BandStructure* (*wband)     (const PMI_Environment& env);
    PMIS_Mobility*      (*wmob)      (const PMI_Environment& env);
    PMIS_Avalanche*     (*wgen)      (const PMI_Environment& env);
    PMIS_Thermal*       (*wthermal)  (const PMI_Environment& env);
    PMIS_Optical*       (*woptical)  (const PMI_Environment& env);
    PMIS_Trap*          (*wtrap)     (const PMI_Environment& env);

    PMI_Environment env = build_PMI_Environment();

    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMIS_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMIS_BasicParameter* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          wband = (PMIS_BandStructure* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          wmob = (PMIS_Mobility* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          wgen = (PMIS_Avalanche* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          wthermal = (PMIS_Thermal* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          woptical = (PMIS_Optical* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          wtrap = (PMIS_Trap* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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



  MaterialInsulator::MaterialInsulator(const SimulationRegion * reg)
      : MaterialBase(reg)
  {
    // open material data base
    std::string _material = FormatMaterialString(material);
    load_material(_material);

    PMI_Environment env = build_PMI_Environment();

    PMII_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMII_BandStructure* (*wband)     (const PMI_Environment& env);
    PMII_Thermal*       (*wthermal)  (const PMI_Environment& env);
    PMII_Optical*       (*woptical)  (const PMI_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMII AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    std::string model_fun_name;

    // init basic parameters for the material
    model_fun_name = "PMII_" + _material + "_BasicParameter_Default";
    wbasic = (PMII_BasicParameter* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMII "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;

    // init band structure model
    model_fun_name = "PMII_" + _material + "_BandStructure_Default";
    wband =  (PMII_BandStructure* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wband) { MESSAGE<<"Open PMII "<< material <<" BandStructure function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Band] = model_fun_name;

    // init Thermal model for lattice temperature equation
    model_fun_name = "PMII_" + _material + "_Thermal_Default";
    wthermal  = (PMII_Thermal* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wthermal) { MESSAGE<<"Open PMII "<< material <<" Thermal function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Thermal] = model_fun_name;

    // init optical data
    model_fun_name = "PMII_" + _material + "_Optical_Default";
    woptical  = (PMII_Optical* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!woptical) { MESSAGE<<"Open PMII "<< material <<" Optical function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Optical] = model_fun_name;

    basic = wbasic(env);
    band  = wband(env);
    thermal  = wthermal(env);
    optical  = woptical(env);
  }


  MaterialInsulator::~MaterialInsulator()
  {
    delete basic;
    delete band;
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
    case Band:
      band->init_node();
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
    case Band:
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
    case Band:
      pmi = band; break;
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
                                  std::vector<Parser::Parameter> & pmi_parameters)
  {
    std::string _material = FormatMaterialString(material);

    PMII_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMII_BandStructure* (*wband)     (const PMI_Environment& env);
    PMII_Thermal*       (*wthermal)  (const PMI_Environment& env);
    PMII_Optical*       (*woptical)  (const PMI_Environment& env);

    PMI_Environment env = build_PMI_Environment();

    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMII_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMII_BasicParameter* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
    case Band     :
      {
        std::string model_fun_name = "PMII_" + _material + "_BandStructure_" + model_name;
        if (active_models[Band] != model_fun_name)
        {
          wband = (PMII_BandStructure* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
          if(!wband) { MESSAGE<<"Open PMII "<< material <<" BandStructure function "<< model_name <<" error!\n"; RECORD(); genius_error(); }
          // delete old PMI object
          delete band;
          // get a new one and do calibrate
          band = wband(env);
          active_models[Band] = model_fun_name;
        }
        if(band->calibrate(pmi_parameters))
        {
          MESSAGE<<"WARNING: PMI "<< material <<" BandStructure Parameter calibrating has mismatch(es)!\n"; RECORD();
        }
        break;
      }
    case Thermal   :
      {
        std::string model_fun_name = "PMII_" + _material + "_Thermal_" + model_name;
        if (active_models[Thermal] != model_fun_name)
        {
          wthermal = (PMII_Thermal* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          woptical = (PMII_Optical* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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



  MaterialConductor::MaterialConductor(const SimulationRegion * reg)
      : MaterialBase(reg)
  {
    // open material data base
    std::string _material = FormatMaterialString(material);
    load_material(_material);

    PMI_Environment env = build_PMI_Environment();

    PMIC_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMIC_Thermal*       (*wthermal)  (const PMI_Environment& env);
    PMIC_Optical*       (*woptical)  (const PMI_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMIC AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    std::string model_fun_name;

    // init basic parameters for the material
    model_fun_name = "PMIC_" + _material + "_BasicParameter_Default";
    wbasic = (PMIC_BasicParameter* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMIC "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;

    // init Thermal model for lattice temperature equation
    model_fun_name = "PMIC_" + _material + "_Thermal_Default";
    wthermal  = (PMIC_Thermal* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wthermal) { MESSAGE<<"Open PMIC "<< material <<" Thermal function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Thermal] = model_fun_name;

    // init optical data
    model_fun_name = "PMIC_" + _material + "_Optical_Default";
    woptical  = (PMIC_Optical* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
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
                                  std::vector<Parser::Parameter> & pmi_parameters)
  {

    std::string _material = FormatMaterialString(material);

    PMIC_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMIC_Thermal*       (*wthermal)  (const PMI_Environment& env);
    PMIC_Optical*       (*woptical)  (const PMI_Environment& env);

    PMI_Environment env = build_PMI_Environment();

    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMIC_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMIC_BasicParameter* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          wthermal = (PMIC_Thermal* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          woptical = (PMIC_Optical* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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


  MaterialVacuum::MaterialVacuum(const SimulationRegion * reg)
      : MaterialBase(reg)
  {
    // open material data base
    std::string _material = FormatMaterialString(material);
    load_material(_material);

    PMI_Environment env = build_PMI_Environment();

    PMIV_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMIV_Thermal*       (*wthermal)  (const PMI_Environment& env);
    PMIV_Optical*       (*woptical)  (const PMI_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMIV AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    std::string model_fun_name;

    // init basic parameters for the material
    model_fun_name = "PMIV_" + _material + "_BasicParameter_Default";
    wbasic = (PMIV_BasicParameter* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMIV "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;

    // init Thermal model for lattice temperature equation
    model_fun_name = "PMIV_" + _material + "_Thermal_Default";
    wthermal  = (PMIV_Thermal* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wthermal) { MESSAGE<<"Open PMIV "<< material <<" Thermal function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Thermal] = model_fun_name;

    // init optical data
    model_fun_name = "PMIV_" + _material + "_Optical_Default";
    woptical  = (PMIV_Optical* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
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
                               std::vector<Parser::Parameter> & pmi_parameters)
  {

    std::string _material = FormatMaterialString(material);

    PMIV_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMIV_Thermal*       (*wthermal)  (const PMI_Environment& env);
    PMIV_Optical*       (*woptical)  (const PMI_Environment& env);

    PMI_Environment env = build_PMI_Environment();

    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMIV_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMIV_BasicParameter* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          wthermal = (PMIV_Thermal* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          woptical = (PMIV_Optical* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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


  MaterialPML::MaterialPML(const SimulationRegion * reg)
      : MaterialBase(reg)
  {

    // open material data base
    std::string _material = FormatMaterialString(material);
    load_material(_material);

    PMI_Environment env = build_PMI_Environment();

    PMIP_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMIP_Thermal*       (*wthermal)  (const PMI_Environment& env);

    //init AD indepedent variable set routine
    set_ad_num = (void* (*) (const unsigned int))LDFUN(dll_file,"set_ad_number");
    if(!set_ad_num) { MESSAGE<<"Open PMIP AD_SET_VARIABLE function error!\n"; RECORD(); genius_error();}

    std::string model_fun_name;

    // init basic parameters for the material
    model_fun_name = "PMIP_" + _material + "_BasicParameter_Default";
    wbasic = (PMIP_BasicParameter* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
    if(!wbasic) { MESSAGE<<"Open PMIP "<< material <<" BasicParameter function "<< "Default" <<" error!\n"; RECORD(); genius_error(); }
    active_models[Basic] = model_fun_name;

    // init Thermal model for lattice temperature equation
    model_fun_name = "PMIP_" + _material + "_Thermal_Default";
    wthermal  = (PMIP_Thermal* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());
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
                            std::vector<Parser::Parameter> & pmi_parameters)
  {

    std::string _material = FormatMaterialString(material);

    PMIP_BasicParameter*(*wbasic)    (const PMI_Environment& env);
    PMIP_Thermal*       (*wthermal)  (const PMI_Environment& env);

    PMI_Environment env = build_PMI_Environment();

    switch ( PMI_Type_string_to_enum(type) )
    {
    case Basic     :
      {
        std::string model_fun_name = "PMIP_" + _material + "_BasicParameter_" + model_name;
        if (active_models[Basic] != model_fun_name)
        {
          wbasic = (PMIP_BasicParameter* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
          wthermal = (PMIP_Thermal* (*) (const PMI_Environment& env))LDFUN( dll_file, model_fun_name.c_str() );
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
