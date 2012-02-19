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


#include <map>
#include <cmath>

#include "enum_solution.h"
#include "physical_unit.h"
using namespace PhysicalUnit;


std::map<std::string, std::string> variable_name_conversion;

static void init_variable_name_conversion_map()
{
  if ( variable_name_conversion.empty() )
  {
    variable_name_conversion["potential"   ] = "potential";
    variable_name_conversion["e.field"     ] = "efield";
    variable_name_conversion["electron"    ] = "electron";
    variable_name_conversion["hole"        ] = "hole";
    variable_name_conversion["temperature" ] = "temperature";
    variable_name_conversion["e.temp"      ] = "elec_temperature";
    variable_name_conversion["h.temp"      ] = "hole_temperature";
    variable_name_conversion["qfn"         ] = "qfn";
    variable_name_conversion["qfp"         ] = "qfp";

    variable_name_conversion["doping.na"   ] = "na";
    variable_name_conversion["doping.nd"   ] = "nd";

    variable_name_conversion["recomb"      ] = "recombination";
    variable_name_conversion["recomb.dir"  ] = "recombination_dir";
    variable_name_conversion["recomb.shr"  ] = "recombination_srh";
    variable_name_conversion["recomb.auger"] = "recombination_auger";

    variable_name_conversion["mole.x"      ] = "mole_x";
    variable_name_conversion["mole.y"      ] = "mole_y";

    variable_name_conversion["ii.gen"             ] = "impact_ionization";
    variable_name_conversion["optical.gen"        ] = "optical_generation";
    variable_name_conversion["particle.gen"       ] = "particle_generation";
  }
}



std::string FormatVariableString(const std::string & var_name)
{
  init_variable_name_conversion_map();

  if( variable_name_conversion.find(var_name) != variable_name_conversion.end() )
    return variable_name_conversion.find(var_name)->second;
  return var_name;
}



//-------------------------------------------------------------------------


std::map<std::string, SolutionVariable> solution_string_to_enum_map;


static void init_solution_string_to_enum_map()
{
  // if the solution_string_to_enum_map not build yet
  if ( solution_string_to_enum_map.empty() )
  {
    solution_string_to_enum_map["potential"   ] = POTENTIAL;
    solution_string_to_enum_map["efield"      ] = E_FIELD;
    solution_string_to_enum_map["electron"    ] = ELECTRON;
    solution_string_to_enum_map["hole"        ] = HOLE;
    solution_string_to_enum_map["temperature" ] = TEMPERATURE;
    solution_string_to_enum_map["elec_temperature"   ] = E_TEMP;
    solution_string_to_enum_map["hole_temperature"   ] = H_TEMP;
    solution_string_to_enum_map["qfn"                ] = QFN;
    solution_string_to_enum_map["qfp"                ] = QFP;

    solution_string_to_enum_map["doping"             ] = DOPING;
    solution_string_to_enum_map["na"                 ] = DOPING_Na;
    solution_string_to_enum_map["nd"                 ] = DOPING_Nd;
    solution_string_to_enum_map["min.carrier"        ] = MIN_CARRIER;
    solution_string_to_enum_map["net.carrier"        ] = NET_CARRIER;
    solution_string_to_enum_map["net.charge"         ] = NET_CHARGE;
    solution_string_to_enum_map["recombination"      ] = RECOMBINATION;
    solution_string_to_enum_map["recombination_dir"  ] = RECOMB_DIR;
    solution_string_to_enum_map["recombination_srh"  ] = RECOMB_SHR;
    solution_string_to_enum_map["recombination_auger"] = RECOMB_AUGER;

    solution_string_to_enum_map["mole_x"      ] = MOLE_X;
    solution_string_to_enum_map["mole_y"      ] = MOLE_Y;

    solution_string_to_enum_map["impact_ionization"  ] = II_GEN;
    solution_string_to_enum_map["optical_generation" ] = OPTICAL_GEN;
    solution_string_to_enum_map["particle_generation"] = PARTICLE_GEN;
  }

}


SolutionVariable solution_string_to_enum ( const std::string & str )
{
  init_solution_string_to_enum_map();

  if ( solution_string_to_enum_map.find ( str ) != solution_string_to_enum_map.end() )
    return solution_string_to_enum_map[str];
  return INVALID_Variable;
}

DataType variable_data_type ( const SolutionVariable v )
{
  switch ( v )
  {
    case POTENTIAL     :
    case ELECTRON      :
    case HOLE          :
    case TEMPERATURE   :
    case E_TEMP        :
    case H_TEMP        :
    case QFN           :
    case QFP           :

    case DOPING        :
    case DOPING_Na     :
    case DOPING_Nd     :
    case MIN_CARRIER   :
    case NET_CARRIER   :
    case NET_CHARGE    :

    case RECOMBINATION :
    case RECOMB_DIR    :
    case RECOMB_SHR    :
    case RECOMB_AUGER  :

    case MOLE_X        :
    case MOLE_Y        :

    case II_GEN        :
    case OPTICAL_GEN   :
    case OPTICAL_HEAT  :
    case PARTICLE_GEN  : return SCALAR;

    case E_FIELD       : return VECTOR;
  }

  return INVALID_DATATYPE;
}


std::string variable_unit_string ( const SolutionVariable v)
{
  switch ( v )
  {
    case POTENTIAL     :  return "V";
    case ELECTRON      :  return "cm^-3";
    case HOLE          :  return "cm^-3";
    case TEMPERATURE   :  return "K";
    case E_TEMP        :  return "K";
    case H_TEMP        :  return "K";
    case QFN           :  return "eV";
    case QFP           :  return "eV";

    case DOPING        :  return "cm^-3";
    case DOPING_Na     :  return "cm^-3";
    case DOPING_Nd     :  return "cm^-3";
    case MIN_CARRIER   :  return "cm^-3";
    case NET_CARRIER   :  return "cm^-3";
    case NET_CHARGE    :  return "cm^-3";

    case RECOMBINATION :  return "cm^-3/s";
    case RECOMB_DIR    :  return "cm^-3/s";
    case RECOMB_SHR    :  return "cm^-3/s";
    case RECOMB_AUGER  :  return "cm^-3/s";

    case MOLE_X        :  return "-";
    case MOLE_Y        :  return "-";

    case II_GEN        :  return "cm^-3/s";
    case OPTICAL_GEN   :  return "cm^-3/s";
    case OPTICAL_HEAT  :  return "eV/s";
    case PARTICLE_GEN  :  return "cm^-3/s";

    case E_FIELD       :  return "V/cm";
  }

  return std::string();
}


double variable_unit ( const SolutionVariable v)
{
  switch ( v )
  {
    case POTENTIAL     :  return V;
    case ELECTRON      :  return pow(cm, -3);
    case HOLE          :  return pow(cm, -3);
    case TEMPERATURE   :  return K;
    case E_TEMP        :  return K;
    case H_TEMP        :  return K;
    case QFN           :  return eV;
    case QFP           :  return eV;

    case DOPING        :  return pow(cm, -3);
    case DOPING_Na     :  return pow(cm, -3);
    case DOPING_Nd     :  return pow(cm, -3);
    case MIN_CARRIER   :  return pow(cm, -3);
    case NET_CARRIER   :  return pow(cm, -3);
    case NET_CHARGE    :  return pow(cm, -3);

    case RECOMBINATION :  return pow(cm, -3)/s;
    case RECOMB_DIR    :  return pow(cm, -3)/s;
    case RECOMB_SHR    :  return pow(cm, -3)/s;
    case RECOMB_AUGER  :  return pow(cm, -3)/s;

    case MOLE_X        :  return 1.0;
    case MOLE_Y        :  return 1.0;

    case II_GEN        :  return pow(cm, -3)/s;
    case OPTICAL_GEN   :  return pow(cm, -3)/s;
    case OPTICAL_HEAT  :  return eV/s;
    case PARTICLE_GEN  :  return pow(cm, -3)/s;

    case E_FIELD       :  return V/cm;
  }

  return 1.0;
}



