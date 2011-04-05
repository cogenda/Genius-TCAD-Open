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


std::map<const std::string, SolutionVariable> solution_string_to_enum_map;


static void init_solution_string_to_enum_map()
{
  // if the solution_string_to_enum_map not build yet
  if ( solution_string_to_enum_map.empty() )
  {
    solution_string_to_enum_map["potential"   ] = POTENTIAL;
    solution_string_to_enum_map["e.field"     ] = E_FIELD;
    solution_string_to_enum_map["electron"    ] = ELECTRON;
    solution_string_to_enum_map["hole"        ] = HOLE;
    solution_string_to_enum_map["temperature" ] = TEMPERATURE;
    solution_string_to_enum_map["e.temp"      ] = E_TEMP;
    solution_string_to_enum_map["h.temp"      ] = H_TEMP;
    solution_string_to_enum_map["qfn"         ] = QFN;
    solution_string_to_enum_map["qfp"         ] = QFP;

    solution_string_to_enum_map["doping"      ] = DOPING;
    solution_string_to_enum_map["doping.na"   ] = DOPING_Na;
    solution_string_to_enum_map["doping.nd"   ] = DOPING_Nd;
    solution_string_to_enum_map["min.carrier" ] = MIN_CARRIER;
    solution_string_to_enum_map["net.carrier" ] = NET_CARRIER;
    solution_string_to_enum_map["net.charge"  ] = NET_CHARGE;

    solution_string_to_enum_map["recomb"      ] = RECOMBINATION;
    solution_string_to_enum_map["recomb.dir"  ] = RECOMB_DIR;
    solution_string_to_enum_map["recomb.shr"  ] = RECOMB_SHR;
    solution_string_to_enum_map["recomb.auger"] = RECOMB_AUGER;

    solution_string_to_enum_map["mole.x"      ] = MOLE_X;
    solution_string_to_enum_map["mole.y"      ] = MOLE_Y;

    solution_string_to_enum_map["ii.gen"      ] = II_GEN;
    solution_string_to_enum_map["optical.gen" ] = OPTICAL_GEN;
    solution_string_to_enum_map["particle.gen"] = PARTICLE_GEN;
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



