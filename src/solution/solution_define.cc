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
#include "enum_solution.h"




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
      solution_string_to_enum_map["doping"      ] = DOPING;
      solution_string_to_enum_map["doping.na"   ] = DOPING_Na;
      solution_string_to_enum_map["doping.nd"   ] = DOPING_Nd;
      solution_string_to_enum_map["min.carrier" ] = MIN_CARRIER;
      solution_string_to_enum_map["net.carrier" ] = NET_CARRIER;
      solution_string_to_enum_map["net.charge"  ] = NET_CHARGE;
      solution_string_to_enum_map["optical.gen" ] = OPTICAL_GEN;
      solution_string_to_enum_map["particle.gen"] = PARTICLE_GEN;
      solution_string_to_enum_map["mole.x"      ] = MOLE_X;
      solution_string_to_enum_map["mole.y"      ] = MOLE_Y;
      solution_string_to_enum_map["qfn"         ] = QFN;
      solution_string_to_enum_map["qfp"         ] = QFP;
 }

}


SolutionVariable solution_string_to_enum(const std::string & str)
{
  init_solution_string_to_enum_map();

  if( solution_string_to_enum_map.find(str)!= solution_string_to_enum_map.end() )
    return solution_string_to_enum_map[str];
  return INVALID_Variable;
}


