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

#ifndef __variable_define_h__
#define __variable_define_h__

#include <string>

#include "enum_data_type.h"
#include "enum_data_location.h"


struct SimulationVariable
{
  SimulationVariable() {}

  SimulationVariable(const std::string &name, DataType data_type, DataLocation data_location,
                     const std::string & unit_string="", unsigned int index=static_cast<unsigned int>(-1),
                     bool valid=true, bool user_defined=false);

  std::string   variable_name;
  DataType      variable_data_type;
  DataLocation  variable_data_location;
  std::string   variable_unit_string;
  double        variable_unit;
  unsigned int  variable_index;
  bool          variable_valid;
  bool          variable_user_defined;
};



#endif
