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

#ifndef __suprem_h__
#define __suprem_h__

#include <map>

#include "stanford.h"

class SupremTIF : public StanfordTIF
{
public:
  /**
   * constructor
   */
  SupremTIF(const std::string & file);

  /**
   * free
   */
  virtual ~SupremTIF()  {}


  /** read suprem file into meta data structure*/
  bool read();

private:

  std::vector<int>       _electrode_info;

  /**
   * mapping from the material int index to material name used in suprem4-gs
   */
  std::map<int, std::string> _material_index_to_string;

  /**
   * mapping from the material name to material int index used in suprem4-gs
   */
  std::map<std::string, int> _material_string_to_index;

  /**
   * mapping from the doping int index to doping name used in suprem4-gs
   */
  std::map<int, std::string> _doping_index_to_string;

  /**
   * mapping from the doping name to doping int index used in suprem4-gs
   */
  std::map<std::string, int> _doping_string_to_index;



  void _init_index_string_map();

};


#endif

