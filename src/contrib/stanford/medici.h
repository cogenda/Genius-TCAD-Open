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

#ifndef __medici_tif_h__
#define __medici_tif_h__

#include "stanford.h"


class MediciTIF : public StanfordTIF
{
public:
  
  MediciTIF() {}

  MediciTIF(const std::string & file);

  virtual ~MediciTIF() {}

  /// read MediciTIF file
  virtual bool read(std::string &err);

  /**
   * export to MEDICI TIF file
   */ 
  void export_tif(const std::string & file) const;

private:

  std::string _version;

  /// map for solution name format
  std::map<std::string, std::string> _solution_name_map;

  std::string _solution_name_format( const std::string &) const;

};


#endif
