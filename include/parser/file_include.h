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

#ifndef _FILE_INCLUDE_H_
#define _FILE_INCLUDE_H_

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace Parser
{
/**
 * process the INCLUDE statement of input file
 * just simply expand the include contex into input file
 */
class FilePreProcess
{
private:

  std::string  _filename;

  /**
   * string buffer
   */
  std::string     _contex;

public:

  FilePreProcess(const char * filename);


  ~FilePreProcess() {}

  /**
   * read contex from input file(s)
   */
  void file_include( const char * filename );

  /**
   * dump file processed by preprocessor
   */
  std::string output();
};

}
#endif
