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

//  $Id: sync_file.cc,v 1.4 2008/07/09 05:58:16 gdiso Exp $

#include "genius_common.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "parallel.h"

/**
 * transport text file to other processor,
 * return local name as filename.processor_id 
 */
const std::string sync_file(const char * filename)
{
#ifdef HAVE_MPI
  // buffer for file text
  std::string   str;
  
  // read only at processor 0
  if (Genius::processor_id() == 0)
  {
    std::ifstream   in(filename);
    // test if file exist
    genius_assert(in.good());
    // read file into string
    std::istreambuf_iterator<char>   beg(in),   end;
    str = std::string(beg,   end);
    
    in.close();
  }

  // Broadcast the string
  Parallel::broadcast(str);
  
  //generate local file name 
  std::string localfilename(filename);
  std::string processor;
  std::stringstream   ss;
  ss << Genius::processor_id();
  ss >> processor;
  localfilename = localfilename + "." + processor;
  
  // write down
  std::ofstream   out(localfilename.c_str());
  out << str;
  out.close();
  
  return localfilename;

#else

  // test if file exist
  std::ifstream   in(filename);
  genius_assert(in.good());
  
  // read file into string
  std::istreambuf_iterator<char>   beg(in),   end;
  std::string str(beg,   end);
  in.close();
  
  //generate local file name 
  std::string localfilename(filename);
  std::string processor;
  std::stringstream   ss;
  ss << Genius::processor_id();
  ss >> processor;
  localfilename = localfilename + "." + processor;
  
  // write down
  std::ofstream   out(localfilename.c_str());
  out << str;
  out.close();
  
  return localfilename;

#endif  

}

