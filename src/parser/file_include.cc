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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


#include "genius_common.h"
#include "genius_env.h"
#include "file_include.h"

using namespace Parser;

FilePreProcess::FilePreProcess(const char * filename): _filename(filename)
{
  file_include( filename );
}


/**
 * write pre-processed file into a temporary file with extersion .pp
 * do following things:
 *   expand include file
 *   remove comment
 *   add file:line info to each line
 */
void FilePreProcess::file_include( const char * filename  )
{
  // open the file
  std::ifstream in(filename);
  if(!in.good())
  {
     std::cerr << "ERROR: file " << filename <<" doesn't exist." << std::endl;
     genius_error();
  }

  std::string str;
  unsigned int line_counter=0;

  while(!in.eof())
  {
    std::getline(in, str);
    ++line_counter;

    // the include statement is at the beginnning of string
    if( str.find(".include")==0 || str.find(".INCLUDE")==0 )
    {
      std::string include, include_file;
      std::stringstream   ss;
      ss << str;
      ss >> include >> include_file;

      // drop "" or <>
      if( include_file.find("\"")==0 && include_file.rfind("\"")==include_file.length()-1 )
         include_file = include_file.substr(1, include_file.length()-2);
      else if( include_file.find("<")==0 && include_file.rfind(">")==include_file.length()-1 )
         include_file = include_file.substr(1, include_file.length()-2);

      file_include(include_file.c_str());
    }
    else
    {
      // write origin file/line info into pre-processed file
      std::stringstream   ss;
      ss << "!\""<<filename<<':'<<line_counter<<"\"";
      std::string file_line_info(ss.str());

      //remove comment begin with '#'
      while(str.rfind("#") < str.length() )
        str = str.substr(0, str.rfind("#"));

      //remove comment begin with '//'
      while(str.rfind("//") < str.length() )
        str = str.substr(0, str.rfind("//"));

      // process continue line
      if(str.rfind("\\") < str.length() )
        _contex += str.substr(0, str.rfind("\\")) + ' ';
      else
        _contex += str + ' ' + file_line_info + '\n';
    }
  }

  in.close();
}


std::string FilePreProcess::output()
{
  // write down
  std::string out_file = _filename + ".pp";
  std::ofstream   out( out_file.c_str() );
  out << _contex;
  out.close();

  return out_file;
}

