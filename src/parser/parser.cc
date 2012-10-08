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

//  $Id: parser.cc,v 1.15 2008/07/09 05:58:16 gdiso Exp $

#include <cstdio>
#include <cstring>
#include "config.h"
#include "parser.h"
#include "log.h"


namespace Parser
{

int InputParser::read_card_file(const char *filename)
{
  InputYY::yyin = fopen(filename, "r");

  if( InputYY::yyin == NULL )
    return 1;

  if( InputYY::yyparse((void*)this) )
    return 1;

  fclose(InputYY::yyin);

  std::list<Card>::iterator it;
  for( it=_card_list.begin(); it!=_card_list.end(); it++ )
  {
    _pattern.check_detail(*it);
  }
  return 0;
}

}

// avoid isatty() problem of Bison 2.3
#define YY_NEVER_INTERACTIVE 1
#ifdef WINDOWS
  #define YY_NO_UNISTD_H 1
extern "C"
{
extern int isatty (int );
}
#endif

namespace InputYY
{
#include "input.yy.c"
#include "input.tab.c"
}

