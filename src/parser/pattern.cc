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

//  $Id: pattern.cc,v 1.10 2008/07/09 05:58:16 gdiso Exp $

#include "genius_env.h"
#include "pattern.h"

#include <algorithm>

namespace Parser
{

  void to_lower(char & x)
  {
    if(isupper(x) )
      x = tolower(x);
  }

  int Pattern::get_pattern(const std::string & filename)
  {
    PatternYY::yyin = fopen(filename.c_str(), "r");

    if( PatternYY::yyin == NULL )
      return 1;

    if( PatternYY::yyparse((void*)this) )
      return 1;

    fclose(PatternYY::yyin);

    return 0;
  }

  int Pattern::ckeck_card (std::string & card_name) const
  {
    std::map< std::string, PatternCard>::const_iterator it;
    //we search card in the pattern map
    if ( _pattern_card_map.count(card_name) ) // exactly march
      return 0;
    else // ambiguous search
    {
      // search if card_name can match at the head of pattern_card
      std::map< std::string, PatternCard>::const_iterator tmp_it;
      unsigned int count = 0;
      for( tmp_it = _pattern_card_map.begin(); tmp_it != _pattern_card_map.end(); tmp_it++)
      {
        std::string pattern_card = (*tmp_it).first;
        if( pattern_card.find(card_name) == 0) { it = tmp_it; count++; }
      }
      // unique name find?
      if ( count == 1 )
      {
        card_name = (*it).first;
        return 0;
      }
    }
    return 1;
  }


  ElemType Pattern::check_parameter_type(const std::string &card_name, std::string &parameter_name) const
  {
    std::map< std::string, PatternCard>::const_iterator it;

    //we search card in the pattern map
    if ( _pattern_card_map.count(card_name) ) // exactly match
      it = _pattern_card_map.find(card_name);
    else // ambiguous search
    {
      // search if card_name can match at the head of pattern_card
      std::map< std::string, PatternCard>::const_iterator tmp_it;
      unsigned int count = 0;
      for( tmp_it = _pattern_card_map.begin(); tmp_it != _pattern_card_map.end(); tmp_it++)
      {
        std::string pattern_card = (*tmp_it).first;
        if( pattern_card.find(card_name) == 0) { it = tmp_it; count++; }
      }
      // unique name find?
      if ( count != 1 )
      {
        return INVALID;
      }
    }

    // get pattern parameter from pattern card
    Parameter q; //empty parameter
    if ( (*it).second.find_parameter(parameter_name,q) )
    {
        return INVALID;
    }
    parameter_name = q.name();
    return q.type();
  }


  int Pattern::check_detail(Parser::Card & c) const
  {
    std::string card_name = c.key();
    std::map< std::string, PatternCard>::const_iterator it;

    //we search card in the pattern map
    if ( _pattern_card_map.count(card_name) ) // exactly march
      it = _pattern_card_map.find(card_name);
    else // ambiguous search
    {
      // search if card_name can match at the head of pattern_card
      std::map< std::string, PatternCard>::const_iterator tmp_it;
      unsigned int count = 0;
      for( tmp_it = _pattern_card_map.begin(); tmp_it != _pattern_card_map.end(); tmp_it++)
      {
        std::string pattern_card = (*tmp_it).first;
        if( pattern_card.find(card_name) == 0) { it = tmp_it; count++; }
      }
      // unique name find?
      if ( count != 1 )
      {
        std::cerr << "ERROR: Card ["<<card_name<< "] at " << c.get_fileline() << " can't be matched in the pattern." << std::endl;
        genius_error();
      }
    }
    // set key to "official" name
    c.set_key( (*it).first );


    //check for all the parameter in vec
    for( unsigned int i=0; i<c.parameter_size(); i++ )
    {
      Parameter p = c.get_parameter(i);
      Parameter q; //empty

      // for user defined parameter, we should not check it
      if(p.do_not_check_me==true) continue;

      // get pattern parameter from pattern card
      if ( (*it).second.find_parameter(p.name(),q) )
      {
        std::cerr << "ERROR: Card ["<<card_name<< "] parameter ["<<p.name()
        << "] at " << c.get_fileline() << " can't be matched in the pattern." << std::endl;
        genius_error();
      }

      // set p to exact name
      p.set_name(q.name());


      // check the ENUM
      if ( q.string_pattern_size() )
      {
        std::string p_sval = p.get_string();

        // change p_sval to lower case
        for_each(p_sval.begin(),p_sval.end(),to_lower);

        //if the p_sval is abbreviated, set to full string
        if( q.string_pattern_match(p_sval) )
        {
          std::cerr << "ERROR: Card ["<<card_name<< "] parameter ["<<p.name()
          << "] at " << c.get_fileline() << " with wrong enum range." << std::endl;
          genius_error();
        }

        //restore (maybe modified) string
        p.set_string(p_sval);
        //set type to ENUM
        p.set_type(ENUM);
      }

      // check the type
      if( p.type() != q.type() )
      {
        std::cerr << "ERROR: Card ["<<card_name<< "] parameter ["<<p.name()
        << "] at "<< c.get_fileline() << " with wrong argument type." << std::endl;
        genius_error();
      }

      //restore p
      c.set_parameter(p,i);
    }

    //rebuild parameter map
    c.rebuild_parameter_map();


    //ok
    return 0;
  }

}



// avoid isatty() problem of Bison 2.3
#define YY_NEVER_INTERACTIVE 1
#ifdef CYGWIN
  #define YY_NO_UNISTD_H 1
extern "C"
{
extern int isatty (int );
}
#endif

namespace PatternYY
{
#include "pattern.yy.c"
#include "pattern.tab.c"
}

