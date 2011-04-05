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

//  $Id: pattern.h,v 1.7 2008/07/09 05:58:16 gdiso Exp $

#ifndef __pattern_h_
#define __pattern_h_

#include "key.h"

namespace Parser
{
  /**
   * The pattern of Card, read from file
   */
  struct PatternCard
  {
    /**
     * keyword
     */
    std::string _key;

    /**
     * description
     */
    std::string _description;

    /**
     * store parameters in map is enough
     */
    std::map< std::string, Parameter> _parameter_map;

    /**
     * find pattern parameter by fuzzy name
     */
    int  find_parameter(const std::string &name, Parameter & p) const
    {
      std::map< std::string, Parameter>::const_iterator it;

      //we search parameter name in the _parameter_map
      if ( _parameter_map.count(name) ) // exactly march
        it = _parameter_map.find(name);
      else // ambiguous search
      {
        // search if parameter_name can match at the head of pattern_parameter_name
        std::map< std::string, Parameter>::const_iterator tmp_it;
        unsigned int count = 0;
        for( tmp_it = _parameter_map.begin(); tmp_it != _parameter_map.end(); tmp_it++)
        {
          std::string pattern_parameter = (*tmp_it).first;
          if( pattern_parameter.find(name) == 0) { it = tmp_it; count++; }
        }

        if ( count != 1 ) { return 1; }
      }

      // get parameter by reference
      p = (*it).second;

      return 0;
    }

    /**
     * clear structure
     */
    void clear()
    {
      _key.clear();
      _description.clear();
      _parameter_map.clear();
    }

  };



  /**
   * The pattern of input deck, read from file by yyparse
   */
  class Pattern
  {
  private:
    /**
     * use map is enough here
     */
    std::map< std::string, PatternCard> _pattern_card_map;

    friend int PatternYY::yyparse(void *);

  public:

    Pattern() {}

    ~Pattern() {}

    /**
     * call yyparse to read pattern file (intput key)
     */
    int get_pattern(const std::string &);

    /**
     * when a card is read from user's input deck, check it with pattern
     * @return 0 if successful.
     *
     * Names (card's and parameter's) may be abbreviated by omitting characters from the
     * end, provided that the abbreviation is unambiguous. We should recognize them and translate to
     * formal names.
     *
     * The enum name is also checked by this rounting
     *
     * For avoiding case problem,
     * the input paser must change card name to upper case and parameter name to lower case!
     * Also enum strings will be converted to lower case.
     */
    int check_detail(Card & c) const;

    /**
     * check if card name matched in pattern, and return "offical" card name
     */
    int ckeck_card (std::string &card_name) const;

    /**
     * check the parameter of card
     * return its type, and return exact parameter name
     */
    ElemType check_parameter_type(const std::string &card_name, std::string &parameter_name) const;

    /**
     * save pattern description to xml file
     * @returns zero if successful
     */
    int save_to_XML(const std::string &fname) const;

    /**
     * read pattern description from xml file, meant to replace get-pattern()
     * @returns zero if successful
     */
    int get_from_XML(const std::string &);
  };

}

#endif // #define __pattern_h_
