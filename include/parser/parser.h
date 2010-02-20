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

//  $Id: parser.h,v 1.8 2008/07/09 05:58:16 gdiso Exp $

#ifndef __parser_h_
#define __parser_h_

// C++ header file
#include <stack>
#include <list>

#include "key.h"
#include "pattern.h"


namespace Parser
{

  /**
   * user's card intepreter
   */
  class InputParser
  {
  public:
    InputParser(Pattern &p) : _pattern(p) {}

    /**
     * read card from user's file
     */
    int read_card_file(const char *filename)
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

    /**
     * card search begin
     */
    void      begin()
    {
      it_card_stack.push(it_card);
      it_card = _card_list.begin();
    }

    /**
     * move to next card
     */
    void      next ()
    { it_card++; }

    /**
     * card search end
     */
    bool      end()
    {
      if( it_card == _card_list.end() )
      {
        it_card = it_card_stack.top();
        it_card_stack.pop();
        return true;
      }
      else
        return false;
    }

    /**
     * safe break the iteration
     */
    void      breakloop ()
    {
        it_card = it_card_stack.top();
        it_card_stack.pop();
    }

    /**
     * return const Card reference
     */
    const Card &  get_current_card()    const
    { return *it_card; }

    /**
     * return writable Card reference
     */
    Card &  get_current_card()
    { return *it_card; }

    /**
     * return true if current card matches card_name
     */
    bool      is_current_card(const std::string & card_name) const
    { return it_card->key() == card_name; }

    /**
     * delete current card, the iterator point to next card
     */
    void      delete_current_card()
    { it_card = _card_list.erase(it_card); }

    /**
     * output card information
     */
    void output()
    {

      std::list<Card>::iterator it;

      for( it=_card_list.begin(); it!=_card_list.end(); it++ )
      {
        it->output();
        std::cout<<std::endl;
      }
    }

    /**
     * @return true if card with card_name can be find, otherwise false
     */
    bool  is_card_exist(const std::string & card_name)
    {
      if ( _card_map.find(card_name) != _card_map.end() )
        return true;
      return false;
    }


  private:
    /**
     * when parser finished, the input cards are stored here as list
     */
    std::list<Card> _card_list;

    /**
     * _card_list iterator
     */
    std::list<Card>::iterator    it_card;

    /**
     * when parser finished, the input cards are stored here as map for some special purpose
     */
    std::multimap<const std::string, Card> _card_map;

    /**
     * the const reference to syntax pattern structure
     */
    const Pattern & _pattern;

    /**
     * map for bool variable defined in input file
     */
    std::map<const std::string, bool>        bool_var;

    /**
     * map for int variable defined in input file
     */
    std::map<const std::string, int>         int_var;

    /**
     * map for double variable defined in input file
     */
    std::map<const std::string, double>      real_var;

    /**
     * map for string variable defined in input file
     */
    std::map<const std::string, std::string> string_var;


    /**
     * A stack for indicating state of searching the card_list
     */
    std::stack<std::list<Card>::iterator> it_card_stack;


    /**
     * let yyparse be the friend of this class, then it can directly access  _card_list and _card_map
     */
    friend int InputYY::yyparse(void *);

    /**
     * let yylex be the friend of this class, then it can directly access  _card_list and _card_map
     */
    friend int InputYY::yylex(void *);


  };

}

#endif // #define __parser_h_

