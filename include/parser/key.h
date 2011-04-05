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

//  $Id: key.h,v 1.13 2008/07/09 05:58:16 gdiso Exp $

# ifndef  __key_h__
# define  __key_h__

#include <vector>
#include <set>
#include <map>
#include <string>
#include <iostream>


#include <cassert>
#include <cstdarg>


#include "genius_common.h"

namespace PatternYY
{
  extern int yyparse(void *);
  extern FILE* yyin;
}


namespace InputYY
{
  extern int yyparse(void *);
  extern int yylex(void *);
  extern FILE* yyin;
}



namespace Parser
{

  enum  ElemType {BOOL, INTEGER, REAL, STRING, ENUM, INVALID};


  /**
   * The parameter class, contains a special parameter (name, type, value)
   */
  class Parameter
  {
  private:
    /**
     * the parameter's name
     */
    std::string  _name;

    /**
     * a short description to the parameter
     */
    std::string _description;

    /**
     * the parameter's type, use enum ElemType to describe
     */
    ElemType     _et;

    /**
     * the parameter's value, stored in structure. for union not support string
     */
    struct Elem
    {
      bool        bval;
      int         ival;
      double      dval;
      std::string sval;
    }
    _elem;

    /**
     * pattern for string value
     */
    std::set<std::string> _string_pattern;

  public:
    // constructors
  Parameter():do_not_check_me(false) {}

    Parameter(const std::string &name, const bool   v): _name(name), _et(BOOL),    do_not_check_me(false)      { _elem.bval = v; }
    Parameter(const std::string &name, const int    v): _name(name), _et(INTEGER), do_not_check_me(false)      { _elem.ival = v; }
    Parameter(const std::string &name, const double v): _name(name), _et(REAL),    do_not_check_me(false)      { _elem.dval = v; }
    Parameter(const std::string &name, const std::string & v): _name(name), _et(STRING),do_not_check_me(false) { _elem.sval = v; }

    /**
     * the const char * constructor is used for define as Parameter p("aaa")
     * without this constructor, the Parameter(const std::string &name, const bool   v)
     * will be choosen by compiler!
     */
    Parameter(const std::string &name, const char * v): _name(name), _et(STRING) { _elem.sval = v; }


    /**
     * set parameter name
     */
    void set_name (const std::string &name) { _name=name; }

    /**
     * @return name
     */
    const std::string name() const {return _name;}

    /**
     * set parameter description
     */
    void set_description (const std::string &desc) { _description = desc; }

    /**
     * @return description
     */
    const std::string description() const { return _description; }

    /**
     * set parameter type
     */
    void set_type (const ElemType type) { _et=type; }

    /**
     * @return the data type of parameter
     */
    ElemType type() const {return _et;}

    /**
     * set bool value
     */
    void    set_bool(const bool v)   {_elem.bval = v; _et = BOOL;}

    /**
     * set integer value
     */
    void    set_int(const int v)    {_elem.ival = v; _et = INTEGER;}

    /**
     * set real value
     */
    void   set_real(const double v)   {_elem.dval = v; _et = REAL;}

    /**
     * set string value
     */
    void   set_string(const std::string & s)  { _elem.sval = s; _et = STRING; }

    /**
     * set ENUM value, we should check if s fit the string pattern
     */
    int    set_enum(const std::string & s)
    {
      // no pattern exist
      if ( _string_pattern.size() == 0 )
      { _elem.sval = s; _et = ENUM; return 0; }

      // s is in the pattern
      if ( _string_pattern.count(s) )
      { _elem.sval = s; _et = ENUM; return 0; }

      // return error
      return 1;
    }

    /**
     * set pattern for string value
     */
    void add_string_pattern (const std::string &pattern)
    {
      _string_pattern.insert(pattern);
    }

    /**
     * the size of patterns for string value
     */
    size_t string_pattern_size ()
    {
      return _string_pattern.size();
    }

    typedef std::set<std::string>::const_iterator StringEnumIterator;
    StringEnumIterator stringPatternBegin() const { return _string_pattern.begin(); }
    StringEnumIterator stringPatternEnd() const { return _string_pattern.end(); }

    /**
     * clear string patterns
     */
    void clear_string_pattern ()   { _string_pattern.clear(); }

    /**
     * test if s match _string_pattern
     * if s is fuzzy matched, set s to exact string
     */
    int string_pattern_match ( std::string &s )
    {
      std::set<std::string>::iterator it;
      //exactly match
      if ( _string_pattern.count(s) ) return 0;
      else
      {
        // search if string s can match at the head of _string_pattern
        std::set<std::string>::iterator tmp_it;
        unsigned int count = 0;
        for( tmp_it = _string_pattern.begin(); tmp_it != _string_pattern.end(); tmp_it++)
        {
          std::string pattern_s = *tmp_it;
          if( pattern_s.find(s) == 0) { it = tmp_it; s = pattern_s; count++; }
        }
        // unique?
        if ( count != 1 ) { return 1; }
      }
      //ok
      return 0;
    }

    /**
     * clear the Parameter
     */
    void clear ()
    {
      _name.clear();
      _description.clear();
      _string_pattern.clear();
      do_not_check_me = false;
    }

    /**
     * for user defined parameter, we should not check it with pattern
     * since it will not exist in the pattern
     */
    bool do_not_check_me;

  public:
    /**
     * get bool value
     */
    bool    get_bool()   const {return _elem.bval;}

    /**
     * get integer value
     */
    int     get_int()    const {return _elem.ival;}

    /**
     * get real value
     */
    double  get_real()   const {return _elem.dval;}

    /**
     * get string value
     */
    std::string  get_string() const {return _elem.sval;}

    /**
     * @return ture if this parameter is user defined
     */
    bool   is_user_defined()
    { return do_not_check_me; }

  };




  /**
   * The card of input deck, contains a keyword and several parameters
   */
  class Card
  {
  public:
    typedef std::multimap< const std::string, Parameter>::const_iterator ConstMapIt;
    typedef std::multimap< const std::string, Parameter>::iterator       MapIt;

    /**
     * Constructor
     */
    Card() {}

    /**
     * Constructor, build the Card from string
     */
    Card(const std::string &description);

    /**
     * insert parameter to  _parameter_map
     */
    void insert(const Parameter &p)
    {
      _parameter_map.insert(std::pair< const std::string, Parameter >(p.name(),p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert bool value parameter to  _parameter_map
     */
    void insert(const std::string &name, const bool value)
    {
      Parameter p(name, value);
      _parameter_map.insert(std::pair< const std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert integer value parameter to  _parameter_map
     */
    void insert(const std::string &name, const int value)
    {
      Parameter p(name, value);
      _parameter_map.insert(std::pair< const std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert float value parameter to  _parameter_map
     */
    void insert(const std::string &name, const double value)
    {
      Parameter p(name, value);
      _parameter_map.insert(std::pair< const std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert string value parameter to  _parameter_map
     */
    void insert(const std::string &name, const std::string & value)
    {
      Parameter p(name, value);
      _parameter_map.insert(std::pair< const std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert const char * value (changed to string) parameter to  _parameter_map
     */
    void insert(const std::string &name, const char * value)
    {
      Parameter p(name, std::string(value));
      _parameter_map.insert(std::pair< const std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }

    /**
     * since parameter with the same name is allowed in miltumap,
     * this function return nth bool value which matches parameter name
     * alias is also supported
     */
    bool    get_n_bool(const char * name , bool default_value, unsigned int idx, unsigned int alias_num, ...)   const;


    /**
     * since parameter with the same name is allowed in miltumap,
     * this function return nth integer value which matches parameter name
     * alias is also supported
     */
    int     get_n_int(const char * name , int default_value, unsigned int idx, unsigned int alias_num, ...)   const;


    /**
     * since parameter with the same name is allowed in miltumap,
     * this function return nth real value which matches parameter name
     * alias is also supported
     */
    double  get_n_real(const char * name , double default_value, unsigned int idx, unsigned int alias_num, ...)   const;


    /**
     * since parameter with the same name is allowed in miltumap,
     * this function return nth string value which matches parameter name
     * alias is also supported
     */
    std::string  get_n_string(const char * name , std::string default_value, unsigned int idx, unsigned int alias_num, ...)   const;

    /**
     * since parameter with the same name is allowed in miltumap,
     * this function return nth lower cased string value which matches parameter name
     * alias is also supported
     */
    std::string  get_n_string_lower_case(const char * name , std::string default_value, unsigned int idx, unsigned int alias_num, ...)   const;


    /**
     * @return first bool value matches name, if no match found, default value will be return.
     * only one alias is allowed.
     */
    bool    get_bool(const char * name , bool default_value, const char * alias=NULL)   const;


    /**
     * @return first int value matches name, if no match found, default value will be return
     */
    int     get_int(const char * name , int default_value, const char * alias=NULL)   const;


    /**
     * @return first real matches name, if no match found, default value will be return
     */
    double  get_real(const char * name , double default_value, const char * alias=NULL)   const;


    /**
     * @return first string matches name, if no match found, default value will be return
     */
    std::string  get_string(const char * name , std::string default_value, const char * alias=NULL)   const;

    /**
     * @return first lower cased string matches name, if no match found, default value will be return
     */
    std::string  get_string_lower_case(const char * name , std::string default_value, const char * alias=NULL)   const;


    /**
     * @return true if parameter exist.
     */
    bool  is_parameter_exist(const std::string & parameter_name) const
      { return _parameter_map.find(parameter_name) != _parameter_map.end() ; }


    /**
     * @return the number of parameters with name parameter_name.
     */
    size_t parameter_count(const std::string & parameter_name) const
      { return _parameter_map.count(parameter_name); }

    /**
     * @return the number of parameters
     */
    size_t parameter_size() const
    {
      genius_assert(_parameter_map.size()==_parameter_vec.size());
      return _parameter_map.size();
    }

    /**
     * @return iterator begin of _parameter_map
     */
    MapIt   parameter_begin()
    { return _parameter_map.begin(); }

    /**
     * @return iterator end of _parameter_map
     */
    MapIt   parameter_end()
    { return _parameter_map.end(); }

    /**
     * @return const iterator begin of _parameter_map
     */
    ConstMapIt   parameter_begin() const
      { return _parameter_map.begin(); }

    /**
     * @return const iterator end of _parameter_map
     */
    ConstMapIt   parameter_end() const
      { return _parameter_map.end(); }

    /**
     * get the parameter by index
     */
    const Parameter & get_parameter(unsigned int idx) const
    {
      return  _parameter_vec[idx];
    }

    /**
     * set the parameter by index
     */
    void set_parameter(const Parameter &p, unsigned int idx)
    {
      _parameter_vec[idx] = p;
    }

    /**
     * rebuild parameter map with vec
     */
    void rebuild_parameter_map()
    {
      _parameter_map.clear();
      for(unsigned int i=0; i<_parameter_vec.size(); ++i)
      {
        const Parameter & p = _parameter_vec[i];
        _parameter_map.insert(std::pair<const std::string, Parameter>(p.name(),p));
      }
    }

    /**
     * @return ture if enum name with value
     */
    bool is_enum_value( const std::string & name,  const std::string &value) const
    {
      return value == get_string(name.c_str() , "");
    }

    /**
     * get keyword
     */
    std::string key()  const
      {return _key;}

    /**
     * set keyword
     */
    void set_key(const std::string key)
    { _key = key; }


    /**
     * clear all the information
     */
    void clear()
    {
      _key.clear();
      _file_line_info.clear();
      _parameter_map.clear();
      _parameter_vec.clear();
    }

    /**
     * get the file line info of this card
     */
    const std::string & get_fileline() const
      { return _file_line_info; }

    /**
     * set the line number of this card
     */
    void set_fileline(const std::string & file_line)
    {  _file_line_info = file_line ; }


    /**
     * get the line number of this card
     */
    int get_lineno() const
      { return _line_number; }

    /**
     * set the line number of this card
     */
    void set_lineno(int line)
    {  _line_number = line ; }


    void output() const
    {
      std::cout<< key() <<std::endl;
      for( unsigned int i=0; i < parameter_size(); i++ )
      {
        Parameter p = get_parameter(i);
        std::cout<< p.name()<<" ";
        switch (p.type())
        {
        case BOOL    : std::cout<< "BOOL   "   << p.get_bool()   <<std::endl; break;
        case INTEGER : std::cout<< "INT    "   << p.get_int()    <<std::endl; break;
        case REAL    : std::cout<< "REAL   "   << p.get_real()   <<std::endl; break;
        case STRING  : std::cout<< "STRING "   << p.get_string() <<std::endl; break;
        case ENUM    : std::cout<< "ENUM   "   << p.get_string() <<std::endl; break;
        case INVALID : std::cout<< "INVALID"   <<std::endl; break;
        }
      }
    }

  private:
    /**
     * keyword
     */
    std::string _key;

    /**
     * the line number of this card in input file
     */
    int _line_number;

    /**
     * the file:line information of this card
     */
    std::string _file_line_info;

    /**
     * store parameters in multi-map,
     * since we may get several parameters with same name
     */
    std::multimap< const std::string, Parameter> _parameter_map;

    /**
     * store parameters in vec
     */
    std::vector<Parameter> _parameter_vec;

  };






}
# endif // #ifndef __key_h__
