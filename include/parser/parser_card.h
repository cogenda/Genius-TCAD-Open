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

# ifndef  __parser_card_h__
# define  __parser_card_h__

#include <vector>
#include <set>
#include <map>
#include <string>
#include <iostream>


#include <cassert>
#include <cstdarg>


#include "genius_common.h"
#include "parser_parameter.h"

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

  /**
   * The card of input deck, contains a keyword and several parameters
   */
  class Card
  {
  public:
    typedef std::multimap< std::string, Parameter>::const_iterator ConstMapIt;
    typedef std::multimap< std::string, Parameter>::iterator       MapIt;

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
      _parameter_map.insert(std::pair< std::string, Parameter >(p.name(),p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert bool value parameter to  _parameter_map
     */
    void insert(const std::string &name, const bool value)
    {
      Parameter p(name, value);
      _parameter_map.insert(std::pair< std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert integer value parameter to  _parameter_map
     */
    void insert(const std::string &name, const int value)
    {
      Parameter p(name, value);
      _parameter_map.insert(std::pair< std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert float value parameter to  _parameter_map
     */
    void insert(const std::string &name, const double value)
    {
      Parameter p(name, value);
      _parameter_map.insert(std::pair< std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert string value parameter to  _parameter_map
     */
    void insert(const std::string &name, const std::string & value)
    {
      Parameter p(name, value);
      _parameter_map.insert(std::pair< std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }

    /**
     * insert const char * value (changed to string) parameter to  _parameter_map
     */
    void insert(const std::string &name, const char * value)
    {
      Parameter p(name, std::string(value));
      _parameter_map.insert(std::pair< std::string, Parameter >(name,p));
      _parameter_vec.push_back(p);
    }


    /**
     * since parameter with the same name is allowed in multimap,
     * this function return nth T value which matches parameter name
     * alias is also supported
     */
    template <typename T>
    T  get_n(const char * name, T default_value, unsigned int idx, unsigned int alias_num, ...)   const;

    /**
     * @return first T value matches name
     */
    template <typename T>
    T  get(const char * name, T default_value, const char * alias=NULL)   const;

    /**
     * since parameter with the same name is allowed in multimap,
     * this function return nth T array which matches parameter name
     * alias is also supported
     */
    template <typename T>
    std::vector<T>  get_n_array(const char * name, unsigned int idx, unsigned int alias_num, ...)   const;

    /**
     * @return first T array matches name
     */
    template <typename T>
    std::vector<T>  get_array(const char * name, const char * alias=NULL)   const;



    /**
     * since parameter with the same name is allowed in multimap,
     * this function return nth bool value which matches parameter name
     * alias is also supported
     */
    bool    get_n_bool(const char * name , bool default_value, unsigned int idx, unsigned int alias_num, ...)   const;


    /**
     * since parameter with the same name is allowed in multimap,
     * this function return nth integer value which matches parameter name
     * alias is also supported
     */
    int     get_n_int(const char * name , int default_value, unsigned int idx, unsigned int alias_num, ...)   const;


    /**
     * since parameter with the same name is allowed in multimap,
     * this function return nth real value which matches parameter name
     * alias is also supported
     */
    double  get_n_real(const char * name , double default_value, unsigned int idx, unsigned int alias_num, ...)   const;


    /**
     * since parameter with the same name is allowed in multimap,
     * this function return nth string value which matches parameter name
     * alias is also supported
     */
    std::string  get_n_string(const char * name , std::string default_value, unsigned int idx, unsigned int alias_num, ...)   const;

    /**
     * since parameter with the same name is allowed in multimap,
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
     * @return true if parameter exist.
     */
    bool  is_parameter_exist(const std::string & parameter_name, const std::string &alias) const
    {
      bool find_p     = _parameter_map.find(parameter_name) != _parameter_map.end() ;
      bool find_alias = _parameter_map.find(alias) != _parameter_map.end() ;
      return (find_p || find_alias);
    }


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
      return value == get_string_lower_case(name.c_str() , "");
    }

    /**
     * get keyword
     */
    const std::string & key()  const
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
    std::multimap< std::string, Parameter> _parameter_map;

    /**
     * store parameters in vec
     */
    std::vector<Parameter> _parameter_vec;

  };






}
# endif
