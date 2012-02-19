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


#include <sstream>

#include "parser_card.h"

namespace Parser
{

  inline void _to_upper(std::string & str)
  {
    for(unsigned int c=0; c<str.length(); c++)
    {
      if(islower(str[c]))
        str[c]=toupper(str[c]);
    }
  }


  inline void _to_lower(std::string & str)
  {
    for(unsigned int c=0; c<str.length(); c++)
    {
      if(isupper(str[c]))
        str[c]=tolower(str[c]);
    }
  }



  //--------------------------------------------------------------------------------------

  Card::Card(const std::string &description)
  {
    // check if first char is '#'
    {
      char first_char;
      std::stringstream ss(description);
      ss >> first_char;
      if(first_char == '#') return;
    }


    std::string str(description);
    // change <>= chars to blank
    {
      std::string filt_elems("<>=");
      std::string::size_type pos = 0;
      while (( pos = str.find_first_of( filt_elems, pos )) != std::string::npos )
        str.replace(pos, 1, " ");
    }


    std::stringstream ss(str);

    // read key
    ss >> _key;
    _to_upper(_key);

    // read parameters
    while(1)
    {
      std::string type, parameter;
      ss >> type;
      if(ss.fail()) break;
      _to_lower(type);

      ss >> parameter;
      _to_lower(parameter);

      if(type == "string")
      {
        std::string string_value;
        ss >> string_value;
        insert(parameter, string_value);
      }

      if(type == "enum")
      {
        std::string string_value;
        ss >> string_value;
        _to_lower(string_value);
        insert(parameter, string_value);
      }

      if(type == "real")
      {
        double real_value;
        ss >> real_value;
        insert(parameter, real_value);
      }

      if(type == "int")
      {
        int int_value;
        ss >> int_value;
        insert(parameter, int_value);
      }

      if(type == "bool")
      {
        std::string string_value;
        ss >> string_value;
        _to_lower(string_value);
        if(string_value=="true" || string_value=="on")
          insert(parameter, true);
        if(string_value=="false" || string_value=="off")
          insert(parameter, false);
      }
    }
  }



  template <typename T>
  std::vector<T>  Card::get_n_array(const char * name, unsigned int idx, unsigned int alias_num, ...)   const
  {
    std::string parameter_name = name;
    std::vector<std::string> alias;
    ConstMapIt it;

    // use va to get alias name and defalut value
    va_list arg;

    // va start from alias_num
    va_start(arg,alias_num);

    // the number of alias should equal to alias_num
    alias.resize(alias_num);
    for(unsigned int i=0;i<alias_num;i++)
      alias[i] = va_arg(arg, const char *);

    // ok
    va_end(arg);

    // if the parameter_name can be find in the parameter multimap?
    it  = _parameter_map.find(parameter_name);

    //parameter_name not match? try alias
    if( it == _parameter_map.end() )
    {
      for(unsigned int i=0;i<alias_num;i++)
        if ( (it = _parameter_map.find(alias[i])) != _parameter_map.end() )
        {
          // format to primary name
          parameter_name = alias[i];
          break;
        }
    }

    // still no find, return empty array
    if ( it == _parameter_map.end() )
      return std::vector<T>();

    // return the idx'th parameter vaule
    it  = _parameter_map.equal_range(parameter_name).first;
    size_t size = _parameter_map.count(name);
    assert( size-1 >= idx );

    for( unsigned int i=0; i<idx; i++,it++ );
    return it->second.get_array<T>();
  }


  template <typename T>
  std::vector<T>  Card::get_array(const char * name, const char * alias)   const
  {
    if( alias == NULL )
      return get_n_array<T>(name, 0, 0);
    return get_n_array<T>(name, 0, 1, alias);
  }


  template <typename T>
  T  Card::get_n(const char * name, T default_value, unsigned int idx,  unsigned int alias_num, ...)   const
  {
    std::string parameter_name = name;
    std::vector<std::string> alias;
    ConstMapIt it;

    // use va to get alias name and defalut value
    va_list arg;

    // va start from alias_num
    va_start(arg,alias_num);

    // the number of alias should equal to alias_num
    alias.resize(alias_num);
    for(unsigned int i=0;i<alias_num;i++)
      alias[i] = va_arg(arg, const char *);

    // ok
    va_end(arg);

    // if the parameter_name can be find in the parameter multimap?
    it  = _parameter_map.find(parameter_name);

    //parameter_name not match? try alias
    if( it == _parameter_map.end() )
    {
      for(unsigned int i=0;i<alias_num;i++)
        if ( (it = _parameter_map.find(alias[i])) != _parameter_map.end() )
        {
          // format to primary name
          parameter_name = alias[i];
          break;
        }
    }

    // still no find, return empty array
    if ( it == _parameter_map.end() )
      return default_value;

    // return the idx'th parameter vaule
    it  = _parameter_map.equal_range(parameter_name).first;
    size_t size = _parameter_map.count(name);
    assert( size-1 >= idx );

    for( unsigned int i=0; i<idx; i++,it++ );
    return it->second.get<T>();
  }


  template <typename T>
  T  Card::get(const char * name, T default_value, const char * alias)   const
    {
      if( alias == NULL )
        return get_n<T>(name , default_value, 0 , 0);
      return get_n<T>(name , default_value, 0 , 1, alias);
    }

  /*-------------------------------------------------------------------------------
   * since parameter with the same name is allowed in miltumap,
   * this function return nth bool value which matches parameter name
   * alias is also supported
   */
  bool    Card::get_n_bool(const char * name , bool default_value, unsigned int idx, unsigned int alias_num, ...)   const
  {
    std::string parameter_name = name;
    std::vector<std::string> alias;
    ConstMapIt it;

    // use va to get alias name and defalut value
    va_list arg;

    // va start from alias_num
    va_start(arg,alias_num);

    // the number of alias should equal to alias_num
    alias.resize(alias_num);
    for(unsigned int i=0;i<alias_num;i++)
      alias[i] = va_arg(arg, const char *);

    // ok
    va_end(arg);

    // if the parameter_name can be find in the parameter multimap?
    it  = _parameter_map.find(parameter_name);

    //parameter_name not match? try alias
    if( it == _parameter_map.end() )
    {
      for(unsigned int i=0;i<alias_num;i++)
        if ( (it = _parameter_map.find(alias[i])) != _parameter_map.end() )
        {
          // format to primary name
          parameter_name = alias[i];
          break;
        }
    }

    // still no find, return default value
    if ( it == _parameter_map.end() )
      return default_value;

    // return the idx'th parameter vaule
    it  = _parameter_map.equal_range(parameter_name).first;
    size_t size = _parameter_map.count(name);
    assert( size-1 >= idx );

    for( unsigned int i=0; i<idx; i++,it++ );

    assert(it->second.type() == BOOL);
    return it->second.get_bool();
  }



  /*-------------------------------------------------------------------------------
   * since parameter with the same name is allowed in miltumap,
   * this function return nth integer value which matches parameter name
   * alias is also supported
   */
  int     Card::get_n_int(const char * name , int default_value, unsigned int idx, unsigned int alias_num, ...)   const
  {
    std::string parameter_name = name;
    std::vector<std::string> alias;
    ConstMapIt it;

    // use va to get alias name and defalut value
    va_list arg;

    // va start from alias_num
    va_start(arg,alias_num);

    // the number of alias should equal to alias_num
    alias.resize(alias_num);
    for(unsigned int i=0;i<alias_num;i++)
      alias[i] = va_arg(arg, const char *);

    // ok
    va_end(arg);

    // if the parameter_name can be find in the parameter multimap?
    it  = _parameter_map.find(parameter_name);

    //parameter_name not match? try alias
    if( it == _parameter_map.end() )
    {
      for(unsigned int i=0;i<alias_num;i++)
        if ( (it = _parameter_map.find(alias[i])) != _parameter_map.end() )
        {
          // format to primary name
          parameter_name = alias[i];
          break;
        }
    }

    // still no find, return default value
    if ( it == _parameter_map.end() )
      return default_value;

    // return the idx'th parameter vaule
    it  = _parameter_map.equal_range(parameter_name).first;
    size_t size = _parameter_map.count(name);
    assert( size-1 >= idx );

    for( unsigned int i=0; i<idx; i++,it++ );

    assert(it->second.type() == INTEGER);
    return it->second.get_int();
  }



  /*-------------------------------------------------------------------------------
   * since parameter with the same name is allowed in miltumap,
   * this function return nth real value which matches parameter name
   * alias is also supported
   */
  double  Card::get_n_real(const char * name , double default_value, unsigned int idx, unsigned int alias_num, ...)   const
  {
    std::string parameter_name = name;
    std::vector<std::string> alias;
    ConstMapIt it;

    // use va to get alias name and defalut value
    va_list arg;

    // va start from alias_num
    va_start(arg,alias_num);

    // the number of alias should equal to alias_num
    alias.resize(alias_num);
    for(unsigned int i=0;i<alias_num;i++)
      alias[i] = va_arg(arg, const char *);

    // ok
    va_end(arg);

    // if the parameter_name can be find in the parameter multimap?
    it  = _parameter_map.find(parameter_name);


    //parameter_name not match? try alias
    if( it == _parameter_map.end() )
    {
      for(unsigned int i=0;i<alias_num;i++)
        if ( (it = _parameter_map.find(alias[i])) != _parameter_map.end() )
        {
          // format to primary name
          parameter_name = alias[i];
          break;
        }
    }


    // still no find, return default value
    if ( it == _parameter_map.end() )
      return default_value;

    // return the idx'th parameter vaule
    it  = _parameter_map.equal_range(parameter_name).first;
    size_t size = _parameter_map.count(name);
    assert( size-1 >= idx );

    for( unsigned int i=0; i<idx; i++,it++ );

    assert(it->second.type() == REAL);
    return it->second.get_real();
  }



  /*-------------------------------------------------------------------------------
   * since parameter with the same name is allowed in miltumap,
   * this function return nth string value which matches parameter name
   * alias is also supported
   */
  std::string  Card::get_n_string(const char * name , std::string default_value, unsigned int idx, unsigned int alias_num, ...)   const
  {
    std::string parameter_name = name;
    std::vector<std::string> alias;
    ConstMapIt it;

    // use va to get alias name and defalut value
    va_list arg;

    // va start from alias_num
    va_start(arg,alias_num);

    // the number of alias should equal to alias_num
    alias.resize(alias_num);
    for(unsigned int i=0;i<alias_num;i++)
      alias[i] = va_arg(arg, const char *);

    // ok
    va_end(arg);

    // if the parameter_name can be find in the parameter multimap?
    it  = _parameter_map.find(parameter_name);

    //parameter_name not match? try alias
    if( it == _parameter_map.end() )
    {
      for(unsigned int i=0;i<alias_num;i++)
        if ( (it = _parameter_map.find(alias[i])) != _parameter_map.end() )
        {
          // format to primary name
          parameter_name = alias[i];
          break;
        }
    }

    // still no find, return default value
    if ( it == _parameter_map.end() )
      return default_value;

    // return the idx'th parameter vaule
    it  = _parameter_map.equal_range(parameter_name).first;
    size_t size = _parameter_map.count(name);
    assert( size-1 >= idx );

    for( unsigned int i=0; i<idx; i++,it++ );

    assert(it->second.type() == STRING || it->second.type() == ENUM);
    return it->second.get_string();
  }




  /*-------------------------------------------------------------------------------
   * since parameter with the same name is allowed in miltumap,
   * this function return nth string value which matches parameter name
   * alias is also supported
   */
  std::string  Card::get_n_string_lower_case(const char * name , std::string default_value, unsigned int idx, unsigned int alias_num, ...)   const
  {
    std::string str;
    std::string parameter_name = name;
    std::vector<std::string> alias;
    ConstMapIt it;

    // use va to get alias name and defalut value
    va_list arg;

    // va start from alias_num
    va_start(arg,alias_num);

    // the number of alias should equal to alias_num
    alias.resize(alias_num);
    for(unsigned int i=0;i<alias_num;i++)
      alias[i] = va_arg(arg, const char *);

    // ok
    va_end(arg);

    // if the parameter_name can be find in the parameter multimap?
    it  = _parameter_map.find(parameter_name);

    //parameter_name not match? try alias
    if( it == _parameter_map.end() )
    {
      for(unsigned int i=0;i<alias_num;i++)
        if ( (it = _parameter_map.find(alias[i])) != _parameter_map.end() )
        {
          // format to primary name
          parameter_name = alias[i];
          break;
        }
    }

    // still no find, return default value
    if ( it == _parameter_map.end() )
      str = default_value;

    // return the idx'th parameter vaule
    it  = _parameter_map.equal_range(parameter_name).first;
    size_t size = _parameter_map.count(name);
    assert( size-1 >= idx );

    for( unsigned int i=0; i<idx; i++,it++ );

    assert(it->second.type() == STRING || it->second.type() == ENUM);
    str = it->second.get_string();

    //convert to lower case
    std::string lower_str;
    for(unsigned int n=0; n<str.size(); ++n)
    {
      char c =  str.at(n);
      if( isalpha(c) )  c = tolower(c);
      lower_str += c;
    }

    return lower_str;

  }





  /*-------------------------------------------------------------------------------
   * @return first bool value matches name, if no match found, default value will be return.
   * only one alias is allowed.
   */
  bool    Card::get_bool(const char * name , bool default_value, const char * alias)   const
  {
    if( alias == NULL )
      return get_n_bool(name , default_value, 0 , 0);
    return get_n_bool(name , default_value, 0 , 1, alias);
  }


  /*-------------------------------------------------------------------------------
   * @return first int value matches name, if no match found, default value will be return
   */
  int     Card::get_int(const char * name , int default_value, const char * alias)   const
  {
    if( alias == NULL )
      return get_n_int(name , default_value, 0 , 0);
    return get_n_int(name , default_value, 0 , 1, alias);
  }


  /*-------------------------------------------------------------------------------
   * @return first int real matches name, if no match found, default value will be return
   */
  double  Card::get_real(const char * name , double default_value, const char * alias)   const
  {
    if( alias == NULL )
      return get_n_real(name , default_value, 0 , 0);
    return get_n_real(name , default_value, 0 , 1, alias);
  }

  /*-------------------------------------------------------------------------------
   * @return first int string matches name, if no match found, default value will be return
   */
  std::string  Card::get_string(const char * name , std::string default_value, const char * alias)   const
  {
    if( alias == NULL )
      return get_n_string(name , default_value, 0 , 0);
    return get_n_string(name , default_value, 0 , 1, alias);
  }

  /*-------------------------------------------------------------------------------
   * @return first int string matches name, if no match found, default value will be return
   */
  std::string  Card::get_string_lower_case(const char * name , std::string default_value, const char * alias)   const
  {
    if( alias == NULL )
      return get_n_string_lower_case(name , default_value, 0 , 0);
    return get_n_string_lower_case(name , default_value, 0 , 1, alias);
  }

  //explicit instantiation

  template bool  Card::get(const char * name, bool default_value, const char * alias)   const;
  template int  Card::get(const char * name, int default_value, const char * alias)   const;
  template double  Card::get(const char * name, double default_value, const char * alias)   const;
  template std::string  Card::get(const char * name, std::string default_value, const char * alias)   const;
  template std::vector<bool>  Card::get_array(const char * name, const char * alias)   const;
  template std::vector<int>  Card::get_array(const char * name, const char * alias)   const;
  template std::vector<double>  Card::get_array(const char * name, const char * alias)   const;
  template std::vector<std::string>  Card::get_array(const char * name, const char * alias)   const;
}



