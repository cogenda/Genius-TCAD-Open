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



#include "parser_parameter.h"

namespace Parser
{

  template <typename T>
  T& Parameter::set (unsigned int idx)
  {
    Value<T>* value = 0;
    if(_values.size() < idx+1)
    {
      _values.resize(idx+1, 0);
      _values[idx] = new Value<T>;
    }
    else
    {
      if(!_values[idx])
        _values[idx] = new Value<T>;
    }
    value = dynamic_cast<Value<T>*>(_values[idx]);
    return  value->set();
  }


  void    Parameter::set_bool(const bool v, unsigned int idx)
{ this->set<bool>(idx)=v; _et = BOOL; }


  void    Parameter::set_int(const int v, unsigned int idx)
  { this->set<int>(idx)=v; _et = INTEGER;}


  void   Parameter::set_real(const double v, unsigned int idx)
  {this->set<double>(idx)=v; _et = REAL;}


  void   Parameter::set_string(const std::string & s, unsigned int idx)
  { this->set<std::string>(idx)=s; _et = STRING; }


  int    Parameter::set_enum(const std::string & s, unsigned int idx)
  {
    // no pattern exist
    if ( _string_pattern.size() == 0 )
    { this->set<std::string>(idx)=s; _et = ENUM; return 0; }

    // s is in the pattern
    if ( _string_pattern.count(s) )
    { this->set<std::string>(idx)=s; _et = ENUM; return 0; }

    // return error
    return 1;
  }


  void    Parameter::set_bool_array(const std::vector<bool> &v)
  {
    for(unsigned int n=0; n<v.size(); ++n)
      this->set<bool>(n)=v[n];
    _et = BOOL;
  }


  void    Parameter::set_int_array(const std::vector<int> & v)
  {
    for(unsigned int n=0; n<v.size(); ++n)
      this->set<int>(n)=v[n];
    _et = INTEGER;
  }


  void   Parameter::set_real_array(const std::vector<double> & v)
  {
    for(unsigned int n=0; n<v.size(); ++n)
      this->set<double>(n)=v[n];
    _et = REAL;
  }


  void   Parameter::set_string_array(const std::vector<std::string> & s)
  {
    for(unsigned int n=0; n<s.size(); ++n)
      this->set<std::string>(n)=s[n];
    _et = STRING;
  }

  int Parameter::string_pattern_match ( std::string &s )
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


  unsigned int Parameter::array_size() const
  { return _values.size(); }


  template <typename T>
  const T & Parameter::get (unsigned int idx) const
    {
      Value<T>* value = dynamic_cast<Value<T>*>(_values[idx]);
      assert(value);
      return  value->get();
    }



  bool    Parameter::get_bool(unsigned int idx)   const
    {return this->get<bool>(idx);}


  int     Parameter::get_int(unsigned int idx)    const
    {return this->get<int>(idx);}


  double  Parameter::get_real(unsigned int idx)   const
    {return this->get<double>(idx);}


  const std::string & Parameter::get_string(unsigned int idx) const
    {return this->get<std::string>(idx);}


  template <typename T>
  std::vector<T> Parameter::get_array() const
  {
    std::vector<T>  array;
    for(unsigned int n=0; n<_values.size(); ++n)
      array.push_back(this->get<T>(n));
    return array;
  }


  void Parameter::clear ()
  {
    _name.clear();
    _description.clear();
    _string_pattern.clear();
    for(unsigned int n=0; n<_values.size(); ++n)
      delete _values[n];
    _values.clear();
    _do_not_check_me = false;
  }


  Parameter & Parameter::operator= ( const Parameter& rhs )
  {
    this->clear();

    _name = rhs._name;
    _description  = rhs._description;
    _et = rhs._et;
    _string_pattern = rhs._string_pattern;
    _do_not_check_me = rhs._do_not_check_me;

    for(unsigned int n=0; n<rhs._values.size(); ++n)
      _values.push_back( rhs._values[n]->clone() );

    return *this;
  }

#if 0


    /**
     * set by template type
     */
    template <typename T>
    T& set (unsigned int idx=0)
    {
      Value<T>* value = 0;
      if(_values.size() < idx+1)
      {
        _values.resize(idx+1, 0);
        _values[idx] = new Value<T>;
      }
      else
      {
        if(!_values[idx])
          _values[idx] = new Value<T>;
      }
      value = dynamic_cast<Value<T>*>(_values[idx]);
      return  value->set();
    }

    /**
     * set bool value
     */
    void    set_bool(const bool v, unsigned int idx=0)
    { this->set<bool>(idx)=v; _et = BOOL; }

    /**
     * set integer value
     */
    void    set_int(const int v, unsigned int idx=0)
    { this->set<int>(idx)=v; _et = INTEGER;}

    /**
     * set real value
     */
    void   set_real(const double v, unsigned int idx=0)
    {this->set<double>(idx)=v; _et = REAL;}

    /**
     * set string value
     */
    void   set_string(const std::string & s, unsigned int idx=0)
    { this->set<std::string>(idx)=s; _et = STRING; }

    /**
     * set ENUM value, we should check if s fit the string pattern
     */
    int    set_enum(const std::string & s, unsigned int idx=0)
    {
      // no pattern exist
      if ( _string_pattern.size() == 0 )
      { this->set<std::string>(idx)=s; _et = ENUM; return 0; }

      // s is in the pattern
      if ( _string_pattern.count(s) )
      { this->set<std::string>(idx)=s; _et = ENUM; return 0; }

      // return error
      return 1;
    }

    /**
     * set bool array
     */
    void    set_bool_array(const std::vector<bool> &v)
    {
      for(unsigned int n=0; n<v.size(); ++n)
        this->set<bool>(n)=v[n];
      _et = BOOL;
    }

    /**
     * set integer array
     */
    void    set_int_array(const std::vector<int> & v)
    {
      for(unsigned int n=0; n<v.size(); ++n)
        this->set<int>(n)=v[n];
      _et = INTEGER;
    }

    /**
     * set real array
     */
    void   set_real_array(const std::vector<double> & v)
    {
      for(unsigned int n=0; n<v.size(); ++n)
        this->set<double>(n)=v[n];
      _et = REAL;
    }

    /**
     * set string array
     */
    void   set_string_array(const std::vector<std::string> & s)
    {
      for(unsigned int n=0; n<s.size(); ++n)
        this->set<std::string>(n)=s[n];
      _et = STRING;
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
     * @return the number of values
     */
    unsigned int array_size() const
      { return _values.size(); }

    /**
     * get by template type
     */
    template <typename T>
    const T & get (unsigned int idx=0) const
      {
        Value<T>* value = dynamic_cast<Value<T>*>(_values[idx]);
        assert(value);
        return  value->get();
      }

    /**
     * get bool value
     */
    bool    get_bool(unsigned int idx=0)   const
      {return this->get<bool>(idx);}

    /**
     * get integer value
     */
    int     get_int(unsigned int idx=0)    const
      {return this->get<int>(idx);}

    /**
     * get real value
     */
    double  get_real(unsigned int idx=0)   const
      {return this->get<double>(idx);}

    /**
     * get string value
     */
    const std::string & get_string(unsigned int idx=0) const
      {return this->get<std::string>(idx);}

    /**
     * get array
     */
    template <typename T>
    std::vector<T> get_array() const
    {
      std::vector<T>  array;
      for(unsigned int n=0; n<_values.size(); ++n)
        array.push_back(this->get<T>(n));
      return array;
    }

    /**
     * clear the Parameter
     */
    void clear ()
    {
      _name.clear();
      _description.clear();
      _string_pattern.clear();
      for(unsigned int n=0; n<_values.size(); ++n)
        delete _values[n];
      _values.clear();
      do_not_check_me = false;
    }

    Parameter & operator= ( const Parameter& rhs )
    {
      this->clear();

      _name = rhs._name;
      _description  = rhs._description;
      _et = rhs._et;
      _string_pattern = rhs._string_pattern;
      do_not_check_me = rhs.do_not_check_me;

      for(unsigned int n=0; n<rhs._values.size(); ++n)
        _values.push_back( rhs._values[n]->clone() );

      return *this;
    }
#endif

  //explicit instantiation

  template bool &         Parameter::set(unsigned int);
  template int &          Parameter::set(unsigned int);
  template double &       Parameter::set(unsigned int);
  template std::string &  Parameter::set(unsigned int);
  template const bool &         Parameter::get(unsigned int)   const;
  template const int &          Parameter::get(unsigned int)   const;
  template const double &       Parameter::get(unsigned int)   const;
  template const std::string &  Parameter::get(unsigned int)   const;
  template std::vector<bool>        Parameter::get_array()   const;
  template std::vector<int>         Parameter::get_array()   const;
  template std::vector<double>      Parameter::get_array()   const;
  template std::vector<std::string> Parameter::get_array()   const;

}



