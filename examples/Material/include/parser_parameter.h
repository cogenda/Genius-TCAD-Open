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

# ifndef  __parser_parameter_h__
# define  __parser_parameter_h__

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <typeinfo>

#include <cassert>
#include <cstdarg>


#include "config.h"




namespace Parser
{

  enum  ElemType {BOOL, INTEGER, REAL, STRING, ENUM, INVALID};


  /**
   * The parameter class, contains a special parameter (name, type, value)
   */
  class Parameter
  {

  public:
    // constructors
    Parameter():_do_not_check_me(false) {}

    Parameter(const Parameter &p)
        :_name(p._name), _description(p._description), _et(p._et),
         _string_pattern(p._string_pattern),_do_not_check_me(p._do_not_check_me)
    {
      for(unsigned int n=0; n<p._values.size(); ++n)
        _values.push_back( p._values[n]->clone() );
    }

    // constructors
    //Parameter(const std::string &name):_name(name), do_not_check_me(false) {}
    Parameter(const std::string &name, bool v):_name(name), _do_not_check_me(false) { this->set_bool(v); }
    Parameter(const std::string &name, int v):_name(name), _do_not_check_me(false) { this->set_int(v); }
    Parameter(const std::string &name, double v):_name(name), _do_not_check_me(false) { this->set_real(v); }
    Parameter(const std::string &name, const std::string s):_name(name), _do_not_check_me(false) { this->set_string(s); }
    Parameter(const std::string &name, const char * s): _name(name),_do_not_check_me(false) { this->set_string(std::string(s)); }

    ~Parameter() { this->clear(); }

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
     * set by template type
     */
    template <typename T>
    T& set (unsigned int idx=0);

    /**
     * set bool value
     */
    void    set_bool(const bool v, unsigned int idx=0);

    /**
     * set integer value
     */
    void    set_int(const int v, unsigned int idx=0);

    /**
     * set real value
     */
    void   set_real(const double v, unsigned int idx=0);

    /**
     * set string value
     */
    void   set_string(const std::string & s, unsigned int idx=0);

    /**
     * set ENUM value, we should check if s fit the string pattern
     */
    int    set_enum(const std::string & s, unsigned int idx=0);

    /**
     * set bool array
     */
    void    set_bool_array(const std::vector<bool> &v);

    /**
     * set integer array
     */
    void    set_int_array(const std::vector<int> & v);

    /**
     * set real array
     */
    void   set_real_array(const std::vector<double> & v);

    /**
     * set string array
     */
    void   set_string_array(const std::vector<std::string> & s);

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
    int string_pattern_match ( std::string &s );


    /**
     * @return the number of values
     */
    unsigned int array_size() const;

    /**
     * get by template type
     */
    template <typename T>
        const T & get (unsigned int idx=0) const;

    /**
     * get bool value
     */
    bool    get_bool(unsigned int idx=0)   const;

    /**
     * get integer value
     */
    int     get_int(unsigned int idx=0)    const;

    /**
     * get real value
     */
    double  get_real(unsigned int idx=0)   const;

    /**
     * get string value
     */
    const std::string & get_string(unsigned int idx=0) const;

    /**
     * get array
     */
    template <typename T>
    std::vector<T> get_array() const;

    /**
     * clear the Parameter
     */
    void clear ();

    Parameter & operator= ( const Parameter& rhs );


    /**
     * set this parameter is user defined
     */
    void  set_user_defined()
    { _do_not_check_me = true; }


    /**
     * @return ture if this parameter is user defined
     */
    bool   is_user_defined() const
    { return _do_not_check_me; }


    typedef std::set<std::string>::const_iterator StringEnumIterator;
    StringEnumIterator stringPatternBegin() const { return _string_pattern.begin(); }
    StringEnumIterator stringPatternEnd() const { return _string_pattern.end(); }

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
     * Abstract definition of a attribute value.
     */
    class ValueBase
    {
    public:

      /**
       * Destructor.
       */
      virtual ~ValueBase() {}

      /**
       * String identifying the type of attribute stored.
       * Must be reimplemented in derived classes.
       */
      virtual std::string type () const = 0;

      /**
       * Clone this value.  Useful in copy-construction.
       * Must be reimplemented in derived classes.
       */
      virtual ValueBase* clone () const = 0;
    };

    /**
     * Concrete definition of a attribute value
     * for a specified type.
     */
    template <typename T>
    class Value : public ValueBase
    {
    public:

      ~Value() {}

      /**
       * @returns a read-only reference to the attribute value.
       */
      const T& get () const { return _value; }

      /**
       * @returns a writeable reference to the attribute value.
       */
      T& set () { return _value; }

      /**
       * String identifying the type of attribute stored.
       */
      std::string type () const { return typeid ( T ).name();  }

      /**
       * Clone this value.  Useful in copy-construction.
       */
      ValueBase* clone () const
      {
        Value<T> *copy = new Value<T>;
        copy->_value = _value;
        return copy;
      }

    private:

      /**
       * Stored attribute value.
       */
      T _value;
    };

    std::vector<ValueBase *> _values;

    /**
     * pattern for string value
     */
    std::set<std::string> _string_pattern;


    /**
     * for user defined parameter, we should not check it with pattern
     * since it will not exist in the pattern
     */
    bool _do_not_check_me;

  };


}
# endif
