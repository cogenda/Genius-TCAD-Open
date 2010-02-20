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



#ifndef __dfise_block_h__
#define __dfise_block_h__

#include <cassert>

#include <map>
#include <vector>
#include <string>


namespace DFISE
{

  /**
   * use this stupid struct to contain int/double/std::string date
   */
  struct TOKEN
  {
    enum TOKEN_TYPE {int_token, float_token, string_token};
    TOKEN_TYPE    token_type;
    void * value;

    void free_value()
    {
      switch(token_type)
      {
      case int_token     : delete (int *)value; break;
      case float_token   : delete (double *)value; break;
      case string_token  : delete (std::string *)value; break;
      }
    }

  };


  /**
   * data structure for DF-ISE "block"
   * (A DFISE file comprises a sequence of blocks)
   */
  class BLOCK
  {
  public:

    /// empty constructor
    BLOCK() {}

    /// constructor with keyword
    BLOCK(const std::string & k):_keyword(k) {}

    /**
     * don't free any TOKEN except call clear()!
     */
    ~BLOCK() {}

    void set_keyword(const std::string & k)
    { _keyword = k; }

    void set_label(const std::string & l)
    { _label = l; }

    void set_index(int i)
    { _index = i; }

    void add_parameter(const std::string & p, const std::vector<TOKEN *> & v)
    { _parameters[p] = v; }

    void add_values(const std::vector<TOKEN *> & new_value)
    {
      for(unsigned int n=0; n<new_value.size(); ++n)
        _values.push_back(new_value[n]);
    }

    void add_sub_block(BLOCK * sub_block)
    { _sub_blocks.push_back(sub_block); }

    void append(const BLOCK & block)
    {
      std::map<std::string, std::vector<TOKEN *> >::const_iterator it = block._parameters.begin();
      for(; it!=block._parameters.end(); ++it)
        _parameters[it->first] = it->second;

      for(unsigned int n=0; n<block._values.size(); ++n)
        _values.push_back(block._values[n]);

      for(unsigned int n=0; n<block._sub_blocks.size(); ++n)
        _sub_blocks.push_back(block._sub_blocks[n]);
    }

    /**
     * free everything
     */
    void clear()
    {
      std::map<std::string, std::vector<TOKEN *> >::iterator it = _parameters.begin();
      for(; it!=_parameters.end(); ++it)
        clear(it->second);
      _parameters.clear();

      clear(_values);

      for(unsigned int n=0; n<_sub_blocks.size(); ++n)
      {
        _sub_blocks[n]->clear();
        delete _sub_blocks[n];
      }
      _sub_blocks.clear();
    }

    void print() const
    {
      std::cout<<_keyword<<std::endl;
      std::cout<<" parameters:"<< _parameters.size() <<std::endl;
      std::cout<<" values:"<< _values.size() <<std::endl;
      std::cout<<" sub blocks:"<< _sub_blocks.size() <<std::endl;
      for(unsigned int n=0; n<_sub_blocks.size(); ++n)
        _sub_blocks[n]->print();
    }


    unsigned int n_sub_blocks() const
    { return _sub_blocks.size(); }

    /**
     * @return sub block by its index
     */
    BLOCK * get_sub_block(unsigned int n)
    { return _sub_blocks[n]; }

    /**
     * @return the first sub block which march the keyword
     */
    BLOCK * get_sub_block(const std::string &k)
    {
      for(unsigned int n=0; n<_sub_blocks.size(); ++n)
        if(_sub_blocks[n]->keyword()==k)
          return _sub_blocks[n];
      return 0;
    }

    const std::string & keyword() const
    { return _keyword; }

    int index() const
      { return _index; }

    const std::string & label() const
      { return _label; }

    /**
     * @return the size of value associated with this parameter
     */
    unsigned int n_values_in_parameter(const std::string & name)
    {
      assert(_parameters.find(name)!=_parameters.end());
      std::vector<TOKEN *> & token = _parameters[name];
      return token.size();
    }

    /**
     * @return the size of values
     */
    unsigned int n_values() const
    {
      return _values.size();
    }

    /**
     * @return the ith pointer to std::string value of parameter name
     */
    std::string get_string_parameter(const std::string & name, unsigned int i)
    {
      assert(_parameters.find(name)!=_parameters.end());
      std::vector<TOKEN *> & token = _parameters[name];
      assert(i<token.size());
      assert(token[i]->token_type == TOKEN::string_token);
      return *(std::string*)token[i]->value;
    }

    /**
     * @return the  ith pointer to int value of parameter name
     */
    int get_int_parameter(const std::string & name, unsigned int i)
    {
      assert(_parameters.find(name)!=_parameters.end());
      std::vector<TOKEN *> & token = _parameters[name];
      assert(i<token.size());
      assert(token[i]->token_type == TOKEN::int_token);
      return *(int*)token[i]->value;
    }

    /**
     * @return the ith pointer to float value of parameter name
     */
    double get_float_parameter(const std::string & name, unsigned int i)
    {
      assert(_parameters.find(name)!=_parameters.end());
      std::vector<TOKEN *> & token = _parameters[name];
      assert(i<token.size());
      assert(token[i]->token_type == TOKEN::int_token || token[i]->token_type == TOKEN::float_token);
      if(token[i]->token_type == TOKEN::int_token)
        return *(int*)token[i]->value;
      return *(double*)token[i]->value;
    }

    std::string get_string_value(unsigned int i)
    {
      assert(i<_values.size());
      assert(_values[i]->token_type == TOKEN::string_token);
      return *(std::string*)_values[i]->value;
    }

    int get_int_value(unsigned int i)
    {
      assert(i<_values.size());
      assert(_values[i]->token_type == TOKEN::int_token);
      return *(int*)_values[i]->value;
    }

    double get_float_value(unsigned int i)
    {
      assert(i<_values.size());
      assert(_values[i]->token_type == TOKEN::int_token || _values[i]->token_type == TOKEN::float_token);
      if(_values[i]->token_type == TOKEN::int_token)
        return *(int*)_values[i]->value;
      return *(double*)_values[i]->value;
    }

  public:

    /**
     * the keyword of the block
     */
    std::string _keyword;

    /**
     * addtional label of the block
     */
    std::string _label;

    /**
     * addtional integer value of the block
     */
    int  _index;

    /**
     * all the parameters in this block
     */
    std::map<std::string, std::vector<TOKEN *> >  _parameters;

    /**
     * all the individual values in this block
     */
    std::vector<TOKEN *>  _values;

    /**
     * the sub blocks
     */
    std::vector<BLOCK *> _sub_blocks;

    /**
     * free tokens and the value in tokens!
     */
    void clear(std::vector<TOKEN *> & tokens)
    {
      for(unsigned int n=0; n<tokens.size(); ++n)
      {
        tokens[n]->free_value();
        delete tokens[n];
      }
      tokens.clear();
    }
  };


}


#endif
