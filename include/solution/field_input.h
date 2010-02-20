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

//  $Id: field_input.h,v 1.2 2008/07/09 05:58:16 gdiso Exp $



#ifndef __field_input_h__
#define __field_input_h__


// C++ inludes
#include <istream>
#include <string>

// Local includes
#include "genius_env.h"
#include "genius_common.h"




/**
 * This class defines an abstract interface for \p solution data input.
 * Specific classes derived from this class actually implement
 * reading various data storage formats.
 */

// ------------------------------------------------------------
// FieldInput class definition
template <class MT>
class FieldInput
{
 protected:

  /**
   * Default constructor. Will set the _obj to NULL, effectively
   * rendering this object useless.
   */
  FieldInput ();
  
  /**
   * Constructor.  Takes a writeable reference to an object.
   * This is the constructor required to read an object.
   */
  FieldInput (MT&);
  
 public:

  /**
   * Destructor.
   */
  virtual ~FieldInput ();
  
  /**
   * This method implements reading data from a specified file.
   */
  virtual void read (const std::string&) = 0;

  
 protected:
  
  /**
   * Returns the object as a writeable reference.
   */
  MT& system ();
  
 
  
 private:
  

  /**
   * A pointer to a non-const object object.
   * This allows us to read the object from file.
   */ 
  MT* _obj;
};



// ------------------------------------------------------------
// FieldInput inline members
template <class MT>
inline
FieldInput<MT>::FieldInput () :
  _obj (NULL)
{
}



template <class MT>
inline
FieldInput<MT>::FieldInput (MT& obj) :
  _obj (&obj)
{
}



template <class MT>
inline
FieldInput<MT>::~FieldInput ()
{
}



template <class MT>
inline
MT& FieldInput<MT>::system ()
{
  if (_obj == NULL) genius_error();
  return *_obj;
}







#endif // #define __field_input_h__
