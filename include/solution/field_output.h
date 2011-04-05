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

//  $Id: field_output.h,v 1.2 2008/07/09 05:58:16 gdiso Exp $


#ifndef __field_output_h__
#define __field_output_h__


// C++ inludes
#include <string>



/**
 * This class defines an abstract interface for \p solution data output.
 * Specific classes derived from this class actually implement
 * writing various storage formats.
 */

// ------------------------------------------------------------
// FieldOutput class definition
template <class MT>
class FieldOutput
{
 protected:

  /**
   * Default constructor. Will set the _obj to NULL, effectively
   * rendering this object useless.
   */
  FieldOutput ();

  /**
   * Constructor.  Takes a reference to a constant object.
   * This constructor will only allow us to write the object.
   */
  FieldOutput (const MT&);


 public:

  /**
   * Destructor.
   */
  virtual ~FieldOutput ();

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string&) = 0;


 protected:


  /**
   * Returns the object as a read-only reference.
   */
  const MT& system() const;


 private:


  /**
   *  A pointer to a constant object.
   * This allows us to write the object to file.
   */
  const MT* const _obj;

};



// ------------------------------------------------------------
// FieldOutput inline members
template <class MT>
inline
FieldOutput<MT>::FieldOutput (const MT& obj) :
  _obj (&obj)
{
}



template <class MT>
inline
FieldOutput<MT>::~FieldOutput ()
{
}




template <class MT>
inline
const MT& FieldOutput<MT>::system () const
{
  return *_obj;
}


#endif // #define __field_output_h__
