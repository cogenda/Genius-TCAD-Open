// $Id: dof_object.h,v 1.2 2008/05/22 04:53:44 gdiso Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef __dof_object_h__
#define __dof_object_h__

// C++ includes


// Local includes
#include "genius_env.h"
#include "genius_common.h"


// Forward declarations
class DofObject;



/**
 * The \p DofObject defines an abstract base class for objects that
 * have degrees of freedom associated with them.  Examples of such
 * objects are the \p Node and \p Elem classes.  This class can
 * not be instantiated, only derived from.
 */
class DofObject
{

protected:

  /**
   * Constructor. Protected so that you can't instantiate one of these.
   */
  DofObject ();

public:

  /**
   * Copy-constructor.
   */
  DofObject (const DofObject&);

  /**
   * Destructor.
   */
  virtual ~DofObject ();

  /**
   * Sets the id to \p invalid_id
   */
  void invalidate_id ();

  /**
   * Sets the processor id to \p invalid_processor_id
   */
  void invalidate_processor_id ();

  /**
   * Invalidates all the indices for this \p DofObject
   */
  void invalidate ();

  /**
   * \returns the \p id for this \p DofObject
   */
  unsigned int id () const;

  /**
   * \returns the \p id for this \p DofObject as a writeable reference.
   */
  unsigned int & set_id ();

  /**
   * Sets the \p id for this \p DofObject
   */
  void set_id (const unsigned int id)
  { this->set_id() = id; }

  /**
   * @returns \p true if this \p DofObject has a valid \p id set,
   * \p false otherwise.
   */
  bool valid_id () const;

  /**
   * @return the first global degree of freedom number of this dof object
   */
  unsigned int global_dof_id() const
  { return _global_dof_id; }

  /**
   * @return the reference of first global degree of freedom number of this dof object
   */
  unsigned int & global_dof_id()
  { return _global_dof_id; }

  /**
   * @return the first local degree of freedom number of this dof object
   */
  unsigned int local_dof_id() const
  { return _local_dof_id; }

  /**
   * @return the reference of first local degree of freedom number of this dof object
   */
  unsigned int & local_dof_id()
  { return _local_dof_id; }

  /**
   * @returns the processor that this element belongs to.
   * To conserve space this is stored as a short integer.
   */
  unsigned short int processor_id () const;

  /**
   * @returns the processor that this element belongs to as a
   * writeable reference.
   */
  unsigned short int & processor_id ();

  /**
   * Sets the \p processor_id for this \p DofObject.
   */
  void processor_id (const unsigned int id);

  /**
   * @returns \p true if this \p DofObject has a valid \p id set,
   * \p false otherwise.
   */
  bool valid_processor_id () const;

  /**
   * @returns true if the DofObject on processor
   */
  bool on_processor () const
  { return _processor_id == Genius::processor_id(); }

  /**
   * @returns true if the DofObject on local
   */
  bool on_local () const
  { return _on_local; }

  /**
   * @returns writable reference to on_local
   */
  bool & on_local ()
  { return _on_local; }

  /**
   * @return true if this object is ghost object
   * a ghost object is on local but its processor_id != Genius::processor_id()
   */
  bool is_ghost () const
  { return _on_local && _processor_id != Genius::processor_id(); }


  /**
   * Implemented in Elem and Node.
   */
  virtual bool operator==(const DofObject& ) const
  { genius_error(); return false; }


  /**
   * An invaild \p id to distinguish an uninitialized \p DofObject
   */
  static const unsigned int invalid_id;

  /**
   * An invalid \p processor_id to distinguish DOFs that have
   * not been assigned to a processor.
   */
  static const unsigned short int invalid_processor_id;


private:


  /**
   * The \p id of the \p DofObject
   */
  unsigned int _id;

  /**
   * The first global degree of freedom number
   */
  unsigned int _global_dof_id;

  /**
   * The first local degree of freedom number
   */
  unsigned int _local_dof_id;

  /**
   * The \p processor_id of the \p DofObject.
   * Degrees of freedom are wholly owned by processors,
   * however they may be duplicated on other processors.
   *
   * This is stored as an unsigned short int since we cannot
   * expect to be solving on 65000+ processors any time soon,
   * can we??
   */
  unsigned short int _processor_id;

  /**
   * this bool value indicate of the DofObject should be hold on local processor
   * ture for both on processor and ghost objects
   */
  bool     _on_local;

};



//------------------------------------------------------
// Inline functions
inline
DofObject::DofObject () :
  _id (invalid_id),
  _global_dof_id(invalid_id),
  _local_dof_id(invalid_id),
  _processor_id (invalid_processor_id),
  _on_local(false)
{
  this->invalidate();
}





inline
DofObject::~DofObject ()
{
}



inline
void DofObject::invalidate_id ()
{
  this->set_id (invalid_id);
}



inline
void DofObject::invalidate_processor_id ()
{
  this->processor_id (invalid_processor_id);
}



inline
void DofObject::invalidate ()
{
  this->invalidate_id ();
  this->invalidate_processor_id ();
}



inline
unsigned int DofObject::id () const
{
  return _id;
}



inline
unsigned int & DofObject::set_id ()
{
  return _id;
}



inline
bool DofObject::valid_id () const
{
  return (DofObject::invalid_id != _id);
}


inline
unsigned short int DofObject::processor_id () const
{
  return _processor_id;
}



inline
unsigned short int & DofObject::processor_id ()
{
  return _processor_id;
}



inline
void DofObject::processor_id (const unsigned int id)
{
  this->processor_id() = id;
}



inline
bool DofObject::valid_processor_id () const
{
  return (DofObject::invalid_processor_id != _processor_id);
}



#endif // #ifndef __dof_object_h__
