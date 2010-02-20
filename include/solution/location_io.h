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


#ifndef __location_io_h__
#define __location_io_h__


// Local includes
#include "genius_common.h"
#include "field_input.h"
#include "field_output.h"
#include "simulation_system.h"

// Forward declarations
class MeshBase;



/**
 * This class write out all the mesh nodes
 */
class LocationIO : public FieldInput<SimulationSystem>,
      public FieldOutput<SimulationSystem>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  LocationIO (SimulationSystem& system);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  LocationIO (const SimulationSystem& system);


  /**
   * read node location
   */
  virtual void read (const std::string&)
  { genius_error(); }
 
  /**
   *write node location
   */
  virtual void write (const std::string&);

  /**
   * set length unit, with respect to SI unit (m)
   */
  void setLengthUnit (const PetscScalar unit)
  {
    genius_assert(unit>0);
    _length_unit = unit;
  }

  /**
   * set if node number should be printed
   */
  void setNumbering (const bool numbering)
  {
    _numbering = numbering;
  }

private:
  PetscScalar _length_unit;
  bool _numbering;
};


// ------------------------------------------------------------
// LocationIO inline members
inline
LocationIO::LocationIO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system),
    _length_unit (1.0),
    _numbering(true)
{
}



inline
LocationIO::LocationIO (const SimulationSystem& system) :
    FieldOutput<SimulationSystem>(system),
    _length_unit (1.0),
    _numbering(true)
{
}

#endif // #define __location_io_h__

 
