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


#ifndef __tif_io_h__
#define __tif_io_h__


// Local includes
#include "genius_common.h"
#include "field_input.h"
#include "field_output.h"
#include "simulation_system.h"


// Forward declarations
class MeshBase;



/**
 * This class implements reading and writing solution data in the
 * TIF (Technology Interchange Format) file.
 * This file format is used by Synopsys Medici(R) and Tsuprem(R).
 */
class TIFIO : public FieldInput<SimulationSystem>,
      public FieldOutput<SimulationSystem>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  TIFIO (SimulationSystem& system);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  TIFIO (const SimulationSystem& system);


  /**
   * This method implements reading a mesh from a specified file
   * in TIF format.
   */
  virtual void read (const std::string& filename);

  /**
   * This method implements writing a mesh to a specified TIF file.
   * NOT implemented yet!
   */
  virtual void write (const std::string& )
  { genius_error(); }


private:
/*
  enum TIF_Record_Type
  {
    h, // TIF file header
    c, // Nodal coordinate
    e, // Edge record
    r, //Region record
    b, //Boundary record
    i, //Interface record
    j, //Interface boundary record
    t, //Triangle record
    s, //Solution record
    n  //Nodal solution record
  };
*/

};



// ------------------------------------------------------------
// TIFIO inline members
inline
TIFIO::TIFIO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system)
{
}



inline
TIFIO::TIFIO (const SimulationSystem& system) :
    FieldOutput<SimulationSystem>(system)
{
}



#endif // #define __tif_io_h__

