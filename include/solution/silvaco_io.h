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


#ifndef __silvaco_io_h__
#define __silvaco_io_h__


// Local includes
#include "genius_common.h"
#include "field_input.h"
#include "field_output.h"
#include "simulation_system.h"
#include "simulation_region.h"

// Forward declarations
class MeshBase;



/**
 * This class implements reading and writing solution data in the
 * Silvaco TIF (very similar to synopsys TIF) file.
 * This file format is used by Silvaco ATLAS(R) and other tools.
 */
class STIFIO : public FieldInput<SimulationSystem>,
               public FieldOutput<SimulationSystem>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  STIFIO (SimulationSystem& system);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  STIFIO (const SimulationSystem& system);


  /**
   * This method implements reading a mesh from a specified file
   * in Silvaco TIF format.
   */
  virtual void read (const std::string& filename);

  /**
   * This method implements writing a mesh to a specified Silvaco TIF file.
   * NOT implemented yet!
   */
  virtual void write (const std::string& )
  { genius_error(); }

};



// ------------------------------------------------------------
// STIFIO inline members
inline
STIFIO::STIFIO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system)
{
}



inline
STIFIO::STIFIO (const SimulationSystem& system) :
    FieldOutput<SimulationSystem>(system)
{
}



#endif // #define __silvaco_io_h__

