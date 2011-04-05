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


#ifndef __tif3d_io_h__
#define __tif3d_io_h__


// Local includes
#include "genius_common.h"
#include "field_input.h"
#include "field_output.h"
#include "simulation_system.h"


// Forward declarations
class MeshBase;



/**
 * This class implements reading and writing solution data in the
 * TIF3D (Technology Interchange Format 3D) file.
 * This file format is used by our GDSII to 3D mesh tool.
 */
class TIF3DIO : public FieldInput<SimulationSystem>,
      public FieldOutput<SimulationSystem>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  TIF3DIO (SimulationSystem& system);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  TIF3DIO (const SimulationSystem& system);


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
  enum TIF3D_Record_Type
  {
    H, // TIF file header
    C, // Nodal coordinate
    F, // Edge record
    T, // tet record
    R, // Region record
    S, // Solution record
    N  // Nodal solution record
  };
*/

};



// ------------------------------------------------------------
// TIF3DIO inline members
inline
TIF3DIO::TIF3DIO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system)
{
}



inline
TIF3DIO::TIF3DIO (const SimulationSystem& system) :
    FieldOutput<SimulationSystem>(system)
{
}



#endif // #define __tif3d_io_h__

