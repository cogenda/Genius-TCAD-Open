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

//  $Id: cgns_io.h,v 1.7 2008/07/09 09:51:29 gdiso Exp $


#ifndef __cgns_io_h__
#define __cgns_io_h__

// C++ includes
#include <map>

// Local includes
#include "genius_common.h"
#include "field_input.h"
#include "field_output.h"
#include "simulation_system.h"


// Forward declarations
class MeshBase;



/**
 * This class implements reading and writing solution data in the
 * CGNS (CFD General Notation System) format.
 * Format description:
 * cf. <a href="http://www.cgns.org/">CGNS home page</a>.
 */

// ------------------------------------------------------------
// CGNSIO class definition
class CGNSIO : public FieldInput<SimulationSystem>,
	       public FieldOutput<SimulationSystem>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  CGNSIO (SimulationSystem& system);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  CGNSIO (const SimulationSystem& system);


  /**
   * This method implements reading a mesh from a specified file
   * in CGNS format.
   */
  virtual void read (const std::string& );

  /**
   * This method implements writing a mesh to a specified cgns file.
   */
  virtual void write (const std::string& );


private:
// some id for cgns file read/write

  /**
   *  cgns file id
   */
  int fn;

  /**
   *  cgns base id
   */
  int B;

  /**
   *  cgns zone id
   */
  int Z;

  /**
   *  cgns coordinate id
   */
  int C;

  /**
   *  cgns section id
   */
  int S;

  /**
   *  cgns boundary condition id
   */
  int BC;

  /**
   *  cgns 1-to-1 connect id
   */
  int I;

  /**
   *  cgns solution id
   */
  int SOL;

  /**
   *  cgns field id
   */
  int F;

  /**
   * store region node's global id
   */
  std::vector< std::vector<int> > region_global_id;

  /**
   * map global id to node's id.
   */
  std::vector< std::map<int, unsigned int > > region_global_id_to_node_id;

  /**
   * region solution
   */
  std::vector< std::map<std::pair<std::string,std::string>, std::vector<double> > > region_solutions;

  /**
   * unit of solution
   */
  std::map< std::string, std::string > region_solution_units;

  /**
   * electrode potential
   */
  std::map< short int, double > electrode_potential;

  /**
   * aux function to sort x by increase order of id
   */
  std::vector<double> & _sort_it (std::vector<double> & x, const std::vector<unsigned int > &id);

};



// ------------------------------------------------------------
// CGNSIO inline members
inline
CGNSIO::CGNSIO (SimulationSystem& system) :
	FieldInput<SimulationSystem> (system),
	FieldOutput<SimulationSystem> (system)
{

}



inline
CGNSIO::CGNSIO (const SimulationSystem& system) :
	FieldOutput<SimulationSystem>(system)
{

}



#endif // #define __cgns_io_h__

