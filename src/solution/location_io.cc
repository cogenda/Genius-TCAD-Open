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

#include <fstream>
#include <iomanip>

#include "mesh_base.h"
#include "location_io.h"

using PhysicalUnit::m;

void LocationIO::write (const std::string& filename)
{
  const MeshBase& mesh = FieldOutput<SimulationSystem>::system().mesh();

  unsigned int n_nodes = mesh.n_nodes();

  // collect mesh nodes, should be executed in parallel
  std::vector<Real> pts;
  mesh.pack_nodes(pts);

  if(Genius::processor_id() == 0)
  {
    std::ofstream out(filename.c_str());
    genius_assert(out.good());

    // set the float number precision
    out.precision(6);

    // set output width and format
    out << std::scientific << std::right;

    // write out nodal data
    for ( unsigned int n=0; n<n_nodes; ++n )
    {
      if (_numbering)
        out << n <<'\t';

      out << std::setw(15) << pts[3*n+0]/m/_length_unit;
      out << std::setw(15) << pts[3*n+1]/m/_length_unit;
      out << std::setw(15) << pts[3*n+2]/m/_length_unit;
      out << std::endl;
    }

    out.close();
  }
}

