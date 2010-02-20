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

#include "mesh.h"
#include "location_io.h"

using PhysicalUnit::m;

void LocationIO::write (const std::string& filename)
{
  if(Genius::processor_id() != 0) return;

  std::ofstream out(filename.c_str());
  genius_assert(out.good());

  // set the float number precision
  out.precision(6);

  // set output width and format
  out << std::scientific << std::right;

  const MeshBase& mesh = FieldOutput<SimulationSystem>::system().mesh();

  MeshBase::const_node_iterator nd = mesh.active_nodes_begin();
  MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();
  // write out nodal data
  for ( ; nd!=nd_end; ++nd )
  {
    const Node * node = (*nd);
    {
      if (_numbering)
        out << node->id()<<'\t';

      out << std::setw(15) << (*node)(0)/m/_length_unit;
      out << std::setw(15) << (*node)(1)/m/_length_unit;
      out << std::setw(15) << (*node)(2)/m/_length_unit;
      out << std::endl;
    }
  }

  out.close();
}

