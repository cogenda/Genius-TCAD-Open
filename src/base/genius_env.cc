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

#include "genius_env.h"
#include "genius_common.h"

// ------------------------------------------------------------
// Genius::GeniusPrivateData data initialization
int  Genius::GeniusPrivateData::_n_processors = 1;
int  Genius::GeniusPrivateData::_processor_id = 0;
std::string Genius::GeniusPrivateData::_input_file;
std::string Genius::GeniusPrivateData::_genius_dir;

bool Genius::init_processors(int *argc, char *** args)
{
  // GENIUS is built on top of PETSC, we should init PETSC first
  PetscInitialize(argc, args, PETSC_NULL, PETSC_NULL);

  // the actual process number
  MPI_Comm_rank (PETSC_COMM_WORLD, &Genius::GeniusPrivateData::_processor_id);
  MPI_Comm_size (PETSC_COMM_WORLD, &Genius::GeniusPrivateData::_n_processors);

  return true;
}

bool Genius::clean_processors()
{
  // end PETSC
  PetscFinalize();
  return true;
}
