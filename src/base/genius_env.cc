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


#include "genius_common.h"
#include "genius_env.h"

#ifdef HAVE_SLEPC
  #include "slepcsys.h"
#endif

// ------------------------------------------------------------
// Genius::GeniusPrivateData data initialization
int  Genius::GeniusPrivateData::_n_processors = 1;
int  Genius::GeniusPrivateData::_processor_id = 0;

#ifdef HAVE_MPI
MPI_Comm Genius::GeniusPrivateData::_comm_world;
MPI_Comm Genius::GeniusPrivateData::_comm_self;
#endif

std::string Genius::GeniusPrivateData::_input_file;
std::string Genius::GeniusPrivateData::_genius_dir;

bool Genius::GeniusPrivateData::_experiment_code=true;

bool Genius::init_processors(int *argc, char *** args)
{
  // GENIUS is built on top of PETSC, we should init PETSC first
#ifdef HAVE_SLEPC
  // if we have slepc, call  SlepcInitialize instead of PetscInitialize
  SlepcInitialize(argc,args,PETSC_NULL,PETSC_NULL);
#else
  PetscInitialize(argc, args, PETSC_NULL, PETSC_NULL);
#endif

#ifdef HAVE_MPI
  // the actual process number
  MPI_Comm_size (PETSC_COMM_WORLD, &Genius::GeniusPrivateData::_n_processors);
  MPI_Comm_rank (PETSC_COMM_WORLD, &Genius::GeniusPrivateData::_processor_id);

  // duplicate an other MPI_Comm for Genius parallel communication
  MPI_Comm_dup( PETSC_COMM_WORLD, &Genius::GeniusPrivateData::_comm_world );
  MPI_Comm_dup( PETSC_COMM_SELF, &Genius::GeniusPrivateData::_comm_self );
#endif

  return true;
}

bool Genius::clean_processors()
{

#ifdef HAVE_MPI
  MPI_Comm_free(&Genius::GeniusPrivateData::_comm_world);
  MPI_Comm_free(&Genius::GeniusPrivateData::_comm_self);
#endif

  // end PETSC
#ifdef  HAVE_SLEPC
  SlepcFinalize();
#else
  PetscFinalize();
#endif


  return true;
}
