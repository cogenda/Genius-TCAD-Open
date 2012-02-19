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


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "genius_common.h"
#include "genius_env.h"
#include "material_define.h"
#include "file_include.h"
#include "sync_file.h"
#include "perf_log.h"
#include "parser.h"
#include "control.h"
#include "parallel.h"



#ifdef WINDOWS
  #include <io.h>      // for windows _access function
#else
  #include <unistd.h>  // for POSIX access function
#endif

static void show_logo();

#if PETSC_VERSION_GE(3,2,0)
static PetscErrorCode genius_error_handler(MPI_Comm comm, int line, const char *func, const char *file, const char *dir,PetscErrorCode n, PetscErrorType p, const char *mess,void *ctx);
#else
static PetscErrorCode genius_error_handler(int line, const char *func, const char *file, const char *dir,PetscErrorCode n, int p, const char *mess,void *ctx);
#endif


// --------------------------------------------------------
// The entrance of GENIUS
int main(int argc, char ** args)
{
  Genius::init_processors(&argc, &args);

  //  PetscExceptionPush(-1);
  PetscPushErrorHandler(genius_error_handler, NULL);

  //show the GENIUS Log
  show_logo();

#ifndef COGENDA_COMMERCIAL_PRODUCT
  if( Genius::n_processors() > 1 )
  {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: Open Source Version does not support multi-processor.\n");
    PetscFinalize();
    exit(0);
  }
#endif

  // record the start time
  PetscLogDouble t_start;
  PetscGetTime(&t_start);

  // count the number of user's input argument
  if(argc<2)
  {
    PetscPrintf(PETSC_COMM_WORLD,"usage: mpirun -n [1-9]+ genius -i card_file [petsc_option]\n");
    PetscFinalize();
    exit(0);
  }


  // test if GENIUS_DIR has been set correctly
  if( getenv("GENIUS_DIR") == NULL )
  {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: User should set entironment variable GENIUS_DIR.\n");
    PetscFinalize();
    exit(0);
  }
  Genius::set_genius_dir(getenv("GENIUS_DIR"));

  // performace log flag
  PetscBool     log_flg;
  PetscOptionsHasName(PETSC_NULL,"-p", &log_flg);
  if(!log_flg)
    perflog.disable_logging();

  // get the name of user input file by PETSC routine
  {
    PetscBool     file_flg;
    char *petsc_arg_buffer = new char[1024];
    PetscOptionsGetString(PETSC_NULL, "-i", petsc_arg_buffer, 1023, &file_flg);
    // no input file? exit...
    if( !file_flg )
    {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR: I want an input file to tell me what to do.\n");
      PetscFinalize();
      exit(0);
    }
    Genius::set_input_file(petsc_arg_buffer);
    delete [] petsc_arg_buffer;
  }

  {
    PetscBool     experiment_code_flg;
    PetscOptionsHasName(PETSC_NULL,"-e", &experiment_code_flg);
    if(experiment_code_flg) Genius::set_experiment_code(false);
  }

  // prepare log system
  std::ofstream logfs;
  if (Genius::processor_id() == 0)
  {
    genius_log.addStream("console", std::cerr.rdbuf());
    std::stringstream log_file;
    log_file << Genius::input_file() << ".log";
    logfs.open(log_file.str().c_str());
    genius_log.addStream("file", logfs.rdbuf());
  }

  MESSAGE<<"Genius boot with " << Genius::n_processors() << " MPI thread.\n\n";  RECORD();

  // test if input file can be opened on processor 0 for read
  if ( Genius::processor_id() == 0 )
  {
#ifdef WINDOWS
    if ( _access( Genius::input_file(),  04 ) == -1 )
#else
    if ( access( Genius::input_file(),  R_OK ) == -1 )
#endif
    {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR: I can't read input file '%s', access failed.\n", Genius::input_file() );
      PetscFinalize();
      exit(0);
    }
  }

  // preprocess include statement of input file
  std::string input_file_pp;
  if (Genius::processor_id() == 0)
  {
    Parser::FilePreProcess * file_preprocess = new Parser::FilePreProcess(Genius::input_file());
    input_file_pp = file_preprocess->output();
    delete file_preprocess;
  }
  Parallel::broadcast(input_file_pp);

  // sync input file to other processor
  const std::string localfile = sync_file(input_file_pp.c_str());

  // read card specification file
  Parser::Pattern pt;
  std::string pattern_file = Genius::genius_dir() +  "/lib/GeniusSyntax.xml";

  // test if pattern file can be open for read
#ifdef WINDOWS
  if ( _access( pattern_file.c_str(),  04 ) == -1 )
#else
  if (  access( pattern_file.c_str(),  R_OK ) == -1 )
#endif
  {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: I can't read pattern file at %s, access failed.\n", pattern_file.c_str() );
    PetscFinalize();
    exit(0);
  }


  if (pt.get_from_XML(pattern_file) )
  {
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: I can't parse pattern file 'GeniusSyntax.xml'.\n" );
    genius_error();
  }

  // parse the input file
  AutoPtr<Parser::InputParser> input = AutoPtr<Parser::InputParser>(new Parser::InputParser(pt));
  if (input->read_card_file(localfile.c_str()) )
  {
    // remove preprocessed file
    if (Genius::processor_id() == 0)
      remove(input_file_pp.c_str());

    remove(localfile.c_str());

    PetscPrintf(PETSC_COMM_WORLD,"ERROR: I can't parse input file.\n");
    PetscFinalize();
    exit(0);
  }

  // remove preprocessed file
  if (Genius::processor_id() == 0)
    remove(input_file_pp.c_str());

  // after that, remove local copy of input file
  remove(localfile.c_str());

  // set material define
  std::string material_file = Genius::genius_dir() +  "/lib/material.def";
  Material::init_material_define(material_file);

  // do solve process here
  AutoPtr<SolverControl>  solve_ctrl = AutoPtr<SolverControl>(new SolverControl());
  solve_ctrl->setDecks(input.get());
  {
    std::stringstream fsol;
    fsol << Genius::input_file() << ".sol";
    solve_ctrl->setSolutionFile(fsol.str().c_str());
  }
  solve_ctrl->mainloop();

  // record the end time
  PetscLogDouble t_end;
  PetscGetTime(&t_end);
  PetscLogDouble elapsed_time = t_end - t_start;

  // change to mm:ss format
  int    min = static_cast<int>(elapsed_time/60);
  double sec = elapsed_time - min*60;
  std::stringstream time_ss;
  time_ss<<"Genius finished. Simulation time is "
  << min <<" min "
  << std::setiosflags(std::ios::fixed) <<  std::setprecision(3)
  << sec<<" second."
  <<" Good bye." << std::endl;

  MESSAGE<<time_ss.str(); RECORD();

  // performace log
  if(log_flg)
  {
    std::string perf_info = perflog.get_log();
    std::vector<char> perf_info_char(perf_info.size());
    std::copy(perf_info.begin(), perf_info.end(), perf_info_char.begin());
    Parallel::allgather(perf_info_char);
    perf_info_char.push_back(0);
    perf_info = std::string(&perf_info_char[0]);
    MESSAGE<<perf_info; RECORD();
  }

  //finish log system
  if (Genius::processor_id() == 0)
  {
    genius_log.removeStream("console");
    genius_log.removeStream("file");
    logfs.close();
  }

  Genius::clean_processors();
  return 0;
}



void show_logo()
{
#ifdef COGENDA_COMMERCIAL_PRODUCT
  // show the GENIUS Log
  // FIXME this logo seems not so pretty...
  PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
  PetscPrintf(PETSC_COMM_WORLD,"*     888888    88888888   88     888  88888   888     888    8888888   *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*   8       8   8          8 8     8     8      8       8    8          *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*  8            8          8  8    8     8      8       8    8          *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*  8            88888888   8   8   8     8      8       8     888888    *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*  8      8888  8          8    8  8     8      8       8           8   *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*   8       8   8          8     8 8     8      8       8           8   *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*     888888    88888888  888     88   88888     8888888     8888888    *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*                                                                       *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*  Parallel Three-Dimensional General Purpose Semiconductor Simulator   *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*                                                                       *\n");
  if(sizeof(PetscScalar)==sizeof(double))
    PetscPrintf(PETSC_COMM_WORLD,"*        Commercial Version %-8s with double precision.          *\n", PACKAGE_VERSION);
  else if(sizeof(PetscScalar)==sizeof(long double))
    PetscPrintf(PETSC_COMM_WORLD,"*        Commercial Version %-8s with long double precision.     *\n", PACKAGE_VERSION);
  PetscPrintf(PETSC_COMM_WORLD,"*                                                                       *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*      Copyright (C) 2007-2010 by Cogenda Pte Ltd.                      *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*                http://www.cogenda.com                                 *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
#else
  PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
  PetscPrintf(PETSC_COMM_WORLD,"*     888888    88888888   88     888  88888   888     888    8888888   *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*   8       8   8          8 8     8     8      8       8    8          *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*  8            8          8  8    8     8      8       8    8          *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*  8            88888888   8   8   8     8      8       8     888888    *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*  8      8888  8          8    8  8     8      8       8           8   *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*   8       8   8          8     8 8     8      8       8           8   *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*     888888    88888888  888     88   88888     8888888     8888888    *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*                                                                       *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*  Parallel Three-Dimensional General Purpose Semiconductor Simulator   *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*                                                                       *\n");
  if(sizeof(PetscScalar)==sizeof(double))
    PetscPrintf(PETSC_COMM_WORLD,"*        Open Source Version %-8s with double precision.         *\n", PACKAGE_VERSION);
  else if(sizeof(PetscScalar)==sizeof(long double))
    PetscPrintf(PETSC_COMM_WORLD,"*        Open Source Version %-8s with long double precision.    *\n", PACKAGE_VERSION);
  PetscPrintf(PETSC_COMM_WORLD,"*                                                                       *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*      Copyright (C) 2007-2010 by Cogenda Pte Ltd.                      *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*                http://www.cogenda.com                                 *\n");
  PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
#endif

  PetscSynchronizedFlush(PETSC_COMM_WORLD);

}


#if PETSC_VERSION_GE(3,2,0)
PetscErrorCode genius_error_handler(MPI_Comm comm, int line, const char *func, const char *file, const char *dir,PetscErrorCode n, PetscErrorType p, const char *mess,void *ctx)
#else
PetscErrorCode genius_error_handler(int line, const char *func, const char *file, const char *dir,PetscErrorCode n, int p, const char *mess,void *ctx)
#endif
{

  MESSAGE << "--------------------- Error Message ------------------------------------" << std::endl;
  MESSAGE << "Fatal Error:";
  if (mess)
    MESSAGE << mess;
  MESSAGE << " at line " << line << " in " << file << std::endl;
  MESSAGE << "------------------------------------------------------------------------" << std::endl << std::endl << std::endl;
  RECORD();

  MPI_Abort(PETSC_COMM_WORLD, -1);
  exit(0);
}

