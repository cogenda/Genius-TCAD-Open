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

#ifndef _genius_env_h_
#define _genius_env_h_

#include "config.h"
#include "genius_petsc.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif


#include <cstdlib>
#include <cstring>
#include <string>


namespace Genius {

  /**
   *  Initialize Petsc
   *  @returns true on success.
   */
  bool init_processors(int *argc, char *** args);

  /**
   * Clean up Petsc
   * @returns true on success.
   */
  bool clean_processors();

  /**
   * @returns the number of processors used in the current simulation.
   */
  unsigned int n_processors();

  /**
   * @return the first processor id, should always be 0
   */
  unsigned int first_processor_id();

  /**
   * @return the last processor id
   */
  unsigned int last_processor_id();

  /**
   * @returns the index of the local processor.
   */
  unsigned int processor_id();

  /**
   * @return true if we are in first processor
   */
  bool is_first_processor();

  /**
   * @return true if we are in last processor
   */
  bool is_last_processor();

#ifdef HAVE_MPI
  /**
   * @return MPI_Comm global communicator
   */
  const MPI_Comm & comm_world();

  /**
   * @return MPI_Comm self communicator
   */
  const MPI_Comm & comm_self();
#endif


  /**
   * @returns the input filename;
   */
  const char * input_file();

  /**
   * Set the input filename
   */
  void set_input_file(const char* fname);

  /**
   * @returns the genius base directory
   */
  std::string genius_dir();

  /**
   * Set the genius base directory
   */
  void set_genius_dir(const std::string &genius_dir);

  /**
   * set flag to enable experiment code
   */
  void set_experiment_code(bool f);

  /**
   * flag to enable experiment code
   */
  bool experiment_code();

  /**
   * current memory usage of this processor in \<virtual memory size, resident Set Size \>
   */
  std::pair<size_t, size_t> memory_size();

  /**
   * Namespaces don't provide private data,
   * so let's take the data we would like
   * private and put it in an obnoxious
   * namespace.  At least that way it is a
   * pain to use, thus discouraging errors.
   */
  class GeniusPrivateData {
  public:
    /**
     * Total number of processors used.
     */
    static int  _n_processors;

    /**
     * The local processor id.
     */
    static int  _processor_id;

#ifdef HAVE_MPI
    /**
     * MPI_Comm global communicator
     */
    static MPI_Comm _comm_world;

    /**
     * MPI_Comm local communicator
     */
    static MPI_Comm _comm_self;

#endif

    /**
     * the user input file.
     */
    static std::string _input_file;

    /**
     * Genius base dir
     */
    static std::string _genius_dir;

    /**
     * flag to enable experiment code
     */
    static bool _experiment_code;

  };
}



// ------------------------------------------------------------
// Genius inline member functions
inline unsigned int Genius::n_processors()
{
  return static_cast<unsigned int>(GeniusPrivateData::_n_processors);
}


inline unsigned int Genius::processor_id()
{
  return static_cast<unsigned int>(GeniusPrivateData::_processor_id);
}


inline unsigned int Genius::first_processor_id()
{
  return static_cast<unsigned int>(0);
}


inline unsigned int Genius::last_processor_id()
{
  return static_cast<unsigned int>(GeniusPrivateData::_n_processors-1);
}


inline bool Genius::is_first_processor()
{
  return GeniusPrivateData::_processor_id == 0;
}


inline bool Genius::is_last_processor()
{
  return GeniusPrivateData::_processor_id == GeniusPrivateData::_n_processors-1;
}


#ifdef HAVE_MPI
inline  const MPI_Comm & Genius::comm_world()
{
  return (GeniusPrivateData::_comm_world);
}

inline  const MPI_Comm & Genius::comm_self()
{
  return (GeniusPrivateData::_comm_self);
}
#endif


inline const char * Genius::input_file()
{
  return GeniusPrivateData::_input_file.c_str();
}

inline void Genius::set_input_file(const char* fname)
{
  GeniusPrivateData::_input_file = std::string(fname);
}

inline std::string Genius::genius_dir()
{
  return GeniusPrivateData::_genius_dir;
}

inline void Genius::set_genius_dir(const std::string &genius_dir)
{
  GeniusPrivateData::_genius_dir = genius_dir;
}

inline void Genius::set_experiment_code(bool f)
{
  GeniusPrivateData::_experiment_code = f;
}

inline bool Genius::experiment_code()
{
  return GeniusPrivateData::_experiment_code;
}


// These are useful macros that behave like functions in the code.
// If you want to make sure you are accessing a section of code just
// stick a here(); in it, for example
#undef genius_here
#define genius_here()     { std::cout << "[" << Genius::processor_id() << "] " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; }

// the stop() macro will stop the code until a SIGCONT signal is received.  This is useful, for example, when
// determining the memory used by a given operation.  A stop() could be inserted before and after a questionable
// operation and the delta memory can be obtained from a ps or top.  This macro only works for serial cases.
#undef genius_stop
#ifdef HAVE_CSIGNAL
#  include <csignal>
#  define genius_stop()     { if (Genius::n_processors() == 0) { genius_here(); std::cout << "Stopping process " << getpid() << "..." << std::endl; std::raise(SIGSTOP); std::cout << "Continuing process " << getpid() << "..." << std::endl; } }
#else
#  define genius_stop()     { if (Genius::n_processors() == 0) { genius_here(); std::cerr << "WARNING:  stop() does not work without the <csignal> header file!" << std::endl; } }
#endif


// Use Petsc error mechanism for assertion tests
#undef genius_assert
#ifdef WINDOWS
#  include <csignal>
#  define genius_assert(a)  {  if (! (a) ) { std::cerr << "[" << Genius::processor_id() << "] " << "Assertion failure at " << __FILE__ << ":" << __LINE__ << std::endl; std::raise(SIGTERM); } }
#else
#  define genius_assert(a)  {  if (! (a) ) { std::cerr << "[" << Genius::processor_id() << "] " << "Assertion failure at " << __FILE__ << ":" << __LINE__ << std::endl; std::abort(); } }
#endif



#endif // #define _genius_env_h_
