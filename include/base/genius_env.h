#ifndef _genius_env_h_
#define _genius_env_h_

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

    /**
     * the user input file.
     */
    static std::string _input_file;

    /**
     * Genius base dir
     */
    static std::string _genius_dir;

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

#endif // #define _genius_env_h_
