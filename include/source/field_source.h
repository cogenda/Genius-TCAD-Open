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


#ifndef __field_source_h__
#define __field_source_h__

#include <vector>
#include <map>
#include "waveform.h"
#include "particle_source.h"
#include "light_source.h"

namespace Parser{
  class InputParser;
  class Card;
}
class LightLenses;
class SimulationSystem;

/**
 * manage all the field sources
 */
class FieldSource
{
public:

  /**
   * constructor
   */
  FieldSource(SimulationSystem & system, Parser::InputParser & decks);

  ~FieldSource();

  /**
   * update the source stimulate to each mesh node
   * when force is true, will update source
   */
  void update(double time, bool force_update_source=false);

  /**
   * update source
   */
  void update_source();

  /**
   * @return applied_to_system flag
   */
  bool applied_to_system() const { return _applied_to_system; }

  /**
   * @return the limited time step
   */
  double limit_dt(double time, double dt, double dt_min) const;

  /**
   * @return true when we have particle incident
   */
  bool is_particle_source_exist() const
    { return _particle_sources.size()>0; }

  /**
   * @return lens system
   */
  LightLenses  * light_lenses() const
  { return _light_lenses; }

  /**
   * @return true when we have light source
   */
  bool is_light_source_exist() const
    { return _light_sources.size()>0; }

  /**
   * if we need serial mesh on each processor
   * @return true when ray tracing solver exist
   */
  bool request_serial_mesh() const;

  /**
   * set the effect waveform
   */
  bool set_effect_waveform(const std::string & str)
  {
    if(_waveforms.find(str)!=_waveforms.end())
    {
      current_waveform = _waveforms[str];
      for(unsigned int n=0; n<_light_sources.size(); ++n)
        _light_sources[n]->set_global_waveform(current_waveform);
      return true;
    }
    else
    {
      current_waveform = 0;
      for(unsigned int n=0; n<_light_sources.size(); ++n)
        _light_sources[n]->set_global_waveform(current_waveform);
      return false;
    }
    return false;
  }

  /**
   * clear effect waveform
   */
  void clear_effect_waveform()  { current_waveform = 0; }

private:

  /**
   * field source is performed on the system
   */
  SimulationSystem & _system;

  /**
  * since reading deck involves stack operate, we can not use const here
  */
  Parser::InputParser & _decks;

  /**
   * add a particle source to this source manager
   */
  void add_particle_source(Particle_Source * ps)
  { _particle_sources.push_back(ps); }

  /**
   * add a light source to this source manager
   */
  void add_light_source(Light_Source * ls)
  { _light_sources.push_back(ls); }


  //void parse_transform(const Parser::Card &c);

  std::vector<Light_Source_From_File *> parse_spectrum_file(const std::string & filename, SimulationSystem & system, const Parser::Card &c);

  /**
   * flag to indicate that the field source has been applied
   */
  bool _applied_to_system;

  /**
   * all the particle sources
   */
  std::vector<Particle_Source *> _particle_sources;

  /**
   * all the light sources
   */
  std::vector<Light_Source *> _light_sources;

  /**
   * lens system
   */
  LightLenses  *_light_lenses;

  /**
   * all the waveforms
   */
  std::map<std::string, Waveform *>  _waveforms;

  /**
   * the effect waveform
   */
  Waveform * current_waveform;

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformUniform(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformSin(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformGauss(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformPulse(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformDoubleExp(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformExpr(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformShell(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformWaveFile(const Parser::Card &c);
};

#endif
