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
#include "parser.h"
#include "waveform.h"
#include "particle_source.h"
#include "light_source.h"


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
   */
  void update(double time);

  /**
   * update system after mesh refine
   */
  void update_system();

  /**
   * @return true when we have particle incident
   */
  bool is_particle_source_exist() const
    { return _particle_sources.size()>0; }

  /**
   * @return true when we have light source
   */
  bool is_light_source_exist() const
    { return _light_sources.size()>0; }

  /**
   * set the effect waveform
   */
  void set_effect_waveform(const std::string & str)
  {
    if(_waveforms.find(str)!=_waveforms.end())
      current_waveform = _waveforms[str];
    else
      current_waveform = 0;
  }

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


  void parse_transform(const Parser::Card &c);

  void parse_spectrum_file(const std::string & filename, SimulationSystem & system, const Parser::Card &c);

  /**
   * all the particle sources
   */
  std::vector<Particle_Source *> _particle_sources;

  /**
   * all the light sources
   */
  std::vector<Light_Source *> _light_sources;

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
  void  SetWaveformGauss(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformPulse(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformExpr(const Parser::Card &c);

  /**
   * private functions for setting each waveform
   */
  void  SetWaveformShell(const Parser::Card &c);
};

#endif
