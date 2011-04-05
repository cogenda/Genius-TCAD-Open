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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "genius_common.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "field_source.h"
#include "solver_specify.h"
#include "log.h"
#include "parallel.h"


using  PhysicalUnit::s;


FieldSource::FieldSource(SimulationSystem & system, Parser::InputParser & decks)
    :_system(system), _decks(decks), current_waveform(0)
{

  // check if any waveform defined
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    if( c.key() == "ENVELOP" )   // It's a ENVELOP (Waveform) card
    {
      if(c.is_parameter_exist("type"))
      {
        if(c.is_enum_value("type","uniform"))            SetWaveformUniform(c);
        else if(c.is_enum_value("type","gaussian"))      SetWaveformGauss(c);
        else if(c.is_enum_value("type","pulse"))         SetWaveformPulse(c);
        else if(c.is_enum_value("type","expression"))    SetWaveformExpr(c);
        else if(c.is_enum_value("type","shell"))         SetWaveformShell(c);
      }
    }
  }

  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    const Parser::Card c = _decks.get_current_card();

    // parse input card to find if PARTICLE source exist
    if(c.key() == "PARTICLE" )
    {
      // PROFILE card should have "type" parameter
      if( !c.is_parameter_exist("profile") )
      {
        MESSAGE<<"ERROR at " << c.get_fileline() <<" PARTICLE: Should have a 'profile' parameter."<<std::endl; RECORD();
        genius_error();
      }

      if( c.is_enum_value("profile", "fromfile2d") || c.is_enum_value("profile", "fromfile3d") )
      {
        Particle_Source * particle_source = new Particle_Source_DataFile(system, c);
        add_particle_source(particle_source);
      }

      if( c.is_enum_value("profile", "analytic") )
      {
        Particle_Source * particle_source = new Particle_Source_Analytic(system, c);
        add_particle_source(particle_source);
      }

      if( c.is_enum_value("profile", "track") )
      {
        Particle_Source * particle_source = new Particle_Source_Track(system, c);
        add_particle_source(particle_source);
      }

    }

    // parse input card to find if light source exist
    if(c.key() == "LIGHT" )
    {
      // read optical wave from spectrum file
      if(c.is_parameter_exist("spectrumfile"))
      {
        parse_spectrum_file(c.get_string("spectrumfile", ""), system, c);
      }
      else
      {
        double wave_length   = c.get_real("wavelength", 0.532, "lambda") * PhysicalUnit::um;   // wave length
        double power         = c.get_real("intensity",  0.0) * PhysicalUnit::J/PhysicalUnit::s/PhysicalUnit::cm/PhysicalUnit::cm;      // incident power;
        double eta           = c.get_real("quan.eff", 1.0);
        bool eta_auto        = !c.is_parameter_exist("quan.eff");

        Light_Source * light_source = new Light_Source_From_File(system, c, "", wave_length, power, eta, eta_auto);
        add_light_source(light_source);
      }
    }

    if(c.key() == "RAYTRACE" )
    {
      Light_Source * light_source = new Light_Source(system);
      add_light_source(light_source);
    }

    if(c.key() == "EMFEM2D")
    {
      Light_Source * light_source = new Light_Source(system);
      add_light_source(light_source);
    }
  }
}


FieldSource::~FieldSource()
{
  // free all the field sources

  std::vector<Particle_Source *>::iterator pit = _particle_sources.begin();
  for(; pit!=_particle_sources.end(); ++pit)
    delete (*pit);

  std::vector<Light_Source *>::iterator lit = _light_sources.begin();
  for(; lit!=_light_sources.end(); ++lit)
    delete (*lit);

  // free waveforms
  std::map<std::string, Waveform *>::iterator it= _waveforms.begin();
  for(; it!=_waveforms.end(); ++it )
    delete it->second;
}


void FieldSource::update(double time)
{
  // fast detect if we need to do something
  if(!SolverSpecify::PatG && !SolverSpecify::OptG) return;

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    if( region->type()== SemiconductorRegion)
    {
      SimulationRegion::processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        FVM_Node * fvm_node = (*it);
        FVM_NodeData * fvm_node_data = fvm_node->node_data();
        genius_assert(fvm_node_data);

        double G=0;

        if(SolverSpecify::PatG)
        {
          std::vector<Particle_Source *>::iterator pit = _particle_sources.begin();
          for(; pit!=_particle_sources.end(); ++pit)
          {
            G += fvm_node_data->PatG()*(*pit)->carrier_generation(time);
          }
        }

        if(SolverSpecify::OptG)
        {
          double A = 1.0;
          if(current_waveform)
            A = current_waveform->waveform(time);
          G += fvm_node_data->OptG()*A;
        }

        fvm_node_data->Field_G() = G;

      }
    }
  }
#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}


void FieldSource::update_system()
{
  std::vector<Particle_Source *>::iterator pit = _particle_sources.begin();
  for(; pit!=_particle_sources.end(); ++pit)
    (*pit)->update_system();

  std::vector<Light_Source *>::iterator lit = _light_sources.begin();
  for(; lit!=_light_sources.end(); ++lit)
    (*lit)->update_system();
}


void FieldSource::parse_spectrum_file(const std::string & filename, SimulationSystem & system, const Parser::Card &c)
{
  // only processor 0 read the spectrum file

  std::vector<Light_Source_From_File * > sources;

  std::ifstream in;
  if(Genius::processor_id()==0)
  {
    in.open(filename.c_str());
    genius_assert(in.good());
  }

  int finished = 0;
  while (! finished )
  {
    std::string fname_ext;
    double wave_length, power, eta;
    int eta_auto;

    if(Genius::processor_id()==0)
    {
      std::string line;

      while (1)
      {
        std::getline(in, line);

        if(in.eof())
        {
          //  reached eof
          in.close();
          finished=1;
          break;
        }

        // skip empty lines
        if(line.size()==0) continue;
        // skip the line begin with '#'
        if(line.find('#')==0) continue;

        // we found a line
        break;
      }

      if (!finished)
      {
        std::stringstream ss;
        ss << line;

        // read wave length and power information from line
        ss >> fname_ext;
        ss >> wave_length;
        ss >> power;

        if(ss.fail())
        {
          MESSAGE<<"ERROR at " << c.get_fileline() <<": error reading spectrum file " << filename << std::endl; RECORD();
          genius_error();
        }

        // test if it contains extra parameter for quantum efficiency
        ss >> eta;
        // if no quantum efficiency parameter find, set eta_auto flag to true
        eta_auto = ss.fail()?1:0;
      }
    }

    Parallel::broadcast(finished);
    if(finished)
      break;

    Parallel::broadcast(fname_ext);
    Parallel::broadcast(wave_length);
    Parallel::broadcast(power);
    Parallel::broadcast(eta);
    Parallel::broadcast(eta_auto);

    wave_length *= PhysicalUnit::um;

    // this is power density (per um wave length)
    power       *= PhysicalUnit::J/PhysicalUnit::s/(PhysicalUnit::cm*PhysicalUnit::cm) / PhysicalUnit::um;

    Light_Source *light_source;
    if (eta_auto)
    {
      light_source = new Light_Source_From_File(system, c, fname_ext, wave_length, power, 1.0, true);
    }
    else
    {
      light_source = new Light_Source_From_File(system, c, fname_ext, wave_length, power, eta, false);
    }

    sources.push_back(dynamic_cast<Light_Source_From_File *>(light_source));

    add_light_source(light_source);
  }

  unsigned int n_sources = sources.size();

  //post process
  genius_assert(n_sources>=2);

  double total_power = 0.0;
  std::vector<double> lambda_distance;
  for(unsigned int n=0; n<n_sources; ++n)
  {
    double ld;
    if(n==0)
      ld = 0.5*(sources[1]->waveLength() -sources[0]->waveLength());
    else if(n==n_sources-1)
      ld = 0.5*(sources[n_sources-1]->waveLength() - sources[n_sources-2]->waveLength());
    else
      ld = 0.5*(sources[n+1]->waveLength() - sources[n-1]->waveLength());

    if(ld<=0)
    {
      MESSAGE << "Wavelength must be strictly increasing." << std::endl; RECORD();
      genius_error();
    }

    lambda_distance.push_back(ld);
    total_power += sources[n]->power()*lambda_distance[n];
  }

  for(unsigned int n=0; n<n_sources; ++n)
  {
    sources[n]->setPower((sources[n]->power()*lambda_distance[n]));
  }

  MESSAGE<<"Read "<<n_sources<<" sources from the spectrum file "<<filename<<".\n"
  <<"Total power intensity is approximate "<<std::setprecision(3)<<std::setw(10)
  <<total_power /( PhysicalUnit::J/PhysicalUnit::s/(PhysicalUnit::cm*PhysicalUnit::cm) )
  <<" W/(cm^2)"
  <<std::endl;
  RECORD();

}



void  FieldSource::SetWaveformUniform(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  assert( label!="" );

  double amplitude = c.get_real("amplitude", 1.0);
  double Tdelay = c.get_real("tdelay",0.0)*s;

  _waveforms[label] = new WaveformUniform(label, Tdelay, amplitude);
}



void FieldSource::SetWaveformGauss(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  assert( label!="" );

  double amplitude = c.get_real("amplitude", 1.0);
  double t0 = c.get_real("t0",0.0)*s;
  double tao = c.get_real("tao",1e-12)*s;

  _waveforms[label] = new WaveformGauss(label, t0, tao, amplitude);
}



void  FieldSource::SetWaveformPulse(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  assert( label!="" );

  double Tdelay = c.get_real("tdelay",0.0)*s;
  double tr = c.get_real("tr",1e-12)*s;
  double tf = c.get_real("tf",1e-12)*s;
  double pw = c.get_real("pw",9e-12)*s;
  double pr = c.get_real("pr",20e-12)*s;
  double Alo = c.get_real("amplitude.low",0.0);
  double Ahi = c.get_real("amplitude.high",1.0);

  _waveforms[label] = new WaveformPulse(label, Tdelay, Alo, Ahi, tr, tf, pw, pr);

}



void  FieldSource::SetWaveformExpr(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  assert( label!="" );

  std::string expr = c.get_string("expression","1.0");

  _waveforms[label] = new WaveformExpression(label, expr);
}




void  FieldSource::SetWaveformShell(const Parser::Card &c)
{

  std::string label = c.get_string("id","");
  assert( label!="" );

  std::string filename = c.get_string("dll","");
  std::string funcname = c.get_string("function","");

  _waveforms[label] = new WaveformShell(label, filename, funcname, s);

}



