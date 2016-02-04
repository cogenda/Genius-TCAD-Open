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
#include "parser.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "field_source.h"
#include "light_lenses.h"
#include "solver_specify.h"
#include "log.h"
#include "parallel.h"



using  PhysicalUnit::s;


FieldSource::FieldSource(SimulationSystem & system, Parser::InputParser & decks)
  :_system(system), _decks(decks), _applied_to_system(false), current_waveform(0)
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
        else if(c.is_enum_value("type","sine"))          SetWaveformSin(c);
        else if(c.is_enum_value("type","gaussian"))      SetWaveformGauss(c);
        else if(c.is_enum_value("type","pulse"))         SetWaveformPulse(c);
        else if(c.is_enum_value("type","dexp"))          SetWaveformDoubleExp(c);
        else if(c.is_enum_value("type","expression"))    SetWaveformExpr(c);
        else if(c.is_enum_value("type","shell"))         SetWaveformShell(c);
        else if(c.is_enum_value("type","wavefile"))      SetWaveformWaveFile(c);
      }
    }
  }

  // if any lens defined?
  _light_lenses = new LightLenses(_decks);

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
      // uniform generation
      if(c.is_parameter_exist("optical.gen") || c.is_parameter_exist("optical.power"))
      {
        Light_Source * light_source = new Light_Source_Uniform(system, c);
        if(c.is_parameter_exist("envelop"))
        {
          std::string waveform = c.get_string("envelop", "");
          if(_waveforms.find(waveform)!=_waveforms.end())
          {
            light_source->set_waveform(_waveforms[waveform]);
          }
          else
          {
            MESSAGE<<"ERROR at " << c.get_fileline() <<" LIGHT: envelop " << waveform << " not be defined."<<std::endl; RECORD();
            genius_error();
          }
        }
        add_light_source(light_source);
        continue;
      }

      // read optical wave from spectrum file
      if(c.is_parameter_exist("spectrumfile"))
      {
        std::vector<Light_Source_From_File *> sources = parse_spectrum_file(c.get_string("spectrumfile", ""), system, c);

        Waveform * w = 0;
        if(c.is_parameter_exist("envelop"))
        {
          std::string waveform = c.get_string("envelop", "");
          if(_waveforms.find(waveform)!=_waveforms.end())
          {
            w = _waveforms[waveform];
          }
          else
          {
            MESSAGE<<"ERROR at " << c.get_fileline() <<" LIGHT: envelop " << waveform << " not be defined."<<std::endl; RECORD();
            genius_error();
          }
        }

        for(unsigned int n=0; n<sources.size(); ++n)
        {
           sources[n]->set_waveform(w);
           add_light_source(sources[n]);
        }
        continue;
      }

      // single spectrum
      {
        double wave_length   = c.get_real("wavelength", 0.532, "lambda") * PhysicalUnit::um;   // wave length
        double power         = c.get_real("intensity",  0.0) * PhysicalUnit::J/PhysicalUnit::s/PhysicalUnit::cm/PhysicalUnit::cm;      // incident power;
        double eta           = c.get_real("quan.eff", 1.0);
        bool eta_auto        = !c.is_parameter_exist("quan.eff");

        Light_Source * light_source = new Light_Source_From_File(system, c, "", wave_length, power, eta, eta_auto);
        if(c.is_parameter_exist("envelop"))
        {
          std::string waveform = c.get_string("envelop", "");
          if(_waveforms.find(waveform)!=_waveforms.end())
          {
            light_source->set_waveform(_waveforms[waveform]);
          }
          else
          {
            MESSAGE<<"ERROR at " << c.get_fileline() <<" LIGHT: envelop " << waveform << " not be defined."<<std::endl; RECORD();
            genius_error();
          }
        }
        add_light_source(light_source);
      }
    }

    if(c.key() == "RAYTRACE" )
    {
      Light_Source * light_source = new Light_Source_RayTracing(system, c);
      if(c.is_parameter_exist("envelop"))
      {
        std::string waveform = c.get_string("envelop", "");
        if(_waveforms.find(waveform)!=_waveforms.end())
        {
          light_source->set_waveform(_waveforms[waveform]);
        }
        else
        {
          MESSAGE<<"ERROR at " << c.get_fileline() <<" LIGHT: envelop " << waveform << " not be defined."<<std::endl; RECORD();
          genius_error();
        }
      }
      add_light_source(light_source);
    }

    if(c.key() == "EMFEM2D")
    {
      Light_Source * light_source = new Light_Source_EMFEM2D(system, c);
      add_light_source(light_source);
    }


    if(c.key() == "XRAYPULSE")
    {
       double dose_rate   = c.get_real("doserate", 0.0) * PhysicalUnit::rad/PhysicalUnit::s;
       Light_Source * light_source = new Light_Source_Xray(system, dose_rate);

       if( c.is_enum_value("waveform", "pulse") )
       {
         double t0          = c.get_real("t0", 0.0) * PhysicalUnit::s;
         double tr          = c.get_real("tr", 10e-9) * PhysicalUnit::s; //10ns
         double tf          = c.get_real("tf", 30e-9) * PhysicalUnit::s; //30ns
         double pw          = c.get_real("pw", 10e-9) * PhysicalUnit::s; //10ns

         Waveform * w = new WaveformPulse("__x.ray__", t0, 0.0, 1.0, tr, tf, pw, 1e30* PhysicalUnit::s);
         _waveforms["__x.ray__"]  = w;
         light_source->set_waveform(w);
       }

       if( c.is_enum_value("waveform", "gaussian") )
       {
         double t0          = c.get_real("t0", 0.0) * PhysicalUnit::s;
         double fwhm        = c.get_real("fwhm", 50e-9) * PhysicalUnit::s; //50ns
         double tao         = fwhm/1.665;
         Waveform * w = new WaveformGauss("__x.ray__", t0, tao, 1.0);
         _waveforms["__x.ray__"] = w;
         light_source->set_waveform(w);
       }

       add_light_source(light_source);
    }
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
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

  //free lens
  delete _light_lenses;
}


void FieldSource::update(double time, bool force_update_source)
{
  // fast detect if we need to do something
  if(!SolverSpecify::PatG && !SolverSpecify::OptG) return;

  if( force_update_source )    this->update_source();

  if( _applied_to_system == false ) this->update_source();


  // clear old particle and optical generation
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    {
      SimulationRegion::processor_node_iterator it = region->on_local_nodes_begin();
      SimulationRegion::processor_node_iterator it_end = region->on_local_nodes_end();
      for(; it!=it_end; ++it)
      {
        FVM_Node * fvm_node = (*it);
        FVM_NodeData * fvm_node_data = fvm_node->node_data();
        fvm_node_data->PatG() = 0.0;
        fvm_node_data->OptG() = 0.0;
      }
    }
  }

  // let particle source update the PatG
  std::vector<Particle_Source *>::iterator pit = _particle_sources.begin();
  for(; pit!=_particle_sources.end(); ++pit)
  {
    (*pit)->carrier_generation(time);
  }

  // let light source update the OptG
  std::vector<Light_Source *>::iterator lit = _light_sources.begin();
  for(; lit!=_light_sources.end(); ++lit)
    (*lit)->carrier_generation(time);


  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    {
      SimulationRegion::processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        FVM_Node * fvm_node = (*it);
        FVM_NodeData * fvm_node_data = fvm_node->node_data();

        double G=0;

        if(SolverSpecify::PatG)
        {
          G += fvm_node_data->PatG();
        }

        if(SolverSpecify::OptG)
        {
          G += fvm_node_data->OptG();
        }

        fvm_node_data->Field_G() = G;
      }
    }
  }
#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}


void FieldSource::update_source()
{
  // clear old particle and optical generation
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    {
      SimulationRegion::processor_node_iterator it = region->on_local_nodes_begin();
      SimulationRegion::processor_node_iterator it_end = region->on_local_nodes_end();
      for(; it!=it_end; ++it)
      {
        FVM_Node * fvm_node = (*it);
        FVM_NodeData * fvm_node_data = fvm_node->node_data();
        fvm_node_data->PatG() = 0.0;
        fvm_node_data->OptG() = 0.0;
      }
    }
  }

  // calculate particle generation
  std::vector<Particle_Source *>::iterator pit = _particle_sources.begin();
  for(; pit!=_particle_sources.end(); ++pit)
    (*pit)->update_source();

  // calculate optical generation
  std::vector<Light_Source *>::iterator lit = _light_sources.begin();
  for(; lit!=_light_sources.end(); ++lit)
    (*lit)->update_source();

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  _applied_to_system = true;
}



double FieldSource::limit_dt(double time, double dt, double dt_min) const
{
  std::vector<Particle_Source *>::const_iterator pit = _particle_sources.begin();
  for(; pit!=_particle_sources.end(); ++pit)
    dt = std::min(dt, (*pit)->limit_dt(time, dt, dt_min));


  std::vector<Light_Source *>::const_iterator lit = _light_sources.begin();
  for(; lit!=_light_sources.end(); ++lit)
    dt = std::min(dt, (*lit)->limit_dt(time, dt, dt_min));

  return dt;
}


bool FieldSource::request_serial_mesh() const
{
  std::vector<Light_Source *>::const_iterator lit = _light_sources.begin();
  for(; lit!=_light_sources.end(); ++lit)
    if( (*lit)->light_source_type() == "light_source_raytracing") return true;
  return false;
}


std::vector<Light_Source_From_File *> FieldSource::parse_spectrum_file(const std::string & filename, SimulationSystem & system, const Parser::Card &c)
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
          MESSAGE<<"ERROR at " << c.get_fileline() <<": Reading spectrum file error." << filename << std::endl; RECORD();
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
  }

  unsigned int n_sources = sources.size();

  //post process
  if(n_sources < 2)
  {
    MESSAGE<<"ERROR at " << c.get_fileline() <<": Spectrum file shoule contain at least 2 spectrums." << filename << std::endl; RECORD();
    genius_error();
  }

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
      MESSAGE<<"ERROR at " << c.get_fileline() <<": Wavelength in spectrum file must be strictly increasing." << filename << std::endl; RECORD();
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

  return sources;

}



void  FieldSource::SetWaveformUniform(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  assert( label!="" );

  double amplitude = c.get_real("amplitude", 1.0);
  double Tdelay = c.get_real("tdelay",0.0)*s;

  _waveforms[label] = new WaveformUniform(label, Tdelay, amplitude);
}


void  FieldSource::SetWaveformSin(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  assert( label!="" );

  double amplitude = c.get_real("amplitude", 1.0);
  double amplitude_offset = c.get_real("amplitude.offset", 0.0);
  double Tdelay = c.get_real("tdelay",0.0)*s;
  double freq = c.get_real("freq",0.0)*1.0/s;
  double alpha= c.get_real("alpha",0.0)*1.0/s;

  _waveforms[label] = new WaveformSin(label, Tdelay, amplitude_offset, amplitude, freq, alpha);
}


void FieldSource::SetWaveformGauss(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  assert( label!="" );

  double amplitude = c.get_real("amplitude", 1.0);
  double t0 = c.get_real("t0",0.0)*s;
  double tao = c.get_real("tau",1e-12)*s;

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
  double pr = c.get_real("pr",tr+tf+pw+pw)*s;
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



void  FieldSource::SetWaveformDoubleExp(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  assert( label!="" );

  double Tdelay = c.get_real("tdelay",0.0)*s;
  double trc = c.get_real("trc",1e-12)*s;
  double tfd = c.get_real("tfd",1e-9)*s;
  double tfc = c.get_real("tfc",1e-12)*s;
  double Alo = c.get_real("amplitude.low",0.0);
  double Ahi = c.get_real("amplitude.high",1.0);

  _waveforms[label] = new WaveformExponential(label, Tdelay, Alo, Ahi, trc, tfd, tfc);
}


void  FieldSource::SetWaveformShell(const Parser::Card &c)
{

  std::string label = c.get_string("id","");
  assert( label!="" );

  std::string filename = c.get_string("dll","");
  std::string funcname = c.get_string("function","");

  _waveforms[label] = new WaveformShell(label, filename, funcname, s);

}


void  FieldSource::SetWaveformWaveFile(const Parser::Card &c)
{

  std::string label = c.get_string("id","");
  assert( label!="" );

  std::string filename = c.get_string("file","");

  _waveforms[label] = new WaveformFile(label, filename);
}





