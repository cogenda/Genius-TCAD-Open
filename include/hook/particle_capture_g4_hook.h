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


#ifndef __particle_capture_g4_hook_h__
#define __particle_capture_g4_hook_h__


#ifdef HAVE_AMS

#include "hook.h"
#include "ams.h"

class DoseRate;

/**
 * interface to G4 particle simulator to get position of trapped particle
 * and RIC enhancement
 */
class ParticleCaptureG4Hook : public Hook
{

public:

  ParticleCaptureG4Hook(SolverBase & solver, const std::string & name, void *);

  virtual ~ParticleCaptureG4Hook();

  /**
   *   This is executed before the initialization of the solver
   */
  virtual void on_init();

  /**
   *   This is executed previously to each solution step.
   */
  virtual void pre_solve();

  /**
   *  This is executed after each solution step.
   */
  virtual void post_solve();

  /**
   *  This is executed after each (nonlinear) iteration
   */
  virtual void post_iteration();

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close();

private:

  /**
   * G4 particle type, only support electron at present
   */
  std::string _particle;

  /**
   * the weight of the particle, a weight 1e4 means one electron takes 1e4 charge
   */
  double _weight;

  /**
   * particle number per cm^2 per second
   */
  double    _particle_flux;

  /**
   * particle source area
   */
  double    _source_area;

  /**
   * each G4 simulation will run this number of particles
   */
  double    _particle_num_per_run;

  double    _lateral_char;

  /**
   * resolution of background mesh
   */
  double    _resolution;

  /**
   * beam particle number, G4 simulator run _beam_num of particles one time
   */
  int    _beam_num;

  /**
   * if the particle flux is constant, do G4 simulation only once
   */
  bool   _const_flux;


  bool   _first_g4_calculation_flag;


  std::string   _dose_rate;

private:

  /// AMS hander

  std::string host;

  int err;

  char *msg;

  AMS_Comm control;

  AMS_Comm data;

  void send_command(const std::string &cmd);

  void load_events();



  /// particle gen and dose rate

  void process_particle_gen();

  DoseRate * dose_rate;

  void process_particle_dose_nodal();

  void process_particle_dose_octree();

  /// particle data


  int EventNumber;

  std::vector<double> ParticleEnergy;

  std::vector<double> ParticleBeginPosition;

  std::vector<double> ParticleEndPosition;

  // the step number of each event
  std::vector<int> EventSteps;

  // end point of each step
  std::vector<double> ParticleStep;

  std::vector<double> ParticleWeight;

  std::vector<double> StepEnergy;


  struct track_t
  {
    Point        start;
    Point        end;
    double       energy;
    double       lateral_char;
  };

  std::vector<track_t>  tracks;

  void build_tracks();
};


#endif


#endif
