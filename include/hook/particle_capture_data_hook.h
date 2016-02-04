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


#ifndef __particle_capture_data_hook_h__
#define __particle_capture_data_hook_h__

#include <vector>

#include "hook.h"

class DoseRate;

/**
 * load G4 particle simulation data to get position of trapped particle
 * and RIC enhancement
 */
class ParticleCaptureDataHook : public Hook
{

public:

  ParticleCaptureDataHook(SolverBase & solver, const std::string & name, void *);

  virtual ~ParticleCaptureDataHook();

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
   * file name
   */
  std::string _track_data_file;


  /**
   * file name
   */
  std::string _particle_data_file;


  double _fraction;

  /**
   * the weight of the particle, a weight 1e4 means one electron takes 1e4 charge
   */
  double _weight;

  /**
   * resolution of background mesh
   */
  double    _resolution;


  bool   _first_g4_calculation_flag;


private:

  /// particle gen and dose rate

  void build_particles();
  bool build_particles_block(std::istream * in, unsigned int block_size);

  void process_particle_gen();



  void build_tracks();
  bool build_tracks_block(std::istream * in, unsigned int block_size);

  void process_particle_dose_octree();





  DoseRate * dose_rate;

  std::vector<double> elem_deposite;

};



#endif
