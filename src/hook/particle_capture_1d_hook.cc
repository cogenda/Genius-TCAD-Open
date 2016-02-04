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

#include <fstream>

#include "solver_base.h"
#include "mesh_base.h"
#include "simulation_system.h"
#include "point_locator_base.h"
#include "particle_capture_1d_hook.h"
#include "parallel.h"



using PhysicalUnit::um;
using PhysicalUnit::mm;
using PhysicalUnit::cm;
using PhysicalUnit::J;
using PhysicalUnit::rad;
using PhysicalUnit::kg;
using PhysicalUnit::s;


ParticleCapture1DHook::ParticleCapture1DHook(SolverBase & solver, const std::string & name, void *param)
  : Hook ( solver, name ), _particle("e-")
{
  Point base_point;
  Point base_dir;
  
  std::string data_file;
  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for(std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
      parm_it != parm_list.end(); parm_it++)
  {

    if(parm_it->name() == "base.point.x" && parm_it->type() == Parser::REAL)
      base_point[0] = parm_it->get_real()*um;
    
    if(parm_it->name() == "base.point.y" && parm_it->type() == Parser::REAL)
      base_point[1] = parm_it->get_real()*um;

    if(parm_it->name() == "base.point.z" && parm_it->type() == Parser::REAL)
      base_point[2] = parm_it->get_real()*um;
    
    
    if(parm_it->name() == "base.dir.x" && parm_it->type() == Parser::REAL)
      base_dir[0] = parm_it->get_real();
    
    if(parm_it->name() == "base.dir.y" && parm_it->type() == Parser::REAL)
      base_dir[1] = parm_it->get_real();

    if(parm_it->name() == "base.dir.z" && parm_it->type() == Parser::REAL)
      base_dir[2] = parm_it->get_real();
    
    if(parm_it->name() == "data.file" && parm_it->type() == Parser::STRING)
      data_file = parm_it->get_string();
    
  }
  
  _base = Plane(base_point, base_dir);
  
  build(data_file);

}


ParticleCapture1DHook::~ParticleCapture1DHook()
{}


void ParticleCapture1DHook::build(const std::string &file)
{
  std::vector<double> depth;
  std::vector<double> gen;
  std::vector<double> dose;
  
  if(Genius::processor_id() == 0)
  {
    std::ifstream fin(file.c_str());
  
    while(!fin.eof())
    {
      double d, g, r;
      fin >>d >> g >> r;
      depth.push_back(d*um);
      gen.push_back(g/(cm*cm*cm)/s);
      dose.push_back(r*rad/s);
    }
  
  }
  
  Parallel::broadcast(depth);
  Parallel::broadcast(gen);
  Parallel::broadcast(dose);
  _gen.set(depth, gen);
  _dose.set(depth, dose);
}




/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void ParticleCapture1DHook::on_init()
{

}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void ParticleCapture1DHook::pre_solve()
{
  double dt = SolverSpecify::dt;
  if(!SolverSpecify::TimeDependent)
    dt = 1*PhysicalUnit::s;

  SimulationSystem & system = _solver.get_system();
  for(unsigned int i=0; i<system.n_regions(); i++)
  {
    SimulationRegion * region = system.region(i);

    if(region->type() == InsulatorRegion && region->material() != "Air")
    {
      SimulationRegion::local_node_iterator it = region->on_local_nodes_begin();
      SimulationRegion::local_node_iterator it_end = region->on_local_nodes_end();
      for(; it!=it_end; ++it)
      {
        FVM_Node * fvm_node = *it;
        Point p = *(fvm_node->root_node());
        double depth = _base.signed_distance(p);
        fvm_node->node_data()->Field_G()  = _gen.evaluate(depth);
        fvm_node->node_data()->DoseRate() = _dose.evaluate(depth);
      }
    }
  }
}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void ParticleCapture1DHook::post_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void ParticleCapture1DHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void ParticleCapture1DHook::on_close()
{

}





#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new ParticleCapture1DHook ( solver, name, fun_data );
  }

}

#endif


