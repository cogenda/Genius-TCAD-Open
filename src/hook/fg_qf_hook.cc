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
#include <ctime>
#include <cstdlib>
#include <iomanip>

#include "mesh_base.h"
#include "solver_base.h"
#include "fg_qf_hook.h"
#include "parallel.h"


#ifdef WINDOWS
  #include <io.h>      // for windows _access function
#else
  #include <unistd.h>  // for POSIX access function
#endif


/*
 * usage: HOOK Load=surface_recombination string<output>=file_name
 */

/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
FloatGateQfHook::FloatGateQfHook ( SolverBase & solver, const std::string & name, void * param)
  : Hook ( solver, name )
{
  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for ( std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
        parm_it != parm_list.end(); parm_it++ )
  {
    
    if(parm_it->name() == "floatgate" && parm_it->type() == Parser::STRING)
      float_gate=parm_it->get_string();
    
    if ( parm_it->name() == "tdelay" && parm_it->type() == Parser::REAL )
      tdelay=parm_it->get_real() * PhysicalUnit::s;
    if ( parm_it->name() == "tr" && parm_it->type() == Parser::REAL )
      tr=parm_it->get_real() * PhysicalUnit::s;
    if ( parm_it->name() == "tf" && parm_it->type() == Parser::REAL )
      tf=parm_it->get_real() * PhysicalUnit::s;
    if ( parm_it->name() == "pw" && parm_it->type() == Parser::REAL )
      pw=parm_it->get_real() * PhysicalUnit::s;
    if ( parm_it->name() == "pr" && parm_it->type() == Parser::REAL )
      pr=parm_it->get_real() * PhysicalUnit::s;
    
    if ( parm_it->name() == "q1" && parm_it->type() == Parser::REAL )
      q1=parm_it->get_real() * PhysicalUnit::C;
    if ( parm_it->name() == "q2" && parm_it->type() == Parser::REAL )
      q2=parm_it->get_real() * PhysicalUnit::C;
  }

}


/*----------------------------------------------------------------------
 * destructor, close file
 */
FloatGateQfHook::~FloatGateQfHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void FloatGateQfHook::on_init()
{}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void FloatGateQfHook::pre_solve()
{
  SimulationSystem & system = _solver.get_system();
  BoundaryConditionCollector * bcs = system.get_bcs();
  for(unsigned int n=0; n<bcs->n_bcs(); n++)
  {
    BoundaryCondition * bc = bcs->get_bc(n);
    if( bc->label() == float_gate && bc->bc_type() == ChargeIntegral )
      bc->scalar("qf") = qwave(SolverSpecify::clock);
  }
}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void FloatGateQfHook::post_solve()
{
 
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void FloatGateQfHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void FloatGateQfHook::on_close()
{

}


//---------------------------------------------------------------------
double FloatGateQfHook::qwave(const double t) const
{
    double _t = t;
    if(_t<tdelay)
      return q1;
    else
    {
      _t-=tdelay;
      while(_t>pr) _t-=pr;
      if(_t<tr)
        return q1+_t*(q2-q1)/tr;
      else if(_t<tr+pw)
        return q2;
      else if(_t<tr+pw+tf)
        return q2-(_t-tr-pw)*(q2-q1)/tf;
      else    return q1;
    }
}



#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new FloatGateQfHook ( solver, name, fun_data );
  }

}

#endif

