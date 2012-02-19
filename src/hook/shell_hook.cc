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

#include <string>
#include <cstdlib>

#include "solver_base.h"
#include "shell_hook.h"


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void ShellHook::on_init()
{
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    std::string exec_cmd = _command + " -init ";
    system(exec_cmd.c_str());
  }
}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void ShellHook::pre_solve()
{
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    std::string exec_cmd = _command + " -pre ";
    system(exec_cmd.c_str());
  }

}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void ShellHook::post_solve()
{
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    std::string exec_cmd = _command + " -post ";
    system(exec_cmd.c_str());
  }
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void ShellHook::post_iteration()
{
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    std::string exec_cmd = _command + " -postit ";
    system(exec_cmd.c_str());
  }
}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void ShellHook::on_close()
{
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    std::string exec_cmd = _command + " -close ";
    system(exec_cmd.c_str());
  }
}



#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new ShellHook(solver, name, fun_data );
  }
}

#endif

