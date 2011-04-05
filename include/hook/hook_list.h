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

//  $Id: hook_list.h,v 1.4 2008/07/09 05:58:16 gdiso Exp $


#ifndef __hook_list_h__
#define __hook_list_h__

#include <map>

#include "hook.h"


/**
 * This is a list of hooks that are executed once at a time
 */
class HookList
{
public:

  HookList() {}

  /**
   * destructor, free all the hooks
   */
  ~HookList()
  {
    this->clear();
  }

  /**
   * add hook to hook list
   */
  void   add_hook(Hook * hook)
  { _hook_list.push_back(hook); }

  /**
   *   This is executed before the initialization of the solver
   */
  void on_init()
  {
    std::deque<Hook *>::iterator it;
    for (it=_hook_list.begin(); it!=_hook_list.end(); ++it)
      (*it)->on_init();
  }

  /**
   *   This is executed previously to each solution step.
   */
  void pre_solve()
  {
    std::deque<Hook *>::iterator it;
    for (it=_hook_list.begin(); it!=_hook_list.end(); ++it)
      (*it)->pre_solve();
  }

  /**
   *  This is executed after each solution step.
   */
  void post_solve()
  {
    std::deque<Hook *>::iterator it;
    for (it=_hook_list.begin(); it!=_hook_list.end(); ++it)
      (*it)->post_solve();
  }


  /**
   *  This is executed before each (nonlinear) iteration
   *  i.e. for analysis the condition number of jacobian matrix
   */
  void pre_iteration()
  {
    std::deque<Hook *>::iterator it;
    for (it=_hook_list.begin(); it!=_hook_list.end(); ++it)
      (*it)->pre_iteration();
  }

  /**
   *  This is executed after each (nonlinear) iteration
   *  i.e. for collecting convergence information or implementing various damping strategy
   */
  void post_iteration()
  {
    std::deque<Hook *>::iterator it;
    for (it=_hook_list.begin(); it!=_hook_list.end(); ++it)
      (*it)->post_iteration();
  }

  /**
   *  This is executed after each (nonlinear) iteration
   *  i.e. for collecting convergence information or implementing various damping strategy
   */
  void post_iteration(void * f, void * x, void * y, void * w, bool & change_y, bool &change_w)
  {
    std::deque<Hook *>::iterator it;
    for (it=_hook_list.begin(); it!=_hook_list.end(); ++it)
      (*it)->post_iteration(f, x, y, w, change_y, change_w);
  }

  /**
   * This is executed after the finalization of the solver
   */
  void on_close()
  {
    std::deque<Hook *>::iterator it;
    for (it=_hook_list.begin(); it!=_hook_list.end(); ++it)
      (*it)->on_close();

    this->clear();
  }

  /**
   * clear all the hooks
   */
  void clear()
  {
    std::deque<Hook *>::iterator it;
    for ( it=_hook_list.begin(); it!=_hook_list.end(); ++it)
      delete (*it);
    _hook_list.clear();
  }

private:

  std::deque<Hook *>  _hook_list;

};



#endif //define __hook_list_h__
