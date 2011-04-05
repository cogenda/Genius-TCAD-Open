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

//  $Id: doping_analytic.h,v 1.7 2008/07/09 05:58:16 gdiso Exp $

#ifndef __doping_analytic_h__
#define __doping_analytic_h__

#include "point_solver.h"
#include "doping_fun.h"
#include "parser.h"

/**
 * this solver compute doping profile by anaytic expression inputted by user
 * since it only involves "local node", it is considered as a PointSolver
 */
class DopingAnalytic : public PointSolver
{
public:
  /**
   *  constructor, do nothing
   */
  DopingAnalytic(SimulationSystem & system, Parser::InputParser & decks)
  : PointSolver(system), _decks(decks)
  {system.record_active_solver(this->solver_type());}

  /**
   * destructor, free the doping functions
   */
  ~DopingAnalytic()
  {
    for(size_t n=0; n<_doping_funs.size(); n++)
      delete _doping_funs[n];
    _doping_funs.clear();
    for(size_t n=0; n<_doping_data.size(); n++)
      delete _doping_data[n].second;
    _doping_data.clear();
    for(std::map<std::string,DopingFunction *>::iterator it=_custom_profile_funs.begin();
        it!=_custom_profile_funs.end(); it++)
      delete it->second;
    _custom_profile_funs.clear();
  }

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type() const
  {return SolverSpecify::DOPING_ANALYTIC;}

  /**
   * we parse input deck for doping profile here
   */
  virtual int create_solver();

  /**
   * virtual function, destroy the solver
   */
  virtual int destroy_solver();

  /**
   * do pre-process before each solve action
   */
  virtual int pre_solve_process(bool =true);

  /**
   * do post-process after each solve action
   */
  virtual int post_solve_process();

  /**
   * assign doping profile to semiconductor region
   */
  virtual int solve();

private:

 /**
  * since reading deck involves stack operate, we can not use const here
  */
 Parser::InputParser & _decks;

 /**
  * parse PROFILE card with uniform doping, build uniform doping function and
  * push into _doping_funs.
  */
 void set_doping_function_uniform(const Parser::Card & c);

 /**
  * parse PROFILE card with analytic doping, build analytic doping function and
  * push into _doping_funs.
  */
 void set_doping_function_analytic(const Parser::Card & c);


 enum Profile_Data_Axes {AXES_X, AXES_Y, AXES_Z, AXES_XY, AXES_XZ, AXES_YZ, AXES_XYZ};

 /**
  * parse PROFILE card with doping file, build doping profile and
  * push into _doping_data.
  */
 void set_doping_function_file(const Parser::Card & c);

 double _do_doping_interp(int i, const Node *node, const std::string &msg=std::string());

 /**
  * @return the acceptor concentration of node
  */
 double doping_Na(const Node * node);

 /**
  * @return the donor concentration of node
  */
 double doping_Nd(const Node * node);

 /**
  * the pointer vector to DopingFunction
  */
 std::vector<DopingFunction * >  _doping_funs;
 std::vector<std::pair<int, InterpolationBase * > > _doping_data;
 std::map<std::string,DopingFunction * >  _custom_profile_funs;

};


#endif // #define __doping_analytic_h__
