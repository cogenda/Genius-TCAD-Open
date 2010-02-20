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

#ifndef __mole_analytic_h__
#define __mole_analytic_h__

#include "point_solver.h"
#include "parser.h"
#include "mole_fun.h"

/**
 * this solver compute mole fraction profile by anaytic expression inputted by user.
 * since it only involves "local node", it is considered as a PointSolver
 */
class MoleAnalytic : public PointSolver
{
public:
  /**
   *  constructor, do nothing
   */
  MoleAnalytic(SimulationSystem & system, Parser::InputParser & decks)
  : PointSolver(system), _decks(decks)
  {system.record_active_solver(this->solver_type());}

  /**
   * destructor, free the mole functions
   */
  ~MoleAnalytic()
  {
    for(size_t n=0; n<_mole_funs.size(); n++)
      delete _mole_funs[n];
    _mole_funs.clear();
  }

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type() const
  {return SolverSpecify::MOLE_ANALYTIC;}

  /**
   * we parse input deck for mole profile here
   */
  virtual int create_solver();

  /**
   * assign mole fraction to nodes in semiconductor region
   */
  virtual int solve();

  /**
   * do nothing
   */
  virtual int destroy_solver()
  {return 0;}


private:

 /**
  * since reading deck involves stack operate, we can not use const here
  */
 Parser::InputParser & _decks;

 /**
  * parse MOLE card, build mole anaytic function and
  * push into _mole_funs.
  */
 void set_mole_function_linear(const Parser::Card & );

 /**
  * parse PROFILE card with doping file, build doping profile and
  * push into _doping_data.
  */
 void set_mole_function_file(const Parser::Card & c);

 /**
  * @return the first mole fraction of given node
  */
 double mole_x(const std::string & region, const Node * node);

 /**
  * @return the second mole fraction of given node
  */
 double mole_y(const std::string & region, const Node * node);

 /**
  * the pointer vector to DopingFunction
  */
 std::vector<MoleFunction * >  _mole_funs;
 std::vector<std::pair<InterpolationBase *, InterpolationBase*> >  _mole_data;

};


#endif // #define __mole_analytic_h__
