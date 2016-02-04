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
#include "parser.h"

class DopingFunction;

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
  DopingAnalytic(SimulationSystem & system, Parser::InputParser & decks);

  /**
   * destructor, free the doping functions
   */
  ~DopingAnalytic();

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
  * parse PROFILE card with linear doping, build linear doping function and
  * push into _doping_funs.
  */
 void set_doping_function_linear(const Parser::Card & c);

 /**
  * parse PROFILE card with analytic doping, build analytic doping function and
  * push into _doping_funs.
  */
 void set_doping_function_analytic(const Parser::Card & c);

 /**
  * parse PROFILE card with analytic2 doping (rectangle mask), build analytic doping function and
  * push into _doping_funs.
  */
 void set_doping_function_analytic2_rec_mask(const Parser::Card & c);

 /**
  * parse PROFILE card with analytic2 doping (poly mask) , build analytic doping function and
  * push into _doping_funs.
  */
 void set_doping_function_analytic2_poly_mask(const Parser::Card & c);


 enum Profile_Data_Axes {AXES_X, AXES_Y, AXES_Z, AXES_XY, AXES_XZ, AXES_YZ, AXES_XYZ};

 /**
  * parse PROFILE card with doping file, build doping profile and
  * push into _doping_data.
  */
 void set_doping_function_file(const Parser::Card & c);


 
 /**
  * doping function struct  
  */ 
 struct DopingFunction_t
 {
   std::string label;
   
   /**
    * doping function 
    */ 
   DopingFunction * df;
   
  /**
   *  the applied regions
   */ 
   std::string region;
 };
 
 /**
  * doping data struct  
  */ 
 struct DopingData_t
 {
   std::string label;
   
   /**
    * doping function 
    */ 
   std::pair<int, InterpolationBase * > df;
      
   /**
    *  the applied regions
    */ 
   std::string region;
 };


 std::vector<DopingFunction_t>  _custom_profile_funs;
 std::vector<DopingData_t>  _custom_profile_data;
 
 void _doping_function_apply(DopingFunction * df, const std::string &region_app);
 void _custom_profile_function_apply(const std::string &ion, DopingFunction * df, const std::string &region_app);
 
 
 void _doping_data_apply(std::pair<int, InterpolationBase * > df, const std::string &region_app);
 void _custom_profile_data_apply(const std::string &ion, std::pair<int, InterpolationBase * > df, const std::string &region_app);
 
  
 double _do_data_interp(std::pair<int, InterpolationBase * > df, const Node *node, const std::string &msg=std::string());

};


#endif // #define __doping_analytic_h__
