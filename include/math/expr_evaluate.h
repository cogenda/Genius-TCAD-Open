#ifndef __expr_evalute_h__
#define __expr_evalute_h__

#include "expr_eval.h"

/**
 * evalute an expression in string
 */
class ExprEvalute
{
public:

  /**
   * general constructor
   */
  ExprEvalute(const std::string & expr);

  double eval(double x, double y, double z, double t);

  double operator () (double x, double y, double z, double t)
  { return eval(x,y,z,t); }

private:

  /**
   * all the variables
   */
  ExprEval::ValueList vlist;

  /**
   * all the functions
   */
  ExprEval::FunctionList flist;

  /**
   * parsed expression
   */
  ExprEval::Expression e;

};

#endif

