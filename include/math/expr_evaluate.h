#ifndef __expr_evalute_h__
#define __expr_evalute_h__

#include "expr_eval.h"


/**
 * evalute an constante expression in string
 */
class ConstanteExprEvalute
{
public:

  /**
   * constante expression constructor
   */
  ConstanteExprEvalute(const std::string & expr);

  /**
     * evalute constante expression
   */
  double eval() ;

  /**
     * evalute constante expression
   */
  double operator () ()
    { return eval(); }

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



/**
 * evalute an expression in string
 */
class ExprEvalute
{
public:

  /**
   * expression with independent variable constructor
   */
  ExprEvalute(const std::string & expr);

  /**
   * evalute an expression, take coordinate(x, y, z) and time as independent variable
   */
  double eval(double x, double y, double z, double t);

  /**
   * evalute an expression, take coordinate(x, y, z) and time as independent variable
   */
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

