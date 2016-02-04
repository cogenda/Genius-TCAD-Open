#include "expr_evaluate.h"
#include "physical_unit.h"

ConstanteExprEvalute::ConstanteExprEvalute(const std::string & expr)
{
  // use default function set
  flist.AddDefaultFunctions();
  e.SetFunctionList(&flist);

  // add Math constant 'E' and 'PI'
  vlist.AddDefaultValues();

  // unit and physical conatant
  vlist.Add("cm",   PhysicalUnit::cm, true);
  vlist.Add("s",    PhysicalUnit::s,  true);
  vlist.Add("V",    PhysicalUnit::V,  true);
  vlist.Add("C",    PhysicalUnit::C,  true);
  vlist.Add("K",    PhysicalUnit::K,  true);

  vlist.Add("m",    PhysicalUnit::m,  true);
  vlist.Add("nm",   PhysicalUnit::nm, true);
  vlist.Add("um",   PhysicalUnit::um, true);
  vlist.Add("mm",   PhysicalUnit::mm, true);
  vlist.Add("J",    PhysicalUnit::J,  true);
  vlist.Add("W",    PhysicalUnit::W,  true);
  vlist.Add("kg",   PhysicalUnit::kg, true);
  vlist.Add("g",    PhysicalUnit::g,  true);
  vlist.Add("eV",   PhysicalUnit::eV, true);
  vlist.Add("ps",   PhysicalUnit::ps, true);
  vlist.Add("A",    PhysicalUnit::A,  true);
  vlist.Add("mA",   PhysicalUnit::mA, true);

  vlist.Add("kb",   PhysicalUnit::kb, true);
  vlist.Add("q",    PhysicalUnit::e,  true);
  vlist.Add("me",   PhysicalUnit::me, true);
  vlist.Add("eps0", PhysicalUnit::eps0, true);
  vlist.Add("mu0",  PhysicalUnit::mu0,  true);
  vlist.Add("h",    PhysicalUnit::h,    true);
  vlist.Add("hbar", PhysicalUnit::hbar, true);

  e.SetValueList(&vlist);

  e.Parse(expr);
}




double ConstanteExprEvalute::eval()
{
  return e.Evaluate();
}


ExprEvalute::ExprEvalute(const std::string & expr)
{
  // use default function set
  flist.AddDefaultFunctions();
  e.SetFunctionList(&flist);

  // add Math constant 'E' and 'PI'
  vlist.AddDefaultValues();

  // user provide variable
  vlist.Add("x");
  vlist.Add("y");
  vlist.Add("z");
  vlist.Add("t");

  // unit and physical conatant
  vlist.Add("cm",   PhysicalUnit::cm, true);
  vlist.Add("s",    PhysicalUnit::s,  true);
  vlist.Add("V",    PhysicalUnit::V,  true);
  vlist.Add("C",    PhysicalUnit::C,  true);
  vlist.Add("K",    PhysicalUnit::K,  true);

  vlist.Add("m",    PhysicalUnit::m,  true);
  vlist.Add("nm",   PhysicalUnit::nm, true);
  vlist.Add("um",   PhysicalUnit::um, true);
  vlist.Add("mm",   PhysicalUnit::mm, true);
  vlist.Add("J",    PhysicalUnit::J,  true);
  vlist.Add("W",    PhysicalUnit::W,  true);
  vlist.Add("kg",   PhysicalUnit::kg, true);
  vlist.Add("g",    PhysicalUnit::g,  true);
  vlist.Add("eV",   PhysicalUnit::eV, true);
  vlist.Add("ps",   PhysicalUnit::ps, true);
  vlist.Add("A",    PhysicalUnit::A,  true);
  vlist.Add("mA",   PhysicalUnit::mA, true);

  vlist.Add("kb",   PhysicalUnit::kb, true);
  vlist.Add("q",    PhysicalUnit::e,  true);
  vlist.Add("me",   PhysicalUnit::me, true);
  vlist.Add("eps0", PhysicalUnit::eps0, true);
  vlist.Add("mu0",  PhysicalUnit::mu0,  true);
  vlist.Add("h",    PhysicalUnit::h,    true);
  vlist.Add("hbar", PhysicalUnit::hbar, true);

  e.SetValueList(&vlist);

  e.Parse(expr);
}



double ExprEvalute::eval(double x, double y, double z, double t)
{
  //assign variable value to the expr
  *(vlist.GetAddress("x")) = x;
  *(vlist.GetAddress("y")) = y;
  *(vlist.GetAddress("z")) = z;
  *(vlist.GetAddress("t")) = t;

  return e.Evaluate();
}

