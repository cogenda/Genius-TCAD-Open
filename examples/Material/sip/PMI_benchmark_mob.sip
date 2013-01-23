%Module PMI_Benckmark_Mob 0

// {{{ %MappedType std::string
%MappedType std::string
{
%TypeHeaderCode
#include <string>
%End

%ConvertFromTypeCode
  if (!sipCpp)
  {
    Py_INCREF (Py_None);
    return Py_None;
  }
  return PyString_FromString(sipCpp->c_str());
%End

%ConvertToTypeCode
  if (sipIsErr == NULL)
    return (PyString_Check(sipPy));

  if (sipPy == Py_None)
  {
     *sipCppPtr = new std::string();
     return 1;
  }
  if (PyString_Check(sipPy))
  {
     *sipCppPtr = new std::string(PyString_AS_STRING(sipPy));
     return 1;
  }
  return 0;
%End

};
// }}}


class PMI_Benckmark_Mob    /NoDefaultCtors/
{
%TypeHeaderCode
#include "PMI_benchmark_mob.h"
%End

  public:
  /**
   * constructor, take material name and the name of mobility model
   */
  PMI_Benckmark_Mob(const std::string &path, const std::string &material, const std::string &model);

  ~PMI_Benckmark_Mob();

  /**
   * set the doping level, with unit cm^-3
   */
  void set_doping(double Na, double Nd);

  /**
   * set the mole fraction
   */
  void set_mole(double mole_x, double mole_y=0.0);

  /**
   * calculate electron mobility with:
   * carrier density (cm^-3)
   * vertical and parallel electrical field (V/cm)
   * temperature (K)
   */
  double mob_electron(double p, double n, double Ep, double Et, double T);

  /**
   * calculate hole mobility with:
   * carrier density (cm^-3)
   * vertical and parallel electrical field (V/cm)
   * temperature (K)
   */
  double mob_hole(double p, double n, double Ep, double Et, double T);

};
