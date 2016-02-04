#ifndef __PMI_benckmark_basic_h__
#define __PMI_benckmark_basic_h__

#include <string>
#include "PMI.h"

class PMI_Benchmark_Basic
{
public:
  /**
   * constructor, take material name and the name of basic parameters
   */
  PMI_Benchmark_Basic(const std::string &path, const std::string &material, const std::string &model="Defalut");

  ~PMI_Benchmark_Basic();

  /**
   * calibrate real parameter
   */
  bool calibrate_real_parameter(const std::string & var_name, double var_value);

  /**
   * calibrate string parameter
   */
  bool calibrate_string_parameter(const std::string & var_name, const std::string &var_value);

  /**
   * set the doping level, with unit cm^-3
   */
  void set_doping(double Na, double Nd);

  /**
   * set the mole fraction
   */
  void set_mole(double mole_x, double mole_y);

  
  /**
   * @return the mass density [g cm^-3] of material
   */
  double Density       (const double Tl) const;

  /**
   * @return the \p relative \p permittivity of material
   */
  double Permittivity  ()                      const;

  /**
   * @return the \p relative \p permeability of material
   */
  double Permeability  ()                      const;

  /**
   * @return the affinity energy [eV] of material
   */
  double Affinity      (const double Tl) const;


private:

  void             * dll_file;

  PMIS_BasicParameter *basic;

  double           cm;
  double           s;
  double           V;
  double           C;
  double           K;

  double           eV;
  double           m;
  double           J;
  double           kg;
  double           g;
};

#endif
