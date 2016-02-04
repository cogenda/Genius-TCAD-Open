#ifndef __PMI_benckmark_thermal_h__
#define __PMI_benckmark_thermal_h__

#include <string>
#include "PMI.h"

class PMI_Benchmark_Thermal
{
public:
  /**
   * constructor, take material name and the name of thermal model
   */
  PMI_Benchmark_Thermal(const std::string &path, const std::string &material, const std::string &model);

  ~PMI_Benchmark_Thermal();

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
   * @return the heat capacity [J/g/K] of the material
   */
  double HeatCapacity  (const double Tl) const;

  /**
   * @return the heat conduction [W/cm/K] of the material
   */
  double HeatConduction(const double Tl) const;

private:

  void             * dll_file;

  PMIS_Thermal     * thermal;

  double           cm;
  double           s;
  double           V;
  double           C;
  double           K;

  double           m;
  double           J;
  double           kg;
  double           g;

};

#endif
