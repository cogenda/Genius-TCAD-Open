class PMI_Benchmark_Mob    /NoDefaultCtors/
{
%TypeHeaderCode
#include "PMI_benchmark_mob.h"
%End

  public:
  /**
   * constructor, take material name and the name of mobility model
   */
  PMI_Benchmark_Mob(const std::string &path, const std::string &material, const std::string &model="Defalut");

  ~PMI_Benchmark_Mob();
  
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
   * set the dmin, with unit nm
   */
  void set_dmin(double d);

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
