class PMI_Benchmark_ImpactIonization    /NoDefaultCtors/
{
%TypeHeaderCode
#include "PMI_benchmark_avalanche.h"
%End

public:
  /**
   * constructor, take material name and the name of thermal model
   */
  PMI_Benchmark_ImpactIonization(const std::string &path, const std::string &material, const std::string &model="Defalut");

  ~PMI_Benchmark_ImpactIonization();

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
   * @return the electron generation rate [1/cm] for DDM simulation
   */
  double  ElecGenRate (const double Tl,const double Ep,const double Eg) const;

  /**
   * @return the hole generation rate [1/cm] for DDM simulation
   */
  double HoleGenRate (const double Tl,const double Ep,const double Eg) const;

};
