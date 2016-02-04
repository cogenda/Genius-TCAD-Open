class PMI_Benchmark_Optical    /NoDefaultCtors/
{
%TypeHeaderCode
#include "PMI_benchmark_optical.h"
%End

  public:
  /**
   * constructor, take material name and the name of mobility model
   */
  PMI_Benchmark_Optical(const std::string &path,const std::string &material, const std::string &model="Defalut");

  ~PMI_Benchmark_Optical();

  /**
   * set the doping level, with unit cm^-3
   */
  void set_doping(double Na, double Nd);

  /**
   * set the mole fraction
   */
  void set_mole(double mole_x, double mole_y=0.0);

  /**
   * refraction index
   */
  double n(double lambda, double T);

  /**
   * refraction index
   */
  double k(double lambda, double T);

  double alpha(double lambda, double T);
};
