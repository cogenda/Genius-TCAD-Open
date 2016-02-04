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
   * refraction index
   */
  double n(double lambda, double T);

  /**
   * refraction index
   */
  double k(double lambda, double T);

  double alpha(double lambda, double T);
};
