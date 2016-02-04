#ifndef __PMI_benckmark_optical_h__
#define __PMI_benckmark_optical_h__

#include <string>
#include "PMI.h"

class PMI_Benchmark_Optical
{
public:
  /**
   * constructor, take material name and the name of mobility model
   */
  PMI_Benchmark_Optical(const std::string &path, const std::string &material, const std::string &model="Defalut");

  ~PMI_Benchmark_Optical();

  /**
   * refraction index
   */
  double n(double lambda, double T);

  /**
   * refraction index
   */
  double k(double lambda, double T);

  /**
   * alpha 1/cm
   */
  double alpha(double lambda, double T);

private:

  void             * dll_file;

  PMIC_Optical     * optical;

  double           cm;
  double           s;
  double           V;
  double           C;
  double           K;

  double           um;

};

#endif
