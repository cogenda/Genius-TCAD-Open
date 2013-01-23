#ifndef __PMI_benckmark_optical_h__
#define __PMI_benckmark_optical_h__

#include <string>
#include "PMI.h"

class PMI_Benckmark_Optical
{
public:
  /**
   * constructor, take material name and the name of mobility model
   */
  PMI_Benckmark_Optical(const std::string &path, const std::string &material, const std::string &model);

  ~PMI_Benckmark_Optical();

  /**
   * set the doping level, with unit cm^-3
   */
  void set_doping(double Na, double Nd);

  /**
   * set the mole fraction
   */
  void set_mole(double mole_x, double mole_y);

  /**
   * refraction index
   */
  double n(double lambda);

  /**
   * refraction index
   */
  double k(double lambda);

private:

  void             * dll_file;

  PMIS_Optical     * optical;

  double           cm;
  double           s;
  double           V;
  double           C;
  double           K;

};

#endif
