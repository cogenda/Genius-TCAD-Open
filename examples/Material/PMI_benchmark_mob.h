#ifndef __PMI_benckmark_mob_h__
#define __PMI_benckmark_mob_h__

#include <string>
#include "PMI.h"

class PMI_Benckmark_Mob
{
public:
  /**
   * constructor, take material name and the name of mobility model
   */
  PMI_Benckmark_Mob(const std::string &material, const std::string &model);

  ~PMI_Benckmark_Mob();

  /**
   * set the doping level, with unit cm^-3
   */
  void set_doping(double Na, double Nd);

  /**
   * set the mole fraction
   */
  void set_mole(double mole_x, double mole_y);

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

private:

  void             * dll_file;

  PMIS_Mobility    * mob;

  double           cm;
  double           s;
  double           V;
  double           C;
  double           K;

};

#endif
