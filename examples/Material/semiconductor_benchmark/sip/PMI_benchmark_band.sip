class PMI_Benchmark_Band    /NoDefaultCtors/
{
%TypeHeaderCode
#include "PMI_benchmark_band.h"
%End

  public:
  /**
   * constructor, take material name and the name of mobility model
   */
  PMI_Benchmark_Band(const std::string &path, const std::string &material, const std::string &model="Defalut");


  ~PMI_Benchmark_Band();

  /**
   * set the doping level, with unit cm^-3
   */
  void set_doping(double Na, double Nd);

  /**
   * set the mole fraction
   */
  void set_mole(double mole_x, double mole_y);

  /**
   * @return band gap of semiconductor
   */
  double   Eg(const double Tl);

  /**
   * @return band gap narrowing due to heavy doping
   */
  double EgNarrow       (const double p, const double n, const double Tl);

  /**
   * @return effective density of states in the conduction band
   */
  double Nc             (const double Tl);

  /**
   * @return effective density of states in the valence band
   */
  double Nv             (const double Tl);


  /**
   * @return intrinsic carrier concentration
   */
  double ni            ( const double Tl);

  /**
   * @return effective intrinsic carrier concentration
   */
  double nie            (const double p, const double n, const double Tl);


  /**
   * @return direct Recombination rate
   */
  double CDIR           (const double Tl);


  /**
   * @return electron lift time in SHR Recombination
   */
  double TAUN           (const double Tl);

  /**
   * @return hole lift time in SHR Recombination
   */
  double TAUP           (const double Tl);


  /**
   * @return electron Auger Recombination rate
   */
  double AUGERN           (const double p, const double n, const double Tl);

  /**
   * @return hole Auger Recombination rate
   */
  double AUGERP           (const double p, const double n, const double Tl);

};
