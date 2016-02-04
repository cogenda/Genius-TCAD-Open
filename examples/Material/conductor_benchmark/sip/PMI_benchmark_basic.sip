class PMI_Benchmark_Basic    /NoDefaultCtors/
{
%TypeHeaderCode
#include "PMI_benchmark_basic.h"
%End

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
  
  /**
   * @return the current density [A/cm^2] under E [V/cm] and T [K]
   */
  double CurrentDensity (const double E, const double Tl) const;
};
