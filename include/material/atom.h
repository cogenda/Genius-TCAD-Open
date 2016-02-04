#ifndef __atom_table_h__
#define __atom_table_h__

#include <string>
#include <vector>


struct Atom
{
  /// atom name, only for converinent
  std::string name;
  /// atom formula
  std::string formula;
  /// atom number
  int         Z;
  /// avarage atom value
  double      atom_value;

  /// isotope
  struct Isotope
  {
    /// name of isotope
    std::string isotope_name;
    /// number of isotope's nucleons
    int         nucleons;
    /// atom mass
    double      atom_value;
    /// fraction of this isotope
    double      fraction;
  };

  std::vector<Isotope>  isotope_info;

  Atom() { Z=0; atom_value=0.0; }

  Atom(const std::string & _name, const std::string & _formula, int _Z, double _atom_value)
      : name(_name), formula(_formula), Z(_Z), atom_value(_atom_value)
  {}

  void add_isotope(const std::string &name, int n, double v, double f)
  {
    Isotope isotope;
    isotope.isotope_name = name;
    isotope.nucleons     = n;
    isotope.atom_value   = v;
    isotope.fraction     = f;
    isotope_info.push_back(isotope);
  }
};



#endif
