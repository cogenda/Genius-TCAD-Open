/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/


#ifndef __gdml_io_h__
#define __gdml_io_h__

// C++ includes
#include <fstream>

// Local includes
#include "genius_common.h"
#include "field_input.h"
#include "field_output.h"
#include "simulation_system.h"


// Forward declarations
class MeshBase;
class Elem;


/**
 * This class implements writing geometry info to
 * Geometry Description Markup Language (GDML) file
 * This file format is used by GEANT4 for build its detector.
 */
class GDMLIO : public FieldInput<SimulationSystem>,
      public FieldOutput<SimulationSystem>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  GDMLIO (SimulationSystem& system);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  GDMLIO (const SimulationSystem& system, bool surface=true);

  /**
   * This method implements reading a geometry from GDML
   */
  virtual void read (const std::string& )
  { genius_error(); }

  /**
   * This method implements writing geometry to GDML
   */
  virtual void write (const std::string& file)
  {
    if(_surface)
      write_gdml_surface(file);
    else
      write_gdml_body(file);
  }

private:

  std::ofstream _out;

  bool _surface;

  void write_gdml_surface(const std::string &file);

  void write_gdml_body(const std::string &file);

  std::map<unsigned int, double> Na, Nd;

  double Na_max, Nd_max;

  void cell_average_doping(const Elem *elem, double &, double &);

private:
  struct Atom
  {
    std::string name;
    std::string formula;
    int         Z;
    /// avarage atom value
    double      atom_value;

    /// isotope
    struct Isotope
    {
      std::string isotope_name;
      int         nucleons;
      double      atom_value;
      double      fraction;
    };

    std::vector<Isotope>  isotope_info;

    Atom() {}

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

  std::map<std::string, Atom> _atoms;

  /**
   * fill atoms map
   */
  void set_atoms();

  double atom_value(const std::string &);

};



// ------------------------------------------------------------
// GDMLIO inline members
inline
GDMLIO::GDMLIO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system)
{
}



inline
GDMLIO::GDMLIO (const SimulationSystem& system, bool surface) :
    FieldOutput<SimulationSystem>(system), _surface(surface)
{
}



#endif // #define __gdml_io_h__

