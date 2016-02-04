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

// c++ include
#include <set>
#include <iomanip>
#include <numeric>

// Local includes
#include "gdml_io.h"
#include "mesh_base.h"
#include "mesh_tools.h"
#include "boundary_info.h"
#include "simulation_system.h"
#include "boundary_condition_collector.h"
#include "material_define.h"
#include "parallel.h"


using PhysicalUnit::g;
using PhysicalUnit::um;
using PhysicalUnit::cm;




void GDMLIO::write_gdml_surface (const std::string& filename)
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  const MeshBase & mesh = system.mesh();

  genius_assert(mesh.mesh_dimension()==3);

  if(Genius::processor_id() == 0)
  {
    _out.open(filename.c_str(), std::ofstream::trunc);
    genius_assert(_out.good());

    //write GDML file header
    _out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    _out << "<gdml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
    << "xsi:noNamespaceSchemaLocation=\"http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd\">" << std::endl << std::endl;

    MeshTools::BoundingBox bbox = MeshTools::bounding_box(mesh);
    double x_c = 0.5*(bbox.second.x() + bbox.first.x())/um;
    double y_c = 0.5*(bbox.second.y() + bbox.first.y())/um;
    double z_c = 0.5*(bbox.second.z() + bbox.first.z())/um;
    double x_len = (bbox.second.x() - bbox.first.x())/um;
    double y_len = (bbox.second.y() - bbox.first.y())/um;
    double z_len = (bbox.second.z() - bbox.first.z())/um;
    //Point center(x_c, y_c, z_c);

    x_len = 2*(std::abs(x_c) + x_len);
    y_len = 2*(std::abs(y_c) + y_len);
    z_len = 2*(std::abs(z_c) + z_len);

    //pick boundary nodes, write into define/position
    std::vector<unsigned int>       el;
    std::vector<unsigned short int> sl;
    std::vector<short int>          il;
    system.mesh().boundary_info->build_side_list (el, sl, il);

    std::set<const Node *> _boundary_nodes;
    for(unsigned int n=0; n<sl.size(); n++)
    {
      const Elem * elem = mesh.elem(el[n]);
      AutoPtr<Elem> boundary_side = elem->build_side(sl[n]);
      for(unsigned int nn=0; nn<boundary_side->n_nodes(); ++nn)
        _boundary_nodes.insert(boundary_side->get_node(nn));
    }

    // begin defind section
    _out << "<define>" << std::endl;
    _out << std::setprecision(8) << std::scientific << std::right;

    for(std::set<const Node *>::iterator it=_boundary_nodes.begin(); it!=_boundary_nodes.end(); ++it)
      {
        //Point p((*it)->x()/um, (*it)->y()/um, (*it)->z()/um);
        _out << "  <position name=\"v" << (*it)->id() << "\" unit=\"micron\" "
        << "x=\"" << (*it)->x()/um <<"\" "
        << "y=\"" << (*it)->y()/um <<"\" "
        << "z=\"" << (*it)->z()/um <<"\"/>"
        <<std::endl;
      }

    _out << "  <position name=\"center\" "
    << "x=\"" << x_c << "\" "
    << "y=\"" << y_c << "\" "
    << "z=\"" << z_c << "\" "
    << "unit=\"micron\"/>" <<std::endl;



    // end defind section
    _out << "</define>" << std::endl << std::endl;

    //begin material section
    _out << "<materials>" << std::endl;
    _out << std::fixed << std::right;

    // find all the atoms we used
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      std::vector<Atom> atoms;
      std::vector<double> fraction;
      region->atom_fraction(atoms, fraction);

      for(unsigned int n=0; n<atoms.size(); ++n)
      {
        const Atom & atom = atoms[n];

        const std::string name = region->name() + '_' + atom.name;
        const std::string formula = atom.formula;
        int Z = atom.Z;
        double atom_value = atom.atom_value;

        // if this element has isotope
        if(atom.isotope_info.size())
        {
          for(unsigned int n=0; n<atom.isotope_info.size(); ++n)
          {
            const Atom::Isotope & isotope =  atom.isotope_info[n];
            _out << "  <isotope name=\"" << isotope.isotope_name <<"\" Z=\"" << Z <<"\" N=\"" << isotope.nucleons <<"\">" << std::endl;
            _out << "    <atom type=\"A\" value=\"" << isotope.atom_value << "\"/>" << std::endl;
            _out << "  </isotope>" <<std::endl;
          }

          _out << "  <element name=\"" << name <<"\" >" <<std::endl;
          for(unsigned int n=0; n<atom.isotope_info.size(); ++n)
          {
            const Atom::Isotope & isotope =  atom.isotope_info[n];
            _out << "    <fraction ref=\"" << isotope.isotope_name <<"\" n=\"" << isotope.fraction <<"\" />" <<std::endl;
          }
          _out << "  </element>" <<std::endl;
        }
        else
        {
          _out << "  <element name=\"" << name << "\" formula=\"" << formula <<"\" Z=\"" << Z << "\">" <<std::endl;
          _out << "  <atom value=\"" << atom_value <<"\" />" << std::endl;
          _out << "  </element>" <<std::endl;
        }
      }
    }


    // hydrogen for Vacuum
    _out << "  <element name=\"" << "Vacuum_Hydrogen" << "\" formula=\"" << "H" <<"\" Z=\"" << 1 << "\">" <<std::endl;
    _out << "  <atom value=\"" << 1.0 <<"\" />" << std::endl;
    _out << "  </element>" <<std::endl;

    // write all the materials
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);

      _out << "  <material name=\"" << region->name() + '_' + Material::FormatMaterialString(region->material()) <<"\" formula=\"" << region->material() <<"\">" <<std::endl;
      _out << "  <D value=\"" << region->get_density()/(g*pow(cm,-3)) <<"\" />" <<std::endl;

      std::vector<Atom> atoms;
      std::vector<double> fraction;
      region->atom_fraction(atoms, fraction);

      double molecule_value= 0;
      for(unsigned int n=0; n<atoms.size(); ++n)
      {
        const Atom & atom = atoms[n];
        molecule_value += fraction[n]*atom.atom_value;
      }

      for(unsigned int n=0; n<atoms.size(); ++n)
      {
        const Atom & atom = atoms[n];
        const std::string & atom_ref = region->name() + '_' + atom.name;
        _out << "  <fraction n=\""<< fraction[n]*atom.atom_value/molecule_value <<"\" ref=\"" << atom_ref <<"\" />" <<std::endl;
      }

      _out << "  </material> " <<std::endl;
    }

    // write material "vacuum", for volume "Word"
    _out << "  <material formula=\"Vacuum\" name=\"vacuum\" >" <<std::endl;
    _out << "  <D value=\"1e-25\" unit=\"g/cm3\"/>" <<std::endl;
    _out << "  <fraction n=\"1.0\" ref=\"Vacuum_Hydrogen\" />" <<std::endl;
    _out << "  </material>" <<std::endl;

    //end material section
    _out << "</materials>" << std::endl << std::endl;


    //begin solids section
    _out << "<solids>" << std::endl;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);

      _out << "  <tessellated name=\"Solid_" << region->name() << "\">" << std::endl;

      for(unsigned int n=0; n<sl.size(); n++)
      {
        const Elem * elem = mesh.elem(el[n]);
        if( elem->subdomain_id() != r) continue;

        AutoPtr<Elem> boundary_side = elem->build_side(sl[n]);

        switch(boundary_side->type())
        {
            case TRI3      :
            case TRI3_FVM  :
            _out << "    <triangular   "
            << "vertex1=\"v" << boundary_side->get_node(0)->id()<<"\" "
            << "vertex2=\"v" << boundary_side->get_node(1)->id()<<"\" "
            << "vertex3=\"v" << boundary_side->get_node(2)->id()<<"\" "
            << "type=\"ABSOLUTE\"/>" << std::endl;
            break;
            case QUAD4 :
            case QUAD4_FVM :
            _out << "    <quadrangular "
            << "vertex1=\"v" << boundary_side->get_node(0)->id()<<"\" "
            << "vertex2=\"v" << boundary_side->get_node(1)->id()<<"\" "
            << "vertex3=\"v" << boundary_side->get_node(2)->id()<<"\" "
            << "vertex4=\"v" << boundary_side->get_node(3)->id()<<"\" "
            << "type=\"ABSOLUTE\"/>" << std::endl;
            break;
            default:
            {
              std::cerr<<"Geometry face can only be triangular or quadrangular"<<std::endl;
              genius_error();
            }
        }
      }

      _out << "  </tessellated>" << std::endl;
    }

    //write solid world
    _out << "  <box name=\"world\" "
         << "x=\"" << x_len <<"\" "
         << "y=\"" << y_len <<"\" "
         << "z=\"" << z_len <<"\" "
         << "lunit=\"micron\" />" <<std::endl;


    //end solids section
    _out << "</solids>" << std::endl << std::endl;


    //begin structure section
    _out << "<structure>" << std::endl;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);

      _out << "  <volume name=\"" << std::string("Logic_") + region->name() << "\">" << std::endl;
      _out << "    <materialref ref=\"" << region->name() + '_' + Material::FormatMaterialString(region->material()) <<"\"/>" << std::endl;
      _out << "    <solidref ref=\"Solid_" << region->name() <<"\"/>" << std::endl;
      _out << "  </volume>" << std::endl;
    }

    _out << "  <volume name=\"World\" >" << std::endl;
    _out << "    <materialref ref=\"vacuum\"/>" << std::endl;
    _out << "    <solidref ref=\"world\"/>" << std::endl;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      _out << "    <physvol name=\"" << region->name() << "\">" << std::endl;
      _out << "       <volumeref ref=\"" << std::string("Logic_") + region->name() << "\"/>" << std::endl;
      _out << "    </physvol>" << std::endl;
    }
    _out << "  </volume>" << std::endl;


    //end structure section
    _out << "</structure>" << std::endl << std::endl;


    //write setup section
    _out << "<setup name=\"GeniusToGeant4\" version=\"1.0\">" << std::endl;
    _out << "  <world ref=\"World\"/>" << std::endl;
    _out << "</setup>" << std::endl << std::endl;

    _out << "</gdml>" << std::endl;

    _out.close();
  }
}







void GDMLIO::set_atoms()
{
  _atoms.clear();

  _atoms["H"  ] = Atom("Hydrogen", "H", 1, 1.00797);

  _atoms["B"  ] = Atom("Boron",    "B", 5, 10.811);
  _atoms["B"  ].add_isotope("B10", 10, 10.0, 0.199);
  _atoms["B"  ].add_isotope("B11", 11, 11.0, 0.801);

  _atoms["C"  ] = Atom("Carbon",   "C", 6, 12.01115);
  _atoms["N"  ] = Atom("Nitrogen", "N", 7, 14.01);
  _atoms["O"  ] = Atom("Oxygen",   "O", 8, 16.0);
  _atoms["F"  ] = Atom("Fluorine", "F", 9, 18.998);

  _atoms["Mg" ] = Atom("Magnesium",  "Mg", 12, 24.3050);
  _atoms["Al" ] = Atom("Aluminum",   "Al", 13, 26.98);
  _atoms["Si" ] = Atom("Silicon",    "Si", 14, 28.086);
  _atoms["P"  ] = Atom("Phosphorus", "P",  15, 30.9738);
  _atoms["S"  ] = Atom("Sulfur",     "S",  16, 32.065);
  _atoms["Cl" ] = Atom("Chlorine",   "Cl", 17, 35.453);

  _atoms["Ti" ] = Atom("Titanium",  "Ti", 22, 47.9);
  _atoms["Fe" ] = Atom("Iron",      "Fe", 26, 55.845);
  _atoms["Co" ] = Atom("Cobalt",    "Co", 27, 58.933);
  _atoms["Ni" ] = Atom("Nickle",    "Ni", 28, 58.69);
  _atoms["Cu" ] = Atom("Copper",    "Cu", 29, 63.54);
  _atoms["Zn" ] = Atom("Zinc",      "Zn", 30, 65.409);
  _atoms["Ga" ] = Atom("Gallium",   "Ga", 31, 69.72);
  _atoms["Ge" ] = Atom("Germanium", "Ge", 32, 72.59);
  _atoms["As" ] = Atom("Arsenic",   "As", 33, 74.922);
  _atoms["Se" ] = Atom("Selenium",  "Se", 34, 78.96);

  
  _atoms["Ag" ] = Atom("Silver",    "Ag", 47, 107.87);
  _atoms["Cd" ] = Atom("Cadmium",   "Cd", 48, 112.40);
  _atoms["In" ] = Atom("Indium",    "In", 49, 114.82);
  _atoms["Sn" ] = Atom("Tin",       "Sn", 50, 118.71);
  _atoms["Sb" ] = Atom("Antimony",  "Sb", 51, 121.75);
  _atoms["Te" ] = Atom("Tellurium", "Te", 52, 127.6);

  _atoms["Hf" ] = Atom("Hafnium",   "Hf", 72, 178.49);
  _atoms["Ta" ] = Atom("Tantalum",  "Ta", 73, 180.95);
  _atoms["W"  ] = Atom("Wolfram",   "W",  74, 183.85);
  _atoms["Au" ] = Atom("Gold",      "Au", 79, 196.97);
  _atoms["Hg" ] = Atom("Mercury",   "Hg", 80, 200.59);
  _atoms["Pb" ] = Atom("Lead",      "Pb", 82, 207.2);

}


double GDMLIO::atom_value(const std::string &name)
{
  if( _atoms.find(name) != _atoms.end() )
    return _atoms.find(name)->second.atom_value;
  else
  {
    for(std::map<std::string, Atom>::const_iterator it = _atoms.begin(); it!=_atoms.end(); ++it)
    {
      if( it->second.name == name )
        return it->second.atom_value;
    }
  }

  return 1.0;
}

