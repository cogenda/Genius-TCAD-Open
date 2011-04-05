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
#include "parallel.h"


using PhysicalUnit::g;
using PhysicalUnit::um;
using PhysicalUnit::cm;




void GDMLIO::write_gdml_surface (const std::string& filename)
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  const MeshBase & mesh = system.mesh();

  genius_assert(mesh.mesh_dimension()==3);

  set_atoms();

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
    std::set<std::string> atom_refs;
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      std::vector<std::string> atoms;
      std::vector<double> fraction;
      region->atom_fraction(atoms, fraction);

      for(unsigned int n=0; n<atoms.size(); ++n)
      {
        genius_assert(_atoms.find(atoms[n])!=_atoms.end());
        atom_refs.insert(atoms[n]);
      }
    }
    // two atoms for Air
    atom_refs.insert("Nitrogen");
    atom_refs.insert("Oxygen");

    // write all the atoms
    for(std::set<std::string>::iterator it=atom_refs.begin(); it!=atom_refs.end(); ++it)
      {
        std::map<std::string, Atom>::iterator atom_it = _atoms.find(*it);
        write_gdml_element(atom_it->first);
      }

    // write all the materials
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);

      _out << "  <material name=\"" << region->name() + "_" + region->material() <<"\" formula=\"" << region->material() <<"\">" <<std::endl;
      _out << "  <D value=\"" << region->get_density(system.T_external())/(g*pow(cm,-3)) <<"\" />" <<std::endl;

      std::vector<std::string> atoms;
      std::vector<double> fraction;
      region->atom_fraction(atoms, fraction);

      double molecule_value= 0;
      for(unsigned int n=0; n<atoms.size(); ++n)
      {
        const std::string & atom_name = atoms[n];
        genius_assert(_atoms.find(atom_name)!=_atoms.end());
        molecule_value += fraction[n]*_atoms[atom_name].atom_value;
      }

      for(unsigned int n=0; n<atoms.size(); ++n)
      {
        const std::string & atom_name = atoms[n];
        _out << "  <fraction n=\""<< fraction[n]*_atoms[atom_name].atom_value/molecule_value <<"\" ref=\"" << atom_name <<"\" />" <<std::endl;
      }

      _out << "  </material> " <<std::endl;
    }

    // write material "air", for volume "Word"
    _out << "  <material formula=\" \" name=\"air\" >" <<std::endl;
    _out << "  <D value=\"1.290\" unit=\"mg/cm3\"/>" <<std::endl;
    _out << "  <fraction n=\"0.7\" ref=\"Nitrogen\" />" <<std::endl;
    _out << "  <fraction n=\"0.3\" ref=\"Oxygen\" />" <<std::endl;
    _out << "  </material>" <<std::endl;

    //end material section
    _out << "</materials>" << std::endl << std::endl;


    //begin solids section
    _out << "<solids>" << std::endl;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);

      _out << "  <tessellated name=\"Region_" << region->name() << "\">" << std::endl;

      for(unsigned int n=0; n<sl.size(); n++)
      {
        const Elem * elem = mesh.elem(el[n]);
        if( elem->subdomain_id() != r) continue;

        AutoPtr<Elem> boundary_side = elem->build_side(sl[n]);

        /*
        Point p1 = *boundary_side->get_node(0);
        Point p2 = *boundary_side->get_node(1);
        Point p3 = *boundary_side->get_node(2);
        Point cross = ((p2 - p1).cross(p3 - p2)).unit();
        Point norm = elem->outside_unit_normal(sl[n]);
        genius_assert( norm*cross > 0 );
        */

        switch(boundary_side->type())
        {
            case TRI3      :
            case TRI3_FVM  :
            _out << "  <triangular   "
            << "vertex1=\"v" << boundary_side->get_node(0)->id()<<"\" "
            << "vertex2=\"v" << boundary_side->get_node(1)->id()<<"\" "
            << "vertex3=\"v" << boundary_side->get_node(2)->id()<<"\" "
            << "type=\"ABSOLUTE\"/>" << std::endl;
            break;
            case QUAD4 :
            case QUAD4_FVM :
            _out << "  <quadrangular "
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

    //write solid world, 2*bbox
    _out << "  <box name=\"world\" "
    << "x=\"" << 2*x_len <<"\" "
    << "y=\"" << 2*y_len <<"\" "
    << "z=\"" << 2*z_len <<"\" "
    << "lunit=\"micron\" />" <<std::endl;

    //end solids section
    _out << "</solids>" << std::endl << std::endl;


    //begin structure section
    _out << "<structure>" << std::endl;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);

      _out << "  <volume name=\"Volume_" << region->name() << "\">" << std::endl;
      _out << "  <materialref ref=\"" << region->name() + "_" + region->material() <<"\"/>" << std::endl;
      _out << "  <solidref ref=\"Region_" << region->name() <<"\"/>" << std::endl;
      _out << "  </volume>" << std::endl;
    }

    _out << "  <volume name=\"World\" >" << std::endl;
    _out << "  <materialref ref=\"air\"/>" << std::endl;
    _out << "  <solidref ref=\"world\"/>" << std::endl;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      _out << "    <physvol>" << std::endl;
      _out << "       <volumeref ref=\"Volume_" << region->name() << "\"/>" << std::endl;
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



//FIXME write all the mesh element into GDML seems over killed too much.
// Geant4 may spend hours for rebuild geometry object
void GDMLIO::write_gdml_body (const std::string& filename)
{
#if 0
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  const MeshBase & mesh = system.mesh();

  genius_assert(mesh.mesh_dimension()==3);
  set_atoms();

  // collent doping information
  const double concentration_scale = pow(cm, -3);
  Na_max = Nd_max = 1.0/concentration_scale;

  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);
    SimulationRegion::const_processor_node_iterator on_processor_nodes_it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator on_processor_nodes_it_end = region->on_processor_nodes_end();
    for(; on_processor_nodes_it!=on_processor_nodes_it_end; ++on_processor_nodes_it)
    {
      const FVM_Node * fvm_node = *on_processor_nodes_it;
      const FVM_NodeData * node_data = fvm_node->node_data();
      const unsigned int id = fvm_node->root_node()->id();

      // if the fvm_node lies on the interface of two material regions,
      // we shall use the node data in the more important region.
      // it just for visualization reason
      if( fvm_node->boundary_id() != BoundaryInfo::invalid_id )
      {
        unsigned int bc_index = system.get_bcs()->get_bc_index_by_bd_id(fvm_node->boundary_id());
        const BoundaryCondition * bc = system.get_bcs()->get_bc(bc_index);
        const FVM_Node * primary_fvm_node = (*bc->region_node_begin(fvm_node->root_node())).second.second;
        node_data = primary_fvm_node->node_data();
      }

      Na[id]  = node_data->Total_Na()/concentration_scale;
      Nd[id]  = node_data->Total_Nd()/concentration_scale;
      Na_max  = std::max(node_data->Total_Na()/concentration_scale, Na_max);
      Nd_max  = std::max(node_data->Total_Nd()/concentration_scale, Nd_max);
    }
  }
  Parallel::gather(0, Na);
  Parallel::gather(0, Nd);

  Parallel::max(Na_max);
  Parallel::max(Nd_max);


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

    // begin define section
    _out << "<define>" << std::endl;
    _out << std::setprecision(8) << std::scientific << std::right;

    // for metal/insulator region, write the whole regions as one element
    // however, semiconductor region should have more info

    //pick boundary nodes, write into define/position
    std::vector<unsigned int>       el;
    std::vector<unsigned short int> sl;
    std::vector<short int>          il;
    system.mesh().boundary_info->build_side_list (el, sl, il);

    std::set<const Node *> _nodes;
    std::vector<const Elem *> semiconductor_cells;

    for(unsigned int n=0; n<sl.size(); n++)
    {
      const Elem * elem = mesh.elem(el[n]);
      AutoPtr<Elem> boundary_side = elem->build_side(sl[n]);
      for(unsigned int nn=0; nn<boundary_side->n_nodes(); ++nn)
        _nodes.insert(boundary_side->get_node(nn));
    }

    MeshBase::const_element_iterator       elem_it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator elem_end = mesh.active_elements_end();
    for (; elem_it != elem_end; ++elem_it)
    {
      const Elem *elem  = (*elem_it);
      const SimulationRegion * region = system.region(elem->subdomain_id());
      if(region->type()!=SemiconductorRegion) continue;
      semiconductor_cells.push_back(elem);

      for(unsigned int n=0; n<elem->n_nodes(); ++n)
      {
        _nodes.insert(elem->get_node(n));
      }
    }


    for(std::set<const Node *>::iterator it=_nodes.begin(); it!=_nodes.end(); ++it)
      {
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
    std::set<std::string> atom_refs;
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      std::vector<std::string> atoms;
      std::vector<double> fraction;
      region->atom_fraction(atoms, fraction);

      for(unsigned int n=0; n<atoms.size(); ++n)
      {
        genius_assert(_atoms.find(atoms[n])!=_atoms.end());
        atom_refs.insert(atoms[n]);
      }
    }
    // atoms for doping
    atom_refs.insert("Boron");
    atom_refs.insert("Phosphorus");

    // two atoms for Air
    atom_refs.insert("Nitrogen");
    atom_refs.insert("Oxygen");

    // write all the atoms
    for(std::set<std::string>::iterator it=atom_refs.begin(); it!=atom_refs.end(); ++it)
      {
        std::map<std::string, Atom>::iterator atom_it = _atoms.find(*it);
        write_gdml_element(atom_it->first);
      }

    // write the materials.  special treatment to semiconductor region
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      double density = region->get_density(system.T_external())/(g*pow(cm,-3));

      if(region->type()!=SemiconductorRegion)
      {
        std::vector<std::string> atoms;
        std::vector<double> fraction;
        region->atom_fraction(atoms, fraction);

        double molecule_value= 0;
        for(unsigned int n=0; n<atoms.size(); ++n)
        {
          const std::string & atom_name = atoms[n];
          genius_assert(_atoms.find(atom_name)!=_atoms.end());
          molecule_value += fraction[n]*_atoms[atom_name].atom_value;
        }

        _out << "  <material name=\"" << region->name() + "_" + region->material() <<"\" formula=\"" << region->material() <<"\">" <<std::endl;
        _out << "  <D value=\"" << region->get_density(system.T_external())/(g*pow(cm,-3)) <<"\" />" <<std::endl;
        for(unsigned int n=0; n<atoms.size(); ++n)
        {
          const std::string & atom_name = atoms[n];
          _out << "  <fraction n=\""<< fraction[n]*_atoms[atom_name].atom_value/molecule_value <<"\" ref=\"" << atom_name <<"\" />" <<std::endl;
        }

        _out << "  </material> " <<std::endl;
      }
      else
      {
        // process doping information
        for(int na=0; na<10; ++na)
          for(int nd=0; nd<10; ++nd)
          {
            std::vector<std::string> atoms;
            std::vector<double> fraction;
            region->atom_fraction(atoms, fraction);

            double molecule_value= 0;
            for(unsigned int n=0; n<atoms.size(); ++n)
            {
              const std::string & atom_name = atoms[n];
              genius_assert(_atoms.find(atom_name)!=_atoms.end());
              molecule_value += fraction[n]*_atoms[atom_name].atom_value;
            }

            double doping_Na = Na_max*pow(0.5, na);
            double doping_Nd = Nd_max*pow(0.5, nd);
            double d = density/molecule_value*6.022142e23; //Avogadro constant
            // insert Na and Nd into atom vector
            atoms.push_back("Boron");
            fraction.push_back(doping_Na/d);
            atoms.push_back("Phosphorus");
            fraction.push_back(doping_Nd/d);

            // current fraction
            double sum_fraction = std::accumulate(fraction.begin(), fraction.end(), 0.0);
            for(unsigned int n=0; n<fraction.size(); ++n)
              fraction[n] /= sum_fraction;
            molecule_value= 0;
            for(unsigned int n=0; n<atoms.size(); ++n)
            {
              const std::string & atom_name = atoms[n];
              genius_assert(_atoms.find(atom_name)!=_atoms.end());
              molecule_value += fraction[n]*_atoms[atom_name].atom_value;
            }

            _out << "  <material name=\"" << region->name() + "_" + region->material()<< "_" << na << "_" << nd <<"\" formula=\"" << region->material() <<"\">" <<std::endl;
            _out << "  <D value=\"" << region->get_density(system.T_external())/(g*pow(cm,-3)) <<"\" />" <<std::endl;
            for(unsigned int n=0; n<atoms.size(); ++n)
            {
              const std::string & atom_name = atoms[n];
              _out << "  <fraction n=\""<< fraction[n]*_atoms[atom_name].atom_value/molecule_value <<"\" ref=\"" << atom_name <<"\" />" <<std::endl;
            }

            _out << "  </material> " <<std::endl;

          }
      }
    }

    // write material "air", for volume "Word"
    _out << "  <material formula=\" \" name=\"air\" >" <<std::endl;
    _out << "  <D value=\"1.290\" unit=\"mg/cm3\"/>" <<std::endl;
    _out << "  <fraction n=\"0.7\" ref=\"Nitrogen\" />" <<std::endl;
    _out << "  <fraction n=\"0.3\" ref=\"Oxygen\" />" <<std::endl;
    _out << "  </material>" <<std::endl;

    //end material section
    _out << "</materials>" << std::endl << std::endl;




    //begin solids section
    _out << "<solids>" << std::endl;

    // write surface except semiconductor region
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      if(region->type()==SemiconductorRegion) continue;

      _out << "  <tessellated name=\"Region_" << region->name() << "\">" << std::endl;

      for(unsigned int n=0; n<sl.size(); n++)
      {
        const Elem * elem = mesh.elem(el[n]);
        if( elem->subdomain_id() != r) continue;
        AutoPtr<Elem> boundary_side = elem->build_side(sl[n]);

        switch(boundary_side->type())
        {
            case TRI3      :
            case TRI3_FVM  :
            _out << "  <triangular   "
            << "vertex1=\"v" << boundary_side->get_node(0)->id()<<"\" "
            << "vertex2=\"v" << boundary_side->get_node(1)->id()<<"\" "
            << "vertex3=\"v" << boundary_side->get_node(2)->id()<<"\" "
            << "type=\"ABSOLUTE\"/>" << std::endl;
            break;
            case QUAD4 :
            case QUAD4_FVM :
            _out << "  <quadrangular "
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

    // write body for semiconductor region
    for (unsigned int e=0; e<semiconductor_cells.size(); ++e)
    {
      const Elem *elem  = semiconductor_cells[e];
      const SimulationRegion * region = system.region(elem->subdomain_id());

      _out << "  <tessellated name=\"Region_"<< region->name()<< "_" << elem->id() << "\">" << std::endl;

      for(unsigned int n=0; n<elem->n_sides(); n++)
      {
        AutoPtr<Elem> boundary_side = elem->build_side(n);

        switch(boundary_side->type())
        {
            case TRI3  :
            _out << "    <triangular   "
            << "vertex1=\"v" << boundary_side->get_node(0)->id()<<"\" "
            << "vertex2=\"v" << boundary_side->get_node(1)->id()<<"\" "
            << "vertex3=\"v" << boundary_side->get_node(2)->id()<<"\" "
            << "type=\"ABSOLUTE\"/>" << std::endl;
            break;
            case QUAD4 :
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

    //write solid world, 2*bbox
    _out << "  <box name=\"world\" "
    << "x=\"" << 2*x_len <<"\" "
    << "y=\"" << 2*y_len <<"\" "
    << "z=\"" << 2*z_len <<"\" "
    << "lunit=\"micron\" />" <<std::endl;

    //end solids section
    _out << "</solids>" << std::endl << std::endl;




    //begin structure section
    _out << "<structure>" << std::endl;

    // write structure except semiconductor region
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      if(region->type()==SemiconductorRegion) continue;
      _out << "  <volume name=\"Volume_" << region->name() << "\">" << std::endl;
      _out << "  <materialref ref=\"" << region->name() + "_" + region->material() <<"\"/>" << std::endl;
      _out << "  <solidref ref=\"Region_" << region->name() <<"\"/>" << std::endl;
      _out << "  </volume>" << std::endl;
    }

    // write structure for semiconductor region, consider doping information
    for (unsigned int e=0; e<semiconductor_cells.size(); ++e)
    {
      const Elem *elem  = semiconductor_cells[e];
      const SimulationRegion * region = system.region(elem->subdomain_id());
      if(region->type()!=SemiconductorRegion) continue;

      double doping_Na, doping_Nd;
      cell_average_doping(elem, doping_Na, doping_Nd);
      doping_Na = std::max(doping_Na, 1.0*concentration_scale);
      doping_Nd = std::max(doping_Nd, 1.0*concentration_scale);
      int na = static_cast<int>(log(doping_Na/Na_max)/log(0.5));
      int nd = static_cast<int>(log(doping_Nd/Nd_max)/log(0.5));
      na = std::max(na, 0);
      na = std::min(na, 9);
      nd = std::max(nd, 0);
      nd = std::min(nd, 9);

      _out << "  <volume name=\"Volume_" << region->name() + "_" << elem->id() << "\">" << std::endl;
      _out << "  <materialref ref=\"" << region->name() + "_" + region->material() << "_" << na << "_" << nd <<"\"/>" << std::endl;
      _out << "  <solidref ref=\"Region_" << region->name() + "_" << elem->id() <<"\"/>" << std::endl;
      _out << "  </volume>" << std::endl;
    }

    _out << "  <volume name=\"World\" >" << std::endl;
    _out << "  <materialref ref=\"air\"/>" << std::endl;
    _out << "  <solidref ref=\"world\"/>" << std::endl;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      if(region->type()==SemiconductorRegion) continue;
      _out << "    <physvol>" << std::endl;
      _out << "       <volumeref ref=\"Volume_" << region->name() << "\"/>" << std::endl;
      _out << "    </physvol>" << std::endl;
    }

    for (unsigned int e=0; e<semiconductor_cells.size(); ++e)
    {
      const Elem *elem  = semiconductor_cells[e];
      const SimulationRegion * region = system.region(elem->subdomain_id());
      if(region->type()!=SemiconductorRegion) continue;
      _out << "    <physvol>" << std::endl;
      _out << "       <volumeref ref=\"Volume_" << region->name() + "_" << elem->id() << "\"/>" << std::endl;
      _out << "    </physvol>" << std::endl;
    }
    _out << "  </volume>" << std::endl;


    //end structure section
    _out << "</structure>" << std::endl << std::endl;




    //write setup section
    _out << "<setup name=\"GeniusToGeant4\" version=\"1.0\">" << std::endl;
    _out << "  <world ref=\"World\"/>" << std::endl;
    _out << "</setup>" << std::endl << std::endl;

    _out << "</gdml>" << std::endl << std::endl;


    //record minimal elem size
    double minimal_size=1e30;
    for ( elem_it = mesh.active_elements_begin(); elem_it != elem_end; ++elem_it)
      minimal_size = std::min(minimal_size, (*elem_it)->hmin());

    _out << "<!-- minimal volume size is " << minimal_size/um << " um -->" << std::endl;

    _out.close();
  }
#endif
}


void GDMLIO::cell_average_doping(const Elem *elem, double &doping_Na, double &doping_Nd)
{
  double _Na = 0;
  double _Nd = 0;
  for(unsigned int n=0; n<elem->n_nodes(); ++n)
  {
    unsigned int node_id = elem->get_node(n)->id();
    _Na += Na.find(node_id)->second;
    _Nd += Nd.find(node_id)->second;
  }

  doping_Na = _Na/elem->n_nodes();
  doping_Nd = _Nd/elem->n_nodes();
}


void GDMLIO::write_gdml_element(const std::string &name)
{
  std::map<std::string, Atom>::const_iterator atom_it = _atoms.find(name);

  const std::string formula = atom_it->second.formula;
  int Z = atom_it->second.Z;
  double atom_value = atom_it->second.atom_value;

  // if this element has isotope
  if(atom_it->second.isotope_info.size())
  {
    for(unsigned int n=0; n<atom_it->second.isotope_info.size(); ++n)
    {
      const Atom::Isotope & isotope =  atom_it->second.isotope_info[n];
      _out << "  <isotope name=\"" << isotope.isotope_name <<"\" Z=\"" << Z <<"\" N=\"" << isotope.nucleons <<"\">" << std::endl;
      _out << "    <atom type=\"A\" value=\"" << isotope.atom_value << "\"/>" << std::endl;
      _out << "  </isotope>" <<std::endl;
    }

    _out << "  <element name=\"" << name <<"\" >" <<std::endl;
    for(unsigned int n=0; n<atom_it->second.isotope_info.size(); ++n)
    {
      const Atom::Isotope & isotope =  atom_it->second.isotope_info[n];
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



void GDMLIO::set_atoms()
{
  _atoms.clear();

  _atoms["Hydrogen"  ] = Atom("Hydrogen", "H", 1, 1.00797);

  _atoms["Boron"     ] = Atom("Boron",    "B", 5, 10.811);
  _atoms["Boron"     ].add_isotope("B10", 10, 10.0, 0.199);
  _atoms["Boron"     ].add_isotope("B11", 11, 11.0, 0.801);

  _atoms["Carbon"    ] = Atom("Carbon",   "C", 6, 12.01115);
  _atoms["Nitrogen"  ] = Atom("Nitrogen", "N", 7, 14.01);
  _atoms["Oxygen"    ] = Atom("Oxygen",   "O", 8, 16.0);

  _atoms["Aluminum"  ] = Atom("Aluminum",   "Al", 13, 26.98);
  _atoms["Silicon"   ] = Atom("Silicon",    "Si", 14, 28.086);
  _atoms["Phosphorus"] = Atom("Phosphorus", "P",  15, 30.9738);

  _atoms["Titanium"  ] = Atom("Titanium",  "Tl", 22, 47.9);
  _atoms["Copper"    ] = Atom("Copper",    "Cu", 29, 63.54);
  _atoms["Gallium"   ] = Atom("Gallium",   "Ga", 31, 69.72);
  _atoms["Germanium" ] = Atom("Germanium", "Ge", 32, 72.59);
  _atoms["Arsenic"   ] = Atom("Arsenic",   "As", 33, 74.922);
  _atoms["Silver"    ] = Atom("Silver",    "Ag", 47, 107.87);
  _atoms["Cadmium"   ] = Atom("Cadmium",   "Cd", 48, 112.40);
  _atoms["Indium"    ] = Atom("Indium",    "In", 49, 114.82);
  _atoms["Antimony"  ] = Atom("Antimony",  "Sb", 51, 121.75);
  _atoms["Tellurium" ] = Atom("Tellurium", "Te", 52, 127.6);

  _atoms["Hafnium"   ] = Atom("Hafnium",   "Hf", 72, 178.49);
  _atoms["Gold"      ] = Atom("Gold",      "Au", 79, 196.97);
  _atoms["Mercury"   ] = Atom("Mercury",   "Hg", 80, 200.59);

}


