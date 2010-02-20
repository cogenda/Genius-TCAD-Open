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

//  $Id: control.cc,v 1.54 2008/07/09 12:56:23 gdiso Exp $

#include "genius_common.h"

#ifdef CYGWIN
    #include <io.h>
#else
    #include <unistd.h>
#endif


#include "parser.h"
#include "control.h"

#include "mesh_generation_tri3.h"
#include "mesh_generation_quad4.h"

// only commercial product support 3D mesh generator
#ifdef COGENDA_COMMERCIAL_PRODUCT
  #include "mesh_generation_hex8.h"
  #include "mesh_generation_prism6.h"
  #include "mesh_generation_tet4.h"
#endif

#include "mesh_tools.h"
#include "mesh_communication.h"
#include "mesh_refinement.h"
#include "mesh_modification.h"
#include "boundary_info.h"
#include "electrical_source.h"
#include "field_source.h"
#include "enum_solution.h"
#include "doping_analytic/doping_analytic.h"
#include "mole_analytic/mole_analytic.h"
#include "poisson/poisson.h"
#include "ddm1/ddm1.h"
#include "ddm2/ddm2.h"
#include "ebm3/ebm3.h"
#include "ddm_ac/ddm_ac.h"
#include "emfem2d/emfem2d.h"
#include "ray_tracing/ray_tracing.h"

#ifndef CYGWIN
 #include "mix1/mix1.h"
 #include "mix2/mix2.h"
 #include "mix3/mix3.h"
#endif

#include "spice_ckt.h"
#include "mixA1/mixA1.h"
#include "mixA2/mixA2.h"
#include "mixA3/mixA3.h"

#include "hall/hall.h"

#include "stress_solver/stress_solver.h"


#include "solver_specify.h"
#include "advanced_model.h"
#include "petsc_type.h"

#include "interpolation_2d_csa.h"
#include "interpolation_3d_qshep.h"
#include "interpolation_3d_nbtet.h"

#include "dlhook.h"
#ifdef CYGWIN
 #include "cv_hook.h"
 #include "rawfile_hook.h"
 #include "gnuplot_hook.h"
 #include "probe_hook.h"
 #include "vtk_hook.h"
#endif

#ifdef ENABLE_VISIT
 #include "visit_hook.h"
#endif

#include "extend_to_3d.h"
#ifdef HAVE_X11
 #include "show_mesh_2d.h"
#endif

#include "MXMLUtil.h"

using PhysicalUnit::A;
using PhysicalUnit::V;
using PhysicalUnit::W;
using PhysicalUnit::C;
using PhysicalUnit::s;


//------------------------------------------------------------------------------
SolverControl::SolverControl()
    : _decks(NULL), _mesh(NULL), _system(NULL)
{
  _dom_solution = mxmlNewXML("1.0");
  mxmlNewElement(_dom_solution, "genius-solutions");
}

SolverControl::~SolverControl()
{
  mxmlDelete(_dom_solution);
}

mxml_node_t* SolverControl::get_dom_solution() const
{
  return _dom_solution;
}

int SolverControl::get_dom_solution_count() const
{
  int cnt=0;
  mxml_node_t* root = mxmlFindElement(_dom_solution, _dom_solution, "genius-solutions", NULL, NULL, MXML_DESCEND_FIRST);
  if (!root)
    return 0;

  for(mxml_node_t* node = mxmlFindElement(root, root, "solution", NULL, NULL, MXML_DESCEND_FIRST);
      node; node = mxmlFindElement(node, root, "solution", NULL, NULL, MXML_NO_DESCEND) )
  {
    cnt++;
  }
  return cnt;
}



void SolverControl::setDecks(Parser::InputParser *input)
{
  _decks = input;
}

int SolverControl::reset_simulation_system()
{
  if (_decks == NULL)
    return 0;

  _mesh = AutoPtr<Mesh>(new Mesh(3));
  _system = AutoPtr<SimulationSystem>(new SimulationSystem(mesh(), decks()));
  return 0;
}

//------------------------------------------------------------------------------
int SolverControl::mainloop()
{
  if (_decks == NULL)
    return 0;

  if (_mesh.get() == NULL || _system.get() == NULL)
    reset_simulation_system();

  // first, we should see if mesh generation card exist
  this->do_mesh();

  // then, we should see if doping profile and/or mole card exist
  this->do_process();

  // from above tow steps, maybe the simulation system has been build.
  // if not, user should use IMPORT command to get an (previous) system into memory.

  // we can begin the main loop here
  for( decks().begin(); !decks().end(); decks().next() )
  {

    Parser::Card c = decks().get_current_card();

    if(c.key() == "MODEL")
      this->set_model ( c );

    if(c.key() == "METHOD")
      this->set_method ( c );

    if(c.key() == "HOOK")
      this->do_hook( c );

    if(c.key() == "SOLVE")
      this->do_solve( c );

    if(c.key() == "EXPORT")
      this->do_export( c );

    if(c.key() == "IMPORT")
      this->do_import( c );

    if(c.key() == "EMFEM2D")
      this->do_em_fem2d_solve( c );

    if(c.key() == "RAYTRACE")
      this->do_ray_trace( c );

    if(c.key() == "NODESET")
      this->set_initial_node_voltage( c );

    if(c.key() == "REFINE.CONFORM")
      this->do_refine_conform( c );

    if(c.key() == "REFINE.HIERARCHICAL")
      this->do_refine_hierarchical( c );

    if(c.key() == "REFINE.UNIFORM")
      this->do_refine_uniform( c );

    if(c.key() == "PMI")
      this->set_physical_model ( c );

    if(c.key() == "ATTACH")
      this->set_electrode_source ( c );

    if(c.key() == "EXTEND")
      this->extend_to_3d( c );

    if(c.key() == "PLOTMESH")
      this->plot_mesh( c );
  }

  return 0;
}


//------------------------------------------------------------------------------
int  SolverControl::do_mesh()
{

  if ( decks().is_card_exist("MESH") )
  {
    // build meshgenerator only on processor 0
    // I am afraid about mesh generator may have different
    // behavior due to float point round-off error
    if (Genius::processor_id() == 0)
    {
      // which mesh generator should we use?
      for( decks().begin(); !decks().end(); decks().next() )
      {
        Parser::Card c = decks().get_current_card();

        if(c.key() == "MESH")
        {
          if( c.is_enum_value("type","s_tri3"))
          {
            meshgen = AutoPtr<MeshGeneratorBase>(new MeshGeneratorTri3(mesh(),decks()));
          }
          else if( c.is_enum_value("type","s_quad4"))
          {
            meshgen = AutoPtr<MeshGeneratorBase>(new MeshGeneratorQuad4(mesh(),decks()));
          }
#ifdef COGENDA_COMMERCIAL_PRODUCT
          else if( c.is_enum_value("type","s_tet4"))
          {
            meshgen = AutoPtr<MeshGeneratorBase>(new MeshGeneratorTet4(mesh(),decks()));
          }
          else if( c.is_enum_value("type","s_prism6"))
          {
            meshgen = AutoPtr<MeshGeneratorBase>(new MeshGeneratorPrism6(mesh(),decks()));
          }
          else if( c.is_enum_value("type","s_hex8"))
          {
            meshgen = AutoPtr<MeshGeneratorBase>(new MeshGeneratorHex8(mesh(),decks()));
          }
#else
          else
          {
            MESSAGE<<"ERROR: 3D mesh generator is not supported by Open Source Version." << std::endl; RECORD();
            genius_error();
          }
#endif
        }
      }
      // ok, generate the mesh
      if ( meshgen->do_mesh() )
      {
        MESSAGE<<"ERROR: Mesh generation failed." << std::endl; RECORD();
        genius_error();
      }
    }

    // since we only build mesh on processor 0,
    // sync mesh to other processors.
    // this procedure also prepare the mesh for using
    MeshCommunication mesh_comm;
    mesh_comm.broadcast(mesh());

    // please note, until here, mesh is still not prepared
    // mesh.is_prepared() will return false

    // now we can build solution system since mesh is done
    system().build_simulation_system();
    system().sync_print_info();
  }

  return 0;
}



int SolverControl::do_process()
{
  // if doping profile card exist
  if ( decks().is_card_exist("DOPING") )
  {
    DopingSolver = AutoPtr<SolverBase>( new DopingAnalytic(system(), decks()) );
    // parse "PROFILE" card
    DopingSolver->create_solver();
    // set doping profile to semiconductor region
    DopingSolver->solve();
    // we will not destroy the doping solver here
  }

  // if mole card exist
  if ( decks().is_card_exist("MOLE") )
  {
    MoleSolver = AutoPtr<SolverBase>( new MoleAnalytic(system(), decks()) );
    // parse "MOLE" card
    MoleSolver->create_solver();
    // set mole fraction to semiconductor region
    MoleSolver->solve();
    // we will not destroy the mole solver here
  }

  // after doping profile and mole is set, we can init system data.
  // i.e. electron and hole initital concentration of semiconductor region

  // note: even no mesh and/or process are done, we can still call it safely
  // although it will do nothing
  system().init_region();

  return 0;
}




//------------------------------------------------------------------------------
int SolverControl::set_method ( const Parser::Card & c )
{
  // reset to default solver parameters
  SolverSpecify::set_default_parameter();

  // set nonlinear solver type
  SolverSpecify::NS = SolverSpecify::nonlinear_solver_type(c.get_string("ns", "linesearch"));

  // set linear solver type
  SolverSpecify::LS = SolverSpecify::linear_solver_type(c.get_string("ls", "bcgs"));

  // set preconditioner type
  SolverSpecify::PC = SolverSpecify::preconditioner_type(c.get_string("pc", "asm"));

  // set Newton damping type
  if(c.is_parameter_exist("damping"))
  {
    if (c.is_enum_value("damping", "no"))        SolverSpecify::Damping = SolverSpecify::DampingNo;
    if (c.is_enum_value("damping", "potential")) SolverSpecify::Damping = SolverSpecify::DampingPotential;
    if (c.is_enum_value("damping", "bankrose"))  SolverSpecify::Damping = SolverSpecify::DampingBankRose;
  }

  //set convergence test
  SolverSpecify::MaxIteration              = c.get_int("maxiteration", 30);

  SolverSpecify::ksp_rtol                  = c.get_real("ksp.rtol", 1e-8);
  SolverSpecify::ksp_atol                  = c.get_real("ksp.atol", 1e-20);
  SolverSpecify::ksp_atol_fnorm            = c.get_real("ksp.atol.fnorm", 1e-7);

  SolverSpecify::relative_toler            = c.get_real("relative.tol", 1e-5);
  SolverSpecify::toler_relax               = c.get_real("toler.relax", 1e4);
  SolverSpecify::poisson_abs_toler         = c.get_real("poisson.tol", 1e-29)*C;
  SolverSpecify::elec_continuity_abs_toler = c.get_real("elec.continuity.tol", 5e-18)*A;
  SolverSpecify::hole_continuity_abs_toler = c.get_real("hole.continuity.tol", 5e-18)*A;
  SolverSpecify::heat_equation_abs_toler   = c.get_real("latt.temp.tol", 1e-11)*W;
  SolverSpecify::elec_energy_abs_toler     = c.get_real("elec.energy.tol", 1e-18)*W;
  SolverSpecify::hole_energy_abs_toler     = c.get_real("hole.energy.tol", 1e-18)*W;
  SolverSpecify::electrode_abs_toler       = c.get_real("electrode.tol", 1e-9)*V;
  SolverSpecify::elec_quantum_abs_toler    = c.get_real("elec.quantum.tol", 1e-29)*C;
  SolverSpecify::hole_quantum_abs_toler    = c.get_real("hole.quantum.tol", 1e-29)*C;

  // set which solver will be used
  if(c.is_parameter_exist("type"))
  {
    if (c.is_enum_value("type", "poisson"))                 SolverSpecify::Solver = SolverSpecify::POISSON;
    else if (c.is_enum_value("type", "ddml1"))              SolverSpecify::Solver = SolverSpecify::DDML1;
    else if (c.is_enum_value("type", "ddml1s"))             SolverSpecify::Solver = SolverSpecify::DDML1MIX;
    else if (c.is_enum_value("type", "ddml1m"))             SolverSpecify::Solver = SolverSpecify::DDML1MIXA;
    else if (c.is_enum_value("type", "hall"))               SolverSpecify::Solver = SolverSpecify::HALLDDML1;
    else if (c.is_enum_value("type", "ddml2"))              SolverSpecify::Solver = SolverSpecify::DDML2;
    else if (c.is_enum_value("type", "ddml2s"))             SolverSpecify::Solver = SolverSpecify::DDML2MIX;
    else if (c.is_enum_value("type", "ddml2m"))             SolverSpecify::Solver = SolverSpecify::DDML2MIXA;
    else if (c.is_enum_value("type", "ebml3"))              SolverSpecify::Solver = SolverSpecify::EBML3;
    else if (c.is_enum_value("type", "ebml3s"))             SolverSpecify::Solver = SolverSpecify::EBML3MIX;
    else if (c.is_enum_value("type", "ebml3m"))             SolverSpecify::Solver = SolverSpecify::EBML3MIXA;
    else if (c.is_enum_value("type", "ddmac"))              SolverSpecify::Solver = SolverSpecify::DDMAC;
  }

  SolverSpecify::ServerPort = static_cast<unsigned short int>(c.get_int("serverport", 1611));

  return 0;
}



/*--------------------------------------------------------
 *  set advanced semiconductor property
 */
int SolverControl::set_model ( const Parser::Card & c )
{
  std::string region_label = c.get_string("region", "");

  AdvancedModel model;

  // advanced mobility model control
  model.ESurface                          = c.get_bool("esurface", false);
  model.HighFieldMobility                 = c.get_bool("highfieldmobility", true, "h.mob");
  model.HighFieldMobilityAD               = c.get_bool("highfieldmobilityad", true, "h.mob.ad");
  model.HighFieldMobilitySelfConsistently = c.get_bool("h.mob.selfconsistent", true);
  if(c.is_parameter_exist("mob.force") || c.is_parameter_exist("mobility.force"))
  {
    if (c.is_enum_value("mob.force", "ej") || c.is_enum_value("mobility.force", "ej"))                 model.Mob_Force = ModelSpecify::EJ;
    if (c.is_enum_value("mob.force", "esimple") || c.is_enum_value("mobility.force", "esimple"))       model.Mob_Force = ModelSpecify::ESimple;
    if (c.is_enum_value("mob.force", "eqf") || c.is_enum_value("mobility.force", "eqf"))               model.Mob_Force = ModelSpecify::EQF;
  }

  //impact ionization model
  model.ImpactIonization = false;
  if(c.is_parameter_exist("impactionization") || c.is_parameter_exist("ii"))
  {
    if(c.is_enum_value("impactionization", "local") || c.is_enum_value("ii", "local"))
    {
      model.ImpactIonization = true;
      if(c.is_parameter_exist("ii.force"))
      {
        model.II_Force = ModelSpecify::GradQf; // default
        if (c.is_enum_value("ii.force", "edotj"))  model.II_Force = ModelSpecify::IIForce_EdotJ;
        if (c.is_enum_value("ii.force", "eside"))  model.II_Force = ModelSpecify::ESide;
        if (c.is_enum_value("ii.force", "evector"))  model.II_Force = ModelSpecify::EVector;
        if (c.is_enum_value("ii.force", "gradqf"))  model.II_Force = ModelSpecify::GradQf;
      }
    }
  }

  model.BandBandTunneling = false;
  if(c.is_parameter_exist("bandbandtunneling") || c.is_parameter_exist("bbt"))
  {
    if(c.is_enum_value("bandbandtunneling", "local") || c.is_enum_value("bbt", "local"))
    {
      model.BandBandTunneling = true;
    }
  }

  // fermi statistics and incomplete ionization
  model.Fermi                 = c.get_bool("fermi", false);
  model.IncompleteIonization  = c.get_bool("incompleteionization", false);

  // charge trapping model
  model.Trap                  = c.get_bool("trap",false);

  //energy balance advanced model
  if(c.is_parameter_exist("eb.level"))
  {
    if (c.is_enum_value("eb.level", "none")) { model.EB_Level = ModelSpecify::NONE;   }
    if (c.is_enum_value("eb.level", "te"))   { model.EB_Level = ModelSpecify::Tn;     }
    if (c.is_enum_value("eb.level", "th"))   { model.EB_Level = ModelSpecify::Tp;     }
    if (c.is_enum_value("eb.level", "tl"))   { model.EB_Level = ModelSpecify::Tl;     }
    if (c.is_enum_value("eb.level", "teth")) { model.EB_Level = ModelSpecify::TnTp;   }
    if (c.is_enum_value("eb.level", "tetl")) { model.EB_Level = ModelSpecify::TnTl;   }
    if (c.is_enum_value("eb.level", "thtl")) { model.EB_Level = ModelSpecify::TpTl;   }
    if (c.is_enum_value("eb.level", "all"))  { model.EB_Level = ModelSpecify::ALL;    }
  }

  if( system().region(region_label)==NULL )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " MODEL: Region " << region_label << " can't be found in mesh regions." << std::endl; RECORD();
    genius_error();
  }


  system().region(region_label)->advanced_model() = model;

  // we should check lattice temperature model, which should be set for all the regions!
  {
    bool temperature = false;
    for(unsigned int n=0; n < system().n_regions(); ++n)
    {
      SimulationRegion * region = system().region(n);
      if (region->get_advanced_model()->enable_Tl()) temperature = true;
    }

    if(temperature == true)
      for(unsigned int n=0; n < system().n_regions(); ++n)
      {
        SimulationRegion * region = system().region(n);
        region->advanced_model().force_temperature_usage();
      }

  }

  return 0;
}



int SolverControl::do_hook( const Parser::Card & c )
{
  std::string name;
  name = c.get_string("load", "");

  std::deque<std::string>::iterator hook_it = SolverSpecify::Hooks.begin();
  for (;hook_it!=SolverSpecify::Hooks.end();hook_it++)
  {
    if (*hook_it==name)
      break;
  }
  if (name.length()>0 && hook_it==SolverSpecify::Hooks.end())
  {
    SolverSpecify::Hooks.push_back(name);

    SolverSpecify::Hook_Parameters.erase(name);
    std::vector<Parser::Parameter> plist;
    for(unsigned int idx=0; idx<c.parameter_size(); idx++)
    {
      Parser::Parameter p = c.get_parameter(idx);
      if ( p.is_user_defined() )
        plist.push_back(p);
    }
    SolverSpecify::Hook_Parameters.insert(std::pair<std::string, std::vector<Parser::Parameter> >(name,plist));
  }

  name = c.get_string("unload", "");
  hook_it = SolverSpecify::Hooks.begin();
  for (;hook_it!=SolverSpecify::Hooks.end();hook_it++)
  {
    if (*hook_it==name)
      break;
  }
  if (name.length()>0 && hook_it!=SolverSpecify::Hooks.end())
  {
    SolverSpecify::Hooks.erase(hook_it);
    SolverSpecify::Hook_Parameters.erase(name);
  }

  return 0;
}


int SolverControl::do_solve( const Parser::Card & c )
{

  // set solution type solver will do
  SolverSpecify::Type = SolverSpecify::INVALID_SolutionType;
  if(c.is_parameter_exist("type"))
    SolverSpecify::Type = SolverSpecify::type_string_to_enum(c.get_string("type", ""));

  if(c.is_parameter_exist("label"))
    SolverSpecify::label = c.get_string("label", "");

  // set more detailed solution parameters
  switch (SolverSpecify::Type)
  {
  case SolverSpecify::EQUILIBRIUM : break;
  case SolverSpecify::STEADYSTATE :
  case SolverSpecify::OP :
    {
      SolverSpecify::OptG        = c.get_bool("optical.gen", false);
      SolverSpecify::PatG        = c.get_bool("particle.gen", false);
      SolverSpecify::RampUpSteps = c.get_int("rampupsteps", 1);
      SolverSpecify::RampUpVStep = c.get_real("rampupvstep", 0.25)*V;
      SolverSpecify::RampUpIStep = c.get_real("rampupistep", 0.1)*A;
      break;
    }
  case SolverSpecify::DCSWEEP     :
    {

      // user should specify vscan OR iscan
      genius_assert(c.is_parameter_exist("vscan") || c.is_parameter_exist("iscan"));
      genius_assert( !(c.is_parameter_exist("vscan") && c.is_parameter_exist("iscan")));

      // clear electrode vector
      SolverSpecify::Electrode_VScan.clear();
      SolverSpecify::Electrode_IScan.clear();

      if(c.is_parameter_exist("vscan"))
      {
        // sweep voltage of an electrode
        if(system().get_circuit()==NULL)
        {
          unsigned int elec_num = c.parameter_count("vscan");
          for(unsigned int n=0; n<elec_num; n++)
          {
            std::string electrode = c.get_n_string("vscan", "", n, 0);
            if( system().get_bcs()->get_bc(electrode) == NULL || system().get_bcs()->get_bc(electrode)->is_electrode() == false )
            {
              MESSAGE<<"ERROR at " <<c.get_fileline()<< " SOLVE: Electrode " << electrode << " can't be found in device structure." << std::endl; RECORD();
              genius_error();
            }
            SolverSpecify::Electrode_VScan.push_back(electrode);
          }

          if( !SolverSpecify::Electrode_VScan.size() )
          {
            MESSAGE<<"ERROR: You must specify at least one electrode for voltage DC scan."<<std::endl; RECORD();
            genius_error();
          }
        }
        // vscan spice voltage source
        else
        {
          std::string spice_vsource = c.get_string_lower_case("vscan", "");
          if( !system().get_circuit()->is_ckt_voltage_source_exist(spice_vsource))
          {
            MESSAGE<<"ERROR at " <<c.get_fileline()<< " SOLVE: VSRC " << spice_vsource << " can't be found in SPICE netlist." << std::endl; RECORD();
            genius_error();
          }
          SolverSpecify::Electrode_VScan.push_back(spice_vsource);
          if( SolverSpecify::Electrode_VScan.size() != 1 )
          {
            MESSAGE<<"ERROR: You must specify one VSRC for voltage DC scan."<<std::endl; RECORD();
            genius_error();
          }
        }
        SolverSpecify::VStart    = c.get_real("vstart", 0.0)*V;
        SolverSpecify::VStep     = c.get_real("vstep", 0.1)*V;
        SolverSpecify::VStepMax  = c.get_real("vstepmax", SolverSpecify::VStep/V)*V;
        SolverSpecify::VStop     = c.get_real("vstop", 5.0)*V;
      }

      if(c.is_parameter_exist("iscan"))
      {
        if(system().get_circuit()==NULL)
        {
          unsigned int elec_num = c.parameter_count("iscan");
          for(unsigned int n=0; n<elec_num; n++)
          {
            std::string electrode = c.get_n_string("iscan", "", n, 0);
            if( system().get_bcs()->get_bc(electrode) == NULL || system().get_bcs()->get_bc(electrode)->is_electrode() == false )
            {
              MESSAGE<<"ERROR at " <<c.get_fileline()<< " SOLVE: Electrode " << electrode << " can't be found in device structure." << std::endl; RECORD();
              genius_error();
            }
            SolverSpecify::Electrode_IScan.push_back(electrode);
          }
          if( !SolverSpecify::Electrode_IScan.size() )
          {
            MESSAGE<<"ERROR: You must specify at least one electrode for current DC scan."<<std::endl; RECORD();
            genius_error();
          }
        }
        // iscan spice current source
        else
        {
          std::string spice_isource = c.get_string_lower_case("iscan", "");
          if( !system().get_circuit()->is_ckt_current_source_exist(spice_isource))
          {
            MESSAGE<<"ERROR at " <<c.get_fileline()<< " SOLVE: ISRC " << spice_isource << " can't be found in SPICE netlist." << std::endl; RECORD();
            genius_error();
          }
          SolverSpecify::Electrode_IScan.push_back(spice_isource);
          if( SolverSpecify::Electrode_IScan.size() != 1 )
          {
            MESSAGE<<"ERROR: You must specify one ISRC for current DC scan."<<std::endl; RECORD();
            genius_error();
          }
        }
        SolverSpecify::IStart    = c.get_real("istart", 0.0)*A;
        SolverSpecify::IStep     = c.get_real("istep", 1e-5)*A;
        SolverSpecify::IStepMax  = c.get_real("istepmax", SolverSpecify::IStep/A)*A;
        SolverSpecify::IStop     = c.get_real("istop", 1e-2)*A;
      }

      SolverSpecify::Predict   = c.get_bool("predict", true);

      SolverSpecify::OptG      = c.get_bool("optical.gen", false);
      SolverSpecify::PatG      = c.get_bool("particle.gen", false);

      // set the waveform of light source
      std::string waveform = c.get_string("waveform", "");
      FieldSource * field_source = system().get_field_source();
      field_source->set_effect_waveform(waveform);

      break;
    }
  case SolverSpecify::TRACE       :
    {
      // user should specify vscan here
      genius_assert(c.is_parameter_exist("vscan"));
      {
        unsigned int elec_num = c.parameter_count("vscan");
        for(unsigned int n=0; n<elec_num; n++)
        {
          std::string electrode = c.get_n_string("vscan", "", n, 0);
          if( system().get_bcs()->get_bc(electrode) == NULL || system().get_bcs()->get_bc(electrode)->is_electrode() == false )
          {
            MESSAGE<<"ERROR at " <<c.get_fileline()<< " SOLVE: Electrode " << electrode << " can't be found in device structure." << std::endl; RECORD();
            genius_error();
          }
          SolverSpecify::Electrode_VScan.push_back(electrode);
        }

        if( !SolverSpecify::Electrode_VScan.size() )
        {
          MESSAGE<<"ERROR: You must specify at least one electrode for IV trace."<<std::endl; RECORD();
          genius_error();
        }
      }
      SolverSpecify::VStart    = c.get_real("vstart", 0.0)*V;
      SolverSpecify::VStep     = c.get_real("vstep", 0.1)*V;
      SolverSpecify::VStepMax  = c.get_real("vstepmax", SolverSpecify::VStep/V)*V;
      SolverSpecify::VStop     = c.get_real("vstop", 5.0)*V;
      SolverSpecify::IStop     = c.get_real("istop", 1.0)*A; //current limit

      SolverSpecify::Predict   = c.get_bool("predict", true);

      SolverSpecify::OptG      = c.get_bool("optical.gen", false);
      SolverSpecify::PatG      = c.get_bool("particle.gen", false);

      // set the waveform of light source
      std::string waveform = c.get_string("waveform", "");
      FieldSource * field_source = system().get_field_source();
      field_source->set_effect_waveform(waveform);

      break;
    }

  case SolverSpecify::ACSWEEP     :
    {
      SolverSpecify::Electrode_ACScan.clear();

      SolverSpecify::FStart    = c.get_real("f.start", 1e6)/s;
      SolverSpecify::FStop     = c.get_real("f.stop", 10e9)/s;
      SolverSpecify::FMultiple = c.get_real("f.multiple", 1.1);
      SolverSpecify::VAC       = c.get_real("vac", 0.0026)*V;

      unsigned int elec_num = c.parameter_count("acscan");
      for(unsigned int n=0; n<elec_num; n++)
      {
        std::string electrode = c.get_n_string("acscan", "", n, 0);
        if( system().get_bcs()->get_bc(electrode) == NULL || system().get_bcs()->get_bc(electrode)->is_electrode() == false )
        {
          MESSAGE<<"ERROR at " <<c.get_fileline()<< " SOLVE: Electrode " << electrode << " can't be found in device structure." << std::endl; RECORD();
          genius_error();
        }
        SolverSpecify::Electrode_ACScan.push_back(electrode);
      }

      if( !SolverSpecify::Electrode_ACScan.size() )
      {
        MESSAGE<<"ERROR: You must specify at least one electrode for AC scan."<<std::endl; RECORD();
        genius_error();
      }

      SolverSpecify::Type = SolverSpecify::ACSWEEP;

      break;
    }

  case SolverSpecify::TRANSIENT  :
    {
      SolverSpecify::AutoStep  = c.get_bool("autostep", true);
      SolverSpecify::Predict   = c.get_bool("predict", true);
      SolverSpecify::UIC       = c.get_bool("uic", false);

      SolverSpecify::TStart    = c.get_real("tstart", 0.0)*s;
      SolverSpecify::TStep     = c.get_real("tstep", 1e-9)*s;
      SolverSpecify::TStepMax  = c.get_real("tstepmax", SolverSpecify::TStep/s)*s;

      SolverSpecify::TStop     = c.get_real("tstop", 1e-6)*s;

      SolverSpecify::TS_rtol   = c.get_real("ts.rtol", 1e-3);
      SolverSpecify::TS_atol   = c.get_real("ts.atol", 1e-4);

      if(c.is_parameter_exist("ts"))
      {
        if (c.is_enum_value("ts", "impliciteuler"))   SolverSpecify::TS_type = SolverSpecify::BDF1;
        if (c.is_enum_value("ts", "bdf1"))            SolverSpecify::TS_type = SolverSpecify::BDF1;
        if (c.is_enum_value("ts", "bdf2"))            SolverSpecify::TS_type = SolverSpecify::BDF2;
      }

      SolverSpecify::OptG      = c.get_bool("optical.gen", false);
      SolverSpecify::PatG      = c.get_bool("particle.gen", false);

      // set the waveform of light source
      std::string waveform = c.get_string("waveform", "");
      FieldSource * field_source = system().get_field_source();
      field_source->set_effect_waveform(waveform);

      break;
    }

  default: break;

  }

  SolverSpecify::out_prefix = c.get_string("out.prefix", "result");

  SolverBase * solver = NULL;

  // call each solver here
  switch (SolverSpecify::Solver)
  {
  case SolverSpecify::POISSON :
    {
      solver = new PoissonSolver(system());
      break;
    }
  case SolverSpecify::DDML1 :
    {
      // if spice circuit exist, use MixA1Solver
      solver = new DDM1Solver(system());
      break;
    }
  case SolverSpecify::DDML1MIXA :
    {
      solver = new MixA1Solver(system());
      break;
    }
  case SolverSpecify::HALLDDML1 :
    {
      solver = new HallSolver(system());
      break;
    }
  case SolverSpecify::DDML2 :
    {
      solver = new DDM2Solver(system());
      break;
    }
  case SolverSpecify::DDML2MIXA :
    {
      solver = new MixA2Solver(system());
      break;
    }
  case SolverSpecify::EBML3 :
    {
      solver = new EBM3Solver(system());
      break;
    }
  case SolverSpecify::EBML3MIXA:
    {
      solver = new MixA3Solver(system());
      break;
    }
  case SolverSpecify::DDMAC :
    {
      solver = new DDMACSolver(system());
      break;
    }
#ifndef CYGWIN
  case SolverSpecify::DDML1MIX :
    {
      solver = new Mix1Solver(system());
      break;
    }
  case SolverSpecify::DDML2MIX :
    {
      solver = new Mix2Solver(system());
      break;
    }
  case SolverSpecify::EBML3MIX :
    {
      solver = new Mix3Solver(system());
      break;
    }
#endif

  default: break;
    MESSAGE<<"ERROR: Selected solver is not supported at present." << std::endl; RECORD();
    break;
  }

  if (solver)
  {
    solver->set_label(SolverSpecify::label);

    // create a solution group;
    mxml_node_t *eGroup = NULL;
    {
      mxml_node_t *eRoot = mxmlFindElement(_dom_solution, _dom_solution, "genius-solutions", NULL, NULL, MXML_DESCEND_FIRST);
      eGroup = mxmlNewElement(eRoot, "solution-group");
      mxml_node_t *eLabel = mxmlNewElement(eGroup, "label");
      mxmlAdd(eLabel, MXML_ADD_AFTER, NULL, MXMLQVariant::makeQVString(solver->label()));

      solver->set_solution_dom_root(eGroup);
    }

    // init (user defined) hook functions here

    if( SolverSpecify::Type == SolverSpecify::DCSWEEP   ||
        SolverSpecify::Type == SolverSpecify::TRANSIENT ||
        SolverSpecify::Type == SolverSpecify::TRACE     ||
        SolverSpecify::Solver == SolverSpecify::DDMAC
      )
    {
#ifdef CYGWIN
      // for windows cygwin system, dynamic link is not supported. we have to use static link.
      // it is not as flexible as unix/linux system().
      // rawfile hook, write electrode IV in SPICE raw file format
      Hook * rawfile_hook =  new RawFileHook(*solver, "rawfile_hook", (void *)Genius::input_file());
      solver->add_hook(rawfile_hook);
      // gnuplot hook, write electrode IV in gnuplot file format
      Hook * gnuplot_hook =  new GnuplotHook(*solver, "gnuplot_hook", (void *)Genius::input_file());
      solver->add_hook(gnuplot_hook);
#else
      Hook * rawfile_hook =  new DllHook(*solver, "rawfile_hook", (void *)(Genius::input_file()));
      solver->add_hook(rawfile_hook);

      Hook * gnuplot_hook =  new DllHook(*solver, "gnuplot_hook", (void *)(Genius::input_file()));
      solver->add_hook(gnuplot_hook);
#endif

    }

#ifdef CYGWIN
    // load static user defined hooks, only support predefined hooks, sigh
    for (std::deque<std::string>::iterator it=SolverSpecify::Hooks.begin();
           it!=SolverSpecify::Hooks.end(); it++)
      {
        Hook * hook=NULL;

        if((*it)=="vtk")
          hook = new VTKHook(*solver, "vtk_hook", (void *)(Genius::input_file()));
        if((*it)=="cv")
          hook = new CVHook (*solver, "cv_hook",  (void *)(Genius::input_file()));
        if((*it)=="probe")
          hook = new ProbeHook (*solver, "probe_hook",  (void *)(Genius::input_file()));

        if(hook) solver->add_hook(hook);
      }

#else
    // dynamic load user defined hooks, stupid win32 platform does not support this function.
    for (std::deque<std::string>::iterator it=SolverSpecify::Hooks.begin();
           it!=SolverSpecify::Hooks.end(); it++)
      {
#ifdef ENABLE_VISIT
        if((*it)=="visit")
          solver->add_hook( new VisitHook (*solver, "visit_hook",  (void *)(Genius::input_file())) );
        else
#endif
          solver->add_hook( new DllHook(*solver, (*it)+"_hook", (void *)(Genius::input_file())) );
      }
#endif

    {
      // always load the control hook. We load it last, such that it is called last
      SolverControlHook * control_hook =  new SolverControlHook(*solver, "control_hook", *this, _fname_solution);
      solver->add_hook(control_hook);
    }

    solver->create_solver();
    solver->solve();
    solver->destroy_solver();

    {
      // if there is a solution in the group, add it to the solution document
      if (mxmlFindElement(eGroup, eGroup, "solution", NULL, NULL, MXML_DESCEND_FIRST)==NULL)
      {
        mxmlDelete(eGroup);
      }
    }

    delete solver;

  }

  return 0;
}




int  SolverControl::set_electrode_source  ( const Parser::Card & c )
{

  std::string electrode = c.get_string("electrode", "");

  if( system().get_bcs()->get_bc(electrode) == NULL || system().get_bcs()->get_bc(electrode)->is_electrode() == false )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " ATTACH: Electrode " << electrode << " can't be found in device structure." << std::endl; RECORD();
    genius_error();
  }

  bool voltage = true;

  if(c.is_parameter_exist("type"))
  {
    if (c.is_enum_value("type", "voltage"))   voltage = true;
    if (c.is_enum_value("type", "current"))   voltage = false;
  }

  if( c.is_parameter_exist("vapp") || c.is_parameter_exist("iapp") )
  {

    std::vector<std::string> source_list;

    // since several vapp or iapp parameters may exist in the card, we should search for all
    for(unsigned int idx=0; idx<c.parameter_size(); idx++)
    {
      Parser::Parameter p = c.get_parameter(idx);
      if ( p.name() == "vapp" && voltage == true)
      {
        if(system().get_sources()->is_vsource_exist(p.get_string()))
          source_list.push_back( p.get_string() );
        else
        {
          MESSAGE<<"ERROR at " <<c.get_fileline()<< " ATTACH: Vapp " << p.get_string() << " hasn't been defined." << std::endl; RECORD();
          genius_error();
        }
      }

      if ( p.name() == "iapp" && voltage == false)
      {
        if(system().get_sources()->is_isource_exist(p.get_string()))
          source_list.push_back( p.get_string() );
        else
        {
          MESSAGE<<"ERROR at " <<c.get_fileline()<< " ATTACH: Iapp " << p.get_string() << " hasn't been defined." << std::endl; RECORD();
          genius_error();
        }
      }
    }

    system().get_sources()->attach_sources_to_electrode(electrode, source_list);
    return 0;
  }

  if( c.is_parameter_exist("vconst") )
  {
    system().get_sources()->attach_voltage_to_electrode(electrode, c.get_real("vconst", 0.0)*V);
    return 0;
  }

  if ( c.is_parameter_exist("iconst") )
  {
    system().get_sources()->attach_current_to_electrode(electrode, c.get_real("iconst", 0.0)*A);
    return 0;
  }

  return 0;
}





int  SolverControl::set_physical_model  ( const Parser::Card & c )
{

  std::string region_label = c.get_string("region", "");
  std::string type         = c.get_string("type", "");
  std::string model        = c.get_string("model", "Default");

  std::vector<Parser::Parameter> pmi_parameters;

  for(unsigned int idx=0; idx<c.parameter_size(); idx++)
  {
    Parser::Parameter p = c.get_parameter(idx);
    // find user defined parameter, which used to calibrate the PMI
    if ( p.is_user_defined() )
    {
      pmi_parameters.push_back(p);
    }
  }

  if( system().region(region_label)==NULL )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " PMI: Region " << region_label << " can't be found in mesh regions." << std::endl; RECORD();
    genius_error();
  }

  if ( (type.length() > 0  ) )
  {
    system().region(region_label)->set_pmi(type, model, pmi_parameters);
    system().get_bcs()->pmi_init_bc(region_label,type);
  }
  else
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " PMI: Must specify the type parameter." << std::endl; RECORD();
    genius_error();
  }

  int print_verbosity = c.get_int("print", 0);
  if ( print_verbosity > 0 )
  {
    // verbose output, let's print out the new material parameters
    MESSAGE << system().region(region_label)->get_pmi_info(type, print_verbosity) << std::endl;
  }

  return 0;

}




int SolverControl::do_export( const Parser::Card & c )
{
  // if export to VTK format is required
  if(c.is_parameter_exist("vtkfile"))
  {
    std::string vtk_filename = c.get_string("vtkfile", "");
    bool ascii = c.get_bool("ascii", false);
    system().export_vtk(vtk_filename, ascii);
  }

  // if export to CGNS format is required
  if(c.is_parameter_exist("cgnsfile"))
  {
    std::string cgns_filename = c.get_string("cgnsfile", "");
    system().export_cgns(cgns_filename);
  }

  // if export GDML surface file is required
  if(c.is_parameter_exist("gdml.surface"))
  {
    std::string gdml_filename = c.get_string("gdml.surface", "");
    system().export_gdml_surface(gdml_filename);
  }

  if(c.is_parameter_exist("gdml.body"))
  {
    std::string gdml_filename = c.get_string("gdml.body", "");
    system().export_gdml_body(gdml_filename);
  }

  // if export boundary condition is required
  if(c.is_parameter_exist("bcinfo"))
  {
    std::string bc_filename = c.get_string("bcinfo", "");
    system().get_bcs()->export_boundary_condition(bc_filename);
  }

  // if export node location is required
  if(c.is_parameter_exist("nodeinfo"))
  {
    std::string node_filename = c.get_string("nodeinfo", "");
    bool numbering = c.get_bool("numbering", true, "");

    if (c.is_enum_value("lunit", "m"))
      system().export_node_location(node_filename, 1.0, numbering );
    else if (c.is_enum_value("lunit", "cm"))
      system().export_node_location(node_filename, 1e-2, numbering );
    else if (c.is_enum_value("lunit", "um"))
      system().export_node_location(node_filename, 1e-6, numbering );
    else if (c.is_enum_value("lunit", "nm"))
      system().export_node_location(node_filename, 1e-9, numbering );
    else
      system().export_node_location(node_filename, 1e-6, numbering );
  }

  return 0;
}



int SolverControl::do_import( const Parser::Card & c )
{
  if(c.is_parameter_exist("cgnsfile"))
  {
    std::string cgns_filename = c.get_string("cgnsfile", "");
#ifdef CYGWIN
    if ( _access( (char *)cgns_filename.c_str(),  04 ) == -1 )
#else
    if (  access( (char *)cgns_filename.c_str(),  R_OK ) == -1 )
#endif
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " IMPORT: CGNSFile " << cgns_filename << " doesn't exist." << std::endl; RECORD();
      genius_error();
    }
    system().import_cgns(cgns_filename);
  }

  if(c.is_parameter_exist("vtkfile"))
  {
    std::string vtk_filename = c.get_string("vtkfile", "");
#ifdef CYGWIN
    if ( _access( (char *)vtk_filename.c_str(),  04 ) == -1 )
#else
    if (  access( (char *)vtk_filename.c_str(),  R_OK ) == -1 )
#endif
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " IMPORT: VTKFile " << vtk_filename << " doesn't exist." << std::endl; RECORD();
      genius_error();
    }
    system().import_vtk(vtk_filename);
  }

  if(c.is_parameter_exist("silvacofile"))
  {
    std::string silvaco_filename = c.get_string("silvacofile", "");
#ifdef CYGWIN
    if ( _access( (char *)silvaco_filename.c_str(),  04 ) == -1 )
#else
    if (  access( (char *)silvaco_filename.c_str(),  R_OK ) == -1 )
#endif
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " IMPORT: Silvaco File " << silvaco_filename << " doesn't exist." << std::endl; RECORD();
      genius_error();
    }
    system().import_silvaco(silvaco_filename);
  }


  if(c.is_parameter_exist("tiffile"))
  {
    std::string tif_filename = c.get_string("tiffile", "");
#ifdef CYGWIN
    if ( _access( (char *)tif_filename.c_str(),  04 ) == -1 )
#else
    if (  access( (char *)tif_filename.c_str(),  R_OK ) == -1 )
#endif
    {
      MESSAGE<<"ERROR at " <<c.get_fileline()<< " IMPORT: TIFFile " << tif_filename << " doesn't exist." << std::endl; RECORD();
      genius_error();
    }
    system().import_tif(tif_filename);
  }


  if(c.is_parameter_exist("isefile"))
  {
    std::string ise_filename = c.get_string("isefile", "");
    system().import_ise(ise_filename);
  }

  return 0;
}




int  SolverControl::do_em_fem2d_solve ( const Parser::Card & c )
{
  EMFEM2DSolver * solver = new EMFEM2DSolver(system(), c);
  solver->create_solver();
  solver->solve();
  solver->destroy_solver();
  delete solver;

  return 0;
}


int SolverControl::do_ray_trace( const Parser::Card & c )
{

  RayTraceSolver * solver = new RayTraceSolver(system(), c);
  solver->create_solver();
  solver->solve();
  solver->destroy_solver();
  delete solver;

  return 0;
}


int  SolverControl::set_initial_node_voltage  ( const Parser::Card & c )
{

  std::string electrode = c.get_string("electrode", "");

  BoundaryCondition * bc = system().get_bcs()->get_bc(electrode);
  if( bc == NULL || bc->is_electrode() == false )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline()<< " NODESET: Electrode " << electrode << " can't be found in device structure." << std::endl; RECORD();
    genius_error();
  }

  bc->ext_circuit()->potential() = c.get_real("v", 0.0);

  return 0;
}




int SolverControl::do_refine_conform(const Parser::Card & c)
{
  // TODO can we refine during the solver solution processing?

  // save previous solution
  AutoPtr<InterpolationBase> interpolator;
  if( mesh().magic_num() < 2008 )
    interpolator = AutoPtr<InterpolationBase>(new Interpolation2D_CSA);
  else
    interpolator = AutoPtr<InterpolationBase>(new Interpolation3D_nbtet);

  if( DopingSolver.get() == NULL )
  {
    system().fill_interpolator(interpolator.get(), "doping.na", InterpolationBase::Asinh);
    system().fill_interpolator(interpolator.get(), "doping.nd", InterpolationBase::Asinh);
  }

  if(system().has_single_compound_semiconductor_region()  && MoleSolver.get() == NULL )
  {
    system().fill_interpolator(interpolator.get(), "mole.x", InterpolationBase::Linear);
  }
  if(system().has_complex_compound_semiconductor_region()  && MoleSolver.get() == NULL )
  {
    system().fill_interpolator(interpolator.get(), "mole.y", InterpolationBase::Linear);
  }

  // fill error vector from system level
  ErrorVector error_per_cell;
  system().estimate_error(c, error_per_cell);

  if (Genius::processor_id() == 0)
  {

    MeshRefinement mesh_refinement(mesh());

    // at least one refine criterion should be exist!
    genius_assert(c.is_parameter_exist("error.fraction") || c.is_parameter_exist("cell.fraction") || c.is_parameter_exist("error.threshold"));

    if(c.is_parameter_exist("error.fraction") )
      mesh_refinement.flag_elements_by_error_fraction (error_per_cell, c.get_real("error.fraction", 0.3), 0.0);

    if(c.is_parameter_exist("cell.fraction") )
      mesh_refinement.flag_elements_by_elem_fraction  (error_per_cell, c.get_real("cell.fraction",  0.3), 0.0);

    if(c.is_parameter_exist("error.threshold") )
      mesh_refinement.flag_elements_by_error_threshold(error_per_cell, c.get_real("error.threshold",0.1), 0.0);

    // if mesh generator exist, we call it to do particular refine
    if( meshgen.get() != NULL )
      meshgen->do_refine(mesh_refinement);
    // else, we have to do general mesh refinement
    else
    {
      genius_assert(mesh().magic_num() != invalid_uint);

      // for 2D, we call triangle
      if( mesh().magic_num() < 2008 )
      {
        MeshGenerator *meshgen = new MeshGeneratorTri3(mesh(), decks());
        meshgen->do_refine(mesh_refinement);
        delete meshgen;
      }
      // for 3D, we still have no idea here.
      else
      {
        MESSAGE<<"ERROR at " <<c.get_fileline()<< " Refine: Genius still can not do 3D conform refine without mesh generator exist." << std::endl; RECORD();
        genius_error();
      }
    }

  }

  // clear the system. however we should reserve mesh information
  system().clear(false);

  // call mesh generator again to make new mesh consistence with FVM request.

  // rebuild the system
  // since we only build mesh on processor 0,
  // sync mesh to other processors.
  // this procedure also prepare the mesh for using
  MeshCommunication mesh_comm;
  mesh_comm.broadcast(mesh());

  // now we can build solution system again
  system().build_simulation_system();
  system().sync_print_info();

  // set doping profile to semiconductor region
  if( DopingSolver.get() != NULL )
    DopingSolver->solve();
  else
  {
    // no doping information?
    system().do_interpolation(interpolator.get(), "doping.na");
    system().do_interpolation(interpolator.get(), "doping.nd");
  }

  // set mole fraction to semiconductor region
  if( MoleSolver.get() != NULL )
    MoleSolver->solve();
  else
  {
    if(system().has_single_compound_semiconductor_region())
      system().do_interpolation(interpolator.get(), "mole.x");
    if(system().has_complex_compound_semiconductor_region())
      system().do_interpolation(interpolator.get(), "mole.y");
  }

  // after doping profile is set, we can init system data.
  system().init_region();

  return 0;

}



int SolverControl::do_refine_hierarchical(const Parser::Card & c)
{

  MESSAGE<<"Hierarchical mesh refinement...\n"<<std::endl; RECORD();

  // save previous solution
  AutoPtr<InterpolationBase> interpolator;
  if( mesh().magic_num() < 2008 )
    interpolator = AutoPtr<InterpolationBase>(new Interpolation2D_CSA);
  else
    interpolator = AutoPtr<InterpolationBase>(new Interpolation3D_nbtet);

  if( DopingSolver.get() == NULL )
  {
    system().fill_interpolator(interpolator.get(), "doping.na", InterpolationBase::Asinh);
    system().fill_interpolator(interpolator.get(), "doping.nd", InterpolationBase::Asinh);
  }

  if(system().has_single_compound_semiconductor_region()  && MoleSolver.get() == NULL )
  {
    system().fill_interpolator(interpolator.get(), "mole.x", InterpolationBase::Linear);
  }
  if(system().has_complex_compound_semiconductor_region()  && MoleSolver.get() == NULL )
  {
    system().fill_interpolator(interpolator.get(), "mole.y", InterpolationBase::Linear);
  }

  // fill error vector from system level
  ErrorVector error_per_cell;
  system().estimate_error(c, error_per_cell);

  if (Genius::processor_id() == 0)
  {

    MeshRefinement mesh_refinement(mesh());

    // at least one refine criterion should be exist!
    genius_assert(c.is_parameter_exist("error.refine.fraction") || c.is_parameter_exist("cell.refine.fraction") || c.is_parameter_exist("error.refine.threshold"));

    if(c.is_parameter_exist("error.refine.fraction") )
      mesh_refinement.flag_elements_by_error_fraction (error_per_cell, c.get_real("error.refine.fraction", 0.3), c.get_real("error.coarsen.fraction", 0.0));

    if(c.is_parameter_exist("cell.refine.fraction") )
      mesh_refinement.flag_elements_by_elem_fraction  (error_per_cell, c.get_real("cell.refine.fraction",  0.3), c.get_real("cell.coarsen.fraction",  0.0));

    if(c.is_parameter_exist("error.refine.threshold") )
      mesh_refinement.flag_elements_by_error_threshold(error_per_cell, c.get_real("error.refine.threshold",0.1), c.get_real("error.coarsen.threshold",0.0));

    // call MeshRefinement class to do FEM refine
    mesh_refinement.refine_and_coarsen_elements ();
  }

  // clear the system(). however we should reserve mesh information
  system().clear(false);

  // call mesh generator again to make new mesh consistence with FVM request.

  // rebuild the system
  // since we only build mesh on processor 0,
  // sync mesh to other processors.
  // this procedure also prepare the mesh for using
  MeshCommunication mesh_comm;
  mesh_comm.broadcast(mesh());

  // now we can build solution system again
  system().build_simulation_system();
  system().sync_print_info();

  // set doping profile to semiconductor region
  if( DopingSolver.get() != NULL )
    DopingSolver->solve();
  else
  {
    // no doping information?
    system().do_interpolation(interpolator.get(), "doping.na");
    system().do_interpolation(interpolator.get(), "doping.nd");
  }

  // set mole fraction to semiconductor region
  if( MoleSolver.get() != NULL )
    MoleSolver->solve();
  else
  {
    if(system().has_single_compound_semiconductor_region())
      system().do_interpolation(interpolator.get(), "mole.x");
    if(system().has_complex_compound_semiconductor_region())
      system().do_interpolation(interpolator.get(), "mole.y");
  }

  // after doping profile is set, we can init system data.
  system().init_region();

  return 0;

}


/*--------------------------------------------------------------------
 * uniform refine the mesh, it is intended to be used only for test!
 */
int SolverControl::do_refine_uniform(const Parser::Card & c)
{

  if (Genius::processor_id() == 0)
  {
    int step =  c.get_int("step", 1);
    MeshRefinement mesh_refinement(mesh());
    mesh_refinement.uniformly_refine(step);
    MeshTools::Modification::flatten(mesh());
  }

  // clear the system. however we should reserve mesh information
  system().clear(false);

  // since we only build mesh on processor 0,
  // sync mesh to other processors.
  // this procedure also prepare the mesh for using
  MeshCommunication mesh_comm;
  mesh_comm.broadcast(mesh());

  // now we can build solution system again
  system().build_simulation_system();
  system().sync_print_info();

  // set doping profile to semiconductor region
  if( DopingSolver.get() != NULL )
    DopingSolver->solve();

  // set mole fraction to semiconductor region
  if( MoleSolver.get() != NULL )
    MoleSolver->solve();

  // after doping profile is set, we can init system data.
  system().init_region();

  return 0;
}


int SolverControl::extend_to_3d ( const Parser::Card & c )
{
  MESSAGE<<"Extend mesh to 3D...\n"<<std::endl; RECORD();
  ExtendTo3D(system(), c)();
  return 0;
}


int SolverControl::plot_mesh(const Parser::Card & c)
{
#ifdef HAVE_X11
  bool inv = c.get_bool("y.inverse", true);
  std::string tiff = c.get_string("tiff.out", "");

  if( system().mesh().mesh_dimension() == 2)
  {
    ShowMesh2D  show_mesh(system(),inv);
    if(Genius::is_first_processor() && !show_mesh.show_mesh_init())
      show_mesh.show_mesh(tiff=="" ? NULL : tiff.c_str());
  }
#endif
  return 0;
}

// ------
SolverControlHook::SolverControlHook(SolverBase & solver, const std::string & name, SolverControl& control, const std::string& fname)
    : Hook(solver, name), _control(control), _fname(fname)
{}

SolverControlHook::~SolverControlHook()
{}

void SolverControlHook::on_init()
{}

void SolverControlHook::on_close()
{}

void SolverControlHook::pre_solve()
{}

void SolverControlHook::post_solve()
{
  if (Genius::processor_id()==0)
  {
    if (!_fname.empty())
    {
      FILE *fout = fopen(_fname.c_str(), "w");
      mxmlSaveFile(_control.get_dom_solution(), fout, MXML_NO_CALLBACK);
      fclose(fout);
    }
  }
}

void SolverControlHook::post_iteration()
{}

