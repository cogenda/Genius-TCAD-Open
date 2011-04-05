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

//  $Id: boundary_condition_collector.cc,v 1.28 2008/07/09 05:58:16 gdiso Exp $
#include <fstream>

#include "mesh_base.h"
#include "boundary_info.h"
#include "fvm_node_info.h"
#include "fvm_node_data.h"
#include "simulation_system.h"
#include "material.h"
#include "boundary_condition_collector.h"
#include "mxml.h"
#include "MXMLUtil.h"


// all the derived boundary conditions
#include "boundary_condition_abs.h"
#include "boundary_condition_ii.h"
#include "boundary_condition_charge.h"
#include "boundary_condition_ir.h"
#include "boundary_condition_is.h"
#include "boundary_condition_ei.h"
#include "boundary_condition_neumann.h"
#include "boundary_condition_gate.h"
#include "boundary_condition_ohmic.h"
#include "boundary_condition_resistance_ohmic.h"
#include "boundary_condition_rr.h"
#include "boundary_condition_schottky.h"
#include "boundary_condition_resistance_schottky.h"
#include "boundary_condition_hetero.h"
#include "boundary_condition_simplegate.h"
#include "boundary_condition_homo.h"
#include "boundary_condition_solderpad.h"

#include "boundary_condition_electrode_interconnect.h"
#include "boundary_condition_charge_integral.h"

//needed physical unit & constant
using PhysicalUnit::J;
using PhysicalUnit::K;
using PhysicalUnit::eps0;
using PhysicalUnit::mu0;
using PhysicalUnit::s;
using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::C;




BoundaryConditionCollector::BoundaryConditionCollector ( SimulationSystem & system, Parser::InputParser & decks )
  : _decks ( decks ), _system ( system ), _mesh ( system.mesh() ), _resistive_metal_mode(system.resistive_metal_mode())
{}



BoundaryCondition * BoundaryConditionCollector::get_bc_nocase(const std::string & label)
{
  std::map<const std::string, BoundaryCondition * >::iterator it = _bc_label_to_bc.begin();
  for(; it != _bc_label_to_bc.end(); ++it)
  {
#ifndef CYGWIN
      if( !strcasecmp(label.c_str(), (*it).first.c_str() ) )
        return (*it).second;
#else
      if(label.size()!=(*it).first.size()) continue;
      bool equal_string = true;
      for(unsigned int n=0; n<label.size(); ++n)
      {
        char n1 =  label.at(n);
        char n2 =  (*it).first.at(n);
        if(n1==n2) continue;
        if( isalpha(n1) && isalpha(n2) && toupper(n1)==toupper(n2)) continue;
        equal_string = false;
        break;
      }

      if(equal_string) return (*it).second;
#endif

  }
  return NULL;
}


const BoundaryCondition * BoundaryConditionCollector::get_bc_nocase(const std::string & label) const
{
  std::map<const std::string, BoundaryCondition * >::const_iterator it = _bc_label_to_bc.begin();
  for(; it != _bc_label_to_bc.end(); ++it)
  {
#ifndef CYGWIN
      if( !strcasecmp(label.c_str(), (*it).first.c_str() ) )
        return (*it).second;
#else
      if(label.size()!=(*it).first.size()) continue;
      bool equal_string = true;
      for(unsigned int n=0; n<label.size(); ++n)
      {
        char n1 =  label.at(n);
        char n2 =  (*it).first.at(n);
        if(n1==n2) continue;
        if( isalpha(n1) && isalpha(n2) && toupper(n1)==toupper(n2)) continue;
        equal_string = false;
        break;
      }

      if(equal_string) return (*it).second;
#endif

  }
  return NULL;
}


void BoundaryConditionCollector::pmi_init_bc ( const std::string &region_label, const std::string &type )
{
  for ( unsigned int i=0; i<_bcs.size(); i++ )
  {
    for ( BoundaryCondition::const_node_iterator node_it = _bcs[i]->nodes_begin();
          node_it != _bcs[i]->nodes_end();
          node_it ++ )
    {

      for ( BoundaryCondition::region_node_iterator rnode_it = _bcs[i]->region_node_begin ( *node_it );
            rnode_it != _bcs[i]->region_node_end ( *node_it );
            rnode_it ++ )
      {
        SimulationRegion *region = ( *rnode_it ).second.first;
        FVM_Node * fvm_node      = ( *rnode_it ).second.second;

        if ( region->name() != region_label ) continue;

        // when node_data is not NULL, this node belongs to current process
        // or at least it is a ghost node
        FVM_NodeData * node_data = fvm_node->node_data();
        if ( node_data == NULL ) continue;

        region->get_material_base()->init_bc_node ( type,_bcs[i]->label(),fvm_node->root_node(),node_data );
      }
    }
  }
}



unsigned int BoundaryConditionCollector::get_bc_index_by_bd_id(short int bd_id) const
{
  std::map<short int, unsigned int>::const_iterator it =  _bd_id_to_bc_index.find(bd_id);
  genius_assert( it != _bd_id_to_bc_index.end() );
  return (*it).second;
}


short int BoundaryConditionCollector::get_bd_id_by_bc_index(unsigned int bc_index) const
{
  std::map<unsigned int, short int>::const_iterator it =  _bc_index_to_bd_id.find(bc_index);
  if( it != _bc_index_to_bd_id.end() )
    return (*it).second;
  return BoundaryInfo::invalid_id;
}


int  BoundaryConditionCollector::bc_setup()
{

  MESSAGE<<"Bulid boundary conditions on all processors..."<<std::endl;  RECORD();

  // the number of boundarys
  unsigned int n_bc = _mesh.boundary_info->n_boundary_ids ();
  _bcs.resize ( n_bc, NULL );

  // set the boundary condition index <--> boundary id map
  build_bc_index_to_bd_id_map ( _mesh.boundary_info->get_boundary_ids() );

  // setup each boundary condition from user's input
  for ( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();

    // the BC card can be classfied into BOUNDARY and CONTACT.
    // BOUNDARY specifies the boundary condition to "FACE"
    // CONTACT specifies the boundary condition to electrode region

    if ( c.key() == "BOUNDARY" )  //  boundary card
    {

      BCType bc_type = BC_string_to_enum ( c.get_string ( "type", "" ) );

      switch ( bc_type )
      {
        case NeumannBoundary             :  { if ( Set_BC_NeumannBoundary ( c ) )              return 1; break;}
        case OhmicContact                :  { if ( Set_BC_OhmicContact ( c ) )                 return 1; break;}
        case IF_Metal_Ohmic              :  { if ( Set_BC_IF_Metal_Ohmic ( c ) )               return 1; break;}
        case SchottkyContact             :  { if ( Set_BC_SchottkyContact ( c ) )              return 1; break;}
        case IF_Metal_Schottky           :  { if ( Set_BC_IF_Metal_Schottky ( c ) )            return 1; break;}
        case GateContact                 :  { if ( Set_BC_GateContact ( c ) )                  return 1; break;}
        case SimpleGateContact           :  { if ( Set_BC_SimpleGateContact ( c ) )            return 1; break;}
        case SolderPad                   :  { if ( Set_BC_Solderpad ( c ) )                    return 1; break;}
        case IF_Insulator_Semiconductor  :  { if ( Set_BC_InsulatorInterface ( c ) )           return 1; break;}
        case AbsorbingBoundary           :  { if ( Set_BC_AbsorbingBoundary ( c ) )            return 1; break;}
        default: genius_error();
      }
    }

    if ( c.key() == "CONTACT" )  //  contact card
    {
      BCType bc_type = BC_string_to_enum ( c.get_string ( "type", "" ) );

      switch ( bc_type )
      {
        case OhmicContact                :  { if ( Set_Electrode_OhmicContact ( c ) )       return 1; break;}
        case SchottkyContact             :  { if ( Set_Electrode_SchottkyContact ( c ) )    return 1; break;}
        case GateContact                 :  { if ( Set_Electrode_GateContact ( c ) )        return 1; break;}
        case ChargedContact              :  { if ( Set_SetFloatMetal ( c ) )                return 1; break;}
        default: genius_error();
      }

    }

  }

  // if we have boundary description in mesh information, set boundary condition by description
  for( unsigned int i=0; i<_bcs.size(); i++)
    if( _bcs[i] == NULL ) // not initialized?
    {
      short int boundary_id = get_bd_id_by_bc_index(i);
      if( boundary_id == BoundaryInfo::invalid_id ) continue;

      std::string description = _mesh.boundary_info->get_description_by_id(boundary_id);
      if(description.empty()) continue;

      Parser::Card c(description);
      if(c.key() != "BOUNDARY") continue;
      BCType bc_type = BC_string_to_enum(c.get_string("type", ""));

      switch(bc_type)
      {
        case NeumannBoundary             :  { if( Set_BC_NeumannBoundary(c) )                  return 1; break;}
        case OhmicContact                :  { if( Set_BC_OhmicContact(c) )                     return 1; break;}
        case IF_Metal_Ohmic              :  { if( Set_BC_IF_Metal_Ohmic ( c ) )                return 1; break;}
        case SchottkyContact             :  { if( Set_BC_SchottkyContact(c) )                  return 1; break;}
        case IF_Metal_Schottky           :  { if( Set_BC_IF_Metal_Schottky ( c ) )             return 1; break;}
        case GateContact                 :  { if( Set_BC_GateContact(c) )                      return 1; break;}
        case SimpleGateContact           :  { if( Set_BC_SimpleGateContact(c) )                return 1; break;}
        case SolderPad                   :  { if( Set_BC_Solderpad(c) )                        return 1; break;}
        case IF_Insulator_Semiconductor  :  { if( Set_BC_InsulatorInterface(c) )               return 1; break;}
        case ChargedContact              :  { if( Set_BC_ChargedContact(c) )                   return 1; break;}
        case AbsorbingBoundary           :  { if( Set_BC_AbsorbingBoundary(c) )                return 1; break;}
        default: break;
      }
    }


  // set remaining bcs
  std::map<short int,  std::pair<unsigned int, unsigned int > > boundary_subdomain_map;
  _mesh.boundary_info->get_subdomains_bd_on(boundary_subdomain_map);
  for ( unsigned int i=0; i<_bcs.size(); i++ )
    if ( _bcs[i] == NULL ) // not initialized?
    {
      short int boundary_id = get_bd_id_by_bc_index ( i );
      std::string label = _mesh.boundary_info->get_label_by_id ( boundary_id );

      const std::pair<unsigned int, unsigned int > & sub_ids = boundary_subdomain_map[boundary_id];
      unsigned int sub_id1 = sub_ids.first, sub_id2=sub_ids.second;

      // on the region outer boundary, should be set to NeumannBC or AbsorbingBoundary
      if ( sub_id2==invalid_uint )
      {
        if ( Material::IsVacuum ( _mesh.subdomain_material ( sub_id1 ) ) )
          _bcs[i] = new AbsorbingBC ( _system, label );
        else
        {
          _bcs[i] = new NeumannBC ( _system, label );
          // the NeumannBC for conductor region is considered as heat sink
          // assign high heat transfer rate to it
          if ( Material::IsConductor ( _mesh.subdomain_material ( sub_id1 ) ) )
            _bcs[i]->Heat_Transfer() = 1e3*J/s/pow ( cm,2 ) /K;
        }
        continue;
      }

      std::string material1 = _mesh.subdomain_material ( sub_id1 );
      std::string material2 = _mesh.subdomain_material ( sub_id2 );


      switch ( determine_bc_by_subdomain ( material1, material2, _resistive_metal_mode ) )
      {
          case IF_Electrode_Semiconductor  :
          {
            _bcs[i] = new OhmicContactBC ( _system, label, true );
            _bcs[i]->build_ext_circuit ( new ExternalCircuit );
            // allow one use electrode label to load the bc
            if ( Material::IsConductor ( material1 ) || Material::IsResistance ( material1 ) )
              _bcs[i]->electrode_label() = _mesh.subdomain_label_by_id ( sub_id1 );
            if ( Material::IsConductor ( material2 ) || Material::IsResistance ( material2 ) )
              _bcs[i]->electrode_label() = _mesh.subdomain_label_by_id ( sub_id2 );
            break;
          }
          case IF_Metal_Semiconductor      :  _bcs[i] = new IF_Metal_OhmicBC ( _system, label );               break;
          case HomoInterface               :  _bcs[i] = new HomoInterfaceBC ( _system, label );                        break;
          case HeteroInterface             :  _bcs[i] = new HeteroInterfaceBC ( _system, label );                      break;
          case IF_Insulator_Semiconductor  :  _bcs[i] = new InsulatorSemiconductorInterfaceBC ( _system, label );      break;
          case IF_Semiconductor_Vacuum     :  _bcs[i] = new NeumannBC ( _system, label );                              break;
          case IF_Insulator_Insulator      :  _bcs[i] = new InsulatorInsulatorInterfaceBC ( _system, label );          break;
          case IF_Insulator_Vacuum         :  _bcs[i] = new NeumannBC ( _system, label );                              break;
          case IF_Electrode_Insulator      :
          {
            // if the electrode has a semicoductor neighbor, set it to ElectrodeInsulatorInterfaceBC
            // else set it to gate
            unsigned int electrode_sbd =  Material::IsConductor ( material1 ) || Material::IsResistance ( material1 ) ? sub_id1 : sub_id2;
            const SimulationRegion * electrode_region = _system.region(electrode_sbd);

            std::vector<unsigned int> neighbor_subdomains;
            _mesh.boundary_info->get_subdomain_neighbors(electrode_sbd, neighbor_subdomains);
            bool electrode_sbd_has_semiconductor_neighbor = false;
            for( unsigned int n=0; n<neighbor_subdomains.size(); ++n )
            {
              if( _system.region(neighbor_subdomains[n])->type() == SemiconductorRegion )
              { electrode_sbd_has_semiconductor_neighbor = true; break; }
            }

            if( electrode_sbd_has_semiconductor_neighbor )
              _bcs[i] = new ElectrodeInsulatorInterfaceBC ( _system, label );
            else
            {
              _bcs[i] = new GateContactBC ( _system, label, true );
              _bcs[i]->electrode_label() = _mesh.subdomain_label_by_id ( electrode_sbd );
              _bcs[i]->Work_Function() = electrode_region->get_affinity(_system.T_external());
              _bcs[i]->build_ext_circuit ( new ExternalCircuit );
            }
            break;
          }
          case IF_Electrode_Electrode      :
          case IF_Electrode_Metal          :
          {
            MESSAGE<<"ERROR: Resistive metal shoule be enabled for this boundary condition!"<<std::endl;  RECORD();
            genius_error();
          }
          case IF_Insulator_Metal          :  _bcs[i] = new ResistanceInsulatorBC ( _system, label );                  break;
          case IF_Metal_Metal              :  _bcs[i] = new ResistanceResistanceBC ( _system, label );                 break;
          case IF_Electrode_Vacuum         :  _bcs[i] = new NeumannBC ( _system, label );                              break;
          case IF_Metal_Vacuum             :  _bcs[i] = new NeumannBC ( _system, label );                              break;
          case IF_PML_Scatter              :  _bcs[i] = 0; break;
          case IF_PML_PML                  :  _bcs[i] = 0; break;
          default: genius_error();
      }

    }



  // after that, set the the _bc_label_to_bc map
  for ( unsigned int i=0; i<_bcs.size(); i++ )
  {
    _bc_label_to_bc[_bcs[i]->label() ] = _bcs[i];

    // allow one use electrode label to load the bc
    if ( ! _bcs[i]->electrode_label().empty() )
      _bc_label_to_bc[_bcs[i]->electrode_label() ] = _bcs[i];
  }

  // set boundary_id to each bc
  for ( unsigned int i=0; i<_bcs.size(); i++ )
  {
    short int boundary_id = this->get_bd_id_by_bc_index ( i );
    if( boundary_id == BoundaryInfo::invalid_id ) continue;

    _bcs[i]->boundary_id() = boundary_id;
  }

  // reset boundary id to the nodes in boundary_info class with a priority order
  std::map<short int, unsigned int>  order;
  for ( unsigned int i=0; i<_bcs.size(); i++ )
  {
    short int boundary_id = this->get_bd_id_by_bc_index ( i );
    if( boundary_id == BoundaryInfo::invalid_id ) continue;

    order.insert ( std::make_pair ( boundary_id, static_cast<unsigned int> ( _bcs[i]->bc_type() ) ) );
  }
  _mesh.boundary_info->build_node_ids_from_priority_order ( order );


  // get node belongs to this boundary
  std::map<short int, std::vector<const Node *> > node_boundary_id_map;
  _mesh.boundary_info->nodes_with_boundary_id ( node_boundary_id_map);

  // get node belongs to which region
  std::map<const Node *, std::set<unsigned int > > node_region_map;
  _mesh.boundary_info->node_in_regions (node_region_map);

  for ( unsigned int i=0; i<_bcs.size(); i++ )
  {
    short int boundary_id = get_bd_id_by_bc_index ( i );
    if( boundary_id == BoundaryInfo::invalid_id ) continue;

    // set regions this bc may have
    {
      const std::pair<unsigned int, unsigned int > & sub_ids = boundary_subdomain_map[boundary_id];
      SimulationRegion * r1 = sub_ids.first  != invalid_uint ? _system.region(sub_ids.first)  : NULL;
      SimulationRegion * r2 = sub_ids.second != invalid_uint ? _system.region(sub_ids.second) : NULL;
      if( r1 && r2 && r1->type() > r2->type() ) std::swap(r1, r2);
      _bcs[i]->set_bc_regions(r1, r2);
    }


    // save nodes which has bd_id into the corresponding bc class
    const std::vector<const Node *> & bd_nodes = node_boundary_id_map[boundary_id];
    for ( unsigned int n=0; n<bd_nodes.size(); n++ )
    {
      const Node * bd_node = bd_nodes[n];
      if( bd_node->processor_id() != Genius::processor_id() ) continue;

      _bcs[i]->add_node ( bd_node );

      // set the corresponding region information and FVM_Node to this node
      const std::set<unsigned int > & regions = node_region_map[bd_node];
      std::set<unsigned int >::const_iterator region_it =  regions.begin();
      for ( ; region_it !=  regions.end(); ++region_it)
      {
        SimulationRegion * region = _system.region ( *region_it );
        FVM_Node * fvm_node;
        // we find a FVM_Node with its root node is bd_node
        if ( ( fvm_node=region->region_fvm_node ( bd_node ) ) != NULL )
          _bcs[i]->insert ( bd_node, region, fvm_node );
      }
    }


    _bcs[i]->set_boundary_id_to_fvm_node();
  }


  // process inter-connect /charge boundary here
  for ( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();
    if ( c.key() == "INTERCONNECT" )  //  InterConnect card
      Set_InterConnect ( c );
    if ( c.key() == "CHARGE" ) //  CHARGE card
      Set_Charge ( c );
  }

  // or if we have extra boundary description in mesh information
  const std::vector<std::string> & extra_boundary_descriptions =  _mesh.boundary_info->extra_descriptions();
  for(unsigned int n=0; n<extra_boundary_descriptions.size(); ++n)
  {
    const std::string & description = extra_boundary_descriptions[n];
    if(description.empty()) continue;

    Parser::Card c(description);
    std::string Identifier = c.get_string ( "id", "" );
    if ( c.key() == "INTERCONNECT" && this->get_bc(Identifier)==NULL )  //  InterConnect card
      Set_InterConnect ( c );
    if ( c.key() == "CHARGE" && this->get_bc(Identifier)==NULL) //  CHARGE card
      Set_Charge ( c );
  }

  MESSAGE<<"Boundary conditions finished.\n"<<std::endl;  RECORD();


  MESSAGE<<"Buliding Geometry Relationship..."; RECORD();
  // let bc set there own data
  for ( unsigned int i=0; i<_bcs.size(); i++ )
    _bcs[i]->prepare_for_use();
  MESSAGE<<"ok.\n"<<std::endl;  RECORD();

  return 0;

}



unsigned int BoundaryConditionCollector::get_bc_from_card(const Parser::Card &c, std::string & Identifier)
{
  if(c.is_parameter_exist("id"))
    Identifier = c.get_string ( "id", "" );
  else if(c.is_parameter_exist("region.first") && c.is_parameter_exist("region.second"))
  {
    std::vector<std::string> ids = _find_interface_by_2_regions(c.get_string ( "region.first", "" ), c.get_string ( "region.second", "" ));
    if(!ids.size())
    {
      MESSAGE<<"ERROR at " <<c.get_fileline() << " Boundary: There is no boundary between region " << c.get_string ( "region.first", "" )
      << " and region " << c.get_string ( "region.second", "" ) <<"."<<std::endl;
      RECORD();
      genius_error();
    }
    if( ids.size() > 1 )
    {
      MESSAGE<<"ERROR at " <<c.get_fileline() << " Boundary: There are more than one boundaries between region " << c.get_string ( "region.first", "" )
      << " and region " << c.get_string ( "region.second", "" ) << ", please specify one by one with its ID." << std::endl;
      RECORD();
      genius_error();
    }
    Identifier = ids[0];
  }

  if( Identifier.empty() )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline() << " Boundary: No boundary specified." <<std::endl; RECORD();
    genius_error();
  }

  // get boundary id by the label
  short int bd_id = _mesh.boundary_info->get_id_by_label ( Identifier );

  if ( bd_id == BoundaryInfo::invalid_id )
  {
    MESSAGE<<"ERROR at " <<c.get_fileline() << " Boundary: ID "<< Identifier <<" can't be found in mesh boundaries." <<std::endl; RECORD();
    genius_error();
  }

  return get_bc_index_by_bd_id ( bd_id );
}


int BoundaryConditionCollector::Set_BC_NeumannBoundary ( const Parser::Card &c )
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new NeumannBC ( _system, Identifier );
  _bcs[bc_index]->T_external()    = c.get_real ( "ext.temp", _system.T_external() /K ) *K;
  _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",0.0 ) *J/s/pow ( cm,2 ) /K;
  _bcs[bc_index]->z_width()       = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
  _bcs[bc_index]->reflection()    = c.get_bool ( "reflection", false );

  return 0;
}



int BoundaryConditionCollector::Set_BC_OhmicContact ( const Parser::Card &c )
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );
  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new OhmicContactBC ( _system, Identifier, false );
  _bcs[bc_index]->T_external() = c.get_real ( "ext.temp", _system.T_external() /K ) *K;
  _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",1e3 ) *J/s/pow ( cm,2 ) /K;
  _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
  _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
  _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
  _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
  _bcs[bc_index]->reflection()    = c.get_bool ( "reflection", true );
  if ( c.is_parameter_exist ( "potential" ) )
    _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;

  if ( c.is_parameter_exist ( "electrode_id" ) )// it is an region interface, not boundary
  {
    _bcs[bc_index]->electrode_label() = c.get_string("electrode_id", "");
    _bcs[bc_index]->set_boundary_type(INTERFACE);
  }
  return 0;
}



int BoundaryConditionCollector::Set_BC_IF_Metal_Ohmic ( const Parser::Card &c )
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );
  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new IF_Metal_OhmicBC ( _system, Identifier );
  if(c.is_parameter_exist("elec.recomb.velocity"))
  {
    double eRecombVelocity = c.get_real ( "elec.recomb.velocity", 2.573e6 )*cm/s;
    if( eRecombVelocity > 1e10*cm/s ) eRecombVelocity = std::numeric_limits<PetscScalar>::infinity();
    _bcs[bc_index]->eRecombVelocity() = eRecombVelocity;
  }
  if(c.is_parameter_exist("hole.recomb.velocity"))
  {
    double hRecombVelocity = c.get_real ( "hole.recomb.velocity", 1.93e6 )*cm/s;
    if( hRecombVelocity > 1e10*cm/s ) hRecombVelocity = std::numeric_limits<PetscScalar>::infinity();
    _bcs[bc_index]->hRecombVelocity() = hRecombVelocity;
  }
  return 0;
}



int BoundaryConditionCollector::Set_BC_SchottkyContact ( const Parser::Card &c )
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new SchottkyContactBC ( _system, Identifier, false );
  _bcs[bc_index]->T_external() = c.get_real ( "ext.temp", _system.T_external() /K ) *K;
  _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",1e3 ) *J/s/pow ( cm,2 ) /K;
  _bcs[bc_index]->Work_Function() = c.get_real ( "workfunction",4.7 ) *V;
  _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
  _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
  _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
  _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
  _bcs[bc_index]->reflection()    = c.get_bool ( "reflection", true );
  if ( c.is_parameter_exist ( "potential" ) )
    _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;

  if ( c.is_parameter_exist ( "electrode_id" ) )// it is an region interface, not boundary
  {
    _bcs[bc_index]->electrode_label() = c.get_string("electrode_id", "");
    _bcs[bc_index]->set_boundary_type(INTERFACE);
  }

  return 0;
}


int BoundaryConditionCollector::Set_BC_IF_Metal_Schottky(const Parser::Card &c)
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new IF_Metal_SchottkyBC ( _system, Identifier );
  return 0;
}


int BoundaryConditionCollector::Set_BC_GateContact ( const Parser::Card &c )
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new GateContactBC ( _system, Identifier, false );
  _bcs[bc_index]->T_external() = c.get_real ( "ext.temp",_system.T_external() /K ) *K;
  _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",1e3 ) *J/s/pow ( cm,2 ) /K;
  _bcs[bc_index]->Work_Function() = c.get_real ( "workfunction",4.17 ) *V;
  _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
  _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
  _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
  _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
  if ( c.is_parameter_exist ( "potential" ) )
    _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;

  if ( c.is_parameter_exist ( "electrode_id" ) )// it is an region interface, not boundary
  {
    _bcs[bc_index]->electrode_label() = c.get_string("electrode_id", "");
    _bcs[bc_index]->set_boundary_type(INTERFACE);
  }

  return 0;
}



int BoundaryConditionCollector::Set_BC_SimpleGateContact ( const Parser::Card &c )
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new SimpleGateContactBC ( _system, Identifier );
  _bcs[bc_index]->T_external() = c.get_real ( "ext.temp",_system.T_external() /K ) *K;
  _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer", 0 ) *J/s/pow ( cm,2 ) /K;
  _bcs[bc_index]->Work_Function() = c.get_real ( "workfunction", 4.17 ) *V;
  _bcs[bc_index]->Thickness() = c.get_real ( "thickness", 1e-2 ) *um;
  _bcs[bc_index]->eps() = c.get_real ( "eps",3.9 ) *eps0;
  _bcs[bc_index]->Qf() = c.get_real ( "qf",1e10 ) /pow ( cm,2 );
  _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
  _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
  _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
  _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
  if ( c.is_parameter_exist ( "potential" ) )
    _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;

  return 0;
}


int BoundaryConditionCollector::Set_BC_Solderpad ( const Parser::Card &c )
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new SolderPadBC ( _system, Identifier );
  _bcs[bc_index]->T_external() = c.get_real ( "ext.temp",_system.T_external() /K ) *K;
  _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer", 0 ) *J/s/pow ( cm,2 ) /K;
  //_bcs[bc_index]->Work_Function() = c.get_real ( "workfunction", 4.17 ) *V;
  _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
  _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
  _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
  _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
  if ( c.is_parameter_exist ( "potential" ) )
    _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;

  return 0;
}


int BoundaryConditionCollector::Set_BC_ChargedContact ( const Parser::Card &c )
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new ChargedContactBC ( _system, Identifier );
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

  return 0;
}


int BoundaryConditionCollector::Set_BC_InsulatorInterface ( const Parser::Card &c )
{
  std::string Identifier;
  unsigned int bc_index = get_bc_from_card ( c, Identifier );

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new InsulatorSemiconductorInterfaceBC ( _system, Identifier );
  _bcs[bc_index]->Qf() = c.get_real ( "qf",1e10 ) /pow ( cm,2 );
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

  return 0;
}









int BoundaryConditionCollector::Set_Electrode_OhmicContact ( const Parser::Card &c )
{
  // get the region label from user input
  std::string region_label = c.get_string ( "id", "" );
  if ( _system.region ( region_label ) == NULL )
  {
    std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" can't be found in mesh structure." <<std::endl;
    genius_error();
  }

  unsigned int subdomain = _mesh.subdomain_id_by_label ( region_label );

  // the region should be electrode with material "Elec"
  std::string region_material = _mesh.subdomain_material(subdomain);
  if(!_resistive_metal_mode)
  {
    if ( !Material::IsConductor(region_material) && !Material::IsResistance(region_material))
    {
      std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" should be conductor or resistance material." <<std::endl;
      genius_error();
    }
  }
  else
  {
    if ( !Material::IsConductor(region_material) )
    {
      std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" should be conductor with material \"Elec\"." <<std::endl;
      genius_error();
    }
  }

  std::set<short int> bd_ids;
  _mesh.boundary_info->get_boundary_ids_by_subdomain ( subdomain, bd_ids );

  // search in the boundary label, which has the format as region1_to_region2,
  // test if region_label matches region1 or region2.

  std::set<short int>::iterator it = bd_ids.begin();

  for ( ; it!=bd_ids.end(); ++it )
  {
    unsigned int sub_id1,sub_id2;

    _mesh.boundary_info->get_subdomains_bd_on ( *it, sub_id1, sub_id2 );

    genius_assert((sub_id1==subdomain) || (sub_id2==subdomain));

    int bc_index = get_bc_index_by_bd_id ( *it );
    std::string label = _mesh.boundary_info->get_label_by_id ( *it );

    // delete previous bc if exist
    if ( _bcs[bc_index]!=NULL ) delete _bcs[bc_index];

    // on the region outer boundary, should be set to NeumannBC
    // at the same time, the NeumannBC here is considered as heat sink
    // with a very high heat transfer rate
    if ( sub_id2==invalid_uint )
    {
      _bcs[bc_index] = new NeumannBC ( _system, label );
      _bcs[bc_index]->Heat_Transfer() = 1e3*J/s/pow ( cm,2 ) /K;
      _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
      continue;
    }

    std::string material1 = _mesh.subdomain_material ( sub_id1 );
    std::string material2 = _mesh.subdomain_material ( sub_id2 );

    switch ( determine_bc_by_subdomain ( material1, material2, _resistive_metal_mode ) )
    {
        case IF_Electrode_Semiconductor  :
        {
          _bcs[bc_index] = new OhmicContactBC ( _system, label, true );
          _bcs[bc_index]->T_external() = c.get_real ( "ext.temp",_system.T_external() /K ) *K;
          _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",1e3 ) *J/s/pow ( cm,2 ) /K;
          _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
          _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
          _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
          _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
          if ( c.is_parameter_exist ( "potential" ) )
            _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;
          _bcs[bc_index]->electrode_label()  = region_label;
          break;
        }
        case IF_Electrode_Insulator      :  _bcs[bc_index] = new ElectrodeInsulatorInterfaceBC ( _system, label );          break;
        case IF_Electrode_Electrode      :  genius_error(); break;
        case IF_Electrode_Vacuum         :  _bcs[bc_index] = new NeumannBC ( _system, label );                              break;
        default: break;
    }
    // override electrode z.width
    _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

  }

  return 0;
}



int BoundaryConditionCollector::Set_Electrode_SchottkyContact ( const Parser::Card &c )
{
  // get the region label from user input
  std::string region_label = c.get_string ( "id","" );
  if ( _system.region ( region_label ) == NULL )
  {
    std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" can't be found in mesh structure." <<std::endl;
    genius_error();
  }

  unsigned int subdomain = _mesh.subdomain_id_by_label ( region_label );

  // the region should be electrode with material "Elec"
  std::string region_material = _mesh.subdomain_material(subdomain);
  if(!_resistive_metal_mode)
  {
    if ( !Material::IsConductor(region_material) && !Material::IsResistance(region_material))
    {
      std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" should be conductor or resistance material." <<std::endl;
      genius_error();
    }
  }
  else
  {
    if ( !Material::IsConductor(region_material) )
    {
      std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" should be conductor with material \"Elec\"." <<std::endl;
      genius_error();
    }
  }

  std::set<short int> bd_ids;
  _mesh.boundary_info->get_boundary_ids_by_subdomain ( subdomain, bd_ids );

  // search in the boundary label, which has the format as region1_to_region2,
  // test if region_label matches region1 or region2.

  std::set<short int>::iterator it = bd_ids.begin();

  for ( ; it!=bd_ids.end(); ++it )
  {
    unsigned int sub_id1,sub_id2;

    _mesh.boundary_info->get_subdomains_bd_on ( *it, sub_id1, sub_id2 );

    genius_assert((sub_id1==subdomain) || (sub_id2==subdomain));

    int bc_index = get_bc_index_by_bd_id ( *it );
    std::string label = _mesh.boundary_info->get_label_by_id ( *it );

    // delete previous bc if exist
    if ( _bcs[bc_index]!=NULL ) delete _bcs[bc_index];

    // on the region outer boundary, should be set to NeumannBC
    // at the same time, the NeumannBC here is considered as heat sink
    // with a very high heat transfer rate
    if ( sub_id2==invalid_uint )
    {
      _bcs[bc_index] = new NeumannBC ( _system, label );
      _bcs[bc_index]->Heat_Transfer() = 1e3*J/s/pow ( cm,2 ) /K;
      _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
      continue;
    }

    std::string material1 = _mesh.subdomain_material ( sub_id1 );
    std::string material2 = _mesh.subdomain_material ( sub_id2 );

    switch ( determine_bc_by_subdomain ( material1, material2, _resistive_metal_mode ) )
    {
        case IF_Electrode_Semiconductor  :
        {
          _bcs[bc_index] = new SchottkyContactBC ( _system, label, true );
          _bcs[bc_index]->T_external() = c.get_real ( "ext.temp",_system.T_external() /K ) *K;
          _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",1e3 ) *J/s/pow ( cm,2 ) /K;
          _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
          _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
          _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
          _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
          if ( c.is_parameter_exist ( "potential" ) )
            _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;
          _bcs[bc_index]->Work_Function() = c.get_real ( "workfunction",4.7 ) *V;
          _bcs[bc_index]->electrode_label()  = region_label;
          break;
        }
        case IF_Electrode_Insulator      :  _bcs[bc_index] = new ElectrodeInsulatorInterfaceBC ( _system, label );          break;
        case IF_Electrode_Electrode      :  genius_error(); break;
        case IF_Electrode_Vacuum         :  _bcs[bc_index] = new NeumannBC ( _system, label );                              break;
        default: break;
    }

    // override electrode z.width
    _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

  }

  return 0;
}



int BoundaryConditionCollector::Set_Electrode_GateContact ( const Parser::Card &c )
{
  // get the region label from user input
  std::string region_label = c.get_string ( "id","" );
  if ( _system.region ( region_label ) == NULL )
  {
    std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" can't be found in mesh structure." <<std::endl;
    genius_error();
  }

  unsigned int subdomain = _mesh.subdomain_id_by_label ( region_label );

  // the region should be electrode with material "Elec"
  std::string region_material = _mesh.subdomain_material(subdomain);
  if(!_resistive_metal_mode)
  {
    if ( !Material::IsConductor(region_material) && !Material::IsResistance(region_material))
    {
      std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" should be conductor or resistance material." <<std::endl;
      genius_error();
    }
  }
  else
  {
    if ( !Material::IsConductor(region_material) )
    {
      std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" should be conductor with material \"Elec\"." <<std::endl;
      genius_error();
    }
  }

  std::set<short int> bd_ids;
  _mesh.boundary_info->get_boundary_ids_by_subdomain ( subdomain, bd_ids );

  // search in the boundary label, which has the format as region1_to_region2,
  // test if region_label matches region1 or region2.

  std::set<short int>::iterator it = bd_ids.begin();

  for ( ; it!=bd_ids.end(); ++it )
  {
    unsigned int sub_id1,sub_id2;

    _mesh.boundary_info->get_subdomains_bd_on ( *it, sub_id1, sub_id2 );

    genius_assert((sub_id1==subdomain) || (sub_id2==subdomain));

    int bc_index = get_bc_index_by_bd_id ( *it );
    std::string label = _mesh.boundary_info->get_label_by_id ( *it );

    // delete previous bc if exist
    if ( _bcs[bc_index]!=NULL ) delete _bcs[bc_index];

    // on the region outer boundary, should be set to NeumannBC
    // at the same time, the NeumannBC here is considered as heat sink
    // with a very high heat transfer rate
    if ( sub_id2==invalid_uint )
    {
      _bcs[bc_index] = new NeumannBC ( _system, label );
      _bcs[bc_index]->Heat_Transfer() = 1e3*J/s/pow ( cm,2 ) /K;
      _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
      continue;
    }

    std::string material1 = _mesh.subdomain_material ( sub_id1 );
    std::string material2 = _mesh.subdomain_material ( sub_id2 );

    switch ( determine_bc_by_subdomain ( material1, material2, _resistive_metal_mode ) )
    {
        case IF_Electrode_Semiconductor  :  genius_error(); break; // we should never reach here
        case IF_Electrode_Insulator      :
        {
          _bcs[bc_index] = new GateContactBC ( _system, label, true );
          _bcs[bc_index]->T_external() = c.get_real ( "ext.temp",_system.T_external() /K ) *K;
          _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",1e3 ) *J/s/pow ( cm,2 ) /K;
          _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
          _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
          _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
          _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
          if ( c.is_parameter_exist ( "potential" ) )
            _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;
          _bcs[bc_index]->Work_Function() = c.get_real ( "workfunction",4.17 ) *V;
          _bcs[bc_index]->electrode_label()  = region_label;
          break;
        }
        case IF_Electrode_Electrode      :  genius_error();                                  break;
        case IF_Electrode_Vacuum         :  _bcs[bc_index] = new NeumannBC ( _system, label );  break;
        default: break;
    }

    // override electrode z.width
    _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

  }

  return 0;
}


int BoundaryConditionCollector::Set_SetFloatMetal ( const Parser::Card &c )
{
  // get the region label from user input
  std::string region_label = c.get_string ( "id","" );
  if ( _system.region ( region_label ) == NULL )
  {
    std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" can't be found in mesh structure." <<std::endl;
    genius_error();
  }

  unsigned int subdomain = _mesh.subdomain_id_by_label ( region_label );

  // the region should be electrode with material "Elec"
  std::string region_material = _mesh.subdomain_material(subdomain);
  if ( !Material::IsConductor(region_material) && !Material::IsResistance(region_material) )
  {
    std::cerr<<"ERROR at " <<c.get_fileline() << " CONTACT: region "<< region_label <<" should be conductor or resistance metal." <<std::endl;
    genius_error();
  }

  // all the pisces of charge boundary
  std::set<BoundaryCondition *> charge_boundaries;

  std::set<short int> bd_ids;
  _mesh.boundary_info->get_boundary_ids_by_subdomain ( subdomain, bd_ids );

  // search in the boundary label, which has the form as region_to_region, and if region_label
  // matches at head or back.

  std::set<short int>::iterator it = bd_ids.begin();
  for ( ; it!=bd_ids.end(); ++it )
  {
    unsigned int sub_id1,sub_id2;

    _mesh.boundary_info->get_subdomains_bd_on ( *it, sub_id1, sub_id2 );

    genius_assert((sub_id1==subdomain) || (sub_id2==subdomain));

    int bc_index = get_bc_index_by_bd_id ( *it );
    std::string label = _mesh.boundary_info->get_label_by_id ( *it );

    // delete previous bc if exist
    if ( _bcs[bc_index]!=NULL ) delete _bcs[bc_index];

    // on the region outer boundary, should not exist
    genius_assert( sub_id2 != invalid_uint );

    std::string material1 = _mesh.subdomain_material ( sub_id1 );
    std::string material2 = _mesh.subdomain_material ( sub_id2 );

    switch ( determine_bc_by_subdomain ( material1, material2, _resistive_metal_mode ) )
    {
        case IF_Electrode_Insulator      :
        case IF_Insulator_Metal     :
        {
          // float metal to insulator interface
          _bcs[bc_index] = new ChargedContactBC ( _system, label );
          _bcs[bc_index]->electrode_label()  = region_label;
          // override electrode z.width
          _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;
          charge_boundaries.insert(_bcs[bc_index]);
          break;
        }
        case IF_Electrode_Semiconductor  :  genius_error(); break; // we should never reach here
        case IF_Electrode_Electrode      :  genius_error(); break;
        case IF_Electrode_Vacuum         :  genius_error(); break;
        default: break;
    }
  }

  {
    // build a new bc for ChargeIntegral
    BoundaryCondition * charge_integral = new ChargeIntegralBC ( _system, region_label+"_ChargeIntegral" );
    charge_integral->set_inter_connect ( charge_boundaries );
    charge_integral->Qf() = c.get_real ( "qf", 0.0 ) *C;

    this->add_bc ( charge_integral );

    // setup each charge bc that should be integral
    std::set<BoundaryCondition * >::iterator it=charge_boundaries.begin();
    for ( ; it != charge_boundaries.end(); ++it )
    {
      ( *it )->set_inter_connect ( charge_boundaries );
      ( *it )->set_inter_connect_hub ( charge_integral );
    }
  }

  return 0;
}


int BoundaryConditionCollector::Set_InterConnect ( const Parser::Card &c )
{
  std::string id = c.get_string ( "id", "" );

  if ( id=="" )
  {
    std::cerr<<"ERROR at " <<c.get_fileline() << " InterConnect: should have an unique ID." <<std::endl;
    genius_error();
  }

  // find all the electrodes inter connect to each other
  std::set<BoundaryCondition *> inter_connect_electrodes;
  for ( unsigned int idx=0; idx<c.parameter_size(); ++idx )
  {
    Parser::Parameter p = c.get_parameter ( idx );
    if ( p.name() == "connectto" )
    {
      std::string electrode_label = p.get_string();
      BoundaryCondition * electrode = this->get_bc ( electrode_label );
      if ( electrode && electrode->is_electrode() )
        inter_connect_electrodes.insert ( electrode );
      else
      {
        std::cerr<<"ERROR at " <<c.get_fileline() << " InterConnect: "<< electrode_label <<" is not an electrode." <<std::endl;
        genius_error();
      }
    }
  }

  // build a new bc for inter connector
  BoundaryCondition * inter_connector = new ElectrodeInterConnectBC ( _system, id );
  inter_connector->set_inter_connect ( inter_connect_electrodes );
  inter_connector->build_ext_circuit ( new ExternalCircuit );
  inter_connector->ext_circuit()->R() = std::max ( c.get_real ( "res", 0.0 ) *V/A, 1e-3*V/A );
  inter_connector->ext_circuit()->C() = c.get_real ( "cap", 0.0 ) *C/V;
  inter_connector->ext_circuit()->L() = c.get_real ( "ind", 0.0 ) *V*s/A;
  if ( c.is_parameter_exist ( "potential" ) )
    inter_connector->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;
  // the default state of inter_connector is float
  if( c.get_bool ( "float", true ) )
    inter_connector->ext_circuit()->set_float();
  else
    inter_connector->ext_circuit()->set_voltage_driven();
  this->add_bc ( inter_connector );


  // setup each electrode belongs to inter connect
  std::set<BoundaryCondition * >::iterator it=inter_connect_electrodes.begin();
  for ( ; it != inter_connect_electrodes.end(); ++it )
  {
    ( *it )->set_inter_connect ( inter_connect_electrodes );
    ( *it )->set_inter_connect_hub ( inter_connector );
    ( *it )->ext_circuit()->R() = std::max ( ( *it )->ext_circuit()->R(), 1e-3*V/A );
  }

  return 0;
}


int BoundaryConditionCollector::Set_Charge(const Parser::Card &c)
{
  std::string id = c.get_string ( "id", "" );

  if ( id=="" )
  {
    std::cerr<<"ERROR at " <<c.get_fileline() << " Charge: should have an unique ID." <<std::endl;
    genius_error();
  }

  // find all the charge boundaries
  std::set<BoundaryCondition *> charge_boundaries;
  for ( unsigned int idx=0; idx<c.parameter_size(); ++idx )
  {
    Parser::Parameter p = c.get_parameter ( idx );
    if ( p.name() == "chargeboundary" )
    {
      std::string charge_boundary_label = p.get_string();
      BoundaryCondition * charge_boundary = this->get_bc ( charge_boundary_label );
      if ( charge_boundary && charge_boundary->bc_type() == ChargedContact )
        charge_boundaries.insert ( charge_boundary );
      else
      {
        std::cerr<<"ERROR at " <<c.get_fileline() << " Charge: "<< charge_boundary_label <<" is not a charge boundary." <<std::endl;
        genius_error();
      }
    }
  }


  // build a new bc for charge integral
  BoundaryCondition * charge_integral = new ChargeIntegralBC ( _system, id );
  charge_integral->set_inter_connect ( charge_boundaries );
  charge_integral->Qf() = c.get_real ( "charge", 0.0, "qf" )*C;
  this->add_bc ( charge_integral );

  // setup each charge boundary belongs to charge integral
  std::set<BoundaryCondition * >::iterator it=charge_boundaries.begin();
  for ( ; it != charge_boundaries.end(); ++it )
  {
    ( *it )->set_inter_connect ( charge_boundaries );
    ( *it )->set_inter_connect_hub ( charge_integral );
  }

  return 0;
}


int BoundaryConditionCollector::Set_BC_AbsorbingBoundary ( const Parser::Card &c )
{
  // get the boundary label from user input
  std::string Identifier = c.get_string ( "id","" );

  // get boundary id by the label
  short int bd_id = _mesh.boundary_info->get_id_by_label ( Identifier );

  if ( bd_id == BoundaryInfo::invalid_id )
  {
    std::cerr<<"ERROR at " <<c.get_fileline() << " Boundary: ID "<< Identifier <<" can't be found in mesh boundaries." <<std::endl;
    genius_error();
  }

  int bc_index = get_bc_index_by_bd_id ( bd_id );

  _bcs[bc_index] = new AbsorbingBC ( _system, Identifier );

  return 0;
}


std::vector<std::string> BoundaryConditionCollector::_find_interface_by_2_regions(const std::string &r1, const std::string &r2) const
{
  unsigned int sub1 =  _mesh.subdomain_id_by_label(r1);
  unsigned int sub2 =  _mesh.subdomain_id_by_label(r2);

  std::vector<std::string> bd_labels;
  std::set<short int> bd_ids = _mesh.boundary_info->get_ids_on_region_interface(sub1, sub2);
  std::set<short int>::iterator it = bd_ids.begin();
  for( ; it != bd_ids.end(); ++it)
    bd_labels.push_back(_mesh.boundary_info->get_label_by_id(*it));
  return bd_labels;
}


void BoundaryConditionCollector::export_boundary_condition ( const std::string& filename ) const
{
  MESSAGE<<"Write Boundary Condition to file "<< filename << "...\n" << std::endl; RECORD();

  if ( Genius::processor_id() == 0 )
  {
    std::ofstream _out ( filename.c_str(), std::ofstream::trunc );

    for ( unsigned int i=0; i<_bcs.size(); i++ )
      _out<<_bcs[i]->boundary_condition_in_string();

    _out.close();
  }
}

mxml_node_t* BoundaryConditionCollector::get_dom_terminal_info() const
{
  mxml_node_t *eTerm = mxmlNewElement(NULL, "terminal-info");

  // record electrode IV information
  // search for all the bcs
  for(unsigned int n=0; n<n_bcs(); n++)
  {
    const BoundaryCondition * bc = get_bc(n);
    // skip bc which is not electrode
    if( !bc->is_electrode() ) continue;
    //
    std::string bc_label = bc->label();
    if(!bc->electrode_label().empty())
      bc_label = bc->electrode_label();
    //
    double volt, vf, curr;
    volt = bc->ext_circuit()->Vapp()/PhysicalUnit::V;
    vf   = bc->ext_circuit()->potential()/PhysicalUnit::V;
    curr = bc->ext_circuit()->current()/PhysicalUnit::A;

    mxml_node_t *eContact = mxmlNewElement(eTerm, "contact");

    mxml_node_t *eLabel = mxmlNewElement(eContact, "label");
    mxmlAdd(eLabel, MXML_ADD_AFTER, 0, MXMLQVariant::makeQVString(bc_label));

    mxml_node_t *eVolt = mxmlNewElement(eContact, "voltage");
    mxmlAdd(eVolt, MXML_ADD_AFTER, 0, MXMLQVariant::makeQVFloat(volt));

    mxml_node_t *eVf   = mxmlNewElement(eContact, "fermi-potential");
    mxmlAdd(eVf, MXML_ADD_AFTER, 0, MXMLQVariant::makeQVFloat(vf));

    mxml_node_t *eCurr = mxmlNewElement(eContact, "current");
    mxmlAdd(eCurr, MXML_ADD_AFTER, 0, MXMLQVariant::makeQVFloat(curr));

  }
  return eTerm;
}
