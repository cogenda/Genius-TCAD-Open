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

#include "mesh.h"
#include "fvm_node_info.h"
#include "fvm_node_data.h"
#include "simulation_system.h"
#include "conductor_region.h"
#include "boundary_info.h"
#include "boundary_condition_collector.h"
#include "mxml.h"
#include "MXMLUtil.h"

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
    : _decks ( decks ), _system ( system ), _mesh ( system.mesh() )
{}


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
        case NeumannBoundary             :  { if ( Set_BC_NeumannBoundary ( c ) )    return 1; break;}
        case OhmicContact                :  { if ( Set_BC_OhmicContact ( c ) )       return 1; break;}
        case SchottkyContact             :  { if ( Set_BC_SchottkyContact ( c ) )    return 1; break;}
        case GateContact                 :  { if ( Set_BC_GateContact ( c ) )        return 1; break;}
        case SimpleGateContact           :  { if ( Set_BC_SimpleGateContact ( c ) )  return 1; break;}
        case IF_Insulator_Semiconductor  :  { if ( Set_BC_InsulatorInterface ( c ) ) return 1; break;}
        case AbsorbingBoundary           :  { if ( Set_BC_AbsorbingBoundary ( c ) )  return 1; break;}
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
      std::string description = _mesh.boundary_info->get_description_by_id(boundary_id);
      if(description.empty()) continue;

      Parser::Card c(description);
      if(c.key() != "BOUNDARY") continue;
      BCType bc_type = BC_string_to_enum(c.get_string("type", ""));

      switch(bc_type)
      {
      case NeumannBoundary             :  { if( Set_BC_NeumannBoundary(c) )    return 1; break;}
      case OhmicContact                :  { if( Set_BC_OhmicContact(c) )       return 1; break;}
      case SchottkyContact             :  { if( Set_BC_SchottkyContact(c) )    return 1; break;}
      case GateContact                 :  { if( Set_BC_GateContact(c) )        return 1; break;}
      case SimpleGateContact           :  { if( Set_BC_SimpleGateContact(c) )  return 1; break;}
      case IF_Insulator_Semiconductor  :  { if( Set_BC_InsulatorInterface(c) ) return 1; break;}
      case AbsorbingBoundary           :  { if( Set_BC_AbsorbingBoundary(c) )  return 1; break;}
      default: genius_error();
      }
    }

  // set remaining bcs
  for ( unsigned int i=0; i<_bcs.size(); i++ )
    if ( _bcs[i] == NULL ) // not initialized?
    {
      short int boundary_id = get_bd_id_by_bc_index ( i );
      std::string label = _mesh.boundary_info->get_label_by_id ( boundary_id );

      unsigned int sub_id1, sub_id2;
      _mesh.boundary_info->get_subdomains_bd_on ( boundary_id, sub_id1, sub_id2 );

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


      switch ( determine_bc_by_subdomain ( material1, material2 ) )
      {
        case IF_Electrode_Semiconductor  :
        {
          _bcs[i] = new OhmicContactBC ( _system, label, true );
          _bcs[i]->build_ext_circuit ( new ExternalCircuit );
          // allow one use electrode label to load the bc
          if ( Material::IsConductor ( material1 ) )
            _bcs[i]->electrode_label() = _mesh.subdomain_label_by_id ( sub_id1 );
          if ( Material::IsConductor ( material2 ) )
            _bcs[i]->electrode_label() = _mesh.subdomain_label_by_id ( sub_id2 );
          break;
        }
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
          unsigned int electrode_sbd =  Material::IsConductor ( material1 ) ? sub_id1 : sub_id2;
          const ConductorSimulationRegion * electrode_region = dynamic_cast<ConductorSimulationRegion *>(_system.region(electrode_sbd));

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
            _bcs[i]->Work_Function() = electrode_region->material()->basic->Affinity(_system.T_external());
            _bcs[i]->build_ext_circuit ( new ExternalCircuit );
          }
          break;
        }
        case IF_Electrode_Electrode      :  genius_error();                                                          break;
        case IF_Electrode_Vacuum         :  _bcs[i] = new NeumannBC ( _system, label );                              break;
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

  // reset boundary id to the nodes in boundary_info class with a priority order

  std::map<short int, unsigned int>  order;
  for ( unsigned int i=0; i<_bcs.size(); i++ )
    { order.insert ( std::make_pair ( get_bd_id_by_bc_index ( i ), static_cast<unsigned int> ( _bcs[i]->bc_type() ) ) ); }

  _mesh.boundary_info->build_node_ids_from_priority_order ( order );


  // map Node * to bd_id
  std::map<const Node *, short int> node_bd_map;


  for ( unsigned int i=0; i<_bcs.size(); i++ )
  {
    short int boundary_id = get_bd_id_by_bc_index ( i );

    // save nodes which has bd_id into the corresponding bc class
    std::vector<const Node *> bd_nodes;
    // this function build the sorted boundary node list with their id
    _mesh.boundary_info->nodes_with_boundary_id ( bd_nodes, boundary_id );

    for ( unsigned int n=0; n<bd_nodes.size(); n++ )
    {
      _bcs[i]->add_node ( bd_nodes[n] );

      // set the corresponding region information and FVM_Node to this node
      for ( unsigned int r=0; r<_system.n_regions(); r++ )
      {
        SimulationRegion * region = _system.region ( r );
        FVM_Node * fvm_node;
        // we find a FVM_Node with its root node is bd_nodes[n]
        if ( ( fvm_node=region->region_fvm_node ( bd_nodes[n] ) ) != NULL )
          _bcs[i]->insert ( bd_nodes[n], region, fvm_node );
      }

      // insert node* into map, for later usage
      node_bd_map.insert ( std::pair<const Node *, short int> ( bd_nodes[n], boundary_id ) );
    }


    // save elements which has bd_id into the corresponding bc class
    std::vector<const Elem *> el;
    std::vector<unsigned int> sl;
    _mesh.boundary_info->active_elem_with_boundary_id ( el, sl, boundary_id );
    for ( unsigned int n=0; n<el.size(); n++ )
    {
      _bcs[i]->add_elem ( el[n], sl[n] );
    }

  }


  // set bd_id (boundary index) to region FVM_Node,
  // then we can easily find boundary_condition_index for each FVM_Node
  for ( unsigned int r=0; r<_system.n_regions(); r++ )
  {
    SimulationRegion * region = _system.region ( r );

    std::map<const Node *, short int>::iterator map_it = node_bd_map.begin();
    for ( ; map_it != node_bd_map.end(); ++map_it )
    {
      FVM_Node * fvm_node = region->region_fvm_node ( ( *map_it ).first );

      if ( fvm_node == NULL ) continue;

      fvm_node->set_boundary_id ( ( *map_it ).second );
    }
  }


  // process inter-connect here
  for ( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();
    if ( c.key() == "INTERCONNECT" )  //  InterConnect card
      Set_InterConnect ( c );
  }

  MESSAGE<<"Boundary conditions finished.\n"<<std::endl;  RECORD();

  return 0;

}



int BoundaryConditionCollector::Set_BC_NeumannBoundary ( const Parser::Card &c )
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

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new OhmicContactBC ( _system, Identifier, false );
  _bcs[bc_index]->T_external() = c.get_real ( "ext.temp", _system.T_external() /K ) *K;
  _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",1e3 ) *J/s/pow ( cm,2 ) /K;
  _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
  _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
  _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
  _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
  if ( c.is_parameter_exist ( "potential" ) )
    _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

  if ( c.is_parameter_exist ( "electrode_id" ) )
    _bcs[bc_index]->electrode_label() = c.get_string("electrode_id", "");

  return 0;
}



int BoundaryConditionCollector::Set_BC_SchottkyContact ( const Parser::Card &c )
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

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new SchottkyContactBC ( _system, Identifier, false );
  _bcs[bc_index]->T_external() = c.get_real ( "ext.temp", _system.T_external() /K ) *K;
  _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",1e3 ) *J/s/pow ( cm,2 ) /K;
  _bcs[bc_index]->Work_Function() = c.get_real ( "workfunction",4.7 ) *V;
  _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
  _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
  _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
  _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
  if ( c.is_parameter_exist ( "potential" ) )
    _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

  if ( c.is_parameter_exist ( "electrode_id" ) )
    _bcs[bc_index]->electrode_label() = c.get_string("electrode_id", "");

  return 0;
}



int BoundaryConditionCollector::Set_BC_GateContact ( const Parser::Card &c )
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

  // build the boundary condition parameter for this boundary
  _bcs[bc_index] = new GateContactBC ( _system, Identifier, false );
  _bcs[bc_index]->T_external() = c.get_real ( "ext.temp",_system.T_external() /K ) *K;
  _bcs[bc_index]->Heat_Transfer() = c.get_real ( "heat.transfer",1e3 ) *J/s/pow ( cm,2 ) /K;
  _bcs[bc_index]->Work_Function() = c.get_real ( "workfunction",4.17 ) *V;
  _bcs[bc_index]->build_ext_circuit ( new ExternalCircuit );
  _bcs[bc_index]->ext_circuit()->R() = c.get_real ( "res",0.0 ) *V/A;
  _bcs[bc_index]->ext_circuit()->C() = c.get_real ( "cap",0.0 ) *C/V;
  _bcs[bc_index]->ext_circuit()->L() = c.get_real ( "ind",0.0 ) *V*s/A;
  if ( c.is_parameter_exist ( "potential" ) )
    _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

  if ( c.is_parameter_exist ( "electrode_id" ) )
    _bcs[bc_index]->electrode_label() = c.get_string("electrode_id", "");

  return 0;
}



int BoundaryConditionCollector::Set_BC_SimpleGateContact ( const Parser::Card &c )
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
  if ( c.is_parameter_exist ( "potential" ) )
    _bcs[bc_index]->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;
  _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

  return 0;
}



int BoundaryConditionCollector::Set_BC_InsulatorInterface ( const Parser::Card &c )
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

  // the region should be electrode with material "conductor"
  std::string region_material = _mesh.subdomain_material(subdomain);
  genius_assert( Material::IsConductor(region_material) );

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

    switch ( determine_bc_by_subdomain ( material1, material2 ) )
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

  // the region should be electrode with material "conductor"
  std::string region_material = _mesh.subdomain_material(subdomain);
  genius_assert( Material::IsConductor(region_material) );

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

    switch ( determine_bc_by_subdomain ( material1, material2 ) )
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

  // the region should be electrode with material "conductor"
  std::string region_material = _mesh.subdomain_material(subdomain);
  genius_assert( Material::IsConductor(region_material) );

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

    switch ( determine_bc_by_subdomain ( material1, material2 ) )
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

  // the region should be electrode with material "conductor"
  std::string region_material = _mesh.subdomain_material(subdomain);
  genius_assert( Material::IsConductor(region_material) );

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

    switch ( determine_bc_by_subdomain ( material1, material2 ) )
    {
      case IF_Electrode_Semiconductor  :  genius_error(); break; // we should never reach here
      case IF_Electrode_Insulator      :
      {
        // float metal to insulator interface
        _bcs[bc_index] = new ChargedContactBC ( _system, label );
        _bcs[bc_index]->Qf() = c.get_real ( "qf",0.0 ) *C;
        _bcs[bc_index]->electrode_label()  = region_label;
        break;
      }
      case IF_Electrode_Electrode      :  genius_error(); break;
      case IF_Electrode_Vacuum         :  genius_error(); break;
      default: break;
    }

    // override electrode z.width
    _bcs[bc_index]->z_width() = c.get_real ( "z.width", _bcs[bc_index]->z_width() /um ) *um;

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
  BoundaryCondition * inter_connector = new InterConnectBC ( _system, id );
  inter_connector->set_inter_connect ( inter_connect_electrodes );
  inter_connector->build_ext_circuit ( new ExternalCircuit );
  inter_connector->ext_circuit()->R() = std::max ( c.get_real ( "res", 0.0 ) *V/A, 1e-3*V/A );
  inter_connector->ext_circuit()->C() = c.get_real ( "cap", 0.0 ) *C/V;
  inter_connector->ext_circuit()->L() = c.get_real ( "ind", 0.0 ) *V*s/A;
  if ( c.is_parameter_exist ( "potential" ) )
    inter_connector->ext_circuit()->potential() = c.get_real ( "potential", 0.0 ) *V;
  // the default state of inter_connector is float
  inter_connector->ext_circuit()->set_float();
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
