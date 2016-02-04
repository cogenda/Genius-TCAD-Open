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

//  $Id: boundary_condition.cc,v 1.22 2008/07/09 05:58:16 gdiso Exp $

#include "mesh_base.h"
#include "boundary_info.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "boundary_condition.h"
#include "material_define.h"

using namespace Material;

using PhysicalUnit::J;
using PhysicalUnit::K;
using PhysicalUnit::eps0;
using PhysicalUnit::mu0;
using PhysicalUnit::s;
using PhysicalUnit::um;
using PhysicalUnit::cm;
using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::C;



// init the static dummy member of BoundaryCondition class
unsigned int BoundaryCondition::_solver_index=0;
PetscScalar BoundaryCondition::_dummy_ = 0;





static std::map<std::string, BCType> bc_name_to_bc_type;


void init_bc_name_to_bc_type_map()
{
  // if the bc_name_to_bc_type not build yet
  if ( bc_name_to_bc_type.empty() )
  {
    bc_name_to_bc_type["neumann"                  ]  = NeumannBoundary;
    bc_name_to_bc_type["ohmiccontact"             ]  = OhmicContact;
    bc_name_to_bc_type["metalohmicinterface"      ]  = IF_Metal_Ohmic;
    bc_name_to_bc_type["schottkycontact"          ]  = SchottkyContact;
    bc_name_to_bc_type["metalschottkyinterface"   ]  = IF_Metal_Schottky;
    bc_name_to_bc_type["gatecontact"              ]  = GateContact;
    bc_name_to_bc_type["polygatecontact"          ]  = PolyGateContact;
    bc_name_to_bc_type["simplegatecontact"        ]  = SimpleGateContact;
    bc_name_to_bc_type["solderpad"                ]  = SolderPad;
    bc_name_to_bc_type["insulatorinterface"       ]  = IF_Insulator_Semiconductor;
    bc_name_to_bc_type["heterojunction"           ]  = HeteroInterface;
    bc_name_to_bc_type["chargedcontact"           ]  = ChargedContact;
    bc_name_to_bc_type["absorbingboundary"        ]  = AbsorbingBoundary;
    bc_name_to_bc_type["sourceboundary"           ]  = SourceBoundary;
    bc_name_to_bc_type["floatmetal"               ]  = ChargedContact;
    bc_name_to_bc_type["fixedpotential"           ]  = FixedPotential;
    bc_name_to_bc_type["chargeemit"               ]  = ChargeEmit;
  }
}



BCType BC_string_to_enum(const std::string & bc_type)
{

  init_bc_name_to_bc_type_map();

  if ( bc_name_to_bc_type.find(bc_type)!= bc_name_to_bc_type.end() )
    return bc_name_to_bc_type[bc_type];
  else
    return INVALID_BC_TYPE;
}


static std::map<BCType, std::string> bc_type_to_bc_name;

void init_bc_type_to_bc_name_map()
{
  // if the bc_type_to_bc_name not build yet
  if ( bc_type_to_bc_name.empty() )
  {
    bc_type_to_bc_name[NeumannBoundary                 ]  = std::string("NeumannBoundary");
    bc_type_to_bc_name[IF_Semiconductor_Vacuum         ]  = std::string("IF_Semiconductor_Vacuum");
    bc_type_to_bc_name[IF_Insulator_Vacuum             ]  = std::string("IF_Insulator_Vacuum");
    bc_type_to_bc_name[IF_Electrode_Vacuum             ]  = std::string("IF_Electrode_Vacuum");
    bc_type_to_bc_name[IF_Metal_Vacuum                 ]  = std::string("IF_Metal_Vacuum");
    bc_type_to_bc_name[IF_PML_Scatter                  ]  = std::string("IF_PML_Scatter");
    bc_type_to_bc_name[IF_PML_PML                      ]  = std::string("IF_PML_PML");
    bc_type_to_bc_name[IF_Electrode_Insulator          ]  = std::string("IF_Electrode_Insulator");
    bc_type_to_bc_name[IF_Insulator_Semiconductor      ]  = std::string("IF_Insulator_Semiconductor");
    bc_type_to_bc_name[IF_Insulator_Insulator          ]  = std::string("IF_Insulator_Insulator");
    bc_type_to_bc_name[IF_Electrode_Electrode          ]  = std::string("IF_Electrode_Electrode");
    bc_type_to_bc_name[IF_Electrode_Metal              ]  = std::string("IF_Electrode_Metal");
    bc_type_to_bc_name[IF_Insulator_Metal              ]  = std::string("IF_Insulator_Metal");
    bc_type_to_bc_name[IF_Metal_Metal                  ]  = std::string("IF_Metal_Metal");
    bc_type_to_bc_name[HomoInterface                   ]  = std::string("HomoInterface");
    bc_type_to_bc_name[HeteroInterface                 ]  = std::string("HeteroInterface");
    bc_type_to_bc_name[IF_Electrode_Semiconductor      ]  = std::string("IF_Electrode_Semiconductor");
    bc_type_to_bc_name[IF_Metal_Semiconductor          ]  = std::string("IF_Metal_Semiconductor");
    bc_type_to_bc_name[OhmicContact                    ]  = std::string("OhmicContact");
    bc_type_to_bc_name[IF_Metal_Ohmic                  ]  = std::string("IF_Metal_Ohmic");
    bc_type_to_bc_name[SchottkyContact                 ]  = std::string("SchottkyContact");
    bc_type_to_bc_name[IF_Metal_Schottky               ]  = std::string("IF_Metal_Schottky");
    bc_type_to_bc_name[SimpleGateContact               ]  = std::string("SimpleGateContact");
    bc_type_to_bc_name[GateContact                     ]  = std::string("GateContact");
    bc_type_to_bc_name[ChargedContact                  ]  = std::string("ChargedContact");
    bc_type_to_bc_name[ChargeIntegral                  ]  = std::string("ChargeIntegral");
    bc_type_to_bc_name[ChargeEmit                      ]  = std::string("ChargeEmit");
    bc_type_to_bc_name[SolderPad                       ]  = std::string("SolderPad");
    bc_type_to_bc_name[InterConnect                    ]  = std::string("InterConnect");
    bc_type_to_bc_name[AbsorbingBoundary               ]  = std::string("AbsorbingBoundary");
    bc_type_to_bc_name[SourceBoundary                  ]  = std::string("SourceBoundary");
    bc_type_to_bc_name[INVALID_BC_TYPE                 ]  = std::string("INVALID_BC_TYPE");
  }
}

std::string BC_enum_to_string(BCType bc_type)
{
  init_bc_type_to_bc_name_map();

  if ( bc_type_to_bc_name.find(bc_type)!= bc_type_to_bc_name.end() )
    return bc_type_to_bc_name[bc_type];
  else
    return std::string("INVALID_BC_TYPE");
}



BCType determine_bc_by_subdomain(const std::string & mat1, const std::string & mat2, bool resistive_metal_mode)
{
  // two materials are both semiconductor, but with different type. it is  Heterojunction
  if( IsSemiconductor(mat1) && IsSemiconductor(mat2) && FormatMaterialString(mat1) != FormatMaterialString(mat2) )
    return HeteroInterface;

  // two materials are both semiconductor and with same type, just a Homojunction
  if( IsSemiconductor(mat1) && IsSemiconductor(mat2) && FormatMaterialString(mat1) == FormatMaterialString(mat2) )
    return HomoInterface;

  //one material is semiconductor
  if( IsSemiconductor(mat1) || IsSemiconductor(mat2) )
  {
    // another material is insulator
    if( IsInsulator(mat1) || IsInsulator(mat2) )
      return IF_Insulator_Semiconductor;

    // another material is vacuum
    if( IsVacuum(mat1) || IsVacuum(mat2) )
      return IF_Semiconductor_Vacuum;

    // another material is conductor (electrode)
    if( IsConductor(mat1) || IsConductor(mat2) )
      return IF_Electrode_Semiconductor;

    // another material is resistance metal
    if( IsResistance(mat1) || IsResistance(mat2) )
    {
      if(resistive_metal_mode)
        return IF_Metal_Semiconductor;
      return IF_Electrode_Semiconductor;
    }
  }

  // two materials are both Insulator
  if( IsInsulator(mat1) && IsInsulator(mat2) )
    return IF_Insulator_Insulator;

  //one material is Insulator
  if( IsInsulator(mat1) || IsInsulator(mat2) )
  {
    // another material is vacuum
    if( IsVacuum(mat1) || IsVacuum(mat2) )
      return IF_Insulator_Vacuum;

    // another material is conductor (electrode)
    if( IsConductor(mat1) || IsConductor(mat2) )
      return IF_Electrode_Insulator;

    // another material is resistance metal
    if( IsResistance(mat1) || IsResistance(mat2) )
    {
      if(resistive_metal_mode)
        return IF_Insulator_Metal;
      return IF_Electrode_Insulator;
    }
  }

  // two materials are both Conductor
  if(IsConductor(mat1) && IsConductor(mat2) )
    return IF_Electrode_Electrode;

  // one material is Conductor
  if(IsConductor(mat1) || IsConductor(mat2))
  {
    if(IsVacuum(mat1) || IsVacuum(mat2))
      return IF_Electrode_Vacuum;

    // another material is resistance metal
    if( IsResistance(mat1) || IsResistance(mat2) )
    {
      if(resistive_metal_mode)
        return IF_Electrode_Metal;
      return IF_Electrode_Electrode;
    }
  }

  // two materials are both resistance metal
  if(IsResistance(mat1) && IsResistance(mat2) )
  {
    if(resistive_metal_mode)
      return IF_Metal_Metal;
    return IF_Electrode_Electrode;
  }

  // one material is resistance metal
  if(IsResistance(mat1) || IsResistance(mat2))
  {
    if(IsVacuum(mat1) || IsVacuum(mat2))
    {
      if(resistive_metal_mode)
        return IF_Metal_Vacuum;
      return IF_Electrode_Vacuum;
    }
  }


  // two materials are both PML
  if(IsPML(mat1) && IsPML(mat2))
    return IF_PML_PML;

  // one material is PML
  if(IsPML(mat1) || IsPML(mat2))
  {
      return IF_PML_Scatter;
  }

  // we should never reach here
  MESSAGE<<"ERROR: The interface type between two materials "<<mat1<<" and "<<mat2<<" can't be determined."<<std::endl;
  genius_error();

  return INVALID_BC_TYPE;
}


bool consistent_bc_by_subdomain(const std::string &mat1, const std::string &mat2, BCType type, bool pisces_compatible_mode)
{
  std::string material1 = mat1;
  std::string material2 = mat2;

  if(Material::material_type(material1) > Material::material_type(material2))
    std::swap(material1, material2);

  switch(type)
  {
    case OhmicContact :
    {
      if(IsSemiconductor(material1) && material2.empty() ) return true;
      if(IsSemiconductor(material1) && IsInsulator(material2) ) return true;
      if(IsSemiconductor(material1) && IsConductor(material2) ) return true;
      if(IsSemiconductor(material1) && IsResistance(material2) ) return true;
      return false;
    }
    case IF_Metal_Ohmic:
    {
      if(IsSemiconductor(material1) && IsResistance(material2) ) return true;
      return false;
    }
    case IF_Metal_Schottky:
    {
      if(IsSemiconductor(material1) && IsResistance(material2) ) return true;
      return false;
    }
    case SchottkyContact:
    {
      if(IsSemiconductor(material1) && material2.empty() ) return true;
      if(IsSemiconductor(material1) && IsInsulator(material2) ) return true;
      if(IsSemiconductor(material1) && IsConductor(material2) ) return true;
      if(IsSemiconductor(material1) && IsResistance(material2) ) return true;
      return false;
    }
    case GateContact:
    {
      if(IsSemiconductor(material1) && material2.empty() ) return true;
      if(IsConductor(material1) && material2.empty() ) return true;
      if(IsResistance(material1) && material2.empty() ) return true;
      if(IsInsulator(material1) && material2.empty() ) return true;
      if(IsInsulator(material1) && IsConductor(material2) ) return true;
      if(IsInsulator(material1) && IsResistance(material2) ) return true;
      if(IsSemiconductor(material1) && IsInsulator(material2) ) return true;
      return false;
    }
    case SimpleGateContact:
    {
      if(IsSemiconductor(material1) && material2.empty() ) return true;
      return false;
    }
    case IF_Insulator_Semiconductor:
    {
      if(IsSemiconductor(material1) && IsInsulator(material2) ) return true;
      return false;
    }
    case SolderPad:
    {
      if(IsResistance(material1) && material2.empty() ) return true;
      if(IsResistance(material1) && IsInsulator(material2) ) return true;
      if(IsResistance(material2) && IsInsulator(material1) ) return true;
      return false;
    }
    case HeteroInterface:
    {
      if(IsSemiconductor(material1) && IsSemiconductor(material2) ) return true;
      return false;
    }
    case ChargeEmit:
    {
      if(IsSemiconductor(material1) && material2.empty() ) return true;
      return false;
    }
  }

  return true;
}



BoundaryCondition::BoundaryCondition(SimulationSystem  & system, const std::string & label)
  : _system(system), _boundary_name(label), _boundary_id(BoundaryInfo::invalid_id), _bc_regions(NULL, NULL), _link_to_spice(false),
    _ext_circuit(NULL), _z_width(system.z_width()),
    _T_Ext(system.T_external()), _inter_connect_hub(0)
{
  for(int i=0; i<4; ++i)
  {
    _global_offset[i] = invalid_uint;
    _local_offset[i]  = invalid_uint;
    _array_offset[i]  = invalid_uint;
  }
}


BoundaryCondition::~BoundaryCondition()
{
  _bd_nodes.clear();
  _bd_fvm_nodes.clear();
  delete _ext_circuit;
}


void BoundaryCondition::insert(const Node * node, SimulationRegion *r, FVM_Node *fn)
{
  //assert(fn->boundary_id() == _boundary_id);
  std::pair<SimulationRegion *, FVM_Node *> mix=std::pair<SimulationRegion *, FVM_Node *>(r, fn);
  (_bd_fvm_nodes[node]).insert(std::pair<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> >(r->type(), mix));
}


std::pair<unsigned int, unsigned int> BoundaryCondition::bc_subdomains() const
{
  std::pair<unsigned int, unsigned int> subdomains;
  subdomains.first = _bc_regions.first ? _bc_regions.first->subdomain_id() : invalid_uint;
  subdomains.second = _bc_regions.second ? _bc_regions.second->subdomain_id() : invalid_uint;
  return subdomains;
}


FVM_Node * BoundaryCondition::get_region_fvm_node(const Node * n, const SimulationRegion * region) const
{
  const_region_node_iterator reg_it     = region_node_begin(n);
  const_region_node_iterator reg_it_end = region_node_end(n);
  for(; reg_it!=reg_it_end; ++reg_it )
    if( (*reg_it).second.first == const_cast<SimulationRegion *>(region) )
      return (*reg_it).second.second;

  return NULL;
}


FVM_Node * BoundaryCondition::get_region_fvm_node(const Node * n, unsigned int subdomain) const
{
  const_region_node_iterator reg_it     = region_node_begin(n);
  const_region_node_iterator reg_it_end = region_node_end(n);
  for(; reg_it!=reg_it_end; ++reg_it )
    if( (*reg_it).second.first->subdomain_id() == subdomain )
      return (*reg_it).second.second;

  return NULL;
}


std::vector<SimulationRegion *>  BoundaryCondition::extra_regions(const Node * n) const
{
  std::vector<SimulationRegion *> regions;

  typedef std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> > RegionNodeMap;
  const RegionNodeMap & region_map = _bd_fvm_nodes.find(n)->second;
  for(RegionNodeMap::const_iterator it = region_map.begin(); it != region_map.end(); it++)
  {
    SimulationRegion * r = it->second.first;
    if( r!= _bc_regions.first && r!= _bc_regions.second )
      regions.push_back(r);
  }

  return regions;
}




bool BoundaryCondition::has_associated_region(const Node * n, SimulationRegionType rt) const
{
  typedef std::map<const Node *, std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> > >::const_iterator It;
  It it = _bd_fvm_nodes.find(n);
  return it->second.find(rt)!=it->second.end();
}

unsigned int BoundaryCondition::n_node_neighbors(const Node * n) const
{
  assert (n->on_processor());

  std::set<const Node *> node_set;

  const_region_node_iterator  reg_it     = region_node_begin(n);
  const_region_node_iterator  reg_it_end = region_node_end(n);
  for(; reg_it!=reg_it_end; ++reg_it )
  {
    const FVM_Node * fvm_node = (*reg_it).second.second;
    FVM_Node::fvm_neighbor_node_iterator it = fvm_node->neighbor_node_begin();
    FVM_Node::fvm_neighbor_node_iterator it_end = fvm_node->neighbor_node_end();
    for(; it!=it_end; ++it)
    {
      const Node * node = (*it).first->root_node();
      if( node->on_processor() )
        node_set.insert(node);
    }
  }
  return  node_set.size();
}


std::vector<const Node *> BoundaryCondition::node_neighbors(const Node * n) const
{
  assert (n->on_processor());

  std::set<const Node *> node_set;

  const_region_node_iterator  reg_it     = region_node_begin(n);
  const_region_node_iterator  reg_it_end = region_node_end(n);
  for(; reg_it!=reg_it_end; ++reg_it )
  {
    const FVM_Node * fvm_node = (*reg_it).second.second;
    FVM_Node::fvm_neighbor_node_iterator it = fvm_node->neighbor_node_begin();
    FVM_Node::fvm_neighbor_node_iterator it_end = fvm_node->neighbor_node_end();
    for(; it!=it_end; ++it)
    {
      const Node * node = (*it).first->root_node();
      if( node->on_processor() )
        node_set.insert(node);
    }
  }

  std::vector<const Node *> neighbor_nodes;
  for(std::set<const Node *>::const_iterator it= node_set.begin(); it != node_set.end(); ++it)
    neighbor_nodes.push_back(*it);
  return  neighbor_nodes;
}



void BoundaryCondition::set_boundary_info_to_fvm_node()
{
  std::map<const Node *, std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> > >::iterator  it = _bd_fvm_nodes.begin();
  for( ; it != _bd_fvm_nodes.end(); ++it)
  {
    std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> > & boundary_fvm_nodes = it->second;
    std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> >::iterator fvm_node_it = boundary_fvm_nodes.begin();
    for(; fvm_node_it != boundary_fvm_nodes.end(); ++fvm_node_it)
    {
      FVM_Node * fvm_node = fvm_node_it->second.second;
      genius_assert(fvm_node->boundary_id() == BoundaryInfo::invalid_id);
      fvm_node->set_bc_type(this->bc_type());
      fvm_node->set_boundary_id(_boundary_id);
    }
  }
}


double BoundaryCondition::scalar(const std::string & name) const
{
  genius_assert(_real_parameters.find(name) != _real_parameters.end());
  return _real_parameters.find(name)->second;
}

double & BoundaryCondition::scalar(const std::string & name)
{
  return _real_parameters[name];
}

bool BoundaryCondition::has_scalar(const std::string & name) const
{
  return _real_parameters.find(name) != _real_parameters.end();
}


std::vector<double> BoundaryCondition::scalar_array(const std::string & name) const
{
  genius_assert(_real_array_parameters.find(name) != _real_array_parameters.end());
  return _real_array_parameters.find(name)->second;
}

std::vector<double> & BoundaryCondition::scalar_array(const std::string & name)
{
  return _real_array_parameters[name];
}

bool BoundaryCondition::has_scalar_array(const std::string & name) const
{
  return _real_array_parameters.find(name) != _real_array_parameters.end();
}


bool BoundaryCondition::flag(const std::string & name) const
{
  genius_assert(_bool_parameters.find(name) != _bool_parameters.end());
  return _bool_parameters.find(name)->second;
}

bool & BoundaryCondition::flag(const std::string & name)
{
  return _bool_parameters[name];
}

bool BoundaryCondition::has_flag(const std::string & name) const
{
  return _bool_parameters.find(name) != _bool_parameters.end();
}


//---------------------------------------------------------------------------------
// constructors for each derived class
//---------------------------------------------------------------------------------
// all the derived boundary conditions
#include "boundary_condition_abs.h"
#include "boundary_condition_ii.h"
#include "boundary_condition_charge.h"
#include "boundary_condition_ir.h"
#include "boundary_condition_is.h"
#include "boundary_condition_ei.h"
#include "boundary_condition_neumann.h"
#include "boundary_condition_gate.h"
#include "boundary_condition_poly_gate.h"
#include "boundary_condition_ohmic.h"
#include "boundary_condition_resistance_ohmic.h"
#include "boundary_condition_rr.h"
#include "boundary_condition_schottky.h"
#include "boundary_condition_resistance_schottky.h"
#include "boundary_condition_hetero.h"
#include "boundary_condition_simplegate.h"
#include "boundary_condition_homo.h"
#include "boundary_condition_solderpad.h"
#include "boundary_condition_iv.h"
#include "boundary_condition_ev.h"
#include "boundary_condition_ee.h"
#include "boundary_condition_charge_emit.h"

#include "boundary_condition_electrode_interconnect.h"
#include "boundary_condition_charge_integral.h"

ElectrodeInterConnectBC::ElectrodeInterConnectBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Inter Connect..."<<std::endl; RECORD();
}


ChargeIntegralBC::ChargeIntegralBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  //std::cout<<"  Initializing \""<< label <<"\" as Charge Integral..."<<std::endl;
  MESSAGE<<"  Initializing \""<< label <<"\" as Charge Integral..."<<std::endl; RECORD();
  scalar("qf")  = 0*C;
  scalar("qf.gen")  = 0*C/s;
  flag("leakage") = false;
}


NeumannBC::NeumannBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Neumann BC..."<<std::endl; RECORD();
  scalar("heat.transfer") = 0.0*J/s/pow(cm,2)/K;
  flag("reflection") = false;
  flag("surface.recombination") = false;
}



OhmicContactBC::OhmicContactBC(SimulationSystem  & system, const std::string & label, bool _interface): BoundaryCondition(system,label)
{
  if(_interface) _bd_type = INTERFACE;
  else _bd_type = BOUNDARY;

  MESSAGE<<"  Initializing \""<< label <<"\" as Ohmic Contact BC..."<<std::endl; RECORD();
  scalar("heat.transfer") = 1e3*J/s/pow(cm,2)/K;
  flag("reflection") = false;
}


IF_Metal_OhmicBC::IF_Metal_OhmicBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Resistance Ohmic Contact BC..."<<std::endl; RECORD();
  scalar("elec.recomb.velocity") = std::numeric_limits<PetscScalar>::infinity();
  scalar("hole.recomb.velocity") = std::numeric_limits<PetscScalar>::infinity();
  //scalar("elec.recomb.velocity") = 2.573e6*cm/s;// for silicon
  //scalar("hole.recomb.velocity") = 1.93e6*cm/s;// for silicon
  _current_flow = 0.0;
}


SchottkyContactBC::SchottkyContactBC(SimulationSystem  & system, const std::string & label, bool _interface): BoundaryCondition(system,label)
{
  if(_interface) _bd_type = INTERFACE;
  else _bd_type = BOUNDARY;

  MESSAGE<<"  Initializing \""<< label <<"\" as Schottky Contact BC..."<<std::endl; RECORD();
  scalar("heat.transfer") = 1e3*J/s/pow(cm,2)/K;
  scalar("workfunction")  = 0.0;
  flag("reflection") = false;
}


IF_Metal_SchottkyBC::IF_Metal_SchottkyBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Resistance Schottky Contact BC..."<<std::endl; RECORD();
  _current_flow = 0.0;
}


GateContactBC::GateContactBC(SimulationSystem  & system, const std::string & label, bool _interface): BoundaryCondition(system,label)
{
  if(_interface) _bd_type = INTERFACE;
  else _bd_type = BOUNDARY;

  MESSAGE<<"  Initializing \""<< label <<"\" as Gate BC..."<<std::endl; RECORD();
  scalar("heat.transfer") = 1e3*J/s/pow(cm,2)/K;
  scalar("workfunction") = 0.0;
}


PolyGateContactBC::PolyGateContactBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as PolyGate BC..."<<std::endl; RECORD();
  scalar("heat.transfer") = 1e3*J/s/pow(cm,2)/K;
  scalar("workfunction") = 0.0;
}


SimpleGateContactBC::SimpleGateContactBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Simple Gate BC..."<<std::endl; RECORD();
  scalar("heat.transfer") = 1e3*J/s/pow(cm,2)/K;
  scalar("workfunction")  = 0.0;
  scalar("thickness")     = 1e-9*cm;
  scalar("eps")           = 3.9*eps0;
  scalar("qf")            = 1e10/pow(cm,2);
}


ResistanceInsulatorBC::ResistanceInsulatorBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Resistance to Insulator BC..."<<std::endl; RECORD();
  _current_flow = 0.0;
}


ResistanceResistanceBC::ResistanceResistanceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Resistance Resistance Interface..."<<std::endl; RECORD();
  _current_flow = 0.0;
}


SolderPadBC::SolderPadBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Solder Pad BC..."<<std::endl; RECORD();
  scalar("heat.transfer") = 1e3*J/s/pow(cm,2)/K;
}

InsulatorSemiconductorInterfaceBC::InsulatorSemiconductorInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Insulator Semiconductor Interface..."<<std::endl; RECORD();
  scalar("qf")  = 1e10/pow(cm,2);
  scalar("J_CBET") = 0.0;
  scalar("J_VBHT") = 0.0;
  scalar("J_VBET") = 0.0;
  flag("surface.recombination") = true;
  _current_flow = 0.0;
}



HomoInterfaceBC::HomoInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Semiconductor Homo Junction..."<<std::endl; RECORD();
}



HeteroInterfaceBC::HeteroInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Semiconductor Hetero Junction..."<<std::endl; RECORD();
  scalar("qf")  = 0.0;
}



ElectrodeInsulatorInterfaceBC::ElectrodeInsulatorInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Electrode Insulator Interface..."<<std::endl; RECORD();
}



InsulatorInsulatorInterfaceBC::InsulatorInsulatorInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Insulator Insulator Interface..."<<std::endl; RECORD();
}


ChargedContactBC::ChargedContactBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Float Metal with Free Charge..."<<std::endl; RECORD();
  _current_flow = 0.0;
}


ChargeEmitBC::ChargeEmitBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Charge Emit boundary..."<<std::endl; RECORD();
  _current_flow = 0.0;
}


AbsorbingBC::AbsorbingBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Absorbing boundary..."<<std::endl; RECORD();
}


InsulatorVacuumInterfaceBC::InsulatorVacuumInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Insulator Vacuum Interface..."<<std::endl; RECORD();
}

ElectrodeVacuumInterfaceBC::ElectrodeVacuumInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Electrode Vacuum Interface..."<<std::endl; RECORD();
}

ElectrodeElectrodeInterfaceBC::ElectrodeElectrodeInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Electrode Electrode Interface..."<<std::endl; RECORD();
}

//---------------------------------------------------------------------------------
// the string which indicates the boundary condition
//---------------------------------------------------------------------------------

std::string NeumannBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
     <<"string<id>="<<this->label()<<" "
     <<"enum<type>=Neumann "
     <<"real<ext.temp>="<<T_external()/K<<" "
     <<"real<heat.transfer>="<<scalar("heat.transfer")/(J/s/pow(cm,2)/K)<<" "
     <<"bool<reflection>="<<( flag("reflection") ? "true" : "false" )<<" "
     <<"bool<surface.recombination>="<<( flag("surface.recombination") ? "true" : "false" )<<" ";

  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;

  ss<<std::endl;

  return ss.str();
}


std::string OhmicContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"enum<type>=OhmicContact ";

  if(this->ext_circuit())
  {
    ss << this->ext_circuit()->format() <<" ";
  }

  if(_bd_type == BOUNDARY)
  {
    ss <<"real<ext.temp>="<<T_external()/K<<" "
       <<"real<heat.transfer>="<<scalar("heat.transfer")/(J/s/pow(cm,2)/K)<<" "
       <<"bool<reflection>="<<( flag("reflection") ? "true" : "false" )<<" ";
  }
  else
  {
    if(!this->electrode_label().empty())
    {
      ss << "string<electrode_id>=" << this->electrode_label() << " ";
    }
  }

  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;
  ss<<std::endl;

  return ss.str();
}


std::string IF_Metal_OhmicBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"enum<type>=resistanceohmiccontact ";

  if( !_infinity_recombination() )
    ss  <<"real<elec.recomb.velocity>="<<scalar("elec.recomb.velocity")/(cm/s)<<" "
        <<"real<hole.recomb.velocity>="<<scalar("hole.recomb.velocity")/(cm/s)<<" ";

  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;

  ss<<std::endl;

  return ss.str();
}



std::string SchottkyContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
     <<"string<id>="<<this->label()<<" "
     <<"enum<type>=SchottkyContact "
     <<"real<workfunction>="<<scalar("workfunction")/V<<" ";

  if(this->ext_circuit())
  {
    ss << this->ext_circuit()->format() <<" ";
  }

  if(_bd_type == BOUNDARY)
  {
    ss <<"real<ext.temp>="<<T_external()/K<<" "
       <<"real<heat.transfer>="<<scalar("heat.transfer")/(J/s/pow(cm,2)/K)<<" "
       <<"bool<reflection>="<<( flag("reflection") ? "true" : "false" )<<" ";
  }
  else
  {
    if(!this->electrode_label().empty())
    {
      ss << "string<electrode_id>=" << this->electrode_label() << " ";
    }
  }

  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;
  ss<<std::endl;

  return ss.str();
}


std::string IF_Metal_SchottkyBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"enum<type>=resistanceschottkycontact ";
  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;

  ss<<std::endl;

  return ss.str();
}



std::string GateContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
     <<"string<id>="<<this->label()<<" "
     <<"enum<type>=GateContact "
     <<"real<workfunction>="<<scalar("workfunction")/V<<" ";

  if(this->ext_circuit())
  {
    ss << this->ext_circuit()->format() <<" ";
  }

  if(_bd_type == BOUNDARY)
  {
    ss <<"real<ext.temp>="<<T_external()/K<<" "
       <<"real<heat.transfer>="<<scalar("heat.transfer")/(J/s/pow(cm,2)/K)<<" ";
  }
  else
  {
    if(!this->electrode_label().empty())
    {
      ss << "string<electrode_id>=" << this->electrode_label() << " ";
    }
  }

  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;
  ss<<std::endl;

  return ss.str();
}


std::string PolyGateContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"enum<type>=PolyGateContact "
      <<"real<workfunction>="<<scalar("workfunction")/V<<" ";

  if(this->ext_circuit())
  {
    ss << this->ext_circuit()->format() <<" ";
  }

  if(_bd_type == BOUNDARY)
  {
    ss <<"real<ext.temp>="<<T_external()/K<<" "
        <<"real<heat.transfer>="<<scalar("heat.transfer")/(J/s/pow(cm,2)/K)<<" ";
  }
  else
  {
    if(!this->electrode_label().empty())
    {
      ss << "string<electrode_id>=" << this->electrode_label() << " ";
    }
  }

  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;
  ss<<std::endl;

  return ss.str();
}


std::string SimpleGateContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"enum<type>=SimpleGateContact "
  <<"real<workfunction>="<<scalar("workfunction")/V<<" "
  <<"real<thickness>="<<scalar("thickness")/um<<" "
  <<"real<eps>="<<scalar("eps")<<" "
  <<"real<qf>="<<scalar("qf")*pow(cm,2)<<" ";

  if(this->ext_circuit())
  {
    ss << this->ext_circuit()->format() <<" ";
  }

  ss<<"real<ext.temp>="<<T_external()/K<<" "
    <<"real<heat.transfer>="<<scalar("heat.transfer")/(J/s/pow(cm,2)/K)<<" ";

  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;

  ss<<std::endl;

  return ss.str();
}


std::string ResistanceInsulatorBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"#BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"  <<the interface between resistance region and insulator region>>"
      <<std::endl;

  return ss.str();
}


std::string ResistanceResistanceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"#BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"  <<the interface between two resistance region>>"
      <<std::endl;

  return ss.str();
}


std::string SolderPadBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss  <<"BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"enum<type>=SolderPad ";

  if(this->ext_circuit())
  {
    ss << this->ext_circuit()->format() <<" ";
  }

  ss  <<"real<ext.temp>="<<T_external()/K<<" "
      <<"real<heat.transfer>="<<scalar("heat.transfer")/(J/s/pow(cm,2)/K)<<" ";

  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;

  ss<<std::endl;

  return ss.str();
}


std::string InsulatorSemiconductorInterfaceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"enum<type>=InsulatorInterface "
  <<"real<qf>="<<scalar("qf")*pow(cm,2)<<" "
  <<"bool<surface.recombination>="<<( flag("surface.recombination") ? "true" : "false" )
  <<std::endl;

  return ss.str();
}


std::string HomoInterfaceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"#BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"  <<the interface between two regions with same semiconductor material>>"
  <<std::endl;

  return ss.str();
}


std::string HeteroInterfaceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"#BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"  <<the interface between two regions with different semiconductor material>>"
  <<std::endl;

  return ss.str();
}


std::string ElectrodeInsulatorInterfaceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"#BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"  <<the interface between electrode and insulator region>>"
  <<std::endl;

  return ss.str();
}


std::string InsulatorInsulatorInterfaceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"#BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"  <<the interface between two insulator region>>"
  <<std::endl;

  return ss.str();
}


std::string ChargedContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
     <<"string<id>="<<this->label()<<" "
     <<"enum<type>=ChargedContact ";

  if(system().mesh().mesh_dimension () == 2)
    ss<<"real<z.width>="<<z_width()/um;

  ss<<std::endl;

  return ss.str();
}


std::string ElectrodeInterConnectBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss  <<"INTERCONNECT "
      <<"string<id>="<<this->label()<<" ";

  if(this->ext_circuit())
  {
    ss << this->ext_circuit()->format() <<" ";
  }

  ss  <<"bool<float>="<<( this->ext_circuit()->is_float() ? "true" : "false" )<<" ";

  const std::vector<BoundaryCondition * > & inter_connect_bcs = this->inter_connect();

  for(unsigned int n=0; n<inter_connect_bcs.size(); n++)
    ss<<"string<connectto>="<<inter_connect_bcs[n]->label()<<" ";

  ss<<std::endl;

  return ss.str();
}


std::string ChargeIntegralBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"CHARGE "
      <<"string<id>="<<this->label()<<" "
      <<"real<charge>="<<this->scalar("qf")/C<<" ";

  const std::vector<BoundaryCondition * > & inter_connect_bcs = this->inter_connect();

  for(unsigned int n=0; n<inter_connect_bcs.size(); n++)
    ss<<"string<chargeboundary>="<<inter_connect_bcs[n]->label()<<" ";

  ss<<std::endl;

  return ss.str();
}


std::string ChargeEmitBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"enum<type>=ChargeEmit "
      <<std::endl;

  return ss.str();
}


std::string AbsorbingBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"enum<type>=AbsorbingBoundary "
  <<std::endl;

  return ss.str();
}


std::string InsulatorVacuumInterfaceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"#BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"  <<the interface between insulator region and vacuum region>>"
      <<std::endl;

  return ss.str();
}


std::string ElectrodeVacuumInterfaceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"#BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"  <<the interface between electrode region and vacuum region>>"
      <<std::endl;

  return ss.str();
}


std::string ElectrodeElectrodeInterfaceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"#BOUNDARY "
      <<"string<id>="<<this->label()<<" "
      <<"  <<the interface between electrode region and electrode region>>"
      <<std::endl;

  return ss.str();
}
