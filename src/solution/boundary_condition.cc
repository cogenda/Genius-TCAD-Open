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


#include "simulation_system.h"
#include "simulation_region.h"
#include "boundary_condition.h"

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
PetscScalar BoundaryCondition::_dummy_ = 0;
bool BoundaryCondition::_bool_dummy_ = false;


static std::map<const std::string, BCType> bc_name_to_bc_type;


void init_bc_name_to_bc_type_map()
{
  // if the bc_name_to_bc_type not build yet
  if ( bc_name_to_bc_type.empty() )
  {
    bc_name_to_bc_type["neumann"           ]  = NeumannBoundary;
    bc_name_to_bc_type["ohmiccontact"      ]  = OhmicContact;
    bc_name_to_bc_type["schottkycontact"   ]  = SchottkyContact;
    bc_name_to_bc_type["gatecontact"       ]  = GateContact;
    bc_name_to_bc_type["simplegatecontact" ]  = SimpleGateContact;
    bc_name_to_bc_type["insulatorinterface"]  = IF_Insulator_Semiconductor;
    bc_name_to_bc_type["heterojunction"    ]  = HeteroInterface;
    bc_name_to_bc_type["chargedcontact"    ]  = ChargedContact;
    bc_name_to_bc_type["absorbingboundary" ]  = AbsorbingBoundary;
    bc_name_to_bc_type["sourceboundary"    ]  = SourceBoundary;
    bc_name_to_bc_type["floatmetal"        ]  = ChargedContact;
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
    bc_type_to_bc_name[NeumannBoundary            ]  = std::string("NeumannBoundary");
    bc_type_to_bc_name[OhmicContact               ]  = std::string("OhmicContact");
    bc_type_to_bc_name[SchottkyContact            ]  = std::string("SchottkyContact");
    bc_type_to_bc_name[GateContact                ]  = std::string("GateContact");
    bc_type_to_bc_name[SimpleGateContact          ]  = std::string("SimpleGateContact");
    bc_type_to_bc_name[IF_Insulator_Semiconductor ]  = std::string("IF_Insulator_Semiconductor");
    bc_type_to_bc_name[HeteroInterface            ]  = std::string("HeteroInterface");
    bc_type_to_bc_name[ChargedContact             ]  = std::string("ChargedContact");
    bc_type_to_bc_name[AbsorbingBoundary          ]  = std::string("AbsorbingBoundary");
    bc_type_to_bc_name[SourceBoundary             ]  = std::string("SourceBoundary");
    bc_type_to_bc_name[ChargedContact             ]  = std::string("ChargedContact");
  }
}

std::string BC_enum_to_string(BCType bc_type)
{
  init_bc_type_to_bc_name_map();

  if ( bc_type_to_bc_name.find(bc_type)!= bc_type_to_bc_name.end() )
    return bc_type_to_bc_name[bc_type];
  else
    return std::string("AutoSet");
}


BCType determine_bc_by_subdomain(const std::string & mat1, const std::string mat2)
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
  }

  // two materials are both Conductor
  if(IsConductor(mat1) && IsConductor(mat2) )
    return IF_Electrode_Electrode;

  // one material is Conductor
  if(IsConductor(mat1) || IsConductor(mat2))
  {
    if(IsVacuum(mat1) || IsVacuum(mat2))
      return IF_Electrode_Vacuum;
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




BoundaryCondition::BoundaryCondition(SimulationSystem  & system, const std::string & label)
    : _system(system), _boundary_name(label), _ext_circuit(NULL), _z_width(system.z_width()), _T_Ext(system.T_external()), _inter_connect_hub(0)
{
  _global_offset = invalid_uint;
  _local_offset  = invalid_uint;
  _array_offset  = invalid_uint;
}


void BoundaryCondition::insert(const Node * node, SimulationRegion *r, FVM_Node *fn)
{
  std::pair<SimulationRegion *, FVM_Node *> mix=std::pair<SimulationRegion *, FVM_Node *>(r, fn);
  (_bd_fvm_nodes[node]).insert(std::pair<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> >(r->type(), mix));
}




//---------------------------------------------------------------------------------
// constructors for each derived class
//---------------------------------------------------------------------------------

InterConnectBC::InterConnectBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Inter Connect..."<<std::endl; RECORD();
}



NeumannBC::NeumannBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Neumann BC..."<<std::endl; RECORD();
  _Heat_Transfer = 0.0;
  _reflection    = false;
}



OhmicContactBC::OhmicContactBC(SimulationSystem  & system, const std::string & label, bool _interface): BoundaryCondition(system,label)
{
  if(_interface) _bd_type = INTERFACE;
  else _bd_type = BOUNDARY;

  MESSAGE<<"  Initializing \""<< label <<"\" as OhmicContact BC..."<<std::endl; RECORD();
  _Heat_Transfer = 1e3*J/s/pow(cm,2)/K;
}



SchottkyContactBC::SchottkyContactBC(SimulationSystem  & system, const std::string & label, bool _interface): BoundaryCondition(system,label)
{
  if(_interface) _bd_type = INTERFACE;
  else _bd_type = BOUNDARY;

  MESSAGE<<"  Initializing \""<< label <<"\" as Schottky BC..."<<std::endl; RECORD();
  _Heat_Transfer = 1e3*J/s/pow(cm,2)/K;
  _WorkFunction = 0.0;
}



GateContactBC::GateContactBC(SimulationSystem  & system, const std::string & label, bool _interface): BoundaryCondition(system,label)
{
  if(_interface) _bd_type = INTERFACE;
  else _bd_type = BOUNDARY;

  MESSAGE<<"  Initializing \""<< label <<"\" as Gate BC..."<<std::endl; RECORD();
  _Heat_Transfer = 1e3*J/s/pow(cm,2)/K;
  _WorkFunction = 0.0;
}



SimpleGateContactBC::SimpleGateContactBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Simple Gate BC..."<<std::endl; RECORD();
  _Heat_Transfer = 0.0;
  _WorkFunction = 0.0;
  _Thickness = 1e-9*cm;
  _eps = 3.9*eps0;
  _Qf  = 1e10/pow(cm,2);
}



InsulatorSemiconductorInterfaceBC::InsulatorSemiconductorInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Insulator Semiconductor Interface..."<<std::endl; RECORD();
  _Qf  = 1e10/pow(cm,2);
}



HomoInterfaceBC::HomoInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Semiconductor Homo Junction..."<<std::endl; RECORD();
}



HeteroInterfaceBC::HeteroInterfaceBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Semiconductor Hetero Junction..."<<std::endl; RECORD();
  _Qf  = 0.0;
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
  _Qf  = 0;
}


AbsorbingBC::AbsorbingBC(SimulationSystem  & system, const std::string & label): BoundaryCondition(system,label)
{
  MESSAGE<<"  Initializing \""<< label <<"\" as Absorbing boundary..."<<std::endl; RECORD();
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
  <<"real<heat.transfer>="<<_Heat_Transfer/(J/s/pow(cm,2)/K)<<" "
  <<"bool<reflection>="<<( _reflection ? "true" : "false" )
  <<std::endl;

  return ss.str();
}


std::string OhmicContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"enum<type>=OhmicContact "
  <<"real<res>="<<this->ext_circuit()->R()/(V/A)<<" "
  <<"real<cap>="<<this->ext_circuit()->C()/(C/V)<<" "
  <<"real<ind>="<<this->ext_circuit()->L()/(V*s/A)<<" "
  <<"real<potential>="<<this->ext_circuit()->potential()/V<<" ";

  if(_bd_type == BOUNDARY)
  {
    ss <<"real<ext.temp>="<<T_external()/K<<" "
       <<"real<heat.transfer>="<<_Heat_Transfer/(J/s/pow(cm,2)/K)<<" ";
  }
  else
  {
    if(!this->electrode_label().empty())
    {
      ss << "string<electrode_id>=" << this->electrode_label() << " ";
    }
  }

  ss<<std::endl;

  return ss.str();
}

std::string SchottkyContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"enum<type>=SchottkyContact "
  <<"real<workfunction>="<<_WorkFunction/V<<" "
  <<"real<res>="<<this->ext_circuit()->R()/(V/A)<<" "
  <<"real<cap>="<<this->ext_circuit()->C()/(C/V)<<" "
  <<"real<ind>="<<this->ext_circuit()->L()/(V*s/A)<<" "
  <<"real<potential>="<<this->ext_circuit()->potential()/V<<" ";

  if(_bd_type == BOUNDARY)
  {
    ss <<"real<ext.temp>="<<T_external()/K<<" "
       <<"real<heat.transfer>="<<_Heat_Transfer/(J/s/pow(cm,2)/K)<<" ";
  }
  else
  {
    if(!this->electrode_label().empty())
    {
      ss << "string<electrode_id>=" << this->electrode_label() << " ";
    }
  }

  ss<<std::endl;

  return ss.str();
}

std::string GateContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"enum<type>=GateContact "
  <<"real<workfunction>="<<_WorkFunction/V<<" "
  <<"real<res>="<<this->ext_circuit()->R()/(V/A)<<" "
  <<"real<cap>="<<this->ext_circuit()->C()/(C/V)<<" "
  <<"real<ind>="<<this->ext_circuit()->L()/(V*s/A)<<" "
  <<"real<potential>="<<this->ext_circuit()->potential()/V<<" ";

  if(_bd_type == BOUNDARY)
  {
    ss <<"real<ext.temp>="<<T_external()/K<<" "
       <<"real<heat.transfer>="<<_Heat_Transfer/(J/s/pow(cm,2)/K)<<" ";
  }
  else
  {
    if(!this->electrode_label().empty())
    {
      ss << "string<electrode_id>=" << this->electrode_label() << " ";
    }
  }

  ss<<std::endl;

  return ss.str();
}

std::string SimpleGateContactBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"enum<type>=SimpleGateContact "
  <<"real<workfunction>="<<_WorkFunction/V<<" "
  <<"real<thickness>="<<_Thickness/um<<" "
  <<"real<eps>="<<_eps<<" "
  <<"real<qf>="<<_Qf*pow(cm,2)<<" "
  <<"real<res>="<<this->ext_circuit()->R()<<" "
  <<"real<cap>="<<this->ext_circuit()->C()<<" "
  <<"real<ind>="<<this->ext_circuit()->L()<<" "
  <<"real<potential>="<<this->ext_circuit()->potential()<<" "
  <<"real<ext.temp>="<<T_external()/K<<" "
  <<"real<heat.transfer>="<<_Heat_Transfer/(J/s/pow(cm,2)/K)<<" "
  <<std::endl;

  return ss.str();
}


std::string InsulatorSemiconductorInterfaceBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"BOUNDARY "
  <<"string<id>="<<this->label()<<" "
  <<"enum<type>=InsulatorInterface "
  <<"real<qf>="<<_Qf*pow(cm,2)<<" "
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
  <<"enum<type>=ChargedContact "
  <<"real<qf>="<<_Qf/C<<" "
  <<std::endl;

  return ss.str();
}


std::string InterConnectBC::boundary_condition_in_string() const
{
  std::stringstream ss;

  ss <<"INTERCONNECT "
  <<"string<id>="<<this->label()<<" "
  <<"real<res>="<<this->ext_circuit()->R()<<" "
  <<"real<cap>="<<this->ext_circuit()->C()<<" "
  <<"real<ind>="<<this->ext_circuit()->L()<<" "
  <<"real<potential>="<<this->ext_circuit()->potential()<<" "
  <<"bool<float>="<<( this->ext_circuit()->is_float() ? "true" : "false" )<<" ";

  const std::vector<BoundaryCondition * > & inter_connect_bcs = this->inter_connect();

  for(unsigned int n=0; n<inter_connect_bcs.size(); n++)
    ss<<"string<connectto>="<<inter_connect_bcs[n]->label()<<" ";

  ss<<std::endl;

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
