%ModuleHeaderCode
#include "enum_region.h"
%End

// {{{ enum  SimulationRegionType
enum  SimulationRegionType {
  SemiconductorRegion=0,
  InsulatorRegion,
  MetalRegion,
  ElectrodeRegion,
  VacuumRegion,
  PMLRegion,
  InvalidRegion
};
// }}}

// {{{ class SimulationRegion
class SimulationRegion /Abstract, NoDefaultCtors/
{
%TypeHeaderCode
#include "simulation_system.h"
#include "simulation_region.h"
%End
public:
  unsigned int n_node() const;
  const std::string& name() const;
  const std::string& material() const;
  std::string type_name() const;
};
// }}}

// {{{ class ExternalCircuit
class ExternalCircuit /NoDefaultCtors/
{
%TypeHeaderCode
#include "external_circuit.h"
%End
public:
  enum DRIVEN{VDRIVEN, IDRIVEN, FLOAT};

  double R() const;
  double C() const;
  double L() const;
  double Vapp() const;
  double Iapp() const;

  double potential() const;
  double current() const;
};
// }}}

%ModuleHeaderCode
#include "boundary_condition.h"
%End
// {{{ enum BCType {
enum BCType
{
  NeumannBoundary             = 0x0001,
  IF_Semiconductor_Vacuum     = 0x0002,
  IF_Insulator_Vacuum         = 0x0003,
  IF_Electrode_Vacuum         = 0x0004,
  IF_Metal_Vacuum        = 0x0005,
  IF_PML_Scatter              = 0x0006,
  IF_PML_PML                  = 0x0007,
  IF_Electrode_Insulator      = 0x0010,
  IF_Insulator_Semiconductor  = 0x0012,
  IF_Insulator_Insulator      = 0x0013,
  IF_Electrode_Electrode      = 0x0014,
  IF_Electrode_Metal          = 0x0015,
  IF_Insulator_Metal          = 0x0016,
  IF_Metal_Metal              = 0x0017,
  HomoInterface               = 0x0018,
  HeteroInterface             = 0x0019,
  IF_Electrode_Semiconductor  = 0x001a,
  IF_Metal_Semiconductor  = 0x001b,
  OhmicContact        = 0x0101,
  IF_Metal_Ohmic      = 0x0102,
  SchottkyContact     = 0x0103,
  IF_Metal_Schottky   = 0x0104,
  SimpleGateContact   = 0x0105,
  GateContact         = 0x0106,
  ChargedContact      = 0x0107,
  ChargeIntegral      = 0x0108,
  SolderPad           = 0x0109,
  InterConnect        = 0x0110,
  AbsorbingBoundary   = 0x1001,
  SourceBoundary      = 0x1002,
  INVALID_BC_TYPE     = 0xffff           // should always be last

};
// }}}

// {{{ class BoundaryCondition
class BoundaryCondition /Abstract, NoDefaultCtors/
{
%TypeHeaderCode
#include "boundary_condition.h"
%End

public:
  const std::string& label() const;
  const std::string& electrode_label() const;
  std::string bc_type_name() const;
  bool is_electrode() const;
  double z_width() const;
  bool is_inter_connect_bc() const;
  bool is_inter_connect_hub() const;

  std::string boundary_condition_in_string() const;

  double T_external() const;
  double Heat_Transfer() const;
  double Work_Function() const;

};
// }}}

// {{{ class SimulationSystem
class SimulationSystem /NoDefaultCtors/
{
%TypeHeaderCode
#include "simulation_system.h"
#include "boundary_condition_collector.h"
%End

public:
  unsigned int n_regions() const;
  SimulationRegion* region(unsigned int n);
  SimulationRegion* region(const std::string &region_lable);

  unsigned int n_bcs() const;
%MethodCode
    sipRes = sipCpp->get_bcs()->n_bcs();
%End

  unsigned int n_electrode_bcs() const;
%MethodCode
    sipRes = sipCpp->get_bcs()->n_electrode_bcs();
%End

  BoundaryCondition* get_bc(unsigned int n);
%MethodCode
    sipRes = sipCpp->get_bcs()->get_bc(a0);
%End

  BoundaryCondition* get_bc(const std::string &label);
%MethodCode
    sipRes = sipCpp->get_bcs()->get_bc(*a0);
%End

};
// }}}



