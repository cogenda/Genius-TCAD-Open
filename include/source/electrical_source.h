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

//  $Id: source.h,v 1.9 2008/07/09 09:10:08 gdiso Exp $

#ifndef __electrical_source_h__
#define __electrical_source_h__

#include <map>

#include "genius_common.h"
#include "physical_unit.h"
#include "vsource.h"
#include "isource.h"


//pre define
namespace Parser{
  class InputParser;
  class Card;
}
class BoundaryConditionCollector;


/**
 * class for manage all the electrical sources
 */
class  ElectricalSource
{
public:
  /**
   * constructor, init all the source defined in input deck
   */
  ElectricalSource(Parser::InputParser & decks);

  /**
   * free all the sources
   */
  ~ElectricalSource();

  /**
   * hold the bcs pointer
   */
  void link_to_bcs(BoundaryConditionCollector * bcs);

  /**
   * attach one source to a electrode bc.
   * delete any source(s) previously attached to this bc
   * @note this function does not set the vapp or iapp
   * for each electrode, user should call update(time)
   * to do it
   */
  void attach_source_to_electrode(const std::string & electrode_label, const std::string & source_list );

  /**
   * attach one or more sources to a electrode bc.
   * delete any source(s) previously attached to this bc
   * @note this function does not set the vapp or iapp
   * for each electrode, user should call update(time)
   * to do it
   */
  void attach_sources_to_electrode(const std::string & electrode_label, const std::vector<std::string> & source_list );


  /**
   * create VDC with given vconst and attach it to a electrode bc.
   */
  void attach_voltage_to_electrode(const std::string & electrode_label, double vconst );


  /**
   * create IDC with given iconst and attach current to a electrode bc.
   */
  void attach_current_to_electrode(const std::string & electrode_label, double iconst );


  /**
   * @return true if source_name can be matched in vsource list
   */
  bool is_vsource_exist(const std::string source_name) const
  { return _vsource_list.find(source_name) != _vsource_list.end(); }

  /**
   * @return true if source_name can be matched in isource list
   */
  bool is_isource_exist(const std::string source_name) const
  { return _isource_list.find(source_name) != _isource_list.end(); }


  /**
   * remove all the source(s) attached with some electrode
   */
  void remove_electrode_source(const std::string & electrode_label);

  /**
   * @return limited time step by max allowd voltage/current change
   */
  double limit_dt(double time, double dt, double dt_min, double v_change, double i_change) const;

  /**
   * update Vapp or Iapp for all the electrode bcs to new time step
   * @note the default vapp/iapp is 0 for all the electrode
   */
  void update(double time);

  /**
   * @return steps for ramping up sources by given voltage and current limiter
   * @param time   gives the intensity of sources
   */
  unsigned int steps_by_limiter(double v_change, double i_change, double time=0.0) const;

  /**
   * ramp up Vapp or Iapp for all the electrode bcs to a*A(time)
   * @param time   gives the intensity of sources
   */
  void rampup(double a, double time=0.0);

  /**
   * force the voltage of electrode to be vapp
   * useful for DC sweep
   */
  void assign_voltage_to(const std::string & electrode_label, double vapp);

  /**
   * force the voltage of some electrodes to be vapp
   * useful for DC sweep
   */
  void assign_voltage_to(const std::vector<std::string> & electrode_labels, double vapp);

  /**
   * force the current of electrode to be iapp
   * useful for DC sweep
   */
  void assign_current_to(const std::string & electrode_label, double iapp);

  /**
   * force the current of some electrodes to be iapp
   * useful for DC sweep
   */
  void assign_current_to(const std::vector<std::string> & electrode_labels, double iapp);

  /**
   * @return the vapp of given bc at time
   */
  double vapp(const std::string &bc, double time) const;

  /**
   * @return the iapp of given bc at time
   */
  double iapp(const std::string &bc, double time) const;

  /**
   * save the source state of each bc
   */
  void save_bc_source_state();


  /**
   * clear the bc information this class holds.
   * when the system is rebuild or a new system imported
   */
  void clear_bc_source_map()
  { _bc_source_map.clear(); }

private:

 /**
  * the bc collector which electrical sources will be assigned to
  */
 BoundaryConditionCollector * _bcs;

 /**
  * counter for implicit defined source
  */
 unsigned int _counter;

 /**
  * all the voltage sources exist
  */
 std::map<std::string, VSource * >  _vsource_list;

 /**
  * all the current sources exist
  */
 std::map<std::string, ISource * >  _isource_list;


 /**
  * the map of electrode bc label to the sources it owns
  */
 std::map<std::string, std::pair<std::vector<VSource *>, std::vector<ISource *> > > _bc_source_map;

 typedef std::map<std::string, std::pair<std::vector<VSource *>, std::vector<ISource *> > >::iterator BIt;

 typedef std::map<std::string, std::pair<std::vector<VSource *>, std::vector<ISource *> > >::const_iterator CBIt;

 /** std::string -- the bc label
  * int          -- cast to ExternalCircuit::DRIVEN
  * double       -- the value of driven force
  */
 std::map<std::string, std::pair<int, double> > _bc_source_state;

 /**
  * @return the vapp state of given bc in _bc_source_state
  */
 double _bc_state_vapp(const std::string &bc) const;

 /**
  * @return the iapp state of given bc in _bc_source_state
  */
 double _bc_state_iapp(const std::string &bc) const;

 /**
  * private functions for setting each source
  */
 void  SetVDC(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetVSIN(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetVEXP(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetVPULSE(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetVUSER(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetVSHELL(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetIDC(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetISIN(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetIEXP(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetIPULSE(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetIUSER(const Parser::Card &c);

 /**
  * private functions for setting each source
  */
 void  SetISHELL(const Parser::Card &c);

};



#endif //#define __electrical_source_h__
