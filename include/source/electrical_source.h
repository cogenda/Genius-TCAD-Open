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

#include "genius_common.h"
#include "parser.h"
#include "physical_unit.h"
#include "vsource.h"
#include "isource.h"


//pre define
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
   * attach one or more sources to a electrode bc.
   * delete any source(s) previously attached to this bc
   * @note this function does not set the vapp or iapp
   * for each electrode, user should call update(time)
   * to do it
   */
  void attach_sources_to_electrode(const std::string & electrode_label, const std::vector<std::string> & source_list );


  /**
   * attach voltage to a electrode bc. delete any source(s) previously attached to this bc
   * the voltage takes effect at once
   */
  void attach_voltage_to_electrode(const std::string & electrode_label, PetscScalar vconst );


  /**
   * attach current to a electrode bc. delete any source(s) previously attached to this bc
   * the current takes effect at once
   */
  void attach_current_to_electrode(const std::string & electrode_label, PetscScalar iconst );


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
  PetscScalar limit_dt(PetscScalar time, PetscScalar dt, PetscScalar v_change, PetscScalar i_change) const;

  /**
   * update Vapp or Iapp for all the electrode bcs to new time step
   * @note the default vapp/iapp is 0 for all the electrode
   */
  void update(PetscScalar time);

  /**
   * force the voltage of electrode to be vapp
   * useful for DC sweep
   */
  void assign_voltage_to(const std::string & electrode_label, PetscScalar vapp);

  /**
   * force the voltage of some electrodes to be vapp
   * useful for DC sweep
   */
  void assign_voltage_to(const std::vector<std::string> & electrode_labels, PetscScalar vapp);

  /**
   * force the current of electrode to be iapp
   * useful for DC sweep
   */
  void assign_current_to(const std::string & electrode_label, PetscScalar iapp);

  /**
   * force the current of some electrodes to be iapp
   * useful for DC sweep
   */
  void assign_current_to(const std::vector<std::string> & electrode_labels, PetscScalar iapp);

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
  * all the voltage sources exist
  */
 std::map<const std::string, VSource * >  _vsource_list;

 /**
  * all the current sources exist
  */
 std::map<const std::string, ISource * >  _isource_list;


 /**
  * the map of electrode bc label to the sources it owns
  */
 std::map<const std::string, std::pair<std::vector<VSource *>, std::vector<ISource *> > > _bc_source_map;

 typedef std::map<const std::string, std::pair<std::vector<VSource *>, std::vector<ISource *> > >::iterator BIt;

 typedef std::map<const std::string, std::pair<std::vector<VSource *>, std::vector<ISource *> > >::const_iterator CBIt;

 PetscScalar _vapp(const std::string &bc, PetscScalar time) const;

 PetscScalar _iapp(const std::string &bc, PetscScalar time) const;

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
