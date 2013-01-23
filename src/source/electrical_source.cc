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

#include "parser.h"
#include "electrical_source.h"
#include "simulation_system.h"
#include "boundary_condition_collector.h"

#ifdef WINDOWS
  #include <Windows.h>
  #undef max
  #undef min
#else
  #include <dlfcn.h>
#endif

// for short
using  PhysicalUnit::s;
using  PhysicalUnit::V;
using  PhysicalUnit::A;


ElectricalSource::ElectricalSource(Parser::InputParser & decks):_bcs(0), _counter(0)
{

  //setup voltage, current source and light source
  MESSAGE<<"Setting each voltage and current source here..."; RECORD();

  for( decks.begin(); !decks.end(); decks.next() )
  {
    Parser::Card c = decks.get_current_card();

    if( c.key() == "VSOURCE" )   // It's a VSOURCE card
    {
      if(c.is_parameter_exist("type"))
      {

        if(c.is_enum_value("type","vdc"))            SetVDC(c);

        else if(c.is_enum_value("type","vsin"))      SetVSIN(c);

        else if(c.is_enum_value("type","vexp"))      SetVEXP(c);

        else if(c.is_enum_value("type","vpulse"))    SetVPULSE(c);

        else if(c.is_enum_value("type","vuser"))     SetVUSER(c);

        else if(c.is_enum_value("type","vshell"))    SetVSHELL(c);

      }
      else
      {
        MESSAGE<<"ERROR at " <<c.get_fileline()<< " VSOURCE: Missing parameter type in this statement." << std::endl; RECORD();
        genius_error();
      }
    }

    if( c.key() == "ISOURCE" )  //it is a isource structure
    {

      if(c.is_parameter_exist("type"))
      {

        if(c.is_enum_value("type","idc"))            SetIDC(c);

        else if(c.is_enum_value("type","isin"))      SetISIN(c);

        else if(c.is_enum_value("type","iexp"))      SetIEXP(c);

        else if(c.is_enum_value("type","ipulse"))    SetIPULSE(c);

        else if(c.is_enum_value("type","iuser"))     SetIUSER(c);

        else if(c.is_enum_value("type","ishell"))    SetISHELL(c);

      }
      else
      {
        MESSAGE<<"ERROR at " <<c.get_fileline()<< " ISOURCE: Missing parameter type in this statement." << std::endl; RECORD();
        genius_error();
      }
    }

  }

  MESSAGE<<"done.\n"<<std::endl;  RECORD();

}



ElectricalSource::~ElectricalSource()
{
  std::map<std::string, VSource * >::iterator it1 = _vsource_list.begin();
  for(; it1!=_vsource_list.end(); ++it1)
    delete (*it1).second;
  _vsource_list.clear();

  std::map<std::string, ISource * >::iterator it2 = _isource_list.begin();
  for(; it2!=_isource_list.end(); ++it2)
    delete (*it2).second;
  _isource_list.clear();

}



void ElectricalSource::link_to_bcs(BoundaryConditionCollector * bcs)
{ _bcs = bcs; }



void ElectricalSource::attach_source_to_electrode(const std::string & electrode_label, const std::string & source )
{
  BoundaryCondition *bc = _bcs->get_bc(electrode_label);

  // bc should not be null
  genius_assert(bc);

  // should be an electrode
  genius_assert( bc->is_electrode() );

  // can we find this electrode already in the table?
  BIt bc_source_it = _bc_source_map.find(electrode_label);

  // if not find, insert a new empty item
  if( bc_source_it == _bc_source_map.end() )
  {
    std::vector<VSource *> vsource_list;
    std::vector<ISource *> isource_list;
    _bc_source_map[electrode_label] = std::pair<std::vector<VSource *>, std::vector<ISource *> >(vsource_list, isource_list);
    bc_source_it = _bc_source_map.find(electrode_label);
    genius_assert( bc_source_it != _bc_source_map.end() );
  }

  // clear old value if exist
  bc_source_it->second.first.clear();
  bc_source_it->second.second.clear();

  // insert the source into map if we find
  if( _vsource_list.find(source) != _vsource_list.end() )
    bc_source_it->second.first.push_back( (*_vsource_list.find(source)).second );

  if( _isource_list.find(source) != _isource_list.end() )
    bc_source_it->second.second.push_back( (*_isource_list.find(source)).second );

  // the electrode should only have one kind sources at a time, either vsource or isource
  // so we should check it
  genius_assert( bc_source_it->second.first.empty() || bc_source_it->second.second.empty() );
}


void ElectricalSource::attach_sources_to_electrode(const std::string & electrode_label, const std::vector<std::string> & source_list )
{
  BoundaryCondition *bc = _bcs->get_bc(electrode_label);

  // bc should not be null
  genius_assert(bc);

  // should be an electrode
  genius_assert( bc->is_electrode() );

  // can we find this electrode already in the table?
  BIt bc_source_it = _bc_source_map.find(electrode_label);

  // if not find, insert a new empty item
  if( bc_source_it == _bc_source_map.end() )
  {
    std::vector<VSource *> vsource_list;
    std::vector<ISource *> isource_list;
    _bc_source_map[electrode_label] = std::pair<std::vector<VSource *>, std::vector<ISource *> >(vsource_list, isource_list);
    bc_source_it = _bc_source_map.find(electrode_label);
    genius_assert( bc_source_it != _bc_source_map.end() );
  }

  // clear old value if exist
  bc_source_it->second.first.clear();
  bc_source_it->second.second.clear();

  // insert the source into map if we find
  for(unsigned int i=0; i<source_list.size(); ++i)
  {
    if( _vsource_list.find(source_list[i]) != _vsource_list.end() )
      bc_source_it->second.first.push_back( (*_vsource_list.find(source_list[i])).second );

    if( _isource_list.find(source_list[i]) != _isource_list.end() )
      bc_source_it->second.second.push_back( (*_isource_list.find(source_list[i])).second );
  }

  // the electrode should only have one kind sources at a time, either vsource or isource
  // so we should check it
  genius_assert( bc_source_it->second.first.empty() || bc_source_it->second.second.empty() );
}


void ElectricalSource::attach_voltage_to_electrode(const std::string & electrode_label, double vconst )
{
  // create an VDC with vconst
  std::stringstream ss;
  ss << "##VCONST##" << _counter++ <<"##";
  std::string vlabel = ss.str();
  _vsource_list[vlabel] = new VDC(vlabel, 0, vconst);

  attach_source_to_electrode(electrode_label, vlabel);
}



void ElectricalSource::attach_current_to_electrode(const std::string & electrode_label, double iconst )
{
  // create an IDC with iconst
  std::stringstream ss;
  ss << "##ICONST##" << _counter++ <<"##";
  std::string ilabel = ss.str();
  _isource_list[ilabel] = new IDC(ilabel, 0, iconst);

  attach_source_to_electrode(electrode_label, ilabel);
}



void ElectricalSource::remove_electrode_source(const std::string & electrode_label)
{
  // can we find this electrode already in the table?
  BIt bc_source_it = _bc_source_map.find(electrode_label);

  // if find, clear all the sources
  if( bc_source_it != _bc_source_map.end() )
  {
    // clear old value if exist
    bc_source_it->second.first.clear();
    bc_source_it->second.second.clear();
  }

}



double ElectricalSource::limit_dt(double time, double dt, double dt_min, double v_change, double i_change) const
{
  double dt_limited = dt;
  double dv = 0;
  double di = 0;

  while( dt_limited>dt_min)
  {
    dv = 0;
    di = 0;

    CBIt it = _bc_source_map.begin();
    for(; it!=_bc_source_map.end(); ++it)
    {
      double bc_dv = 0;
      for(unsigned int i=0; i<(*it).second.first.size(); ++i)
        bc_dv += (*it).second.first[i]->dv_max(time, dt_limited, dt_min);

      double bc_di = 0;
      for(unsigned int i=0; i<(*it).second.second.size(); ++i)
        bc_di += (*it).second.second[i]->di_max(time, dt_limited, dt_min);

      if( std::abs(bc_dv) > std::abs(dv) )
        dv = bc_dv;
      if( std::abs(bc_di) > std::abs(di) )
        di = bc_di;
    }


    if(std::abs(dv) > v_change || std::abs(di) > i_change)
      dt_limited *= 0.9;//std::min(v_change/std::abs(dv), i_change/std::abs(di));
    else
      break;
  }

  CBIt it = _bc_source_map.begin();
  for(; it!=_bc_source_map.end(); ++it)
  {
    for(unsigned int i=0; i<(*it).second.first.size(); ++i)
      (*it).second.first[i]->dt_critial_limit(time, dt_limited, dt_min);

    for(unsigned int i=0; i<(*it).second.second.size(); ++i)
      (*it).second.second[i]->dt_critial_limit(time, dt_limited, dt_min);
  }

  return dt_limited;
}



void ElectricalSource::update(double time)
{
  BIt it = _bc_source_map.begin();

  for(; it!=_bc_source_map.end(); ++it)
  {
    // get the bc
    BoundaryCondition *bc = _bcs->get_bc( (*it).first );

    // bc should not be null
    genius_assert(bc);

    // should be an electrode
    genius_assert( bc->is_electrode() );

    // we find voltage source
    if( !(*it).second.first.empty() )
    {
      double vapp = 0;
      for(unsigned int i=0; i<(*it).second.first.size(); ++i)
        vapp += (*it).second.first[i]->vapp(time);

      bc->ext_circuit()->Vapp() = vapp;
      bc->ext_circuit()->set_voltage_driven();
    }

    // or we find current source
    if( !(*it).second.second.empty() )
    {
      double iapp = 0;
      for(unsigned int i=0; i<(*it).second.second.size(); ++i)
        iapp += (*it).second.second[i]->iapp(time);

      bc->ext_circuit()->Iapp() = iapp;
      bc->ext_circuit()->set_current_driven();
    }
  }
}




unsigned int ElectricalSource::steps_by_limiter(double v_change, double i_change, double time) const
{
  double vapp_max = 0;
  double iapp_max = 0;

  CBIt it = _bc_source_map.begin();
  for(; it!=_bc_source_map.end(); ++it)
  {
    BoundaryCondition *bc = _bcs->get_bc( it->first );

    double v_0 = _bc_state_vapp(bc->label());
    double i_0 = _bc_state_iapp(bc->label());

    double v_t = this->vapp(it->first, time);
    double i_t = this->iapp(it->first, time);

    vapp_max = std::max(vapp_max, std::abs(v_t - v_0));
    iapp_max = std::max(iapp_max, std::abs(i_t - i_0));
  }

  return static_cast<unsigned int>(ceil(std::max(1.0, std::max(vapp_max/v_change, iapp_max/i_change))));
}


void ElectricalSource::save_bc_source_state()
{
  _bc_source_state.clear();

  CBIt it = _bc_source_map.begin();
  for(; it!=_bc_source_map.end(); ++it)
  {
    // get the bc
    const BoundaryCondition *bc = _bcs->get_bc( it->first );
    // bc should not be null
    genius_assert(bc);
    // should be an electrode
    genius_assert( bc->is_electrode() );

    switch( bc->ext_circuit()->driven_state() )
    {
      case ExternalCircuit::VDRIVEN : _bc_source_state[ bc->label() ] = std::make_pair( static_cast<int>(ExternalCircuit::VDRIVEN) , bc->ext_circuit()->Vapp() ); break;
      case ExternalCircuit::IDRIVEN : _bc_source_state[ bc->label() ] = std::make_pair( static_cast<int>(ExternalCircuit::IDRIVEN) , bc->ext_circuit()->Iapp() ); break;
    }
  }
}


double ElectricalSource::_bc_state_vapp(const std::string &bc) const
{
  if(_bc_source_state.find(bc) != _bc_source_state.end())
  {
    const std::pair<int, double> & source_state = _bc_source_state.find(bc)->second;
    if( static_cast<ExternalCircuit::DRIVEN>(source_state.first) == ExternalCircuit::VDRIVEN )
    {
      return  source_state.second;
    }
  }
  return 0.0;
}


double ElectricalSource::_bc_state_iapp(const std::string &bc) const
{
  if(_bc_source_state.find(bc) != _bc_source_state.end())
  {
    const std::pair<int, double> & source_state = _bc_source_state.find(bc)->second;
    if( static_cast<ExternalCircuit::DRIVEN>(source_state.first) == ExternalCircuit::IDRIVEN )
    {
      return  source_state.second;
    }
  }
  return 0.0;
}


void ElectricalSource::rampup(double a, double time)
{
  BIt it = _bc_source_map.begin();
  for(; it!=_bc_source_map.end(); ++it)
  {
    // get the bc
    BoundaryCondition *bc = _bcs->get_bc( it->first );

    // bc should not be null
    genius_assert(bc);

    // should be an electrode
    genius_assert( bc->is_electrode() );

    // we find voltage source
    if( !it->second.first.empty() )
    {
      double vapp_0 = _bc_state_vapp(bc->label());
      double vapp = 0;
      for(unsigned int i=0; i<it->second.first.size(); ++i)
        vapp += it->second.first[i]->vapp(time);

      bc->ext_circuit()->Vapp() = vapp_0 + a*(vapp-vapp_0);
      bc->ext_circuit()->set_voltage_driven();
    }

    // or we find current source
    if( !it->second.second.empty() )
    {
      double iapp_0 = _bc_state_iapp(bc->label());
      double iapp = 0;
      for(unsigned int i=0; i<it->second.second.size(); ++i)
        iapp += it->second.second[i]->iapp(time);

      bc->ext_circuit()->Iapp() = iapp_0 + a*(iapp-iapp_0);
      bc->ext_circuit()->set_current_driven();
    }
  }
}




double ElectricalSource::vapp(const std::string &bc, double time) const
{
  CBIt it = _bc_source_map.find(bc);
  double vapp = 0;
  for(unsigned int i=0; i<(*it).second.first.size(); ++i)
    vapp += (*it).second.first[i]->vapp(time);
  return vapp;
}


double ElectricalSource::iapp(const std::string &bc, double time) const
{
  CBIt it = _bc_source_map.find(bc);
  double iapp = 0;
  for(unsigned int i=0; i<(*it).second.second.size(); ++i)
    iapp += (*it).second.second[i]->iapp(time);
  return iapp;
}


void  ElectricalSource::assign_voltage_to(const std::string & electrode_label, double vapp)
{
  // get the bc
  BoundaryCondition *bc = _bcs->get_bc( electrode_label );

  // bc should not be null
  genius_assert( bc );

  // should be an electrode
  genius_assert( bc->is_electrode() );

  bc->ext_circuit()->Vapp() = vapp;

  bc->ext_circuit()->set_voltage_driven();
}


void ElectricalSource::assign_voltage_to(const std::vector<std::string> & electrode_labels, double vapp)
{

  for(unsigned int n=0; n<electrode_labels.size(); ++n)
    assign_voltage_to(electrode_labels[n], vapp);

}



void ElectricalSource::assign_current_to(const std::string & electrode_label, double iapp)
{
  // get the bc
  BoundaryCondition *bc = _bcs->get_bc( electrode_label );

  // bc should not be null
  genius_assert( bc );

  // should be an electrode
  genius_assert( bc->is_electrode() );

  bc->ext_circuit()->Iapp() = iapp;

  bc->ext_circuit()->set_current_driven();
}



void ElectricalSource::assign_current_to(const std::vector<std::string> & electrode_labels, double vapp)
{

  for(unsigned int n=0; n<electrode_labels.size(); ++n)
    assign_current_to(electrode_labels[n], vapp);

}



void  ElectricalSource::SetVDC(const Parser::Card &c)
{

  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  double v = c.get_real("vconst",0.0)*V;
  double Tdelay = c.get_real("tdelay",0.0)*s;

  _vsource_list[label] = new VDC(label, Tdelay, v);
/*
  MESSAGE<<"\nVDC    : "<<label<<"\t"
  <<"td="<<Tdelay/s<<" "
  <<"v="<<v/V<<"\n";
  RECORD();
*/
}



void  ElectricalSource::SetVSIN(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  double v0 = c.get_real("vconst",0.0)*V;
  double vamp = c.get_real("vamp",0.0)*V;
  double Tdelay = c.get_real("tdelay",0.0)*s;
  double freq = c.get_real("freq",0.0)*1.0/s;
  double alpha= c.get_real("alpha",0.0)*1.0/s;

  _vsource_list[label] = new VSIN(label, Tdelay, v0, vamp, freq, alpha);

/*
  MESSAGE<<"\nVSIN   : "<<label<<"\t"
  <<"td="<<Tdelay/s<<" "
  <<"v0="<<v0/V<<" "
  <<"vamp="<<vamp/V<<" "
  <<"freq="<<freq*s<<" "
  <<"alpha="<<alpha*s<<"\n";
  RECORD();
*/
}


void ElectricalSource::SetVEXP(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  double Tdelay = c.get_real("tdelay",0.0)*s;
  double trc = c.get_real("trc",0.0)*s;
  double tfd = c.get_real("tfd",0.0)*s;
  double tfc = c.get_real("tfc",0.0)*s;

  double v1, v2;
  if ( c.is_parameter_exist("v1") && c.is_parameter_exist("v2"))
  {
    v1 = c.get_real("v1",0.0)*V;
    v2 = c.get_real("v2",0.0)*V;
  }
  else
  {
    v1 = c.get_real("vlo",0.0)*V;
    v2 = c.get_real("vhi",0.0)*V;
  }

  _vsource_list[label] = new VEXP(label, Tdelay, v1, v2, trc, tfd, tfc);

/*
  MESSAGE<<"\nVEXP   : "<<label<<"\t"
  <<"td="<<Tdelay/s<<" "
  <<"trc="<<trc/s<<" "
  <<"tfd="<<tfd/s<<" "
  <<"tfc="<<tfc/s<<"\n\t\t"
  <<"vlo="<<vlo/V<<" "
  <<"vhi="<<vhi/V<<"\n";
  RECORD();
*/
}



void  ElectricalSource::SetVPULSE(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  double Tdelay = c.get_real("tdelay",0.0)*s;
  double tr = c.get_real("tr",1e-9)*s;
  double tf = c.get_real("tf",1e-9)*s;
  double pw = c.get_real("pw",5e-7)*s;
  double pr = c.get_real("pr",1e-6)*s;

  double v1, v2;
  if ( c.is_parameter_exist("v1") && c.is_parameter_exist("v2"))
  {
    v1 = c.get_real("v1",0.0)*V;
    v2 = c.get_real("v2",1.0)*V;
  }
  else
  {
    v1 = c.get_real("vlo",0.0)*V;
    v2 = c.get_real("vhi",1.0)*V;
  }

  _vsource_list[label] = new VPULSE(label, Tdelay, v1, v2, tr, tf, pw, pr);

/*
  MESSAGE<<"\nVPULSE : "<<label<<"\t"
  <<"td="<<Tdelay/s<<" "
  <<"tr="<<tr/s<<" "
  <<"tf="<<tf/s<<" "
  <<"pw="<<pw/s<<" "
  <<"pr="<<pr/s<<"\n\t\t"
  <<"vlo="<<vlo/V<<" "
  <<"vhi="<<vhi/V<<"\n";
  RECORD();
*/
}



void  ElectricalSource::SetVUSER(const Parser::Card &c)
{

  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  std::string expr = c.get_string("expression","1.0*V");

  _vsource_list[label] = new VUSER(label, expr);
}




void  ElectricalSource::SetVSHELL(const Parser::Card &c)
{

  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  std::string filename = c.get_string("dll","");
  std::string funcname = c.get_string("func","");

#ifdef WINDOWS

  HINSTANCE hInstLibrary = LoadLibrary(filename.c_str());

  void *fp = GetProcAddress(hInstLibrary, funcname.c_str());

  _vsource_list[label] = new VSHELL(label, hInstLibrary, fp, s, V);

#else

#ifdef RTLD_DEEPBIND
  void * dp = dlopen(filename.c_str(), RTLD_LAZY|RTLD_DEEPBIND);
#else
  void * dp = dlopen(filename.c_str(), RTLD_LAZY);
#endif
  genius_assert(dp);

  void *fp = dlsym(dp,funcname.c_str());
  genius_assert(fp);

  _vsource_list[label] = new VSHELL(label, dp, fp, s, V);

#endif

/*
  MESSAGE<<"\nVSHELL : "<<label<<" load from "<<filename<<"\n";
  RECORD();
*/



}





void  ElectricalSource::SetIDC(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  double i = c.get_real("iconst",0.0)*A;
  double Tdelay = c.get_real("tdelay",0.0)*s;

  _isource_list[label] = new IDC(label, Tdelay, i);

/*
  MESSAGE<<"\nIDC    : "<<label<<"\t"
  <<"td="<<Tdelay/s<<" "
  <<"i="<<i/A<<"\n";
  RECORD();
*/
}



void  ElectricalSource::SetISIN(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  double i0 = c.get_real("iconst",0.0)*A;
  double iamp = c.get_real("iamp",0.0)*A;
  double Tdelay = c.get_real("tdelay",0.0)*s;
  double freq = c.get_real("freq",0.0)*1.0/s;
  double alpha= c.get_real("alpha",0.0)*1.0/s;

  _isource_list[label] = new ISIN(label, Tdelay, i0, iamp, freq, alpha);

/*
  MESSAGE<<"\nISIN   : "<<label<<"\t"
  <<"td="<<Tdelay/s<<" "
  <<"iamp="<<iamp/A<<" "
  <<"freq="<<freq*s<<" "
  <<"alpha="<<alpha*s<<"\n";
  RECORD();
*/
}



void  ElectricalSource::SetIEXP(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  double Tdelay = c.get_real("tdelay",0.0)*s;
  double trc = c.get_real("trc",0.0)*s;
  double tfd = c.get_real("tfd",0.0)*s;
  double tfc = c.get_real("tfc",0.0)*s;

  double i1, i2;
  if ( c.is_parameter_exist("i1") && c.is_parameter_exist("i2"))
  {
    i1 = c.get_real("i1",0.0)*A;
    i2 = c.get_real("i2",0.0)*A;
  }
  else
  {
    i1 = c.get_real("ilo",0.0)*A;
    i2 = c.get_real("ihi",0.0)*A;
  }

  _isource_list[label] = new IEXP(label, Tdelay, i1, i2, trc, tfd, tfc);

/*
  MESSAGE<<"\nIEXP   : "<<label<<"\t"
  <<"td="<<Tdelay/s<<" "
  <<"trc="<<trc/s<<" "
  <<"tfd="<<tfd/s<<" "
  <<"tfc="<<tfc/s<<"\n\t\t"
  <<"ilo="<<ilo/A<<" "
  <<"ihi="<<ihi/A<<"\n";
  RECORD();
*/
}



void  ElectricalSource::SetIPULSE(const Parser::Card &c)
{
  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  double Tdelay = c.get_real("tdelay",0.0)*s;
  double tr = c.get_real("tr",1e-9)*s;
  double tf = c.get_real("tf",1e-9)*s;
  double pw = c.get_real("pw",5e-7)*s;
  double pr = c.get_real("pr",1e-6)*s;

  double i1, i2;
  if ( c.is_parameter_exist("v1") && c.is_parameter_exist("v2"))
  {
    i1 = c.get_real("i1",0.0)*A;
    i2 = c.get_real("i2",1.0)*A;
  }
  else
  {
    i1 = c.get_real("ilo",0.0)*A;
    i2 = c.get_real("ihi",1.0)*A;
  }

  _isource_list[label] = new IPULSE(label, Tdelay, i1, i2, tr, tf, pw, pr);

/*
  MESSAGE<<"\nIPULSE : "<<label<<"\t"
  <<"td="<<Tdelay/s<<" "
  <<"tr="<<tr/s<<" "
  <<"tf="<<tf/s<<" "
  <<"pw="<<pw/s<<" "
  <<"pr="<<pr/s<<"\n\t\t"
  <<"ilo="<<ilo/A<<" "
  <<"ihi="<<ihi/A<<"\n";
  RECORD();
*/
}


void  ElectricalSource::SetIUSER(const Parser::Card &c)
{

  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  std::string expr = c.get_string("expression","1.0*mA");

  _isource_list[label] = new IUSER(label, expr);
}



void  ElectricalSource::SetISHELL(const Parser::Card &c)
{

  std::string label = c.get_string("id","");
  genius_assert( label!="" );

  std::string filename = c.get_string("dll","");
  std::string funcname = c.get_string("func","");

#ifdef WINDOWS

  HINSTANCE hInstLibrary = LoadLibrary(filename.c_str());

  void *fp = GetProcAddress(hInstLibrary, funcname.c_str());

  _isource_list[label] = new ISHELL(label, hInstLibrary, fp, s, A);

#else
#ifdef RTLD_DEEPBIND
  void * dp = dlopen(filename.c_str(), RTLD_LAZY|RTLD_DEEPBIND);
#else
  void * dp = dlopen(filename.c_str(), RTLD_LAZY);
#endif
  genius_assert(dp);


  void *fp = dlsym(dp,funcname.c_str());
  genius_assert(fp);

  _isource_list[label] = new ISHELL(label, dp, fp, s, A);

#endif



/*
  MESSAGE<<"\nISHELL : "<<label<<" load from "<<filename<<"\n";
  RECORD();
*/

}




