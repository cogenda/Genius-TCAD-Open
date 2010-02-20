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


#include "electrical_source.h"
#include "simulation_system.h"
#include "boundary_condition_collector.h"

// for short
using  PhysicalUnit::s;
using  PhysicalUnit::V;
using  PhysicalUnit::A;


ElectricalSource::ElectricalSource(Parser::InputParser & decks)
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
    }

  }

  MESSAGE<<"done.\n"<<std::endl;  RECORD();

}



ElectricalSource::~ElectricalSource()
{
  std::map<const std::string, VSource * >::iterator it1 = _vsource_list.begin();
  for(; it1!=_vsource_list.end(); ++it1)
    delete (*it1).second;
  _vsource_list.clear();

  std::map<const std::string, ISource * >::iterator it2 = _isource_list.begin();
  for(; it2!=_isource_list.end(); ++it2)
    delete (*it2).second;
  _isource_list.clear();

}



void ElectricalSource::link_to_bcs(BoundaryConditionCollector * bcs)
{ _bcs = bcs; }



void ElectricalSource::attach_sources_to_electrode(const std::string & electrode_label, const std::vector<std::string> & source_list )
{
  BoundaryCondition *bc = _bcs->get_bc(electrode_label);

  // bc should not be null
  genius_assert(bc);

  // should be an electrode
  genius_assert( bc->is_electrode() );

  // can we find this electrode already in the table?
  BIt it = _bc_source_map.find(electrode_label);

  // if not find, insert a new empty item
  if( it == _bc_source_map.end() )
  {
    std::vector<VSource *> vsource_list;
    std::vector<ISource *> isource_list;
    _bc_source_map[electrode_label] = std::pair<std::vector<VSource *>, std::vector<ISource *> >(vsource_list, isource_list);
    it = _bc_source_map.find(electrode_label);
    genius_assert( it != _bc_source_map.end() );
  }

  // clear old value if exist
  (*it).second.first.clear();
  (*it).second.second.clear();

  // insert the source into map if we find
  for(unsigned int i=0; i<source_list.size(); ++i)
  {
    if( _vsource_list.find(source_list[i]) != _vsource_list.end() )
      (*it).second.first.push_back( (*_vsource_list.find(source_list[i])).second );

    if( _isource_list.find(source_list[i]) != _isource_list.end() )
      (*it).second.second.push_back( (*_isource_list.find(source_list[i])).second );
  }

  // the electrode should only have one kind sources at a time, either vsource or isource
  // so we should check it
  genius_assert( (*it).second.first.empty() || (*it).second.second.empty() );
}


void ElectricalSource::attach_voltage_to_electrode(const std::string & electrode_label, PetscScalar vconst )
{
  BoundaryCondition *bc = _bcs->get_bc(electrode_label);

  // bc should not be null
  genius_assert(bc);

  // should be an electrode
  genius_assert( bc->is_electrode() );

  // remove old sources
  remove_electrode_source(electrode_label);

  assign_voltage_to(electrode_label, vconst);
}



void ElectricalSource::attach_current_to_electrode(const std::string & electrode_label, PetscScalar iconst )
{
  BoundaryCondition *bc = _bcs->get_bc(electrode_label);

  // bc should not be null
  genius_assert(bc);

  // should be an electrode
  genius_assert( bc->is_electrode() );

  // remove old sources
  remove_electrode_source(electrode_label);

  assign_current_to(electrode_label, iconst);
}



void ElectricalSource::remove_electrode_source(const std::string & electrode_label)
{
  // can we find this electrode already in the table?
  BIt it = _bc_source_map.find(electrode_label);

  // if find, clear all the sources
  if( it != _bc_source_map.end() )
  {
    // clear old value if exist
    (*it).second.first.clear();
    (*it).second.second.clear();
  }

}




void ElectricalSource::update(PetscScalar time)
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



void  ElectricalSource::assign_voltage_to(const std::string & electrode_label, PetscScalar vapp)
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


void ElectricalSource::assign_voltage_to(const std::vector<std::string> & electrode_labels, PetscScalar vapp)
{

  for(unsigned int n=0; n<electrode_labels.size(); ++n)
    assign_voltage_to(electrode_labels[n], vapp);

}



void ElectricalSource::assign_current_to(const std::string & electrode_label, PetscScalar iapp)
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



void ElectricalSource::assign_current_to(const std::vector<std::string> & electrode_labels, PetscScalar vapp)
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

  double tr = c.get_real("tr",1e-12)*s;

  double tf = c.get_real("tf",1e-12)*s;

  double pw = c.get_real("pw",0.0)*s;

  double pr = c.get_real("pr",0.0)*s;

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

#ifdef CYGWIN

  HINSTANCE hInstLibrary = LoadLibrary(filename.c_str());

  void *fp = GetProcAddress(hInstLibrary, funcname.c_str());

  _vsource_list[label] = new VSHELL(label, hInstLibrary, fp, s, V);

#else

  void * dp = dlopen(filename.c_str(), RTLD_LAZY);
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
    i1 = c.get_real("i1",0.0)*V;
    i2 = c.get_real("i2",0.0)*V;
  }
  else
  {
    i1 = c.get_real("ilo",0.0)*V;
    i2 = c.get_real("ihi",0.0)*V;
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

  double tr = c.get_real("tr",1e-12)*s;

  double tf = c.get_real("tf",1e-12)*s;

  double pw = c.get_real("pw",0.0)*s;

  double pr = c.get_real("pr",0.0)*s;

  double i1, i2;
  if ( c.is_parameter_exist("v1") && c.is_parameter_exist("v2"))
  {
    i1 = c.get_real("i1",0.0)*V;
    i2 = c.get_real("i2",0.0)*V;
  }
  else
  {
    i1 = c.get_real("ilo",0.0)*V;
    i2 = c.get_real("ihi",0.0)*V;
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

#ifdef CYGWIN

  HINSTANCE hInstLibrary = LoadLibrary(filename.c_str());

  void *fp = GetProcAddress(hInstLibrary, funcname.c_str());

  _isource_list[label] = new ISHELL(label, hInstLibrary, fp, s, A);

#else
  void * dp = dlopen(filename.c_str(), RTLD_LAZY);
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




