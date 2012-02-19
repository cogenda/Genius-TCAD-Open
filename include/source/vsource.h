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

//  $Id: vsource.h,v 1.3 2008/07/09 05:58:16 gdiso Exp $


//sevel types of voltage sources are defined here
//reference: spice voltage source model

#ifndef __vsource_h__
#define __vsource_h__

#include <string>
#include <cmath>

#include "config.h"
#include "expr_evaluate.h"

#ifdef WINDOWS
  class HINSTANCE__; // Forward or never
  typedef HINSTANCE__* HINSTANCE;
#endif




/**
 * virtual base class for voltage source
 */
class VSource
{
private:

  /**
   * the unique label for a voltage source
   */
  std::string _label;

public:

  /**
   * constructor
   */
  VSource(const std::string & s): _label(s) {}

  /**
   * have nothing to do
   */
  virtual ~VSource(){}

  /**
   * virtual function, @return vapp for a given time step
   */
  virtual double vapp(double t)=0;

  /**
   * @return const reference of label
   */
  const std::string & label() const
  { return _label;}

  /**
   * @return writable reference to label
   */
  std::string & label()
  { return _label;}

};



/**
 * Direct Current source
 */
class VDC: public VSource
{
private:

  /**
   * time delay
   */
  double td;

  /**
   * voltage amplitude
   */
  double Vdc;

public:

  /**
   * constructor
   */
  VDC(const std::string & s, double t1,double v1)
  : VSource(s), td(t1),Vdc(v1)
  {}

  /**
   * destructor have nothing to do
   */
  ~VDC(){}

  /**
   * @return constant vapp when t>td, else 0
   */
  double vapp(double t)
  { return t>=td? Vdc:0;}

};



/**
 * sine wave voltage source
 */
class VSIN: public VSource
{
private:
  /**
   * time delay
   */
  double td;

  /**
   * bias dc voltage
   */
  double V0;

  /**
   * sine wave amplitude
   */
  double Vamp;

  /**
   * sine wave frequency
   */
  double fre;

  /**
   * attenuation parameter, should be >=0
   */
  double alpha;

public:
  /**
   * constructor
   */
  VSIN(const std::string & s, double t1,double v0,double v1,double f1,double a1=0)
  : VSource(s), td(t1),V0(v0),Vamp(v1),fre(f1),alpha(a1)
  {}

  /**
   * destructor have nothing to do
   */
  ~VSIN(){}

  /**
   * @return voltage with sine wave when t>td, else return bias dc voltage
   */
  double vapp(double t)
  { return t>=td ? V0+Vamp*exp(-alpha*(t-td))*sin(2*3.14159265358979323846*fre*(t-td)) : V0; }
};



/**
 * pulsed voltage source
 */
class VPULSE: public VSource
{
private:
  /**
   * time delay
   */
  double td;

  /**
   * time for raising edge
   */
  double tr;

  /**
   * time for falling edge
   */
  double tf;

  /**
   * pulse width
   */
  double pw;

  /**
   * pulse repeat period
   */
  double pr;

  /**
   * the initial voltage level
   */
  double V1;

  /**
   * the pulsed voltage level
   */
  double V2;

public:

  /**
   * constructor
   */
  VPULSE(const std::string & s, double t1,double v1,double v2,double t2,double t3,double t4, double t5)
  :VSource(s),td(t1),tr(t2),tf(t3),pw(t4),pr(t5),V1(v1),V2(v2)
  {}

  /**
   * destructor have nothing to do
   */
  ~VPULSE(){}

  /**
   * @return waveform of pulsed source
   */
  double vapp(double t)
  {
    if(t<td)
      return V1;
    else
    {
      t-=td;
      while(t>pr) t-=pr;
      if(t<tr)
        return V1+t*(V2-V1)/tr;
      else if(t<tr+pw)
        return V2;
      else if(t<tr+pw+tf)
        return V2-(t-tr-pw)*(V2-V1)/tf;
      else    return V1;
    }
  }

};



/**
 * dual exponential source
 */
class VEXP: public VSource
{
private:
  /**
   * time delay
   */
  double td;

  /**
   * time constant for raising edge
   */
  double trc;

  /**
   * time delay for falling
   */
  double tfd;

  /**
   * time constant for falling edge
   */
  double tfc;

  /**
   * the initial voltage level
   */
  double V1;

  /**
   * the pulsed voltage level
   */
  double V2;

public:
  /**
   * constructor
   */
  VEXP(const std::string & s, double t1,double v1,double v2,double t2,double t3,double t4)
  :VSource(s), td(t1),trc(t2),tfd(t3),tfc(t4),V1(v1),V2(v2)
  {}

  /**
   * destructor have nothing to do
   */
  ~VEXP(){}

  /**
   * @return waveform of exponential source
   */
  double vapp(double t)
  {
    if(t<=td)
      return V1;
    else if(t<=tfd)
      return V1+(V2-V1)*(1-exp(-(t-td)/trc));
    else
      return V1+(V2-V1)*(1-exp(-(t-td)/trc))+(V1-V2)*(1-exp(-(t-tfd)/tfc));

  }

};


/**
 * user defind source from analytic express
 */
class VUSER : public VSource
{
private:

  ExprEvalute expr_evaluator;

public:

 VUSER(const std::string &s, const std::string & expr) : VSource(s), expr_evaluator(expr)
 {}

 /**
  * @return waveform of exponential source
  */
  double vapp(double t)
  {
    return expr_evaluator(0,0,0,t);
  }

};



/**
 * A shell for user provide source from dll file
 */
class VSHELL: public VSource
{
private:
  /**
   * the pointer to dll file
   */
#ifdef WINDOWS
   HINSTANCE                  dll;
#else
   void                      *dll;
#endif

  /**
   * the pointer to function in the dll file
   */
  double (*Vapp_Shell)(double);

  /**
   * time scale parameter
   */
  double scale_t;

  /**
   * voltage scale parameter
   */
  double scale_V;

public:

  /**
   * constructor, hold the pointer to dll file and function
   */
#ifdef WINDOWS
  VSHELL(const std::string & s, HINSTANCE dp, void * fp, double s_t, double s_V);
#else
  VSHELL(const std::string & s, void * dp, void * fp, double s_t, double s_V);
#endif


  /**
   * destructor, free the pointer
   */
  ~VSHELL();

  /**
   * call Vapp_Shell to get the user provide value
   */
  double vapp(double t)
  {
     return scale_V*Vapp_Shell(t/scale_t);
  }

};



#endif //#ifndef __vsource_h__

