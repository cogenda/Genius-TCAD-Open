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
#include <map>
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
  virtual double vapp(const double t) const=0;

  /**
   * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
   */
  virtual double dv_max(const double t, const double dt, const double dt_min) const
  {
    double v = vapp(t);
    double dv = vapp(t+dt) - v;
    for(double clock=t; clock < t+dt; clock+=dt_min)
      if( std::abs(dv) < std::abs(vapp(clock) - v))
        dv = vapp(clock) - v;
    return dv;
  }


  /**
   * limit the dt by wave form
   */
  virtual void dt_critial_limit(const double t, double & dt, const double dt_min)  const {}

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
  double vapp(const double t) const
  { return t>=td? Vdc:0;}


  /**
   * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
   */
  double dv_max(const double t, const double dt, const double dt_min) const
  {
    if(t<td && t+dt >=td) return Vdc-0.0;
    return 0.0;
  }

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
  double vapp(const double t) const
  { return t>=td ? V0+Vamp*exp(-alpha*(t-td))*sin(2*3.14159265358979323846*fre*(t-td)) : V0; }


  /**
   * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
   */
  double dv_max(const double t, const double dt, const double dt_min) const
  {
    int cycles = fre*(t-td);
    double t1 = t;
    double t2 = t+dt;
    double ta1 = cycles/fre + 0.25/fre + td;
    double ta2 = cycles/fre + 0.75/fre + td;
    double ta3 = (cycles+1)/fre + 0.25/fre + td;
    double ta4 = (cycles+1)/fre + 0.75/fre + td;

    std::map<double, double> v_wave;
    v_wave[t1 ] =  vapp(t1);
    v_wave[t2 ] =  vapp(t2);
    v_wave[ta1] =  vapp(ta1);
    v_wave[ta2] =  vapp(ta2);
    v_wave[ta3] =  vapp(ta3);
    v_wave[ta4] =  vapp(ta4);

    std::map<double, double>::iterator it1 = v_wave.find(t1);
    std::map<double, double>::iterator it2 = v_wave.find(t2);

    double v1 = it1->second;
    double dv = it2->second - v1;
    while(it1++ != it2)
    {
      if( std::abs(dv) < std::abs(it1->second - v1) )
        dv = it1->second - v1;
    }
    return dv;
  }


  void dt_critial_limit(const double t, double & dt, const double dt_min) const
  {
    int cycles = fre*(t-td);

    double ta1 = cycles/fre + 0.25/fre + td;
    double ta2 = cycles/fre + 0.75/fre + td;

    if( t<ta1 && t+dt>ta1 && ta1-t>=dt_min ) { dt = ta1-t; return; }
    if( t<ta2 && t+dt>ta2 && ta2-t>=dt_min ) { dt = ta2-t; return; }
  }


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
  double vapp(const double t) const
  {
    double _t = t;
    if(_t<td)
      return V1;
    else
    {
      _t-=td;
      while(_t>pr) _t-=pr;
      if(_t<tr)
        return V1+_t*(V2-V1)/tr;
      else if(_t<tr+pw)
        return V2;
      else if(_t<tr+pw+tf)
        return V2-(_t-tr-pw)*(V2-V1)/tf;
      else    return V1;
    }
  }


  /**
   * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
   */
  double dv_max(const double t, const double dt, const double dt_min) const
  {
    int cycles = (t-td)/pr;

    double ta1 = 0  + cycles*pr + td;
    double ta2 = tr + cycles*pr + td;
    double ta3 = tr + pw + cycles*pr + td;
    double ta4 = tr + pw + tf + cycles*pr + td;
    double ta5 = 0  + (cycles+1)*pr + td;
    double ta6 = tr + (cycles+1)*pr + td;
    double ta7 = tr + pw + (cycles+1)*pr + td;
    double ta8 = tr + pw + tf + (cycles+1)*pr + td;

    double t1  = t;
    double t2  = t+dt;

    std::map<double, double> v_wave;
    v_wave[t1 ] =  vapp(t1);
    v_wave[t2 ] =  vapp(t2);
    v_wave[ta1] =  vapp(ta1);
    v_wave[ta2] =  vapp(ta2);
    v_wave[ta3] =  vapp(ta3);
    v_wave[ta4] =  vapp(ta4);
    v_wave[ta5] =  vapp(ta5);
    v_wave[ta6] =  vapp(ta6);
    v_wave[ta7] =  vapp(ta7);
    v_wave[ta8] =  vapp(ta8);

    std::map<double, double>::iterator it1 = v_wave.find(t1);
    std::map<double, double>::iterator it2 = v_wave.find(t2);

    double v1 = it1->second;
    double dv = it2->second - v1;
    while(it1++ != it2)
    {
      if( std::abs(dv) < std::abs(it1->second - v1) )
        dv = it1->second - v1;
    }
    return dv;
  }


  void dt_critial_limit(const double t, double & dt, const double dt_min) const
  {
    int cycles = (t-td)/pr;

    double ta2 = tr + cycles*pr + td;
    double ta3 = tr + pw + cycles*pr + td;
    double ta4 = tr + pw + tf + cycles*pr + td;
    double ta5 = 0  + (cycles+1)*pr + td;

    if( t<ta2 && t+dt>ta2 && ta2-t>=dt_min ) { dt = ta2-t; return; }
    if( t<ta3 && t+dt>ta3 && ta3-t>=dt_min ) { dt = ta3-t; return; }
    if( t<ta4 && t+dt>ta4 && ta4-t>=dt_min ) { dt = ta4-t; return; }
    if( t<ta5 && t+dt>ta5 && ta5-t>=dt_min ) { dt = ta5-t; return; }

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
  double vapp(const double t) const
  {
    if(t<=td)
      return V1;
    else if(t<=tfd)
      return V1+(V2-V1)*(1-exp(-(t-td)/trc));
    else
      return V1+(V2-V1)*(1-exp(-(t-td)/trc))+(V1-V2)*(1-exp(-(t-tfd)/tfc));

  }


  /**
   * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
   */
  double dv_max(const double t, const double dt, const double dt_min) const
  {
    double t1 = (t-td);
    double t2 = (t+dt-td);
    double ta1 = tfd;

    std::map<double, double> v_wave;
    v_wave[t1 ] =  vapp(t1);
    v_wave[t2 ] =  vapp(t2);
    v_wave[ta1] =  vapp(ta1);

    std::map<double, double>::iterator it1 = v_wave.find(t1);
    std::map<double, double>::iterator it2 = v_wave.find(t2);

    double v1 = it1->second;
    double dv = it2->second - v1;
    while(it1++ != it2)
    {
      if( std::abs(dv) < std::abs(it1->second - v1) )
        dv = it1->second - v1;
    }
    return dv;
  }


  void dt_critial_limit(const double t, double & dt, const double dt_min) const
  {
    double ta1 = td;
    double ta2 = tfd;

    if( t<ta1 && t+dt>ta1 && ta1-t>=dt_min ) { dt = ta1-t; return; }
    if( t<ta2 && t+dt>ta2 && ta2-t>=dt_min ) { dt = ta2-t; return; }
  }

};



/**
 * single gaussion pulse
 */
class VGAUSSION : public VSource
{
private:

  double _vamp;
  double _t0;
  double _tau;

public:

  VGAUSSION(const std::string &s, double vamp, double t0, double tau) : VSource(s), _vamp(vamp), _t0(t0), _tau(tau)
  {}

 /**
  * @return waveform of exponential source
  */
  double vapp(const double t) const
  {
    double q = (t-_t0)*(t-_t0)/(_tau*_tau);
    if(q > 30) return 0.0;
    return _vamp*exp(-q);
  }

 /**
  * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
  */
  double dv_max(const double t, const double dt, const double dt_min) const
  {
    if(t+dt < _t0) return  vapp(t+dt) - vapp(t);
    if(t    > _t0) return  vapp(t) - vapp(t+dt);
    return std::max( _vamp - vapp(t), _vamp - vapp(t+dt) );
  }
  
  void dt_critial_limit(const double t, double & dt, const double dt_min) const
  {
    if( t<_t0 && t+dt>_t0 && _t0-t>=dt_min ) { dt = _t0-t; return; }
  }
  
};


/**
 * double expotential pulse, as a description of HEMP pulse
 */
class VDEXP : public VSource
{
private:

  double _vamp;
  double _h;
  double _alpha;
  double _beta;
  double _tp;

public:

  VDEXP(const std::string &s, double vamp, double h, double alpha, double beta) : VSource(s), _vamp(vamp), _h(h), _alpha(alpha), _beta(beta)
  {
    _tp = (log(_alpha)-log(_beta)) / (_alpha - _beta);
  }

 /**
  * @return waveform of exponential source
  */
  double vapp(const double t) const
  {
    double q1 = _alpha*t;
    double q2 = _beta*t;
    if(q1 > 30 && q2 > 30) return 0.0;
    return _h*_vamp*(exp(-q1) - exp(-q2));
  }

 /**
  * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
  */
  double dv_max(const double t, const double dt, const double dt_min) const
  {
    if(t+dt < _tp) return  vapp(t+dt) - vapp(t);
    if(t    > _tp) return  vapp(t) - vapp(t+dt);
    return std::max( _vamp - vapp(t), _vamp - vapp(t+dt) );
  }
};


/**
 * user defind source from analytic express
 */
class VUSER : public VSource
{
private:

  mutable ExprEvalute expr_evaluator;

public:

 VUSER(const std::string &s, const std::string & expr) : VSource(s), expr_evaluator(expr)
 {}

 /**
  * @return waveform of exponential source
  */
  double vapp(const double t) const
  {
    return expr_evaluator(0,0,0,t);
  }

 /**
  * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
  */
  virtual double dv_max(const double t, const double dt, const double dt_min) const
  {
    return 0.0;
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
  double vapp(double t) const
  {
     return scale_V*Vapp_Shell(t/scale_t);
  }


 /**
  * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
  */
  virtual double dv_max(const double t, const double dt, const double dt_min) const
  {
    return 0.0;
  }

};


/**
 * Voltage mixer Current source
 */
class VMIX: public VSource
{
private:

  const VSource * _v1;

  const VSource * _v2;

public:

  /**
   * constructor
   */
  VMIX(const std::string & s, const VSource *v1, const VSource *v2)
   : VSource(s), _v1(v1), _v2(v2)
  {}

  /**
   * destructor have nothing to do
   */
  ~VMIX()  {}

  /**
   * @return constant vapp when t>td, else 0
   */
  double vapp(const double t) const
  { return _v1->vapp(t) * _v2->vapp(t) ;}


  /**
   * @return the max value of vapp(t+delta_t)-vapp(t), delta_t in [0, dt]
   */
  double dv_max(const double t, const double dt, const double dt_min) const
  {
    double dv_max1 = _v1->dv_max(t, dt, dt_min);
    double dv_max2 = _v2->dv_max(t, dt, dt_min);
    return dv_max1*dv_max2;
  }

};



#endif //#ifndef __vsource_h__

