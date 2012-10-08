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

//  $Id: isource.h,v 1.3 2008/07/09 05:58:16 gdiso Exp $


//sevel types of current sources are defined here
//reference: spice current source model


#ifndef __isource_h__
#define __isource_h__

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
 * virtual base class for current source
 */
class ISource
{
private:

  /**
   * the unique label for a current source
   */
  std::string _label;

public:

  /**
   * constructor
   */
  ISource(const std::string & s): _label(s) {}

  /**
   * have nothing to do
   */
  virtual ~ISource(){}

  /**
   * virtual function, @return current for a given time step
   */
  virtual double iapp(double t)=0;

  /**
   * @return the max value of iapp(t+delta_t)-iapp(t), delta_t in [0, dt]
   */
  virtual double di_max(double t, double dt, double dt_min)
  {
    double i = iapp(t);
    double di = iapp(t+dt) - i;
    for(double clock=t; clock < t+dt; clock+=dt_min)
      if( std::abs(di) < std::abs(iapp(clock) - i))
        di = iapp(clock) - i;
    return di;
  }

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
class IDC: public ISource
{
private:
  /**
   * time delay
   */
  double td;

  /**
   * constant current amplitude
   */
  double Idc;

public:
  /**
   * constructor
   */
  IDC(const std::string & s, double t1, double I1)
  :ISource(s), td(t1),Idc(I1)
  {}

  /**
   * destructor have nothing to do
   */
  ~IDC(){}

  /**
   * @return constant current when t>td, else 0
   */
  double iapp(double t)
  { return t>=td? Idc:0;}


  /**
   * @return the max value of iapp(t+delta_t)-iapp(t), delta_t in [0, dt]
   */
  double di_max(double t, double dt, double dt_min)
  {
    if(t<td && t+dt >=td) return Idc-0.0;
    return 0.0;
  }
};


/**
 * sine wave current source
 */
class ISIN: public ISource
{
private:
  /**
   * time delay
   */
  double td;

  /**
   * bias dc current
   */
  double I0;

  /**
   * sine wave amplitude
   */
  double Iamp;

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
  ISIN(const std::string & s, double t1, double I, double Ia, double f1,double a1=0)
  : ISource(s),td(t1),I0(I),Iamp(Ia),fre(f1),alpha(a1)
  {}

  /**
   * destructor have nothing to do
   */
  ~ISIN(){}

  /**
   * @return current with sine wave when t>td, else return dc current of I0
   */
  double iapp(double t)
  { return t>=td? Iamp*exp(-alpha*(t-td))*sin(2*3.14159265358979323846*fre*(t-td)):0;}


  /**
   * @return the max value of iapp(t+delta_t)-iapp(t), delta_t in [0, dt]
   */
  double di_max(double t, double dt, double dt_min)
  {
    int cycles = fre*(t-td);
    double t1 = t;
    double t2 = t+dt;
    double ta1 = cycles/fre + 0.25/fre + td;
    double ta2 = cycles/fre + 0.75/fre + td;
    double ta3 = (cycles+1)/fre + 0.25/fre + td;
    double ta4 = (cycles+1)/fre + 0.75/fre + td;

    std::map<double, double> i_wave;
    i_wave[t1 ] =  iapp(t1);
    i_wave[t2 ] =  iapp(t2);
    i_wave[ta1] =  iapp(ta1);
    i_wave[ta2] =  iapp(ta2);
    i_wave[ta3] =  iapp(ta3);
    i_wave[ta4] =  iapp(ta4);

    std::map<double, double>::iterator it1 = i_wave.find(t1);
    std::map<double, double>::iterator it2 = i_wave.find(t2);

    double i1 = it1->second;
    double di = it2->second - i1;
    while(it1++ != it2)
    {
      if( std::abs(di) < std::abs(it1->second - i1) )
        di = it1->second - i1;
    }
    return di;
  }
};



/**
 * pulsed current source
 */
class IPULSE: public ISource
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
   * the initial current level
   */
  double I1;

  /**
   * the pulsed current level
   */
  double I2;

public:
  /**
   * constructor
   */
  IPULSE(const std::string & s, double t1,double i1,double i2,double t2,double t3,double t4, double t5)
  :ISource(s),td(t1),tr(t2),tf(t3),pw(t4),pr(t5),I1(i1),I2(i2)
  {}

  /**
   * destructor have nothing to do
   */
  ~IPULSE(){}


  /**
   * @return waveform of pulsed source
   */
  double iapp(double t)
  {
    if(t<td)
      return I1;
    else
    {
      t-=td;
      while(t>pr) t-=pr;
      if(t<tr)
        return I1+t*(I2-I1)/tr;
      else if(t<tr+pw)
        return I2;
      else if(t<tr+pw+tf)
        return I2-(t-tr-pw)*(I2-I1)/tf;
      else    return I1;
    }
  }

  /**
   * @return the max value of iapp(t+delta_t)-iapp(t), delta_t in [0, dt]
   */
  double di_max(double t, double dt, double dt_min)
  {
    int cycles = (t-td)/pr;
    double t1  = t;
    double t2  = t+dt;
    double ta1 = 0  + cycles*pr + td;
    double ta2 = tr + cycles*pr + td;
    double ta3 = tr + pw + cycles*pr + td;
    double ta4 = tr + pw + tf + cycles*pr + td;
    double ta5 = 0  + (cycles+1)*pr + td;
    double ta6 = tr + (cycles+1)*pr + td;
    double ta7 = tr + pw + (cycles+1)*pr + td;
    double ta8 = tr + pw + tf + (cycles+1)*pr + td;

    std::map<double, double> i_wave;
    i_wave[t1 ] =  iapp(t1);
    i_wave[t2 ] =  iapp(t2);
    i_wave[ta1] =  iapp(ta1);
    i_wave[ta2] =  iapp(ta2);
    i_wave[ta3] =  iapp(ta3);
    i_wave[ta4] =  iapp(ta4);
    i_wave[ta5] =  iapp(ta5);
    i_wave[ta6] =  iapp(ta6);
    i_wave[ta7] =  iapp(ta7);
    i_wave[ta8] =  iapp(ta8);

    std::map<double, double>::iterator it1 = i_wave.find(t1);
    std::map<double, double>::iterator it2 = i_wave.find(t2);

    double i1 = it1->second;
    double di = it2->second - i1;
    while(it1++ != it2)
    {
      if( std::abs(di) < std::abs(it1->second - i1) )
        di = it1->second - i1;
    }
    return di;
  }
};



/**
 * dual exponential source, useful for cap charge/discharge
 */
class IEXP: public ISource
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
   * the initial current level
   */
  double I1;

  /**
   * the pulsed current level
   */
  double I2;

public:
  /**
   * constructor
   */
  IEXP(const std::string & s, double t1,double i1,double i2,double t2,double t3,double t4)
  :ISource(s),td(t1),trc(t2),tfd(t3),tfc(t4),I1(i1),I2(i2)
  {}


  /**
   * destructor have nothing to do
   */
  ~IEXP(){}


  /**
   * @return waveform of exponential current source
   */
  double iapp(double t)
  {
    if(t<=td)
      return I1;
    else if(t<=tfd)
      return I1+(I2-I1)*(1-exp(-(t-td)/trc));
    else
      return I1+(I2-I1)*(1-exp(-(t-td)/trc))+(I1-I2)*(1-exp(-(t-tfd)/tfc));
  }


  /**
   * @return the max value of iapp(t+delta_t)-iapp(t), delta_t in [0, dt]
   */
  double di_max(double t, double dt, double dt_min)
  {
    double t1 = (t-td);
    double t2 = (t+dt-td);
    double ta1 = tfd;

    std::map<double, double> i_wave;
    i_wave[t1 ] =  iapp(t1);
    i_wave[t2 ] =  iapp(t2);
    i_wave[ta1] =  iapp(ta1);

    std::map<double, double>::iterator it1 = i_wave.find(t1);
    std::map<double, double>::iterator it2 = i_wave.find(t2);

    double i1 = it1->second;
    double di = it2->second - i1;
    while(it1++ != it2)
    {
      if( std::abs(di) < std::abs(it1->second - i1) )
        di = it1->second - i1;
    }
    return di;
  }

};



/**
 * user defind source from analytic express
 */
class IUSER : public ISource
{
private:

  ExprEvalute expr_evaluator;

public:

 IUSER(const std::string &s, const std::string & expr) : ISource(s), expr_evaluator(expr)
 {}

 /**
  * @return waveform of exponential source
  */
  double iapp(double t)
  {
    return expr_evaluator(0,0,0,t);
  }

};


/**
 * A shell for user provide source from dll file
 */
class ISHELL: public ISource
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
  double (*Iapp_Shell)(double);

  /**
   * time scale parameter
   */
  double scale_t;

  /**
   * current scale parameter, scale to A
   */
  double scale_A;

public:

  /**
   * constructor
   */
#ifdef WINDOWS
  ISHELL(const std::string & s, HINSTANCE dp, void * fp, double s_t, double s_A);
#else
  ISHELL(const std::string & s, void * dp, void * fp, double s_t, double s_A);
#endif


  /**
   * destructor, free the pointer
   */
  ~ISHELL();


  /**
   * call Iapp_Shell to get the user provide current value
   */
  double iapp(double t)
  {
     return scale_A*Iapp_Shell(t/scale_t);
  }

};


#endif //#ifndef __isource_h__

