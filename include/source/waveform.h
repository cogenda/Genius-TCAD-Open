#ifndef __waveform_h_
#define __waveform_h_

#include <string>
#include <cmath>

#include "config.h"
#include "expr_evaluate.h"

#ifdef CYGWIN
  #include <Windows.h>
  #undef max
  #undef min
#else
  #include <dlfcn.h>
#endif


/**
 * basic class of Waveform
 */
class Waveform
{
public:
   Waveform(const std::string & s): _label(s) {}

   ~Waveform() {};

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


  virtual double waveform(double )=0;

private:

  /**
   * the unique label for a Waveform
   */
  std::string _label;

};




/**
 * Uniform Waveform
 */
class UniformWaveform: public Waveform
{
private:

  /**
   * time delay
   */
  double td;

  /**
   * wave amplitude
   */
  double amplitude;

public:

  /**
   * constructor
   */
  UniformWaveform(const std::string & s, double t1,double v1)
  : Waveform(s), td(t1), amplitude(v1)
  {}

  /**
   * destructor have nothing to do
   */
  ~UniformWaveform(){}

  /**
   * @return constant amplitude when t>td, else 0
   */
  double waveform(double t)
  { return t>=td? amplitude:0.0;}

};




/**
 * pulsed waveform
 */
class PulseWaveform: public Waveform
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
   * the low amplitude of pulse
   */
  double amplitude_low;

  /**
   * the high amplitude of pulse
   */
  double amplitude_high;

public:

  /**
   * constructor
   */
  PulseWaveform(const std::string & s, double t1,double v1,double v2,double t2,double t3,double t4, double t5)
  :Waveform(s),td(t1),tr(t2),tf(t3),pw(t4),pr(t5),amplitude_low(v1),amplitude_high(v2)
  {}

  /**
   * destructor have nothing to do
   */
  ~PulseWaveform(){}

  /**
   * @return waveform of pulsed source
   */
  double waveform(double t)
  {
    if(t<td)
      return amplitude_low;
    else
    {
      t-=td;
      while(t>pr) t-=pr;
      if(t<tr)
        return amplitude_low+t*(amplitude_high-amplitude_low)/tr;
      else if(t<tr+pw)
        return amplitude_high;
      else if(t<tr+pw+tf)
        return amplitude_high-(t-tr-pw)*(amplitude_high-amplitude_low)/tf;
      else    return amplitude_low;
    }
  }

};



/**
 * gauss waveform
 */
class GaussWaveform: public Waveform
{
private:
  /**
   * pulse center time
   */
  double _t0;

  /**
   * gauss constant
   */
  double _tao;

  /**
   * the amplitude of gauss
   */
  double _amplitude;


public:
  /**
   * constructor
   */
  GaussWaveform(const std::string & s, double t0, double tao, double A)
  :Waveform(s), _t0(t0), _tao(tao), _amplitude(A)
  {}

  /**
   * destructor have nothing to do
   */
  ~GaussWaveform() {}

  /**
   * @return waveform of gauss source
   */
  double waveform(double t)
  {
    return _amplitude*exp(-(t-_t0)*(t-_t0)/(2.0*_tao*_tao));
  }

};


/**
 * user defind waveform from analytic express
 */

class ExprWaveform : public Waveform
{
private:

  ExprEvalute expr_evaluator;

public:

 ExprWaveform(const std::string & s, const std::string & expr)
 :Waveform(s),expr_evaluator(expr)
 {}

 /**
  * @return user defined waveform
  */
  double waveform(double t)
  {
    return expr_evaluator(0,0,0,t);
  }

};



/**
 * A shell for user provide waveform from dll file
 */
class ShellWaveform: public Waveform
{
private:
  /**
   * the pointer to dll file
   */
#ifdef CYGWIN
   HINSTANCE                  dll;
#else
   void                      *dll;
#endif

  /**
   * the pointer to function in the dll file
   */
  double (*Waveform_Shell)(double);

  /**
   * time scale parameter
   */
  double scale_t;

public:

  /**
   * constructor, hold the pointer to dll file and function
   */
#ifdef CYGWIN
  ShellWaveform(const std::string & s, HINSTANCE dp, void * fp, double s_t):Waveform(s)
#else
  ShellWaveform(const std::string & s, void * dp, void * fp, double s_t):Waveform(s)
#endif
  {
     dll = dp;
     Waveform_Shell = (double (*)(double)) fp;
     assert(Waveform_Shell);
     scale_t = s_t;
  }

  /**
   * destructor, free the pointer
   */
  ~ShellWaveform()
  {
     //delete Vapp_Shell;
#ifdef CYGWIN
    FreeLibrary(dll);
#else
    if ( dll )
      dlclose( dll );
#endif
  }

  /**
   * call Waveform_Shell to get the user provide value
   */
  double waveform(double t)
  {
     return Waveform_Shell(t/scale_t);
  }

};


#endif

