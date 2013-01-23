#ifndef __waveform_h_
#define __waveform_h_

#include <string>
#include <cmath>

#include "config.h"
#include "expr_evaluate.h"

#ifdef WINDOWS
class HINSTANCE__; // Forward or never
typedef HINSTANCE__* HINSTANCE;
#endif

class MonotCubicInterpolator;

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
class WaveformUniform: public Waveform
{
private:

  /**
   * time delay
   */
  double _td;

  /**
   * wave amplitude
   */
  double _amplitude;

public:

  /**
   * constructor
   */
  WaveformUniform(const std::string & s, double td,double amplitude)
  : Waveform(s), _td(td), _amplitude(amplitude)
  {}

  /**
   * destructor have nothing to do
   */
  ~WaveformUniform(){}

  /**
   * @return constant amplitude when t>td, else 0
   */
  double waveform(double t)
{ return t>=_td? _amplitude:0.0;}

};



/**
 * sine wave form
 */
class WaveformSin: public Waveform
{
private:
  /**
   * time delay
   */
  double _td;

  /**
   * bias amplitude
   */
  double _amplitude_bias;

  /**
   * sine wave amplitude
   */
  double _amplitude;

  /**
   * sine wave frequency
   */
  double _frequency;

  /**
   * attenuation parameter, should be >=0
   */
  double _alpha;

public:
  /**
   * constructor
   */
  WaveformSin(const std::string & s, double td,double a0,double a,double f,double alpha=0)
  : Waveform(s), _td(td),_amplitude_bias(a0),_amplitude(a),_frequency(f),_alpha(alpha)
  {}

  /**
     * destructor have nothing to do
   */
  ~WaveformSin(){}

  /**
   * @return sine wave when t>td, else return bias amplitude
   */
  double waveform(double t)
  {
    double amp = 0.0;
    if (t>=_td )
      amp = _amplitude_bias+_amplitude*exp(-_alpha*(t-_td))*sin(2*3.14159265358979323846*_frequency*(t-_td));
    else
      amp  = _amplitude_bias;
    return std::max(0.0, amp);
  }
};



/**
 * pulsed waveform
 */
class WaveformPulse: public Waveform
{
private:
  /**
   * time delay
   */
  double _td;

  /**
   * time for raising edge
   */
  double _tr;

  /**
   * time for falling edge
   */
  double _tf;

  /**
   * pulse width
   */
  double _pw;

  /**
   * pulse repeat period
   */
  double _pr;

  /**
   * the low amplitude of pulse
   */
  double _amplitude_low;

  /**
   * the high amplitude of pulse
   */
  double _amplitude_high;

public:

  /**
   * constructor
   */
  WaveformPulse(const std::string & s, double td,double a1,double a2,double tr,double tf,double pw, double pr)
  :Waveform(s),_td(td),_tr(tr),_tf(tf),_pw(pw),_pr(pr),_amplitude_low(a1),_amplitude_high(a2)
  {}

  /**
   * destructor have nothing to do
   */
  ~WaveformPulse(){}

  /**
   * @return waveform of pulsed source
   */
  double waveform(double t)
  {
    if(t<_td)
      return _amplitude_low;
    else
    {
      t-=_td;
      while(t>_pr) t-=_pr;
      if(t<_tr)
        return _amplitude_low+t*(_amplitude_high-_amplitude_low)/_tr;
      else if(t<_tr+_pw)
        return _amplitude_high;
      else if(t<_tr+_pw+_tf)
        return _amplitude_high-(t-_tr-_pw)*(_amplitude_high-_amplitude_low)/_tf;
      else    return _amplitude_low;
    }
  }

};



/**
 * double exponential waveform
 */
class WaveformExponential: public Waveform
{
private:
  /**
   * time delay
   */
  double _td;

  /**
   * time constant for raising edge
   */
  double _trc;

  /**
   * time delay for falling
   */
  double _tfd;

  /**
   * time constant for falling edge
   */
  double _tfc;

  /**
   * the low amplitude of pulse
   */
  double _amplitude_low;

  /**
   * the high amplitude of pulse
   */
  double _amplitude_high;

public:
  /**
   * constructor
   */
  WaveformExponential(const std::string & s, double td,double a1,double a2,double trc,double tfd,double tfc)
  :Waveform(s), _td(td),_trc(trc),_tfd(tfd),_tfc(tfc),_amplitude_low(a1),_amplitude_high(a2)
  {}

  /**
   * destructor have nothing to do
   */
  ~WaveformExponential(){}

  /**
     * @return waveform of exponential source
   */
  double waveform(double t)
  {
    if(t<=_td)
      return _amplitude_low;

    double exp1 = ( (t-_td)/_trc > 32) ? 0 : exp(-(t-_td)/_trc);
    double exp2 = ( (t-_tfd)/_tfc > 32) ? 0 : exp(-(t-_tfd)/_tfc);

    if(t<=_tfd)
      return _amplitude_low+(_amplitude_high-_amplitude_low)*(1-exp1);
    else
      return _amplitude_low+(_amplitude_high-_amplitude_low)*(1-exp1)+(_amplitude_low-_amplitude_high)*(1-exp2);

  }

};



/**
 * gauss waveform
 */
class WaveformGauss: public Waveform
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
  WaveformGauss(const std::string & s, double t0, double tao, double A)
      :Waveform(s), _t0(t0), _tao(tao), _amplitude(A)
  {}

  /**
   * destructor have nothing to do
   */
  ~WaveformGauss() {}

  /**
   * @return waveform of gauss source
   */
  double waveform(double t)
  {
    return _amplitude*exp(-(t-_t0)*(t-_t0)/(_tao*_tao));
  }

};


/**
 * user defind waveform from analytic express
 */

class WaveformExpression : public Waveform
{
private:

  ExprEvalute expr_evaluator;

public:

  WaveformExpression(const std::string & s, const std::string & expr)
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
class WaveformShell: public Waveform
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
  double (*Waveform_Shell)(double);

  /**
   * time scale parameter
   */
  double scale_t;

public:

  /**
   * constructor, hold the pointer to dll file and function
   */
  WaveformShell(const std::string & s, const std::string & dll_file, const std::string & function, double s_t);

  /**
   * destructor, free the pointer
   */
  ~WaveformShell();

  /**
   * call Waveform_Shell to get the user provide value
   */
  double waveform(double t)
  {
    return Waveform_Shell(t/scale_t);
  }

};


/**
 * wave form from file
 */
class WaveformFile: public Waveform
{
private:

  std::vector<double> _time;

  std::vector<double> _wave;

  MonotCubicInterpolator * _int;

public:

  /**
   * constructor, read the file
   */
  WaveformFile(const std::string & s, const std::string & file);

  /**
   * destructor
   */
  ~WaveformFile();

  /**
   * call Waveform_Shell to get the user provide value
   */
  double waveform(double t);
};


#endif

