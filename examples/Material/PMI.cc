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


#include <cmath>
#include <iomanip>
#include <fstream>

#include "point.h"
#include "fvm_node_data.h"

#include "PMI.h"

using namespace adtl;

/**
 * aux function return node coordinate.
 */
void PMI_Server::ReadCoordinate (PetscScalar& x, PetscScalar& y, PetscScalar& z) const
{
  if(pp_point)
  {
    x = (*pp_point)->x();
    y = (*pp_point)->y();
    z = (*pp_point)->z();
  }
  else
  {
    x = 0.0;
    y = 0.0;
    z = 0.0;
  }
}

/**
 * aux function return current time.
 */
PetscScalar PMI_Server::ReadTime () const
{
  if( p_clock )
    return *p_clock;
  return 0.0;
}


/**
 * check iff given variable eixst
 */
bool PMI_Server::HasVariable(const std::string &v, DataType t) const
{
  if ( !pp_variables || (*pp_variables)->find(v) == (*pp_variables)->end() ) return false;
  return ((*pp_variables)->find(v)->second.variable_data_type == t);
}


/**
 * @return index of given variable
 */
unsigned int PMI_Server::VariableIndex(const std::string &v) const
{
  if( !pp_variables || (*pp_variables)->find(v) == (*pp_variables)->end() ) return invalid_uint;
  return (*pp_variables)->find(v)->second.variable_index;
}


/**
 * aux function return scalar value of given variable.
 */
PetscScalar PMI_Server::ReadRealVariable (const unsigned int v) const
{
  if( pp_node_data )
    return (*pp_node_data)->data<Real>(v);
  return 0.0;
}


/**
 * aux function return scalar value of given variable.
 */
PetscScalar PMI_Server::ReadRealVariable (const std::string & v) const
{
  if( pp_node_data )
    return (*pp_node_data)->data<Real>(v);
  return 0.0;
}


#ifdef   __CALIBRATE__

/**
 * set numeric parameter value by its name.
 * @return 0 for success, 1 for variable no find.
 */
int PMI_Server::calibrate_real_parameter(const std::string & var_name, PetscScalar var_value)
{
  std::map<std::string, PARA>::iterator it = parameter_map.find(var_name);
  if( it != parameter_map.end() && it->second.type==PARA::Real)
  {
    *((PetscScalar*)it->second.value) = var_value*(it->second.unit_in_real);
    return 0;
  }

  _calibrate_error_info += "  unrecognized parameter " + var_name + "\n";

  return 1;
}

/**
 * set string parameter value by its name.
 * @return 0 for success, 1 for variable no find.
 */
int PMI_Server::calibrate_string_parameter(const std::string & var_name, const std::string &var_value)
{
  std::map<std::string, PARA>::iterator it = parameter_map.find(var_name);
  if( it != parameter_map.end() && it->second.type==PARA::String)
  {
    *((std::string*)it->second.value) = var_value;
    return 0;
  }

  _calibrate_error_info += "  unrecognized parameter " + var_name + "\n";

  return 1;
}

/**
 * set numeric and string parameters value by its name.
 */
int PMI_Server::calibrate(std::vector<Parser::Parameter> & pmi_parameters)
{
  int ierr = 0;
  std::vector<Parser::Parameter>::const_iterator it = pmi_parameters.begin();
  for(; it != pmi_parameters.end(); ++it)
  {
    switch (it->type())
    {
    case (Parser::REAL):
      ierr += this->calibrate_real_parameter(it->name(), it->get_real());
      break;
    case Parser::STRING :
      ierr += this->calibrate_string_parameter(it->name(), it->get_string());
      break;
    default:
      break;
    }
  }

  this->post_calibrate_process();

  return ierr;
}

#endif

/**
 * an interface for main code to access the parameter information in the material database
 */
std::map<std::string, PARA > & PMI_Server::get_parameter_info()
{
  return parameter_map;
}

/**
 * an interface for main code to access the library information in the material database
 */
const std::string & PMI_Server::get_PMI_info()
{
  return PMI_Info;
}

const std::string & PMI_Server::get_parameter_string(const int verbosity)
{
  std::stringstream output;

  // total screen width, and column width
  int wd = 120, wd_name=10, wd_unit=30, wd_val= 12, wd_sep=3;

  // desc text skip
  int tskip = wd_name + wd_unit + wd_val + wd_sep;
  // desc text width
  int twd = wd - tskip;

  output << std::setw(wd_name) << "Name"
         << std::setw(wd_val) << "Value"
         << std::setw(wd_unit) << "Unit" << "   "
         << std::setw(twd) << std::left << "Description" << std::right << std::endl;

  for( std::map<std::string, PARA>::const_iterator it = parameter_map.begin();
      it != parameter_map.end() ; ++it )
  {
    output << std::setw(wd_name) << it->second.name ;
    if ( it->second.type == PARA::String )
    {
      const std::string value = *((std::string*)it->second.value);
      output << std::setw(wd_val) << value;
    }
    else if ( it->second.type == PARA::Real )
    {
      const PetscScalar value =  *((PetscScalar*)it->second.value) / (it->second.unit_in_real);
      output << std::setw(wd_val) << value;
    }
    output << std::setw(wd_unit) << it->second.unit_in_string;
    output << std::setw(wd_sep) << "";

    int l = 0;
    for (std::string text = it->second.brief_intro; text.length()>0; text.erase(0,twd))
    {
      if (l>0)
        output << std::endl << std::setw(tskip) << "";

      output << text.substr(0,twd);
      l++;
    }
    output << std::right;
    output << std::endl;
  }
  _param_string = output.str();

  return _param_string;
}

/**
 * constructor, link PMI object to material class
 * also set the physical constants
 */
PMI_Server::PMI_Server(const PMI_Environment &env)
  : pp_variables(env.pp_variables), pp_point(env.pp_point), pp_node_data(env.pp_node_data), p_clock(env.p_clock)
{

  m  = env.m;
  s  = env.s;
  V  = env.V;
  C  = env.C;
  K  = env.K;

  cm = 1e-2*m;
  um = 1e-4*cm;
  nm = 1e-3*um;
  J  = C*V;
  W  = J/s;
  kg = J/(m*m)*s*s;
  g  = 1e-3*kg;
  eV = 1.602176462e-19*J;
  ps = 1e-12*s;
  A  = C/s;
  mA = 1e-3*A;
  Ohm= V/A;

  kb   = 1.3806503e-23*J/K;
  e    = 1.602176462e-19*C;
  me   = 9.10938188e-31*kg;
  eps0 = 8.854187818e-12*C/V/m;
  mu0  = 12.56637061e-7*std::pow(s,2)/C*V/m;
  h    = 6.62606876e-34*J*s;
  hbar = 1.054571596e-34*J*s;

  pi   = 3.14159265358979323846;
}


/*****************************************************************************
 *               Physical Model Interface for Semiconductor
 ****************************************************************************/

void PMIS_Server::SetFakeDopingEnvironment(PetscScalar Na, PetscScalar Nd)
{
  _Na = Na;
  _Nd = Nd;
}

void PMIS_Server::SetFakeMoleEnvironment(PetscScalar mole_x, PetscScalar mole_y)
{
  _mole_x = mole_x;
  _mole_y = mole_y;
}

void PMIS_Server::SetFakeDminEnvironment(PetscScalar dmin)
{
  _dmin = dmin;
}

void PMIS_Server::SetFakeStrainEnvironment(const TensorValue<PetscScalar> &T)
{
  _strain = T;
}


/**
 * aux function return first mole function of current node.
 */
PetscScalar PMIS_Server::ReadxMoleFraction () const
{
  if(pp_node_data) return (*pp_node_data)->mole_x();
  return _mole_x;
}

/**
 * aux function return first mole function with upper and lower bind of current node.
 */
PetscScalar PMIS_Server::ReadxMoleFraction (const PetscScalar mole_xmin, const PetscScalar mole_xmax) const
{
  if(pp_node_data)
  {
    PetscScalar mole_x=(*pp_node_data)->mole_x();
    if( mole_x < mole_xmin ) return mole_xmin;
    if( mole_x > mole_xmax ) return mole_xmax;
    return mole_x;
  }
  return _mole_x;
}

/**
 * aux function return second mole function of current node.
 */
PetscScalar PMIS_Server::ReadyMoleFraction () const
{
  if(pp_node_data) return (*pp_node_data)->mole_y();
  return _mole_y;
}

/**
 * aux function return second mole function with upper and lower bind of current node.
 */
PetscScalar PMIS_Server::ReadyMoleFraction (const PetscScalar mole_ymin, const PetscScalar mole_ymax) const
{
  if(pp_node_data)
  {
    PetscScalar mole_y=(*pp_node_data)->mole_y();
    if( mole_y < mole_ymin ) return mole_ymin;
    if( mole_y > mole_ymax ) return mole_ymax;
    return mole_y;
  }
  return _mole_y;
}

/**
 * aux function return total Acceptor concentration of current node
 */
PetscScalar PMIS_Server::ReadDopingNa () const
{
  if(pp_node_data)  return (*pp_node_data)->Total_Na();
  return _Na;
}

/**
 * aux function return total Donor concentration of current node
 */
PetscScalar PMIS_Server::ReadDopingNd () const
{
  if(pp_node_data) return (*pp_node_data)->Total_Nd();
  return _Nd;
}

/**
 * aux function return minimal distance to surface
 */
PetscScalar PMIS_Server::ReadDmin () const
{
  if(pp_node_data) return (*pp_node_data)->dmin();
  return _dmin;
}


/**
 * aux function return strain tensor
 */
TensorValue<PetscScalar> PMIS_Server::ReadStrain() const
{
  if(pp_node_data) return (*pp_node_data)->strain();
  return _strain;
}




/*****************************************************************************
 *               Physical Model Interface for Optical
 ****************************************************************************/


void RefractionSplineInterp::_build(const std::vector<PetscScalar> &t,
                                    const std::vector<PetscScalar> &y, std::vector<PetscScalar> &ypp)
{
  int size = t.size();

  std::vector<PetscScalar> a(3*size);
  std::vector<PetscScalar> b(size);

  //
  //  Set up the first equation. the second derivative at the left endpoint should be 0.0
  //
  {
    b[0] = 0.0; //second derivative
    a[1+0*3] = 1.0;
    a[0+1*3] = 0.0;
  }

  //
  //  Set up the intermediate equations.
  //
  for (int i = 1; i < size-1; i++ )
  {
    b[i] = ( y[i+1] - y[i] ) / ( t[i+1] - t[i] ) - ( y[i] - y[i-1] ) / ( t[i] - t[i-1] );
    a[2+(i-1)*3] = ( t[i] - t[i-1] ) / 6.0;
    a[1+ i   *3] = ( t[i+1] - t[i-1] ) / 3.0;
    a[0+(i+1)*3] = ( t[i+1] - t[i] ) / 6.0;
  }
  //
  //  Set up the last equation. the second derivative at the right endpoint should be 0.0
  //
  {
    b[size-1] = 0.0; //second derivative
    a[2+(size-2)*3] = 0.0;
    a[1+(size-1)*3] = 1.0;
  }

  //
  //  Solve the linear system.
  //
  if ( size == 2 )
  {
    ypp.resize(size);

    ypp[0] = 0.0;
    ypp[1] = 0.0;
  }
  else
  {
    ypp.resize(size);

    for (int i = 0; i < size; i++ )
    {
      ypp[i] = b[i];
    }

    for (int  i = 1; i < size; i++ )
    {
      double xmult = a[2+(i-1)*3] / a[1+(i-1)*3];
      a[1+i*3] = a[1+i*3] - xmult * a[0+i*3];
      ypp[i] = ypp[i] - xmult * ypp[i-1];
    }

    ypp[size-1] = ypp[size-1] / a[1+(size-1)*3];
    for (int  i = size-2; 0 <= i; i-- )
    {
      ypp[i] = ( ypp[i] - a[0+(i+1)*3] * ypp[i+1] ) / a[1+i*3];
    }
  }
}


PetscScalar RefractionSplineInterp::_eval(PetscScalar tval,
                                         const std::vector<PetscScalar> &t, const std::vector<PetscScalar> &y, const std::vector<PetscScalar> &ypp) const
{
  if(tval < t.front() ) return y.front();
  if(tval > t.back() )  return y.back();


  int size = t.size();
  double yval = 0.0;
  //
  //  Determine the interval [ double(I), double(I+1) ] that contains TVAL.
  //  Values below double[0] or above double[N-1] use extrapolation.
  //
  int begin=0, end=size-1;
  while(end-begin>2)
  {
    int mid = (begin+end)/2;
    if( tval < t[mid] )
      end = mid;
    else
      begin = mid;
  }

  int ival = size - 2;
  for (int i = begin; i < end; i++ )
  {
    if ( tval < t[i+1] )
    {
      ival = i;
      break;
    }
  }

  //
  //  In the interval I, the polynomial is in terms of a normalized
  //  coordinate between 0 and 1.
  //
  double dt = tval - t[ival];
  double h = t[ival+1] - t[ival];

  yval = y[ival]+ dt * ( ( y[ival+1] - y[ival] ) / h
      - ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * h
      + dt * ( 0.5 * ypp[ival] + dt * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0 * h ) ) ) );



  return  yval;
}



/**
 * when refraction_data_file is not empty, read from it
 */
void PMIS_Optical ::post_calibrate_process()
{
  if(_refraction_data_file.empty()) return;

  std::ifstream in(_refraction_data_file.c_str());
  if(! in.good()) return;

  _wave_table.clear();
  while(!in.eof())
  {
    std::string line;
    std::getline(in, line);

    if(in.fail()) break;

    std::stringstream ss(line);
    RefractionItem item;
    ss >> item.wavelength >> item.RefractionIndexRe >> item.RefractionIndexIm;
    _wave_table.push_back(item);
  }

  in.close();
}


const std::string& PMIS_Optical ::get_parameter_string(const int verbosity)
{
  std::stringstream output;

  // total screen width, and column width
  int wd = 120, wd_name=10, wd_unit=30, wd_val= 12, wd_sep=3;

  // desc text skip
  int tskip = wd_name + wd_unit + wd_val + wd_sep;
  // desc text width
  int twd = wd - tskip;

  output << std::setw(wd_name) << "lambda (um)"
      << std::setw(wd_val) << "  n  "
      << std::setw(wd_val) << "  k  " << "     "
      << std::setw(twd) << std::left << " Description " << std::right << std::endl;


  for( unsigned int n=0; n<_wave_table.size(); ++n)
  {
    output << std::setw(wd_name) << _wave_table[n].wavelength
        << std::setw(wd_val) << _wave_table[n].RefractionIndexRe
        << std::setw(wd_val) << _wave_table[n].RefractionIndexIm << std::endl;
  }

  _param_string = output.str();

  return _param_string;
}



/**
 * when refraction_data_file is not empty, read from it
 */
void PMII_Optical ::post_calibrate_process()
{
  if(_refraction_data_file.empty()) return;

  std::ifstream in(_refraction_data_file.c_str());
  if(! in.good()) return;

  _wave_table.clear();
  while(!in.eof())
  {
    std::string line;
    std::getline(in, line);

    if(in.fail()) break;

    std::stringstream ss(line);
    RefractionItem item;
    ss >> item.wavelength >> item.RefractionIndexRe >> item.RefractionIndexIm;
    _wave_table.push_back(item);
  }

  in.close();
}

/**
 * when refraction_data_file is not empty, read from it
 */
void PMIC_Optical ::post_calibrate_process()
{
  if(_refraction_data_file.empty()) return;

  std::ifstream in(_refraction_data_file.c_str());
  if(! in.good()) return;

  _wave_table.clear();
  while(!in.eof())
  {
    std::string line;
    std::getline(in, line);

    if(in.fail()) break;

    std::stringstream ss(line);
    RefractionItem item;
    ss >> item.wavelength >> item.RefractionIndexRe >> item.RefractionIndexIm;
    _wave_table.push_back(item);
  }

  in.close();
}






