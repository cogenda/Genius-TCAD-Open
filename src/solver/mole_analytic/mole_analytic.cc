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

#include "mesh_base.h"
#include "mole_analytic/mole_analytic.h"
#include "interpolation_2d_csa.h"
//#include "interpolation_3d_qshep.h"
#include "interpolation_3d_nbtet.h"
#include "parallel.h"

using PhysicalUnit::cm;
using PhysicalUnit::um;

/*------------------------------------------------------------------
 * we parse input deck for mole fraction here
 */
int MoleAnalytic::create_solver()
{
  // search again all the decks...
  // have any way with better efficency?
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    const Parser::Card c = _decks.get_current_card();
    if(c.key() == "MOLE" )
    {
      if(c.is_enum_value("type","file"))
        set_mole_function_file(c);
      else
        set_mole_function_linear(c);
    }
  }

  return 0;

}



/*------------------------------------------------------------------
 * assign mole fraction profile to semiconductor region
 */
int MoleAnalytic::solve()
{
  //search for all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    //we only process semiconductor region
    if( region->type() != SemiconductorRegion ) continue;

    SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
    SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = *node_it;
      FVM_NodeData * node_data = fvm_node->node_data();
      node_data->mole_x() = this->mole_x( region->name(), fvm_node->root_node() );
      node_data->mole_y() = this->mole_y( region->name(), fvm_node->root_node() );
    }
  }

  return 0;
}



//-----------------------------------------------------------------------
// private functions
void MoleAnalytic::set_mole_function_file(const Parser::Card & c)
{
  std::string fname = c.get_string("file", "");
  int skip_line      = c.get_int("skipline", 0);

  VectorValue<double> translate(c.get_real("translate.x", 0.0) * um,
                                  c.get_real("translate.y", 0.0) * um,
                                  c.get_real("translate.z", 0.0) * um);
  TensorValue<double> transform(c.get_real("transform.xx", 1.0), c.get_real("transform.xy", 0.0), c.get_real("transform.xz", 0.0),
                                c.get_real("transform.yx", 0.0), c.get_real("transform.yy", 1.0), c.get_real("transform.yz", 0.0),
                                c.get_real("transform.zx", 0.0), c.get_real("transform.zy", 0.0), c.get_real("transform.zz", 1.0));

  bool is_3D_mesh = _system.mesh().magic_num() < 2008 ? false : true;
  PetscScalar LUnit;
  if (c.is_enum_value("lunit", "m"))
    LUnit = 1.0;
  else if (c.is_enum_value("lunit", "cm"))
    LUnit = 1e-2;
  else if (c.is_enum_value("lunit", "um"))
    LUnit = 1e-6;
  else if (c.is_enum_value("lunit", "nm"))
    LUnit = 1e-9;
  else
    LUnit = 1e-6;

  InterpolationBase *interpolatorX, *interpolatorY;
  if(is_3D_mesh)
  {
    interpolatorX = new Interpolation3D_nbtet; // 3D mesh
    interpolatorY = new Interpolation3D_nbtet; // 3D mesh
  }
  else
  {
    interpolatorX = new Interpolation2D_CSA; // 2D mesh
    interpolatorY = new Interpolation2D_CSA; // 2D mesh
  }

  interpolatorX->set_interpolation_type(0, InterpolationBase::Linear);
  interpolatorY->set_interpolation_type(0, InterpolationBase::Linear);

  Point p;
  double moleX, moleY;

  int hasY = 0;
  if(Genius::processor_id()==0)
  {
    std::ifstream in(fname.c_str());
    if(!in.good())
    {
      MESSAGE<<"ERROR at " << c.get_fileline() <<": can not open mole profile file " << fname << std::endl; RECORD();
      genius_error();
    }

    MESSAGE<<"Loading mole profile from " << fname << ".\n"; RECORD();

    int i;
    for (i=0; i<skip_line && !in.eof(); i++)
    {
      std::string line;
      std::getline(in, line);
    }

    while(!in.eof() && in.good())
    {
      std::string buf;
      std::getline(in, buf);
      std::istringstream line(buf);

      line >> p[0];// read x location
      line >> p[1];// read y location
      if(is_3D_mesh)
        line >> p[2];// read z location

      p *= PhysicalUnit::m*LUnit;     // scale to length unit
      p = transform*p + translate; //do transform & translate

      line >> moleX;

      if(!line.fail())
      {
        interpolatorX->add_scatter_data(p, 0, moleX);

        line >> moleY;
        if(!line.fail())
        {
          interpolatorY->add_scatter_data(p, 0, moleY);
          hasY = 1;
        }
      }
      i++;
    }
    if(!in.eof())
    {
      MESSAGE<<"ERROR at " << c.get_fileline() <<": error reading mole profile file " << fname;
      MESSAGE<<" at line " << i << "." << std::endl; RECORD();
      genius_error();
    }
    in.close();
  }
  Parallel::broadcast(hasY, 0);

  if (hasY)
  {
    interpolatorX->broadcast(0);
    interpolatorX->setup(0);
    interpolatorY->broadcast(0);
    interpolatorY->setup(0);
    _mole_data.push_back(std::pair<InterpolationBase*, InterpolationBase*>(interpolatorX, interpolatorY));
  }
  else
  {
    interpolatorX->broadcast(0);
    interpolatorX->setup(0);
    _mole_data.push_back(std::pair<InterpolationBase*, InterpolationBase*>(interpolatorX, (InterpolationBase*)NULL));
    delete interpolatorY;
  }
}

void MoleAnalytic::set_mole_function_linear(const Parser::Card & c)
{
  std::string region = c.get_string("region", "");

  std::pair<Point, Point> bbox = _system.region(region)->region_bound_box();

  // get the bound box
  double xmin = c.get_real("x.min", (bbox.first) (0)/um, "x.left")*um;
  double xmax = c.get_real("x.max", (bbox.second)(0)/um, "x.right")*um;
  double ymin = c.get_real("y.min", (bbox.first) (1)/um, "y.top")*um;
  double ymax = c.get_real("y.max", (bbox.second)(1)/um, "y.bottom")*um;
  double zmin = c.get_real("z.min", (bbox.first) (2)/um, "z.front")*um;
  double zmax = c.get_real("z.max", (bbox.second)(2)/um, "z.back")*um;

  genius_assert(xmin<=xmax);
  genius_assert(ymin<=ymax);
  genius_assert(zmin<=zmax);

  // the begin and end location of mole fraction
  double x_base, x_end;

  // read the mole fraction and its slope from user's input
  double x_mole       =  c.get_real("x.mole", 0.0);
  double x_mole_slope =  c.get_real("x.mole.slope", 0.0)/um;

  // read the slope direction
  Direction x_dir     =  Y_Direction;
  if( c.is_enum_value("x.mole.grad", "x.linear") )
    x_dir = X_Direction;
  if( c.is_enum_value("x.mole.grad", "z.linear") )
    x_dir = Z_Direction;

  // determine the begin and end location
  switch(x_dir)
  {
  case X_Direction: x_base = xmin; x_end = xmax; break;
  case Y_Direction: x_base = ymin; x_end = ymax; break;
  case Z_Direction: x_base = zmin; x_end = zmax; break;
  default : genius_error();
  }

  // if user input x.mole.end instead of x.mole.slope, convert it
  if(c.is_parameter_exist("x.mole.end") && !c.is_parameter_exist("x.mole.slope"))
  {
    double x_mole_end   =  c.get_real("x.mole.end", 0.0);
    x_mole_slope = (x_mole_end - x_mole)/(x_end - x_base);
  }

  //build linear mole function for first mole fraction
  MoleFunction * mf = new LinearMoleFunction(region, x_base, x_dir, x_mole, x_mole_slope, true);

  //push it into _mole_funs vector
  _mole_funs.push_back(mf);


  // for complex compound material, process second mole fraction
  if(c.is_parameter_exist("y.mole"))
  {
    // the begin and end location of mole fraction
    double y_base, y_end;

    // read the mole fraction and its slope from user's input
    double y_mole       =  c.get_real("y.mole", 0.0);
    double y_mole_slope =  c.get_real("y.mole.slope", 0.0)/um;

    // read the slope direction
    Direction y_dir     =  Y_Direction;
    if( c.is_enum_value("y.mole.grad", "x.linear") )
      y_dir = X_Direction;
    if( c.is_enum_value("y.mole.grad", "z.linear") )
      y_dir = Z_Direction;

    // determine the begin and end location
    switch(x_dir)
    {
    case X_Direction: y_base = xmin; y_end = xmax; break;
    case Y_Direction: y_base = ymin; y_end = ymax; break;
    case Z_Direction: y_base = zmin; y_end = zmax; break;
    default : genius_error();
    }

    // if user input x.mole.end instead of x.mole.slope, convert it
    if(c.is_parameter_exist("y.mole.end") && !c.is_parameter_exist("y.mole.slope"))
    {
      double y_mole_end   =  c.get_real("y.mole.end", 0.0);
      y_mole_slope = (y_mole_end - y_mole)/(y_end - y_base);
    }

    //build linear mole function for first mole fraction
    MoleFunction * mf = new LinearMoleFunction(region, y_base, y_dir, y_mole, y_mole_slope, false);

    //push it into _mole_funs vector
    _mole_funs.push_back(mf);

  }
}



double MoleAnalytic::mole_x(const std::string & region, const Node * node)
{
  double mole = 0.0;

  // use the first Mole profile that specified an X composition
  for(size_t i=0; i<_mole_funs.size(); i++)
  {
    if(_mole_funs[i]->region() != region) continue;

    if(_mole_funs[i]->mole_x())
      return  _mole_funs[i]->profile((*node)(0),(*node)(1),(*node)(2));
  }

  // then check mole profile file
  for(size_t i=0; i<_mole_data.size(); i++)
  {
    InterpolationBase * interp = _mole_data[i].first;
    if (!interp)
      continue;
    double mx = interp->get_interpolated_value(*node, 0);
    mx = mx<0.0 ? 0.0 : mx;
    mx = mx>1.0 ? 1.0 : mx;
#if defined(HAVE_FENV_H) && defined(DEBUG)
    if(fetestexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW))
    {
      MESSAGE<< "Warning: problem in interpolating x mole fraction at ";
      MESSAGE<< (*node)(0)/um << "\t";
      MESSAGE<< (*node)(1)/um << "\t";
      MESSAGE<< (*node)(2)/um << " " << mx << " , ignored.\n";
      RECORD();
      feclearexcept(FE_ALL_EXCEPT);
      return mole;
    }
#endif
    return mx;
  }

  return mole;
}



double MoleAnalytic::mole_y(const std::string & region, const Node * node)
{
  double mole = 0.0;

  for(size_t i=0; i<_mole_funs.size(); i++)
  {
    if(_mole_funs[i]->region() != region) continue;

    if(_mole_funs[i]->mole_y())
      return  _mole_funs[i]->profile((*node)(0),(*node)(1),(*node)(2));
  }

  // then check mole profile file
  for(size_t i=0; i<_mole_data.size(); i++)
  {
    InterpolationBase * interp = _mole_data[i].second;
    if (!interp)
      continue;
    double mx = interp->get_interpolated_value(*node, 0);
    mx = mx<0.0 ? 0.0 : mx;
    mx = mx>1.0 ? 1.0 : mx;
#if defined(HAVE_FENV_H) && defined(DEBUG)
    if(fetestexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW))
    {
      MESSAGE<< "Warning: problem in interpolating x mole fraction at ";
      MESSAGE<< (*node)(0)/um << "\t";
      MESSAGE<< (*node)(1)/um << "\t";
      MESSAGE<< (*node)(2)/um << " " << mx << " , ignored.\n";
      RECORD();
      feclearexcept(FE_ALL_EXCEPT);
      return mole;
    }
#endif
    return mx;
  }

  return mole;
}

