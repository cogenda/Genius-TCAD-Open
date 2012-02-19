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

//  $Id: doping_analytic.cc,v 1.10 2008/07/09 05:58:16 gdiso Exp $


#include "polygon.h"
#include "mesh_base.h"
#include "doping_analytic/doping_fun.h"
#include "doping_analytic/doping_analytic.h"
#include "semiconductor_region.h"
#include "interpolation_1d_linear.h"
#include "interpolation_1d_spline.h"
#include "interpolation_2d_csa.h"
#include "interpolation_2d_nn.h"
//#include "interpolation_3d_qshep.h"
#include "interpolation_3d_nbtet.h"

using PhysicalUnit::cm;
using PhysicalUnit::um;


DopingAnalytic::DopingAnalytic(SimulationSystem & system, Parser::InputParser & decks)
  : PointSolver(system), _decks(decks)
{system.record_active_solver(this->solver_type());}


DopingAnalytic::~DopingAnalytic()
{
  for(size_t n=0; n<_doping_funs.size(); n++)
    delete _doping_funs[n];
  _doping_funs.clear();
  for(size_t n=0; n<_doping_data.size(); n++)
    delete _doping_data[n].second;
  _doping_data.clear();
  for(std::map<std::string,DopingFunction *>::iterator it=_custom_profile_funs.begin();
      it!=_custom_profile_funs.end(); it++)
    delete it->second;
  _custom_profile_funs.clear();
}


/*------------------------------------------------------------------
 * we parse input deck for doping profile here
 */
int DopingAnalytic::create_solver()
{
  // search again all the decks...
  // have any way with better efficency?
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    const Parser::Card c = _decks.get_current_card();
    if(c.key() == "PROFILE.DOPING" || c.key() == "PROFILE")
    {
      // PROFILE card should have "type" parameter
      if( !c.is_parameter_exist("type") )
      {
        MESSAGE<<"ERROR at " << c.get_fileline() <<" PROFILE: Should have a 'type' parameter."<<std::endl; RECORD();
        genius_error();
      }

      // set doping function here.
      if(c.is_enum_value("type","uniform"))
        set_doping_function_uniform(c);

      if(c.is_enum_value("type","analytic"))
        set_doping_function_analytic(c);

      if(c.is_enum_value("type","analytic2"))
      {
        if(c.is_parameter_exist("mask.polygon"))
          set_doping_function_analytic2_poly_mask(c);
        else
          set_doping_function_analytic2_rec_mask(c);
      }
      if(c.is_enum_value("type","file"))
        set_doping_function_file(c);
    }
  }

  return 0;

}

int DopingAnalytic::destroy_solver()
{
  return 0;
}

int DopingAnalytic::pre_solve_process(bool)
{
  return 0;
}

int DopingAnalytic::post_solve_process()
{
  return 0;
}

/*------------------------------------------------------------------
 * assign doping profile to semiconductor region
 */
int DopingAnalytic::solve()
{
  //search for all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    //we only process semiconductor region
    if( region->type() != SemiconductorRegion ) continue;
    SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<SemiconductorSimulationRegion *>(region);
    // prepare region custom defined variable
    std::map<std::string, std::pair<unsigned int, int> > ion_map;
    for (std::map<std::string,DopingFunction *>::iterator it = _custom_profile_funs.begin();
         it!=_custom_profile_funs.end(); it++)
    {
      const std::string & name = it->first;
      unsigned int ion_index = region->add_variable(SimulationVariable(name, SCALAR, POINT_CENTER, "cm^-3", invalid_uint, true, true));
      int ion_type = semiconductor_region->material()->band->IonType(name);
      ion_map.insert(std::make_pair(name, std::make_pair(ion_index, ion_type)));
    }

    SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
    SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = *node_it;
      FVM_NodeData * node_data = fvm_node->node_data();
      genius_assert(node_data!=NULL);

      const Node * node = fvm_node->root_node();
      node_data->Na() = doping_Na( node );
      node_data->Nd() = doping_Nd( node );
      // fill custom defined variable
      for (std::map<std::string,DopingFunction *>::iterator it = _custom_profile_funs.begin();
           it!=_custom_profile_funs.end(); it++)
      {
        const std::string & name = it->first;
        double d = it->second->profile((*node)(0),(*node)(1),(*node)(2));
        node_data->data<Real>(ion_map[name].first) = d;
        if(ion_map[name].second < 0 ) node_data->Na() += d;
        if(ion_map[name].second > 0 ) node_data->Nd() += d;
      }
    }

  }

  return 0;
}


//-----------------------------------------------------------------------
// private functions

void DopingAnalytic::set_doping_function_file(const Parser::Card & c)
{
  std::string fname = c.get_string("file", "");
  int skip_line     = c.get_int("skipline", 0);

  VectorValue<double> translate;
  TensorValue<double> transform;
  if( c.is_parameter_exist("translate") )
  {
    std::vector<double> dummy = c.get_array<double>("translate");
    translate = VectorValue<double>(&dummy[0])*um;
  }
  else
  {
    translate = VectorValue<double>(c.get_real("translate.x", 0.0)*um,
                                    c.get_real("translate.y", 0.0)*um,
                                    c.get_real("translate.z", 0.0)*um);
  }

  if( c.is_parameter_exist("transform") )
  {
    std::vector<double> dummy = c.get_array<double>("transform");
    transform = TensorValue<double>(&dummy[0]);
  }
  else
  {
    transform = TensorValue<double>(c.get_real("transform.xx", 1.0), c.get_real("transform.xy", 0.0), c.get_real("transform.xz", 0.0),
                                    c.get_real("transform.yx", 0.0), c.get_real("transform.yy", 1.0), c.get_real("transform.yz", 0.0),
                                    c.get_real("transform.zx", 0.0), c.get_real("transform.zy", 0.0), c.get_real("transform.zz", 1.0));
  }

  double ion  = 0;
  if(c.is_enum_value("ion","donor"))
    ion = 1.0;
  else if(c.is_enum_value("ion","acceptor"))
    ion = -1.0;

  bool is_3D_mesh = (_system.mesh().mesh_dimension() == 3);
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

  int axes=AXES_XY;  // default axes for 2D mesh
  if(is_3D_mesh)
    axes=AXES_XYZ; // default axes for 3D mesh

  if (c.is_enum_value("axes", "x"))
    axes=AXES_X;
  else if (c.is_enum_value("axes", "y"))
    axes=AXES_Y;
  else if (c.is_enum_value("axes", "z"))
    axes=AXES_Z;
  else if (c.is_enum_value("axes", "xy"))
    axes=AXES_XY;
  else if (c.is_enum_value("axes", "xz"))
    axes=AXES_XZ;
  else if (c.is_enum_value("axes", "yz"))
    axes=AXES_YZ;
  else if (c.is_enum_value("axes", "xyz"))
    axes=AXES_XYZ;

  InterpolationBase * interpolator;
  switch(axes)
  {
  case AXES_X:
  case AXES_Y:
  case AXES_Z:
    interpolator = new Interpolation1D_Linear; // 1D profile
    break;
  case AXES_XY:
  case AXES_XZ:
  case AXES_YZ:
    interpolator = new Interpolation2D_NN; // 2D profile
    //interpolator = new Interpolation2D_CSA; // 2D profile
    break;
  case AXES_XYZ:
    interpolator = new Interpolation3D_nbtet; // 3D profile
  }

  interpolator->set_interpolation_type(0, InterpolationBase::Linear);

  Point p;
  double doping;

  if(Genius::processor_id()==0)
  {
    std::ifstream in(fname.c_str());
    if(!in.good())
    {
      MESSAGE<<"ERROR at " << c.get_fileline() <<": can not open doping file " << fname << std::endl; RECORD();
      genius_error();
    }

    MESSAGE<<"Loading doping profile from " << fname << ".\n"; RECORD();

    int i;
    for (i=0; i<skip_line && !in.eof(); i++)
    {
      std::string line;
      std::getline(in, line);
    }

    while(!in.eof() && in.good())
    {
      switch(axes)
      {
      case AXES_X:
        in >> p[0]; break; // read the only coordinate
      case AXES_Y:
        in >> p[1]; break;
      case AXES_Z:
        in >> p[2]; break;
      case AXES_XY:
        in >> p[0] >> p[1]; break; // read two coordinates
      case AXES_XZ:
        in >> p[0] >> p[2]; break; // read two coordinates
        break;
      case AXES_YZ:
        in >> p[1] >> p[2]; break; // read two coordinates
      case AXES_XYZ:
        in >> p[0] >> p[1] >> p[2]; break; // read three coordinates
      }

      p *= PhysicalUnit::m*LUnit;     // scale to length unit
      p = transform*p + translate; //do transform & translate

      Point p1;
      switch(axes)
      {
      case AXES_X:
        p1[0] = p[0]; break;
      case AXES_Y:
        p1[0] = p[1]; break;
      case AXES_Z:
        p1[0] = p[2]; break;
      case AXES_XY:
        p1[0] = p[0]; p1[1] = p[1]; break;
      case AXES_XZ:
        p1[0] = p[0]; p1[1] = p[2]; break;
      case AXES_YZ:
        p1[0] = p[1]; p1[1] = p[2]; break;
      case AXES_XYZ:
        p1[0] = p[0]; p1[1] = p[1]; p1[2] = p[2]; break;
      }

      in >> doping;
      doping *= ion;

      if(!in.fail())
      {
        interpolator->add_scatter_data(p1, 0, doping);
      }
      i++;
    }
    if(!in.eof())
    {
      MESSAGE<<"ERROR at " << c.get_fileline() <<": error reading doping file " << fname;
      MESSAGE<<" at line " << i << "." << std::endl; RECORD();
      genius_error();
    }
    in.close();
  }

  interpolator->broadcast(0);
  interpolator->setup(0);
  _doping_data.push_back(std::pair<int, InterpolationBase * >(axes, interpolator));
}

void DopingAnalytic::set_doping_function_uniform(const Parser::Card & c)
{

  // get the doping profile bond box
  double xmin = c.get_real("x.min", 0.0, "x.left")*um;
  double xmax = c.get_real("x.max", 0.0, "x.right")*um;
  double ymin = c.get_real("y.min", 0.0, "y.top")*um;
  double ymax = c.get_real("y.max", 0.0, "y.bottom")*um;
  double zmin = c.get_real("z.min", 0.0, "z.front")*um;
  double zmax = c.get_real("z.max", 0.0, "z.back")*um;

  if(xmin>xmax || ymin>ymax || zmin>zmax)
  {
    MESSAGE<<"ERROR at " << c.get_fileline() <<" PROFILE: Profile has incorrect XYZ bound."<<std::endl; RECORD();
    genius_error();
  }

  // the peak value of doping profile
  double peak = c.get_real("n.peak",0.0)/pow(cm,3);
  genius_assert(peak>=0.0);

  // the ion type
  genius_assert(c.is_parameter_exist("ion"));

  // explicit specified ion
  if(c.is_enum_value("ion","donor") || c.is_enum_value("ion","acceptor"))
  {
    double ion  = 0;
    if(c.is_enum_value("ion","donor"))
      ion = 1.0;
    else if(c.is_enum_value("ion","acceptor"))
      ion = -1.0;

    //build uniform doping function
    DopingFunction * df = new UniformDopingFunction(ion,xmin,xmax,ymin,ymax,zmin,zmax,peak);

    //push it into _doping_funs vector
    _doping_funs.push_back(df);
  }

  // custom specified, the N or P ion will be determined by region material vai PMI
  if(c.is_enum_value("ion","custom"))
  {
    genius_assert(c.is_parameter_exist("id"));
    std::string name = c.get_string("id", "custom_profile");
    //build uniform doping function
    DopingFunction * df = new UniformDopingFunction(1.0,xmin,xmax,ymin,ymax,zmin,zmax,peak);
    _custom_profile_funs.insert(std::pair<std::string,DopingFunction*>(name,df));
  }
}



void DopingAnalytic::set_doping_function_analytic(const Parser::Card & c)
{

  // get the doping profile mask box
  double xmin = c.get_real("x.min", 0.0, "x.left")*um;
  double xmax = c.get_real("x.max", 0.0, "x.right")*um;
  double ymin = c.get_real("y.min", 0.0, "y.top")*um;
  double ymax = c.get_real("y.max", 0.0, "y.bottom")*um;
  double zmin = c.get_real("z.min", 0.0, "z.front")*um;
  double zmax = c.get_real("z.max", 0.0, "z.back")*um;

  if(xmin>xmax || ymin>ymax || zmin>zmax)
  {
    MESSAGE<<"ERROR at " << c.get_fileline() <<" PROFILE: Porfile has incorrect XYZ bound."<<std::endl; RECORD();
    genius_error();
  }

  // the peak value of doping profile
  double peak = c.get_real("n.peak",0.0)/pow(cm,3);
  genius_assert(peak>=0.0);

  // the ion type
  genius_assert(c.is_parameter_exist("ion"));
  double ion  = 0;
  if(c.is_enum_value("ion","donor"))
    ion = 1.0;
  else if(c.is_enum_value("ion","acceptor"))
    ion = -1.0;

  // we should get characteristic length in every direction
  double YCHAR, XCHAR, ZCHAR;
  // the depth of junction in y direction (norm to the silicon)
  double YJUNC;
  // aux variable, the background doping concentration
  double dop=0;
  // the centre of doping profile
  double slice_x = 0.5*(xmin+xmax);
  double slice_z = 0.5*(zmin+zmax);

  //characteristic length in y direction (norm to the silicon surface)
  // however, there are several ways...
  if(c.is_parameter_exist("y.char"))          // directly
    YCHAR = c.get_real("y.char",0.25)*um;
  else if(c.is_parameter_exist("y.junction")) // by "y.junction"
  {
    //Junction is an absolute location.
    YJUNC = c.get_real("y.junction",0.0)*um;

    //get concentration of background doping profile
    for(size_t i=0;i<_doping_funs.size();i++)
      dop += _doping_funs[i]->profile(slice_x,YJUNC,slice_z);

    //Can we even find a junction?
    genius_assert(dop*ion<0.0);
    genius_assert(peak>0.0);

    //Now convert junction depth into char. length
    YCHAR = (YJUNC-ymin)/sqrt(log(fabs(peak/dop)));
  }

  //characteristic length in x direction (parallel to the silicon surface)
  if(c.is_parameter_exist("x.char"))
    XCHAR = c.get_real("x.char",YCHAR/um)*um;
  else if(c.is_parameter_exist("xy.ratio"))
    XCHAR = YCHAR*c.get_real("xy.ratio",1.0);
  else
    XCHAR = YCHAR;

  //characteristic length in z direction (also parallel to the silicon surface)
  if(c.is_parameter_exist("z.char"))
    ZCHAR = c.get_real("z.char",YCHAR/um)*um;
  else if(c.is_parameter_exist("zy.ratio"))
    ZCHAR = YCHAR*c.get_real("zy.ratio",1.0);
  else
    ZCHAR = YCHAR;

  genius_assert(XCHAR>=0 && YCHAR>=0 && ZCHAR>=0);

  // we can get peak value by dose, when YCHAR is known
  if(c.is_parameter_exist("dose"))
  {
    double dose = c.get_real("dose",0.0)/pow(cm,2);
    genius_assert(dose>0.0);
    peak=dose/(YCHAR*sqrt(M_PI));
  }

  // check if we should use erfc in x or z direction
  bool erfcx = c.get_bool("x.erfc",false);
  bool erfcy = c.get_bool("y.erfc",false);
  bool erfcz = c.get_bool("z.erfc",false);

  // explicit specified ion
  if(c.is_enum_value("ion","donor") || c.is_enum_value("ion","acceptor"))
  {
    //build analytic doping function
    DopingFunction * df = new AnalyticDopingFunction(ion,xmin,xmax,ymin,ymax,zmin,zmax,peak,XCHAR,YCHAR,ZCHAR,erfcx,erfcy,erfcz);

    //push it into _doping_funs vector
    _doping_funs.push_back(df);
  }

  // custom specified, the N or P ion will be determined by region material vai PMI
  if(c.is_enum_value("ion","custom"))
  {
    genius_assert(c.is_parameter_exist("id"));
    std::string name = c.get_string("id", "custom_profile");
    //build analytic doping function
    DopingFunction * df = new AnalyticDopingFunction(1.0,xmin,xmax,ymin,ymax,zmin,zmax,peak,XCHAR,YCHAR,ZCHAR,erfcx,erfcy,erfcz);
    _custom_profile_funs.insert(std::pair<std::string,DopingFunction*>(name,df));
  }

}


void DopingAnalytic::set_doping_function_analytic2_poly_mask(const Parser::Card & c)
{

  // get the doping profile mask polygon
  std::vector<Point> poly;

  std::vector<double> coods = c.get_array<double>("mask.polygon");
  if( coods.size() % 3 != 0 )
  {
    MESSAGE<<"ERROR at " << c.get_fileline() <<" PROFILE: Mask polygon has incorrect coordinates."<<std::endl; RECORD();
    genius_error();
  }
  if( coods.size() < 9 )
  {
    MESSAGE<<"ERROR at " << c.get_fileline() <<" PROFILE: Mask polygon should have at least 3 vertex."<<std::endl; RECORD();
    genius_error();
  }

  for(unsigned int i=0; i<coods.size();)
  {
    poly.push_back(Point(coods[i++]*um, coods[i++]*um, coods[i++]*um));
  }

  Polygon polygon(poly);
  if(!polygon.valid())
  {
    MESSAGE<<"ERROR at " << c.get_fileline() <<" PROFILE: Mask polygon should coplaner."<<std::endl; RECORD();
    genius_error();
  }


  // the peak value of doping profile
  double peak = c.get_real("n.peak",0.0)/pow(cm,3);
  genius_assert(peak>=0.0);

  // the ion type
  genius_assert(c.is_parameter_exist("ion"));
  double ion  = 0;
  if(c.is_enum_value("ion","donor"))
    ion = 1.0;
  else if(c.is_enum_value("ion","acceptor"))
    ion = -1.0;

  double rmin = c.get_real("implant.rmin", 0.0)*um;
  double rmax = c.get_real("implant.rmax", 0.0)*um;

  double theta = c.get_real("implant.theta", 0.0);
  double phi = c.get_real("implant.phi", 0.0);

  double char_depth = c.get_real("depth.char", 0.1)*um;
  double char_lateral = c.get_real("lateral.char", 0.1)*um;

  double resolution_factor = c.get_real("resolution", 4.0);

  // explicit specified ion
  if(c.is_enum_value("ion","donor") || c.is_enum_value("ion","acceptor"))
  {
    //build analytic doping function
    DopingFunction * df = new PolyMaskDopingFunction(ion,poly,theta,phi,rmin,rmax,peak,char_depth,char_lateral,resolution_factor);

    //push it into _doping_funs vector
    _doping_funs.push_back(df);
  }

  // custom specified, the N or P ion will be determined by region material vai PMI
  if(c.is_enum_value("ion","custom"))
  {
    genius_assert(c.is_parameter_exist("id"));
    std::string name = c.get_string("id", "custom_profile");
    //build analytic doping function
    DopingFunction * df = new PolyMaskDopingFunction(1.0,poly,theta,phi,rmin,rmax,peak,char_depth,char_lateral,resolution_factor);
    _custom_profile_funs.insert(std::pair<std::string,DopingFunction*>(name,df));
  }

}


void DopingAnalytic::set_doping_function_analytic2_rec_mask(const Parser::Card & c)
{

  // get the doping profile mask box
  double xmin = c.get_real("mask.xmin", 0.0, "x.min")*um;
  double xmax = c.get_real("mask.xmax", 0.0, "x.max")*um;
  double ymin = c.get_real("mask.ymin", 0.0, "y.min")*um;
  double ymax = c.get_real("mask.ymax", 0.0, "y.max")*um;
  double zmin = c.get_real("mask.zmin", 0.0, "z.min")*um;
  double zmax = c.get_real("mask.zmax", 0.0, "z.max")*um;

  if(xmin>xmax || ymin>ymax || zmin>zmax)
  {
    MESSAGE<<"ERROR at " << c.get_fileline() <<" PROFILE: Mask has incorrect XYZ bound."<<std::endl; RECORD();
    genius_error();
  }

  if(c.is_parameter_exist("mask.x"))
    xmin = xmax = c.get_real("mask.x", 0.0)*um;
  if(c.is_parameter_exist("mask.y"))
    ymin = ymax = c.get_real("mask.y", 0.0)*um;
  if(c.is_parameter_exist("mask.z"))
    zmin = zmax = c.get_real("mask.z", 0.0)*um;

  // the peak value of doping profile
  double peak = c.get_real("n.peak",0.0)/pow(cm,3);
  genius_assert(peak>=0.0);

  // the ion type
  genius_assert(c.is_parameter_exist("ion"));
  double ion  = 0;
  if(c.is_enum_value("ion","donor"))
    ion = 1.0;
  else if(c.is_enum_value("ion","acceptor"))
    ion = -1.0;

  double rmin = c.get_real("implant.rmin", 0.0)*um;
  double rmax = c.get_real("implant.rmax", 0.0)*um;

  double theta = c.get_real("implant.theta", 0.0);
  double phi = c.get_real("implant.phi", 0.0);

  double char_depth = c.get_real("depth.char", 0.1)*um;
  double char_lateral = c.get_real("lateral.char", 0.1)*um;

  double resolution_factor = c.get_real("resolution", 4.0);

  // explicit specified ion
  if(c.is_enum_value("ion","donor") || c.is_enum_value("ion","acceptor"))
  {
    //build analytic doping function
    DopingFunction * df = new RecMaskDopingFunction(ion,xmin,xmax,ymin,ymax,zmin,zmax,theta,phi,rmin,rmax,peak,char_depth,char_lateral,resolution_factor);

    //push it into _doping_funs vector
    _doping_funs.push_back(df);
  }

  // custom specified, the N or P ion will be determined by region material vai PMI
  if(c.is_enum_value("ion","custom"))
  {
    genius_assert(c.is_parameter_exist("id"));
    std::string name = c.get_string("id", "custom_profile");
    //build analytic doping function
    DopingFunction * df = new RecMaskDopingFunction(1.0,xmin,xmax,ymin,ymax,zmin,zmax,theta,phi,rmin,rmax,peak,char_depth,char_lateral,resolution_factor);
    _custom_profile_funs.insert(std::pair<std::string,DopingFunction*>(name,df));
  }

}


double DopingAnalytic::_do_doping_interp(int i, const Node *node, const std::string &msg)
{
  int axes = _doping_data[i].first;
  InterpolationBase * interp = _doping_data[i].second;

  Point p;
  switch(axes)
  {
    case AXES_X:
      p[0]=(*node)(0);
      break;
    case AXES_Y:
      p[0]=(*node)(1);
      break;
    case AXES_Z:
      p[0]=(*node)(2);
      break;
    case AXES_XY:
      p[0]=(*node)(0);
      p[1]=(*node)(1);
      break;
    case AXES_XZ:
      p[0]=(*node)(0);
      p[1]=(*node)(2);
      break;
    case AXES_YZ:
      p[0]=(*node)(1);
      p[1]=(*node)(2);
      break;
    case AXES_XYZ:
      p[0]=(*node)(0);
      p[1]=(*node)(1);
      p[2]=(*node)(2);
      break;
  }
  double d = interp->get_interpolated_value(p, 0);
#if defined(HAVE_FENV_H) && defined(DEBUG)
  if(fetestexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW))
  {
    MESSAGE<< "Warning: problem in interpolating " << msg << " at ";
    MESSAGE<< (*node)(0)/um << "\t";
    MESSAGE<< (*node)(1)/um << "\t";
    MESSAGE<< (*node)(2)/um << " " << d << " , ignored.\n";
    RECORD();
    feclearexcept(FE_ALL_EXCEPT);
  }
#endif
  return d;
}

double DopingAnalytic::doping_Na(const Node * node)
{
  double dop = 0.0;

  //only add negative value
  for(size_t i=0; i<_doping_funs.size(); i++)
  {
    double d = _doping_funs[i]->profile((*node)(0),(*node)(1),(*node)(2));
    dop += d < 0.0 ? d: 0.0;
  }
  double unit = 1.0/std::pow(PhysicalUnit::cm,3.0);
  for(size_t i=0; i<_doping_data.size(); i++)
  {
    double d = unit * _do_doping_interp(i, node, "Na");
    dop += d < 0.0 ? d: 0.0;
  }

  return std::abs(dop);
}



double DopingAnalytic::doping_Nd(const Node * node)
{
  double dop = 0.0;

  //only add positive value
  for(size_t i=0; i<_doping_funs.size(); i++)
  {
    double d = _doping_funs[i]->profile((*node)(0),(*node)(1),(*node)(2));
    dop += d > 0.0 ? d: 0.0;
  }
  double unit = 1.0/std::pow(PhysicalUnit::cm,3.0);
  for(size_t i=0; i<_doping_data.size(); i++)
  {
    double d = unit * _do_doping_interp(i, node, "Nd");
    dop += d > 0.0 ? d: 0.0;
  }

  return std::abs(dop);
}
