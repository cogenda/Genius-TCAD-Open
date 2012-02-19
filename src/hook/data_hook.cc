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
#include <time.h>
#include <string>
#include <cstdlib>
#include <iomanip>

#include "mesh_base.h"
#include "solver_base.h"
#include "data_hook.h"
#include "mxml.h"
#include "MXMLUtil.h"

using PhysicalUnit::um;

/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
DataHook::DataHook ( SolverBase & solver, const std::string & name, void * param )
  : Hook ( solver, name ), _count(0), _output_prefix ( SolverSpecify::out_prefix )
{

  const SimulationSystem &system = get_solver().get_system();
  _is_2d = system.mesh().mesh_dimension() == 2;


  const std::vector<Parser::Parameter> & parm_list = * ( ( std::vector<Parser::Parameter> * ) param );
  for ( std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
        parm_it != parm_list.end(); parm_it++ )
  {
    if ( parm_it->name() == "variable" )
    {
      _variable_name.push_back ( parm_it->get_string() );
    }

    if ( parm_it->name() == "region" )
    {
      std::string region = parm_it->get_string();
      if(system.region ( region ))
        _region.push_back( region );
      else
        if( Genius::is_first_processor() )
          std::cerr<<"DataHook: Invalid given region "<< region  <<  ", ignored." << std::endl;
    }

    if ( parm_it->name() == "regions" )
    {
      std::vector<std::string> regions = parm_it->get_array<std::string>();
      for(unsigned int n=0; n<regions.size(); ++n)
      {
        std::string region = regions[n];
        if(system.region ( region ))
          _region.push_back( region );
        else
          if( Genius::is_first_processor() )
            std::cerr<<"DataHook: Invalid given region "<< region  <<  ", ignored." << std::endl;
      }
    }

    if ( parm_it->name() == "output" )
    {
      _output_prefix = parm_it->get_string();
    }
  }

  // when no region specified by user, record all
  if( _region.empty() )
  {
    for ( unsigned int r=0; r<system.n_regions(); ++r )
    {
      const SimulationRegion * region = system.region ( r );
      _region.push_back( region->name() );
    }
  }
}


/*----------------------------------------------------------------------
 * destructor, close file
 */
DataHook::~DataHook()
{

}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void DataHook::on_init()
{

}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void DataHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void DataHook::post_solve()
{
  std::ostringstream data_filename;
  data_filename << _output_prefix << ( this->_count++ ) << ".fd";

  if ( !Genius::processor_id() )
    _out.open ( data_filename.str().c_str() );

  // prepare the file head
  // only root processor do this
  if ( !Genius::processor_id() )
  {
    // get simulation time
    time_t simulation_time;
    time ( &simulation_time );

    // write file head
    _out << "# Title: Field Data File Created by Genius TCAD Simulation" << std::endl;
    _out << "# Date: " << ctime ( &simulation_time ) << std::endl;

    // write variables
    _out << "# Variables: " << std::endl;
    _out << '#' <<'\t' << "x" << " [um]"<< std::endl;
    _out << '#' <<'\t' << "y" << " [um]"<< std::endl;
    _out << '#' <<'\t' << "z" << " [um]"<< std::endl;
    if(_is_2d)
      _out << '#' <<'\t' << "cell_area" << " [um^2]"<< std::endl;
    else
      _out << '#' <<'\t' << "cell_volumn" << " [um^3]"<< std::endl;
    for ( unsigned int n=0; n<_variable_name.size(); ++n )
    {
      SolutionVariable variable = solution_string_to_enum ( FormatVariableString(_variable_name[n]) );
      _out << '#' <<'\t' << _variable_name[n] << " [" << variable_unit_string ( variable ) << "]" <<std::endl;
    }

    _out<<std::endl;

    // set the float number precision
    _out.precision ( 6 );

    // set output width and format
    _out<< std::scientific << std::right;

  }

  const SimulationSystem &system = get_solver().get_system();
  const MeshBase & mesh = system.mesh();

  // allgather node location, parallel code
  std::vector<Real> pts;
  mesh.pack_nodes(pts);

  for ( unsigned int r=0; r<_region.size(); ++r )
  {
    const SimulationRegion * region = system.region ( _region[r] );
    genius_assert(region);
    if ( !Genius::processor_id() )
      _out << "# Region: " << region->name() << std::endl;

    // data buffer
    std::vector< std::vector<Real> > values;

    // sync field data between all the processor, must executed in parallel!
    for ( unsigned int n=0; n<_variable_name.size(); ++n )
    {
      SimulationVariable v;
      bool find = region->get_variable(FormatVariableString(_variable_name[n]), POINT_CENTER, v)  ;
      if(!find) continue;

      switch ( v.variable_data_type )
      {
        case  SCALAR :
        {
          std::vector<Real> value;
          region->get_variable_data<Real>( FormatVariableString(_variable_name[n]), POINT_CENTER, value );
          values.push_back ( value );
          break;
        }
        default: break;
      }
    }

    // allgather region node, parallel code
    std::vector<unsigned int> nodes;
    region->region_node(nodes);


    // save data
    if ( !Genius::processor_id() )
    {
      for (unsigned int n=0 ; n<nodes.size(); ++n )
      {
        unsigned int node_id = nodes[n];
        _out << std::setw ( 15 ) << pts[3*node_id+0]/um;
        _out << std::setw ( 15 ) << pts[3*node_id+1]/um;
        _out << std::setw ( 15 ) << pts[3*node_id+2]/um;
        for (unsigned int v=0; v<values.size(); ++v)
        {
          const std::vector<Real> & value = values[v];
          _out << std::setw ( 15 ) << value[n];
        }
        _out << std::endl;
      }
    }


  }



  // close file
  if ( !Genius::processor_id() )
    _out.close();
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void DataHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void DataHook::on_close()
{}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new DataHook ( solver, name, fun_data );
  }

}

#endif

