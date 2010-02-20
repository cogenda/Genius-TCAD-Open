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



// system include
#include <string>
#include <cstdlib>


// genius include
#include "solver_base.h"
#include "solver_specify.h"
#include "visit_hook.h"
#include "parallel.h"

/******************************************************************************
 *  VisitHook static member
 *****************************************************************************/

VisitState   VisitHook::visit_state = UNKNOWN;
SolverBase*  VisitHook::_p_solver   = 0;
std::vector<std::pair<std::string, std::string> >  VisitHook::_variables;
std::vector< std::vector<double> > VisitHook::_values;
unsigned int VisitHook::_n_values   = 0;
std::map<unsigned int, int> VisitHook::global_id_to_index;



/******************************************************************************
 *  Visit parallel aux functions
 *****************************************************************************/

static int visit_broadcast_int_callback(int *value, int sender)
{return MPI_Bcast(value, 1, MPI_INT, sender, PETSC_COMM_WORLD);}

static int visit_broadcast_string_callback(char *str, int len, int sender)
{return MPI_Bcast(str, len, MPI_CHAR, sender, PETSC_COMM_WORLD);}

#define VISIT_COMMAND_PROCESS 0
#define VISIT_COMMAND_SUCCESS 1
#define VISIT_COMMAND_FAILURE 2

/* Helper function for ProcessVisItCommand */
static void BroadcastSlaveCommand(int *command)
{
  MPI_Bcast(command, 1, MPI_INT, 0, PETSC_COMM_WORLD);
}

/* Callback involved in command communication. */
static void SlaveProcessCallback()
{
  int command = VISIT_COMMAND_PROCESS;
  BroadcastSlaveCommand(&command);
}

/* aux function for run VisItProcessEngineCommand() on all processor. */
static int ProcessVisItCommand(void)
{
  int command;
  if (Genius::processor_id()==0)
  {
    int success = VisItProcessEngineCommand();

    if (success)
    {
      command = VISIT_COMMAND_SUCCESS;
      BroadcastSlaveCommand(&command);
      return 1;
    }
    else
    {
      command = VISIT_COMMAND_FAILURE;
      BroadcastSlaveCommand(&command);
      return 0;
    }
  }
  else
  {
    /* Note: only through the SlaveProcessCallback callback
     * above can the rank 0 process send a VISIT_COMMAND_PROCESS
     * instruction to the non-rank 0 processes. */
    while (1)
    {
      BroadcastSlaveCommand(&command);
      switch (command)
      {
      case VISIT_COMMAND_PROCESS:
        VisItProcessEngineCommand();
        break;
      case VISIT_COMMAND_SUCCESS:
        return 1;
      case VISIT_COMMAND_FAILURE:
        return 0;
      }
    }
  }
}



/******************************************************************************
 *  VisitHook member function
 *****************************************************************************/



/*----------------------------------------------------------------------
 * constructor, connect to visit
 */
VisitHook::VisitHook(SolverBase & _solver, const std::string & name, void * dummy)
    : Hook(_solver, name)
{
  _p_solver = &_solver;

  if(SolverSpecify::Type != SolverSpecify::TRANSIENT && SolverSpecify::Type != SolverSpecify::DCSWEEP ) return;

  connect_to_visit();
}



/*----------------------------------------------------------------------
 * destructor, close the connection to visit
 */
VisitHook::~VisitHook()
{}



/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void VisitHook::on_init()
{
  if(SolverSpecify::Type != SolverSpecify::TRANSIENT && SolverSpecify::Type != SolverSpecify::DCSWEEP ) return;

  // record electrode IV (and time) information
  {
    // if transient simulation, we need to record time
    if ( SolverSpecify::Type == SolverSpecify::TRANSIENT )
      _variables.push_back( std::pair<std::string, std::string>("time", "time") );

    // record electrode IV information
    const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
    // search for all the bcs
    for(unsigned int n=0; n<bcs->n_bcs(); n++)
    {
      const BoundaryCondition * bc = bcs->get_bc(n);
      // skip bc which is not electrode
      if( !bc->is_electrode() ) continue;
      //
      _variables.push_back( std::pair<std::string, std::string>(bc->label() + "_Vapp", "voltage") );
      _variables.push_back( std::pair<std::string, std::string>(bc->label() + "_potential", "voltage")  );
      _variables.push_back( std::pair<std::string, std::string>(bc->label() + "_current", "current") );
    }

    _values.resize( _variables.size() );
  }
}




/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void VisitHook::pre_solve()
{}




/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void VisitHook::post_solve()
{
  if(SolverSpecify::Type != SolverSpecify::TRANSIENT && SolverSpecify::Type != SolverSpecify::DCSWEEP ) return;

  // record electrode IV (and time) information
  {
    // variable counter
    unsigned int i=0;

    // if transient simulation, we need to record time
    if (SolverSpecify::Type == SolverSpecify::TRANSIENT)
      _values[i++].push_back( SolverSpecify::clock/PhysicalUnit::s );

    // search for all the bc
    const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
    for(unsigned int n=0; n<bcs->n_bcs(); n++)
    {
      const BoundaryCondition * bc = bcs->get_bc(n);
      // skip bc which is not electrode
      if( !bc->is_electrode() ) continue;

      _values[i++].push_back( bc->ext_circuit()->Vapp()/PhysicalUnit::V );
      _values[i++].push_back( bc->ext_circuit()->potential()/PhysicalUnit::V );
      _values[i++].push_back( bc->ext_circuit()->current()/PhysicalUnit::A );
    }

    if(i)  _n_values++;
  }


  // recive information from visit and answer it!
  comm_to_visit();

}




/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void VisitHook::post_iteration()
{}




/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void VisitHook::on_close()
{

  if(visit_state  != CONNECTED) return;
  std::cout << "Please exit visit to continue..." << std::endl;
  // wait for visit until it exist
  disconnect_to_visit();
}



int VisitHook::visit_cell_type(const Elem *elem)
{
  int celltype;

  switch(elem->type())
  {
  case EDGE2:
  case EDGE2_FVM:
    celltype = VISIT_CELL_BEAM;
    break;
  case TRI3:
  case TRI3_FVM:
    celltype = VISIT_CELL_TRI;
    break;// 3
  case QUAD4:
  case QUAD4_FVM:
    celltype = VISIT_CELL_QUAD;
    break;// 5
  case TET4:
  case TET4_FVM:
    celltype = VISIT_CELL_TET;
    break;// 8
  case HEX8:
  case HEX8_FVM:
    celltype = VISIT_CELL_HEX;
    break;// 10
  case PRISM6:
  case PRISM6_FVM:
    celltype = VISIT_CELL_WEDGE;
    break;// 13
  case PYRAMID5:
  case PYRAMID5_FVM:
    celltype = VISIT_CELL_PYR;
    break;// 16
  default:
    {
      std::cerr<<"element type "<<elem->type()<<" not implemented"<<std::endl;
      genius_error();
    }
  }

  return celltype;
}



/******************************************************************************
 * ACK routine
 *****************************************************************************/


void VisitHook::connect_to_visit()
{
  /* this should be run in parallel */
  parallel_only();

  /* Initialize environment variables. */
  VisItSetupEnvironment();

  /* Install callback functions for global communication. */
  VisItSetBroadcastIntFunction(visit_broadcast_int_callback);
  VisItSetBroadcastStringFunction(visit_broadcast_string_callback);

  /* Tell libsim whether the simulation is parallel. */
  VisItSetParallel(Genius::n_processors() > 1);
  VisItSetParallelRank(Genius::processor_id());

  /* Write out .sim file that VisIt uses to connect. Only do it on processor 0 */
  if(Genius::processor_id()==0)
  {
    VisItInitializeSocketAndDumpSimFile("genius",
                                        "Genius real time display",
                                        "",
                                        NULL,
                                        NULL,
                                        "genius.sim1");


    // we fork a child process to run visit
    pid_t pid = fork();

    // child process
    if (pid == (pid_t) 0)
    {
      // get visit executable file
      char *s = getenv("VISIT");
      if (s)
      { 
        std::string visit_dir(s);
        std::string visit_exec = visit_dir + "/../../bin/" + "visit";

        // exec visit and try to link with genius
        int err = execl(visit_exec.c_str(), "visit", "-o", "genius.sim1", 0);
      }
      exit(0);
    }
    else if (pid < (pid_t) 0)
    {
      // The fork failed.
      genius_error();
    }
  }

  // wait until connected to visit
  do
  {
    /* Get input from VisIt */
    int visit_input;
    if(Genius::processor_id()==0)
      visit_input = VisItDetectInput(1, -1);

    MPI_Bcast(&visit_input, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    switch (visit_input)
    {
    case 0 :/* There was no input from VisIt. Continue waiting . */
      {
        break;
      }

    case 1 :/* VisIt is trying to connect to sim. */
      {
        if(VisItAttemptToCompleteConnection())
        {
          visit_state  = CONNECTED;
	  std::cout << "VisIt connected! \n"<<std::endl;
          VisItSetSlaveProcessCallback(SlaveProcessCallback);
          //VisItSetCommandCallback(ControlCommandCallback);
        }
        else
        {
          std::cout << "VisIt NOT connected! \n" << VisItGetLastError() << std::endl;
          visit_state  = DISCONNECTED;
        }

        break;
      }

    case 2 :/* VisIt wants to tell the engine something. */
      {
        if(!ProcessVisItCommand())
        {
          /* Disconnect on an error or closed connection. */
          VisItDisconnect();
          visit_state  = DISCONNECTED;
          remove("genius.sim1");
        }
        break;
      }
    default: break;
    }
  }
  while( visit_state == UNKNOWN );

}



void VisitHook::comm_to_visit()
{

  if(visit_state  != CONNECTED) return;

  do
  {

    /* Get input from VisIt */
    int visit_input;
    if(Genius::processor_id()==0)
      visit_input = VisItDetectInput(0, -1);
    MPI_Bcast(&visit_input, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    switch (visit_input)
    {
    case 0 :/* There was no input from VisIt. Send simulation data to visit */
      {

        /**
         * we should use this function to tell visit that simulation updated
         */
        VisItTimeStepChanged();
        VisItUpdatePlots();
        return;
      }

    case 2 :/* VisIt wants to tell the engine something. */
      {
        if(!ProcessVisItCommand())
        {
          /* Disconnect on an error or closed connection. */
          VisItDisconnect();
          visit_state  = DISCONNECTED;
          if(Genius::processor_id()==0)
            remove("genius.sim1");
          return;
        }
        break;
      }
    default: break;
    }
  }
  while(1);
}





void VisitHook::disconnect_to_visit()
{

  do
  {
    /* Get input from VisIt */
    int visit_input;
    if(Genius::processor_id()==0)
      visit_input = VisItDetectInput(0, -1);
    MPI_Bcast(&visit_input, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    switch (visit_input)
    {
    case 0 :/* There was no input from VisIt. */
      {
        break;
      }

    case 2 :/* VisIt wants to tell the engine something. */
      {
        if(!ProcessVisItCommand())
        {
          /* Disconnect on an error or closed connection. */
          VisItDisconnect();
          visit_state  = UNKNOWN;
          if(Genius::processor_id()==0)
            remove("genius.sim1");
          return;
        }
        break;
      }
    default: break;
    }
  }
  while(1);
}




/******************************************************************************
 * Data access routine
 *****************************************************************************/


/**
 * Data access function for metadata, only processor 0 process it!
 */
VisIt_SimulationMetaData * VisitHook::get_meta_data(void)
{

  /* Create a metadata object with no variables. */
  size_t sz = sizeof(VisIt_SimulationMetaData);
  VisIt_SimulationMetaData *md = (VisIt_SimulationMetaData *)malloc(sz);
  memset(md, 0, sz);

  /* Set the simulation state. */
  md->currentMode =  VISIT_SIMMODE_RUNNING;

  /* transient simulation */
  if(SolverSpecify::TimeDependent == true)
  {
    md->currentCycle = SolverSpecify::T_Cycles;
    md->currentTime = SolverSpecify::clock/PhysicalUnit::s;
  }
  else /* dc simulation */
  {
    md->currentCycle = SolverSpecify::DC_Cycles;
    md->currentTime = 0;
  }


  /* Allocate enough room for mesh in the metadata. */
  md->numMeshes = 1;
  sz = sizeof(VisIt_MeshMetaData) * md->numMeshes;
  md->meshes = (VisIt_MeshMetaData *)malloc(sz);
  memset(md->meshes, 0, sz);

  /* Set mesh properties.*/
  md->meshes[0].name = strdup("Genius Mesh3D");
  md->meshes[0].meshType = VISIT_MESHTYPE_UNSTRUCTURED;
  md->meshes[0].topologicalDimension = 3;
  md->meshes[0].spatialDimension = 3;
  md->meshes[0].numBlocks = Genius::n_processors();
  md->meshes[0].blockTitle = strdup("Domains");
  md->meshes[0].blockPieceName = strdup("domain");
  md->meshes[0].numGroups = 0;
  md->meshes[0].units = strdup("um");
  md->meshes[0].xLabel = strdup("Width");
  md->meshes[0].yLabel = strdup("Height");
  md->meshes[0].zLabel = strdup("Depth");


  /* Add some scalar variables. */
  md->numScalars = 3;
  sz = sizeof(VisIt_ScalarMetaData) * md->numScalars;
  md->scalars = (VisIt_ScalarMetaData *)malloc(sz);
  memset(md->scalars, 0, sz);

  /* Add potential nodal variable on mesh3d. */
  md->scalars[0].name = strdup("potential");
  md->scalars[0].meshName = strdup("Genius Mesh3D");
  md->scalars[0].centering = VISIT_VARCENTERING_NODE;

  /* Add electron density nodal variable on mesh3d. */
  md->scalars[1].name = strdup("electron_density");
  md->scalars[1].meshName = strdup("Genius Mesh3D");
  md->scalars[1].centering = VISIT_VARCENTERING_NODE;

  /* Add hole density nodal variable on mesh3d. */
  md->scalars[2].name = strdup("hole_density");
  md->scalars[2].meshName = strdup("Genius Mesh3D");
  md->scalars[2].centering = VISIT_VARCENTERING_NODE;

  /* Add curve variable. */

  if ( SolverSpecify::Type == SolverSpecify::TRANSIENT )
    md->numCurves = _variables.size()-1;
  else
    md->numCurves = _variables.size()/3*2;

  sz = sizeof(VisIt_CurveMetaData) * md->numCurves;
  md->curves = (VisIt_CurveMetaData *)malloc(sz);
  memset(md->curves, 0, sz);

  unsigned int i=0;

  // search for all the bc
  const BoundaryConditionCollector * bcs = _p_solver->get_system().get_bcs();
  for(unsigned int n=0; n<bcs->n_bcs(); n++)
  {
    const BoundaryCondition * bc = bcs->get_bc(n);
    // skip bc which is not electrode
    if( !bc->is_electrode() ) continue;

    if (SolverSpecify::Type == SolverSpecify::TRANSIENT)
    {
      // plot Vapp vs Time
      std::string name = "Transient " + bc->label() + " Vapp vs Time";
      md->curves[3*i].name   = strdup(name.c_str());
      md->curves[3*i].xUnits = strdup("second");
      md->curves[3*i].xLabel = strdup("time");
      md->curves[3*i].yUnits = strdup("V");
      md->curves[3*i].yLabel = strdup("voltage");

      // plot Potential vs Time
      name = "Transient " + bc->label() + " Potential vs Time";
      md->curves[3*i+1].name   = strdup(name.c_str());
      md->curves[3*i+1].xUnits = strdup("second");
      md->curves[3*i+1].xLabel = strdup("time");
      md->curves[3*i+1].yUnits = strdup("V");
      md->curves[3*i+1].yLabel = strdup("potential");

      // plot Current vs Time
      name = "Transient " + bc->label() + " Current vs Time";
      md->curves[3*i+2].name   = strdup(name.c_str());
      md->curves[3*i+2].xUnits = strdup("second");
      md->curves[3*i+2].xLabel = strdup("time");
      md->curves[3*i+2].yUnits = strdup("A");
      md->curves[3*i+2].yLabel = strdup("current");
    }

    if (SolverSpecify::Type == SolverSpecify::DCSWEEP)
    {
      // plot  Current vs Vapp
      std::string name = "DC " + bc->label() + " Current vs Vapp";
      md->curves[2*i].name   = strdup(name.c_str());
      md->curves[2*i].xUnits = strdup("V");
      md->curves[2*i].xLabel = strdup("voltage");
      md->curves[2*i].yUnits = strdup("A");
      md->curves[2*i].yLabel = strdup("current");

      // plot Current vs Potential
      name = "DC " + bc->label() + " Current vs Potential";
      md->curves[2*i+1].name   = strdup(name.c_str());
      md->curves[2*i+1].xUnits = strdup("V");
      md->curves[2*i+1].xLabel = strdup("potential");
      md->curves[2*i+1].yUnits = strdup("A");
      md->curves[2*i+1].yLabel = strdup("current");
    }

    i++;
  }


  return md;
}



VisIt_CurveData * VisitHook::get_electrode_iv(const std::string & visit_cv_name)
{

  size_t sz = sizeof(VisIt_CurveData);
  VisIt_CurveData *curve = (VisIt_CurveData*)malloc(sz);
  memset(curve, 0, sz);

  float *x = NULL;
  float *y = NULL;
  x = (float*)malloc(_n_values * sizeof(float));
  y = (float*)malloc(_n_values * sizeof(float));

  // search for all the bc
  unsigned int value_index=0;
  const BoundaryConditionCollector * bcs = _p_solver->get_system().get_bcs();
  for(unsigned int n=0; n<bcs->n_bcs(); n++)
  {
    const BoundaryCondition * bc = bcs->get_bc(n);
    // skip bc which is not electrode
    if( !bc->is_electrode() ) continue;

    if (SolverSpecify::Type == SolverSpecify::TRANSIENT)
    {
      // plot Vapp vs Time
      std::string name = "Transient " + bc->label() + " Vapp vs Time";
      if( name == visit_cv_name )
      {
        for(unsigned int i = 0; i < _n_values; ++i)
        {
          x[i] = _values[0][i];
          y[i] = _values[3*value_index+0][i];
        }
      }

      // plot Potential vs Time
      name = "Transient " + bc->label() + " Potential vs Time";
      if( name == visit_cv_name )
      {
        for(unsigned int i = 0; i < _n_values; ++i)
        {
          x[i] = _values[0][i];
          y[i] = _values[3*value_index+1][i];
        }
      }

      // plot Current vs Time
      name = "Transient " + bc->label() + " Current vs Time";
      if( name == visit_cv_name )
      {
        for(unsigned int i = 0; i < _n_values; ++i)
        {
          x[i] = _values[0][i];
          y[i] = _values[3*value_index+2][i];
        }
      }

    }


    if (SolverSpecify::Type == SolverSpecify::DCSWEEP)
    {
      // plot  Current vs Vapp
      std::string name = "DC " + bc->label() + " Current vs Vapp";
      if( name == visit_cv_name )
      {
        for(unsigned int i = 0; i < _n_values; ++i)
        {
          x[i] = _values[3*value_index+0][i];
          y[i] = _values[3*value_index+2][i];
        }
      }

      // plot  Current vs Potential
      name = "DC " + bc->label() + " Current vs Potential";
      if( name == visit_cv_name )
      {
        for(unsigned int i = 0; i < _n_values; ++i)
        {
          x[i] = _values[3*value_index+1][i];
          y[i] = _values[3*value_index+2][i];
        }
      }
    }

    value_index++;
  }

  // Give the arrays to VisIt. VisIt will free them since VISIT_OWNER_VISIT flag is set.
  curve->len = _n_values;
  curve->x = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_VISIT, x);
  curve->y = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_VISIT, y);

  return curve;
}



VisIt_MeshData *  VisitHook::get_mesh(int domain,  const std::string &name)
{

  VisIt_MeshData *mesh = NULL;
  size_t sz = sizeof(VisIt_MeshData);

  genius_assert( domain < Genius::n_processors());
  genius_assert( name == "Genius Mesh3D" );

  const MeshBase & _mesh = _p_solver->get_system().mesh();

  static std::vector<float> x;
  static std::vector<float> y;
  static std::vector<float> z;
  static std::vector<int>   connectivity;

  /* Allocate VisIt_MeshData. */
  mesh = (VisIt_MeshData *)malloc(sz);
  memset(mesh, 0, sz);

  /* Make VisIt_MeshData contain a VisIt_UnstructuredMesh. */
  sz = sizeof(VisIt_UnstructuredMesh);
  mesh->umesh = (VisIt_UnstructuredMesh *)malloc(sz);memset(mesh->umesh, 0, sz);

  /* Tell VisIt which mesh object to use. */
  mesh->meshType = VISIT_MESHTYPE_UNSTRUCTURED;
  /* Set the mesh¡¯s number of dimensions. */
  mesh->umesh->ndims = 3;

  /* Set the number of nodes and zones in the mesh domain. */
  mesh->umesh->nnodes = _mesh.n_nodes();
  mesh->umesh->nzones = _mesh.n_elem();

  /* Set the indices for the first and last real zones. */
  mesh->umesh->firstRealZone = 0;
  mesh->umesh->lastRealZone  = _mesh.n_elem()-1;

  // clear old data
  x.clear();
  y.clear();
  z.clear();
  connectivity.clear();
  global_id_to_index.clear();

  //get node location
  MeshBase::const_node_iterator node_it = _mesh.active_nodes_begin();
  MeshBase::const_node_iterator node_it_end = _mesh.active_nodes_end();
  for(int i=0; node_it!=node_it_end; ++node_it)
  {
    const Node * node = (*node_it);
    x.push_back((*node)(0)/PhysicalUnit::um);
    y.push_back((*node)(1)/PhysicalUnit::um);
    z.push_back((*node)(2)/PhysicalUnit::um);
    global_id_to_index[node->id()] = i++;
  }

  // get cell connectivity information
  MeshBase::const_element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator elem_it_end = _mesh.active_elements_end();
  for(; elem_it!=elem_it_end; ++elem_it)
  {
    if( (*elem_it)->processor_id() != Genius::processor_id() ) continue;

    std::vector<unsigned int> conn;
    (*elem_it)->connectivity(0, VTK, conn);

    connectivity.push_back(visit_cell_type(*elem_it));
    for(unsigned int i=0; i<conn.size(); ++i)
      connectivity.push_back(global_id_to_index[conn[i]]);
  }

  // fill node location and cell connectivity information to visit

  /* Let VisIt use simulation¡¯s copy of the mesh coordinates. */
  mesh->umesh->xcoords = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_SIM, &x[0]);
  mesh->umesh->ycoords = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_SIM, &y[0]);
  mesh->umesh->zcoords = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_SIM, &z[0]);

  /* Let VisIt use the simulation¡¯s copy of the connectivity. */
  mesh->umesh->connectivity = VisIt_CreateDataArrayFromInt(VISIT_OWNER_SIM, &connectivity[0]);
  mesh->umesh->connectivityLen = connectivity.size();


  return mesh;
}



VisIt_ScalarData * VisitHook::get_scalar(int domain, const std::string &name)
{

  size_t sz = sizeof(VisIt_ScalarData);
  VisIt_ScalarData *scalar = (VisIt_ScalarData*)malloc(sz);
  memset(scalar, 0, sz);

  genius_assert( domain < Genius::n_processors() );
  genius_assert( !global_id_to_index.empty() );

  static std::vector<double> solution;
  solution.resize(global_id_to_index.size());

  for( unsigned int r=0; r<_p_solver->get_system().n_regions(); r++)
  {
    const SimulationRegion * region = _p_solver->get_system().region(r);

    SimulationRegion::const_node_iterator node_it     = region->nodes_begin();
    SimulationRegion::const_node_iterator node_it_end = region->nodes_end();
    for(int i=0; node_it!=node_it_end; ++node_it)
    {
      const Node * node = (*node_it).second->root_node();
      const FVM_NodeData * node_data = (*node_it).second->node_data();

      if( !node->on_local() ) continue;

      int index = global_id_to_index[node->id()];

      if(name == "potential")
        solution[index] = (node_data->psi()/PhysicalUnit::V);
      if(name == "electron_density")
        solution[index] = (node_data->n()/std::pow(PhysicalUnit::cm,-3));
      if(name == "hole_density")
        solution[index] = (node_data->p()/std::pow(PhysicalUnit::cm,-3));
    }
  }

  scalar->len  = solution.size();
  scalar->data = VisIt_CreateDataArrayFromDouble(VISIT_OWNER_SIM, &solution[0]);

  return  scalar;

}





#ifndef CYGWIN

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new VisitHook(solver, name, fun_data );
  }
}

#endif




/******************************************************************************
 *  VISIT INTERFACE FUNCTIONS
 *****************************************************************************/

VisIt_SimulationMetaData *VisItGetMetaData(void)
{

  return VisitHook::get_meta_data();
}


VisIt_MeshData *VisItGetMesh(int domain, const char *name)
{
  return VisitHook::get_mesh(domain, name);
}


VisIt_CurveData *VisItGetCurve(const char *name)
{
  return VisitHook::get_electrode_iv(name);
}


VisIt_ScalarData *VisItGetScalar(int domain, const char *name)
{
  return VisitHook::get_scalar(domain, name);
}

VisIt_DomainList *VisItGetDomainList(void)
{
  size_t sz = sizeof(VisIt_DomainList);
  VisIt_DomainList *dl = (VisIt_DomainList*)malloc(sz);
  memset(dl, 0, sz);

  /* Get number of processors and rank from MPI. */
  int np   = Genius::n_processors();
  int rank = Genius::processor_id();

  dl->nTotalDomains = np;
  dl->nMyDomains = 1;
  dl->myDomains = VisIt_CreateDataArrayFromInt(VISIT_OWNER_SIM, &rank);
  return dl;
}


VisIt_SimulationCallback visitCallbacks =
  {
    &VisItGetMetaData,
    &VisItGetMesh,
    NULL,               /* GetMaterial */
    NULL,               /* GetSpecies */
    &VisItGetScalar,    /* GetScalar */
    &VisItGetCurve,     /* GetCurve */
    NULL,               /* GetMixedScalar */
    &VisItGetDomainList /* GetDomainList */
  };



