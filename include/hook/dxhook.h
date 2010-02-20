// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: dxhook.h,v 1.1 2008/06/09 05:52:30 gdiso Exp $

#ifndef DXHOOK_H
#define DXHOOK_H

#ifdef USE_PTHREADS
#include <pthread.h>
#endif
#ifdef USE_SSL
#include <SSL/sockets.h>
#endif
#include <src/dvector.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class FieldGen {
public:
  virtual ~FieldGen() {}
  virtual void init(int ndof,TextHashTable* options,char *name)=0;
  virtual int n()=0;
  virtual void field(int j,string &name,vector<int> &rank)=0;
  virtual void values(int j,vector<double> &in,vector<double> &out)=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class FieldGenList : public vector<FieldGen *> {
public:
  ~FieldGenList();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This hooks is in charge of sending node coordinates, element
    conectivities and results to the DX client. */ 
class dx_hook : public Hook {
#ifdef USE_SSL
private:
  /// Table of options
  TextHashTableFilter *options;
  /// The srvr_root is created first and then a srvr is established. 
  Socket *srvr_root,*srvr;
  /// The FEM mesh
  Mesh *mesh;
  /// The dofmap
  Dofmap *dofmap;
  /// Auxiliary variables
  int step_cntr, steps, ierr, dx_auto_combine;
  /// A list of fields
  FieldGenList field_gen_list;
  /// The name of the state file to be read
  string state_file;
  /** the record in the file (many states may be stored
      in the same file */
  int record;
  /// Integer parameters
  int ndim,nnod,ndof,nu;
  /** Flags whether the mesh changes and 
      in which file it is read */
  string dx_node_coordinates;
  double coef0, coef;
  int read_coords;
  int dx_do_make_command;
  dvector<double> x0;
  // The step solicited by DX
  int dx_step;

  /// Flags reading states from files
  int dx_read_state_from_file;
  /// Flags error condition when can't reading state
  /// from file. 
  int dx_stop_on_bad_file;
  /// Flags sending coordinates each step or using cache
  int dx_cache_coords;

#ifdef USE_PTHREADS
  enum connection_state_t {
    not_launched, not_connected, connected} connection_state_m,
    connection_state_master;
  void re_launch_connection();
  pthread_t thread;
  connection_state_t connection_state();
  void set_connection_state(connection_state_t s);
#endif

public:
  dx_hook();
  ~dx_hook() { delete options; }
  void init(Mesh &mesh,Dofmap &dofmap,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
  void *wait_connection();
  virtual Vec state()=0;
  virtual TimeData *time_data()=0;
  typedef int (dx_hook::*build_state_fun_t)(double *);
  int build_state_from_state(double *);
  int build_state_from_file(double *);
  void send_state(int step,build_state_fun_t bf);
#else
public:
  void init(Mesh &mesh,Dofmap &dofmap,const char *name) {
    PETSCFEM_ERROR0("Hook unavailable. Code not compiled with sockets library!!");  
    PetscFinalize();
    exit(0);
  }
#endif
};

#endif
