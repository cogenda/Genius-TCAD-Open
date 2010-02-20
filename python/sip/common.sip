
// {{{ namespace Genius 
namespace Genius {

%TypeHeaderCode
#include "log.h"
#include <iostream>
#include <string>
#include <map>
#include <cstdlib>
#include "genius_env.h"
#include "genius_common.h"
int argc;
char **argv;

#ifndef WIN32
#include <ext/stdio_filebuf.h>
typedef __gnu_cxx::stdio_filebuf<char> cfilebuf_t;
#else
#include <fstream>
typedef std::basic_filebuf<char> cfilebuf_t;
#endif
std::map<std::string, cfilebuf_t*> filebufs;

int _genius_error=0;

static PetscErrorCode genius_error_handler(int line, const char *func, const char *file, const char *dir,PetscErrorCode n,int p, const char *mess,void *ctx)
{

  MESSAGE << "--------------------- Error Message ------------------------------------" << std::endl;
  MESSAGE << "Fatal Error:";
  if (mess)
    MESSAGE << mess;
  MESSAGE << " at line " << line << " in " << file << std::endl;
  MESSAGE << "------------------------------------------------------------------------" << std::endl << std::endl << std::endl;
  RECORD();

  _genius_error=1;
#ifdef WIN32
  _sleep(86400);
#else
  sleep(86400);
#endif
  exit(-1);
}

%End

  bool init_processors(SIP_PYLIST);
%MethodCode
    argc = PyList_Size(a0);
    argv = new char*[argc];
    for (int i=0; i<argc; i++)
    {
      PyObject *strObj = PyList_GetItem(a0, i);
      argv[i] = PyString_AsString( strObj );
    }

    Py_BEGIN_ALLOW_THREADS
    sipRes = Genius::init_processors(&argc, &argv);

    PetscPushErrorHandler(genius_error_handler, NULL);
    
    Py_END_ALLOW_THREADS
%End

  bool clean_processors();
%MethodCode
    Py_BEGIN_ALLOW_THREADS
    sipRes = Genius::clean_processors();
    Py_END_ALLOW_THREADS
    delete argv;
%End

  unsigned int n_processors();
  unsigned int processor_id();
  const char * input_file();
  void set_input_file(const char* fname);
  std::string genius_dir();
  void set_genius_dir(const std::string &genius_dir);

  void add_log_file(const std::string &name, SIP_PYOBJECT);
%MethodCode
    if (PyFile_Check(a1))
    {
#ifndef WIN32
      cfilebuf_t  *ofbuf =  new cfilebuf_t(PyFile_AsFile(a1), std::ios_base::out);
#else
      cfilebuf_t  *ofbuf =  new cfilebuf_t(stderr);
#endif
      filebufs.insert(std::pair<std::string, cfilebuf_t*>(*a0, ofbuf));
      genius_log.addStream(*a0, ofbuf);
    }
%End

  void add_log_filename(const std::string &name, const std::string &fname);
%MethodCode
    genius_log.addStream(*a0, *a1);
%End

  void remove_log_file(const std::string name);
%MethodCode
  genius_log.removeStream(*a0);
  std::map<std::string, cfilebuf_t*>::iterator it = filebufs.find(*a0);
  if (it != filebufs.end())
  {
    delete it->second;
    filebufs.erase(it);
  }
%End

  void log(const char *msg /in/);
%MethodCode
    Py_BEGIN_ALLOW_THREADS
    MESSAGE << a0;
    RECORD();
    Py_END_ALLOW_THREADS
%End

  int error();
%MethodCode
    sipRes = _genius_error;
%End
};

// }}}

// {{{ namespace Parallel
namespace Parallel
{
%TypeHeaderCode
#include "parallel.h"
#include <string>
#include <iostream>
%End

SIP_PYOBJECT broadcast(SIP_PYOBJECT, const unsigned int root);
%MethodCode
  int type = 0;
  bool iamroot = false;

  if (Genius::n_processors() == 1)
  {
    std::string str = PyString_AsString(a0);
    sipRes = PyString_FromString(str.c_str());
  }
  else
  {

    if (Genius::processor_id() == a1)
    {
      iamroot = true;
      if (PyString_Check(a0)) type = 1;        // string
      else if (PyInt_Check(a0)) type = 2;      // int
      else if (PyFloat_Check(a0)) type = 3;    // float
    }
    Parallel::broadcast(type, a1);

    switch(type)
    {
    case 1: // string
    { 
      std::string str;
      if (iamroot) str = PyString_AsString(a0);
      Parallel::broadcast(str, a1);
      sipRes = PyString_FromString(str.c_str()); 
      break;
    }
    case 2: // integer
    {
      long val;
      if (iamroot) val = PyInt_AsLong(a0);
      Parallel::broadcast(val, a1);
      sipRes = PyInt_FromLong(val);
      break;
    }
    case 3: // float
    {
      double val;
      if (iamroot) val = PyFloat_AsDouble(a0);
      Parallel::broadcast(val, a1);
      sipRes = PyFloat_FromDouble(val);
      break;
    }
    default:
      PyErr_SetString(PyExc_TypeError, "Type not supported");
      sipIsErr = 1;
    }
  }
%End

};
// }}}


