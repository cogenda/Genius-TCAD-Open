%Module(name=genius, version=0)

%Include types.sip

%Include common.sip
%Include simulation_system.sip
%Include material.sip
%Include control.sip

%ModuleCode
#include <cstring>
#include "control.h"

void sipDelayedDtors(const sipDelayedDtor *dd_list)
{
  const sipDelayedDtor *dd = dd_list;
  for(;dd; dd=dd->dd_next)
  {
    void* ptr;
    ptr = dd->dd_ptr;

    //FIXME: SIP documents doesn't seem to be accurate here
    //need to study further
    //if (dd->dd_isderived)
    //  ptr = ((sipSimpleWrapper*)ptr)->u.cppPtr;

    if (strncmp(dd->dd_name, "SolverControl", 20)==0)
    {
      delete reinterpret_cast<SolverControl *>(ptr);
    }
  }
}

%End

