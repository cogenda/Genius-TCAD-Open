// {{{ %MappedType std::string
%MappedType std::string
{
%TypeHeaderCode
#include <string>
%End

%ConvertFromTypeCode
  if (!sipCpp)
  {
    Py_INCREF (Py_None);
    return Py_None;
  }
  return PyString_FromString(sipCpp->c_str());
%End

%ConvertToTypeCode
  if (sipIsErr == NULL)
    return (PyString_Check(sipPy));

  if (sipPy == Py_None)
  {
     *sipCppPtr = new std::string();
     return 1;
  }
  if (PyString_Check(sipPy))
  {
     *sipCppPtr = new std::string(PyString_AS_STRING(sipPy));
     return 1;
  }
  return 0;
%End

};
// }}}


