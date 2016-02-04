#include <vector>
#include <cstdlib>

#include "CogendaHDF5.h"

#ifdef HAVE_HDF5

namespace CogendaHDF5
{

template<>
std::string getAttribute<std::string>(hid_t parent_grp, const std::string& obj, const std::string& name, const std::string& defVal)
{
  hid_t attr, atype;
  attr = H5Aopen_by_name(parent_grp, obj.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
  if (attr<0)
    return defVal;

  char* buf;
  atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, H5T_VARIABLE);

  H5Aread(attr, atype, &buf);
  std::string str(buf);

  free(buf);
  H5Tclose(atype);
  H5Aclose(attr);

  return str;
}

template<>
void setAttribute<std::string>(hid_t obj, const std::string& name, const std::string& val)
{
  hid_t attr, atype, aspace;

  atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, H5T_VARIABLE);

  aspace = H5Screate(H5S_SCALAR);

  attr = H5Acreate_by_name(obj, ".", name.c_str(), atype, aspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  const char* buf = val.c_str();
  H5Awrite(attr, atype, &buf);

  H5Tclose(atype);
  H5Sclose(aspace);
  H5Aclose(attr);
}

typedef std::vector<std::string> str_list;
typedef const std::vector<std::string> const_str_list;

template<>
str_list getAttribute<str_list>(hid_t parent_grp, const std::string& obj, const std::string& name, const_str_list& defVal)
{
  hid_t attr, atype, aspace;
  str_list lst;

  attr = H5Aopen_by_name(parent_grp, obj.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
  if (attr<0)
    return defVal;

  atype = H5Aget_type(attr);
  if (H5Tis_variable_str(atype))
  {
    // variable-length string

    aspace = H5Aget_space(attr);
    hsize_t n = H5Sget_simple_extent_npoints(aspace);

    std::vector<char*> ptrs;
    ptrs.resize(n);

    H5Aread(attr, atype, &ptrs[0]);
    for (std::vector<char*>::const_iterator i=ptrs.begin(); i!=ptrs.end(); i++)
      lst.push_back(std::string(*i));

    H5Dvlen_reclaim(atype, aspace, H5P_DEFAULT, &ptrs[0]);
  }
  else
  {
    // fixed-length string

    aspace = H5Aget_space(attr);
    hsize_t n = H5Sget_simple_extent_npoints(aspace);
    size_t len = H5Tget_size(atype);

    std::vector<char> buf(n*len);
    H5Aread(attr, atype, &buf[0]);

    for (hsize_t i=0; i<n; i++)
    {
      lst.push_back(std::string(&buf[i*n]));
    }
  }

  H5Sclose(aspace);
  H5Tclose(atype);
  H5Aclose(attr);

  return lst;
}

template<>
void setAttribute<str_list>(hid_t obj, const std::string& name, const_str_list& lst)
{
  hid_t attr, atype, aspace;

  atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, H5T_VARIABLE);

  hsize_t dim = lst.size();
  aspace = H5Screate_simple(1, &dim, &dim);

  attr = H5Acreate_by_name(obj, ".", name.c_str(), atype, aspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  std::vector<const char*> ptrs;
  for (const_str_list::const_iterator i=lst.begin(); i!=lst.end(); i++)
      ptrs.push_back(i->c_str());

  H5Awrite(attr, atype, &ptrs[0]);

  H5Tclose(atype);
  H5Sclose(aspace);
  H5Aclose(attr);

}

template<typename T>
std::vector<T> getAttribute(hid_t parent_grp, const std::string& obj, const std::string& name, const std::vector<T>& defVal=T())
{
  hid_t attr, atype, aspace;
  attr = H5Aopen_by_name(parent_grp, obj.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
  if (attr<0)
    return defVal;

  T data;
  atype = H5Tcopy(getHDF5MemType<T>());

  aspace = H5Aget_space(attr);
  hsize_t n = H5Sget_simple_extent_npoints(aspace);

  std::vector<T> lst;
  H5Aread(attr, atype, &lst[0]);

  H5Tclose(atype);
  H5Sclose(aspace);
  H5Aclose(attr);

  return data;
}

template<typename T>
void setAttribute(hid_t obj, const std::string& name, const std::vector<T>& lst)
{
  hid_t attr, atype, aspace;

  atype = H5Tcopy(getHDF5MemType<T>());

  hsize_t dim = lst.size();
  aspace = H5Screate_simple(1, &dim, &dim);

  attr = H5Acreate_by_name(obj, ".", name.c_str(), atype, aspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Awrite(attr, atype, &lst[0]);

  H5Tclose(atype);
  H5Sclose(aspace);
  H5Aclose(attr);
}

template<>
LookupTable getAttribute(hid_t parent_grp, const std::string &obj, const std::string &name, const LookupTable& defVal)
{
  return LookupTable(getAttribute<str_list>(parent_grp, obj, name, defVal()));
}

template<>
void setAttribute(hid_t obj, const std::string &name, const LookupTable& val)
{
  setAttribute<str_list>(obj, name, val());
}

template<typename T>
T getAttribute(hid_t parent_grp, const std::string& obj, const std::string& name, const T& defVal)
{
  hid_t attr, atype;
  attr = H5Aopen_by_name(parent_grp, obj.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
  if (attr<0)
    return defVal;

  T data;
  atype = H5Tcopy(getHDF5MemType<T>());

  H5Aread(attr, atype, &data);

  H5Tclose(atype);
  H5Aclose(attr);

  return data;
}

template<typename T>
void setAttribute(hid_t obj, const std::string& name, const T& val)
{
  hid_t attr, atype, aspace;

  atype = H5Tcopy(getHDF5MemType<T>());
  aspace = H5Screate(H5S_SCALAR);

  attr = H5Acreate_by_name(obj, ".", name.c_str(), atype, aspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Awrite(attr, atype, &val);

  H5Tclose(atype);
  H5Sclose(aspace);
  H5Aclose(attr);
}

template float          getAttribute(hid_t, const std::string&, const std::string&, const float&);
template double         getAttribute(hid_t, const std::string&, const std::string&, const double&);
template char           getAttribute(hid_t, const std::string&, const std::string&, const char&);
template unsigned char  getAttribute(hid_t, const std::string&, const std::string&, const unsigned char&);
template int            getAttribute(hid_t, const std::string&, const std::string&, const int&);
template unsigned int   getAttribute(hid_t, const std::string&, const std::string&, const unsigned int&);
template long           getAttribute(hid_t, const std::string&, const std::string&, const long&);
template unsigned long  getAttribute(hid_t, const std::string&, const std::string&, const unsigned long&);
template long long           getAttribute(hid_t, const std::string&, const std::string&, const long long&);
template unsigned long long  getAttribute(hid_t, const std::string&, const std::string&, const unsigned long long&);

template void setAttribute(hid_t obj, const std::string& name, const float& val);
template void setAttribute(hid_t obj, const std::string& name, const double& val);
template void setAttribute(hid_t obj, const std::string& name, const char& val);
template void setAttribute(hid_t obj, const std::string& name, const unsigned char& val);
template void setAttribute(hid_t obj, const std::string& name, const int& val);
template void setAttribute(hid_t obj, const std::string& name, const unsigned int& val);
template void setAttribute(hid_t obj, const std::string& name, const long& val);
template void setAttribute(hid_t obj, const std::string& name, const unsigned long& val);
template void setAttribute(hid_t obj, const std::string& name, const long long& val);
template void setAttribute(hid_t obj, const std::string& name, const unsigned long long& val);



}

#endif


