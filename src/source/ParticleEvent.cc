#include <cstring>
#include <iostream>


#include "ParticleEvent.h"
#include "CogendaHDF5.h"

#ifdef HAVE_HDF5

using namespace CogendaHDF5::GSeat;


const std::string Event::fullname = "com.cogenda.hdf5.gseat.ParticleEvent";
const std::string Event::version  = "0.2";



Step::Dataset::Dataset()
{
  Record f = {{0.0, 0.0, 0.0}, 0, 0.0, 0.0};
  _filler = f;
  _mkMemDataType();
  _mkDiskDataType();
}

void Step::Dataset::_mkMemDataType()
{
  hid_t dtEnd;

  _memDataType = H5Tcreate(H5T_COMPOUND, sizeof(Record));

  const hsize_t dims[] = {3};
  dtEnd = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dims);


  H5Tinsert(_memDataType, "End",           HOFFSET(Record, EndPoint), dtEnd);
  H5Tinsert(_memDataType, "Region",        HOFFSET(Record, Region),   H5T_NATIVE_UINT);
  //H5Tinsert(_memDataType, "Material",      HOFFSET(Record, Material), H5T_NATIVE_UINT);
  //H5Tinsert(_memDataType, "Density",       HOFFSET(Record, Density),  H5T_NATIVE_DOUBLE);
  H5Tinsert(_memDataType, "EnergyDeposit", HOFFSET(Record, EnergyDeposit), H5T_NATIVE_DOUBLE);
  H5Tinsert(_memDataType, "ehPairRadius",  HOFFSET(Record, ehPairRadius),  H5T_NATIVE_DOUBLE);

  H5Tclose(dtEnd);
}

void Step::Dataset::_mkDiskDataType()
{
  hid_t dtEnd;

  _diskDataType = H5Tcreate(H5T_COMPOUND, sizeof(Record));

  const hsize_t dims[] = {3};
  dtEnd = H5Tarray_create(H5T_IEEE_F64LE, 1, dims);

  H5Tinsert(_diskDataType, "End",           HOFFSET(Record, EndPoint), dtEnd);
  H5Tinsert(_diskDataType, "Region",        HOFFSET(Record, Region),   H5T_STD_U32LE);
  //H5Tinsert(_diskDataType, "Material",      HOFFSET(Record, Material), H5T_STD_U32LE);
  //H5Tinsert(_diskDataType, "Density",       HOFFSET(Record, Density),  H5T_IEEE_F64LE);
  H5Tinsert(_diskDataType, "EnergyDeposit", HOFFSET(Record, EnergyDeposit), H5T_IEEE_F64LE);
  H5Tinsert(_diskDataType, "ehPairRadius",  HOFFSET(Record, ehPairRadius),  H5T_IEEE_F64LE);
  H5Tpack(_diskDataType);

  H5Tclose(dtEnd);
}

Step::Dataset::~Dataset()
{
  H5Tclose(_memDataType);
  H5Tclose(_diskDataType);
}

Track::Dataset::Dataset()
{
  _mkMemDataType();
  _mkDiskDataType();

  Record f = {0, "", 0, 0.0, 0.0, 0, {0.0, 0.0, 0.0}, 0, 0};
  _filler = f;
}

void Track::Dataset::_mkMemDataType()
{
  hid_t dtLabel, dtPoint;

  _memDataType = H5Tcreate(H5T_COMPOUND, sizeof(Record));

  dtLabel = H5Tcopy (H5T_C_S1);
  H5Tset_size (dtLabel, lenStringLabel+1);  // reserve one more byte in memory, for NULL terminator.

  const hsize_t dims[] = {3};
  dtPoint = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dims);

  H5Tinsert(_memDataType, "ID",            HOFFSET(Record, ID),            H5T_NATIVE_ULONG);
  H5Tinsert(_memDataType, "Particle",      HOFFSET(Record, Particle),      dtLabel);
  H5Tinsert(_memDataType, "ParentID",      HOFFSET(Record, ParentID),      H5T_NATIVE_ULONG);
  H5Tinsert(_memDataType, "InitialEnergy", HOFFSET(Record, InitialEnergy), H5T_NATIVE_DOUBLE);
  H5Tinsert(_memDataType, "EnergyDeposit", HOFFSET(Record, EnergyDeposit), H5T_NATIVE_DOUBLE);
  H5Tinsert(_memDataType, "Charge",        HOFFSET(Record, Charge),        H5T_NATIVE_INT);
  H5Tinsert(_memDataType, "Start",         HOFFSET(Record, StartPoint),    dtPoint);
  H5Tinsert(_memDataType, "StepOffset",    HOFFSET(Record, StepOffset),    H5T_NATIVE_HSIZE);
  H5Tinsert(_memDataType, "NumSteps",      HOFFSET(Record, NumSteps),      H5T_NATIVE_HSIZE);

  H5Tclose(dtLabel);
  H5Tclose(dtPoint);
}

void Track::Dataset::_mkDiskDataType()
{
  hid_t dtLabel, dtPoint;

  _diskDataType = H5Tcreate(H5T_COMPOUND, sizeof(Record));

  dtLabel = H5Tcopy (H5T_C_S1);
  H5Tset_size (dtLabel, lenStringLabel);

  const hsize_t dims[] = {3};
  dtPoint = H5Tarray_create(H5T_IEEE_F64LE, 1, dims);

  H5Tinsert(_diskDataType, "ID",            HOFFSET(Record, ID),            H5T_STD_U64LE);
  H5Tinsert(_diskDataType, "Particle",      HOFFSET(Record, Particle),      dtLabel);
  H5Tinsert(_diskDataType, "ParentID",      HOFFSET(Record, ParentID),      H5T_STD_U64LE);
  H5Tinsert(_diskDataType, "InitialEnergy", HOFFSET(Record, InitialEnergy), H5T_IEEE_F64LE);
  H5Tinsert(_diskDataType, "EnergyDeposit", HOFFSET(Record, EnergyDeposit), H5T_IEEE_F64LE);
  H5Tinsert(_diskDataType, "Charge",        HOFFSET(Record, Charge),        H5T_STD_I32LE);
  H5Tinsert(_diskDataType, "Start",         HOFFSET(Record, StartPoint),    dtPoint);
  H5Tinsert(_diskDataType, "StepOffset",    HOFFSET(Record, StepOffset),    H5T_STD_U64LE);
  H5Tinsert(_diskDataType, "NumSteps",      HOFFSET(Record, NumSteps),      H5T_STD_U64LE);

  H5Tpack(_diskDataType);

  H5Tclose(dtLabel);
  H5Tclose(dtPoint);
}


Track::Dataset::~Dataset()
{
  H5Tclose(_memDataType);
}


bool Event::readData(hid_t parent_grp, const std::string &name)
{
  if (fullname!=CogendaHDF5::getAttribute<std::string>(parent_grp, name, ATTR_fullname)) return false;

  hid_t grp;

  grp = H5Gopen(parent_grp, name.c_str(), H5P_DEFAULT);
  if (grp<0)
    return false;

  bool ret = true;

  setID(CogendaHDF5::getAttribute(grp, ".", "ID", 0));
  setEnergyDeposit(CogendaHDF5::getAttribute(grp, ".", "EnergyDeposit", 0.0));

  ret = ret && _tracks.readData(grp, "Tracks");
  ret = ret && _steps.readData(grp, "Steps");

  H5Gclose(grp);
  return ret;
}

bool Event::writeData(hid_t parent_grp, const std::string &name)
{
  hid_t grp;

  grp = H5Gcreate(parent_grp, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (grp<0)
    return false;

  CogendaHDF5::setAttribute(grp, ATTR_fullname, fullname);
  CogendaHDF5::setAttribute(grp, ATTR_version,  version);
  CogendaHDF5::setAttribute(grp, "ID", ID());
  CogendaHDF5::setAttribute(grp, "EnergyDeposit", EnergyDeposit());

  bool ret = true;

  ret = ret && _tracks.writeData(grp, "Tracks");
  ret = ret && _steps.writeData(grp, "Steps");

  H5Gclose(grp);
  return ret;
}

#endif

