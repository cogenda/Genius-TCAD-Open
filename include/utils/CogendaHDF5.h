#ifndef COGENDAHDF5_H
#define COGENDAHDF5_H

#include "config.h"

#ifdef HAVE_HDF5

#include <string>
#include <vector>
#include <map>
#include "hdf5.h"

/**
 * Everything related to HDF5 is in this name space
 */

namespace CogendaHDF5
{

  static const std::string ATTR_fullname = "meta-fullname";
  static const std::string ATTR_version  = "meta-version";

  /**
   * @returns the HDF5 type constant (for in-memory type) corresponding to the template type.
   */
  template<typename T> inline hid_t getHDF5MemType();
  template<> inline hid_t getHDF5MemType<double> () { return H5T_NATIVE_DOUBLE; }
  template<> inline hid_t getHDF5MemType<float> () { return H5T_NATIVE_FLOAT; }
  template<> inline hid_t getHDF5MemType<char> () { return H5T_NATIVE_CHAR; }
  template<> inline hid_t getHDF5MemType<unsigned char> () { return H5T_NATIVE_UCHAR; }
  template<> inline hid_t getHDF5MemType<int> () { return H5T_NATIVE_INT; }
  template<> inline hid_t getHDF5MemType<unsigned int> () { return H5T_NATIVE_UINT; }
  template<> inline hid_t getHDF5MemType<long> () { return H5T_NATIVE_LONG; }
  template<> inline hid_t getHDF5MemType<unsigned long> () { return H5T_NATIVE_ULONG; }
  template<> inline hid_t getHDF5MemType<long long> () { return H5T_NATIVE_LLONG; }
  template<> inline hid_t getHDF5MemType<unsigned long long> () { return H5T_NATIVE_ULLONG; }

  /**
   * @returns the attribute \p name of the object \p <parent_grp>/<obj> .
   * @param parent_grp  id of the parent group containing the object
   * @param obj         name of the object
   * @param name        name of the attribute
   * @param defVal      default value if the attribute does not exist
   */
  template<typename T>
  T getAttribute(hid_t parent_grp, const std::string& obj, const std::string& name, const T& defVal=T());

  /**
   * set the attribute \p name of the object \p obj
   * @param obj   id of the object
   * @param name  name of the attribute
   * @param val   value
   */
  template<typename T>
  void setAttribute(hid_t obj, const std::string& name, const T& val);


  /**
   * Attribute that contains a string Lookup table
   */
  class LookupTable
  {
  public:
    /**
     * Constructor. build the LUP from a list of strings.
     */
    LookupTable(const std::vector<std::string>& vals=std::vector<std::string>())
      : _lst(vals)
    {
      int i=0;
      for (std::vector<std::string>::const_iterator it=_lst.begin(); it!=_lst.end(); it++)
        _map.insert(std::make_pair(*it, i++));
    }

    /**
     * @returns the string at index i
     */
    inline std::string operator[](int i) const { return _lst[i]; }

    /**
     * @returns the entire table as string list
     */
    inline std::vector<std::string> operator()() const { return _lst; }

    /**
     * @returns the index for the string val; -1 if not found.
     */
    inline int find(const std::string& val) const
    {
      std::map<std::string, int>::const_iterator it= _map.find(val);
      if (it==_map.end())
        return -1;
      return it->second;
    }

    /**
     * append an item the LUP
     */
    inline int append(const std::string& val)
    {
      if(_map.find(val) != _map.end())
        return _map.find(val)->second;

      int i=_lst.size();
      _lst.push_back(val);
      _map.insert(std::make_pair(val, i));
      return i;
    }

    /**
     * append an item in stream-style
     */
    LookupTable& operator<<(const std::string& val)
    {
      append(val);
      return *this;
    }

  private:
    std::vector<std::string> _lst;
    std::map<std::string, int> _map;
  };

  template <class RecordStruct>
  class RecordDataset
  {
  public:

    virtual ~RecordDataset() {}

    /**
     * @returns the number of data records, i.e. number of tracks
     */
    inline size_t size() const { return _data.size(); }

    /**
     * @returns the \p i-th data record
     */
    inline const RecordStruct& operator[](size_t i) const {return _data[i]; }
    /**
     * @returns the \p i-th data record
     */
    inline RecordStruct& operator[](size_t i) {return _data[i]; }

    /**
     * @returns the last data record
     */
    inline const RecordStruct& back() const { return _data.back(); }
    /**
     * @returns the last data record
     */
    inline RecordStruct& back() { return _data.back(); }

    /**
     * Append a record, and returns the newly added, empty record.
     */
    RecordStruct& append() { _data.push_back(filler()); return _data.back(); }

    /**
     * read data from the hdf5 file.
     * @param grp   id of the group containing the dataset. Tracks will be read from <grp>/dName.
     * @returns true if success, false otherwise
     */
    bool readData(hid_t grp, const std::string& dName)
    {
      hid_t dset, dspace;

      dset = H5Dopen(grp, dName.c_str(), H5P_DEFAULT);

      if (dset<0)    return false;

      dspace = H5Dget_space(dset);
      if (dspace<0)
      {
        H5Dclose(dset);
        return false;
      }
      hsize_t nelmts = H5Sget_simple_extent_npoints(dspace);

      _data.resize(nelmts, filler());

      H5Dread(dset, memDataType(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &_data[0]);

      H5Sclose(dspace);
      H5Dclose(dset);
      return true;
    }

    /**
     * write data to the hdf5 file.
     * @param grp   id of the group containing the dataset. Tracks will be read from <grp>/dName.
     * @returns true if success, false otherwise
     */
    bool writeData(hid_t grp, const std::string& dName)
    {
      hid_t dset, dspace, plist;

      const hsize_t dims[1] = {size()};
      dspace = H5Screate_simple(1, dims, NULL);

      plist = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_chunk(plist, 1, dims);
      H5Pset_shuffle(plist);
      H5Pset_deflate(plist, 6);

      dset = H5Dcreate(grp, dName.c_str(), diskDataType(), dspace, H5P_DEFAULT, plist, H5P_DEFAULT);
      if (dset<0)
      {
        H5Sclose(dspace);
        H5Pclose(plist);
        return false;
      }

      H5Dwrite(dset, memDataType(), H5S_ALL, dspace, H5P_DEFAULT, &_data[0]);

      H5Dclose(dset);
      H5Sclose(dspace);
      H5Pclose(plist);

      return true;
    }

  protected:
    virtual hid_t memDataType() const = 0;
    virtual hid_t diskDataType() const = 0;
    virtual RecordStruct filler() const = 0;

  private:
    std::vector<RecordStruct> _data;
  };

}

#endif //HAVE_HDF5

#endif // COGENDAHDF5_H
