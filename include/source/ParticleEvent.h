#ifndef PARTICLEEVENT_H
#define PARTICLEEVENT_H

#include <limits>
#include <cstring>
#include <vector>

#include "CogendaHDF5.h"

#ifdef HAVE_HDF5


namespace CogendaHDF5
{

namespace GSeat
{

static const size_t lenStringLabel = 32;

/**
 * Interface class for users to access the properties of a Step.
 * No data record is actually stored here. We only keep a reference
 * to the corresponding record in DatasetStep.
 * User typically discard a Step object after using it.
 */
class Step
{
public:
  /**
   * data for one step.
   */
  struct Record
  {
    double       EndPoint[3];
    unsigned int Region;        /// Region ID (see Event)
    //unsigned int Material;      /// Material ID (see Event)
    //double       Density;
    double       EnergyDeposit;
    double       ehPairRadius;
  };
  /**
   * Container holding the data of all steps in an Event.
   * Provides \p size() and \p [] operator for accessing records.
   * Typical user should not have to use this class.
   */
  class Dataset : public RecordDataset<Record>
  {
  public:
    Dataset();
    virtual ~Dataset();

  protected:
    virtual hid_t memDataType() const {return _memDataType;}
    virtual hid_t diskDataType() const {return _diskDataType;}
    virtual Record filler() const { return _filler; }

  private:
    void _mkMemDataType();
    void _mkDiskDataType();

    hid_t _memDataType, _diskDataType;
    Record _filler;
  };

public:
  /**
   * Constructor. User typically do not call this directly.
   * Instead, one should obtain an instance from one of the
   * iterator/appender provided by the Event class.
   */
  Step(Record &rec) : _rec(rec) {}

  /**
   * return the coordinates of the end point
   */
  inline std::vector<double> EndPoint() const
  {
    std::vector<double> r;
    r.push_back(_rec.EndPoint[0]);
    r.push_back(_rec.EndPoint[1]);
    r.push_back(_rec.EndPoint[2]);
    return r;
  }
  /**
   * return the individual coordinate of the end point
   */
  inline double EndPointX() const { return _rec.EndPoint[0]; }
  inline double EndPointY() const { return _rec.EndPoint[1]; }
  inline double EndPointZ() const { return _rec.EndPoint[2]; }
  /**
   * set the end point coordinates
   */
  inline void setEndPoint(double x, double y, double z)
  {
    _rec.EndPoint[0] = x;
    _rec.EndPoint[1] = y;
    _rec.EndPoint[2] = z;
  }

  /**
   * get region id
   */
  inline unsigned int Region() const { return _rec.Region; }
  /**
   * set region id
   */
  inline void setRegion(unsigned int rgn) { _rec.Region = rgn ; }

  /**
   * get energy deposit in this step
   */
  inline double EnergyDeposit() const { return _rec.EnergyDeposit; }
  /**
   * set energy deposit in this step
   */
  inline void setEnergyDeposit(double e) { _rec.EnergyDeposit = e; }

  inline double ehPairRadius() const { return _rec.ehPairRadius; }
  inline void setehPairRadius(double r) { _rec.ehPairRadius = r; }

private:
  Record &_rec;
};

/**
 * Interface class for users to access the properties of a Track.
 * No data record is actually stored here. We only keep a reference
 * to the corresponding record in DatasetTrack.
 * User typically discard a Track object after using it.
 */
class Track
{
public:
  /**
   * data for one track.
   */
  struct Record
  {
    unsigned long ID;
    char      Particle[lenStringLabel+1];
    unsigned long ParentID;
    double    InitialEnergy;
    double    EnergyDeposit;
    int       Charge;
    double    StartPoint[3];
    hsize_t   StepOffset;
    hsize_t   NumSteps;
  };
  /**
   * Container holding the data of all tracks in an Event.
   * Provides \p size() and \p [] operator for accessing records.
   * Typical user should not have to use this class.
   */
  class Dataset : public RecordDataset<Record>
  {
  public:
    Dataset();
    ~Dataset();
  protected:
    virtual hid_t memDataType() const {return _memDataType;}
    virtual hid_t diskDataType() const {return _diskDataType;}
    virtual Record filler() const { return _filler; }
  private:
    void _mkMemDataType();
    void _mkDiskDataType();

    hid_t _memDataType, _diskDataType;
    Record _filler;
  };

public:
  /**
   * Constructor. User typically do not call this directly.
   * Instead, one should obtain an instance from one of the
   * iterator/appender provided by the Event class.
   *
   * Since steps of all tracks are stored consequtively,
   * one must be careful with the StepOffset and NumSteps properties.
   * Steps for a track are stored in the index range [ StepOffset ... StepOffset+NumSteps-1 ]
   */
  Track(Record &rec) : _rec(rec) {}

  inline unsigned long ID() const { return _rec.ID; }
  inline void setID(unsigned long id) { _rec.ID = id; }

  inline std::string Particle() const { return std::string(_rec.Particle); }
  inline void setParticle(const std::string& pat) { strncpy(_rec.Particle, pat.c_str(), lenStringLabel+1); }

  inline unsigned long ParentID() const { return _rec.ParentID; }
  inline void setParentID(unsigned long pid) { _rec.ParentID = pid; }

  inline double InitialEnergy() const { return _rec.InitialEnergy; }
  inline void setInitialEnergy(double e) { _rec.InitialEnergy = e; }

  inline double EnergyDeposit() const { return _rec.EnergyDeposit; }
  inline void setEnergyDeposit(double e) { _rec.EnergyDeposit = e; }

  inline int Charge() const { return _rec.Charge; }
  inline void setCharge(int c) { _rec.Charge = c; }

  inline std::vector<double> StartPoint() const
  {
    std::vector<double> r;
    r.push_back(_rec.StartPoint[0]);
    r.push_back(_rec.StartPoint[1]);
    r.push_back(_rec.StartPoint[2]);
    return r;
  }
  inline double StartPointX() const { return _rec.StartPoint[0]; }
  inline double StartPointY() const { return _rec.StartPoint[1]; }
  inline double StartPointZ() const { return _rec.StartPoint[2]; }
  inline void setStartPoint(double x, double y, double z)
  {
    _rec.StartPoint[0] = x;
    _rec.StartPoint[1] = y;
    _rec.StartPoint[2] = z;
  }

  inline size_t StepOffset() const { return _rec.StepOffset; }
  inline void setStepOffset(size_t offset) { _rec.StepOffset = offset; }

  inline size_t NumSteps() const { return _rec.NumSteps; }
  inline void setNumSteps(size_t n) { _rec.NumSteps = n; }

private:
  Record &_rec;
};

/**
 * A particle event, consists of many tracks, one for each particle.
 * The 0th track is the root track, for recording the incidient position.
 * Each track consists of many steps.
 * Users should access the tracks/steps with one of the iterators.
 * The tracks/steps should be added with the appenders.
 */
class Event
{
public:
  /**
   * constructor.
   */
  Event(const CogendaHDF5::LookupTable &region_lut, const std::vector<std::string> &materials)
  :  _regions(region_lut), _materials(materials),_energyDeposit(0)
  {}


  /**
   * to be written to meta-fullname attribute
   */
  static const std::string fullname;

  /**
   * to be written to meta-version attribute
   */
  static const std::string version;

  /**
   * Iterator for accessing tracks.
   * User can use Event::track_begin() to obtain an iterator.
   */
  class TrackIter
  {
  protected:
    friend class Event;
    TrackIter(Event& evt, size_t i=0)
      : _base(evt), _offset(i)
    {}

  public:
    /**
     * @returns a copy of Track instance for the current track.
     * The copy points to the actual track record in the dataset,
     * so modifying the copy changes the dataset.
     */
    Track operator*()
    {
      return Track(_base._tracks[_offset]);
    }

    /**
     * move to the next track.
     */
    inline TrackIter& operator++()
    {
      _offset++;
      if (_offset>=_base._tracks.size())
        _offset=std::numeric_limits<size_t>::max();
      return *this;
    }

    /**
     * move to the previous track.
     */
    inline TrackIter& operator--()
    {
      if (_offset==0)
        _offset=std::numeric_limits<size_t>::max();
      else
        _offset--;
      return *this;
    }

    /**
     * check if two iterators equal. TODO: we should also check if they belong to the same Event.
     */
    inline bool operator==(const TrackIter& o) const { return _offset==o._offset; }

    /**
     * check if two iterators differ.
     */
    inline bool operator!=(const TrackIter& o) const { return !(o==*this); }
  private:
    Event& _base;     /// underlying event object
    size_t _offset;   /// offset in the Event._tracks container.
  };

  /**
   * Iterator for accessing steps of a track.
   */
  class StepIter
  {
  protected:
    friend class Event;
    StepIter(Event& evt, const Track& track, size_t i=0)
      : _base(evt), _track(track), _offset(i) {}

  public:
    /**
     * @returns a copy of Step instance for the current step.
     * The copy points to the actual step record in the dataset,
     * so modifying the copy changes the dataset.
     */
    Step operator*()
    {
      size_t i = _track.StepOffset()+_offset;
      return Step(_base._steps[i]);
    }

    /**
     * move to the next track.
     */
    inline StepIter& operator++()
    {
      _offset++;
      if (_offset>=_track.NumSteps())
        _offset=std::numeric_limits<size_t>::max();
      return *this;
    }

    /**
     * move to the previous track.
     */
    inline StepIter& operator--()
    {
      if (_offset==0)
        _offset=std::numeric_limits<size_t>::max();
      else
        _offset--;

      return *this;
    }

    /**
     * check if two iterators equal. TODO: we should also check if they belong to the same Event and same Track.
     */
    inline bool operator==(const StepIter& o) const { return _offset==o._offset; }

    /**
     * check if two iterators differ.
     */
    inline bool operator!=(const StepIter& o) const { return !(o==*this); }
  private:
    Event& _base;       /// underlying event object
    Track _track;       /// underlying track object
    size_t _offset;     /// Track.StepOffset + _offset is the actual Event._steps container
  };


  inline Track appendTrack()
  {
    Track track(_tracks.append());
    track.setStepOffset(_steps.size());
    return track;
  }

  inline Step appendStep()
  {
    Track track(_tracks.back());
    track.setNumSteps(track.NumSteps()+1);
    return Step(_steps.append());
  }

  /**
   * @returns an iterator, pointing at the first track
   */
  TrackIter track_begin() { return TrackIter(*this, 0); }
  /**
   * @returns an iterator that signals the end-condition
   */
  TrackIter track_end() { return TrackIter(*this, std::numeric_limits<size_t>::max()); }

  /**
   * @returns an iterator, pointing at the last track
   */
  TrackIter track_rbegin() { return TrackIter(*this, _tracks.size()-1); }

  /**
   * @returns an iterator pointing at the first step of the give track
   * @param   track   Track object
   */
  StepIter step_begin(const Track& track) { return StepIter(*this, track, 0); }
  /**
   * @returns an iterator that signals the end-condition
   */
  StepIter step_end(const Track& track)   { return StepIter(*this, track, std::numeric_limits<size_t>::max()); }


  /**
   * read data from the hdf5 file.
   * @param parent_grp   id of the group containing the Event group. Event will be read from <parent_grp>/<name>.
   * @param name         name of the Event group.
   * @returns true if success, false otherwise
   */
  bool readData(hid_t parent_grp, const std::string& name);
  /**
   * write data to the hdf5 file.
   * @param parent_grp   id of the group containing the Event group. Event will be written to <parent_grp>/<name>.
   * @param name         name of the Event group.
   * @returns true if success, false otherwise
   */
  bool writeData(hid_t parent_grp, const std::string& name);

  /**
   * getter for event id
   */
  inline unsigned long ID() const { return _ID; }
  /**
   * setter for event id
   */
  inline void setID(unsigned long id) { _ID = id; }

  /**
   * getter for energy deposit
   */
  inline double EnergyDeposit() const { return _energyDeposit; }
  /**
   * setter for energy deposit
   */
  inline void setEnergyDeposit(double e) { _energyDeposit = e; }

  /**
   * convert region name to its id
   */
  inline int requestRegion(const std::string& name)
  {
    return _regions.find(name);
  }

  inline std::string region(int region_id) const
  { return _regions[region_id]; }

  inline std::string material(int material_id) const
  { return _materials[material_id]; }

  const LookupTable& regions() const {return _regions;}

  const std::vector<std::string>& materials() const {return _materials; }

private:
  Track::Dataset _tracks;
  Step::Dataset _steps;
  const LookupTable & _regions;
  const std::vector<std::string> & _materials;

  unsigned long _ID;
  double _energyDeposit;

};

}
}

#endif // HAVE_HDF5

#endif // PARTICLEEVENT_H
