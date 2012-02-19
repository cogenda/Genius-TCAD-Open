// $Id: parallel.h,v 1.7 2008/06/21 12:42:19 gdiso Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#ifndef __parallel_h__
#define __parallel_h__

// System includes
#include <string>
#include <vector>
#include <map>
#include <set>
#include <complex>

// Local includes
#include "genius_common.h" // for Real
#include "genius_env.h"
#include "vector_value.h"
#include "tensor_value.h"
#include "perf_log.h"


// Macro to identify and debug functions which should only be called in
// parallel on every processor at once

#undef parallel_only
#ifdef DEBUG
  #define parallel_only() { genius_assert(Parallel::verify(std::string(__FILE__))); assert(Parallel::verify(__LINE__)); }
#else
  #define parallel_only() { }
#endif

/**
 * The Parallel namespace is for wrapper functions
 * for common general parallel synchronization tasks.
 *
 * For MPI 1.1 compatibility, temporary buffers are used
 * instead of MPI 2's MPI_IN_PLACE
 */
namespace Parallel
{
#ifdef HAVE_MPI
  //-------------------------------------------------------------------
  /**
   * Data types for communication
   */
  typedef MPI_Datatype data_type;
  typedef MPI_Datatype data_type_with_int;

  /**
   * Templated function to return the appropriate MPI datatype
   * for use with built-in C types
   */
  template <typename T>
  inline data_type datatype();

  template <typename T>
  inline data_type datatype_with_int();

  /**
   * Request object for non-blocking I/O
   */
  typedef MPI_Request request;

  /**
   * Default message tag id
   */
  const int any_tag=MPI_ANY_TAG;

  /**
   * Accept from any source
   */
  const int any_source=MPI_ANY_SOURCE;

  /**
   * MPI error string
   */
  inline std::string error_string(int error_code)
  {
    char error[MPI_MAX_ERROR_STRING];
    int length;
    MPI_Error_string(error_code, error, &length);
    return std::string(error);
  }

#else
  // These shouldn't be needed
  typedef unsigned int data_type;
  typedef unsigned int request;

  const int any_tag=-1;
  const int any_source=0;
#endif // HAVE_MPI



  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Status struct.  Allows the source and size
   * of the message to be determined.
   */
  class MPIStatus
  {
  public:
    MPIStatus () {}

#ifndef HAVE_MPI

    int source () const
      { return 0; }

#else

    MPIStatus (const MPI_Status &mpi_status,
            const MPI_Datatype &data_type) :
        _status(mpi_status),
        _datatype(data_type)
    {}

    int source () const
    {
      return _status.MPI_SOURCE;
    }

    unsigned int size () const
    {
      int msg_size;
      MPI_Get_count (const_cast<MPI_Status*>(&_status), _datatype, &msg_size);
      assert (msg_size >= 0);
      return msg_size;
    }

  private:

    MPI_Status   _status;
    MPI_Datatype _datatype;
#endif

  };


  //-------------------------------------------------------------------
  /**
   * Pause execution until all processors reach a certain point.
   */
  inline void barrier ()
  {
#ifdef HAVE_MPI
    MPI_Barrier (Genius::comm_world());
#endif
    return;
  }

  //-------------------------------------------------------------------
  /**
   * Verify that a local variable has the same value on all processors
   */
  template <typename T>
  inline bool verify(const T &r);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the minimum of it's values
   * on all processors
   */
  template <typename T>
  inline void min(T &r);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the minimum
   * of it's values on all processors
   */
  template <typename T>
  inline void min(std::vector<T> &r);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the maximum of it's values
   * on all processors
   */
  template <typename T>
  inline void max(T &r);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the maximum
   * of it's values on all processors
   */
  template <typename T>
  inline void max(std::vector<T> &r);

  //-------------------------------------------------------------------
  /**
   * Take a local variable, replace it with the maximum of its values on
   * all processors, along with a flag indicating whether the global
   * maximum resides on the processor
   */
  template <typename T>
  inline void max_loc(T &r, unsigned int &loc);

  //-------------------------------------------------------------------
  /**
   * Take a local variable, replace it with the minimum of its values on
   * all processors, along with a flag indicating whether the global
   * minimum resides on the processor
   */
  template <typename T>
  inline void min_loc(T &r, unsigned int &loc);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the sum of it's values
   * on all processors
   */
  template <typename T>
  inline void sum(T &r);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the sum of
   * it's values on all processors
   */
  template <typename T>
  inline void sum(std::vector<T> &r);

  //-------------------------------------------------------------------
  /**
   * Blocking-send vector to one processor.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
                    std::vector<T> &buf,
                    const int tag=0);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-send vector to one processor.
   */
  template <typename T>
  inline void isend (const unsigned int dest_processor_id,
                     std::vector<T> &buf,
                     request &r,
                     const int tag=0);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-send vector to one processor with user-defined type.
   */
  template <typename T>
  inline void isend (const unsigned int dest_processor_id,
                     std::vector<T> &buf,
                     data_type &type,
                     request &r,
                     const int tag=0);

  //-------------------------------------------------------------------
  /**
   * Blocking-receive vector from one processor.
   */
  template <typename T>
  inline MPIStatus recv (const unsigned int src_processor_id,
                      std::vector<T> &buf,
                      const int tag=any_tag);

  //-------------------------------------------------------------------
  /**
   * Blocking-receive vector from one processor with user-defined type
   */
  template <typename T>
  inline MPIStatus recv (const unsigned int src_processor_id,
                      std::vector<T> &buf,
                      data_type &type,
                      const int tag=any_tag);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive vector from one processor.
   */
  template <typename T>
  inline void irecv (const unsigned int src_processor_id,
                     std::vector<T> &buf,
                     request &r,
                     const int tag=any_tag);

  //-------------------------------------------------------------------
  /**
   * Wait for a non-blocking send or receive to finish
   */
  inline void wait (request &r);

  //-------------------------------------------------------------------
  /**
   * Wait for a non-blocking send or receive to finish
   */
  inline void wait (std::vector<request> &r);

  //-------------------------------------------------------------------
  /**
   * Send vector send to one processor while simultaneously receiving
   * another vector recv from a (potentially different) processor.
   */
  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           T &send,
                           const unsigned int source_processor_id,
                           T &recv);

  //-------------------------------------------------------------------
  /**
   * Send vector send to one processor while simultaneously receiving
   * another vector recv from a (potentially different) processor using
   * a user-specified MPI Dataype.
   */
  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           T &send,
                           const unsigned int source_processor_id,
                           T &recv,
                           data_type &type);

  //-------------------------------------------------------------------
  /**
   * Take a vector of length n_processors, and on processor root_id fill in
   * recv[processor_id] = the value of send on processor processor_id
   */
  template <typename T>
  inline void gather(const unsigned int root_id,
                     T send,
                     std::vector<T> &recv);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and expand it on processor root_id
   * to include values from all processors
   */
  template <typename T>
  inline void gather(const unsigned int root_id, std::vector<T> &r);


  //-------------------------------------------------------------------
  /**
   * Take a distributed set
   */
  template <typename T>
  inline void gather(const unsigned int root_id, std::set<T> &r);

  //-------------------------------------------------------------------
  /**
   * gather a distributed map
   */
  template <typename T1, typename T2>
  inline void gather(const unsigned int root_id, std::map<T1, T2> &map);

  //-------------------------------------------------------------------
  /**
   * Take a vector of length n_processors, and fill in recv[processor_id] = the
   * value of send on that processor
   */
  template <typename T>
  inline void allgather(T send,
                        std::vector<T> &recv);


  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and expand it to include
   * values from all processors
   */
  template <typename T>
  inline void allgather(std::vector<T> &r);


  //-------------------------------------------------------------------
  /**
   * Take a set of local variables and expand it to include
   * values from all processors
   */
  template <typename T>
  inline void allgather(std::set<T> &r);

  //-------------------------------------------------------------------
  /**
   * allgather a distributed map
   */
  template <typename T1, typename T2>
  inline void allgather(std::map<T1, T2> &map);

  //-------------------------------------------------------------------
  /**
   * Effectively transposes the input vector across all processors.
   * The jth entry on processor i is replaced with the ith entry
   * from processor j.
   */
  template <typename T>
  inline void alltoall(std::vector<T> &r);



  //-------------------------------------------------------------------
  /**
   * Take a local value and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.
   */
  template <typename T>
  inline void broadcast(T &data, const unsigned int root_id=0);



  //-------------------------------------------------------------------
  /**
   * Take a local vector and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.
   */
  template <typename T>
  inline void broadcast(std::vector<T> &data, const unsigned int root_id=0);

  // gcc appears to need an additional declaration to make sure it
  // uses the right definition below
  template <typename T>
  inline void broadcast(std::vector<std::complex<T> > &data, const unsigned int root_id=0);

  template <typename T>
  inline void broadcast(std::set<T> &set, const unsigned int root_id=0);

  template <typename T1, typename T2>
  inline void broadcast(std::map<T1, T2> &map, const unsigned int root_id=0);

  //-----------------------------------------------------------------------
  // Parallel members

  // Internal helper function to create vector<something_useable> from
  // vector<bool> for compatibility with MPI bitwise operations
  template <typename T>
  inline void pack_vector_bool(const std::vector<bool> &in,
                               std::vector<T> &out)
  {
    unsigned int data_bits = 8*sizeof(T);
    unsigned int in_size = in.size();
    unsigned int out_size = in_size/data_bits + (in_size%data_bits?1:0);
    out.clear();
    out.resize(out_size);
    for (unsigned int i=0; i != in_size; ++i)
    {
      unsigned int index = i/data_bits;
      unsigned int offset = i%data_bits;
      out[index] += (in[i]?1:0) << offset;
    }
  }

  // Internal helper function to create vector<something_useable> from
  // vector<bool> for compatibility with MPI byte operations
  template <typename T>
  inline void unpack_vector_bool(const std::vector<T> &in,
                                 std::vector<bool> &out)
  {
    unsigned int data_bits = 8*sizeof(T);
    // We need the output vector to already be properly sized
    unsigned int in_size = in.size();
    unsigned int out_size = out.size();
    genius_assert(out_size/data_bits + (out_size%data_bits?1:0) == in_size);

    for (unsigned int i=0; i != out_size; ++i)
    {
      unsigned int index = i/data_bits;
      unsigned int offset = i%data_bits;
      out[i] = in[index] << (data_bits-1-offset) >> (data_bits-1);
    }
  }

#ifdef HAVE_MPI
  template<>
  inline MPI_Datatype datatype<char>() { return MPI_CHAR; }

  template<>
  inline MPI_Datatype datatype<unsigned char>() { return MPI_UNSIGNED_CHAR; }

  template<>
  inline MPI_Datatype datatype<short int>() { return MPI_SHORT; }

  template<>
  inline MPI_Datatype datatype<unsigned short int>() { return MPI_UNSIGNED_SHORT; }

  template<>
  inline MPI_Datatype datatype<int>() { return MPI_INT; }

  template<>
  inline MPI_Datatype datatype<unsigned int>() { return MPI_UNSIGNED; }

  template<>
  inline MPI_Datatype datatype<long>() { return MPI_LONG; }

  template<>
  inline MPI_Datatype datatype<long long>() { return MPI_LONG_LONG; }

  template<>
  inline MPI_Datatype datatype<unsigned long>() { return MPI_UNSIGNED_LONG; }

  template<>
  inline MPI_Datatype datatype<unsigned long long>() { return MPI_UNSIGNED_LONG_LONG; }

  template<>
  inline MPI_Datatype datatype<float>() { return MPI_FLOAT; }

  template<>
  inline MPI_Datatype datatype<double>() { return MPI_DOUBLE; }

  template<>
  inline MPI_Datatype datatype<long double>() { return MPI_LONG_DOUBLE; }

  template<>
  inline MPI_Datatype datatype_with_int<short int>() { return MPI_SHORT_INT; }

  template<>
  inline MPI_Datatype datatype_with_int<int>() { return MPI_2INT; }

  template<>
  inline MPI_Datatype datatype_with_int<long>() { return MPI_LONG_INT; }

  template<>
  inline MPI_Datatype datatype_with_int<long long>() { return MPI_LONG_LONG_INT; }

  template<>
  inline MPI_Datatype datatype_with_int<float>() { return MPI_FLOAT_INT; }

  template<>
  inline MPI_Datatype datatype_with_int<double>() { return MPI_DOUBLE_INT; }

  template<>
  inline MPI_Datatype datatype_with_int<long double>() { return MPI_LONG_DOUBLE_INT; }

  template <typename T>
  inline bool verify(const T &r)
  {
    if (Genius::n_processors() > 1)
    {
      T tempmin = r, tempmax = r;
      Parallel::min(tempmin);
      Parallel::max(tempmax);
      bool verified = (r == tempmin) &&
                      (r == tempmax);
      Parallel::min(verified);
      return verified;
    }
    return true;
  }


  template <>
  inline bool verify(const std::string & r)
  {
    if (Genius::n_processors() > 1)
    {
      // Cannot use <char> since MPI_MIN is not
      // strictly defined for chars!
      std::vector<short int> temp; temp.reserve(r.size());
      for (unsigned int i=0; i != r.size(); ++i)
        temp.push_back(r[i]);
      return Parallel::verify(temp);
    }
    return true;
  }


  template <typename T>
  inline void min(T &r)
  {
    START_LOG("min()", "Parallel");

    if (Genius::n_processors() > 1)
    {
      T temp;
      MPI_Allreduce (&r,
                     &temp,
                     1,
                     datatype<T>(),
                     MPI_MIN,
                     Genius::comm_world());
      r = temp;
    }

    STOP_LOG("min()", "Parallel");
  }


  template <>
  inline void min(bool &r)
  {
    START_LOG("min()", "Parallel");

    if (Genius::n_processors() > 1)
    {
      unsigned int tempsend = r;
      unsigned int temp;
      MPI_Allreduce (&tempsend,
                     &temp,
                     1,
                     datatype<unsigned int>(),
                     MPI_MIN,
                     Genius::comm_world());
      r = temp;
    }

    STOP_LOG("min()", "Parallel");
  }


  template <typename T>
  inline void min(std::vector<T> &r)
  {
    START_LOG("min()", "Parallel");

    if (Genius::n_processors() > 1)
    {
      std::vector<T> temp(r.size());
      MPI_Allreduce (r.size() ? &r[0] : NULL,
                     temp.size() ? &temp[0] : NULL,
                     r.size(),
                     datatype<T>(),
                     MPI_MIN,
                     Genius::comm_world());
      r = temp;
    }

    STOP_LOG("min()", "Parallel");
  }


  template <>
  inline void min(std::vector<bool> &r)
  {
    START_LOG("min()", "Parallel");

    if (Genius::n_processors() > 1)
    {
      std::vector<unsigned int> rchar;
      pack_vector_bool(r, rchar);
      std::vector<unsigned int> temp(rchar.size());
      MPI_Allreduce (rchar.size() ? &rchar[0] : NULL,
                     temp.size() ? &temp[0] : NULL,
                     rchar.size(),
                     datatype<unsigned int>(),
                     MPI_BAND,
                     Genius::comm_world());
      unpack_vector_bool(temp, r);
    }

    STOP_LOG("min()", "Parallel");
  }


  template <typename T>
  inline void max(T &r)
  {
    START_LOG("max()", "Parallel");

    if (Genius::n_processors() > 1)
    {
      T temp;
      MPI_Allreduce (&r,
                     &temp,
                     1,
                     datatype<T>(),
                     MPI_MAX,
                     Genius::comm_world());
      r = temp;
    }

    STOP_LOG("max()", "Parallel");
  }


  template <>
  inline void max(bool &r)
  {
    START_LOG("max()", "Parallel");

    if (Genius::n_processors() > 1)
    {
      unsigned int tempsend = r;
      unsigned int temp;
      MPI_Allreduce (&tempsend,
                     &temp,
                     1,
                     datatype<unsigned int>(),
                     MPI_MAX,
                     Genius::comm_world());
      r = temp;
    }

    STOP_LOG("max()", "Parallel");
  }


  template <typename T>
  inline void max(std::vector<T> &r)
  {
    START_LOG("max()", "Parallel");

    if (Genius::n_processors() > 1)
    {
      std::vector<T> temp(r.size());
      MPI_Allreduce (r.size() ? &r[0] : NULL,
                     temp.size() ? &temp[0] : NULL,
                     r.size(),
                     datatype<T>(),
                     MPI_MAX,
                     Genius::comm_world());
      r = temp;
    }

    STOP_LOG("max()", "Parallel");
  }


  template <>
  inline void max(std::vector<bool> &r)
  {
    START_LOG("max()", "Parallel");

    if (Genius::n_processors() > 1)
    {
      std::vector<unsigned int> rchar;
      pack_vector_bool(r, rchar);
      std::vector<unsigned int> temp(rchar.size());
      MPI_Allreduce (rchar.size() ? &rchar[0] : NULL,
                     temp.size() ? &temp[0] : NULL,
                     rchar.size(),
                     datatype<unsigned int>(),
                     MPI_BOR,
                     Genius::comm_world());
      unpack_vector_bool(temp, r);
    }

    STOP_LOG("max()", "Parallel");
  }


  template <typename T>
  inline void min_loc(T &r, unsigned int &loc)
  {
    if (Genius::n_processors() > 1)
    {
      struct
      {
        T val;
        unsigned int rank;
      } in, out;

      in.val  = r;
      in.rank = Genius::processor_id();

      MPI_Allreduce (&in,
                     &out,
                     1,
                     datatype_with_int<T>(),
                     MPI_MINLOC,
                     Genius::comm_world());
      r   = out.val;
      loc = out.rank;
    }
    else
      loc = Genius::processor_id();
  }


  template <typename T>
  inline void max_loc(T &r, unsigned int &loc)
  {
    if (Genius::n_processors() > 1)
    {
      struct
      {
        T val;
        unsigned int rank;
      } in, out;

      in.val  = r;
      in.rank = Genius::processor_id();

      MPI_Allreduce (&in,
                     &out,
                     1,
                     datatype_with_int<T>(),
                     MPI_MAXLOC,
                     Genius::comm_world());
      r   = out.val;
      loc = out.rank;
    }
    else
      loc = Genius::processor_id();
  }


  template <typename T>
  inline void sum(T &r)
  {
    if (Genius::n_processors() > 1)
    {
      T temp = r;
      MPI_Allreduce (&temp,
                     &r,
                     1,
                     datatype<T>(),
                     MPI_SUM,
                     Genius::comm_world());
    }
  }


  template <typename T>
  inline void sum(std::vector<T> &r)
  {
    if (Genius::n_processors() > 1 && !r.empty())
    {
      std::vector<T> temp(r);
      MPI_Allreduce (&temp[0],
                     &r[0],
                     r.size(),
                     datatype<T>(),
                     MPI_SUM,
                     Genius::comm_world());
    }
  }


  template <typename T>
  inline void sum(std::complex<T> &r)
  {
    if (Genius::n_processors() > 1)
    {
      std::complex<T> temp(r);
      MPI_Allreduce (&temp,
                     &r,
                     2,
                     datatype<T>(),
                     MPI_SUM,
                     Genius::comm_world());
    }
  }


  template <typename T>
  inline void sum(std::vector<std::complex<T> > &r)
  {
    if (Genius::n_processors() > 1 && !r.empty())
    {
      std::vector<std::complex<T> > temp(r);
      MPI_Allreduce (&temp[0],
                     &r[0],
                     r.size() * 2,
                     datatype<T>(),
                     MPI_SUM,
                     Genius::comm_world());
    }
  }



  template <typename T>
  inline void send (const unsigned int dest_processor_id,
                    std::vector<T> &buf,
                    const int tag)
  {
    START_LOG("send()", "Parallel");

    const int ierr =
      MPI_Send (buf.empty() ? NULL : &buf[0],
                buf.size(),
                datatype<T>(),
                dest_processor_id,
                tag,
                Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("send()", "Parallel");
  }


  template <typename T>
  inline void send (const unsigned int dest_processor_id,
                    std::vector<std::complex<T> > &buf,
                    const int tag)
  {
    START_LOG("send()", "Parallel");

    const int ierr =
      MPI_Send (buf.empty() ? NULL : &buf[0],
                buf.size() * 2,
                datatype<T>(),
                dest_processor_id,
                tag,
                Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("send()", "Parallel");
  }



  template <typename T>
  inline void isend (const unsigned int dest_processor_id,
                     std::vector<T> &buf,
                     request &r,
                     const int tag)
  {
    START_LOG("isend()", "Parallel");

    const int ierr =
      MPI_Isend (buf.empty() ? NULL : &buf[0],
                 buf.size(),
                 datatype<T>(),
                 dest_processor_id,
                 tag,
                 Genius::comm_world(),
                 &r);
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("isend()", "Parallel");
  }


  template <typename T>
  inline void isend (const unsigned int dest_processor_id,
                     std::vector<std::complex<T> > &buf,
                     request &r,
                     const int tag)
  {
    START_LOG("isend()", "Parallel");

    const int ierr =
      MPI_Isend (buf.empty() ? NULL : &buf[0],
                 buf.size() * 2,
                 datatype<T>(),
                 dest_processor_id,
                 tag,
                 Genius::comm_world(),
                 &r);
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("isend()", "Parallel");
  }



  template <typename T>
  inline void isend (const unsigned int dest_processor_id,
                     std::vector<T> &buf,
                     MPI_Datatype &type,
                     request &r,
                     const int tag)
  {
    START_LOG("isend()", "Parallel");

    const int ierr =
      MPI_Isend (buf.empty() ? NULL : &buf[0],
                 buf.size(),
                 type,
                 dest_processor_id,
                 tag,
                 Genius::comm_world(),
                 &r);
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("isend()", "Parallel");
  }



  template <typename T>
  inline MPIStatus recv (const unsigned int src_processor_id,
                      std::vector<T> &buf,
                      const int tag)
  {
    START_LOG("recv()", "Parallel");

    MPI_Status status;

    const int ierr =
      MPI_Recv (buf.empty() ? NULL : &buf[0],
                buf.size(),
                datatype<T>(),
                src_processor_id,
                tag,
                Genius::comm_world(),
                &status);
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("recv()", "Parallel");

    return MPIStatus(status, datatype<T>());
  }



  template <typename T>
  inline MPIStatus recv (const unsigned int src_processor_id,
                      std::vector<T> &buf,
                      MPI_Datatype &type,
                      const int tag)
  {
    START_LOG("recv()", "Parallel");

    MPI_Status status;

    const int ierr =
      MPI_Recv (buf.empty() ? NULL : &buf[0],
                buf.size(),
                type,
                src_processor_id,
                tag,
                Genius::comm_world(),
                &status);
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("recv()", "Parallel");

    return MPIStatus(status, type);
  }


  template <typename T>
  inline MPIStatus recv (const unsigned int src_processor_id,
                      std::vector<std::complex<T> > &buf,
                      const int tag)
  {
    START_LOG("recv()", "Parallel");

    MPI_Status status;

    const int ierr =
      MPI_Recv (buf.empty() ? NULL : &buf[0],
                buf.size() * 2,
                datatype<T>(),
                src_processor_id,
                tag,
                Genius::comm_world(),
                &status);
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("recv()", "Parallel");

    return MPIStatus(status, datatype<T>());
  }



  template <typename T>
  inline void irecv (const unsigned int src_processor_id,
                     std::vector<T> &buf,
                     request &r,
                     const int tag)
  {
    START_LOG("irecv()", "Parallel");

    const int ierr =
      MPI_Irecv (buf.empty() ? NULL : &buf[0],
                 buf.size(),
                 datatype<T>(),
                 src_processor_id,
                 tag,
                 Genius::comm_world(),
                 &r);
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("irecv()", "Parallel");
  }


  template <typename T>
  inline void irecv (const unsigned int src_processor_id,
                     std::vector<std::complex<T> > &buf,
                     request &r,
                     const int tag)
  {
    START_LOG("irecv()", "Parallel");

    const int ierr =
      MPI_Irecv (buf.empty() ? NULL : &buf[0],
                 buf.size() * 2,
                 datatype<T>(),
                 src_processor_id,
                 tag,
                 Genius::comm_world(),
                 &r);
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("irecv()", "Parallel");
  }



  inline void wait (request &r)
  {
    START_LOG("wait()", "Parallel");

    MPI_Wait (&r, MPI_STATUS_IGNORE);

    STOP_LOG("wait()", "Parallel");
  }



  inline void wait (std::vector<request> &r)
  {
    START_LOG("wait()", "Parallel");

    MPI_Waitall (r.size(), r.empty() ? NULL : &r[0], MPI_STATUSES_IGNORE);

    STOP_LOG("wait()", "Parallel");
  }



  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           T &send,
                           const unsigned int source_processor_id,
                           T &recv)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == Genius::processor_id() &&
        source_processor_id == Genius::processor_id())
    {
      recv = send;
      STOP_LOG("send_receive()", "Parallel");
      return;
    }

    MPI_Status status;
    MPI_Sendrecv(&send, 1, datatype<T>(),
                 dest_processor_id, 0,
                 &recv, 1, datatype<T>(),
                 source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    STOP_LOG("send_receive()", "Parallel");
  }


  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           std::complex<T> &send,
                           const unsigned int source_processor_id,
                           std::complex<T> &recv)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == Genius::processor_id() &&
        source_processor_id == Genius::processor_id())
    {
      recv = send;
      STOP_LOG("send_receive()", "Parallel");
      return;
    }

    MPI_Status status;
    MPI_Sendrecv(&send, 2, datatype<T>(),
                 dest_processor_id, 0,
                 &recv, 2, datatype<T>(),
                 source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    STOP_LOG("send_receive()", "Parallel");
  }



  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           std::vector<T> &send,
                           const unsigned int source_processor_id,
                           std::vector<T> &recv)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == Genius::processor_id() &&
        source_processor_id == Genius::processor_id())
    {
      recv = send;
      STOP_LOG("send_receive()", "Parallel");
      return;
    }

    // Trade buffer sizes first
    unsigned int sendsize = send.size(), recvsize;
    MPI_Status status;
    MPI_Sendrecv(&sendsize, 1, datatype<unsigned int>(),
                 dest_processor_id, 0,
                 &recvsize, 1, datatype<unsigned int>(),
                 source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    recv.resize(recvsize);

    MPI_Sendrecv(sendsize ? &send[0] : NULL, sendsize, datatype<T>(),
                 dest_processor_id, 0,
                 recvsize ? &recv[0] : NULL, recvsize, datatype<T>(),
                 source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    STOP_LOG("send_receive()", "Parallel");
  }


  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           std::vector<std::complex<T> > &send,
                           const unsigned int source_processor_id,
                           std::vector<std::complex<T> > &recv)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == Genius::processor_id() &&
        source_processor_id == Genius::processor_id())
    {
      recv = send;
      STOP_LOG("send_receive()", "Parallel");
      return;
    }

    // Trade buffer sizes first
    unsigned int sendsize = send.size(), recvsize;
    MPI_Status status;
    MPI_Sendrecv(&sendsize, 1, datatype<unsigned int>(),
                 dest_processor_id, 0,
                 &recvsize, 1, datatype<unsigned int>(),
                 source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    recv.resize(recvsize);

    MPI_Sendrecv(sendsize ? &send[0] : NULL, sendsize * 2, datatype<T>(),
                 dest_processor_id, 0,
                 recvsize ? &recv[0] : NULL, recvsize * 2, datatype<T>(),
                 source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    STOP_LOG("send_receive()", "Parallel");
  }



  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           std::vector<T> &send,
                           const unsigned int source_processor_id,
                           std::vector<T> &recv,
                           MPI_Datatype &type)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == Genius::processor_id() &&
        source_processor_id == Genius::processor_id())
    {
      recv = send;
      STOP_LOG("send_receive()", "Parallel");
      return;
    }

    // Trade buffer sizes first
    unsigned int sendsize = send.size(), recvsize;
    MPI_Status status;
    MPI_Sendrecv(&sendsize, 1, datatype<unsigned int>(),
                 dest_processor_id, 0,
                 &recvsize, 1, datatype<unsigned int>(),
                 source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    recv.resize(recvsize);

    MPI_Sendrecv(sendsize ? &send[0] : NULL, sendsize, type,
                 dest_processor_id, 0,
                 recvsize ? &recv[0] : NULL, recvsize, type,
                 source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    STOP_LOG("send_receive()", "Parallel");
  }



  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           std::vector<std::vector<T> > &send,
                           const unsigned int source_processor_id,
                           std::vector<std::vector<T> > &recv)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == Genius::processor_id() &&
        source_processor_id == Genius::processor_id())
    {
      recv = send;
      STOP_LOG("send_receive()", "Parallel");
      return;
    }

    // Trade outer buffer sizes first
    unsigned int sendsize = send.size(), recvsize;
    MPI_Status status;
    MPI_Sendrecv(&sendsize, 1, datatype<unsigned int>(),
                 dest_processor_id, 0,
                 &recvsize, 1, datatype<unsigned int>(),
                 source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    recv.resize(recvsize);

    // Trade inner buffer sizes next
    std::vector<unsigned int> sendsizes(sendsize), recvsizes(recvsize);
    unsigned int sendsizesum = 0, recvsizesum = 0;
    for (unsigned int i = 0; i != sendsize; ++i)
    {
      sendsizes[i] = send[i].size();
      sendsizesum += sendsizes[i];
    }

    MPI_Sendrecv(sendsize ? &sendsizes[0] : NULL, sendsize,
                 datatype<unsigned int>(), dest_processor_id, 0,
                 recvsize ? &recvsizes[0] : NULL, recvsize,
                 datatype<unsigned int>(), source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    for (unsigned int i = 0; i != recvsize; ++i)
    {
      recvsizesum += recvsizes[i];
      recv[i].resize(recvsizes[i]);
    }

    // Build temporary buffers third
    // We can't do multiple Sendrecv calls instead because send.size() may
    // differ on different processors
    std::vector<T> senddata(sendsizesum), recvdata(recvsizesum);

    // Fill the temporary send buffer
    typename std::vector<T>::iterator out = senddata.begin();
    for (unsigned int i = 0; i != sendsize; ++i)
    {
      out = std::copy(send[i].begin(), send[i].end(), out);
    }
    genius_assert(out == senddata.end());

    MPI_Sendrecv(sendsizesum ? &senddata[0] : NULL, sendsizesum,
                 datatype<T>(), dest_processor_id, 0,
                 recvsizesum ? &recvdata[0] : NULL, recvsizesum,
                 datatype<T>(), source_processor_id, 0,
                 Genius::comm_world(),
                 &status);

    // Empty the temporary recv buffer
    typename std::vector<T>::iterator in = recvdata.begin();
    for (unsigned int i = 0; i != recvsize; ++i)
    {
      std::copy(in, in + recvsizes[i], recv[i].begin());
      in += recvsizes[i];
    }
    genius_assert(in == recvdata.end());

    STOP_LOG("send_receive()", "Parallel");
  }


  template <typename T>
  inline void gather(const unsigned int root_id,
                     T send,
                     std::vector<T> &recv)
  {
    genius_assert(root_id < Genius::n_processors());

    if (Genius::processor_id() == root_id)
      recv.resize(Genius::n_processors());

    if (Genius::n_processors() > 1)
    {
      MPI_Gather(&send,
                 1,
                 datatype<T>(),
                 recv.empty() ? NULL : &recv[0],
                 1,
                 datatype<T>(),
                 root_id,
                 Genius::comm_world());
    }
    else
      recv[0] = send;
  }



  template <typename T>
  inline void gather(const unsigned int root_id,
                     std::complex<T> send,
                     std::vector<std::complex<T> > &recv)
  {
    genius_assert(root_id < Genius::n_processors());

    if (Genius::processor_id() == root_id)
      recv.resize(Genius::n_processors());

    if (Genius::n_processors() > 1)
    {
      MPI_Gather(&send,
                 2,
                 datatype<T>(),
                 recv.empty() ? NULL : &recv[0],
                 2,
                 datatype<T>(),
                 root_id,
                 Genius::comm_world());
    }
    else
      recv[0] = send;
  }



  /**
   * This function provides a convenient method
   * for combining vectors from each processor into one
   * contiguous chunk on one processor.  This handles the
   * case where the lengths of the vectors may vary.
   * Specifically, this function transforms this:
   \verbatim
    Processor 0: [ ... N_0 ]
    Processor 1: [ ....... N_1 ]
      ...
    Processor M: [ .. N_M]
   \endverbatim
   *
   * into this:
   *
   \verbatim
   [ [ ... N_0 ] [ ....... N_1 ] ... [ .. N_M] ]
   \endverbatim
   *
   * on processor root_id. This function is collective and therefore
   * must be called by all processors.
   */
  template <typename T>
  inline void gather(const unsigned int root_id,
                     std::vector<T> &r)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id()==root_id);
      return;
    }

    std::vector<int>
    sendlengths  (Genius::n_processors(), 0),
    displacements(Genius::n_processors(), 0);

    const int mysize = r.size();
    Parallel::allgather(mysize, sendlengths);

    START_LOG("gather()", "Parallel");

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0;
    for (unsigned int i=0; i != Genius::n_processors(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += sendlengths[i];
    }

    // Check for quick return
    if (globalsize == 0)
    {
      STOP_LOG("gather()", "Parallel");
      return;
    }

    // copy the input buffer
    std::vector<T> r_src = r;

    // now resize it to hold the global data
    // on the receiving processor
    if (root_id == Genius::processor_id())
      r.resize(globalsize);

    // and get the data from the remote processors
    const int ierr =
      MPI_Gatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
                   r.empty() ? NULL :  &r[0], &sendlengths[0],
                   &displacements[0], datatype<T>(),
                   root_id,
                   Genius::comm_world());

    assert (ierr == MPI_SUCCESS);

    STOP_LOG("gather()", "Parallel");
  }


  template <typename T>
  inline void gather(const unsigned int root_id,
                     std::vector<std::complex<T> > &r)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id()==root_id);
      return;
    }

    std::vector<int>
    sendlengths  (Genius::n_processors(), 0),
    displacements(Genius::n_processors(), 0);

    const int mysize = r.size() * 2;
    Parallel::allgather(mysize, sendlengths);

    START_LOG("gather()", "Parallel");

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0;
    for (unsigned int i=0; i != Genius::n_processors(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += sendlengths[i];
    }

    // Check for quick return
    if (globalsize == 0)
    {
      STOP_LOG("gather()", "Parallel");
      return;
    }

    // Make temporary buffers for the input/output data
    std::vector<std::complex<T> > r_src(r);

    // now resize r to hold the global data
    // on the receiving processor
    if (root_id == Genius::processor_id())
      r.resize(globalsize/2);

    // and get the data from the remote processors
    const int ierr =
      MPI_Gatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
                   r.empty() ? NULL : &r[0], &sendlengths[0],
                   &displacements[0], datatype<T>(),
                   root_id, Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("gather()", "Parallel");
  }


  /**
   * gather a distributed set
   */
  template <typename T>
  inline void gather(const unsigned int root_id,
                     std::set<T> &set)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id()==root_id);
      return;
    }

    START_LOG("gather()", "Parallel");

    std::vector<T> key;
    key.reserve(set.size());

    typename std::set<T>::iterator set_it = set.begin();
    typename std::set<T>::iterator set_it_end = set.end();
    for(; set_it!=set_it_end; ++set_it)
    {
      key.push_back(*set_it);
    }

    gather(root_id, key);


    if(Genius::processor_id()==root_id)
    {
      for(unsigned int n=0; n<key.size(); n++)
        set.insert(key[n]);
    }

    STOP_LOG("gather()", "Parallel");
  }

  /**
   * gather a distributed map
   */
  template <typename T1, typename T2>
  inline void gather(const unsigned int root_id,
                     std::map<T1, T2> &map)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id()==root_id);
      return;
    }

    START_LOG("gather()", "Parallel");

    std::vector<T1> key;
    std::vector<T2> value;
    key.reserve(map.size());
    value.reserve(map.size());

    typename std::map<T1, T2>::iterator map_it = map.begin();
    typename std::map<T1, T2>::iterator map_it_end = map.end();
    for(; map_it!=map_it_end; ++map_it)
    {
      key.push_back((*map_it).first);
      value.push_back((*map_it).second);
    }

    gather(root_id, key);
    gather(root_id, value);


    if(Genius::processor_id()==root_id)
    {
      genius_assert(key.size()==value.size());
      for(unsigned int n=0; n<key.size(); n++)
        map.insert(std::pair<T1, T2>(key[n], value[n]) );
    }

    STOP_LOG("gather()", "Parallel");
  }



  template <typename T>
  inline void allgather(T send,
                        std::vector<T> &recv)
  {
    START_LOG ("allgather()","Parallel");
    //genius_assert(recv.size() == Genius::n_processors());
    recv.resize(Genius::n_processors());

    if (Genius::n_processors() > 1)
    {
      MPI_Allgather (&send,
                     1,
                     datatype<T>(),
                     &recv[0],
                     1,
                     datatype<T>(),
                     Genius::comm_world());
    }
    else
      recv[0] = send;

    STOP_LOG ("allgather()","Parallel");
  }



  template <typename T>
  inline void allgather(std::complex<T> send,
                        std::vector<std::complex<T> > &recv)
  {
    START_LOG ("allgather()","Parallel");
    //genius_assert(recv.size() == Genius::n_processors());
    recv.resize(Genius::n_processors());

    if (Genius::n_processors() > 1)
    {
      MPI_Allgather (&send,
                     2,
                     datatype<T>(),
                     &recv[0],
                     2,
                     datatype<T>(),
                     Genius::comm_world());
    }
    else
      recv[0] = send;

    STOP_LOG ("allgather()","Parallel");
  }



  /**
   * This function provides a convenient method
   * for combining vectors from each processor into one
   * contiguous chunk.  This handles the case where the
   * lengths of the vectors may vary.  Specifically, this
   * function transforms this:
   \verbatim
    Processor 0: [ ... N_0 ]
    Processor 1: [ ....... N_1 ]
      ...
    Processor M: [ .. N_M]
   \endverbatim
   *
   * into this:
   *
   \verbatim
   [ [ ... N_0 ] [ ....... N_1 ] ... [ .. N_M] ]
   \endverbatim
   *
   * on each processor. This function is collective and therefore
   * must be called by all processors.
   */
  template <typename T>
  inline void allgather(std::vector<T> &r)
  {
    if (Genius::n_processors() == 1)
      return;

    START_LOG("allgather()", "Parallel");

    std::vector<int>
    sendlengths  (Genius::n_processors(), 0),
    displacements(Genius::n_processors(), 0);

    const int mysize = r.size();
    Parallel::allgather(mysize, sendlengths);

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0;
    for (unsigned int i=0; i != Genius::n_processors(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += sendlengths[i];
    }

    // Check for quick return
    if (globalsize == 0)
    {
      STOP_LOG("allgather()", "Parallel");
      return;
    }

    // copy the input buffer
    std::vector<T> r_src(r);

    // now resize it to hold the global data
    r.resize(globalsize);

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
    const int ierr =
      MPI_Allgatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
                      r.empty()     ? NULL : &r[0],     &sendlengths[0],
                      &displacements[0], datatype<T>(), Genius::comm_world());

    assert (ierr == MPI_SUCCESS);

    STOP_LOG("allgather()", "Parallel");
  }


  template <typename T>
  inline void allgather(std::vector<std::complex<T> > &r)
  {
    if (Genius::n_processors() == 1)
      return;

    START_LOG("allgather()", "Parallel");

    std::vector<int>
    sendlengths  (Genius::n_processors(), 0),
    displacements(Genius::n_processors(), 0);

    const int mysize = r.size() * 2;
    Parallel::allgather(mysize, sendlengths);

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0;
    for (unsigned int i=0; i != Genius::n_processors(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += sendlengths[i];
    }

    // Check for quick return
    if (globalsize == 0)
    {
      STOP_LOG("allgather()", "Parallel");
      return;
    }

    // copy the input buffer
    std::vector<std::complex<T> > r_src(r);

    // now resize it to hold the global data
    r.resize(globalsize);

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
    const int ierr =
      MPI_Allgatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
                      r.empty()     ? NULL : &r[0],     &sendlengths[0],
                      &displacements[0], datatype<T>(),
                      Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("allgather()", "Parallel");
  }

  template <typename T>
  inline void allgather(std::vector< VectorValue<T> > &r)
  {
    if (Genius::n_processors() == 1)
      return;

    START_LOG("allgather()", "Parallel");

    std::vector<int>
        sendlengths  (Genius::n_processors(), 0),
    displacements(Genius::n_processors(), 0);

    const int mysize = r.size() * 3;
    Parallel::allgather(mysize, sendlengths);

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0;
    for (unsigned int i=0; i != Genius::n_processors(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += sendlengths[i];
    }

    // Check for quick return
    if (globalsize == 0)
    {
      STOP_LOG("allgather()", "Parallel");
      return;
    }

    // set the input/output buffer
    std::vector< T > r_src(mysize), r_dst(globalsize);
    for(unsigned int n=0; n<r.size(); ++n)
    {
      r_src[3*n+0] = r[n][0];
      r_src[3*n+1] = r[n][1];
      r_src[3*n+2] = r[n][2];
    }

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
    const int ierr =
        MPI_Allgatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
                        r.empty()     ? NULL : &r_dst[0], &sendlengths[0],
                        &displacements[0], datatype<T>(),
                        Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    r.clear();
    for(unsigned int n=0; n<r_dst.size(); )
    {
      r.push_back( VectorValue<T>(r_dst[n++], r_dst[n++], r_dst[n++]) );
    }

    STOP_LOG("allgather()", "Parallel");
  }

  template <typename T>
  inline void allgather(std::vector< TensorValue<T> > &r)
  {
    if (Genius::n_processors() == 1)
      return;

    START_LOG("allgather()", "Parallel");

    std::vector<int>
        sendlengths  (Genius::n_processors(), 0),
    displacements(Genius::n_processors(), 0);

    const int mysize = r.size() * 9;
    Parallel::allgather(mysize, sendlengths);

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0;
    for (unsigned int i=0; i != Genius::n_processors(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += sendlengths[i];
    }

    // Check for quick return
    if (globalsize == 0)
    {
      STOP_LOG("allgather()", "Parallel");
      return;
    }

    // set the input/output buffer
    std::vector< T > r_src(mysize), r_dst(globalsize);
    for(unsigned int n=0; n<r.size(); ++n)
    {
      r_src[9*n+0] = r[n][0];  r_src[9*n+1] = r[n][1];  r_src[9*n+2] = r[n][2];
      r_src[9*n+3] = r[n][3];  r_src[9*n+4] = r[n][4];  r_src[9*n+5] = r[n][5];
      r_src[9*n+6] = r[n][6];  r_src[9*n+7] = r[n][7];  r_src[9*n+8] = r[n][8];
    }

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
    const int ierr =
        MPI_Allgatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
                        r.empty()     ? NULL : &r_dst[0], &sendlengths[0],
                        &displacements[0], datatype<T>(),
                        Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    r.clear();
    for(unsigned int n=0; n<r_dst.size(); )
    {
      r.push_back( TensorValue<T>(r_dst[n++], r_dst[n++], r_dst[n++],
                                  r_dst[n++], r_dst[n++], r_dst[n++],
                                  r_dst[n++], r_dst[n++], r_dst[n++]
                                 ) );
    }

    STOP_LOG("allgather()", "Parallel");
  }



  /**
   * allgather a distributed set
   */
  template <typename T>
  inline void allgather(std::set<T> &set)
  {
    if (Genius::n_processors() == 1)
      return;

    START_LOG("allgather()", "Parallel");

    std::vector<T> key;
    key.reserve(set.size());

    typename  std::set<T>::iterator set_it = set.begin();
    typename  std::set<T>::iterator set_it_end = set.end();
    for(; set_it!=set_it_end; ++set_it)
    {
      key.push_back(*set_it);
    }

    allgather(key);

    for(unsigned int n=0; n<key.size(); ++n)
      set.insert(key[n]);

    STOP_LOG("allgather()", "Parallel");
  }


  /**
   * allgather a distributed map
   */
  template <typename T1, typename T2>
  inline void allgather(std::map<T1, T2> &map)
  {
    if (Genius::n_processors() == 1)
      return;

    START_LOG("allgather()", "Parallel");

    std::vector<T1> key;
    std::vector<T2> value;
    key.reserve(map.size());
    value.reserve(map.size());

    typename  std::map<T1, T2>::iterator map_it = map.begin();
    typename  std::map<T1, T2>::iterator map_it_end = map.end();
    for(; map_it!=map_it_end; ++map_it)
    {
      key.push_back((*map_it).first);
      value.push_back((*map_it).second);
    }

    allgather(key);
    allgather(value);
    genius_assert(key.size()==value.size());

    for(unsigned int n=0; n<key.size(); ++n)
      map.insert(std::pair<T1, T2>(key[n], value[n]) );

    STOP_LOG("allgather()", "Parallel");
  }


  /**
   * Replaces the input buffer with the result of MPI_Alltoall.
   * The vector size must be of te form N*n_procs, where N is
   * the number of elements to be sent/received from each
   * processor.
   */
  template <typename T>
  inline void alltoall(std::vector<T> &buf)
  {
    if (Genius::n_processors() == 1)
      return;

    START_LOG("alltoall()", "Parallel");

    // the per-processor size.  this is the same for all
    // processors using MPI_Alltoall, could be variable
    // using MPI_Alltoallv
    const unsigned int size_per_proc =
      buf.size()/Genius::n_processors();

    assert (buf.size() % Genius::n_processors() == 0);

    std::vector<T> tmp(buf);

    const int ierr =
      MPI_Alltoall (tmp.empty() ? NULL : &tmp[0],
                    size_per_proc,
                    datatype<T>(),
                    buf.empty() ? NULL : &buf[0],
                    size_per_proc,
                    datatype<T>(),
                    Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("alltoall()", "Parallel");
  }



  template <typename T>
  inline void broadcast (T &data, const unsigned int root_id)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id() == root_id);
      return;
    }

    START_LOG("broadcast()", "Parallel");

    // Spread data to remote processors.
    const int ierr = MPI_Bcast (&data, 1, datatype<T>(), root_id, Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }


  template <typename T>
  inline void broadcast (std::complex<T> &data, const unsigned int root_id)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id() == root_id);
      return;
    }

    START_LOG("broadcast()", "Parallel");

    // Spread data to remote processors.
    const int ierr = MPI_Bcast (&data, 2, datatype<T>(), root_id, Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }


  template <>
  inline void broadcast (std::string &data, const unsigned int root_id)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id() == root_id);
      return;
    }

    START_LOG("broadcast()", "Parallel");

    unsigned int data_size = data.size();
    Parallel::broadcast(data_size, root_id);

    // empty string?
    if(data_size==0)
    {
      data.clear();
      STOP_LOG("broadcast()", "Parallel");
      return;
    }

    std::vector<char> data_c(data_size);
    std::string orig(data);

    if (Genius::processor_id() == root_id)
      for(unsigned int i=0; i<data.size(); i++)
        data_c[i] = data[i];

    Parallel::broadcast (data_c,root_id);

    data.clear(); data.reserve(data_c.size());
    for(unsigned int i=0; i<data_c.size(); i++)
      data.push_back(data_c[i]);

    if (Genius::processor_id() == root_id)
      genius_assert(data == orig);

    STOP_LOG("broadcast()", "Parallel");
  }



  template <typename T>
  inline void broadcast (std::vector<T> &data, const unsigned int root_id)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id() == root_id);
      return;
    }

    START_LOG("broadcast()", "Parallel");

    unsigned int data_size = data.size();
    Parallel::broadcast(data_size, root_id);

    // empty
    if(data_size==0)
    {
      data.clear();
      STOP_LOG("broadcast()", "Parallel");
      return;
    }
    data.resize(data_size);


    // and get the data from the remote processors.
    const int ierr = MPI_Bcast (&data[0], data_size, datatype<T>(),root_id, Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }

  template <>
  inline void broadcast (std::vector<std::string> &data, const unsigned int root_id)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id() == root_id);
      return;
    }

    START_LOG("broadcast()", "Parallel");

    unsigned int data_size = data.size();
    Parallel::broadcast(data_size, root_id);

    // empty
    if(data_size==0)
    {
      data.clear();
      STOP_LOG("broadcast()", "Parallel");
      return;
    }
    data.resize(data_size);

    // pack to char array
    std::vector<char> data_c;
    std::vector<unsigned int> data_length;
    if (Genius::processor_id() == root_id)
      for(unsigned int i=0; i<data.size(); i++)
      {
        const std::string string = data[i];
        data_length.push_back(string.size());
        for(unsigned int n=0; n<string.size(); n++)
          data_c.push_back(string[n]);
      }

    Parallel::broadcast(data_c, root_id);
    Parallel::broadcast(data_length, root_id);

    // unpack
    if (Genius::processor_id() != root_id)
    {
      unsigned int data_c_index = 0;
      for(unsigned int i=0; i<data.size(); i++)
      {
        data[i].clear();
        data[i].reserve(data_length[i]);
        for(unsigned int n=0; n<data_length[i]; n++)
          data[i].push_back(data_c[data_c_index++]);
      }
    }

    STOP_LOG("broadcast()", "Parallel");
  }


  template <>
  inline void broadcast (std::vector<bool> &data, const unsigned int root_id)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id() == root_id);
      return;
    }

    START_LOG("broadcast()", "Parallel");

    unsigned int data_size = data.size();
    Parallel::broadcast(data_size, root_id);

    // empty
    if(data_size==0)
    {
      data.clear();
      STOP_LOG("broadcast()", "Parallel");
      return;
    }

    data.resize(data_size);

    std::vector<unsigned int> int_data;
    pack_vector_bool(data, int_data);

    // and get the data from the remote processors.
    const int ierr =  MPI_Bcast (&int_data[0], int_data.size(), datatype<unsigned int>(), root_id, Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    unpack_vector_bool(int_data, data);

    STOP_LOG("broadcast()", "Parallel");
  }


  template <typename T>
  inline void broadcast (std::vector<std::complex<T> > &data,
                         const unsigned int root_id)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id() == root_id);
      return;
    }

    START_LOG("broadcast()", "Parallel");

    unsigned int data_size = data.size();
    Parallel::broadcast(data_size, root_id);

    // empty
    if(data_size==0)
    {
      data.clear();
      STOP_LOG("broadcast()", "Parallel");
      return;
    }

    data.resize(data_size);


    // and get the data from the remote processors.
    const int ierr =  MPI_Bcast (&data[0], data_size * 2, datatype<T>(), root_id, Genius::comm_world());
    assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }


  template <typename T>
  inline void broadcast (std::set<T> &set, const unsigned int root_id)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id() == root_id);
      return;
    }

    START_LOG("broadcast()", "Parallel");

    std::vector<T> key;

    if(Genius::processor_id() == root_id)
    {
      typename  std::set<T>::iterator set_it = set.begin();
      for(; set_it!=set.end(); ++set_it)
      {
        key.push_back(*set_it);
      }
    }

    unsigned int data_size = key.size();
    Parallel::broadcast(data_size, root_id);

    // empty
    if(data_size==0)
    {
      set.clear();
      STOP_LOG("broadcast()", "Parallel");
      return;
    }

    key.resize(data_size);

    genius_assert( MPI_Bcast (&key[0], data_size, datatype<T>(), root_id, Genius::comm_world())== MPI_SUCCESS );

    if(Genius::processor_id() != root_id)
    {
      for(unsigned int n=0; n<key.size(); ++n)
        set.insert(key[n]);
    }

    STOP_LOG("broadcast()", "Parallel");
  }


  template <typename T1, typename T2>
  inline void broadcast (std::map<T1, T2> &map, const unsigned int root_id)
  {
    if (Genius::n_processors() == 1)
    {
      assert (Genius::processor_id() == root_id);
      return;
    }

    START_LOG("broadcast()", "Parallel");

    std::vector<T1> key;
    std::vector<T2> value;

    if(Genius::processor_id() == root_id)
    {
      typename  std::map<T1, T2>::iterator map_it = map.begin();
      for(; map_it!=map.end(); ++map_it)
      {
        key.push_back((*map_it).first);
        value.push_back((*map_it).second);
      }
    }
    broadcast(key, root_id);
    broadcast(value, root_id);

    genius_assert(key.size()==value.size());

    if(Genius::processor_id() != root_id)
    {
      for(unsigned int n=0; n<key.size(); ++n)
        map[key[n]] = value[n];
    }

    STOP_LOG("broadcast()", "Parallel");
  }




#else // HAVE_MPI

  template <typename T>
  inline bool verify(const T &) { return true; }

  template <typename T>
  inline void min(T &) {}

  template <typename T>
  inline void min(std::vector<T> &) {}

  template <typename T>
  inline void max(T &) {}

  template <typename T>
  inline void max(std::vector<T> &) {}

  template <typename T>
  inline void min_loc(T &, unsigned int &loc) { loc=0; }

  template <typename T>
  inline void max_loc(T &, unsigned int &loc) { loc=0; }

  template <typename T>
  inline void sum(T &) {}

  template <typename T>
  inline void sum(std::vector<T> &) {}

  // Blocking sends don't make sense on one processor
  template <typename T>
  inline void send (const unsigned int,
                    std::vector<T> &,
                    const unsigned int) { genius_error(); }

  template <typename T>
  inline void isend (const unsigned int,
                     std::vector<T> &,
                     request &,
                     const int) {}

  // Blocking receives don't make sense on one processor
  template <typename T>
  inline MPIStatus recv (const unsigned int,
                      std::vector<T> &,
                      const int) { genius_error(); return MPIStatus(); }

  template <typename T>
  inline void irecv (const unsigned int,
                     std::vector<T> &,
                     request &,
                     const int) {}

  inline void wait (request &) {}

  inline void wait (std::vector<request> &) {}

  template <typename T>
  inline void send_receive(const unsigned int send_tgt,
                           T &send,
                           const unsigned int recv_source,
                           T &recv)
  {
    assert (send_tgt == recv_source);
    recv = send;
  }

  template <typename T>
  inline void gather(const unsigned int root_id,
                     T send,
                     std::vector<T> &recv)
  {
    assert (!root_id);
    recv.resize(1);
    recv[0] = send;
  }

  template <typename T>
  inline void gather(const unsigned int, std::vector<T> &) {}


  template <typename T>
  inline void gather(const unsigned int, std::set<T> &) {}


  template <typename T1, typename T2>
  inline void gather(const unsigned int, std::map<T1, T2> &) {}


  template <typename T>
  inline void allgather(T send,
                        std::vector<T> &recv)
  {
    recv.resize(1);
    recv[0] = send;
  }

  template <typename T>
  inline void allgather(std::vector<T> &) {}

  template <typename T>
  inline void allgather(std::set<T> &) {}

  template <typename T1, typename T2>
  inline void allgather(std::map<T1, T2> &) {}

  template <typename T>
  inline void alltoall(std::vector<T> &) {}

  template <typename T>
  inline void broadcast (T &, const unsigned int) {}

  template <typename T>
  inline void broadcast (std::vector<T> &, const unsigned int) {}

  template <typename T>
  inline void broadcast (std::set<T> &, const unsigned int) {}

  template <typename T1, typename T2>
  inline void broadcast(std::map<T1, T2> &, const unsigned int) {}

#endif // HAVE_MPI

}

#endif // #define __parallel_h__
