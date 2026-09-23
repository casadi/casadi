/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_STATS_RECORDER_HPP
#define CASADI_STATS_RECORDER_HPP

#include "shared_object.hpp"
#include "printable.hpp"
#include "generic_type.hpp"

namespace casadi {
  // Forward declaration
  class StatsRecorderInternal;

  /** \brief Records the statistics of an evaluation as a call tree

      Pass it first: f(S, x), f.call(S, arg). It is reset when the call starts; afterwards
      it holds one node per call made (solvers, their oracle calls, the calls of a map, each
      with its iterate index, ...) as a CBOR stream, the same format generated code writes
      with codegen option "stats".

      Options:
      - "size": the fixed buffer, in bytes (default 131072). Records that do not fit are
        dropped and truncated() becomes true; nbytes() then says what size would have done. */
  class CASADI_EXPORT StatsRecorder
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<StatsRecorder>) {
  public:
    /// An empty recording of the default size
    StatsRecorder();

    /// An empty recording
    explicit StatsRecorder(const Dict& opts);

    /// Load a stream written by save, or elsewhere, e.g. the buffer of generated code
    static StatsRecorder load(const std::string& filename);

    /** \brief A stream written elsewhere, e.g. the buffer of generated code

        In Python, pass a bytes object. */
    static StatsRecorder from_bytes(const std::string& bytes);

    /// Readable name of the public class
    static std::string type_name() {return "StatsRecorder";}

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    /// Bytes the recording needed (more than the buffer holds if truncated)
    casadi_int nbytes() const;

    /// Were records dropped because the buffer was full?
    bool truncated() const;

    /** \brief Write the stream to a file, for load

        Not to a .json file: the call tree as JSON is export_json. */
    void save(const std::string& filename) const;

    /** \brief The records of an interleaved (multithreaded) stream, reordered

        Each call's records become contiguous, depth first, map iterates in order.
        The result describes the same call tree. */
    StatsRecorder deinterleave() const;

    /// The call tree as JSON: {"version": 1, "calls": to_native()}
    std::string to_json() const;

    /// Write the call tree as JSON (to_json) to a file; one-way, load takes the stream only
    void export_json(const std::string& filename) const;

    /** \brief The call tree as native types

        A list of root calls, each a dict with name, id, mem, flag and stats, and its calls
        under children, directly or by section: pre, iterations and post for a call with
        iteration rows (a solver), iterations (one per iterate) for a map. */
    GenericType to_native() const;

    /** \brief A stat of the last call of the function named fname

        key is a stats entry (e.g. "return_status"), or "iterations" for the
        per-iteration columns. */
    GenericType get_stat(const std::string& fname, const std::string& key) const;

#ifndef SWIG
    /** \brief  Create from node */
    static StatsRecorder create(StatsRecorderInternal* node);

    StatsRecorderInternal* get() const;

    /// Access a member function or object
    const StatsRecorderInternal* operator->() const;

    /// Access a member function or object
    StatsRecorderInternal* operator->();

  private:
    explicit StatsRecorder(StatsRecorderInternal* node);
#endif // SWIG
  };

} // namespace casadi

#endif // CASADI_STATS_RECORDER_HPP
