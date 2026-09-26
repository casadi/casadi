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

      Pass it as first argument: f(S, x), f.call(S, arg). It is cleared when the call starts.
      Records are kept as a CBOR stream, the same format generated code writes with codegen
      option "stats".

      Options:
      - "size": buffer size in bytes (default 131072). Records that do not fit are dropped;
        see truncated() and nbytes(). */
  class CASADI_EXPORT StatsRecorder
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<StatsRecorder>) {
  public:
    /// Default constructor
    StatsRecorder();

    /// Constructor with options
    explicit StatsRecorder(const Dict& opts);

    /// Load a stream written by save or by generated code
    static StatsRecorder load(const std::string& filename);

    /** \brief Wrap a raw stream, e.g. the buffer of generated code

        In Python, pass a bytes object. */
    static StatsRecorder from_bytes(const std::string& bytes);

    /// Readable name of the public class
    static std::string type_name() {return "StatsRecorder";}

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    /// Drop all records (a call does this itself)
    void clear();

    /// Bytes needed, including those dropped
    casadi_int nbytes() const;

    /// Were records dropped because the buffer was full?
    bool truncated() const;

    /** \brief Write the stream to a file, for load

        For JSON, use export_json. */
    void save(const std::string& filename) const;

    /** \brief Reorder an interleaved (multithreaded) stream

        Same call tree, with each node's records contiguous, depth first. */
    StatsRecorder deinterleave() const;

    /// The call tree as JSON: {"version": 1, "calls": to_native()}
    std::string to_json() const;

    /// Write to_json to a file (not loadable)
    void export_json(const std::string& filename) const;

    /** \brief The call tree as native types

        A list of root calls, each a dict with name, id, mem, flag and stats. Nested calls go
        under children, or under sections (e.g. pre, post) and iterations when present. */
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
