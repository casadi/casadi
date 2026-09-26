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

    /** \brief Store where each node ends, after the last write: queries skip subtrees

        Optional; queries give the same results without. */
    void index();

    /** \brief Next call after cursor (-1: first) at any depth below scope (-1: root)

        Nodes are stream offsets, the same as in generated code; -1 if none.
        A pattern containing ':' matches the full id (e.g. "#5:solver"), otherwise the
        function name; empty matches any. */
    casadi_int find_function(const std::string& pattern, casadi_int scope=-1,
      casadi_int cursor=-1) const;

    /// As find_function, among the direct children of parent
    casadi_int select_function(casadi_int parent, const std::string& pattern,
      casadi_int cursor=-1) const;

    /// Next iteration of parent after cursor (-1: first) with index (-1: any), or -1
    casadi_int select_iteration(casadi_int parent, casadi_int index,
      casadi_int cursor=-1) const;

    /// Iteration of parent recorded last, or -1
    casadi_int select_last_iteration(casadi_int parent) const;

    /// Next section of parent after cursor (-1: first) named name (empty: any), or -1
    casadi_int select_section(casadi_int parent, const std::string& name,
      casadi_int cursor=-1) const;

    /** \brief Value of key of node

        key is a stats entry (e.g. "return_status"), an iteration field (e.g. "obj"),
        or reserved: id, name, mem and flag of a call, name of a section, index of an
        iteration. */
    GenericType get_stat(casadi_int node, const std::string& key) const;

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
