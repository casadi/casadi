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


#ifndef CASADI_STATS_RECORDER_INTERNAL_HPP
#define CASADI_STATS_RECORDER_INTERNAL_HPP

#include "stats_recorder.hpp"
#include "shared_object.hpp"
#include "runtime/casadi_runtime.hpp"

#include <map>
#include <ostream>
#include <string>
#include <vector>

/// \cond INTERNAL
namespace casadi {

  /** \brief A node of a decoded stream: a call, a section or an iteration */
  struct StatsNode {
    enum Kind {CALL, SECTION, ITERATION};
    Kind kind = CALL;
    // name: function name for a call, section name for a section
    std::string id, name;
    GenericType mem;
    casadi_int index = -1, flag = 0;
    bool ended = false;
    casadi_int parent = -1;
    // In stream order: own records (index >= 0) and child nodes (-1-index)
    std::vector<casadi_int> events;
    casadi_int begin = -1, end = -1;
    Dict stats;
    // Declared iteration fields (call) or field values (iteration)
    bool has_fields = false;
    std::vector<std::string> fields;
    std::map<casadi_int, GenericType> field_values;
  };

  /** \brief A decoded record, split off its reference for rewriting */
  struct StatsRecord {
    casadi_int n, kind, ref;  // as casadi_stats_read_head reads them
    std::string tail;  // the other n-2 items
  };

  /** \brief A decoded stream */
  struct StatsTree {
    std::vector<StatsNode> nodes;
    std::vector<StatsRecord> records;
    std::vector<casadi_int> roots;
    bool truncated = false;
  };

  /** \brief Internal class for StatsRecorder */
  class CASADI_EXPORT StatsRecorderInternal : public SharedObjectInternal {
  public:
    /// Constructor
    explicit StatsRecorderInternal(casadi_int size);

    std::string class_name() const override {return "StatsRecorderInternal";}
    void disp(std::ostream& stream, bool more) const override;

    /// Drop all records
    void clear();

    /// Bytes needed, including those dropped
    casadi_int nbytes() const;

    /// The recorded bytes
    std::string bytes() const;

    /// Replace the stream
    void set_bytes(const std::string& b);

    /// Decode the recording
    StatsTree decode() const;

    /// Iteration fields of a call as columns, plus "iter"
    static Dict iteration_columns(const StatsTree& t, const StatsNode& c);

    /// Append a decoded record with a new reference (deinterleave)
    casadi_int put_record(const StatsRecord& r, casadi_int ref);

    /// Were records dropped?
    bool truncated() const;

    /// Size of the buffer
    casadi_int size() const { return size_; }

    /// The sink evaluations write to
    casadi_stats_sink* sink() { return &s_; }
    const casadi_stats_sink* sink() const { return &s_; }

  private:
    casadi_int size_;

    std::vector<unsigned char> buf_;
    // Sink on buf_; needed > size_ once truncated
    casadi_stats_sink s_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_STATS_RECORDER_INTERNAL_HPP
