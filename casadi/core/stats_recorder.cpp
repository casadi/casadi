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


#include "stats_recorder_internal.hpp"
#include "casadi_misc.hpp"
#include "filesystem_impl.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

namespace casadi {

  // ---------- StatsRecorder ----------

  StatsRecorder::StatsRecorder() : StatsRecorder(Dict()) {
  }

  StatsRecorder::StatsRecorder(StatsRecorderInternal* node) {
    own(node);
  }

  StatsRecorder::StatsRecorder(const Dict& opts) {
    casadi_int size = 1 << 17;
    for (auto&& op : opts) {
      if (op.first=="size") {
        size = static_cast<casadi_int>(op.second.to_double());
      } else {
        casadi_error("Unknown StatsRecorder option: " + op.first);
      }
    }
    casadi_assert(size>0, "StatsRecorder size must be positive");
    own(new StatsRecorderInternal(size));
  }

namespace {
    // A .json file (any case) is for the call tree (export_json), not the stream
    bool is_json(const std::string& filename) {
      if (filename.size() < 5) return false;
      std::string ext = filename.substr(filename.size() - 5);
      for (auto& ch : ext) ch = std::tolower(static_cast<unsigned char>(ch));
      return ext==".json";
    }
}  // namespace

  StatsRecorder StatsRecorder::load(const std::string& filename) {
    casadi_assert(!is_json(filename),
      "Cannot load '" + filename + "': JSON is an export, load the stream (e.g. .cbor)");
    auto is = Filesystem::ifstream_ptr(filename, std::ios::binary);
    std::stringstream ss;
    ss << is->rdbuf();
    return from_bytes(ss.str());
  }

  StatsRecorder StatsRecorder::from_bytes(const std::string& bytes) {
    StatsRecorder ret = create(new StatsRecorderInternal(std::max<casadi_int>(bytes.size(), 1)));
    ret->set_bytes(bytes);
    return ret;
  }

  StatsRecorder StatsRecorder::create(StatsRecorderInternal* node) {
    return StatsRecorder(node);
  }

  bool StatsRecorder::test_cast(const SharedObjectInternal* ptr) {
    return dynamic_cast<const StatsRecorderInternal*>(ptr)!=nullptr;
  }

  StatsRecorderInternal* StatsRecorder::get() const {
    return static_cast<StatsRecorderInternal*>(SharedObject::get());
  }

  const StatsRecorderInternal* StatsRecorder::operator->() const {
    return static_cast<const StatsRecorderInternal*>(SharedObject::operator->());
  }

  StatsRecorderInternal* StatsRecorder::operator->() {
    return static_cast<StatsRecorderInternal*>(SharedObject::operator->());
  }

  casadi_int StatsRecorder::nbytes() const {
    return (*this)->nbytes();
  }

  bool StatsRecorder::truncated() const {
    return (*this)->truncated() || (*this)->decode().truncated;
  }

  void StatsRecorder::save(const std::string& filename) const {
    casadi_assert(!is_json(filename),
      "Cannot save '" + filename + "': JSON is export_json, save the stream (e.g. .cbor)");
    auto os = Filesystem::ofstream_ptr(filename, std::ios::binary);
    *os << (*this)->bytes();
  }

  // ---------- Sink ----------

  StatsRecorderInternal::StatsRecorderInternal(casadi_int size)
      : size_(size), buf_(size), needed_(0) {
    s_.data = this;
    s_.reserve = reserve;
  }

  void StatsRecorderInternal::disp(std::ostream& stream, bool more) const {
    stream << "StatsRecorder(" << nbytes() << " of " << size_ << " bytes)";
  }

  unsigned char* StatsRecorderInternal::reserve(casadi_stats_sink* s, casadi_int len,
      casadi_int* pos) {
    // Concurrent writers get disjoint windows from one atomic add
    auto* me = static_cast<StatsRecorderInternal*>(s->data);
    *pos = me->needed_.fetch_add(len);
    if (*pos + len <= me->size_) return me->buf_.data() + *pos;
    // Recording stops; a CBOR break marks the cut
    if (*pos < me->size_) me->buf_[*pos] = 0xff;
    return nullptr;
  }

  void StatsRecorderInternal::reset() {
    needed_ = 0;
  }

  casadi_int StatsRecorderInternal::nbytes() const {
    return needed_;
  }

  bool StatsRecorderInternal::truncated() const {
    return needed_ > size_;
  }

  std::string StatsRecorderInternal::bytes() const {
    return std::string(reinterpret_cast<const char*>(buf_.data()),
                       std::min<casadi_int>(needed_, size_));
  }

  void StatsRecorderInternal::set_bytes(const std::string& b) {
    casadi_assert_dev(static_cast<casadi_int>(b.size()) <= size_);
    std::copy(b.begin(), b.end(), buf_.begin());
    needed_ = b.size();
  }

  casadi_int StatsRecorderInternal::put_record(const StatsRecord& r, casadi_int ref) {
    casadi_int pos;
    unsigned char* p = casadi_stats_reserve_record(&s_, r.n, static_cast<casadi_stats_kind>(r.kind),
      ref, r.tail.size(), &pos);
    if (p) std::copy(r.tail.begin(), r.tail.end(), p);
    return pos;
  }

  // ---------- Decoding ----------

namespace {

    // Arrays of one kind become the matching native vector
    GenericType homogenize(const std::vector<GenericType>& e) {
      if (e.empty()) return std::vector<double>();
      bool all_int = true, all_num = true, all_str = true;
      for (auto&& ek : e) {
        all_int = all_int && ek.is_int();
        all_num = all_num && (ek.is_int() || ek.is_double());
        all_str = all_str && ek.is_string();
      }
      if (all_int) {
        std::vector<casadi_int> r;
        for (auto&& ek : e) r.push_back(ek.as_int());
        return r;
      }
      if (all_num) {
        std::vector<double> r;
        for (auto&& ek : e) r.push_back(ek.to_double());
        return r;
      }
      if (all_str) {
        std::vector<std::string> r;
        for (auto&& ek : e) r.push_back(ek.as_string());
        return r;
      }
      return e;
    }

    // An item of what casadi_cbor.hpp writes, as a native type: the position after it, or null
    const unsigned char* read_item(const unsigned char* p, const unsigned char* end,
        GenericType& v) {
      const unsigned char *q, *t;
      casadi_int i, n;
      double d;
      int b;
      if ((q = casadi_cbor_read_int(p, end, &i))) {
        v = i;
      } else if ((q = casadi_cbor_read_real(p, end, &d))) {
        v = d;
      } else if ((q = casadi_cbor_read_text(p, end, &t, &n))) {
        v = std::string(reinterpret_cast<const char*>(t), n);
      } else if ((q = casadi_cbor_read_bool(p, end, &b))) {
        v = static_cast<bool>(b);
      } else if ((q = casadi_cbor_read_null(p, end))) {
        v = GenericType();
      } else if ((q = casadi_cbor_read_array(p, end, &n))) {
        // Each element takes a byte at least
        if (n > end - q) return nullptr;
        std::vector<GenericType> e(n);
        for (auto&& ek : e) {
          q = read_item(q, end, ek);
          if (!q) return nullptr;
        }
        v = homogenize(e);
      }
      return q;
    }

}  // namespace

  StatsTree StatsRecorderInternal::decode() const {
    std::string b = bytes();
    const unsigned char* s = reinterpret_cast<const unsigned char*>(b.data());
    const unsigned char* end = s + b.size();
    StatsTree t;
    // Node by the address of the record that created it
    std::map<casadi_int, casadi_int> by_offset;
    std::vector<std::map<std::string, std::vector<GenericType> > > appended;
    // Record by record, up to the end of the last complete one
    const unsigned char* p = s;
    while (p < end) {
      casadi_int start = p - s;
      StatsRecord rec;
      const unsigned char* q = casadi_stats_read_head(p, end, &rec.n, &rec.kind, &rec.ref);
      // Each item takes a byte at least
      if (!q || rec.n - 2 > end - q) break;
      const unsigned char* tail = q;
      std::vector<GenericType> rest(rec.n - 2);
      for (auto&& e : rest) {
        q = read_item(q, end, e);
        if (!q) break;
      }
      if (!q) break;
      rec.tail = std::string(reinterpret_cast<const char*>(tail), q - tail);
      casadi_int ri = t.records.size();
      // The node the record references, if any
      auto it = by_offset.find(rec.ref);
      casadi_int owner = it==by_offset.end() ? -1 : it->second;
      if (rec.ref>=0 && owner<0) break;
      if (rec.kind==CASADI_STATS_BEGIN_CALL || rec.kind==CASADI_STATS_BEGIN_SECTION
          || rec.kind==CASADI_STATS_BEGIN_ITERATION) {
        StatsNode c;
        if (rec.kind==CASADI_STATS_BEGIN_CALL) {
          // [1, parent|null, id, mem|null]
          if (rest.size()!=2 || !rest[0].is_string()) break;
          if (!rest[1].is_int() && rest[1].getType()!=OT_NULL) break;
          c.id = rest[0].as_string();
          // After the first colon, as casadi_stats_id_is
          c.name = c.id.substr(c.id.find(':') + 1);
          c.mem = rest[1];
        } else if (rec.kind==CASADI_STATS_BEGIN_SECTION) {
          // [5, parent, name]
          if (rest.size()!=1 || !rest[0].is_string() || owner<0) break;
          c.kind = StatsNode::SECTION;
          c.name = rest[0].as_string();
        } else {
          // [7, parent, index]
          if (rest.size()!=1 || !rest[0].is_int() || owner<0) break;
          c.kind = StatsNode::ITERATION;
          c.index = rest[0].as_int();
        }
        c.begin = ri;
        c.parent = owner;
        casadi_int ci = t.nodes.size();
        if (owner>=0) {
          t.nodes[owner].events.push_back(-1 - ci);
        } else {
          t.roots.push_back(ci);
        }
        by_offset[start] = ci;
        t.nodes.push_back(c);
        appended.emplace_back();
      } else {
        if (owner<0) break;
        StatsNode& c = t.nodes[owner];
        if (rec.kind==CASADI_STATS_END_CALL) {
          // [2, call, flag]
          if (rest.size()!=1 || !rest[0].is_int() || c.kind!=StatsNode::CALL) break;
          c.ended = true;
          c.flag = rest[0].to_int();
          c.end = ri;
        } else if (rec.kind==CASADI_STATS_END_SCOPE) {
          // [6, node]
          if (!rest.empty() || c.kind==StatsNode::CALL) break;
          c.ended = true;
          c.end = ri;
        } else {
          if (rec.kind==CASADI_STATS_SET || rec.kind==CASADI_STATS_APPEND) {
            // [3|4, node, key, value]
            if (rest.size()!=2 || !rest[0].is_string()) break;
            if (rec.kind==CASADI_STATS_SET) {
              c.stats[rest[0].as_string()] = rest[1];
            } else {
              appended[owner][rest[0].as_string()].push_back(rest[1]);
            }
          } else if (rec.kind==CASADI_STATS_DECLARE_FIELDS) {
            // [8, call, [names]]; an empty array decodes as a numeric vector
            if (rest.size()!=1) break;
            if (rest[0].is_string_vector()) {
              c.fields = rest[0].as_string_vector();
            } else if (!rest[0].is_double_vector() || !rest[0].as_double_vector().empty()) {
              break;
            }
            c.has_fields = true;
          } else if (rec.kind==CASADI_STATS_SET_FIELD) {
            // [9, iteration, field, value]
            if (rest.size()!=2 || !rest[0].is_int() || rest[0].as_int()<0) break;
            if (c.kind!=StatsNode::ITERATION) break;
            c.field_values[rest[0].as_int()] = rest[1];
          } else {
            break;
          }
          c.events.push_back(ri);
        }
      }
      t.records.push_back(rec);
      p = q;
    }
    t.truncated = p < end;

    for (size_t k=0; k<t.nodes.size(); ++k) {
      StatsNode& c = t.nodes[k];
      for (auto&& e : appended[k]) c.stats[e.first] = e.second;
      // The fields of the iterations as columns, and their indices as "iter"
      if (c.has_fields) c.stats["iterations"] = iteration_columns(t, c);
    }
    return t;
  }

  Dict StatsRecorderInternal::iteration_columns(const StatsTree& t, const StatsNode& c) {
    std::vector<const StatsNode*> its;
    for (casadi_int e : c.events) {
      if (e<0 && t.nodes[-1-e].kind==StatsNode::ITERATION) its.push_back(&t.nodes[-1-e]);
    }
    std::stable_sort(its.begin(), its.end(),
      [](const StatsNode* a, const StatsNode* b) { return a->index < b->index; });
    // Fields beyond the declared ones are named by their index
    std::vector<std::string> names = c.fields;
    for (const StatsNode* it : its) {
      for (auto&& v : it->field_values) {
        for (casadi_int j=names.size(); j<=v.first; ++j) names.push_back("field" + str(j));
      }
    }
    Dict cols;
    std::vector<casadi_int> iter;
    for (const StatsNode* it : its) iter.push_back(it->index);
    cols["iter"] = iter;
    for (casadi_int j=0; j<static_cast<casadi_int>(names.size()); ++j) {
      // One kind for a column, where its values allow: a missing value is NaN if numeric
      std::vector<GenericType> col;
      bool all_int = true, all_num = true, all_str = true, all_bool = true, missing = false;
      for (const StatsNode* it : its) {
        auto v = it->field_values.find(j);
        if (v==it->field_values.end()) {
          missing = true;
          col.push_back(GenericType());
          continue;
        }
        col.push_back(v->second);
        all_int = all_int && v->second.is_int();
        all_num = all_num && (v->second.is_int() || v->second.is_double());
        all_str = all_str && v->second.is_string();
        all_bool = all_bool && v->second.is_bool();
      }
      if (all_int && !missing) {
        std::vector<casadi_int> r;
        for (auto&& e : col) r.push_back(e.as_int());
        cols[names[j]] = r;
      } else if (all_num) {
        std::vector<double> r;
        for (auto&& e : col) r.push_back(e.getType()==OT_NULL ? nan : e.to_double());
        cols[names[j]] = r;
      } else if (all_str && !missing) {
        std::vector<std::string> r;
        for (auto&& e : col) r.push_back(e.as_string());
        cols[names[j]] = r;
      } else if (all_bool && !missing) {
        std::vector<bool> r;
        for (auto&& e : col) r.push_back(e.as_bool());
        cols[names[j]] = r;
      } else {
        cols[names[j]] = col;
      }
    }
    return cols;
  }

namespace {

    // Events in presentation order: iterations by index
    std::vector<casadi_int> ordered_events(const StatsTree& t, const StatsNode& c) {
      std::vector<casadi_int> ev = c.events, slots, its;
      for (size_t k=0; k<ev.size(); ++k) {
        if (ev[k]<0 && t.nodes[-1-ev[k]].kind==StatsNode::ITERATION) {
          slots.push_back(k);
          its.push_back(ev[k]);
        }
      }
      std::stable_sort(its.begin(), its.end(), [&](casadi_int a, casadi_int b) {
        return t.nodes[-1-a].index < t.nodes[-1-b].index;
      });
      for (size_t k=0; k<slots.size(); ++k) ev[slots[k]] = its[k];
      return ev;
    }

    Dict native(const StatsTree& t, casadi_int ci);

    // The nodes below a node: calls under children, sections under their name, iterations by
    // index under iterations
    void native_contents(const StatsTree& t, const StatsNode& c, Dict& d) {
      std::vector<Dict> children;
      std::vector<Dict> its;
      std::map<std::string, std::vector<Dict> > secs;
      for (casadi_int e : c.events) {
        if (e>=0) continue;
        const StatsNode& ch = t.nodes[-1-e];
        if (ch.kind==StatsNode::CALL) {
          children.push_back(native(t, -1-e));
        } else if (ch.kind==StatsNode::SECTION) {
          secs[ch.name].push_back(native(t, -1-e));
        } else if (ch.index>=0) {
          // A cut stream may miss iterations: those stay empty
          if (ch.index>=static_cast<casadi_int>(its.size())) {
            its.resize(ch.index + 1, Dict{{"children", std::vector<Dict>()}});
          }
          its[ch.index] = native(t, -1-e);
        }
      }
      for (auto&& s : secs) {
        if (s.second.size()==1) {
          d[s.first] = s.second[0];
        } else {
          d[s.first] = s.second;
        }
      }
      if (!its.empty()) d["iterations"] = its;
      if (!children.empty() || c.kind!=StatsNode::CALL || (its.empty() && secs.empty())) {
        d["children"] = children;
      }
    }

    Dict native(const StatsTree& t, casadi_int ci) {
      const StatsNode& c = t.nodes[ci];
      Dict d;
      if (c.kind==StatsNode::CALL) {
        d["name"] = c.name;
        d["id"] = c.id;
        d["mem"] = c.mem;
        d["flag"] = c.ended ? GenericType(c.flag) : GenericType();
        d["stats"] = c.stats;
      } else if (!c.stats.empty()) {
        d["stats"] = c.stats;
      }
      native_contents(t, c, d);
      return d;
    }

    void json_string(std::ostream& os, const std::string& s) {
      os << '"';
      for (unsigned char ch : s) {
        if (ch=='"' || ch=='\\') {
          os << '\\' << ch;
        } else if (ch < 0x20) {
          os << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(ch)
             << std::dec;
        } else {
          os << ch;
        }
      }
      os << '"';
    }

    void json_number(std::ostream& os, double v) {
      if (std::isfinite(v)) {
        os << std::setprecision(std::numeric_limits<double>::max_digits10) << v;
      } else {
        os << "null";
      }
    }

    void json(std::ostream& os, const GenericType& v) {
      if (v.is_bool()) {
        os << (v.as_bool() ? "true" : "false");
      } else if (v.is_int()) {
        os << v.as_int();
      } else if (v.is_double()) {
        json_number(os, v.as_double());
      } else if (v.is_string()) {
        json_string(os, v.as_string());
      } else if (v.is_dict()) {
        os << '{';
        bool first = true;
        for (auto&& e : v.as_dict()) {
          if (!first) os << ',';
          first = false;
          json_string(os, e.first);
          os << ':';
          json(os, e.second);
        }
        os << '}';
      } else if (v.is_dict_vector() || v.is_vector() || v.is_int_vector()
                 || v.is_double_vector() || v.is_string_vector() || v.is_bool_vector()) {
        std::vector<GenericType> e;
        if (v.is_dict_vector()) {
          for (auto&& d : v.as_dict_vector()) e.push_back(d);
        } else if (v.is_vector()) {
          e = v.as_vector();
        } else if (v.is_int_vector()) {
          for (auto x : v.as_int_vector()) e.push_back(x);
        } else if (v.is_double_vector()) {
          for (auto x : v.as_double_vector()) e.push_back(x);
        } else if (v.is_string_vector()) {
          for (auto&& x : v.as_string_vector()) e.push_back(x);
        } else {
          for (bool x : v.as_bool_vector()) e.push_back(x);
        }
        os << '[';
        for (size_t k=0; k<e.size(); ++k) {
          if (k) os << ',';
          json(os, e[k]);
        }
        os << ']';
      } else {
        os << "null";
      }
    }

    void deinterleave_node(const StatsTree& t, casadi_int ci, casadi_int parent,
        StatsRecorderInternal* out) {
      const StatsNode& c = t.nodes[ci];
      casadi_int me = out->put_record(t.records[c.begin], parent);
      for (casadi_int e : ordered_events(t, c)) {
        if (e<0) {
          deinterleave_node(t, -1-e, me, out);
        } else {
          out->put_record(t.records[e], me);
        }
      }
      if (c.end>=0) out->put_record(t.records[c.end], me);
    }

}  // namespace

  GenericType StatsRecorder::to_native() const {
    StatsTree t = (*this)->decode();
    std::vector<Dict> roots;
    for (casadi_int r : t.roots) roots.push_back(native(t, r));
    return roots;
  }

  std::string StatsRecorder::to_json() const {
    std::stringstream ss;
    json(ss, Dict{{"version", 1}, {"calls", to_native()}});
    return ss.str();
  }

  void StatsRecorder::export_json(const std::string& filename) const {
    auto os = Filesystem::ofstream_ptr(filename);
    *os << to_json();
  }

  GenericType StatsRecorder::get_stat(const std::string& fname,
      const std::string& key) const {
    StatsTree t = (*this)->decode();
    // The last call by that name; which call to pick when there are several is still open
    for (auto it = t.nodes.rbegin(); it != t.nodes.rend(); ++it) {
      if (it->kind!=StatsNode::CALL || it->name!=fname) continue;
      auto s = it->stats.find(key);
      casadi_assert(s!=it->stats.end(),
        "No stat '" + key + "' for the last call of '" + fname + "'");
      return s->second;
    }
    casadi_error("No call of '" + fname + "' recorded");
  }

  StatsRecorder StatsRecorder::deinterleave() const {
    StatsTree t = (*this)->decode();
    // References may need more bytes once reordered: room for 8 more per record
    StatsRecorder ret = create(new StatsRecorderInternal(
      std::max<casadi_int>((*this)->size(), (*this)->nbytes() + 8 * t.records.size())));
    for (casadi_int r : t.roots) deinterleave_node(t, r, -1, ret.get());
    return ret;
  }

} // namespace casadi
