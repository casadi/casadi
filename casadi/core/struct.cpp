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


#define CASADI_STRUCT_CPP
#include "struct.hpp"
#include "mx_node.hpp"
#include <algorithm>
#include <set>

namespace casadi {

namespace {
    bool is_slice(const GenericType& key) {
      return key.is_string() && key.as_string().find(':')!=std::string::npos;
    }

    bool is_index(const GenericType& key) {
      return key.is_int() || (key.is_double() && key.as_double()==std::floor(key.as_double()));
    }

    casadi_int canonical(casadi_int i, casadi_int n, bool ind1) {
      if (ind1) {
        casadi_assert(i>=1 && i<=n, "Index " + str(i) + " out of range [1, " + str(n) + "].");
        return i-1;
      }
      casadi_assert(i>=-n && i<n, "Index " + str(i) + " out of range [" + str(-n) + ", "
        + str(n) + ").");
      return i<0 ? i+n : i;
    }

    Slice parse_slice(const std::string& s) {
      std::vector<std::string> part(1);
      for (char c : s) {
        if (c==':') {
          part.emplace_back();
        } else {
          part.back() += c;
        }
      }
      casadi_assert(part.size()<=3, "Invalid slice '" + s + "'.");
      part.resize(3);
      casadi_int step = part[2].empty() ? 1 : std::stoll(part[2]);
      casadi_int start = part[0].empty() ? std::numeric_limits<casadi_int>::min()
        : std::stoll(part[0]);
      casadi_int stop = part[1].empty() ? std::numeric_limits<casadi_int>::max()
        : std::stoll(part[1]);
      return Slice(start, stop, step);
    }

    bool is_all(const Slice& s) {
      return (s.start==0 || s.start==std::numeric_limits<casadi_int>::min())
        && s.stop==std::numeric_limits<casadi_int>::max() && s.step==1;
    }

    // Positions selected by a path element along a dimension of length n
    bool selection(const GenericType& key, casadi_int n, bool ind1,
        std::vector<casadi_int>& sel) {
      sel.clear();
      if (is_index(key)) {
        sel.push_back(canonical(key.to_int(), n, ind1));
      } else if (key.is_int_vector()) {
        for (casadi_int i : key.as_int_vector()) sel.push_back(canonical(i, n, ind1));
      } else if (key.is_double_vector()) {
        for (double i : key.as_double_vector()) sel.push_back(canonical(i, n, ind1));
      } else if (is_slice(key)) {
        sel = parse_slice(key.as_string()).all(n);
      } else {
        return false;
      }
      return true;
    }

    // Stack the nonzeros of each part
    template<class T>
    T veccat_nz(const std::vector<T>& v) {
      std::vector<T> ret = {T(0, 1)};
      for (auto&& e : v) ret.push_back(T::sparsity_cast(e, Sparsity::dense(e.nnz(), 1)));
      return vertcat(ret);
    }

    // Readable path, 0-based like labels and symbol names
    std::string path_str(const std::vector<GenericType>& p, bool ind1=false) {
      std::stringstream ss;
      for (casadi_int k=0; k<p.size(); ++k) {
        if (p[k].is_string() && !is_slice(p[k])) {
          ss << (k==0 ? "" : ".") << p[k].as_string();
        } else if (ind1 && is_index(p[k])) {
          ss << "[" << p[k].to_int()-1 << "]";
        } else if (ind1 && (p[k].is_int_vector() || p[k].is_double_vector())) {
          std::vector<casadi_int> v;
          if (p[k].is_int_vector()) {
            v = p[k].as_int_vector();
          } else {
            for (double d : p[k].as_double_vector()) v.push_back(static_cast<casadi_int>(d));
          }
          for (casadi_int& i : v) i--;
          ss << "[" << str(v) << "]";
        } else {
          ss << "[" << p[k] << "]";
        }
      }
      return ss.str();
    }

    casadi_int nonzero_count(const DM& x) {
      return std::count_if(x->begin(), x->end(), [](double v) { return v!=0;});
    }
    casadi_int nonzero_count(const SX& x) {
      return std::count_if(x->begin(), x->end(), [](const SXElem& v) { return !v.is_zero();});
    }
    casadi_int nonzero_count(const MX&) {
      return 0;
    }

    // Nonzeros of a column, cast to a sparsity
    template<class M>
    M nz_select(const M& x, const Sparsity& sp, const std::vector<casadi_int>& nz) {
      M r;
      x.get_nz(r, false, IM(nz));
      return M::sparsity_cast(r, sp);
    }
    MX nz_select(const MX& x, const Sparsity& sp, const std::vector<casadi_int>& nz) {
      return x->get_nzref(sp, nz);
    }

    // Flat indices of a matrix at an offset
    IM nz_range(casadi_int offset, const Sparsity& sp) {
      return IM(sp, range(offset, offset+sp.nnz()));
    }
}  // namespace

  Struct::Struct() : is_leaf_(false), symmetric_(false), nnz_(0) {
  }

  Struct::Struct(const std::vector<std::string>& names) : Struct() {
    for (auto&& n : names) add(n);
  }

  Struct Struct::leaf(const Sparsity& sp, bool symmetric) {
    casadi_assert(!symmetric || sp.is_square(),
      "Symmetric matrix must be square, got " + sp.dim() + ".");
    Struct ret;
    ret.is_leaf_ = true;
    ret.symmetric_ = symmetric;
    ret.sp_ = symmetric ? sp + sp.T() : sp;
    ret.nnz_ = ret.stored().nnz();
    return ret;
  }

  Struct Struct::matrix(const Struct& rows, const Struct& cols, bool symmetric) {
    Struct ret = leaf(Sparsity::dense(rows.nnz(), cols.nnz()), symmetric);
    if (!rows.is_leaf()) ret.rows_ = std::make_shared<Struct>(rows);
    if (!cols.is_leaf()) ret.cols_ = std::make_shared<Struct>(cols);
    return ret;
  }

  Struct Struct::matrix(const Struct& rows, casadi_int ncol) {
    return matrix(rows, leaf(Sparsity::dense(ncol, 1)));
  }

  Sparsity Struct::stored() const {
    return symmetric_ ? Sparsity::triu(sp_) : sp_;
  }

  void Struct::add(const std::string& name, const Struct& s,
      const std::vector<casadi_int>& repeat) {
    casadi_assert(!is_leaf_, "Cannot add entries to a matrix.");
    casadi_assert(!name.empty() && name.find(':')==std::string::npos,
      "Invalid entry name '" + name + "'.");
    casadi_assert(find(name)<0, "Duplicate entry '" + name + "'.");
    Entry e;
    e.name = name;
    e.repeat = repeat;
    e.s = std::make_shared<Struct>(s);
    e.n = repeat.empty() ? 1 : repeat[0];
    e.size = s.nnz();
    for (casadi_int i=0; i<repeat.size(); ++i) {
      casadi_assert(repeat[i]>=0, "Negative repetition for '" + name + "'.");
      if (i>0) e.size *= repeat[i];
    }
    groups_.push_back({static_cast<casadi_int>(entries_.size())});
    entries_.push_back(e);
    update();
  }

  void Struct::add(const std::string& name, const Sparsity& sp,
      const std::vector<casadi_int>& repeat) {
    add(name, leaf(sp), repeat);
  }

  void Struct::add(const std::string& name, casadi_int nrow, casadi_int ncol,
      const std::vector<casadi_int>& repeat) {
    add(name, Sparsity::dense(nrow, ncol), repeat);
  }

  void Struct::interleave(const std::vector<std::string>& names) {
    std::vector<casadi_int> merged;
    for (auto&& n : names) {
      merged.push_back(find(n));
      casadi_assert(merged.back()>=0, "No entry '" + n + "'.");
    }
    std::vector<std::vector<casadi_int> > groups;
    bool placed = false;
    for (auto&& g : groups_) {
      if (std::find(merged.begin(), merged.end(), g.front())==merged.end()) {
        groups.push_back(g);
      } else {
        casadi_assert(g.size()==1, "Entry '" + entries_[g.front()].name
          + "' is already interleaved.");
        if (!placed) groups.push_back(merged);
        placed = true;
      }
    }
    groups_ = groups;
    update();
  }

  void Struct::update() {
    start_.clear();
    block_entry_.clear();
    block_rep_.clear();
    nnz_ = 0;
    for (auto&& g : groups_) {
      casadi_int n = 0;
      for (casadi_int e : g) {
        n = std::max(n, entries_[e].n);
        entries_[e].start.resize(entries_[e].n);
      }
      for (casadi_int i=0; i<n; ++i) {
        for (casadi_int e : g) {
          if (i>=entries_[e].n) continue;
          entries_[e].start[i] = nnz_;
          start_.push_back(nnz_);
          block_entry_.push_back(e);
          block_rep_.push_back(i);
          nnz_ += entries_[e].size;
        }
      }
    }
  }

  casadi_int Struct::find(const std::string& name) const {
    for (casadi_int e=0; e<entries_.size(); ++e) {
      if (entries_[e].name==name) return e;
    }
    return -1;
  }

  Matrix<casadi_int> Struct::index(const std::vector<GenericType>& path, bool ind1) const {
    IM ret = get<IM>(path, ind1, nz_range);
    if (ind1) for (casadi_int& i : ret.nonzeros()) i++;
    return ret;
  }

  template<class T>
  T Struct::get(const std::vector<GenericType>& path, bool ind1, const Leaf<T>& leaf) const {
    try {
      return get(0, path, 0, ind1, leaf);
    } catch (std::exception& e) {
      casadi_error("Cannot index " + get_str() + " with path " + path_str(path, ind1) + ":\n"
        + e.what());
    }
  }

  template<class T>
  T Struct::get(casadi_int offset, const std::vector<GenericType>& p, casadi_int k,
      bool ind1, const Leaf<T>& leaf) const {
    if (is_leaf_) return get_leaf(offset, p, k, ind1, leaf);
    std::vector<std::string> sel;
    if (k==p.size()) {
      std::vector<T> ret;
      for_each_leaf([&](casadi_int i, const Sparsity& sp, const std::vector<GenericType>&) {
        ret.push_back(leaf(offset+i, sp));
      });
      return veccat_nz(ret);
    } else if (p[k].is_string() && !is_slice(p[k])) {
      casadi_int e = find(p[k].as_string());
      casadi_assert(e>=0, "No entry '" + p[k].as_string() + "', expected one of "
        + str(names()) + ".");
      std::vector<casadi_int> rep;
      return get(e, offset, p, k+1, rep, ind1, leaf);
    } else if (p[k].is_string_vector()) {
      sel = p[k].as_string_vector();
    } else if (is_slice(p[k]) && is_all(parse_slice(p[k].as_string()))) {
      sel = names();
    } else {
      casadi_error("Expected an entry name, list of names or ':', got " + str(p[k]) + ".");
    }
    std::vector<T> ret;
    std::vector<GenericType> q = p;
    for (auto&& n : sel) {
      q[k] = n;
      ret.push_back(get(offset, q, k, ind1, leaf));
    }
    return veccat_nz(ret);
  }

  template<class T>
  T Struct::get(casadi_int e, casadi_int offset, const std::vector<GenericType>& p,
      casadi_int k, std::vector<casadi_int>& rep, bool ind1, const Leaf<T>& leaf) const {
    const Entry& en = entries_[e];
    if (rep.size()<en.repeat.size()) {
      // Repetitions not addressed are all selected
      std::vector<casadi_int> sel;
      if (k<p.size() && selection(p[k], en.repeat[rep.size()], ind1, sel)) {
        k++;
      } else {
        sel = range(en.repeat[rep.size()]);
      }
      std::vector<T> ret;
      for (casadi_int i : sel) {
        rep.push_back(i);
        ret.push_back(get(e, offset, p, k, rep, ind1, leaf));
        rep.pop_back();
      }
      return horzcat(ret);
    }
    offset += en.start.at(rep.empty() ? 0 : rep[0]);
    casadi_int lin = 0;
    for (casadi_int d=1; d<rep.size(); ++d) lin = lin*en.repeat[d] + rep[d];
    return en.s->get(offset + lin*en.s->nnz(), p, k, ind1, leaf);
  }

  template<class T>
  T Struct::get_leaf(casadi_int offset, const std::vector<GenericType>& p, casadi_int k,
      bool ind1, const Leaf<T>& leaf) const {
    T ret = leaf(offset, stored());
    if (symmetric_) ret = triu2symm(ret);
    casadi_int nk = p.size()-k;
    casadi_assert(nk<=2, "Expected at most 2 indices for a matrix, got " + str(nk) + ".");
    if (nk==0) return ret;
    T r;
    if (nk==1) {
      ret.get(r, false, dim_index(p[k], sp_.size2()==1 ? rows_ : nullptr, sp_.numel(), ind1));
    } else {
      ret.get(r, false, dim_index(p[k], rows_, sp_.size1(), ind1),
        dim_index(p[k+1], cols_, sp_.size2(), ind1));
    }
    return r;
  }

  IM Struct::dim_index(const GenericType& key, const std::shared_ptr<const Struct>& labels,
      casadi_int n, bool ind1) {
    std::vector<casadi_int> sel;
    if (!selection(key, n, ind1, sel)) {
      casadi_assert(labels, "Matrix has no labels, cannot index with " + str(key) + ".");
      sel = labels->index(std::vector<GenericType>{key}).nonzeros();
    }
    return IM(sel);
  }

  std::vector<GenericType> Struct::path(casadi_int i, bool ind1) const {
    std::vector<GenericType> ret;
    path(canonical(i, nnz_, ind1), ret);
    if (ind1) {
      for (auto&& e : ret) if (e.is_int()) e = e.as_int() + 1;
    }
    return ret;
  }

  void Struct::path(casadi_int i, std::vector<GenericType>& p) const {
    if (is_leaf_) {
      Sparsity sp = stored();
      if (!sp_.is_scalar()) p.push_back(sp.row(i) + sp.get_col().at(i)*sp.size1());
      return;
    }
    casadi_int b = std::upper_bound(start_.begin(), start_.end(), i) - start_.begin() - 1;
    const Entry& e = entries_[block_entry_[b]];
    p.push_back(e.name);
    i -= start_[b];
    casadi_int lin = i / e.s->nnz();
    i %= e.s->nnz();
    std::vector<casadi_int> rep(e.repeat.size());
    for (casadi_int d=rep.size()-1; d>0; --d) {
      rep[d] = lin % e.repeat[d];
      lin /= e.repeat[d];
    }
    if (!rep.empty()) rep[0] = block_rep_[b];
    for (casadi_int r : rep) p.push_back(r);
    e.s->path(i, p);
  }

  std::vector<std::string> Struct::labels() const {
    std::vector<std::string> ret(nnz_);
    for (casadi_int i=0; i<nnz_; ++i) ret[i] = path_str(path(i, false));
    return ret;
  }

  std::vector<std::string> Struct::names() const {
    std::vector<std::string> ret;
    for (auto&& e : entries_) ret.push_back(e.name);
    return ret;
  }

  std::vector<casadi_int> Struct::repeat(const std::string& name) const {
    casadi_int e = find(name);
    casadi_assert(e>=0, "No entry '" + name + "'.");
    return entries_[e].repeat;
  }

  Struct Struct::child(const std::string& name) const {
    casadi_int e = find(name);
    casadi_assert(e>=0, "No entry '" + name + "'.");
    return *entries_[e].s;
  }

  const Sparsity& Struct::sparsity() const {
    casadi_assert(is_leaf_, "Not a matrix.");
    return sp_;
  }

  void Struct::disp(std::ostream& stream, bool more) const {
    if (is_leaf_) {
      if (symmetric_) stream << "symm(";
      if (rows_) {
        stream << rows_->get_str() << "x" << (cols_ ? cols_->get_str() : str(sp_.size2()));
      } else {
        stream << sp_.dim(!sp_.is_dense());
      }
      if (symmetric_) stream << ")";
      return;
    }
    stream << "{";
    for (casadi_int e=0; e<entries_.size(); ++e) {
      if (e>0) stream << ", ";
      stream << entries_[e].name;
      if (!entries_[e].repeat.empty()) stream << str(entries_[e].repeat);
      stream << ": " << entries_[e].s->get_str();
    }
    stream << "}";
  }

  void Struct::for_each_leaf(const LeafVisitor& f) const {
    std::vector<GenericType> p;
    for_each_leaf(0, p, f);
  }

  void Struct::for_each_leaf(casadi_int offset, std::vector<GenericType>& p,
      const LeafVisitor& f) const {
    if (is_leaf_) return f(offset, stored(), p);
    for (casadi_int b=0; b<start_.size(); ++b) {
      const Entry& e = entries_[block_entry_[b]];
      casadi_int cs = e.s->nnz();
      if (cs==0) continue;
      for (casadi_int lin=0; lin<e.size/cs; ++lin) {
        casadi_int np = p.size();
        p.push_back(e.name);
        if (!e.repeat.empty()) {
          p.push_back(block_rep_[b]);
          p.resize(np+1+e.repeat.size());
          for (casadi_int d=e.repeat.size()-1, r=lin; d>0; --d) {
            p[np+1+d] = r % e.repeat[d];
            r /= e.repeat[d];
          }
        }
        e.s->for_each_leaf(offset + start_[b] + lin*cs, p, f);
        p.resize(np);
      }
    }
  }

  std::string Struct::slice_str(const Slice& s) {
    std::stringstream ss;
    if (s.start!=0 && s.start!=std::numeric_limits<casadi_int>::min()) ss << s.start;
    ss << ":";
    if (s.stop!=std::numeric_limits<casadi_int>::max()) ss << s.stop;
    if (s.step!=1) ss << ":" << s.step;
    return ss.str();
  }

  template<class M>
  StructValue<M>::StructValue(const Struct& s, const M& data) : s_(s), cached_(false) {
    M flat;
    if (data.is_scalar()) {
      flat = M(Sparsity::dense(s.nnz(), 1), densify(data));
    } else {
      casadi_assert(data.numel()==s.nnz() && (data.is_vector() || data.is_empty()),
        "Expected a vector of length " + str(s.nnz()) + ", got " + data.dim() + ".");
      flat = densify(vec(data));
    }
    s.for_each_leaf([&](casadi_int i, const Sparsity& sp, const std::vector<GenericType>&) {
      offset_.push_back(i);
      leaves_.push_back(nz_select(flat, sp, range(i, i+sp.nnz())));
    });
  }

  template<class M>
  StructValue<M> StructValue<M>::sym(const Struct& s, const std::string& prefix) {
    StructValue<M> ret(s, 0);
    casadi_int k = 0;
    s.for_each_leaf([&](casadi_int, const Sparsity& sp, const std::vector<GenericType>& p) {
      std::string name = prefix;
      for (auto&& e : p) {
        if (!name.empty()) name += "_";
        name += e.is_string() ? e.as_string() : str(e.to_int());
      }
      ret.leaves_[k++] = M::sym(name, sp);
    });
    return ret;
  }

  template<class M>
  casadi_int StructValue<M>::leaf(casadi_int offset) const {
    return std::lower_bound(offset_.begin(), offset_.end(), offset) - offset_.begin();
  }

  template<class M>
  const M& StructValue<M>::cat() const {
    if (!cached_) cat_ = veccat_nz(leaves_);
    cached_ = true;
    return cat_;
  }

  template<class M>
  M StructValue<M>::get(const std::vector<GenericType>& path, bool ind1) const {
    return s_.get<M>(path, ind1, [this](casadi_int offset, const Sparsity& sp) {
      return sp.nnz()==0 ? M::zeros(sp) : leaves_.at(leaf(offset));
    });
  }

  template<class M>
  void StructValue<M>::set(const std::vector<GenericType>& path, const M& value, bool ind1) {
    assign(s_.get<IM>(path, ind1, nz_range), value, path_str(path, ind1));
  }

  template<class M>
  void StructValue<M>::set_nz(const IM& ind, const M& value, bool ind1) {
    std::vector<casadi_int> nz = ind.nonzeros();
    for (casadi_int& i : nz) i = canonical(i, s_.nnz(), ind1);
    assign(IM(ind.sparsity(), nz), value, "the indexed part");
  }

  template<class M>
  void StructValue<M>::assign(const IM& ind, const M& value, const std::string& what) {
    M v = value;
    if (v.is_scalar()) {
      v = M(ind.sparsity(), densify(v));
    } else if (v.size1()==1 && v.size2()==ind.size2() && ind.size1()>1) {
      v = repmat(v, ind.size1(), 1);
    } else if (v.size1()==ind.size1() && v.size2()>0 && v.size2()<ind.size2()
        && ind.size2() % v.size2()==0) {
      v = repmat(v, 1, ind.size2()/v.size2());
    } else if (v.size()!=ind.size() && v.is_vector() && ind.is_vector()
        && v.numel()==ind.numel()) {
      v = v.T();
    }
    casadi_assert(v.size()==ind.size(), "Cannot assign " + value.dim(false) + " to "
      + what + ", which is " + ind.dim(false) + ".");
    if (v.sparsity()!=ind.sparsity()) {
      M p = project(v, ind.sparsity());
      casadi_assert(nonzero_count(v)==nonzero_count(p), "Cannot assign nonzeros outside the "
        "sparsity pattern of " + what + ".");
      v = p;
    }
    v = M::sparsity_cast(v, Sparsity::dense(v.nnz(), 1));
    // Positions and local nonzeros per matrix; symmetric matrices appear twice: last one wins
    const std::vector<casadi_int>& nz = ind.nonzeros();
    std::map<casadi_int, std::pair<std::vector<casadi_int>, std::vector<casadi_int> > > parts;
    std::set<casadi_int> seen;
    for (casadi_int k=nz.size()-1; k>=0; --k) {
      if (!seen.insert(nz[k]).second) continue;
      casadi_int l = std::upper_bound(offset_.begin(), offset_.end(), nz[k]) - offset_.begin() - 1;
      parts[l].first.push_back(k);
      parts[l].second.push_back(nz[k] - offset_[l]);
    }
    for (auto&& e : parts) {
      std::vector<casadi_int>& pos = e.second.first;
      std::vector<casadi_int>& loc = e.second.second;
      std::reverse(pos.begin(), pos.end());
      std::reverse(loc.begin(), loc.end());
      M& l = leaves_[e.first];
      if (loc==range(l.nnz())) {
        l = nz_select(v, l.sparsity(), pos);
      } else {
        l.set_nz(nz_select(v, Sparsity::dense(pos.size(), 1), pos), false, IM(loc));
      }
    }
    cached_ = false;
  }

  template<class M>
  std::string StructValue<M>::type_name() {
    return "Struct" + M::type_name();
  }

  template<class M>
  void StructValue<M>::disp(std::ostream& stream, bool more) const {
    stream << s_.get_str() << ": " << cat().get_str(more);
  }

  typedef std::vector<GenericType> Path;
  template CASADI_EXPORT IM Struct::get(const Path&, bool, const Leaf<IM>&) const;
  template CASADI_EXPORT DM Struct::get(const Path&, bool, const Leaf<DM>&) const;
  template CASADI_EXPORT SX Struct::get(const Path&, bool, const Leaf<SX>&) const;
  template CASADI_EXPORT MX Struct::get(const Path&, bool, const Leaf<MX>&) const;
  template class CASADI_EXPORT StructValue<DM>;
  template class CASADI_EXPORT StructValue<SX>;
  template class CASADI_EXPORT StructValue<MX>;

} // namespace casadi
