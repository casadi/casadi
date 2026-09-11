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

#include "shf_spline_impl.hpp"
#include "casadi_misc.hpp"
#include "interpolant_impl.hpp"
#include "serializing_stream.hpp"
#include "casadi_low.hpp"
#include <limits>

namespace casadi {

  casadi_int SHFSplineFunction::get_coeff_size(casadi_int m,
      const std::vector<casadi_int>& offset) {
    casadi_int ret = m;
    for (casadi_int i=0;i<offset.size()-1;++i) ret *= offset[i+1]-offset[i];
    return ret;
  }

  void SHFSplineFunction::prepare(casadi_int m, const std::vector<casadi_int>& offset,
      casadi_int& coeffs_size, std::vector<casadi_int>& coeffs_dims,
      std::vector<casadi_int>& strides) {
    casadi_int n_dims = offset.size()-1;
    coeffs_size = get_coeff_size(m, offset);
    // The coefficients are the look-up table itself: one per grid point, no fitting.
    coeffs_dims.resize(n_dims+1);
    coeffs_dims[0] = m;
    for (casadi_int i=0;i<n_dims;++i) coeffs_dims[i+1] = offset[i+1]-offset[i];
    strides.resize(n_dims);
    strides[0] = m;
    for (casadi_int i=0;i<n_dims-1;++i) strides[i+1] = strides[i]*coeffs_dims[i+1];
  }

  void SHFSplineFunction::derive(const std::vector<double>& grid,
      const std::vector<casadi_int>& offset, std::vector<double>& inv_h,
      std::vector<casadi_int>& width, std::vector<double>& min_h) {
    casadi_int n_dims = offset.size()-1;
    inv_h.clear();
    width.resize(n_dims);
    min_h.assign(n_dims, std::numeric_limits<double>::infinity());
    for (casadi_int r=0;r<n_dims;++r) {
      casadi_int ng = offset[r+1]-offset[r];
      casadi_assert(ng>=2, "Axis " + str(r) + " needs at least two grid points.");
      const double* g = get_ptr(grid)+offset[r];
      for (casadi_int i=0;i<ng-1;++i) {
        double h = g[i+1]-g[i];
        casadi_assert(h>0, "Grid must be strictly increasing.");
        inv_h.push_back(1.0/h);
        if (h<min_h[r]) min_h[r] = h;
      }
      // Padded stencil width; an axis with only two grid points is never smoothed.
      width[r] = ng>=3 ? 3 : 2;
    }
  }

  std::vector<casadi_int> SHFSplineFunction::multi_index(casadi_int n_dims, casadi_int p) {
    // Canonical order: the leading axis counts down, so order 1 enumerates as
    // e_0, e_1, ... and jacobian columns land in axis order.
    if (n_dims==1) return std::vector<casadi_int>(1, p);
    std::vector<casadi_int> ret;
    for (casadi_int i=p; i>=0; --i) {
      std::vector<casadi_int> tail = multi_index(n_dims-1, p-i);
      for (size_t b=0; b<tail.size()/(n_dims-1); ++b) {
        ret.push_back(i);
        for (casadi_int r=0; r<n_dims-1; ++r) ret.push_back(tail[b*(n_dims-1)+r]);
      }
    }
    return ret;
  }

  std::vector<double> SHFSplineFunction::bernstein_ctrl(casadi_int k, casadi_int order) {
    // Control points of s_k: b_i = 0 for i <= k, (2i-n)/n beyond (Lemma 3.1); row p
    // holds those of s_k^(p), the p-th forward difference scaled by n!/(n-p)!, with
    // its max(0, k+1-p) leading zeros dropped so the kernel blends only the tail
    casadi_int n = 2*k+1;
    std::vector<double> b(n+1);
    for (casadi_int i=0; i<=n; ++i) b[i] = i<=k ? 0.0 : (2.0*i-n)/n;
    std::vector<double> ctrl((order+1)*(n+1), 0.0);
    double f = 1;
    for (casadi_int p=0; p<=order; ++p) {
      if (p>0) {
        for (casadi_int i=0; i<=n-p; ++i) b[i] = b[i+1]-b[i];
        f *= n-p+1;
      }
      casadi_int z = std::max(static_cast<casadi_int>(0), k+1-p);
      for (casadi_int i=z; i<=n-p; ++i) ctrl[p*(n+1)+i-z] = f*b[i];
    }
    return ctrl;
  }

  size_t SHFSplineFunction::n_iw(casadi_int n_dims, casadi_int nb) {
    // starts, per-axis weight offsets, and (nb==1 path) a compacted w_offset
    return n_dims + n_dims*nb + (n_dims+1);
  }

  size_t SHFSplineFunction::n_w(casadi_int n_dims, casadi_int k, casadi_int order,
      casadi_int nb) {
    // per-axis weights of order 0..order, one accumulator set per recursion
    // level, and the de Casteljau scratch
    size_t acc = (n_dims+1)*nb;
    if (nb==1 && static_cast<size_t>(3*n_dims) > acc) acc = 3*n_dims;
    return n_dims*(order+1)*3 + acc + k+1;
  }

  SHFSplineFunction::SHFSplineFunction(const std::string& name,
      const std::vector<double>& grid, const std::vector<casadi_int>& offset,
      const std::vector<double>& values, casadi_int k, casadi_int m, casadi_int order)
      : FunctionInternal(name), grid_(grid), offset_(offset), k_(k), m_(m),
        values_(values), parametric_(false), order_(order) {
    init_derived_members();
    casadi_assert(static_cast<casadi_int>(values_.size())==coeffs_size_,
      "Expected " + str(coeffs_size_) + " table values, got " + str(values_.size()) + ".");
  }

  SHFSplineFunction::SHFSplineFunction(const std::string& name,
      const std::vector<double>& grid, const std::vector<casadi_int>& offset,
      casadi_int k, casadi_int m, casadi_int order)
      : FunctionInternal(name), grid_(grid), offset_(offset), k_(k), m_(m),
        parametric_(true), order_(order) {
    init_derived_members();
  }

  void SHFSplineFunction::init_derived_members() {
    casadi_assert(offset_.size()>=2, "Grid must have at least one axis.");
    casadi_assert(k_>=1, "k must be at least 1, got " + str(k_) + ".");
    casadi_assert(m_>=1, "m must be at least 1, got " + str(m_) + ".");
    multi_ = multi_index(n_dims(), order_);
    prepare(m_, offset_, coeffs_size_, coeffs_dims_, strides_);
    derive(grid_, offset_, inv_h_, width_, min_h_);
    ctrl_ = bernstein_ctrl(k_, order_);
  }

  SHFSplineFunction::~SHFSplineFunction() {
    clear_mem();
  }

  const Options SHFSplineFunction::options_
  = {{&FunctionInternal::options_},
     {{"lookup_mode",
       {OT_STRINGVECTOR,
        "Specifies, for each grid dimension, the lookup algorithm used to find the "
        "correct index. 'linear' uses a forward linear search. 'exact' uses "
        "a comparator function optimized for uniformly distributed data "
        "(requires equally spaced grid points). 'binary' uses a binary search. "
        "'auto' (default) uses 'linear' for small grids and 'binary' for large."}}
     }
  };

  void SHFSplineFunction::init(const Dict& opts) {
    FunctionInternal::init(opts);

    for (auto&& op : opts) {
      if (op.first=="lookup_mode") {
        lookup_modes_ = op.second;
      }
    }

    lookup_mode_ = Interpolant::interpret_lookup_mode(lookup_modes_, grid_, offset_,
      std::vector<casadi_int>(), std::vector<casadi_int>());

    alloc_iw(n_iw(n_dims(), nb()));
    alloc_w(n_w(n_dims(), k_, order_, nb()));
  }

  Sparsity SHFSplineFunction::get_sparsity_in(casadi_int i) {
    if (i==0) return Sparsity::dense(n_dims());
    if (parametric_ && i==arg_c()) return Sparsity::dense(coeffs_size_);
    if (i==arg_eps()) return Sparsity::dense(n_dims());
    casadi_assert_dev(false);
    return Sparsity();
  }

  Sparsity SHFSplineFunction::get_sparsity_out(casadi_int i) {
    return Sparsity::dense(m_, nb());
  }

  std::string SHFSplineFunction::get_name_in(casadi_int i) {
    if (i==0) return "x";
    if (parametric_ && i==arg_c()) return "C";
    if (i==arg_eps()) return "eps";
    casadi_assert_dev(false);
    return "";
  }

  std::string SHFSplineFunction::get_name_out(casadi_int i) {
    return "f";
  }

  int SHFSplineFunction::eval(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const {
    if (!res[0]) return 0;
    const double* eps = arg[arg_eps()];
    // Neighbouring epsilon-balls must stay disjoint, otherwise the partition of
    // unity -- and with it the value -- is silently wrong.
    if (eps) {
      for (casadi_int r=0; r<n_dims(); ++r) {
        casadi_assert(eps[r]>=0 && eps[r]<0.5*min_h_[r], "epsilon (" + str(eps[r])
          + ") for axis " + str(r) + " must lie in [0, " + str(0.5*min_h_[r]) + "), half "
          "the smallest grid spacing, so that the epsilon-neighbourhoods of neighbouring "
          "grid points stay disjoint.");
      }
    }
    casadi_clear(res[0], m_*nb());
    casadi_shf_eval_multi(res[0], n_dims(), get_ptr(grid_), get_ptr(offset_),
      get_ptr(inv_h_), get_ptr(width_), get_ptr(strides_),
      parametric_ ? arg[arg_c()] : get_ptr(values_), m_, arg[0],
      eps, 1, k_, get_ptr(ctrl_), get_ptr(multi_), nb(), order_,
      get_ptr(lookup_mode_), iw, w);
    return 0;
  }

  void SHFSplineFunction::codegen_body(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_SHF_EVAL);
    g.add_auxiliary(CodeGenerator::AUX_CLEAR);
    std::string coeffs = parametric_ ? g.arg(arg_c()) : g.constant(values_);
    g << "if (" << g.res(0) << ") {\n";
    g << g.clear(g.res(0), m_*nb()) << "\n";
    g << "CASADI_PREFIX(shf_eval_multi)(" << g.res(0) << ","
      << n_dims() << "," << g.constant(grid_) << "," << g.constant(offset_) << ","
      << g.constant(inv_h_) << "," << g.constant(width_) << ","
      << g.constant(strides_) << "," << coeffs << "," << m_ << ","
      << g.arg(0) << "," << g.arg(arg_eps()) << ",1," << k_ << ","
      << g.constant(ctrl_) << "," << g.constant(multi_) << "," << nb() << "," << order_ << ","
      << g.constant(lookup_mode_) << ", iw, w);\n";
    g << "}\n";
  }

  Function SHFSplineFunction::next_order() const {
    // Every mixed partial of the next order comes from ONE traversal of the
    // coefficient stencil, so raising the derivative order costs one kernel
    // call regardless of the order or the dimension.
    std::string fname = name_ + "_der";
    Function f;
    if (!incache(fname, f)) {
      Dict opts;
      opts["lookup_mode"] = lookup_modes_;
      if (parametric_) {
        f = Function::create(new SHFSplineFunction(fname, grid_, offset_, k_, m_,
          order_+1), opts);
      } else {
        f = Function::create(new SHFSplineFunction(fname, grid_, offset_, values_,
          k_, m_, order_+1), opts);
      }
      tocache(f);
    }
    return f;
  }

  DM SHFSplineFunction::selector(casadi_int s) const {
    casadi_int nd = n_dims();
    std::vector<casadi_int> up = multi_index(nd, order_+1);
    casadi_int nup = up.size()/nd;
    DM S = DM::zeros(nup, nb());
    for (casadi_int b=0; b<nb(); ++b) {
      // locate multi_[b] + e_s among the next-order multi-indices
      for (casadi_int a=0; a<nup; ++a) {
        bool hit = true;
        for (casadi_int r=0; r<nd; ++r) {
          casadi_int want = multi_[b*nd+r] + (r==s ? 1 : 0);
          if (up[a*nd+r]!=want) { hit = false; break; }
        }
        if (hit) { S(a, b) = 1; break; }
      }
    }
    return S;
  }

  std::vector<MX> SHFSplineFunction::jac_cols(const std::vector<MX>& arg) const {
    MX up = next_order()(arg)[0];
    std::vector<MX> cols(n_dims());
    for (casadi_int s=0; s<n_dims(); ++s) cols[s] = mtimes(up, MX(selector(s)));
    return cols;
  }

  Function SHFSplineFunction::get_forward(casadi_int nfwd, const std::string& name,
      const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const Dict& opts) const {
    std::vector<MX> arg = mx_in();
    std::vector<MX> cols = jac_cols(arg);
    std::vector<MX> in = arg;
    in.push_back(MX::sym(inames[n_in_], Sparsity(size_out(0))));
    for (casadi_int i=0; i<n_in_; ++i) {
      in.push_back(MX::sym(inames[n_in_+1+i], repmat(sparsity_in(i), 1, nfwd)));
    }
    // Only x carries a seed
    std::vector<MX> seed = horzsplit(in[n_in_+1]);
    std::vector<MX> sens(nfwd);
    for (casadi_int d=0; d<nfwd; ++d) {
      MX r = cols[0]*seed[d](0);
      for (casadi_int s=1; s<n_dims(); ++s) r = r + cols[s]*seed[d](s);
      sens[d] = r;
    }
    Dict wopts = opts;
    wopts["always_inline"] = true;
    return Function(name, in, {horzcat(sens)}, inames, onames, wopts);
  }

  Function SHFSplineFunction::get_reverse(casadi_int nadj, const std::string& name,
      const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const Dict& opts) const {
    std::vector<MX> arg = mx_in();
    std::vector<MX> cols = jac_cols(arg);
    std::vector<MX> in = arg;
    in.push_back(MX::sym(inames[n_in_], Sparsity(size_out(0))));
    in.push_back(MX::sym(inames[n_in_+1], repmat(sparsity_out(0), 1, nadj)));
    std::vector<MX> seed = horzsplit(in[n_in_+1], nb());
    std::vector<MX> sens(nadj);
    for (casadi_int d=0; d<nadj; ++d) {
      std::vector<MX> parts(n_dims());
      for (casadi_int s=0; s<n_dims(); ++s) parts[s] = dot(cols[s], seed[d]);
      sens[d] = vertcat(parts);
    }
    // Only x gets a sensitivity; the others are structurally empty
    std::vector<MX> out(n_in_);
    out[0] = horzcat(sens);
    for (casadi_int j=1; j<n_in_; ++j) out[j] = MX(size1_in(j), size2_in(j)*nadj);
    Dict wopts = opts;
    wopts["always_inline"] = true;
    return Function(name, in, out, inames, onames, wopts);
  }

  Function SHFSplineFunction::get_jacobian(const std::string& name,
      const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const Dict& opts) const {
    std::vector<MX> arg = mx_in();
    // Only x is differentiated; the other blocks are structurally empty
    std::vector<MX> cols = jac_cols(arg);
    for (auto& c : cols) c = vec(c);
    std::vector<MX> jac_out(n_in_);
    jac_out[0] = horzcat(cols);
    for (casadi_int j=1; j<n_in_; ++j) jac_out[j] = MX(m_*nb(), nnz_in(j));
    std::vector<MX> jac_in = arg;
    jac_in.push_back(mx_out(0));
    Dict wopts = opts;
    wopts["always_inline"] = true;
    return Function(name, jac_in, jac_out, inames, onames, wopts);
  }

  void SHFSplineFunction::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);
    s.version("SHFSplineFunction", 1);
    s.pack("SHFSplineFunction::grid", grid_);
    s.pack("SHFSplineFunction::offset", offset_);
    s.pack("SHFSplineFunction::k", k_);
    s.pack("SHFSplineFunction::m", m_);
    s.pack("SHFSplineFunction::values", values_);
    s.pack("SHFSplineFunction::parametric", parametric_);
    s.pack("SHFSplineFunction::order", order_);
    s.pack("SHFSplineFunction::lookup_modes", lookup_modes_);
    s.pack("SHFSplineFunction::lookup_mode", lookup_mode_);
  }

  SHFSplineFunction::SHFSplineFunction(DeserializingStream& s) : FunctionInternal(s) {
    s.version("SHFSplineFunction", 1);
    s.unpack("SHFSplineFunction::grid", grid_);
    s.unpack("SHFSplineFunction::offset", offset_);
    s.unpack("SHFSplineFunction::k", k_);
    s.unpack("SHFSplineFunction::m", m_);
    s.unpack("SHFSplineFunction::values", values_);
    s.unpack("SHFSplineFunction::parametric", parametric_);
    s.unpack("SHFSplineFunction::order", order_);
    s.unpack("SHFSplineFunction::lookup_modes", lookup_modes_);
    s.unpack("SHFSplineFunction::lookup_mode", lookup_mode_);
    init_derived_members();
  }

  ProtoFunction* SHFSplineFunction::deserialize(DeserializingStream& s) {
    return new SHFSplineFunction(s);
  }

  Function shf_spline(const std::string& name,
      const std::vector< std::vector<double> >& grid, const std::vector<double>& values,
      casadi_int k, casadi_int m, const Dict& opts) {
    std::vector<casadi_int> offset;
    std::vector<double> stacked;
    Interpolant::stack_grid(grid, offset, stacked);
    return Function::create(new SHFSplineFunction(name, stacked, offset, values, k, m), opts);
  }

  Function shf_spline(const std::string& name,
      const std::vector< std::vector<double> >& grid,
      casadi_int k, casadi_int m, const Dict& opts) {
    std::vector<casadi_int> offset;
    std::vector<double> stacked;
    Interpolant::stack_grid(grid, offset, stacked);
    return Function::create(new SHFSplineFunction(name, stacked, offset, k, m), opts);
  }

} // namespace casadi
