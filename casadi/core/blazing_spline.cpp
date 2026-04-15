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


#include "blazing_spline_impl.hpp"
#include "interpolant_impl.hpp"
#include "bspline_impl.hpp"
#include "casadi_misc.hpp"
#include "serializer.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

namespace casadi {

  static std::vector<casadi_int> knot_offsets(const std::vector<casadi_int>& knot_dims) {
    std::vector<casadi_int> offsets(knot_dims.size() + 1);
    offsets[0] = 0;
    for (size_t i = 0; i < knot_dims.size(); ++i)
      offsets[i + 1] = offsets[i] + knot_dims[i];
    return offsets;
  }

  template<typename M>
  static M compute_inv(const M& K, const std::vector<casadi_int>& offsets) {
    casadi_int nd = offsets.size() - 1;
    std::vector<M> inv_parts;
    for (casadi_int d = 0; d < nd; ++d) {
      casadi_int off = offsets[d];
      casadi_int n_k = offsets[d + 1] - off;
      M t = K(Slice(off, off + n_k));
      for (casadi_int span = 1; span <= 3; ++span) {
        M inv_span;
        if (span < n_k) {
          M diff = t(Slice(span, n_k)) - t(Slice(0, n_k - span));
          diff = vertcat(diff, M::zeros(span, 1));
          inv_span = if_else_zero(diff, 1.0 / (diff + 1e-100));
        } else {
          inv_span = M::zeros(n_k, 1);
        }
        inv_parts.push_back(inv_span);
      }
    }
    return vertcat(inv_parts);
  }

  // MX version of BSplineCommon::derivative_coeff for parametric knots.
  // Builds the bidiagonal transformation matrix T from symbolic knots
  // and applies it to the coefficient tensor along axis i.
  static MX derivative_coeff_mx(casadi_int i,
      const std::vector<MX>& knots_per_dim,
      const std::vector<casadi_int>& degree,
      const std::vector<casadi_int>& coeffs_dims,
      const MX& coeffs) {
    casadi_int n_dims = degree.size();
    casadi_int n_k = knots_per_dim[i].size1();
    casadi_int n = n_k - degree[i] - 1;

    MX K_i = knots_per_dim[i];
    MX delta_knots = K_i(range(1+degree[i], n_k-1))
          - K_i(range(1, n_k-degree[i]-1));
    MX d = static_cast<double>(degree[i]) / delta_knots;    // length n-1

    // T = diag(-d) + upper_band(+d) is a scaled finite-difference operator.
    // Apply via slice-subtract + broadcast-multiply — no T, no kron, no densify.
    std::vector<casadi_int> coeffs_dims_new = coeffs_dims;
    coeffs_dims_new[i+1] = n - 1;

    casadi_int L = 1, R = 1;
    for (casadi_int k=0; k<=i; ++k) L *= coeffs_dims[k];
    for (casadi_int k=i+2; k<(casadi_int)coeffs_dims.size(); ++k) R *= coeffs_dims[k];
    casadi_int K = coeffs_dims[i+1];
    casadi_int Kp = n - 1;

    MX M_coeffs = reshape(coeffs, L*K, R);
    MX top = M_coeffs(Slice(L,   L*K),    Slice());
    MX bot = M_coeffs(Slice(0, L*(K-1)),  Slice());
    MX diffed = top - bot;

    std::vector<casadi_int> dims{L, Kp, R};
    std::vector<casadi_int> a{-1, -2, -3};
    std::vector<casadi_int> b{-2};
    std::vector<casadi_int> c{-1, -2, -3};
    return MX::einstein(vec(diffed), d,
      dims, std::vector<casadi_int>{Kp}, dims,
      a, b, c);
  }

  Function blazing_spline(const std::string& name,
      const std::vector< std::vector<double> >& knots,
      const Dict& opts) {
    return Function::create(new BlazingSplineFunction(name, knots, 0), opts);
  }

  Function blazing_spline(const std::string& name,
      const std::vector<casadi_int>& knot_dims,
      const Dict& opts) {
    bool precompute_coeff = false, precompute_grid = false;
    auto it = opts.find("precompute_coeff");
    if (it != opts.end()) precompute_coeff = it->second;
    it = opts.find("precompute_grid");
    if (it != opts.end()) precompute_grid = it->second;
    bool use_inv = precompute_grid;
    Function F_inner = Function::create(
      new BlazingSplineFunction(name, knot_dims, 0,
        precompute_coeff, precompute_grid, use_inv), opts);
    if (!use_inv) return F_inner;

    // Wrap: user sees (x, C, knots), wrapper computes inv and calls inner
    std::vector<casadi_int> offsets = knot_offsets(knot_dims);
    casadi_int nd = knot_dims.size();
    casadi_int nk = offsets.back();
    MX x = MX::sym("x", nd);
    MX C = MX::sym("C", F_inner.size_in(1));
    MX knots = MX::sym("knots", nk);
    MX inv = compute_inv(knots, offsets);
    std::vector<MX> ret = F_inner(std::vector<MX>{x, C, knots, inv});
    return Function(name, {x, C, knots}, ret,
                    {"x", "C", "knots"}, F_inner.name_out(),
                    {{"always_inline", true}});
  }

  size_t BlazingSplineFunction::get_n_in() {
    casadi_int n = 2; // x, C
    if (precompute_coeff_) n += diff_order_; // dC, ddC
    if (has_parametric_knots()) n += 1; // knots
    if (inv_input_) n += 1; // inv
    return n;
  }
  size_t BlazingSplineFunction::get_n_out() {
    return diff_order_+1;
  }

  bool BlazingSplineFunction::get_diff_in(casadi_int i) {
    return i==0;
  }

  Sparsity BlazingSplineFunction::get_sparsity_in(casadi_int i) {
    if (i==0) {
      return Sparsity::dense(ndim());
    } else if (i==1) {
      return Sparsity::dense(nc_);
    } else if (has_parametric_knots() && i==arg_knots()) {
      return Sparsity::dense(knots_offset_.back());
    } else if (inv_input_ && i==arg_inv()) {
      return Sparsity::dense(3 * knots_offset_.back());
    } else if (precompute_coeff_ && i==2+has_parametric_knots()+inv_input_) {
      return Sparsity::dense(ndc_);
    } else if (precompute_coeff_ && i==3+has_parametric_knots()+inv_input_) {
      return Sparsity::dense(nddc_);
    } else {
      casadi_assert_dev(false);
      return Sparsity();
    }
  }
  Sparsity BlazingSplineFunction::get_sparsity_out(casadi_int i) {
    if (i==0) {
      return Sparsity::dense(1, 1);
    } else if (i==1) {
      return Sparsity::dense(1, ndim());
    } else if (i==2) {
      return Sparsity::dense(ndim(), ndim());
    } else {
      casadi_assert_dev(false);
      return Sparsity();
    }
  }

  std::string BlazingSplineFunction::get_name_in(casadi_int i) {
    if (i==0) {
      return "x";
    } else if (i==1) {
      return "C";
    } else if (has_parametric_knots() && i==arg_knots()) {
      return "knots";
    } else if (inv_input_ && i==arg_inv()) {
      return "inv";
    } else if (precompute_coeff_ && i==2+has_parametric_knots()+inv_input_) {
      return "dC";
    } else if (precompute_coeff_ && i==3+has_parametric_knots()+inv_input_) {
      return "ddC";
    } else {
      casadi_assert_dev(false);
      return "";
    }
  }
  std::string BlazingSplineFunction::get_name_out(casadi_int i) {
    if (i==0) {
      return "f";
    } else if (i==1) {
      return "g";
    } else if (i==2) {
      return "h";
    } else {
      casadi_assert_dev(false);
      return "";
    }
  }

  BlazingSplineFunction::BlazingSplineFunction(const std::string& name,
    const std::vector< std::vector<double> >& knots,
    casadi_int diff_order,
    bool precompute_coeff,
    bool precompute_grid) : FunctionInternal(name), diff_order_(diff_order),
    precompute_coeff_(precompute_coeff), precompute_grid_(precompute_grid),
    knots_(knots) {

    init_derived_members();

    casadi_assert(knots.size()>=1, "blazing_spline only defined for 1D-5D");
    casadi_assert(knots.size()<=5, "blazing_spline only defined for 1D-5D");
  }

  BlazingSplineFunction::BlazingSplineFunction(const std::string& name,
    const std::vector<casadi_int>& knot_dims,
    casadi_int diff_order,
    bool precompute_coeff,
    bool precompute_grid,
    bool inv_input) : FunctionInternal(name), diff_order_(diff_order),
    precompute_coeff_(precompute_coeff), precompute_grid_(precompute_grid),
    inv_input_(inv_input) {
    // knots_ left empty → has_parametric_knots() returns true
    // Build knots_offset_ from dimension sizes
    knots_offset_.resize(knot_dims.size()+1);
    knots_offset_[0] = 0;
    for (size_t i=0; i<knot_dims.size(); ++i) {
      knots_offset_[i+1] = knots_offset_[i] + knot_dims[i];
    }

    init_derived_members();

    casadi_assert(knot_dims.size()>=1, "blazing_spline only defined for 1D-5D");
    casadi_assert(knot_dims.size()<=5, "blazing_spline only defined for 1D-5D");
  }

  void BlazingSplineFunction::init_derived_members() {
    // For non-parametric knots, stack the grid to compute offsets
    if (!has_parametric_knots()) {
      Interpolant::stack_grid(knots_, knots_offset_, knots_stacked_);
    }

    casadi_int nd = ndim();

    // Compute coefficient tensor size
    nc_ = 1;
    for (casadi_int i=0; i<nd; ++i) {
      nc_ *= (knots_offset_[i+1] - knots_offset_[i]) - 4;
    }

    // Compute derivative coefficient tensor size
    ndc_ = 0;
    for (casadi_int k=0;k<nd;++k) {
      casadi_int ndc = 1;
      for (casadi_int i=0;i<nd;++i) {
        ndc *= (knots_offset_[i+1] - knots_offset_[i]) - 4 - (i==k);
      }
      ndc_+= ndc;
    }

    nddc_ = 0;
    for (casadi_int k=0;k<nd;++k) {
      for (casadi_int kk=0;kk<nd;++kk) {
        casadi_int ndc = 1;
        for (casadi_int i=0;i<nd;++i) {
          ndc *= (knots_offset_[i+1] - knots_offset_[i]) - 4 - (i==k)-(i==kk);
        }
        // We only need the triangular part
        if (kk>=k) {
          nddc_+= ndc;
        }
      }
    }

    // Precompute reciprocal knot spans (only when knot values are known)
    if (!has_parametric_knots()) {
      DM inv_dm = compute_inv(DM(knots_stacked_), knots_offset_);
      knots_inv_ = inv_dm.nonzeros();
    }
  }

  const Options BlazingSplineFunction::options_
  = {{&FunctionInternal::options_},
     {{"precompute_coeff",
       {OT_BOOL,
        "If true, derivative evaluation requires precomputed derivative "
        "coefficient tensors (dC, ddC) as function inputs. Only supported "
        "up to 3D. Default: true for fixed knots, false for parametric knots."}},
      {"precompute_grid",
       {OT_BOOL,
        "If true, precompute reciprocal knot spans to replace runtime "
        "divisions with multiplications. For parametric knots, inv is "
        "computed symbolically from the knots input. Default: false."}},
      {"lookup_mode",
       {OT_STRINGVECTOR,
        "Specifies, for each grid dimension, the lookup algorithm used to find the "
        "correct index. 'linear' uses a forward linear search. 'exact' uses "
        "a comparator function optimized for uniformly distributed data "
        "(requires equally spaced knots). 'binary' uses a binary search. "
        "'auto' (default) uses 'linear' for small grids and 'binary' for large."}}
     }
  };

  void BlazingSplineFunction::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="precompute_coeff") {
        precompute_coeff_ = op.second;
      } else if (op.first=="precompute_grid") {
        precompute_grid_ = op.second;
      } else if (op.first=="lookup_mode") {
        lookup_modes_ = op.second;
      }
    }

    casadi_int n_dims = ndim();

    if (precompute_coeff_) {
      casadi_assert(n_dims<=3,
        "blazing_spline with precompute_coeff=true only supports up to 3D. "
        "Use precompute_coeff=false for 4D/5D.");
    }

    // Arrays for holding inputs and outputs
    alloc_iw(4*n_dims+2);
    alloc_w(n_dims+1);
  }

  BlazingSplineFunction::~BlazingSplineFunction() {
    clear_mem();
  }

  void BlazingSplineFunction::codegen_body(CodeGenerator& g) const {
    casadi_int nd = ndim();
    switch (nd) {
      case 1: g.add_auxiliary(CodeGenerator::AUX_BLAZING_1D_BOOR_EVAL); break;
      case 2: g.add_auxiliary(CodeGenerator::AUX_BLAZING_2D_BOOR_EVAL); break;
      case 3: g.add_auxiliary(CodeGenerator::AUX_BLAZING_3D_BOOR_EVAL); break;
      case 4: g.add_auxiliary(CodeGenerator::AUX_BLAZING_4D_BOOR_EVAL); break;
      case 5: g.add_auxiliary(CodeGenerator::AUX_BLAZING_5D_BOOR_EVAL); break;
      default: casadi_assert_dev(false);
    }
    g.add_include("simde/x86/avx2.h");
    g.add_include("simde/x86/fma.h");

    std::string knots_offset = g.constant(knots_offset_);
    std::string knots_stacked = has_parametric_knots() ?
      g.arg(arg_knots()) : g.constant(knots_stacked_);
    std::string knots_inv;
    if (inv_input_) {
      knots_inv = g.arg(arg_inv());
    } else if (!has_parametric_knots() && precompute_grid_) {
      knots_inv = g.constant(knots_inv_);
    } else {
      knots_inv = "0";
    }

    std::vector<casadi_int> degree(nd, 3);
    std::vector<casadi_int> mode =
      Interpolant::interpret_lookup_mode(
        lookup_modes_, knots_stacked_, knots_offset_, degree, degree);

    std::string fun_name = "casadi_blazing_" + str(nd) + "d_boor_eval";
    std::string f_ptr = "res[0]";
    std::string J_ptr = (diff_order_>=1) ? "res[1]" : "0";
    std::string H_ptr = (diff_order_>=2) ? "res[2]" : "0";

    std::string dc_ptr = "0", ddc_ptr = "0";
    if (precompute_coeff_) {
      casadi_int dc_idx = 2 + has_parametric_knots() + inv_input_;
      if (diff_order_>=1) dc_ptr = g.arg(dc_idx);
      if (diff_order_>=2) ddc_ptr = g.arg(dc_idx+1);
    }

    g << fun_name + "(" + f_ptr + ", " + J_ptr + ", " + H_ptr + ", " +
          knots_stacked + ", " +
          knots_inv + ", " +
          knots_offset + ", " +
          "arg[1], " + dc_ptr + ", " + ddc_ptr + ", " +
          "arg[0], " +
          g.constant(mode) + ", " +
          "iw, w);\n";
  }

  bool BlazingSplineFunction::has_jacobian() const {
    return diff_order_<2;
  }

  Function BlazingSplineFunction::get_jacobian(const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    casadi_int N = ndim();
    bool parametric = has_parametric_knots();
    casadi_int nk = parametric ? knots_offset_.back() : 0;
    casadi_int n_inv = 3 * nk;

    MX x = MX::sym("x", N);
    MX C = MX::sym("C", nc_);
    MX knots_sym;
    if (parametric) knots_sym = MX::sym("knots", nk);
    MX inv_sym;
    if (inv_input_) inv_sym = MX::sym("inv", n_inv);

    Dict Jopts = combine(jacobian_options_, der_options_);
    Jopts = combine(opts, Jopts);
    Jopts = combine(Jopts, generate_options("jacobian"));
    Jopts["derivative_of"] = self();

    std::string fJname = name_ + "_der";

    // --- Synthesize dC/ddC tensors (coeff mode only) ---
    std::vector<casadi_int> coeffs_dims(N+1);
    coeffs_dims[0] = 1;
    for (casadi_int i=0; i<N; ++i) {
      coeffs_dims[i+1] = knots_offset_[i+1]-knots_offset_[i]-4;
    }
    std::vector<casadi_int> degree(N, 3);

    std::vector<MX> dCv;
    MX dC, ddC;
    // Per-dim degree after one derivative
    std::vector< std::vector<casadi_int> > degree_d(N);
    // Numeric derivative knots (non-parametric path); filled by derivative_coeff
    std::vector< std::vector< std::vector<double> > > knots_d_num(N);
    // Parametric per-dim knot vectors
    std::vector<MX> K_per_dim;

    if (precompute_coeff_) {
      if (parametric) {
        K_per_dim.resize(N);
        for (casadi_int i=0; i<N; ++i) {
          casadi_int off = knots_offset_[i];
          K_per_dim[i] = knots_sym(Slice(off, off + (knots_offset_[i+1]-off)));
        }
      }
      for (casadi_int i=0; i<N; ++i) {
        if (parametric) {
          dCv.push_back(derivative_coeff_mx(i, K_per_dim, degree, coeffs_dims, C));
        } else {
          dCv.push_back(BSplineCommon::derivative_coeff(
            i, knots_, degree, coeffs_dims, C, knots_d_num[i], degree_d[i]));
        }
        degree_d[i].assign(N, 3);
        degree_d[i][i] = 2;
      }
      dC = vertcat(dCv);

      if (diff_order_>=1) {
        // ddC ordering (must match runtime layout expected by 2d/3d_boor_eval):
        //   diagonals  (i,i) for i in [0,N),
        //   off-diags  (0,1) for N==2;  (0,1), (1,2), (2,0) for N==3.
        std::vector<std::pair<casadi_int, casadi_int>> dd_pairs;
        for (casadi_int i=0; i<N; ++i) dd_pairs.emplace_back(i, i);
        if (N==2) {
          dd_pairs.emplace_back(0, 1);
        } else if (N==3) {
          dd_pairs.emplace_back(0, 1);
          dd_pairs.emplace_back(1, 2);
          dd_pairs.emplace_back(2, 0);
        }

        std::vector<MX> parts;
        parts.reserve(dd_pairs.size());
        std::vector< std::vector<double> > knots_dummy;
        std::vector<casadi_int> degree_dummy;
        for (auto& p : dd_pairs) {
          casadi_int di = p.first, dj = p.second;
          std::vector<casadi_int> cd = coeffs_dims;
          cd[di+1] -= 1;
          if (parametric) {
            std::vector<MX> Kd(N);
            for (casadi_int k=0; k<N; ++k) {
              casadi_int n_ki = knots_offset_[k+1]-knots_offset_[k];
              Kd[k] = (k==di) ? K_per_dim[k](Slice(1, n_ki-1)) : K_per_dim[k];
            }
            parts.push_back(derivative_coeff_mx(dj, Kd, degree_d[di], cd, dCv[di]));
          } else {
            parts.push_back(BSplineCommon::derivative_coeff(
              dj, knots_d_num[di], degree_d[di], cd, dCv[di],
              knots_dummy, degree_dummy));
          }
        }
        ddC = vertcat(parts);
      }
    }

    // --- Create child function fJ (diff_order_+1) ---
    Function fJ;
    if (!incache(fJname, fJ)) {
      if (parametric) {
        std::vector<casadi_int> kdims(N);
        for (casadi_int i=0; i<N; ++i)
          kdims[i] = knots_offset_[i+1]-knots_offset_[i];
        fJ = Function::create(
          new BlazingSplineFunction(fJname, kdims, diff_order_+1,
            precompute_coeff_, precompute_grid_, /*inv_input=*/precompute_grid_), Jopts);
      } else {
        fJ = Function::create(
          new BlazingSplineFunction(fJname, knots_, diff_order_+1,
            precompute_coeff_, precompute_grid_), Jopts);
      }
      tocache(fJ);
    }

    // --- Child inputs: [x, C, [knots], [inv], [dC], [ddC]] ---
    std::vector<MX> in_child = {x, C};
    if (parametric) in_child.push_back(knots_sym);
    if (precompute_grid_ && parametric) {
      MX inv_mx = inv_input_ ? inv_sym : compute_inv(knots_sym, knots_offset_);
      in_child.push_back(inv_mx);
    }
    if (precompute_coeff_) {
      in_child.push_back(dC);
      if (diff_order_ >= 1) in_child.push_back(ddC);
    }

    std::vector<MX> ret = fJ(in_child);

    // --- User-facing jacobian inputs (mirror original function inputs) ---
    std::vector<MX> jac_in = {x, C};
    std::vector<casadi_int> in_sizes = {N, nc_};
    if (parametric)     { jac_in.push_back(knots_sym);        in_sizes.push_back(nk); }
    if (inv_input_)     { jac_in.push_back(inv_sym);          in_sizes.push_back(n_inv); }
    if (precompute_coeff_ && diff_order_>=1) {
      jac_in.push_back(MX(1, ndc_));                          in_sizes.push_back(ndc_);
    }

    // --- Jacobian outputs: for each orig output k, for each in_user, a block ---
    std::vector<MX> jac_out;
    for (casadi_int k=0; k<=diff_order_; ++k) {
      casadi_int nrows = 1;
      for (casadi_int j=0; j<k; ++j) nrows *= N;
      for (size_t j=0; j<jac_in.size(); ++j) {
        jac_out.push_back(j==0 ? ret[k+1] : MX(nrows, in_sizes[j]));
      }
    }

    // --- Append adjoint seeds (one per original output) ---
    for (casadi_int k=0; k<=diff_order_; ++k) {
      if (k==0) jac_in.push_back(MX(1, 1));
      else if (k==1) jac_in.push_back(MX(1, N));
      else if (k==2) jac_in.push_back(MX(N, N));
    }

    return Function(name, jac_in, jac_out, inames, onames, {{"always_inline", true}});
  }

  void BlazingSplineFunction::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);

    s.version("BlazingSplineFunction", 2);
    s.pack("BlazingSplineFunction::diff_order", diff_order_);
    s.pack("BlazingSplineFunction::precompute_coeff", precompute_coeff_);
    s.pack("BlazingSplineFunction::precompute_grid", precompute_grid_);
    s.pack("BlazingSplineFunction::knots", knots_);
    s.pack("BlazingSplineFunction::lookup_modes", lookup_modes_);
    s.pack("BlazingSplineFunction::parametric_knots", has_parametric_knots());
    if (has_parametric_knots()) {
      s.pack("BlazingSplineFunction::knots_offset", knots_offset_);
      s.pack("BlazingSplineFunction::inv_input", inv_input_);
    }
  }

  BlazingSplineFunction::BlazingSplineFunction(DeserializingStream & s) : FunctionInternal(s) {
    int v = s.version("BlazingSplineFunction", 1, 2);
    s.unpack("BlazingSplineFunction::diff_order", diff_order_);
    if (v>=2) {
      s.unpack("BlazingSplineFunction::precompute_coeff", precompute_coeff_);
      s.unpack("BlazingSplineFunction::precompute_grid", precompute_grid_);
    } else {
      precompute_coeff_ = true;
      precompute_grid_ = false;
    }
    s.unpack("BlazingSplineFunction::knots", knots_);
    if (v>=2) {
      s.unpack("BlazingSplineFunction::lookup_modes", lookup_modes_);
      bool parametric;
      s.unpack("BlazingSplineFunction::parametric_knots", parametric);
      if (parametric) {
        s.unpack("BlazingSplineFunction::knots_offset", knots_offset_);
        s.unpack("BlazingSplineFunction::inv_input", inv_input_);
      }
    }
    init_derived_members();
  }

  ProtoFunction* BlazingSplineFunction::deserialize(DeserializingStream& s) {
    return new BlazingSplineFunction(s);
  }

  class BlazingSplineIncrementalSerializer {
    public:

    BlazingSplineIncrementalSerializer() : serializer(ss) {
    }

    std::string generate_id(const std::vector<MX>& a) {
      ref.insert(ref.end(), a.begin(), a.end());
      if (a.empty()) return "";

      std::vector<MX> ordered = Function::order(a);
      // First serialize may introduce unknown dependencies (e.g. sparsity)
      // and hence definitions
      // Subsequent serialization will have references instead.
      // In order to still get a match with a later common subexpression,
      // make sure that all dependencies are already defined.
      serializer.pack(ordered);
      ss.str("");
      ss.clear();
      serializer.pack(ordered);
      std::string ret = ss.str();
      ss.str("");
      ss.clear();
      return ret;
    }

    private:
      std::stringstream ss;
      // List of references to keep alive
      std::vector<MX> ref;
      SerializingStream serializer;
  };

  void BlazingSplineFunction::merge(const std::vector<MX>& arg,
      std::vector<MX>& subs_from,
      std::vector<MX>& subs_to) const {

    Function base = self();
    for (casadi_int i=0;i<diff_order_;++i) {
      base = base->derivative_of_;
    }

    // Sort graph
    Function f("f", {}, arg, {{"allow_free", true}, {"max_io", 0}});


    std::unordered_map<std::string, std::vector<MX> > targets0;
    std::unordered_map<std::string, std::vector<MX> > targets1;
    std::vector<MX> targets2;

    BlazingSplineIncrementalSerializer ss;
    std::string key;

    // Loop over instructions
    for (int k=0; k<f.n_instructions(); ++k) {
      MX e = f.instruction_MX(k);
      if (e.is_call()) {
        Function fun = e.which_function();

        // Check if the function is a BlazingSplineFunction
        if (fun.class_name()=="BlazingSplineFunction") {
          key = ss.generate_id(e->dep_);
          // Which derivative level?
          if (fun==base) {
            targets0[key].push_back(e);
          } else if (!fun->derivative_of_.is_null() &&
                     fun->derivative_of_==base) {
            targets1[key].push_back(e);
          } else if (!fun->derivative_of_.is_null() &&
                    !fun->derivative_of_->derivative_of_.is_null() &&
                    fun->derivative_of_->derivative_of_==base) {
            targets2.push_back(e);
          }
        }
      }
    }

    // Loop over second order targets, targets2
    for (const auto& e : targets2) {

      // Compute key that matches targets1
      // Precompute: strip last arg (ddC) to match targets1's (x, C, dC)
      // NPC: all levels share the same deps (x, C), use directly
      key = precompute_coeff_ ?
        ss.generate_id(vector_init(e->dep_)) :
        ss.generate_id(e->dep_);

      // Loop over all matching target1 entries
      for (const auto& ee : targets1[key]) {
        // Mark all matches for substitution
        subs_from.push_back(ee);
        // Substitute with self
        subs_to.push_back(e);
      }

      // Compute key that matches targets0
      // Precompute coeff: strip two args (ddC, dC) to match targets0's (x, C)
      // Parametric grid: strip inv to match targets0's (x, C, K)
      // NPC: same deps already match
      if (precompute_coeff_) {
        key = ss.generate_id(vector_init(vector_init(e->dep_)));
      }

      // Loop over all matching target0 entries
      for (const auto& ee : targets0[key]) {
        // Mark all matches for substitution
        subs_from.push_back(ee);
        // Substitute with self
        subs_to.push_back(e);
      }
    }

    // Loop over first order targets, targets1
    for (const auto& ee : targets1) {
      for (const auto& e : ee.second) {
        // Compute key that matches targets0
        // Precompute coeff: strip last arg (dC) to match targets0's (x, C)
        // NPC/grid: all levels share the same deps, use directly
        key = precompute_coeff_ ?
          ss.generate_id(vector_init(e->dep_)) :
          ss.generate_id(e->dep_);

        // Loop over all matching target0 entries
        for (const auto& ee : targets0[key]) {
          // Mark all matches for substitution
          subs_from.push_back(ee);
          // Substitute with self
          subs_to.push_back(e);
        }
      }
    }

  }


} // namespace casadi
