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

#include "shf_interpolant.hpp"
#include <limits>

namespace casadi {

  extern "C"
  int CASADI_INTERPOLANT_SHF_EXPORT
  casadi_register_interpolant_shf(Interpolant::Plugin* plugin) {
    plugin->creator = SHFInterpolant::creator;
    plugin->name = "shf";
    plugin->doc = SHFInterpolant::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &SHFInterpolant::options_;
    plugin->deserialize = &SHFInterpolant::deserialize;
    plugin->exposed.do_inline = nullptr;
    return 0;
  }

  extern "C"
  void CASADI_INTERPOLANT_SHF_EXPORT casadi_load_interpolant_shf() {
    Interpolant::registerPlugin(casadi_register_interpolant_shf);
  }

  const Options SHFInterpolant::options_
  = {{&Interpolant::options_},
     {{"smoothness_order",
       {OT_INT,
        "Smoothness order k: the approximation is C^k. Default: 3."}},
      {"epsilon",
       {OT_DOUBLEVECTOR,
        "Absolute half-width of the smoothing balls around the interior grid points, "
        "one value shared by all axes or one per axis. "
        "Default: 0.1 times the smallest grid spacing of each axis."}},
      {"epsilon_parametric",
       {OT_BOOL,
        "Take epsilon as an additional (non-differentiable) input 'eps' of length ndim "
        "instead of fixing it at construction. Default: false."}}
     }
  };

  SHFInterpolant::~SHFInterpolant() {
    clear_mem();
  }

  SHFInterpolant::
  SHFInterpolant(const std::string& name,
                 const std::vector<double>& grid,
                 const std::vector<casadi_int>& offset,
                 const std::vector<double>& values,
                 casadi_int m)
                 : Interpolant(name, grid, offset, values, m) {
    epsilon_parametric_ = false;
  }

  void SHFInterpolant::init(const Dict& opts) {
    smoothness_order_ = 3;
    epsilon_.clear();
    epsilon_parametric_ = false;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="smoothness_order") {
        smoothness_order_ = op.second;
      } else if (op.first=="epsilon") {
        epsilon_ = op.second;
      } else if (op.first=="epsilon_parametric") {
        epsilon_parametric_ = op.second;
      }
    }

    // Call the base class initializer
    Interpolant::init(opts);

    casadi_assert(!has_parametric_grid(), "Parametric grid not supported");

    std::vector< std::vector<double> > grid;
    for (casadi_int k=0;k<ndim_;++k) {
      grid.push_back(std::vector<double>(grid_.begin()+offset_[k], grid_.begin()+offset_[k+1]));
    }

    // Smallest spacing of each axis: the epsilon-balls around neighbouring grid
    // points must stay disjoint, otherwise the partition of unity -- and with it
    // the value -- is silently wrong
    std::vector<double> min_h(ndim_, std::numeric_limits<double>::infinity());
    for (casadi_int k=0;k<ndim_;++k) {
      for (size_t i=0;i+1<grid[k].size();++i) {
        min_h[k] = std::min(min_h[k], grid[k][i+1]-grid[k][i]);
      }
    }

    if (epsilon_parametric_) {
      casadi_assert(epsilon_.empty(),
        "Options 'epsilon' and 'epsilon_parametric' are mutually exclusive.");
    } else if (epsilon_.empty()) {
      // A tenth of the smallest spacing of each axis
      for (casadi_int k=0;k<ndim_;++k) epsilon_.push_back(0.1*min_h[k]);
    } else {
      casadi_assert(epsilon_.size()==1 || epsilon_.size()==static_cast<size_t>(ndim_),
        "Option 'epsilon' must hold one value or one per axis (" + str(ndim_) + "), got "
        + str(epsilon_.size()) + ".");
      for (casadi_int k=0;k<ndim_;++k) {
        double e = epsilon_[epsilon_.size()==1 ? 0 : k];
        casadi_assert(e>=0, "epsilon must be non-negative, got " + str(e) + ".");
        casadi_assert(e<0.5*min_h[k], "epsilon (" + str(e) + ") must be smaller than "
          "half the smallest grid spacing of axis " + str(k) + " (" + str(0.5*min_h[k])
          + ") so that the epsilon-neighbourhoods of neighbouring grid points stay "
          "disjoint.");
      }
    }

    Dict opts_shf;
    opts_shf["lookup_mode"] = lookup_modes_;
    Function F = has_parametric_values() ?
      shf_spline(name_ + "_shf", grid, smoothness_order_, m_, opts_shf) :
      shf_spline(name_ + "_shf", grid, values_, smoothness_order_, m_, opts_shf);

    MX x = MX::sym("x", ndim_, batch_x_);
    std::vector<MX> args = {x};
    std::vector<std::string> names = {"x"};
    MX coeff;
    if (has_parametric_values()) {
      coeff = MX::sym("c", coeff_size());
      args.push_back(coeff);
      names.push_back("c");
    }
    MX eps;
    if (epsilon_parametric_) {
      eps = MX::sym("eps", ndim_);
      args.push_back(eps);
      names.push_back("eps");
    } else {
      eps = DM(epsilon_);
    }

    // One call per column of x
    std::vector<MX> cols = horzsplit(x);
    for (auto& c : cols) {
      std::vector<MX> a = {c};
      if (has_parametric_values()) a.push_back(coeff);
      a.push_back(eps);
      c = F(a)[0];
    }

    // Only x is differentiable
    std::vector<bool> is_diff(args.size(), false);
    is_diff[0] = true;
    Dict opts_wrapper;
    opts_wrapper["is_diff_in"] = is_diff;
    opts_wrapper["is_diff_in"] = is_diff;
    S_ = Function("wrapper", args, {horzcat(cols)}, names, {"f"}, opts_wrapper);

    alloc_w(S_.sz_w());
    alloc_iw(S_.sz_iw());
    alloc_arg(S_.sz_arg());
    alloc_res(S_.sz_res());
  }

  Sparsity SHFInterpolant::get_sparsity_in(casadi_int i) {
    if (epsilon_parametric_ && i==arg_epsilon()) return Sparsity::dense(ndim_);
    return Interpolant::get_sparsity_in(i);
  }

  std::string SHFInterpolant::get_name_in(casadi_int i) {
    if (epsilon_parametric_ && i==arg_epsilon()) return "eps";
    return Interpolant::get_name_in(i);
  }

  void SHFInterpolant::find(
    std::map<FunctionInternal*, std::pair<Function, size_t > >& all_fun,
      casadi_int max_depth) const {
    // Call to base class
    FunctionInternal::find(all_fun, max_depth);
    add_embedded(all_fun, S_, max_depth);
  }

  int SHFInterpolant::eval(const double** arg, double** res,
                           casadi_int* iw, double* w, void* mem) const {
    setup(mem, arg, res, iw, w);
    scoped_checkout<Function> m(S_);
    return S_(arg, res, iw, w, m);
  }

  void SHFInterpolant::eval_mx(const MXVector& arg, MXVector& res,
      bool always_inline, bool never_inline) const {
    // An opaque call node would hide the seeds from the AD of the enclosing graph:
    // second derivatives through it evaluate J*(seed of seed), a structural zero
    // in a Hessian, at the cost of a full extra kernel call.
    if (never_inline || never_inline_) {
      FunctionInternal::eval_mx(arg, res, always_inline, never_inline);
    } else {
      S_->eval_mx(arg, res, true, false);
    }
  }

  void SHFInterpolant::codegen_body(CodeGenerator& g) const {
    S_->codegen_body(g);
  }

  void SHFInterpolant::codegen_declarations(CodeGenerator& g) const {
    S_->codegen_declarations(g);
  }

  Function SHFInterpolant::
  get_jacobian(const std::string& name,
               const std::vector<std::string>& inames,
               const std::vector<std::string>& onames,
               const Dict& opts) const {
    return S_->get_jacobian(name, inames, onames, opts);
  }

  Function SHFInterpolant::
  get_forward(casadi_int nfwd, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    return S_->get_forward(nfwd, name, inames, onames, opts);
  }

  Function SHFInterpolant::
  get_reverse(casadi_int nadj, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    return S_->get_reverse(nadj, name, inames, onames, opts);
  }

  SHFInterpolant::SHFInterpolant(DeserializingStream& s) : Interpolant(s) {
    s.version("SHFInterpolant", 1);
    s.unpack("SHFInterpolant::s", S_);
    s.unpack("SHFInterpolant::epsilon_parametric", epsilon_parametric_);
  }

  void SHFInterpolant::serialize_body(SerializingStream &s) const {
    Interpolant::serialize_body(s);
    s.version("SHFInterpolant", 1);
    s.pack("SHFInterpolant::s", S_);
    s.pack("SHFInterpolant::epsilon_parametric", epsilon_parametric_);
  }

} // namespace casadi
