/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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


#include "bspline_interpolant.hpp"
#include "casadi/core/bspline.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_INTERPOLANT_BSPLINE_EXPORT
  casadi_register_interpolant_bspline(Interpolant::Plugin* plugin) {
    plugin->creator = BSplineInterpolant::creator;
    plugin->name = "bspline";
    plugin->doc = BSplineInterpolant::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &BSplineInterpolant::options_;
    plugin->deserialize = &BSplineInterpolant::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_INTERPOLANT_BSPLINE_EXPORT casadi_load_interpolant_bspline() {
    Interpolant::registerPlugin(casadi_register_interpolant_bspline);
  }

  const Options BSplineInterpolant::options_
  = {{&Interpolant::options_},
      {{"degree",
       {OT_INTVECTOR,
        "Sets, for each grid dimension, the degree of the spline."}},
       {"linear_solver",
        {OT_STRING,
         "Solver used for constructing the coefficient tensor."}},
       {"algorithm",
        {OT_STRING,
         "Algorithm used for fitting the data: 'not_a_knot' (default, same as Matlab),"
        " 'smooth_linear'."}},
       {"smooth_linear_frac",
        {OT_DOUBLE,
         "When 'smooth_linear' algorithm is active, determines sharpness between"
         " 0 (sharp, as linear interpolation) and 0.5 (smooth)."
         "Default value is 0.1."}}
     }
  };

  BSplineInterpolant::~BSplineInterpolant() {
    clear_mem();
  }

  BSplineInterpolant::
  BSplineInterpolant(const string& name,
                    const std::vector<double>& grid,
                    const std::vector<casadi_int>& offset,
                    const vector<double>& values,
                    casadi_int m)
                    : Interpolant(name, grid, offset, values, m) {

  }

  std::vector<double> BSplineInterpolant::not_a_knot(const std::vector<double>& x, casadi_int k) {
    std::vector<double> ret;
    if (k%2) {
      casadi_int m = (k-1)/2;
      casadi_assert(x.size()>=2*m+2, "Need more data points");
      for (casadi_int i=0;i<k+1;++i) ret.push_back(x[0]);
      for (casadi_int i=0;i<x.size()-2*m-2;++i) ret.push_back(x[m+1+i]);
      for (casadi_int i=0;i<k+1;++i) ret.push_back(x[x.size()-1]);
    } else {
      casadi_error("Not implemented");
      //for (casadi_int i=0;i<k+1;++i) ret.push_back(x[0]);
      //for (casadi_int i=0;i<x.size()-2*m-2;++i) ret.push_back(x[m+1+i]);
      //for (casadi_int i=0;i<k+1;++i) ret.push_back(x[x.size()-1]);
    }
    return ret;
  }

  void BSplineInterpolant::init(const Dict& opts) {

    degree_  = std::vector<casadi_int>(offset_.size()-1, 3);

    linear_solver_ = "lsqr";
    algorithm_ = ALG_NOT_A_KNOT;
    smooth_linear_frac_ = 0.1;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="degree") {
        degree_ = op.second;
      } else if (op.first=="linear_solver") {
        linear_solver_ = op.second.to_string();
      } else if (op.first=="algorithm") {
        std::string alg = op.second.to_string();
        if (alg=="not_a_knot") {
          algorithm_ = ALG_NOT_A_KNOT;
        } else if (alg=="smooth_linear") {
          algorithm_ = ALG_SMOOTH_LINEAR;
        } else {
          casadi_error("Algorithm option invalid: " + get_options().info("algorithm"));
        }
      } else if (op.first=="smooth_linear_frac") {
        smooth_linear_frac_ = op.second;
        casadi_assert(smooth_linear_frac_>0 && smooth_linear_frac_<0.5,
          "smooth_linear_frac must be in ]0,0.5[");
      }
    }

    casadi_assert_dev(degree_.size()==offset_.size()-1);

    // Call the base class initializer
    Interpolant::init(opts);

    casadi_assert(!has_parametric_grid(), "Parametric grid not supported");

    MX x = MX::sym("x", ndim_);

    if (has_parametric_values()) {
      MX coeff = MX::sym("coeff", coeff_size());

      MX e = construct_graph(x, coeff, opts);

      S_ = Function("wrapper", {x, coeff}, {e});
    } else {
      MX e = construct_graph(x, DM(values_), opts);
      S_ = Function("wrapper", {x}, {e});
    }

    alloc_w(S_.sz_w());
    alloc_iw(S_.sz_iw());
    alloc_arg(S_.sz_arg());
    alloc_res(S_.sz_res());
  }

  std::vector<double> BSplineInterpolant::greville_points(const std::vector<double>& x,
                                                          casadi_int deg) {
    casadi_int dim = x.size()-deg-1;
    std::vector<double> ret(dim);
    for (casadi_int i = 0; i < dim; ++i) {
      ret[i] = 0;
      for (casadi_int j = 0; j < deg; j++) {
        ret[i] += x[i+1+j];
      }
      ret[i] = ret[i] / deg;
    }
    return ret;
  }

  int BSplineInterpolant::eval(const double** arg, double** res,
                                casadi_int* iw, double* w, void* mem) const {
    scoped_checkout<Function> m(S_);
    return S_(arg, res, iw, w, m);
  }

  void BSplineInterpolant::codegen_body(CodeGenerator& g) const {
    S_->codegen_body(g);
  }

  void BSplineInterpolant::codegen_declarations(CodeGenerator& g) const {
    S_->codegen_declarations(g);
  }

  Function BSplineInterpolant::
  get_jacobian(const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {
    return S_->get_jacobian(name, inames, onames, opts);
  }

  BSplineInterpolant::BSplineInterpolant(DeserializingStream& s) : Interpolant(s) {
    s.version("BSplineInterpolant", 1);
    s.unpack("BSplineInterpolant::s", S_);
  }

  void BSplineInterpolant::serialize_body(SerializingStream &s) const {
    Interpolant::serialize_body(s);

    s.version("BSplineInterpolant", 1);
    s.pack("BSplineInterpolant::s", S_);
  }

} // namespace casadi
