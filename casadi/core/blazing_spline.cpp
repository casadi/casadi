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

#include <fstream>
#include <iostream>
#include <sstream>

namespace casadi {

  Function blazing_spline(const std::string& name,
      const std::vector< std::vector<double> >& knots,
      const Dict& opts) {
    return Function::create(new BlazingSplineFunction(name, knots, 0), opts);
  }

  size_t BlazingSplineFunction::get_n_in() {
    return 2+diff_order_;
  }
  size_t BlazingSplineFunction::get_n_out() {
    return diff_order_+1;
  }

  bool BlazingSplineFunction::get_diff_in(casadi_int i) {
    return i==0;
  }

  Sparsity BlazingSplineFunction::get_sparsity_in(casadi_int i) {
    if (i==0) {
      return Sparsity::dense(knots_.size());
    } else if (i==1) {
      return Sparsity::dense(nc_);
    } else if (i==2) {
      return Sparsity::dense(ndc_);
    } else {
      casadi_assert_dev(false);
      return Sparsity();
    }
  }
  Sparsity BlazingSplineFunction::get_sparsity_out(casadi_int i) {
    if (i==0) {
      return Sparsity::dense(1, 1);
    } else if (i==1) {
      return Sparsity::dense(1, knots_.size());
    } else if (i=2) {
      return Sparsity::dense(knots_.size(), knots_.size());
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
    } else if (i==2) {
      return "dC";
    } else if (i==3) {
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
    casadi_int diff_order) : FunctionInternal(name), knots_(knots), diff_order_(diff_order) {
    
    // Compute coefficient tensor size
    nc_ = 1;
    for (const auto& e : knots) {
      nc_ *= e.size()-4;
    }

    // Compute coefficient tensor size
    ndc_ = 0;
    for (casadi_int k=0;k<knots.size();++k) {
      casadi_int ndc = 1;
      for (casadi_int i=0;i<knots.size();++i) {
        ndc *= knots[i].size() - 4 - (i==k);
      }
      ndc_+= ndc;
    }


    std::vector<casadi_int> offset;
    std::vector<double> stacked;
    Interpolant::stack_grid(knots, knots_offset_, knots_stacked_);

  }

  const Options BlazingSplineFunction::options_
  = {{&FunctionInternal::options_},
     {{"buffered",
      {OT_BOOL,
        "Buffer the calls, user does not need to "}},
       {"jac",
      {OT_STRING,
        "Function body for Jacobian"}},
      {"hess",
       {OT_STRING,
        "Function body for Hessian"}}
     }
  };

  void BlazingSplineFunction::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="buffered") {
        buffered_ = op.second;
      } else if (op.first=="jac") {
        jac_body_ = op.second.to_string();
      } else if (op.first=="hess") {
        hess_body_ = op.second.to_string();
      }
    }

    casadi_int n_dims = knots_.size();

    // Arrays for holding inputs and outputs
    alloc_iw(4*n_dims+2+1000);
    alloc_w(4*(n_dims-1)+2*3+1+n_dims+1+1000);
  }

  BlazingSplineFunction::~BlazingSplineFunction() {
    clear_mem();
  }

  void BlazingSplineFunction::codegen_body(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_BLAZING_SPLINE);
    g.add_include("simde/x86/avx2.h");
    g.add_include("simde/x86/fma.h");

    std::string knots_offset = g.constant(knots_offset_);
    std::string knots_stacked = g.constant(knots_stacked_);

    std::vector<std::string> lookup_mode;
    std::vector<casadi_int> degree(3, knots_.size());
    std::vector<casadi_int> mode = Interpolant::interpret_lookup_mode(lookup_mode, knots_stacked_, knots_offset_, degree, degree);


    if (diff_order_==0) {
      // Codegen function body
      g << "casadi_blazing_nd_boor_eval_dim3(res[0], 0, " +
            knots_stacked + ", " +
            knots_offset + ", " +
            "arg[1], 0, " + 
            "arg[0], " + 
            g.constant(mode) + ", " +
            "iw, w);\n";
    } else if (diff_order_==1) {
      // Codegen function body
      g << "casadi_blazing_nd_boor_eval_dim3(res[0], res[1], " +
            knots_stacked + ", " +
            knots_offset + ", " +
            "arg[1], arg[2], " + 
            "arg[0], " + 
            g.constant(mode) + ", " +
            "iw, w);\n";
    }
  }

  bool BlazingSplineFunction::has_jacobian() const {
    return true;
  }

  Function BlazingSplineFunction::get_jacobian(const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {


    MX C = MX::sym("C", nc_);
    MX x = MX::sym("x", knots_.size());

    // Compute coefficient tensor size
    casadi_int ndc_ = 1;
    for (const auto& e : knots_) {
      ndc_ *= e.size()-5;
    }

    std::vector<casadi_int> coeffs_dims(knots_.size()+1);
    coeffs_dims[0] = 1;
    for (casadi_int i=0; i<knots_.size(); ++i) {
      coeffs_dims[i+1] = knots_[i].size()-4;
    }

    std::vector<casadi_int> offset;
    std::vector<double> stacked;
    Interpolant::stack_grid(knots_, offset, stacked);

    std::vector<casadi_int> degree(knots_.size(), 3);
    MX dC0 = BSplineCommon::derivative_coeff(0, stacked, offset, degree, coeffs_dims, C);
    MX dC1 = BSplineCommon::derivative_coeff(1, stacked, offset, degree, coeffs_dims, C);
    MX dC2 = BSplineCommon::derivative_coeff(2, stacked, offset, degree, coeffs_dims, C);
    MX dC = vertcat(dC0,dC1,dC2);

    Dict Jopts = combine(jacobian_options_, der_options_);
    Jopts = combine(opts, Jopts);
    Jopts = combine(Jopts, generate_options("jacobian"));

    Function fJ = Function::create(new BlazingSplineFunction(name, knots_, diff_order_+1), Jopts);

    MX out_f = MX(1, 1);

    std::vector<MX> ret = fJ(std::vector<MX>{x, C, dC});

    Function jac(name,
      {x, C, out_f},
      {ret[1], MX(1, nc_)},
      inames,
      onames,
      {{"always_inline", true}});

    return jac;
  }

} // namespace casadi
