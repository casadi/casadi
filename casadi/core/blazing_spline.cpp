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
    } else if (i==3) {
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
      return Sparsity::dense(1, knots_.size());
    } else if (i==2) {
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

    init_derived_members();

    casadi_assert(knots.size()>=1, "blazing_spline only defined for 1D/2D/3D");
    casadi_assert(knots.size()<=3, "blazing_spline only defined for 1D/2D/3D");
  }

  void BlazingSplineFunction::init_derived_members() {
    // Compute coefficient tensor size
    nc_ = 1;
    for (const auto& e : knots_) {
      nc_ *= e.size()-4;
    }

    // Compute coefficient tensor size
    ndc_ = 0;
    for (casadi_int k=0;k<knots_.size();++k) {
      casadi_int ndc = 1;
      for (casadi_int i=0;i<knots_.size();++i) {
        ndc *= knots_[i].size() - 4 - (i==k);
      }
      ndc_+= ndc;
    }

    nddc_ = 0;
    for (casadi_int k=0;k<knots_.size();++k) {
      for (casadi_int kk=0;kk<knots_.size();++kk) {
        casadi_int ndc = 1;
        for (casadi_int i=0;i<knots_.size();++i) {
          ndc *= knots_[i].size() - 4 - (i==k)-(i==kk);
        }
        // We only need the triangular part
        if (kk>=k) {
          nddc_+= ndc;
        }
      }
    }


    std::vector<casadi_int> offset;
    std::vector<double> stacked;
    Interpolant::stack_grid(knots_, knots_offset_, knots_stacked_);
  }

  const Options BlazingSplineFunction::options_
  = {{&FunctionInternal::options_},
     {
     }
  };

  void BlazingSplineFunction::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
    }

    casadi_int n_dims = knots_.size();

    // Arrays for holding inputs and outputs
    // TODO
    alloc_iw(4*n_dims+2+1000);
    alloc_w(4*(n_dims-1)+2*3+1+n_dims+1+1000);
  }

  BlazingSplineFunction::~BlazingSplineFunction() {
    clear_mem();
  }

  void BlazingSplineFunction::codegen_body(CodeGenerator& g) const {
    if (knots_.size() == 3) {
      g.add_auxiliary(CodeGenerator::AUX_BLAZING_3D_BOOR_EVAL);
    }
    if (knots_.size() == 2) {
      g.add_auxiliary(CodeGenerator::AUX_BLAZING_2D_BOOR_EVAL);
    }
    if (knots_.size() == 1) {
      g.add_auxiliary(CodeGenerator::AUX_BLAZING_1D_BOOR_EVAL);
    }
    g.add_include("simde/x86/avx2.h");
    g.add_include("simde/x86/fma.h");

    std::string knots_offset = g.constant(knots_offset_);
    std::string knots_stacked = g.constant(knots_stacked_);

    std::vector<std::string> lookup_mode;
    std::vector<casadi_int> degree(3, knots_.size());
    std::vector<casadi_int> mode = Interpolant::interpret_lookup_mode(lookup_mode, knots_stacked_, knots_offset_, degree, degree);

    std::string fun_name = "casadi_blazing_" + str(knots_.size()) + "d_boor_eval";

    if (diff_order_==0) {
      // Codegen function body
      g << fun_name + "(res[0], 0, 0, " +
            knots_stacked + ", " +
            knots_offset + ", " +
            "arg[1], 0, 0, " + 
            "arg[0], " + 
            g.constant(mode) + ", " +
            "iw, w);\n";
    } else if (diff_order_==1) {
      // Codegen function body
      g << fun_name + "(res[0], res[1], 0, " +
            knots_stacked + ", " +
            knots_offset + ", " +
            "arg[1], arg[2], 0, " + 
            "arg[0], " + 
            g.constant(mode) + ", " +
            "iw, w);\n";
    } else if (diff_order_==2) {
      // Codegen function body
      g << fun_name + "(res[0], res[1], res[2], " +
            knots_stacked + ", " +
            knots_offset + ", " +
            "arg[1], arg[2], arg[3], " + 
            "arg[0], " + 
            g.constant(mode) + ", " +
            "iw, w);\n";
    } else {
      casadi_assert_dev(false);
    }
  }

  bool BlazingSplineFunction::has_jacobian() const {
    return diff_order_<2;
  }

  Function BlazingSplineFunction::get_jacobian(const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    size_t N = knots_.size();
    MX C = MX::sym("C", nc_);
    MX x = MX::sym("x", N);

    std::vector<casadi_int> coeffs_dims(N+1);
    coeffs_dims[0] = 1;
    for (casadi_int i=0; i<N; ++i) {
      coeffs_dims[i+1] = knots_[i].size()-4;
    }

    std::vector<casadi_int> degree(N, 3);
    std::vector< std::vector< std::vector<double> > > knots_d(N);
    std::vector< std::vector< casadi_int> > degree_d(N);

    std::vector<MX> dCv;
    for (size_t i=0;i<N;++i) {
      dCv.push_back(BSplineCommon::derivative_coeff(i, knots_stacked_, knots_offset_, degree, coeffs_dims, C, knots_d[i], degree_d[i]));
    }

    MX dC = vertcat(dCv);

    Dict Jopts = combine(jacobian_options_, der_options_);
    Jopts = combine(opts, Jopts);
    Jopts = combine(Jopts, generate_options("jacobian"));
    Jopts["derivative_of"] = self();

    std::string fJname = name_ + "_der";

    Function fJ;
    if (!incache(fJname, fJ)) {
      fJ = Function::create(new BlazingSplineFunction(fJname, knots_, diff_order_+1), Jopts);
      // Save in cache
      tocache(fJ);
    }

    if (diff_order_==0) {
      std::vector<MX> ret = fJ(std::vector<MX>{x, C, dC});

      Function jac(name,
        {x, C, MX(1, 1)},
        {ret[1], MX(1, nc_)},
        inames,
        onames,
        {{"always_inline", true}});

      return jac;
    } else {
      MX ddC;
      if (N==3) {
        int size0 = knots_[0].size();
        int size1 = knots_[1].size();
        int size2 = knots_[2].size();
        std::vector<casadi_int> offset;
        std::vector<double> stacked;
        std::vector< std::vector<double> > knots_dummy;
        std::vector< casadi_int> degree_dummy;
        Interpolant::stack_grid(knots_d[0], offset, stacked);
        MX ddC00 = BSplineCommon::derivative_coeff(0, stacked, offset, degree_d[0],
          {1, size0-5, size1-4, size2-4}, dCv[0], knots_dummy, degree_dummy);
        Interpolant::stack_grid(knots_d[1], offset, stacked);
        MX ddC11 = BSplineCommon::derivative_coeff(1, stacked, offset, degree_d[1],
          {1, size0-4, size1-5, size2-4}, dCv[1], knots_dummy, degree_dummy);
        Interpolant::stack_grid(knots_d[2], offset, stacked);
        MX ddC22 = BSplineCommon::derivative_coeff(2, stacked, offset, degree_d[2],
          {1, size0-4, size1-4, size2-5}, dCv[2], knots_dummy, degree_dummy);
        Interpolant::stack_grid(knots_d[0], offset, stacked);
        MX ddC01 = BSplineCommon::derivative_coeff(1, stacked, offset, degree_d[0],
          {1, size0-5, size1-4, size2-4}, dCv[0], knots_dummy, degree_dummy);
        Interpolant::stack_grid(knots_d[1], offset, stacked);
        MX ddC12 = BSplineCommon::derivative_coeff(2, stacked, offset, degree_d[1],
          {1, size0-4, size1-5, size2-4}, dCv[1], knots_dummy, degree_dummy);
        Interpolant::stack_grid(knots_d[2], offset, stacked);
        MX ddC20 = BSplineCommon::derivative_coeff(0, stacked, offset, degree_d[2],
          {1, size0-4, size1-4, size2-5}, dCv[2], knots_dummy, degree_dummy);
        ddC = vertcat(ddC00,ddC11,ddC22,ddC01,ddC12,ddC20);
      } else if (N==2) {
        int size0 = knots_[0].size();
        int size1 = knots_[1].size();
        std::vector<casadi_int> offset;
        std::vector<double> stacked;
        std::vector< std::vector<double> > knots_dummy;
        std::vector< casadi_int> degree_dummy;
        Interpolant::stack_grid(knots_d[0], offset, stacked);
        MX ddC00 = BSplineCommon::derivative_coeff(0, stacked, offset, degree_d[0],
          {1, size0-5, size1-4}, dCv[0], knots_dummy, degree_dummy);
        Interpolant::stack_grid(knots_d[1], offset, stacked);
        MX ddC11 = BSplineCommon::derivative_coeff(1, stacked, offset, degree_d[1],
          {1, size0-4, size1-5}, dCv[1], knots_dummy, degree_dummy);
        Interpolant::stack_grid(knots_d[0], offset, stacked);
        MX ddC01 = BSplineCommon::derivative_coeff(1, stacked, offset, degree_d[0],
          {1, size0-5, size1-4}, dCv[0], knots_dummy, degree_dummy);
        ddC = vertcat(ddC00,ddC11,ddC01);
      } else if (N==1) {
        int size0 = knots_[0].size();
        std::vector<casadi_int> offset;
        std::vector<double> stacked;
        std::vector< std::vector<double> > knots_dummy;
        std::vector< casadi_int> degree_dummy;
        Interpolant::stack_grid(knots_d[0], offset, stacked);
        ddC = BSplineCommon::derivative_coeff(0, stacked, offset, degree_d[0],
          {1, size0-5}, dCv[0], knots_dummy, degree_dummy);
      }

      std::vector<MX> ret = fJ(std::vector<MX>{x, C, dC, ddC});

      Function jac(name,
        {x, C, MX(1, ndc_), MX(1, 1), MX(1, N)},
        {ret[1], MX(1, nc_), MX(1, ndc_), ret[2], MX(N, nc_), MX(N, ndc_)},
        inames,
        onames,
        {{"always_inline", true}});

      return jac;
    }
  }

  void BlazingSplineFunction::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);

    s.version("BlazingSplineFunction", 1);
    s.pack("BlazingSplineFunction::diff_order", diff_order_);
    s.pack("BlazingSplineFunction::knots", knots_);
  }

  BlazingSplineFunction::BlazingSplineFunction(DeserializingStream & s) : FunctionInternal(s) {
    s.version("BlazingSplineFunction", 1);
    s.unpack("BlazingSplineFunction::diff_order", diff_order_);
    s.unpack("BlazingSplineFunction::knots", knots_);
    init_derived_members();
  }

  ProtoFunction* BlazingSplineFunction::deserialize(DeserializingStream& s) {
    return new BlazingSplineFunction(s);
  }


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

    StringSerializer ss;
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
          } else if (fun->derivative_of_==base) {
            targets1[key].push_back(e);
          } else if (fun->derivative_of_->derivative_of_==base) {
            targets2.push_back(e);
          }
        }
      }
    }

    // Loop over second order targets, targets2
    for (const auto& e : targets2) {

      // Compute key of all but last args
      key = ss.generate_id(vector_init(e->dep_));

      // Loop over all matching target1 entries
      for (const auto& ee : targets1[key]) {
        // Mark all matches for substitution
        subs_from.push_back(ee);
        // Substitute with self
        subs_to.push_back(e);
      }

      // Compute key of all but last args
      key = ss.generate_id(vector_init(vector_init(e->dep_)), 0);

      // Loop over all matching target1 entries
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
        // Compute key of all but last args
        key = ss.generate_id(vector_init(e->dep_), 0);

        // Loop over all matching target1 entries
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
