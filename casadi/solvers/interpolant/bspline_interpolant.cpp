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
    return 0;
  }

  extern "C"
  void CASADI_INTERPOLANT_BSPLINE_EXPORT casadi_load_interpolant_bspline() {
    Interpolant::registerPlugin(casadi_register_interpolant_bspline);
  }

  Options BSplineInterpolant::options_
  = {{&Interpolant::options_},
      {{"degree",
       {OT_INTVECTOR,
        "Sets, for each grid dimenion, the degree of the spline."}},
       {"linear_solver",
        {OT_STRING,
         "Solver used for constructing the coefficient tensor."}}
     }
  };

  BSplineInterpolant::
  BSplineInterpolant(const string& name,
                    const std::vector<double>& grid,
                    const std::vector<int>& offset,
                    const vector<double>& values)
                    : Interpolant(name, grid, offset, values) {

  }

  BSplineInterpolant::~BSplineInterpolant() {
  }

  std::vector<double> meshgrid(const std::vector< std::vector<double> >& grid) {
    std::vector<int> cnts(grid.size()+1, 0);
    std::vector<int> sizes(grid.size(), 0);
    for (int k=0;k<grid.size();++k) sizes[k]= grid[k].size();

    int total_iter = 1;
    for (int k=0;k<grid.size();++k) total_iter*= sizes[k];

    int n_dims = grid.size();

    std::vector<double> ret(total_iter*n_dims);
    for (int i=0;i<total_iter;++i) {

      for (int j=0;j<grid.size();++j) {
        ret[i*n_dims+j] = grid[j][cnts[j]];
      }

      cnts[0]++;
      int j = 0;
      while (j<n_dims && cnts[j]==sizes[j]) {
        cnts[j] = 0;
        j++;
        cnts[j]++;
      }

    }

    return ret;
  }

  std::vector<double> not_a_knot(const std::vector<double>& x, int k) {
    std::vector<double> ret;
    if (k%2) {
      int m = (k-1)/2;
      for (int i=0;i<k+1;++i) ret.push_back(x[0]);
      for (int i=0;i<x.size()-2*m-2;++i) ret.push_back(x[m+1+i]);
      for (int i=0;i<k+1;++i) ret.push_back(x[x.size()-1]);
    } else {
      casadi_error("Not implemented");
      //for (int i=0;i<k+1;++i) ret.push_back(x[0]);
      //for (int i=0;i<x.size()-2*m-2;++i) ret.push_back(x[m+1+i]);
      //for (int i=0;i<k+1;++i) ret.push_back(x[x.size()-1]);
    }
    return ret;
  }

  void BSplineInterpolant::init(const Dict& opts) {
    // Call the base class initializer
    Interpolant::init(opts);

    degree_  = std::vector<int>(offset_.size()-1, 3);

    linear_solver_ = "lsqr";

    // Read options
    for (auto&& op : opts) {
      if (op.first=="degree") {
        degree_ = op.second;
      } else if (op.first=="linear_solver") {
        linear_solver_ = op.second.to_string();
      }
    }

    casadi_assert(degree_.size()==offset_.size()-1);

    std::vector< std::vector<double> > knots;
    std::vector< std::vector<double> > grid;
    for (int k=0;k<degree_.size();++k) {
      std::vector<double> local_grid(grid_.begin()+offset_[k], grid_.begin()+offset_[k+1]);
      grid.push_back(local_grid);
      knots.push_back(not_a_knot(local_grid, degree_[k]));
    }

    Dict opts_dual;
    opts_dual["ad_weight_sp"] = 0;

    Function B = Function::bspline_dual("spline", knots, meshgrid(grid), degree_, 1, false,
      opts_dual);

    Function Jf = B.jacobian_old(0, 0);

    MX C = MX::sym("C", B.size_in(0));

    MX Js = Jf(std::vector<MX>{C})[0];
    Function temp = Function("J", {C}, {Js});
    DM J = temp(std::vector<DM>{0})[0];

    casadi_assert(J.size1()==J.size2());

    DM C_opt = solve(J, DM(values_), linear_solver_);

    double fit = static_cast<double>(norm_1(mtimes(J, C_opt) - DM(values_)));

    userOut() << "Lookup table fitting error: " << fit << std::endl;

    S_ = Function::bspline("spline", knots, C_opt.nonzeros(), degree_, 1);

    alloc_w(S_->sz_w(), true);
    alloc_iw(S_->sz_iw(), true);

  }

  void BSplineInterpolant::eval(void* mem, const double** arg, double** res,
                               int* iw, double* w) const {
    S_->eval(mem, arg, res, iw, w);
  }

  void BSplineInterpolant::generateBody(CodeGenerator& g) const {
    S_->generateBody(g);
  }

  Function BSplineInterpolant::
  get_jacobian2(const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {
    return S_->get_jacobian2(name, inames, onames, opts);
  }

} // namespace casadi
