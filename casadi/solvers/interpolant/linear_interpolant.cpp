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


#include "linear_interpolant.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_INTERPOLANT_LINEAR_EXPORT
  casadi_register_interpolant_linear(Interpolant::Plugin* plugin) {
    plugin->creator = LinearInterpolant::creator;
    plugin->name = "linear";
    plugin->doc = LinearInterpolant::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &LinearInterpolant::options_;
    return 0;
  }

  extern "C"
  void CASADI_INTERPOLANT_LINEAR_EXPORT casadi_load_interpolant_linear() {
    Interpolant::registerPlugin(casadi_register_interpolant_linear);
  }

  Options LinearInterpolant::options_
  = {{&Interpolant::options_},
     {{"lookup_mode",
       {OT_STRINGVECTOR,
        "Sets, for each grid dimenion, the lookup algorithm used to find the correct index. "
        "'linear' uses a for-loop + break; "
        "'exact' uses floored division (only for uniform grids)."}}
     }
  };

  LinearInterpolant::
  LinearInterpolant(const string& name,
                    const std::vector<double>& grid,
                    const std::vector<int>& offset,
                    const vector<double>& values)
                    : Interpolant(name, grid, offset, values) {
  }

  LinearInterpolant::~LinearInterpolant() {
  }

  void LinearInterpolant::init(const Dict& opts) {
    // Call the base class initializer
    Interpolant::init(opts);

    lookup_mode_ = std::vector<int>(offset_.size()-1, 0);

    std::vector<std::string> lookup_mode;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="lookup_mode") {
        lookup_mode = op.second;
      }
    }

    if (!lookup_mode.empty()) {
      casadi_assert(lookup_mode.size()==offset_.size()-1);
      for (int i=0;i<offset_.size()-1;++i) {
        if (lookup_mode[i]=="linear") {
          lookup_mode_[i] = 0;
        } else if (lookup_mode[i]=="exact") {
          lookup_mode_[i] = 1;
          std::vector<double> grid(
              grid_.begin()+offset_[i],
              grid_.begin()+offset_[i+1]);
          casadi_assert(isIncreasing(grid) && isEquallySpaced(grid));
        } else {
          casadi_error("Unknown lookup_mode option '" + lookup_mode[i] + ". "
                       "Allowed values: linear, exact.");
        }
      }
    }



    // Needed by casadi_interpn
    alloc_w(ndim_, true);
    alloc_iw(2*ndim_, true);
  }

  void LinearInterpolant::eval(void* mem, const double** arg, double** res,
                               int* iw, double* w) const {
    if (res[0]) {
      res[0][0] = casadi_interpn(ndim_, get_ptr(grid_), get_ptr(offset_),
                                 get_ptr(values_), arg[0], get_ptr(lookup_mode_), iw, w);
    }
  }

  void LinearInterpolant::generateBody(CodeGenerator& g) const {
    g << "  if (res[0]) {\n"
      << "    res[0][0] = " << g.interpn(ndim_, g.constant(grid_), g.constant(offset_),
      g.constant(values_), "arg[0]", g.constant(lookup_mode_), "iw", "w") << "\n"
      << "  }\n";
  }

  Function LinearInterpolant::
  getFullJacobian(const std::string& name,
                  const std::vector<std::string>& i_names,
                  const std::vector<std::string>& o_names,
                  const Dict& opts) {
    Function ret;
    ret.assignNode(new LinearInterpolantJac(name));
    ret->construct(opts);
    return ret;
  }

  void LinearInterpolantJac::init(const Dict& opts) {
    // Call the base class initializer
    FunctionInternal::init(opts);

    // Needed by casadi_interpn
    auto m = derivative_of_.get<LinearInterpolant>();
    alloc_w(2*m->ndim_, true);
    alloc_iw(2*m->ndim_, true);
  }

  void LinearInterpolantJac::eval(void* mem, const double** arg, double** res,
                               int* iw, double* w) const {
    auto m = derivative_of_.get<LinearInterpolant>();
    casadi_interpn_grad(res[0], m->ndim_, get_ptr(m->grid_), get_ptr(m->offset_),
                        get_ptr(m->values_), arg[0], get_ptr(m->lookup_mode_), iw, w);
  }


  void LinearInterpolantJac::generateBody(CodeGenerator& g) const {

    auto m = derivative_of_.get<LinearInterpolant>();

    g << "  " << g.interpn_grad("res[0]", m->ndim_,
      g.constant(m->grid_), g.constant(m->offset_), g.constant(m->values_),
      "arg[0]", g.constant(m->lookup_mode_), "iw", "w") << "\n";
  }

} // namespace casadi
