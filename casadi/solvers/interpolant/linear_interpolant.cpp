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
    plugin->version = 30;
    return 0;
  }

  extern "C"
  void CASADI_INTERPOLANT_LINEAR_EXPORT casadi_load_interpolant_linear() {
    Interpolant::registerPlugin(casadi_register_interpolant_linear);
  }

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

    // Needed by casadi_interpn
    alloc_w(ndim_, true);
    alloc_iw(2*ndim_, true);
  }

  void LinearInterpolant::eval(void* mem, const double** arg, double** res,
                               int* iw, double* w) const {
    if (res[0]) {
      res[0][0] = casadi_interpn(ndim_, get_ptr(grid_), get_ptr(offset_),
                                 get_ptr(values_), arg[0], iw, w);
    }
  }

  Function LinearInterpolant::
  getFullJacobian(const std::string& name, const Dict& opts) {
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
                        get_ptr(m->values_), arg[0], iw, w);
  }

} // namespace casadi
