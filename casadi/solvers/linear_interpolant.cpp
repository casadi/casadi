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

    // Temporary memory
    alloc_w(ndim_, true); // alpha
    alloc_iw(ndim_, true); // index
    alloc_iw(ndim_, true); // corner
  }

  void LinearInterpolant::eval(void* mem, const double** arg, double** res,
                               int* iw, double* w) const {
    if (res[0]) {
      res[0][0] = casadi_interpn_linear(ndim_, get_ptr(grid_), get_ptr(offset_),
                                        get_ptr(values_), arg[0], iw, w);
    }
  }

} // namespace casadi
