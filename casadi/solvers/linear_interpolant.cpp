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
    // Quick return?
    if (!res[0]) return;

    // Work vectors
    double* alpha = w; w += ndim_;
    int* index = iw; iw += ndim_;
    int* corner = iw; iw += ndim_;

    // Left index and fraction of interval
    for (int i=0; i<ndim_; ++i) {
      // input
      double x = arg[0] ? arg[0][i] : 0;
      // Grid
      const double* g = get_ptr(grid_) + offset_[i];
      int ng = offset_[i+1]-offset_[i];
      // Find left index
      int j = index[i] = low(x, g, ng);
      // Get interpolation/extrapolation alpha
      alpha[i] = (x-g[j])/(g[j+1]-g[j]);
    }
    // Return value
    double ret = 0;

    // Loop over all corners, add contribution to output
    casadi_fill(corner, ndim_, 0);
    do {
      // Get weight and value for corner
      double w=1;
      int ld=1; // leading dimension
      int ind=0;
      for (int i=0; i<ndim_; ++i) {
        if (corner[i]) {
          w *= alpha[i];
        } else {
          w *= 1-alpha[i];
        }
        ind += (index[i]+corner[i])*ld;
        ld *= offset_[i+1]-offset_[i];
      }

      // Add contribution to return value
      ret += w*values_.at(ind);
    } while (flip(corner, ndim_));

    // Return interpolation
    res[0][0] = ret;
  }

} // namespace casadi
