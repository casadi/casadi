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
                    const vector<vector<double> >& grid,
                    const vector<double>& values) : Interpolant(name, grid, values) {
  }

  LinearInterpolant::~LinearInterpolant() {
  }

  void LinearInterpolant::init(const Dict& opts) {
    // Call the base class initializer
    Interpolant::init(opts);

    // Temporary memory
    alloc_w(grid_.size(), true); // alpha
    alloc_iw(grid_.size(), true); // index
    alloc_iw(grid_.size(), true); // corner
  }

  void LinearInterpolant::eval(void* mem, const double** arg, double** res,
                               int* iw, double* w) const {
    // Quick return?
    if (!res[0]) return;

    // Number of Dimensions
    int ndim = grid_.size();

    // Work vectors
    double* alpha = w; w += ndim;
    int* index = iw; iw += ndim;
    int* corner = iw; iw += ndim;

    // Left index and fraction of interval
    for (int i=0; i<ndim; ++i) {
      // input
      double x = arg[0] ? arg[0][i] : 0;
      // Grid
      const double* g = get_ptr(grid_[i]);
      int ng = grid_[i].size();
      // Find left index
      int j = index[i] = low(x, g, ng);
      // Get interpolation/extrapolation alpha
      alpha[i] = (x-g[j])/(g[j+1]-g[j]);
    }
    // Return value
    double ret = 0;

    // Loop over all corners, add contribution to output
    casadi_fill(corner, ndim, 0);
    do {
      // Get weight and value for corner
      double w=1;
      int ld=1; // leading dimension
      int ind=0;
      for (int i=0; i<ndim; ++i) {
        if (corner[i]) {
          w *= alpha[i];
        } else {
          w *= 1-alpha[i];
        }
        ind += (index[i]+corner[i])*ld;
        ld *= grid_[i].size();
      }

      // Add contribution to return value
      ret += w*values_.at(ind);
    } while (flip(corner, ndim));

    // Return interpolation
    res[0][0] = ret;
  }

} // namespace casadi
