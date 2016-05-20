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
    // Number of Dimensions
    int ndim = grid_.size();
    // Left index and fraction of interval
    double* alpha = w; w += ndim;
    int* index = iw; iw += ndim;
    for (int i=0; i<ndim; ++i) {
      // Number of gridpoints
      int ng = grid_[i].size();
      // Get x-value
      double x = arg[0] ? arg[0][i] : 0;
      // Find left index
      int ind = -1;
      for (int j=0; j<ng; ++j, ++ind) {
        if (x < grid_[i][j]) break;
      }
      // Confine to [0, ng-2]
      index[i] = ind = max(0, min(ind, ng-2));
      // Get interpolation/extrapolation alpha
      alpha[i] = (x-grid_[i][ind])/(grid_[i][ind+1]-grid_[i][ind]);
    }
    // Loop over all corners, add contribution to output
    double ret = 0;
    int* corner = iw; iw += ndim;
    casadi_fill(corner, ndim, 0);
    while (true) {
      // Get weight and value for corner
      double w=1;
      int ld=1; // leading dimension
      int ind=0;
      for (int i=0; i<ndim; ++i) {
        w *= corner[i] ? alpha[i] : 1-alpha[i];
        ind += (index[i]+corner[i])*ld;
        ld *= grid_[i].size();
      }

      // Add contribution to return value
      ret += w*values_.at(ind);

      // Find next corner
      bool done = true;
      for (int i=0; i<ndim; ++i) {
        if (corner[i]) {
          corner[i]=0;
        } else {
          corner[i]=1;
          done = false;
          break;
        }
      }
      // Has all corners been visited?
      if (done) break;
    }

    // Return interpolation
    casadi_copy(&ret, 1, res[0]);
  }

} // namespace casadi
