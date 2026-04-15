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


#ifndef CASADI_BSPLINE_IMPL_HPP
#define CASADI_BSPLINE_IMPL_HPP

#include "bspline.hpp"
#include "interpolant_impl.hpp"
#include "casadi_low.hpp"

namespace casadi {

    template<class M>
    M BSplineCommon::derivative_coeff(casadi_int i,
        const std::vector< std::vector<double> >& knots,
        const std::vector<casadi_int>& degree,
        const std::vector<casadi_int>& coeffs_dims, const M& coeffs,
        std::vector< std::vector<double> >& new_knots,
        std::vector<casadi_int>& new_degree) {
      casadi_int n_dims = degree.size();

      casadi_int n_knots = knots[i].size();
      casadi_int n = n_knots-degree[i]-1;
      DM K = knots[i];
      DM delta_knots = K(range(1+degree[i], n_knots-1))
            - K(range(1, n_knots-degree[i]-1));
      DM d = degree[i]/delta_knots;     // length n-1

      std::vector<casadi_int> coeffs_dims_new = coeffs_dims;
      coeffs_dims_new[i+1] = n-1;

      // T = diag(-d) + upper_band(+d) is a scaled finite-difference operator.
      // Apply via slice-subtract + broadcast-multiply — no T, no kron, no densify,
      // no permutation mapping baked into generated code.
      casadi_int L = 1, R = 1;
      for (casadi_int k=0; k<=i; ++k)                                          L *= coeffs_dims[k];
      for (casadi_int k=i+2; k<(casadi_int)coeffs_dims.size(); ++k)            R *= coeffs_dims[k];
      casadi_int K_sz = coeffs_dims[i+1];
      casadi_int Kp = n-1;

      M M_coeffs = reshape(coeffs, L*K_sz, R);
      M top = M_coeffs(Slice(L,   L*K_sz),    Slice());
      M bot = M_coeffs(Slice(0, L*(K_sz-1)),  Slice());
      M diffed = top - bot;

      std::vector<casadi_int> dims{L, Kp, R};
      std::vector<casadi_int> a{-1, -2, -3};
      std::vector<casadi_int> b{-2};
      std::vector<casadi_int> c{-1, -2, -3};
      M coeff_matrix = einstein(vec(diffed), M(d),
        dims, std::vector<casadi_int>{Kp}, dims,
        a, b, c);

      new_knots.clear();
      new_degree.clear();
      for (casadi_int k=0;k<degree.size();++k) {
        if (i==k) {
          new_knots.push_back(
            std::vector<double>(knots[k].begin()+1, knots[k].end()-1));
          new_degree.push_back(degree[k]-1);
        } else {
          new_knots.push_back(knots[k]);
          new_degree.push_back(degree[k]);
        }
      }

      // Return the flat vector
      return coeff_matrix;
    }
    
    template<class T>
    MX BSplineCommon::jac(const MX& x, const T& coeffs) const {
    casadi_int n_dims = degree_.size();
    std::vector<MX> parts;

    Dict opts;
    std::vector<std::string> lookup_mode;
    for (auto e : lookup_mode_) lookup_mode.push_back(Low::lookup_mode_from_enum(e));
    opts["lookup_mode"] = lookup_mode;

    // Unflatten knots
    std::vector< std::vector<double> > knots_unflat(n_dims);
    for (casadi_int k=0;k<n_dims;++k) {
      knots_unflat[k] = std::vector<double>(
        get_ptr(knots_)+offset_[k], get_ptr(knots_)+offset_[k+1]);
    }

    // Loop over dimensions
    for (casadi_int k=0;k<n_dims;++k) {
      std::vector< std::vector<double> > knots;
      std::vector< casadi_int> degree;
      T dC = derivative_coeff(k, knots_unflat, degree_, coeffs_dims_, coeffs, knots, degree);
      MX d = MX::bspline(x, dC, knots, degree, m_, opts);
      parts.push_back(d);
    }

    return horzcat(parts);
  }

} // namespace casadi

#endif // CASADI_BSPLINE_IMPL_HPP