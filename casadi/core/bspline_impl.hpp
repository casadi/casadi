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


#ifndef CASADI_BSPLINE_IMPL_HPP
#define CASADI_BSPLINE_IMPL_HPP

#include "bspline.hpp"
#include "interpolant_impl.hpp"
#include "casadi_low.hpp"

namespace casadi {

    template<class M>
    M BSplineCommon::derivative_coeff(casadi_int i, const M& coeffs) const {
      casadi_int n_dims = degree_.size();

      casadi_int n_knots = offset_[i+1]-offset_[i];
      casadi_int n = n_knots-degree_[i]-1;
      DM knots = std::vector<double>(get_ptr(knots_)+offset_[i], get_ptr(knots_)+offset_[i+1]);
      DM delta_knots = knots(range(1+degree_[i], n_knots-1))
            - knots(range(1, n_knots-degree_[i]-1));
      Sparsity sp_diag = vertsplit(Sparsity::diag(n), {0, n-1, n})[0];
      Sparsity sp_band = vertsplit(Sparsity::band(n, -1), {0, n-1, n})[0];
      DM delta_knots_inv = 1/delta_knots;
      DM T = DM(sp_diag, -delta_knots_inv) + DM(sp_band, delta_knots_inv);
      T*= degree_[i];

      std::vector<casadi_int> coeffs_dims_new = coeffs_dims_;
      coeffs_dims_new[i+1] = T.size1();

      // Apply transformation T on axis i

      // Bring axis i to the back
      std::vector<casadi_int> order = range(n_dims+1);
      std::swap(order.back(), order[i+1]);
      std::vector<casadi_int> mapping = tensor_permute_mapping(coeffs_dims_, order);
      M coeff_matrix = coeffs.nz(mapping); // NOLINT(cppcoreguidelines-slicing)

      // Cast as matrix
      coeff_matrix = reshape(coeff_matrix, -1, T.size2());

      // Apply the transformation matrix from the right
      coeff_matrix = mtimes(coeff_matrix, T.T());

      // Bring axis i back to the original place
      mapping = tensor_permute_mapping(permute(coeffs_dims_new, order), order);
      coeff_matrix = coeff_matrix.nz(mapping); // NOLINT(cppcoreguidelines-slicing)

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
    
    // Loop over dimensions
    for (casadi_int k=0;k<n_dims;++k) {
      std::vector< std::vector<double> > knots;
      std::vector< casadi_int> degree;
      for (casadi_int i=0;i<degree_.size();++i) {
        if (i==k) {
          knots.push_back(
            std::vector<double>(get_ptr(knots_)+offset_[i]+1, get_ptr(knots_)+offset_[i+1]-1));
          degree.push_back(degree_[i]-1);
        } else {
          knots.push_back(
            std::vector<double>(get_ptr(knots_)+offset_[i], get_ptr(knots_)+offset_[i+1]));
          degree.push_back(degree_[i]);
        }
      }
      MX d = MX::bspline(x, derivative_coeff(k, coeffs), knots, degree, m_, opts);
      parts.push_back(d);
    }

    return horzcat(parts);
  }

} // namespace casadi

#endif // CASADI_BSPLINE_IMPL_HPP
