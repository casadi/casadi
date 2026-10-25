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


#ifndef CASADI_INTERPOLANT_HPP
#define CASADI_INTERPOLANT_HPP

#include "function.hpp"

namespace casadi {


  /** \defgroup main_interpolant Title
      \par

   * An interpolant function for lookup table data
   *
   * \param[in] name label for the resulting Function
   * \param[in] solver name of the plugin
   * \param[in] grid collection of 1D grids whose outer product
   *            defines the full N-D rectangular grid
   * \param[in] values flattened vector of all values
   *            for all gridpoints
   *
   * Examples in 1D and 2D; the values of a 2D grid are flattened in column-major order.
   *
   * \lang{python}
   * \doctest
   * >>> xgrid = np.linspace(1, 6, 6)
   * >>> V = [-1, -1, -2, -3, 0, 2]
   * >>> LUT = interpolant("LUT", "bspline", [xgrid], V)
   * >>> print(LUT(2.5))
   * -1.35
   * >>> xgrid = np.linspace(-5, 5, 11)
   * >>> ygrid = np.linspace(-4, 4, 9)
   * >>> X, Y = np.meshgrid(xgrid, ygrid, indexing="ij")
   * >>> R = np.sqrt(5*X**2 + Y**2) + 1
   * >>> V = np.sin(R)/R
   * >>> LUT = interpolant("LUT", "bspline", [xgrid, ygrid], V.ravel(order="F"))
   * >>> print(LUT([0.5, 1]))
   * 0.245507
   * \enddoctest
   * \endlang
   *
   * \lang{matlab}
   * \doctest
   * >> xgrid = 1:6;
   * >> V = [-1 -1 -2 -3 0 2];
   * >> LUT = interpolant('LUT', 'bspline', {xgrid}, V);
   * >> LUT(2.5)
   * -1.35
   * >> xgrid = -5:5;
   * >> ygrid = -4:4;
   * >> [X, Y] = ndgrid(xgrid, ygrid);
   * >> R = sqrt(5*X.^2 + Y.^2) + 1;
   * >> V = sin(R)./R;
   * >> LUT = interpolant('LUT', 'bspline', {xgrid, ygrid}, V(:));
   * >> LUT([0.5 1])
   * 0.245507
   * \enddoctest
   * \endlang
   *
   * \lang{julia}
   * \doctest
   * julia> xgrid = collect(1.0:6.0);
   * julia> V = [-1, -1, -2, -3, 0, 2];
   * julia> LUT = ca.interpolant("LUT", "bspline", [xgrid], V);
   * julia> LUT(2.5)
   * DM(-1.35)
   * julia> xgrid = collect(-5.0:5.0);
   * julia> ygrid = collect(-4.0:4.0);
   * julia> R = [sqrt(5x^2 + y^2) + 1 for x in xgrid, y in ygrid];
   * julia> V = sin.(R) ./ R;
   * julia> LUT = ca.interpolant("LUT", "bspline", [xgrid, ygrid], vec(V));
   * julia> LUT([0.5, 1])
   * DM(0.245507)
   * \enddoctest
   * \endlang
   *
   *  \generalsection{Interpolant}
   *  \pluginssection{Interpolant}
   * \author Joel Andersson
   * \date 2016

      \identifier{21p} */

  /** \defgroup interpolant Title
  * @copydoc main_interpolant
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_interpolant
  * \endif
  */
  ///@{
  CASADI_EXPORT Function interpolant(const std::string& name,
                                     const std::string& solver,
                                     const std::vector<std::vector<double> >& grid,
                                     const std::vector<double>& values,
                                     const Dict& opts=Dict());
  /** \brief Parametric variant of interpolant
   * 
   * The resulting function will have an additional argument for the coefficients
   * 
   * By default, derivatives wrt the coefficients are not supported (zero).
   * Some interpolant plugins may support the `inline=true` which enables correct derivatives
   * 
   *
      \identifier{1p3} */
  CASADI_EXPORT Function interpolant(const std::string& name,
                                     const std::string& solver,
                                     const std::vector<std::vector<double> >& grid,
                                     casadi_int m=1,
                                     const Dict& opts=Dict());
  ///@}

  /** \brief Parametric variant of interpolant
   * 
   * The resulting function will have additional arguments for the grid and coefficients
   *
   * By default, derivatives wrt the coefficients are not supported (zero).
   * Some interpolant plugins may support the `inline=true` which enables correct derivatives
   *
      \identifier{1p4} */
  CASADI_EXPORT Function interpolant(const std::string& name,
                                     const std::string& solver,
                                     const std::vector<casadi_int>& grid_dims,
                                     casadi_int m=1,
                                     const Dict& opts=Dict());

  /** \brief Parametric variant of interpolant
   * 
   * The resulting function will have an additional argument for the grid
   * 
   * By default, derivatives wrt the coefficients are not supported (zero).
   * Some interpolant plugins may support the `inline=true` which enables correct derivatives

      \identifier{1p5} */
  CASADI_EXPORT Function interpolant(const std::string& name,
                                     const std::string& solver,
                                     const std::vector<casadi_int>& grid_dims,
                                     const std::vector<double>& values,
                                     const Dict& opts=Dict());

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_interpolant(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_interpolant(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_interpolant(const std::string& name);

  /** @} */

} // namespace casadi

#endif // CASADI_INTERPOLANT_HPP
