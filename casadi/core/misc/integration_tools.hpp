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


#ifndef CASADI_INTEGRATION_TOOLS_HPP
#define CASADI_INTEGRATION_TOOLS_HPP

#include <vector>
#include <algorithm>
#include "casadi/core/casadi_exception.hpp"
#include "casadi/core/options_functionality.hpp"
#include "casadi/core/mx/mx.hpp"

namespace casadi {

  ///@{
  /** \brief Obtain collocation points of specific order and scheme
  \param scheme  'radau' or 'legendre'
  **/
  CASADI_CORE_EXPORT
    std::vector<double> collocationPoints(int order, const std::string& scheme="radau");
#ifndef SWIG
  CASADI_CORE_EXPORT
    std::vector<long double> collocationPointsL(int order, const std::string& scheme="radau");
#endif // SWIG
  ///@}

  /** \brief Obtain collocation interpolating matrices
  \param tau_root  location of collocation points, as obtained from collocationPoints
  \param[out] C interpolating coefficients to obtain derivatives
      Length: order+1, order + 1

      dX/dt @collPoint(j) ~ Sum_i C[j][i]*X@collPoint(i)

  \param[out] D interpolating coefficients to obtain end state
      Length: order+1
  */
#ifndef SWIG
  CASADI_CORE_EXPORT
    void collocationInterpolators(const std::vector<double> & tau_root,
                                  std::vector< std::vector<double> > &C,
                                  std::vector< double > &D);
#else // SWIG
  CASADI_CORE_EXPORT
    void collocationInterpolators(const std::vector<double> & tau_root,
                                  std::vector< std::vector<double> > &OUTPUT,
                                  std::vector< double > &OUTPUT);
#endif // SWIG

#ifndef SWIG
extern const long double legendre_points1[2];
extern const long double legendre_points2[3];
extern const long double legendre_points3[4];
extern const long double legendre_points4[5];
extern const long double legendre_points5[6];
extern const long double legendre_points6[7];
extern const long double legendre_points7[8];
extern const long double legendre_points8[9];
extern const long double legendre_points9[10];
extern const long double* legendre_points[10];

// Radau collocation points
extern const long double radau_points1[2];
extern const long double radau_points2[3];
extern const long double radau_points3[4];
extern const long double radau_points4[5];
extern const long double radau_points5[6];
extern const long double radau_points6[7];
extern const long double radau_points7[8];
extern const long double radau_points8[9];
extern const long double radau_points9[10];
extern const long double* radau_points[10];

extern const long double** collocation_points[2];
#endif // SWIG

  // Type of collocation points
  enum CollocationPoints {LEGENDRE, RADAU};

  /** \brief Construct an explicit Runge-Kutta integrator
  * \param f dynamical system
  * \copydoc scheme_DAEInput
  * \copydoc scheme_DAEOutput
  * \param tf    Integration end time
  * \param order Order of integration
  * \param ne    Number of times the \e RK primitive is repeated over the integration interval
  */
  CASADI_CORE_EXPORT Function explicitRK(Function& f, const MX &tf=1,
                                                int order=4, int ne = 1);

  /** \brief Construct an implicit Runge-Kutta integrator
  * \param f dynamical system
  * \copydoc scheme_DAEInput
  * \copydoc scheme_DAEOutput
  * \param tf    Integration end time
  * \param order Order of integration
  * \param scheme Collocation scheme, as excepted by collocationPoints function.
  * \param ne    Number of times the \e RK primitive is repeated over the integration interval
  */
  CASADI_CORE_EXPORT
  Function implicitRK(Function& f, const std::string& impl,
                      const Dictionary& dict = Dictionary(), const MX &tf=1, int order=4,
                      const std::string& scheme="radau", int ne = 1);

} // namespace casadi

#endif // CASADI_INTEGRATION_TOOLS_HPP
