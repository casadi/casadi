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


#ifndef CASADI_INTEGRATION_TOOLS_HPP
#define CASADI_INTEGRATION_TOOLS_HPP

#include "casadi/core/function.hpp"

namespace casadi {

  ///@{
  /** \brief Obtain collocation points of specific order and scheme

  \param order Which order (1 to 9 supported)
  \param scheme  'radau' or 'legendre'
  *
  \identifier{1so} */
  CASADI_EXPORT
    std::vector<double> collocation_points(casadi_int order, const std::string& scheme="radau");
#ifndef SWIG
  CASADI_EXPORT
    std::vector<long double> collocation_pointsL(casadi_int order,
      const std::string& scheme="radau");
#endif // SWIG
  ///@}

  /** \brief Obtain collocation interpolating matrices
  
  A collocation method poses a polynomial Pi that interpolates exactly through
  an initial state (0,X_0) and helper states at collocation points (tau_j,X:collPoint(j)).

  This function computes the linear mapping between dPi/dt and coefficients Z=[X_0 X:collPoints].

  \param tau  location of collocation points, as obtained from collocation_points
  \param[out] C interpolating coefficients to obtain derivatives.
      Length: order+1, order+1

    \verbatim
      dPi/dt @Z_j = (1/h) Sum_i C[j][i]*Z_i,
    \endverbatim

    with h the length of the integration interval.

  \param[out] D interpolating coefficients to obtain end state.
      Length: order+1

    \verbatim
      Pi @X_f = Sum_i D[i]*Z_i
    \endverbatim

  \identifier{1sp} */
  CASADI_EXPORT void
  collocation_interpolators(const std::vector<double> & tau,
                            std::vector< std::vector<double> > &SWIG_OUTPUT(C),
                            std::vector< double > &SWIG_OUTPUT(D));

  /** \brief Obtain collocation interpolating matrices
  
  A collocation method poses a polynomial Pi that interpolates exactly through
  an initial state (0,X_0) and helper states at collocation points (tau_j,Xc_j)
  with j=1..degree.

  This function computes the linear mapping between dPi/dt and coefficients Z=[X_0 Xc].

  \param tau  location of collocation points (length: degree), as obtained from collocation_points
  \param[out] C interpolating coefficients to obtain derivatives.
      Size: (degree+1)-by-degree

    You may find the slopes of Pi at the collocation points as
    \verbatim
      dPi/dt @ Xc = (1/h) Z*C,
    \endverbatim

    with h the length of the integration interval.

  \param[out] D interpolating coefficients to obtain end state.
      Size: (degree+1)-by-1

    You may find the end point of Pi as
    \verbatim
      Pi @X_f = Z*D
    \endverbatim

  \param[out] B quadrature coefficients
      Size: degree-by-1

    Given quadrature righ-hand-sides 'quad' evaluated at the collocation points,
    you may find the integrated quadratures as
    \verbatim
      q = quad*B*h
    \endverbatim

  \identifier{1sq} */
  CASADI_EXPORT void
  collocation_coeff(const std::vector<double> & tau,
                            DM &SWIG_OUTPUT(C),
                            DM &SWIG_OUTPUT(D),
                            DM &SWIG_OUTPUT(B));

  // Type of collocation points
  enum CollocationPoints {LEGENDRE, RADAU};

  /** \brief Construct an explicit Runge-Kutta integrator

   * The constructed function has three inputs,
   * corresponding to initial state (x0), parameter (p) and integration time (h)
   * and one output, corresponding to final state (xf).
   *
   * \param f     ODE function with two inputs (x and p) and one output (xdot)
   * \param N     Number of integrator steps
   * \param order Order of interpolating polynomials

  \identifier{1sr} */
  CASADI_EXPORT Function simpleRK(Function f, casadi_int N=10, casadi_int order=4);

  /** \brief Construct an implicit Runge-Kutta integrator using a collocation scheme

   * The constructed function has three inputs,
   * corresponding to initial state (x0), parameter (p) and integration time (h)
   * and one output, corresponding to final state (xf).
   *
   * \param f      ODE function with two inputs (x and p) and one output (xdot)
   * \param N      Number of integrator steps
   * \param order  Order of interpolating polynomials
   * \param scheme Collocation scheme, as excepted by collocationPoints function.
   * \param solver Solver plugin
   * \param solver_options Options to be passed to the solver plugin

  \identifier{1ss} */
  CASADI_EXPORT
  Function simpleIRK(Function f, casadi_int N=10, casadi_int order=4,
                      const std::string& scheme="radau",
                      const std::string& solver="newton",
                      const Dict& solver_options = Dict());

  /** \brief Simplified wrapper for the Integrator class

  \identifier{1st} */
  CASADI_EXPORT
  Function simpleIntegrator(Function f, const std::string& integrator="cvodes",
                              const Dict& integrator_options = Dict());


  /** \brief Reduce index
   *
   *
   * Index reduction leads to a new set of variables and equations.
   *
   * In the process, a set of constraints (algebraic equations or derivatives) a.k.a invariants is constructed
   * that are invariant to the problem: whenever an initial point satisfies these constraints,
   * the boundary-value-problem outcome will keep satisfying those constraints automatically,
   * even though they are *not* part of the reduced DAE.
   *
   * For any practical numerical integration method, there will be numerical drift away from satisfaction of those constraints.
   * In other words, you will see the value of invariants slowly moving away from original zero.
   * 
   * A classic mitigation technique is Baumgarte stabilization: you add these invariants to the reduced DAE
   * as a correction term that acts in a way to make small (numerical) perturbations to the invariants decay to the origin
   * as a dampened linear system.
   * 
   * in which a certain set of constraints (algebraic equations or derivatives) has been dropped in favour of 
   * 
   * 
   * 
   * 
   * \param dae Expression dictionary describing the DAE
   *   
   *   Each value must be a dense column vector.
   * 
   *   keys:
   *      - x_impl:  symbol for implicit differential states
   *      - dx_impl: symbol for implicit differential state derivatives
   *      - z:       symbol for algebraic variables
   *      - alg:     expression for algebraic equations
   *      - t:       symbol for time
   *      - p:       symbol for parameters
   * \param opts Option dictionary
   * 
   * 
   *   'baumgarte_pole': double
   *      Poles (inverse time constants) of the Baumgarte invariant correction term.
   *      Must be <0 to dampen out perturbations
   *      0 (default) amounts to no correction.
   *      Corresponds to -gamma of equation (1.5) in
   *      Ascher, Uri M., Hongsheng Chin, and Sebastian Reich. "Stabilization of DAEs and invariant manifolds." Numerische Mathematik 67.2 (1994): 131-149.
   * 
   * \param stats Statistics
   * 
   * \return Expression dictionary describing the reduced DAE
   *   
   *   In addition the fields allowed in the input DAE, the following keys occur:
   *
   *      - x:   symbol for explicit differential states
   *      - ode: expression for right-hand-side of explicit differential states
   *      - I:   expression for invariants
   *

      \identifier{23h} */
  /// @{
  CASADI_EXPORT
  MXDict dae_reduce_index(const MXDict& dae, Dict& SWIG_OUTPUT(stats), const Dict& opts={});
  CASADI_EXPORT
  SXDict dae_reduce_index(const SXDict& dae, Dict& SWIG_OUTPUT(stats), const Dict& opts={});
  /// @}


  /** \brief Turn a reduced DAE into a semi explicit form suitable for CasADi integrator
  * 
  * \param[in]  dae Original (unreduced) DAE structure
  * \param[in]  dae_red Reduced DAE (see dae_reduce_index)
  * \param[out] state_to_orig A mapping of integrator (semi explicit) states to states of the original DAE
  * \param[out] phi A function to compute the invariants of the reduced DAE
  *             Inputs:
  *               - x and z: (semi explicit) integrator states; typically integrator outputs xf and zf
  *               - p: parameters
  *               - t: time
  * \return Semi explicit DAE dictionary, suitable to pass to a CasADi integrator
  * 
  * \sa dae_reduce_index

  \identifier{1su} */
  /// @{
  CASADI_EXPORT
  MXDict dae_map_semi_expl(const MXDict& dae, const MXDict& dae_red,
    Function& SWIG_OUTPUT(state_to_orig), Function& SWIG_OUTPUT(phi));
  CASADI_EXPORT
  SXDict dae_map_semi_expl(const SXDict& dae, const SXDict& dae_red,
    Function& SWIG_OUTPUT(state_to_orig), Function& SWIG_OUTPUT(phi));
  /// @}

/** \brief Obtain a generator Function for producing consistent initial guesses of a reduced DAE
  * 
  * \param[in] dae Original (unreduced) DAE structure
  * \param[in] dae_red Reduced DAE (see dae_reduce_index) 
  * \param[in] init_solver NLP solver plugin name for nlpsol used to construct an initial guess
  * \param[in] init_strength Influence the nature of the NLP
  *               Structure with keys x_impl, dx_impl, z corresponding to inputs of init_gen
  *               Each key maps to a DM that should match the variable size corresponding to that key.
  *               For each variable the meaning of the corresponding DM value is as follows:
  *                When >=0, indicates that the provided initial guess is used in a quadratic penalty (value used as weight)
  *                When -1, indicates that the provided initial guess must be observed (simple bound on variable)
  * \param[in] init_solver_options NLP solver options to be passed to nlpsol
  *
  * \return init_gen A function to generate a consistent initial guess that can
  *             be used to pass to an integrator constructed from a semi explict reduced DAE
  *             Inputs:
  *               - x_impl, dx_impl, z: initial guesses in the original DAE space
  *               - p: parameters
  *               - t: time
  *             Outputs:
  *               - x0, z0: (semi explicit) integrator states and algebraic variables;
  *                         typically used as input for integrators
  * 
  * \sa dae_reduce_index

  \identifier{1sv} */
  /// @{
  CASADI_EXPORT
  Function dae_init_gen(const MXDict& dae, const MXDict& dae_red,
    const std::string& init_solver, const DMDict& init_strength=DMDict(),
    const Dict& init_solver_options=Dict());
  CASADI_EXPORT
  Function dae_init_gen(const SXDict& dae, const SXDict& dae_red,
    const std::string& init_solver, const DMDict& init_strength=DMDict(),
    const Dict& init_solver_options=Dict());
  /// @}

} // namespace casadi

#endif // CASADI_INTEGRATION_TOOLS_HPP
