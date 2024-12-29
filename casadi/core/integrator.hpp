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


#ifndef CASADI_INTEGRATOR_HPP
#define CASADI_INTEGRATOR_HPP

#include "function.hpp"
#include "linsol.hpp"
#include "rootfinder.hpp"
#include "casadi_enum.hpp"

namespace casadi {

  /** \defgroup main_integrator Title
      \par

      Create an ODE/DAE integrator
      Solves an initial value problem (IVP) coupled to a terminal value problem
      with differential equation given as an implicit ODE coupled to an algebraic
      equation and a set of quadratures:

      \verbatim
      Initial conditions at t=t0
      x(t0)  = x0
      q(t0)  = 0

      Forward integration from t=t0 to t=tf
      der(x) = function(x, z, p, t)                  Forward ODE
      0 = fz(x, z, p, t)                  Forward algebraic equations
      der(q) = fq(x, z, p, t)                  Forward quadratures

      Terminal conditions at t=tf
      rx(tf)  = rx0
      rq(tf)  = 0

      Backward integration from t=tf to t=t0
      der(rx) = gx(rx, rz, rp, x, z, p, t)        Backward ODE
      0 = gz(rx, rz, rp, x, z, p, t)        Backward algebraic equations
      der(rq) = gq(rx, rz, rp, x, z, p, t)        Backward quadratures

      where we assume that both the forward and backwards integrations are index-1
      (i.e. dfz/dz, dgz/drz are invertible) and furthermore that
      gx, gz and gq have a linear dependency on rx, rz and rp.

      \endverbatim

      \generalsection{Integrator}
      \pluginssection{Integrator}

      \author Joel Andersson
      \date 2011-2015

      \identifier{21k} */
  /** \defgroup integrator Title
  * @copydoc main_integrator
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_integrator
  * \endif
  */
  ///@{
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const SXDict& dae, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const SXDict& dae, double t0, const std::vector<double>& tout, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, double t0, const std::vector<double>& tout, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, double t0, const std::vector<double>& tout, const Dict& opts=Dict());
#ifndef SWIGMATLAB
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const SXDict& dae, double t0, double tf, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, double t0, double tf, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, double t0, double tf, const Dict& opts=Dict());
#endif // SWIGMATLAB
  ///@}

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_integrator(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_integrator(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_integrator(const std::string& name);

  /** \brief Get input scheme of integrators

      \identifier{7b} */
  CASADI_EXPORT std::vector<std::string> integrator_in();

  /** \brief Get integrator output scheme of integrators

      \identifier{7c} */
  CASADI_EXPORT std::vector<std::string> integrator_out();

  /** \brief Get integrator input scheme name by index

      \identifier{7d} */
  CASADI_EXPORT std::string integrator_in(casadi_int ind);

  /** \brief Get output scheme name by index

      \identifier{7e} */
  CASADI_EXPORT std::string integrator_out(casadi_int ind);

  /** \brief Get the number of integrator inputs

      \identifier{7f} */
  CASADI_EXPORT casadi_int integrator_n_in();

  /** \brief Get the number of integrator outputs

      \identifier{7g} */
  CASADI_EXPORT casadi_int integrator_n_out();

  /** \brief Get input scheme of a DAE function

      \identifier{25p} */
  CASADI_EXPORT std::vector<std::string> dyn_in();

  /** \brief Get output scheme of a DAE function

      \identifier{25q} */
  CASADI_EXPORT std::vector<std::string> dyn_out();

  /** \brief Get input scheme of a DAE function by index

      \identifier{25r} */
  CASADI_EXPORT std::string dyn_in(casadi_int ind);

  /** \brief Get output scheme of a DAE function by index

      \identifier{25s} */
  CASADI_EXPORT std::string dyn_out(casadi_int ind);

  /** \brief Get the number of inputs for a DAE function

      \identifier{25t} */
  CASADI_EXPORT casadi_int dyn_n_in();

  /** \brief Get the number of outputs for a DAE function

      \identifier{25u} */
  CASADI_EXPORT casadi_int dyn_n_out();

  /** \brief Get input scheme of an event transition function

      \identifier{2b4} */
  CASADI_EXPORT std::vector<std::string> event_in();

  /** \brief Get output scheme of an event transition functions

      \identifier{2b5} */
  CASADI_EXPORT std::vector<std::string> event_out();
  /** @} */

#ifndef SWIG
/// Inputs of the symbolic representation of the DAE
enum DynIn {
  DYN_T,
  DYN_X,
  DYN_Z,
  DYN_P,
  DYN_U,
  DYN_NUM_IN};

/// Outputs of the symbolic representation of the DAE
enum DynOut {
  DYN_ODE,
  DYN_ALG,
  DYN_QUAD,
  DYN_ZERO,
  DYN_NUM_OUT};

/// Inputs of an event transition function
enum EventIn {
  EVENT_INDEX,
  EVENT_T,
  EVENT_X,
  EVENT_Z,
  EVENT_P,
  EVENT_U,
  EVENT_NUM_IN};

/// Outputs of an event transition function
enum EventOut {
  EVENT_POST_X,
  EVENT_POST_Z,
  EVENT_NUM_OUT};

/// Input arguments of an integrator
enum IntegratorInput {
  /// Differential state at the initial time
  INTEGRATOR_X0,
  /// Initial guess for the algebraic variable at the initial time
  INTEGRATOR_Z0,
  /// Parameters
  INTEGRATOR_P,
  /// Piecewise constant control, a new control interval starts at each output time
  INTEGRATOR_U,
  /// Adjoint seeds corresponding to the states at the output times
  INTEGRATOR_ADJ_XF,
  /// Adjoint seeds corresponding to the algebraic variables at the output times
  INTEGRATOR_ADJ_ZF,
  /// Adjoint seeds corresponding to the quadratures at the output times
  INTEGRATOR_ADJ_QF,
  /// Number of input arguments of an integrator
  INTEGRATOR_NUM_IN
};

/// Output arguments of an integrator
enum IntegratorOutput {
  /// Differential state at all output times
  INTEGRATOR_XF,
  /// Algebraic variable at all output times
  INTEGRATOR_ZF,
  /// Quadrature state at all output times
  INTEGRATOR_QF,
  /// Adjoint sensitivities corresponding to the initial state
  INTEGRATOR_ADJ_X0,
  /// Adjoint sensitivities corresponding to the algebraic variable guess
  INTEGRATOR_ADJ_Z0,
  /// Adjoint sensitivities corresponding to the parameter vector
  INTEGRATOR_ADJ_P,
  /// Adjoint sensitivities corresponding to the control vector
  INTEGRATOR_ADJ_U,
  /// Number of output arguments of an integrator
  INTEGRATOR_NUM_OUT
};

///@{
/// Number of entries in enums
template<> struct enum_traits<DynIn> {
  static const size_t n_enum = DYN_NUM_IN;
};
template<> struct enum_traits<DynOut> {
  static const size_t n_enum = DYN_NUM_OUT;
};
template<> struct enum_traits<EventIn> {
  static const size_t n_enum = EVENT_NUM_IN;
};
template<> struct enum_traits<EventOut> {
  static const size_t n_enum = EVENT_NUM_OUT;
};
///@}

///@{
/// Convert to string
CASADI_EXPORT std::string to_string(DynIn v);
CASADI_EXPORT std::string to_string(DynOut v);
CASADI_EXPORT std::string to_string(EventIn v);
CASADI_EXPORT std::string to_string(EventOut v);
///@}

#endif // SWIG

} // namespace casadi

#endif // CASADI_INTEGRATOR_HPP
