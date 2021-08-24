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


#ifndef CASADI_SIMULATOR_HPP
#define CASADI_SIMULATOR_HPP

#include "function.hpp"
#include "linsol.hpp"
#include "rootfinder.hpp"

namespace casadi {

  /** \defgroup main_simulator
      Create an ODE/DAE simulator
      Solves an initial value problem (IVP) coupled to a terminal value problem
      with differential equation given as an implicit ODE coupled to an algebraic
      equation and a set of quadratures:

      \verbatim
      Time grid: [t0, t1, ..., tN]

      Initial conditions for forward integration
      x(t0)  = x0
      qk(tk)  = 0, k = 0, ..., N-1

      Sequential forward integration from t=tk to t=t{k+1}
      der(x) = function(t, x, z, p, uk)   Forward ODE
      0 = fz(t, x, z, p, uk)              Forward algebraic equations
      der(qk) = fq(t, x, z, p, uk)        Forward quadratures
      yk = fy(t, x, z, p, uk)             Forward output equations

      Terminal conditions for backward integration
      rx(tf)  = rx0
      rqk(tk{k+1})  = 0, k = N-1, ..., 0

      Sequential backward integration from t=t{k+1} to t=tk
      der(rx) = gx(t, rx, rz, rp, ruk, x, z, p, uk)    Backward ODE
      0 = gz(t, rx, rz, rp, ruk, x, z, p, uk)          Backward algebraic equations
      der(rq) = gq(t, rx, rz, rp, ruk, x, z, p, uk)    Backward quadratures
      yk = gy(t, rx, rz, rp, ruk, x, z, p, uk)         Backward output equations

      where we assume that both the forward and backwards integrations are index-1
      (i.e. dfz/dz, dgz/drz are invertible) and furthermore that
      gx, gz, gq and gy have a linear dependency on rx, rz, rp and ru.

      \endverbatim

      \generalsection{Simulator}
      \pluginssection{Simulator}

      \author Joel Andersson
      \date 2011-2021
  */
  /** \defgroup simulator
  * @copydoc main_simulator
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_simulator
  * \endif
  */
  ///@{
  CASADI_EXPORT Function simulator(const std::string& name, const std::string& solver,
    const SXDict& dae, const std::vector<double>& grid, const Dict& opts=Dict());
  CASADI_EXPORT Function simulator(const std::string& name, const std::string& solver,
    const MXDict& dae, const std::vector<double>& grid, const Dict& opts=Dict());
  CASADI_EXPORT Function simulator(const std::string& name, const std::string& solver,
    const Function& dae, const std::vector<double>& grid, const Dict& opts=Dict());
  ///@}

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_simulator(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_simulator(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_simulator(const std::string& name);

  /** \brief Get input scheme of simulators */
  CASADI_EXPORT std::vector<std::string> simulator_in();

  /** \brief Get simulator output scheme of simulators */
  CASADI_EXPORT std::vector<std::string> simulator_out();

  /** \brief Get simulator input scheme name by index */
  CASADI_EXPORT std::string simulator_in(casadi_int ind);

  /** \brief Get output scheme name by index */
  CASADI_EXPORT std::string simulator_out(casadi_int ind);

  /** \brief Get the number of simulator inputs */
  CASADI_EXPORT casadi_int simulator_n_in();

  /** \brief Get the number of simulator outputs */
  CASADI_EXPORT casadi_int simulator_n_out();
  /** @} */

#ifndef SWIG
/// Inputs of the symbolic representation of the DAE
enum DynIn {
  DYN_T,
  DYN_X,
  DYN_Z,
  DYN_P,
  DYN_RX,
  DYN_RZ,
  DYN_RP,
  DYN_NUM_IN};

/// Inputs of the symbolic representation of the DAE
enum DynOut {
  DYN_ODE,
  DYN_ALG,
  DYN_Y,
  DYN_QUAD,
  DYN_RODE,
  DYN_RALG,
  DYN_RY,
  DYN_RQUAD,
  DYN_NUM_OUT};

/// Input arguments of an ODE/DAE function
enum DAE2Input {
  /// Differential state
  DAE2_X,
  /// Algebraic state
  DAE2_Z,
  /// Parameter
  DAE2_P,
  /// Explicit time dependence
  DAE2_T,
  /// Number of arguments
  DAE2_NUM_IN
};

/// Output arguments of an DAE function
enum DAE2Output {
  /// Right hand side of the implicit ODE
  DAE2_ODE,
  /// Right hand side of algebraic equations
  DAE2_ALG,
  /// Right hand side of quadratures equations
  DAE2_QUAD,
  /// Number of arguments
  DAE2_NUM_OUT
};

/// Input arguments of an ODE/DAE backward integration function
enum RDAE2Input {
  /// Backward differential state
  RDAE2_RX,
  /// Backward algebraic state
  RDAE2_RZ,
  /// Backward  parameter vector
  RDAE2_RP,
  /// Forward differential state
  RDAE2_X,
  /// Forward algebraic state
  RDAE2_Z,
  /// Parameter vector
  RDAE2_P,
  /// Explicit time dependence
  RDAE2_T,
  /// Number of arguments
  RDAE2_NUM_IN
};

/// Output arguments of an ODE/DAE backward integration function
enum RDAE2Output {
  /// Right hand side of ODE
  RDAE2_ODE,
  /// Right hand side of algebraic equations
  RDAE2_ALG,
  /// Right hand side of quadratures
  RDAE2_QUAD,
  /// Number of arguments
  RDAE2_NUM_OUT
};

/// Input arguments of an simulator
enum SimulatorInput {
  /// Differential state at the initial time
  SIMULATOR_X0,
  /// Parameters
  SIMULATOR_P,
  /// Initial guess for the algebraic variable
  SIMULATOR_Z0,
  /// Backward differential state at the final time
  SIMULATOR_RX0,
  /// Backward parameter vector
  SIMULATOR_RP,
  /// Initial guess for the backwards algebraic variable
  SIMULATOR_RZ0,
  /// Number of input arguments of an simulator
  SIMULATOR_NUM_IN
};

/// Output arguments of an simulator
enum SimulatorOutput {
  /// Differential state at the final time
  SIMULATOR_XF,
  /// Outputs
  SIMULATOR_Y,
  /// Quadratures
  SIMULATOR_QF,
  /// Algebraic variable at the final time
  SIMULATOR_ZF,
  /// Backward differential state at the initial time
  SIMULATOR_RXF,
  /// Backward outputs
  SIMULATOR_RY,
  /// Backward quadratures
  SIMULATOR_RQF,
  /// Backward algebraic variable at the initial time
  SIMULATOR_RZF,
  /// Number of output arguments of an simulator
  SIMULATOR_NUM_OUT
};
#endif // SWIG

} // namespace casadi

#endif // CASADI_SIMULATOR_HPP
