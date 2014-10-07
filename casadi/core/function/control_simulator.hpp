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


#ifndef CASADI_CONTROLSIMULATOR_HPP
#define CASADI_CONTROLSIMULATOR_HPP

#include "integrator.hpp"

namespace casadi {

  /// Input arguments of an ODE/DAE function [controldaeIn]
  enum ControlledDAEInput {
    /// Global physical time. (1-by-1) [t]
    CONTROL_DAE_T,
    /// State vector (dimension nx-by-1).
    /// Should have the same amount of
    /// non-zeros as DAEOutput:DAE_RES [x]
    CONTROL_DAE_X,
    /// Algebraic state vector (dimension np-by-1). [z]
    CONTROL_DAE_Z,
    /// Parameter vector (dimension np-by-1). [p]
    CONTROL_DAE_P,
    /// Control vector (dimension nu-by-1). [u]
    CONTROL_DAE_U,
    /// Control vector, linearly interpolated (dimension nu-by-1). [u_interp]
    CONTROL_DAE_U_INTERP,
    /// State vector (dimension nx-by-1) at the last major time-step [x_major]
    CONTROL_DAE_X_MAJOR,
    /// Time at start of control interval (1-by-1) [t0]
    CONTROL_DAE_T0,
    /// Time at end of control interval (1-by-1) [tf]
    CONTROL_DAE_TF,
    /// Number of arguments.
    CONTROL_DAE_NUM_IN
  };

  /// Input arguments of a control simulator [controlsimulatorIn]
  enum ControlSimulatorInput {
    /// Differential or algebraic state at t0  (dimension nx-by-1) [x0]
    CONTROLSIMULATOR_X0,
    /// Parameters that are fixed over the entire horizon  (dimension np-by-1) [p]
    CONTROLSIMULATOR_P,
    /// Parameters that change over the integration intervals (dimension nu-by-(ns-1)) [u]
    CONTROLSIMULATOR_U,
    /// Number of input arguments of a piecewise simulator
    CONTROLSIMULATOR_NUM_IN};

  // Forward declaration of internal class
  class ControlSimulatorInternal;

  /** \brief Piecewise Simulation class

      A ControlSimulator can be seen as a chain of Simulators whereby some parameters change
      from one Simulator to the next.

      These changing parameters can typically be interpreted as "controls" in the context
      of dynamic optimization.

      We discriminate between the following time steps:
      * Major time-steps. These are the time steps provided by the supplied grid.
        Controls are constant inbetween major time-steps\n
      * Minor time-steps. These are time steps linearly interpolated from one major time-step
        to the next. The option 'nf' regulates how many minor time-steps are taken.\n
      * Integration time-steps. Time steps that the supplied integrator might choose to
        integrate the continuous dynamics.
        They are not important what ControlSimulator is concerned.\n

      np  Number of parameters
      nu  Number of controls
      ns  The number of major grid points, as supplied in the constructor
      nf  The number of minor grid points per major interval

      \author Joris Gillis
      \date 2011
  */

  class CASADI_CORE_EXPORT ControlSimulator : public Function {
  public:

    /// Default constructor
    ControlSimulator();

    /** \brief Creates a piecewise simulator
     * \param ffcn Continuous time dynamics, an casadi::Function with the following mapping:
     * \copydoc scheme_ControlledDAEInput
     * \copydoc scheme_DAEOutput
     *
     *
     * \param output_fcn output function which maps ControlledDAEInput or DAEInput to n outputs.
     * \copydoc scheme_DAEInput
     * \copydoc scheme_ControlledDAEInput
     * \param grid the major time grid
     */
    ControlSimulator(const Function& dae, const Function& output_fcn,
                     const std::vector<double>& grid);
    ControlSimulator(const Function& dae, const Function& output_fcn,
                     const Matrix<double>& grid);

    /// Output function equal to the state
    ControlSimulator(const Function& dae, const std::vector<double>& grid);
    ControlSimulator(const Function& dae, const Matrix<double>& grid);

    /// Access functions of the node.
    ControlSimulatorInternal* operator->();

    /// Const access functions of the node.
    const ControlSimulatorInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /** Get the (minor) time grid
     *  The length is (ns-1)*nf + 1
     */
    std::vector<double> getMinorT() const;

    /** \brief Get the controls, sampled on the minor timescale.
     * Number of rows is (ns-1)*nf
     */
    Matrix<double> getMinorU() const;


    /** \brief Get the index i such that <tt>gridminor[i] == gridmajor</tt>
     */
    std::vector< int > getMajorIndex() const;

  };

} // namespace casadi

#endif // CASADI_CONTROLSIMULATOR_HPP
