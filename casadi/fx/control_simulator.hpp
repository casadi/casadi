/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef CONTROLSIMULATOR_HPP
#define CONTROLSIMULATOR_HPP

#include "integrator.hpp"

namespace CasADi{

/// Input arguments of an ODE/DAE function
enum ControlledDAEInput{
  /** Time. (1-by-1) */
  CONTROL_DAE_T,
  /** State vector (dimension nx-by-1). Should have same amount of non-zeros as DAEOutput:DAE_RES */
  CONTROL_DAE_Y,
  /** Parameter vector (dimension np-by-1). */
  CONTROL_DAE_P,
  /** Control vector (dimension nu-by-1). */
  CONTROL_DAE_U,
  /** State derivative vector (dimension nx-by-1). Should have same amount of non-zeros as DAEOutput:DAE_RES */
  CONTROL_DAE_YDOT,
  /** State vector (dimension nx-by-1) at the last major time-step */
  CONTROL_DAE_Y_MAJOR,
  /** Number of arguments. */
  CONTROL_DAE_NUM_IN
};

/// Input arguments of an integrator
enum ControlSimulatorInput{
  /** Differential or algebraic state at t0  (dimension nx-by-1) */
  CONTROLSIMULATOR_X0, 
  /** Parameters that are fixed over the entire horizon  (dimension np-by-1) */
  CONTROLSIMULATOR_P, 
  /** Parameters that change over the integration intervals (dimension (ns-1)-by-nu) */
  CONTROLSIMULATOR_U, 
  /** State derivative at t0  (dimension nx-by-1)
  * Only relevant for implicit integrators.
  */
  CONTROLSIMULATOR_XP0, 
  /** Number of input arguments of a piecewise simulator */
  CONTROLSIMULATOR_NUM_IN};
  
// Forward declaration of internal class
class ControlSimulatorInternal;

/** \brief Piecewise Simulation class
  A ControlSimulator can be seen as a chain of Simulators whereby some parameters change from one Simulator to the next.
  
  These changing parameters can typically be interpreted as "controls" in the context of dynamic optimization.
  
  
  We discriminate between the following time steps:
   * Major time-steps. These are the time steps provided by the supplied grid. Controls are constant inbetween major time-steps\n
   * Minor time-steps. These are time steps linearly interpolated from one major time-step to the next. The option 'nf' regulates how many minor time-steps are taken.\n
   * Integration time-steps. Time steps that the supplied integrator might choose to integrate the continous dynamics. They are not important what ControlSimulator is concerned.\n

  np  Number of parameters
  nu  Number of controls
  ns  The number of major grid points, as supplied in the constructor
  nf  The number of minor grid points per coarse grid point interval
  
  \author Joris Gillis 
  \date 2011
*/

class ControlSimulator : public FX{
public:

  /// Default constructor 
  ControlSimulator();
  
    /** \brief Creates a piecewise simulator
    * \param ffcn Continuous time dynamics, an CasADi::FX with the folowing mapping:
    * \copydoc scheme_ControlledDAEInput
    * \copydoc scheme_DAEOutput
    *
    *
    * \param output_fcn output function which maps ControlledDAEInput to n outputs.
    * \param grid  the major time grid
    * \copydoc scheme_DAEInput
    * 
    */
  ControlSimulator(const FX& dae, const FX& output_fcn, const std::vector<double>& grid);
  
  /// Output function equal to the state
  ControlSimulator(const FX& dae, const std::vector<double>& grid);

  /// Access functions of the node.
  ControlSimulatorInternal* operator->();

  /// Const access functions of the node.
  const ControlSimulatorInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /** Get the (minor) time grid
  *  The length is (ns-1)*nf + 1
  */
 	std::vector<double> getMinorT() const; 
 	
  /** \brief Get the controls, sampled on the fine timescale.
  * Number of rows is (ns-1)*nf
  */
 	Matrix<double> getMinorU() const; 
 	
 	
 	/** \brief Get the index i such that gridminor[i] == gridmajor 
  */
 	std::vector< int > getMajorIndex() const; 
 	
};
  
} // namespace CasADi

#endif //CONTROLSIMULATOR_HPP
