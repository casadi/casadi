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

/// Input arguments of an integrator
enum ControlSimulatorInput{
  /** Differential or algebraic state at t0  (dimension nx-by-1) */
  PW_SIMULATOR_X0, 
  /** Parameters that are fixed over the entire horizon  (dimension np-by-1) */
  PW_SIMULATOR_P, 
  /** Parameters that change over the integration intervals (dimension (ns-1)-by-nv) */
  PW_SIMULATOR_V, 
  /** State derivative at t0  (dimension nx-by-1)
  * Only relevant for implicit integrators.
  */
  PW_SIMULATOR_XP0, 
  /** Number of input arguments of a piecewise simulator */
  PW_SIMULATOR_NUM_IN};
  
// Forward declaration of internal class
class ControlSimulatorInternal;

/** \brief Piecewise Simulation class
  A "piecewise simulator" can be seen as a chain of Simulators whereby some parameters change from one Simulator to the next.
  
  These changing parameters can typically be intrepreted as "controls" in the context of dynamic optimization.
  
  The ControlSimulator starts from a suplied continuous dynamics function, dae.
  
  We discriminate between the following tim steps:
   * Coarse time-steps. These are the time steps provided by the supplied grid.
   * Fine time-steps. These are time steps linearly interpolated from one coarse time-step to the next. The option 'nf' regulates how many fine time-steps are taken.
   * Integration time-steps. Time steps that the supplied integrator might choose to integrate the continous dynamics. They are not important what ControlSimulator is concerned.
  
  We divide the set of parameters dae.input(DAE) into static and non-static, i.e. parameters that are fixed over the entire time horizon and parameters that change at each coarse time steps. This division is carried out by an integer scalar, option 'np', that denotes the number of static parameters.
  
  _______________________
  |                     |
  |  np static params   |
  |_____________________|   ns in total
  |                     |   
  |  nv changing params |
  |_____________________|
  

  np  Number of parameters that do not change over the entire intergation time
  nv  Number of parameters that change during the integration time, but not on a fine level
  ns  The number of (coarse) grid points, as supplied in the constructor
  nf  The number of fine grid points per coarse grid point interval
  
  
  \author Joris Gillis 
  \date 2011
*/

class ControlSimulator : public FX{
public:

  /// Default constructor 
  ControlSimulator();
  
    /** \brief Creates a piecewise simulator
    * \param ffcn Continuous time dynamics, an CasADi::FX with the folowing mapping:
    * \copydoc scheme_DAEInput
    * \copydoc scheme_DAEOutput
    * Important notes:
    *  - In the above table, INTEGRATOR_P input is not really of shape (np x 1), but rather ( (np+nv) x 1 ).
    *  - The first np entries of the INTEGRATOR_P input are interpreted as parameters that are fixed on the whole domain. The remainder are interpreted as parameters that change at each coarse time-step. 
    *
    * \param output_fcn output function which maps to n outputs.
    * \param grid  the coarse time grid
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
  
  /** Get the (fine-scaled) time grid
  *  The length is (ns-1)*nf + 1
  */
 	std::vector<double> getGrid() const; 
 	
  /** \brief Get the parameters that change on a coarse time scale, sampled on the fine timescale.
  * Number of rows is (ns-1)*nf
  */
 	Matrix<double> getVFine() const; 
 	
};
  
} // namespace CasADi

#endif //CONTROLSIMULATOR_HPP
