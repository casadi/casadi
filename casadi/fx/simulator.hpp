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

#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include "integrator.hpp"

namespace CasADi{

/// Input arguments of an integrator
enum SimulatorInput{
  /** Differential or algebraic state at t0  (dimension nx-by-1) */
  SIMULATOR_X0, 
  /** Parameters p  (dimension np-by-1) */
  SIMULATOR_P,  
  /** State derivative at t0  (dimension nx-by-1)
  * Only relevant for implicit intergators.
  * This input may be changed during an IDASIntegrator::evaluate()
  */
  SIMULATOR_XP0, 
  /** Parameters that change over the integration intervals (dimension (ns-1)-by-nv) */
  SIMULATOR_V, 
  /** Number of input arguments of an integrator */
  SIMULATOR_NUM_IN};
  
// Forward declaration of internal class
class SimulatorInternal;

/** \brief Integrator class
  An "simulator" integrates an IVP, stopping at a (fixed) number of grid points and 
  evaluates a set of output functions at these points.
  The internal stepsizes of the integrator need not coincide with the gridpoints.
  
  
  Simulator is an CasADi::FX mapping from CasADi::IntegratorInput to n. \\
  
  The output function needs to be a mapping from CasADi::DAEInput to n. The default output has n=1 and the output is the (flattened) differential state for each time step:  number of rows equals ns*nf
  
  np  Number of parameters that do not change over the entire intergation time
  nv  Number of parameters that change during the integration time, but not on a fine level
  ns  The number of (coarse) grid points, as suplied in the constructor
  nf  The number of fine grid points per coarse grid point interval
  
  nf = 1
  ns = 5
  
  |||||
  
  nf = 3
  ns = 5
  
  |..|..|..|..|
  
  
  \author Joel Andersson 
  \date 2010
*/

class Simulator : public FX{
public:

  /// Default constructor 
  Simulator();
  
  /** \brief Constructor
  * \param output_fcn output function which maps to n outputs.
  * \copydoc scheme_DAEInput
  * 
  */
  Simulator(const Integrator& integrator, const FX& output_fcn, const std::vector<double>& grid);
  
  /// Output function equal to the state
  Simulator(const Integrator& integrator, const std::vector<double>& grid);

  /// Access functions of the node.
  SimulatorInternal* operator->();

  /// Const access functions of the node.
  const SimulatorInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Get the (fine-scaled) time grid used by the integrator
  std::vector<double> getGrid() const;
};
  
} // namespace CasADi

#endif //SIMULATOR_HPP
