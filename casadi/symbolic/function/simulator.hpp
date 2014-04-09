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

namespace casadi{

// Forward declaration of internal class
class SimulatorInternal;

/** \brief Integrator class

  An "simulator" integrates an IVP, stopping at a (fixed) number of grid points and
  evaluates a set of output functions at these points.
  The internal stepsizes of the integrator need not coincide with the gridpoints.


  Simulator is an casadi::Function mapping from casadi::IntegratorInput to n. \\

  The output function needs to be a mapping from casadi::DAEInput to n. The default output has n=1 and the output is the (vectorized) differential state for each time step.

  \author Joel Andersson
  \date 2010
*/

class CASADI_SYMBOLIC_EXPORT Simulator : public Function{
public:

  /// Default constructor
  Simulator();

  /** \brief Constructor
  * \param output_fcn output function which maps to n outputs.
  * \copydoc scheme_DAEInput
  *
  */
  Simulator(const Integrator& integrator, const Function& output_fcn, const std::vector<double>& grid);
  Simulator(const Integrator& integrator, const Function& output_fcn, const Matrix<double>& grid);

  /// Output function equal to the state
  Simulator(const Integrator& integrator, const std::vector<double>& grid);
  Simulator(const Integrator& integrator, const Matrix<double>& grid);

  /// Access functions of the node.
  SimulatorInternal* operator->();

  /// Const access functions of the node.
  const SimulatorInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
};

} // namespace casadi

#endif //SIMULATOR_HPP
