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

/// Indices of the inputs of the output function
enum OutputInput{OUTPUT_T, OUTPUT_X, OUTPUT_P, OUTPUT_NUM_IN};

/// Indices of the inputs of the function
enum SimulatorInput{SIMULATOR_X0, SIMULATOR_P, SIMULATOR_NUM_IN};

/// Indices of the outputs of the function
//enum SimulatorOutput{SIMULATOR_Y, SIMULATOR_XF, SIMULATOR_NUM_OUT};

/** \brief  Forward declaration of internal class */
class SimulatorInternal;

/** \brief Integrator class
  An "simulator" integrates an IVP, stopping at a (fixed) number of grid points and 
  evaluates a set of output functions at these points.
  The internal stepsizes of the integrator need not coincide with the gridpoints.

  inputs:
  0: Vector of output times (dimension ngrid)
  1: State at t0  (dimension consistent with integrator)
  2: Parameters to be passed to integrator (dimension consistent with integrator)

  outputs:
  Same as the output function with the different times as the second dimension
  
  \author Joel Andersson 
  \date 2010

*/
class Simulator : public FX{
public:
  
  /** \brief  Constructor */
  Simulator(const Integrator& integrator, const FX& output_fcn, const vector<double>& grid);
  
  /** \brief  Output function equal to the state */
  Simulator(const Integrator& integrator, const vector<double>& grid);

  /** \brief  Default constructor */
  Simulator();

  /** \brief  Copy constructors */
  Simulator(const Simulator& ref);

  //@{
  /** \brief  Access functions of the node.
  A regular user is not supposed to use these methods*/
  SimulatorInternal* operator->();
  const SimulatorInternal* operator->() const;
  //@}

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

};
  
} // namespace CasADi

#endif //SIMULATOR_HPP
