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

#ifndef SIMULATOR_INTERNAL_HPP
#define SIMULATOR_INTERNAL_HPP

#include "simulator.hpp"
#include "fx_internal.hpp"

namespace CasADi{

/** \brief Simulator data storage classs
  \author Joel Andersson 
  \date 2010
*/
class SimulatorInternal : public FXInternal{
public:
  
  /** \brief  Constructor */
  SimulatorInternal(const Integrator& integrator, const FX& output_fcn, const std::vector<double>& grid);
  
  /** \brief  Destructor */
  virtual ~SimulatorInternal();
  
  /** \brief  Clone */
  virtual SimulatorInternal* clone() const{ return new SimulatorInternal(deepcopy(integrator_),deepcopy(output_fcn_),grid_);}

  /** \brief  initialize */
  virtual void init();

  /** \brief  Integrate */
  virtual void evaluate();
  
  Integrator integrator_;
  FX output_fcn_;
  
  std::vector<double> grid_;
  
  std::vector< Matrix<double> > states_;
};
  
} // namespace CasADi

#endif // SIMULATOR_INTERNAL_HPP
