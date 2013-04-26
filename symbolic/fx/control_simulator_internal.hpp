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

#ifndef CONTROLSIMULATOR_INTERNAL_HPP
#define CONTROLSIMULATOR_INTERNAL_HPP

#include "control_simulator.hpp"
#include "simulator.hpp"
#include "fx_internal.hpp"

namespace CasADi{

/** \brief ControlSimulator data storage classs
  \author Joel Andersson 
  \date 2010
*/
class ControlSimulatorInternal : public FXInternal{
public:
  
  /** \brief  Constructor */
  ControlSimulatorInternal(const FX& dae, const FX& output_fcn, const std::vector<double>& gridc);
  
  /** \brief  Destructor */
  virtual ~ControlSimulatorInternal();
  
  /** \brief  Clone */
  virtual ControlSimulatorInternal* clone() const{ return new ControlSimulatorInternal(deepcopy(dae_),deepcopy(output_fcn_),gridc_);}

  /** \brief  initialize */
  virtual void init();

  /** \brief  Integrate */
  virtual void evaluate(int nfdir, int nadir);
  
  /** \brief  Update the number of sensitivity directions during or after initialization */
  virtual void updateNumSens(bool recursive);
  
  /// Get the parameters that change on a coarse time scale, sampled on the fine timescale
         Matrix<double> getVFine() const; 
         
         /** \brief Get the index i such that gridfine[i] == gridcoarse 
  */
         std::vector< int > getCoarseIndex() const; 
         

  Integrator integrator_;
  FX dae_;
  FX control_dae_;
  Simulator simulator_;
  
  // original output function
  FX orig_output_fcn_;
  
  // adapted output function
  FX output_fcn_;
  
  /** \brief The hart of this class, a casadi of simulator calls */
  FX all_output_;
  
  /** grid */
  std::vector<double> grid_;
  
  /** Coarse grid */
  std::vector<double> gridc_;
  
  /** The local non-dimensional time grid */
  std::vector<double> gridlocal_;

  /** \brief Number of states */
  int ny_;
  
  /** \brief Number of static parameters */
  int np_;
  
  /** \brief Number of controls */
  int nu_;
  
  /** \brief Number of interpolated controls */
  int nu_interp_;
  
  /** \brief Number of coarse time steps */
  int ns_;
  
  /** \brief Number of fine-grained time steps */
  int nf_;
  
};
  
} // namespace CasADi

#endif // CONTROLSIMULATOR_INTERNAL_HPP
