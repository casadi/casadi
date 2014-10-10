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


#ifndef CASADI_CONTROLSIMULATOR_INTERNAL_HPP
#define CASADI_CONTROLSIMULATOR_INTERNAL_HPP

#include "control_simulator.hpp"
#include "simulator.hpp"
#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief ControlSimulator data storage class
      \author Joel Andersson
      \date 2010
  */
  class CASADI_CORE_EXPORT ControlSimulatorInternal : public FunctionInternal {
  public:

    /** \brief  Constructor */
    ControlSimulatorInternal(const Function& dae, const Function& output_fcn,
                             const std::vector<double>& gridc);

    /** \brief  Destructor */
    virtual ~ControlSimulatorInternal();

    /** \brief  Clone */
    virtual ControlSimulatorInternal* clone() const {
        return new ControlSimulatorInternal(deepcopy(dae_), deepcopy(output_fcn_), gridc_);}

    /** \brief  initialize */
    virtual void init();

    /** \brief  Integrate */
    virtual void evaluate();

    /// Get the parameters that change on a coarse time scale, sampled on the fine timescale
    Matrix<double> getVFine() const;

    /** \brief Get the index i such that <tt>gridfine[i] == gridcoarse</tt>
     */
    std::vector< int > getCoarseIndex() const;


    Integrator integrator_;
    Function dae_;
    Function control_dae_;
    Simulator simulator_;

    // original output function
    Function orig_output_fcn_;

    // adapted output function
    Function output_fcn_;

    /** \brief The hart of this class, a casadi of simulator calls */
    Function all_output_;

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

} // namespace casadi
/// \endcond

#endif // CASADI_CONTROLSIMULATOR_INTERNAL_HPP
