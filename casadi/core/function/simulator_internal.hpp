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


#ifndef CASADI_SIMULATOR_INTERNAL_HPP
#define CASADI_SIMULATOR_INTERNAL_HPP

#include "simulator.hpp"
#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief Simulator data storage class
      \author Joel Andersson
      \date 2010
  */
  class CASADI_CORE_EXPORT SimulatorInternal : public FunctionInternal {
  public:

    /** \brief  Constructor */
    SimulatorInternal(const Integrator& integrator, const Function& output_fcn,
                      const std::vector<double>& grid);

    /** \brief  Destructor */
    virtual ~SimulatorInternal();

    /** \brief  Clone */
    virtual SimulatorInternal* clone() const { return new SimulatorInternal(*this);}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  initialize */
    virtual void init();

    /** \brief  Integrate */
    virtual void evaluate();

    // Integrator instance
    Integrator integrator_;

    // Output function to be evaluated at each grid point
    Function output_fcn_;

    // Time grid
    std::vector<double> grid_;

    // Iterators to current outputs
    std::vector<std::vector<double>::iterator> output_its_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_SIMULATOR_INTERNAL_HPP
