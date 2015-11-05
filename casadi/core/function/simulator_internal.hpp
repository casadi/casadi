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

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief Simulator data storage class
      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT SimulatorInternal : public FunctionInternal {
  public:

    /** \brief  Constructor */
    SimulatorInternal(const std::string& name, const Function& integrator);

    /** \brief  Destructor */
    virtual ~SimulatorInternal();

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() const { return IVPSOL_NUM_IN;}
    virtual size_t get_n_out() const { return IVPSOL_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int ind) const {
      return integrator_.sparsity_in(ind);
    }
    virtual Sparsity get_sparsity_out(int ind) const {
      return repmat(integrator_.sparsity_out(ind), 1, grid_.size());
    }
    /// @}

    /** \brief  initialize */
    virtual void init();

    /** \brief  Integrate */
    virtual void evalD(void* mem, const double** arg, double** res, int* iw, double* w);

    // Ivpsol instance
    Function integrator_;

    // Time grid
    std::vector<double> grid_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_SIMULATOR_INTERNAL_HPP
