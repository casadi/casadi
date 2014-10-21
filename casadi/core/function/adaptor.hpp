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


#ifndef CASADI_ADAPTOR_HPP
#define CASADI_ADAPTOR_HPP

#include "../function/function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief A helper class for a Solver that delegates work to another solver
   *
   *
   *    \author Joris Gillis
   *    \date 2014
   */
  template< class Derived, class Solver>
  class Adaptor {
  public:

    // Solver name of target solver
    static std::string solvername() { return Solver::shortname() + "_solver"; }

    // Options name of target solver
    static std::string optionsname() { return solvername() + "_options"; }

    /// Add options that are common to all Adaptor classes
    void addOptions();
  };

  template< class Derived, class Solver>
  void Adaptor<Derived, Solver>::addOptions() {
    Derived* this_ = static_cast<Derived*>(this);
    this_->addOption(solvername(),      OT_STRING,      GenericType(),
                     "Name of solver.");
    this_->addOption(optionsname(),     OT_DICTIONARY,  GenericType(),
                     "Options to be passed to solver.");
  }

} // namespace casadi

/// \endcond

#endif // CASADI_ADAPTOR_HPP
