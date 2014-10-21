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

    /// Initialize
    void init();

    /// Copy the necessary options to the target solver
    void setTargetOptions();

    // Solver name of target solver
    static std::string solvername() { return Solver::shortname() + "_solver"; }

    // Options name of target solver
    static std::string optionsname() { return solvername() + "_options"; }

    /// \brief Get the name of the target solver
    std::string targetName();

    /// Add options that are common to all Adaptor classes
    void addOptions();
  };

  template< class Derived, class Solver>
  void Adaptor<Derived, Solver>::addOptions() {
    Derived* this_ = static_cast<Derived*>(this);

    // TODO(@jgillis): Fix option descriptions
    this_->addOption("target",          OT_BOOLEAN,     false,
                     "Options to be passed to the target solver.");
    this_->addOption("target_options",  OT_DICTIONARY,  GenericType(),
                     "Options to be passed to the target solver.");
    this_->addOption(solvername(),      OT_STRING,      GenericType(),
                     "User-defined DPLE solver class.");
    this_->addOption(optionsname(),     OT_DICTIONARY,  GenericType(),
                     "Options to be passed to the DPLE solver.");
  }

  template< class Derived, class Solver>
  void Adaptor<Derived, Solver>::init() {
    Derived* this_ = static_cast<Derived*>(this);

    if (this_->getOption("target")) {
      if (this_->hasSetOption("target_options")) {
        this_->setOption(this_->getOption("target_options"));
      }
    }
  }

  template< class Derived, class Solver>
  void Adaptor<Derived, Solver>::setTargetOptions() {
    Derived* this_ = static_cast<Derived*>(this);

    if (this_->hasSetOption(optionsname())) {
      this_->solver_.setOption(this_->getOption(optionsname()));
    }
    if (!this_->getOption("target")) {
      if (this_->solver_.hasOption("target_options")) {
        if (this_->hasSetOption("target_options")) {
          this_->solver_.setOption("target_options", this_->getOption("target_options"));
        }
      } else {
        if (this_->hasSetOption("target_options")) {
          this_->solver_.setOption(this_->getOption("target_options"));
        }
      }
    }
  }

  template< class Derived, class Solver>
  std::string Adaptor<Derived, Solver>::targetName() {
    return static_cast<Derived*>(this)->getOption(solvername());
  }

} // namespace casadi

/// \endcond

#endif // CASADI_ADAPTOR_HPP
