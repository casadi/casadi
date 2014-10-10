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

  class AdaptorBase {
    public:
      virtual std::string prefix() = 0;

      bool target_reached_;

  };

  /** \brief A helper class for a Solver that delegates work to another solver
  *
  *
  *    \author Joris Gillis
  *    \date 2014
  */
  template< class Derived, class Solver>
  class Adaptor : public virtual AdaptorBase {
    public:

      /// Initialize
      void init();

      /// Copy the necessary options to the target solver
      void setTargetOptions();

      /// Instance of the target solver
      Solver solver_;

      /// Get the prefix of the target solver
      virtual std::string prefix();

      /// \brief Get the name of the target solver
      std::string targetName();

      /// Add options that are common to all Adaptor classes
      void addOptions();

      /// Load the plugin needed by the target solver
      static void adaptorLoader(const std::string& name);

  };

#ifndef SWIG
template< class Derived, class Solver>
void Adaptor<Derived, Solver>::adaptorLoader(const std::string& name) {
  Solver::loadPlugin(name);
}

template< class Derived, class Solver>
void Adaptor<Derived, Solver>::addOptions() {
  Derived* d = static_cast<Derived*>(this);

  d->addOption("target", OT_BOOLEAN, false, "Options to be passed to the DPLE solver.");
  d->addOption("target_options",         OT_DICTIONARY,   GenericType(),
              "Options to be passed to the DPLE solver.");

  d->addOption(prefix() + "_solver",            OT_STRING, GenericType(),
              "User-defined DPLE solver class.");
  d->addOption(prefix() + "_solver_options",    OT_DICTIONARY,   GenericType(),
              "Options to be passed to the DPLE solver.");
}

template< class Derived, class Solver>
void Adaptor<Derived, Solver>::init() {
  Derived* d = static_cast<Derived*>(this);

  if (d->getOption("target")) {
    if (d->hasSetOption("target_options")) {
      d->setOption(d->getOption("target_options"));
    }
  }
}

template< class Derived, class Solver>
void Adaptor<Derived, Solver>::setTargetOptions() {
  Derived* d = static_cast<Derived*>(this);

  if (d->hasSetOption(prefix() + "_solver_options")) {
    solver_.setOption(d->getOption(prefix() + "_solver_options"));
  }
  if (!d->getOption("target")) {
    if (solver_.hasOption("target_options")) {
      if (d->hasSetOption("target_options")) {
        solver_.setOption("target_options", d->getOption("target_options"));
      }
    } else {
      if (d->hasSetOption("target_options")) {
        solver_.setOption(d->getOption("target_options"));
      }
    }
  }
}

template< class Derived, class Solver>
std::string Adaptor<Derived, Solver>::prefix() {
  std::string p = Solver::infix();
  return p.substr(0, p.size()-6); // solver
}

template< class Derived, class Solver>
std::string Adaptor<Derived, Solver>::targetName() {
  return static_cast<Derived*>(this)->getOption(prefix()+"_solver");
}
#endif

} // namespace casadi

/// \endcond

#endif // CASADI_ADAPTOR_HPP
