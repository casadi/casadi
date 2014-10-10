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


#ifndef CASADI_WRAPPER_HPP
#define CASADI_WRAPPER_HPP

#include "../function/function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief A helper class for a Function that wrap another Function
  *
  *
  *    \author Joris Gillis
  *    \date 2014
  */
  template< class Derived>
  class Wrapper  {
    public:
      /// Check the dimensions of the internal function after initialization
      void checkDimensions();
      /// Evaluate the internal function and make it external
      void evaluate();
    protected:
      /// The internal function that is being wrapped
      Function f_;
  };

#ifndef SWIG
template< class Derived>
void Wrapper<Derived>::checkDimensions() {

  Derived* d = static_cast<Derived*>(this);

  // Check number of inputs/outputs
  casadi_assert(d->getNumInputs()==f_.getNumInputs());
  casadi_assert(d->getNumOutputs()==f_.getNumOutputs());

  // Check sparsities of inputs/outputs
  for (int i=0;i< d->getNumInputs();++i) {
    casadi_assert_message(d->input(i).sparsity()==f_.input(i).sparsity(),
      "Sparsity mismatch for input " << i << ":" <<
      d->input(i).dimString() << " <-> " << f_.input(i).dimString() << ".");
  }
  for (int i=0;i< d->getNumOutputs();++i) {
    casadi_assert_message(d->output(i).sparsity()==f_.output(i).sparsity(),
      "Sparsity mismatch for output " << i << ":" <<
      d->output(i).dimString() << " <-> " << f_.output(i).dimString() << ".");
  }

}

template< class Derived>
void Wrapper<Derived>::evaluate() {
  Derived* d = static_cast<Derived*>(this);

  // Copy the inputs from external to internal
  for (int i=0;i< d->getNumInputs();++i) {
    std::copy(d->input(i).begin(), d->input(i).end(), f_.input(i).begin());
  }

  // Evaluate the internal function
  f_.evaluate();

  // Copy the outputs from internal to external
  for (int i=0;i< d->getNumOutputs();++i) {
    std::copy(f_.output(i).begin(), f_.output(i).end(), d->output(i).begin());
  }
}
#endif
} // namespace casadi

/// \endcond

#endif // CASADI_WRAPPER_HPP
