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


#ifndef CASADI_DENSE_IO_HPP
#define CASADI_DENSE_IO_HPP

#include "../function/function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief A helper class for Functions that work on dense Inputs/Outputs
  *
  * 
  *    \author Joris Gillis
  *    \date 2014
  */
  template< class Derived>
  class DenseIO  {
    public:
      /// Initialize
      void init();
      /// Read the sparse inputs into the dense inputs
      void readInputs();
      /// Write the dense outputs to the sparse outputs
      void writeOutputs();

      DMatrix & inputD(int i);
      DMatrix & outputD(int i);
      const DMatrix & inputD(int i) const;
      const DMatrix & outputD(int i) const;
    protected:
      /// The dense io interface
      std::vector<DMatrix> dense_inputs_;
      /// The dense io interface
      std::vector<DMatrix> dense_outputs_;
  };

#ifndef SWIG
template< class Derived>
void DenseIO<Derived>::init() {

  Derived* d = static_cast<Derived*>(this);

  dense_inputs_.resize(d->getNumInputs());
  dense_outputs_.resize(d->getNumOutputs());

  for (int i=0;i< d->getNumInputs();++i) {
    if (!d->input(i).isDense()) {
      dense_inputs_[i] = dense(d->input(i));
    }
  }
  for (int i=0;i< d->getNumOutputs();++i) {
    if (!d->output(i).isDense()) {
      dense_outputs_[i] = dense(d->output(i));
    }
  }

}

template< class Derived>
DMatrix & DenseIO<Derived>::inputD(int i) {
  Derived* d = static_cast<Derived*>(this);

  if (d->input(i).isDense()) {
    return d->input(i);
  } else {
    return dense_inputs_[i];
  }
}

template< class Derived>
const DMatrix & DenseIO<Derived>::inputD(int i) const {
  Derived* d = static_cast<Derived*>(this);

  if (d->input(i).isDense()) {
    return d->input(i);
  } else {
    return dense_inputs_[i];
  }
}

template< class Derived>
DMatrix & DenseIO<Derived>::outputD(int i) {
  Derived* d = static_cast<Derived*>(this);

  if (d->output(i).isDense()) {
    return d->output(i);
  } else {
    return dense_outputs_[i];
  }
}

template< class Derived>
const DMatrix & DenseIO<Derived>::outputD(int i) const {
  Derived* d = static_cast<Derived*>(this);

  if (d->output(i).isDense()) {
    return d->output(i);
  } else {
    return dense_outputs_[i];
  }
}

template< class Derived>
void DenseIO<Derived>::readInputs() {

  Derived* d = static_cast<Derived*>(this);

  for (int i=0;i< d->getNumInputs();++i) {
    if (!d->input(i).isDense()) {
      inputD(i).set(d->input(i));
    }
  }

}

template< class Derived>
void DenseIO<Derived>::writeOutputs() {

  Derived* d = static_cast<Derived*>(this);

  for (int i=0;i< d->getNumOutputs();++i) {
    if (!d->output(i).isDense()) {
      d->setOutput(outputD(i), i);
    }
  }

}

#endif
} // namespace casadi

/// \endcond

#endif // CASADI_DENSE_IO_HPP
