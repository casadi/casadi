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


#ifndef CASADI_PARALLELIZER_INTERNAL_HPP
#define CASADI_PARALLELIZER_INTERNAL_HPP

#include <vector>
#include "parallelizer.hpp"
#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief  Internal node class for Parallelizer
      \author Joel Andersson
      \date 2010
  */
  class CASADI_CORE_EXPORT ParallelizerInternal : public FunctionInternal {
    friend class Parallelizer;

  protected:
    /// Constructor
    explicit ParallelizerInternal(const std::vector<Function>& funcs);

  public:
    /// clone
    virtual ParallelizerInternal* clone() const {
      ParallelizerInternal* ret = new ParallelizerInternal(*this);
      for (std::vector<Function>::iterator it=ret->funcs_.begin(); it!=ret->funcs_.end(); ++it) {
        it->makeUnique();
      }
      return ret;
    }

    /// Destructor
    virtual ~ParallelizerInternal();

    /// Evaluate the all the tasks
    virtual void evaluate();

    /// Evaluate a single task
    virtual void evaluateTask(int task);

    /// Reset the sparsity propagation
    virtual void spInit(bool use_fwd);

    /// Propagate the sparsity pattern through a set of directional derivatives forward or backward
    virtual void spEvaluate(bool use_fwd);

    /** \brief Propagate the sparsity pattern through a set of directional derivatives
     * forward or backward, one task only
     */
    void spEvaluateTask(bool use_fwd, int task);

    /// Is the class able to propagate seeds through the algorithm?
    virtual bool spCanEvaluate(bool fwd) { return true;}

    /// Generate a function that calculates nfwd forward derivatives and nadj adjoint derivatives
    virtual Function getDerivative(int nfwd, int nadj);

    /// Generate a function that calculates a Jacobian function
    virtual Function getJacobian(int iind, int oind, bool compact, bool symmetric);

    /// Initialize
    virtual void init();

    /// Generate the sparsity of a Jacobian block
    virtual Sparsity getJacSparsity(int iind, int oind, bool symmetric);

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /// Functions
    std::vector<Function> funcs_;

    /// Input argument indices
    std::vector<int> inind_;

    /// Output argument indices
    std::vector<int> outind_;

    /// Is a function a copy of another
    std::vector<int> copy_of_;

    /// Parallelization modes
    enum Mode {SERIAL, OPENMP, MPI};

    /// Mode
    Mode mode_;
  };



} // namespace casadi

/// \endcond
#endif // CASADI_PARALLELIZER_INTERNAL_HPP

