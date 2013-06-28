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

#ifndef PARALLELIZER_INTERNAL_HPP
#define PARALLELIZER_INTERNAL_HPP

#include <vector>
#include "parallelizer.hpp"
#include "fx_internal.hpp"

namespace CasADi{
 
  /** \brief  Internal node class for Parallelizer
  \author Joel Andersson 
  \date 2010
*/
class ParallelizerInternal : public FXInternal{
  friend class Parallelizer;
  
  protected:
    /// Constructor
    explicit ParallelizerInternal(const std::vector<FX>& funcs);

  public:
    /// clone
    virtual ParallelizerInternal* clone() const{ 
      ParallelizerInternal* ret = new ParallelizerInternal(*this);
      for(std::vector<FX>::iterator it=ret->funcs_.begin(); it!=ret->funcs_.end(); ++it){
        it->makeUnique();
      }
      return ret;
    }
    
    /// Destructor
    virtual ~ParallelizerInternal();
    
    /// Evaluate the all the tasks
    virtual void evaluate(int nfdir, int nadir);

    /// Evaluate a single task
    virtual void evaluateTask(int task, int nfdir, int nadir);

    /// Reset the sparsity propagation
    virtual void spInit(bool use_fwd);
    
    /// Propagate the sparsity pattern through a set of directional derivatives forward or backward
    virtual void spEvaluate(bool use_fwd);
    
    /// Propagate the sparsity pattern through a set of directional derivatives forward or backward, one task only
    void spEvaluateTask(bool use_fwd, int task);

    /// Is the class able to propate seeds through the algorithm?
    virtual bool spCanEvaluate(bool fwd){ return true;}
    
    /// Generate a function that calculates nfwd forward derivatives and nadj adjoint derivatives
    virtual FX getDerivative(int nfwd, int nadj);
    
    /// Generate a function that calculates a Jacobian function
    virtual FX getJacobian(int iind, int oind, bool compact, bool symmetric);

    /// Initialize
    virtual void init();

    /// Generate the sparsity of a Jacobian block
    virtual CRSSparsity getJacSparsity(int iind, int oind, bool symmetric);

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /** \brief  Update the number of sensitivity directions during or after initialization */
    virtual void updateNumSens(bool recursive);
    
    /// Functions
    std::vector<FX> funcs_;
    
    /// Input argument indices
    std::vector<int> inind_;
    
    /// Output argument indices
    std::vector<int> outind_;
    
    /// Is a function a copy of another
    std::vector<int> copy_of_;
    
    /// Parallelization modes
    enum Mode{SERIAL,OPENMP,MPI};
    
    /// Mode
    Mode mode_;
};



} // namespace CasADi


#endif // PARALLELIZER_INTERNAL_HPP

