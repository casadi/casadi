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
    /// Destructor
    virtual ~ParallelizerInternal();
    
    /// Evaluate the all the tasks
    virtual void evaluate(int nfdir, int nadir);

    /// Evaluate a single task
    virtual void evaluateTask(int task, int nfdir, int nadir);

    /// Initialize
    virtual void init();

    /// Functions
    std::vector<FX> funcs_;
    
    /// Input argument indices
    std::vector<int> inind_;
    
    /// Output argument indices
    std::vector<int> outind_;
    
    /// Parallelization modes
    enum Mode{SERIAL,OPENMP,MPI};
    
    /// Mode
    Mode mode_;
    
    /// Save corrected input values after evaluation
    bool save_corrected_input_;
};



} // namespace CasADi


#endif // PARALLELIZER_INTERNAL_HPP

