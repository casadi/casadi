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

#ifndef VERTSPLIT_HPP
#define VERTSPLIT_HPP

#include "multiple_output.hpp"
#include <map>
#include <stack>

namespace CasADi{
  /** \brief Vertical split, x0,x1,...,[y0,y1,...] -> x0+y0, x1+y1,...
      \author Joel Andersson
      \date 2013
  */
  class Vertsplit : public MultipleOutput{
  public:

    /// Constructor
    Vertsplit(const std::vector<MX>& x, const MX& y);

    /** \brief  Number of outputs */
    virtual int getNumOutputs() const{ return ndep()-1; }
        
    /** \brief  Get the sparsity of output oind */
    virtual const CRSSparsity& sparsity(int oind) const{ return dep(oind).sparsity();}

    /// Clone function
    virtual Vertsplit* clone() const;
      
    /// Destructor
    virtual ~Vertsplit(){}
    
    /// Evaluate the function numerically
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (SX)
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /// Propagate sparsity
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;
    
    /// Evaluate the function (template)
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);    
    
    /** \brief Get the operation */
    virtual int getOp() const{ return OP_VERTSPLIT;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    //virtual int numInplace() const{ return getNumOutputs();}
  };

} // namespace CasADi

#endif // VERTSPLIT_HPP
