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

#ifndef NONLINEAR_SOLVE_HPP
#define NONLINEAR_SOLVE_HPP

#include "mx_node.hpp"
#include "../fx/implicit_function.hpp"

namespace CasADi{
  /** \brief An MX atomic for a nonlinear solve
      \author Joel Andersson
      \date 2013
  */
  class NonlinearSolve : public MXNode{
  public:
    
    /** \brief  Constructor */
    NonlinearSolve(const std::vector<MX>& x, const ImplicitFunction& implicit_function);

    /** \brief  Destructor */
    virtual ~NonlinearSolve(){}

    /** \brief  Clone function */
    virtual NonlinearSolve* clone() const{ return new NonlinearSolve(*this);}

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;

    /// Evaluate the function numerically
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    //    /** \brief  Propagate sparsity */
    //    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);
    
    /** \brief Get the operation */
    virtual int getOp() const{ return OP_NONLINEAR_SOLVE;}

    /** \brief  Get function reference */
    virtual FX& getFunction(){ return implicit_function_;}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /// Implicit function (may be shared between multiple nodes)
    ImplicitFunction implicit_function_;
  };


} // namespace CasADi


#endif // NONLINEAR_SOLVE_HPP
