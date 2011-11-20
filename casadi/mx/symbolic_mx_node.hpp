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

#ifndef SYMBOLIC_MATRIX_HPP
#define SYMBOLIC_MATRIX_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents a symbolic MX
  \author Joel Andersson 
  \date 2010
  A regular user is not supposed to work with this Node class.
  This user can call MX(name,n,m) directly.
*/
class SymbolicMatrix : public MXNode{
  public:

    /** \brief  Constructors */
    explicit SymbolicMatrix(const std::string& name, int n=1, int m=1);

    /** \brief  Constructors */
    explicit SymbolicMatrix(const std::string& name, const CRSSparsity & sp);
    
    /** \brief  Clone function */
    virtual SymbolicMatrix* clone() const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief  Evaluate the function numerically */
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(const DMatrixPtrV& input, DMatrixPtrV& output);

    /** \brief  Is symbolic */
    virtual bool isSymbolic() const;
    
    /** \brief  Get the name */
    virtual const std::string& getName() const;
    
    /// Symbolic evaluation (matrix graph)
    virtual MX eval(const std::vector<MX>& x){return MX::create(this);}

    /// Partial derivatives
    virtual std::vector<MX> partial(const std::vector<MX>& x);

  protected:
    // Name of the varible
    std::string name_;
};

} // namespace CasADi


#endif // SYMBOLIC_MATRIX_HPP
