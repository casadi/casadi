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

#ifndef INVERSE_MAPPING_HPP
#define INVERSE_MAPPING_HPP

#include "mx_node.hpp"
#include "multiple_output.hpp"
#include <map>

namespace CasADi{
  /** \brief Maps non-zero elements
  \author Joel Andersson
  \date 2011
*/
class InverseMapping : public MultipleOutput{
  public:

    /// Constructor
    InverseMapping(const MX& dep, const std::vector<CRSSparsity>& sp, const std::vector<int>& nzind, const std::vector<int>& depind);

    /// Clone function
    virtual InverseMapping* clone() const;
      
    /// Destructor
    virtual ~InverseMapping(){}
    
    /// Evaluate the function numerically
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (SX)
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /// Propagate sparsity
    virtual void propagateSparsity(const DMatrixPtrV& input, DMatrixPtrV& output);

    /** \brief  Number of outputs */
    virtual int getNumOutputs() const;
    
    /** \brief  Get the sparsity of output oind */
    virtual const CRSSparsity& sparsity(int oind);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /// Sparsity of the outputs
    std::vector<CRSSparsity> sp_;
    
    /// Mapping from the output non-zero to the dependency nonzero index
    std::vector<int> nzind_;

    /// Mapping from the output non-zero index of the dependency index
    std::vector<int> depind_;
    
};

} // namespace CasADi

#endif // INVERSE_MAPPING_HPP
