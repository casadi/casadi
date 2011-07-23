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

#ifndef MAPPING_HPP
#define MAPPING_HPP

#include "mx_node.hpp"
#include <map>

namespace CasADi{
  /** \brief Maps non-zero elements
  \author Joel Andersson
  \date 2011
*/
class Mapping : public MXNode{
  public:

    /// Constructor
    Mapping(const CRSSparsity& sp);

    /// Clone function
    Mapping* clone() const;
      
    /// Destructor
    virtual ~Mapping(){}
    
    /// Evaluate the function numerically
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (SX)
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens);
    
    /// Print
    virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;
    
    /// Is a mapping matrix
    virtual bool isMapping() const{return true;}
    
    /// Add a dependency (index given)
    virtual void addDependency(int depind, const std::vector<int>& nz_d, const std::vector<int>& nz);
    
    /// Add a dependency
    virtual void addDependency(const MX& d, const std::vector<int>& nz_d, const std::vector<int>& nz);
    
    /// Add a dependency
    virtual void addDependency(const MX& d, const std::vector<int>& nz_d);
    
    /// Symbolic forward sensitivities
    virtual MX adFwd(const std::vector<MX>& jx);
    
    /// Check if the mapping is ready
    bool isReady() const;
    
    /// Mapping from the output non-zero to the dependency nonzero index
    Matrix<int> nzmap_;

    /// Mapping from the output non-zero index of the dependency index
    std::vector<int> depind_;

    /// Map to locate the dependencies
    std::map<const MXNode*, int> depmap_;
};

} // namespace CasADi

#endif // MAPPING_HPP
