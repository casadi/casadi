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

#ifndef MX_NODE_HPP
#define MX_NODE_HPP

#include "mx.hpp"
#include <vector>

namespace CasADi{
  typedef double* Dptr;
  typedef std::vector<double*> VDptr;
  typedef std::vector<std::vector<double* > > VVDptr;
  
  
/** \brief Node class for MX objects
    \author Joel Andersson 
    \date 2010
    Internal class.
*/
class MXNode : public SharedObjectNode{
  friend class MX;
  friend class MXFunctionInternal;
  
  public:
    /// Constructor
    MXNode();
  
    /** \brief  Destructor: better _not_ to make this destructor virtual? In the SX class, a virtual destructor leads to stack overflow due to recursive calling */
    virtual ~MXNode();

    /** \brief  Clone function */
    virtual MXNode* clone() const = 0;

    /** \brief  Print description */
    virtual void print(std::ostream &stream) const;

    /** \brief  Evaluate the function */
    virtual void evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj) = 0;
                          
    /** \brief  Get the name */
    virtual const std::string& getName() const;
    
    /** \brief  Check if symbolic */
    virtual bool isSymbolic() const;

    /** \brief  Check if constant */
    virtual bool isConstant() const;

    /** \brief  dependencies - functions that have to be evaluated before this one */
    const MX& dep(int ind=0) const;
    
    /** \brief  Number of dependencies */
    int ndep() const;

    /// Get the sparsity
    const CRSSparsity& sparsity() const;

  protected:
    
    /// Set the sparsity
    void setSparsity(const CRSSparsity& sparsity);
    
    /// Set unary dependency
    void setDependencies(const MX& dep);
    
    /// Set binary dependencies
    void setDependencies(const MX& dep1, const MX& dep2);
    
    /// Set ternary dependencies
    void setDependencies(const MX& dep1, const MX& dep2, const MX& dep3);
    
    /// Set multiple dependencies
    void setDependencies(const std::vector<MX>& dep);
        
    /// Number of elements
    int numel() const;
    
    /// Get size
    int size() const;
    
    /// Get size
    int size1() const;
    
    /// Get size
    int size2() const;
    
    /** \brief  dependencies - functions that have to be evaluated before this one */
    std::vector<MX> dep_;
    
    /** \brief  The sparsity pattern */
    CRSSparsity sparsity_;
};

} // namespace CasADi


#endif // MX_NODE_HPP
