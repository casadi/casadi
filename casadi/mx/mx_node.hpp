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

/* Developer nodes:
  
  This MX classes need to get the ability to generate whole Jacobian matrices, not just directional derivatives. There exist two principle ways of
  doing this:
  
  1. The current approach:
     Generate the Jacobian by seeding in all directions, using curtis powel reed scaling, the number of directions can be reduced.
     - Problems:
       * Need a way to calculate the Jacobian sparsity pattern, detecting it completely black box is expensiveunsafe and potentially unsafe.
       * No easy way to generate second and higher order sensitivities
       
  2. The source code approach:
     Generate a new tree for the derivatives, as for the SXFunction class.
     - Advantages:
      * Arbitrary order derivatives
      * Easier to parallelize efficiently
     - Problems:
      * More implementation than approach 1.
*/


namespace CasADi{
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

    /** \brief  Print expression */
    virtual void print(std::ostream &stream, const std::vector<std::string>& args) const=0;
    
    /** \brief  Print expression */
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
    
    /** \brief Generate all the partial derivatives for the node - needed for symbolic AD algorithm */
    virtual std::vector<MX> partial() const;

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
