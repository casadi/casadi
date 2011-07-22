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

#ifndef JACOBIAN_REFERENCE_HPP
#define JACOBIAN_REFERENCE_HPP

#include "mx_node.hpp"
#include "../fx/fx.hpp"
#include <map>

namespace CasADi{
  /** \brief Maps non-zero elements
  \author Joel Andersson
  \date 2011
*/
class JacobianReference : public MXNode{
  public:

    /// Default constructor
    JacobianReference(const MX& output, int iind);

    /// Clone function
    JacobianReference* clone() const;
      
    /// Destructor
    virtual ~JacobianReference(){}
    
    /// Evaluate the function and store the result in the node
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /// Print
    virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

    /** \brief  Check if jacobian reference */
    virtual bool isJacobian() const{return true;}

    /** \brief  Get function input */
    virtual int getFunctionInput() const{ return iind_;}

    /** \brief  Get function output */
    virtual int getFunctionOutput() const{ return dep(0)->getFunctionOutput();}

    /** \brief  Get function reference */
    virtual FX& getFunction();

    // Input index
    int iind_;
};

} // namespace CasADi

#endif // JACOBIAN_REFERENCE_HPP
