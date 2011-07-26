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

#ifndef MULTIPLE_OUTPUT_HPP
#define MULTIPLE_OUTPUT_HPP

#include "mx_node.hpp"
#include "../fx/fx.hpp"
#include <set>

namespace CasADi{

/// Forward declaration
class OutputNode;
  
/** 
  \author Joel Andersson 
  \date 2010
*/
class MultipleOutput : public MXNode{
  friend class OutputNode;
  public:

    /** \brief  Constructor */
    MultipleOutput();
 
    /** \brief  Destructor */
    virtual ~MultipleOutput();
 
    /** \brief  Number of outputs */
    virtual int getNumOutputs() const=0;
    
    /** \brief  Get the sparsity of output oind */
    virtual const CRSSparsity& sparsity(int oind)=0;

    /** \brief  Check if a multiple output node */
    virtual bool isMultipleOutput() const{return true;}

};

class OutputNode : public MXNode{
  public:
  
    /** \brief  Constructor */
    OutputNode(const MX& parent, int oind);

    /** \brief  Destructor */
    virtual ~OutputNode();

    /** \brief  Clone function */
    virtual OutputNode* clone() const{ return new OutputNode(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /** \brief  Print */
    virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

    /** \brief Is the node nonlinear */
    virtual bool isNonLinear(){return true;} 

    /** \brief  Check if evaluation output */
    virtual bool isOutputNode() const{return true;}
    
    /** \brief  Get function input */
    virtual int getFunctionInput() const{ return -1;}

    /** \brief  Get function output */
    virtual int getFunctionOutput() const{ return oind_;}

    /** \brief  Get the jacobian of an function evaluation with respect to the iind-th argument */
    virtual MX jac(int iind);
    
    /** \brief  Output index */
    int oind_;
};








} // namespace CasADi

#endif // MULTIPLE_OUTPUT_HPP
