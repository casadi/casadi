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

#ifndef MX_FUNCTION_INTERNAL_HPP
#define MX_FUNCTION_INTERNAL_HPP

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "mx_function.hpp"
#include "x_function_internal.hpp"
#include "../mx/mx_node.hpp"

namespace CasADi{

  /** \brief  Internal node class for MXFunction
      \author Joel Andersson 
      \date 2010
  */
  class MXFunctionInternal : public XFunctionInternal<MXFunction,MXFunctionInternal,MX,MXNode>{
    friend class MXFunction;
  
  public:

    /** \brief  Multiple input, multiple output constructor, only to be accessed from MXFunction, therefore protected */
    MXFunctionInternal(const std::vector<MX>& input, const std::vector<MX>& output);

    /** \brief  Make a deep copy */
    virtual MXFunctionInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);
    
    /** \brief  Destructor */
    virtual ~MXFunctionInternal();

    /** \brief  Evaluate the algorithm */
    virtual void evaluate(int nfdir, int nadir);

    /** \brief  Print description */
    virtual void print(std::ostream &stream) const;

    /** \brief  Initialize */
    virtual void init();

    /** \brief  Update the number of sensitivity directions during or after initialization */
    virtual void updateNumSens(bool recursive);

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(std::ostream &stream, const std::string& type, CodeGenerator& gen) const;
    
    /** \brief Generate code for the body of the C function */
    virtual void generateBody(std::ostream &stream, const std::string& type, CodeGenerator& gen) const;

    /** \brief Extract the residual function G and the modified function Z out of an expression (see Albersmeyer2010 paper) */
    void generateLiftingFunctions(MXFunction& vdef_fcn, MXFunction& vinit_fcn);

    /** \brief Generate a function that calculates a Jacobian function by operator overloading */
    virtual FX getNumericJacobian(int iind, int oind, bool compact, bool symmetric);
    
    /** \brief  An elemenent of the algorithm, namely an MX node */
    typedef MXAlgEl AlgEl;

    /** \brief  All the runtime elements in the order of evaluation */
    std::vector<AlgEl> algorithm_;

    /** \brief  Working vector for numeric calculation */
    std::vector<FunctionIO> work_;
  
    /** \brief  Temporary vectors needed for the evaluation (integer) */
    std::vector<int> itmp_;
    
    /** \brief  Temporary vectors needed for the evaluation (real) */
    std::vector<double> rtmp_;

    /** \brief  "Tape" with spilled variables */
    std::vector<std::pair<std::pair<int,int>,DMatrix> > tape_;
    
    /// Free variables
    std::vector<MX> free_vars_;
        
    /** \brief Evaluate symbolically, SX type*/
    virtual void evalSXsparse(const std::vector<SXMatrix>& input, std::vector<SXMatrix>& output, 
                        const std::vector<std::vector<SXMatrix> >& fwdSeed, std::vector<std::vector<SXMatrix> >& fwdSens, 
                        const std::vector<std::vector<SXMatrix> >& adjSeed, std::vector<std::vector<SXMatrix> >& adjSens);
                        
    /** \brief Evaluate symbolically, MX type */
    virtual void evalMX(const std::vector<MX>& input, std::vector<MX>& output, 
                        const std::vector<std::vector<MX> >& fwdSeed, std::vector<std::vector<MX> >& fwdSens, 
                        const std::vector<std::vector<MX> >& adjSeed, std::vector<std::vector<MX> >& adjSens);

    /** \brief Expand the matrix valued graph into a scalar valued graph */
    SXFunction expand(const std::vector<SXMatrix>& inputv );
    
    // Update pointers to a particular element
    void updatePointers(const AlgEl& el, int nfdir, int nadir);
    
    // Vectors to hold pointers during evaluation
    DMatrixPtrV mx_input_;
    DMatrixPtrV mx_output_;
    DMatrixPtrVV mx_fwdSeed_;
    DMatrixPtrVV mx_fwdSens_;
    DMatrixPtrVV mx_adjSeed_;
    DMatrixPtrVV mx_adjSens_;

    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<MX> symbolicInput() const{ return inputv_;}

    /// Propagate a sparsity pattern through the algorithm
    virtual void spEvaluate(bool fwd);

    /// Is the class able to propate seeds through the algorithm?
    virtual bool spCanEvaluate(bool fwd){ return true;}

    /// Reset the sparsity propagation
    virtual void spInit(bool fwd);
    
    /// Print work vector
    void printWork(int nfdir=0, int nadir=0, std::ostream &stream=std::cout);
    
    /// Print tape
    void printTape(std::ostream &stream=std::cout);
    
    /// Allocate tape
    void allocTape();
    
    // print an element of an algorithm
    void print(std::ostream &stream, const AlgEl& el) const;
    
  };

} // namespace CasADi


#endif // MX_FUNCTION_INTERNAL_HPP

