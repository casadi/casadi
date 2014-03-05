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

#ifndef SPLIT_HPP
#define SPLIT_HPP

#include "multiple_output.hpp"
#include <map>
#include <stack>

namespace CasADi{

  /** \brief Split: Split into multiple expressions splitting the nonzeros
      \author Joel Andersson
      \date 2014
  */
  class Split : public MultipleOutput{
  public:
    /// Constructor
    Split(const MX& x, const std::vector<int>& offset);

    /// Destructor
    virtual ~Split() = 0;

    /** \brief  Number of outputs */
    virtual int getNumOutputs() const{ return output_sparsity_.size(); }
        
    /** \brief  Get the sparsity of output oind */
    virtual const Sparsity& sparsity(int oind) const{ return output_sparsity_.at(oind);}

    /// Evaluate the function numerically
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp, std::vector<SXElement>& rtmp);

    /// Propagate sparsity
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;

    /// Evaluate the function (template)
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp, std::vector<T>& rtmp);    

    // Sparsity pattern of the outputs
    std::vector<int> offset_;
    std::vector<Sparsity> output_sparsity_;
  };

  /** \brief Horizontal split, x -> x0, x1,...
      \author Joel Andersson
      \date 2013
  */
  class Horzsplit : public Split{
  public:
    
    /// Constructor
    Horzsplit(const MX& x, const std::vector<int>& offset);

    /// Destructor
    virtual ~Horzsplit(){}

    /// Clone function
    virtual Horzsplit* clone() const{ return new Horzsplit(*this);}
          
    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;
        
    /** \brief Get the operation */
    virtual int getOp() const{ return OP_HORZSPLIT;}    
  };

  /** \brief Vertical split of vectors, x -> x0, x1,...
      \author Joel Andersson
      \date 2014
  */
  class Vertsplit : public Split{
  public:
    
    /// Constructor
    Vertsplit(const MX& x, const std::vector<int>& offset);

    /// Destructor
    virtual ~Vertsplit(){}

    /// Clone function
    virtual Vertsplit* clone() const{ return new Vertsplit(*this);}
          
    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;
        
    /** \brief Get the operation */
    virtual int getOp() const{ return OP_VERTSPLIT;}
  };

} // namespace CasADi

#endif // SPLIT_HPP
