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
#include <stack>

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
    virtual Mapping* clone() const;
      
    /// Destructor
    virtual ~Mapping(){}
    
    /// Evaluate the function numerically
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (SX)
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /// Propagate sparsity
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;
    
    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;

    /// Assign/add nonzeros
    virtual void assign(const MX& d, const std::vector<int>& inz, const std::vector<int>& onz, bool add=false);

    /// Assign/add nonzeros, outputs sequential
    virtual void assign(const MX& d, const std::vector<int>& inz, bool add=false);

    /// Initialize
    virtual void init();

    /// Map to locate the dependencies
    std::map<const MXNode*, int> depmap_;

    /// Input nonzero and dependency index
    struct OutputNZ{
      int inz, iind;
    };
    
    /// Operation sequence
    typedef std::vector<std::pair<int,int> > IOMap;
    
    /* Operations sorted by output nonzero - always available
    *  The outer vector is size size()
    *  The inner vector lists elements to be summed
    */
    std::vector<std::vector<OutputNZ> > output_sorted_;
    
    /// Operations sorted by input and output index and then by output nonzero (this is the runtime)
    std::vector<std::vector<IOMap> > index_output_sorted_;

    /// Evaluate a block given the data vectors
    template<typename T>
    void evaluateBlock(int iind, int oind, const std::vector<T>& idata, std::vector<T>& odata, bool fwd) const;

    /// Evaluate the function (template)
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);
    
    /// Compare output index
    static bool outputSmaller(const std::pair<int,int>& el1, const std::pair<int,int>& el2){return el1.second<el2.second;}
    
    /// Check equality for output index
    static bool outputEqual(const std::pair<int,int>& el1, const std::pair<int,int>& el2){return el1.second==el2.second;}
    
    /// Construct the IMatrix that maps from the iind'th input to the output 
    Matrix<int> mapping(int iind=0) const;
    
    /// Get mapping from the output non-zero index of the dependency index
    std::vector<int> getDepInd() const;
    
    /// Check if the mapping is in fact an identity mapping
    bool isIdentity() const;
    
    /** \brief Get the operation */
    virtual int getOp() const{ return OP_MAPPING;}    
};

} // namespace CasADi

#endif // MAPPING_HPP
