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
    virtual Mapping* clone() const;
      
    /// Destructor
    virtual ~Mapping(){}
    
    /// Evaluate the function numerically
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (SX)
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /// Propagate sparsity
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;
    
    /// Is a mapping matrix
    virtual bool isMapping() const{return true;}

    /// Assign nonzeros (mapping matrix)
    virtual void assign(const MX& d, const IOMap& iomap);

    /// Assign nonzeros (index given) -> change to multiple indices?
    virtual void assignIndex(int depind, const IOMap& iomap);

    /// Initialize
    virtual void init();

    /// Check if the mapping is ready
    bool isReady() const;
    
    /// Mapping from the output non-zero to the dependency nonzero index
    std::vector<int> nzind_;

    /// Mapping from the output non-zero index of the dependency index
    std::vector<int> depind_;

    /// Map to locate the dependencies
    std::map<const MXNode*, int> depmap_;

    /// Assignment operations (runtime)
    std::vector<std::vector<IOMap> > assignments_;

    /// Addition operations (runtime)
    std::vector<std::vector<IOMap> > additions_;
    
    /// Unsorted operations
    struct Unsorted{
      // Constructor
      Unsorted(int inz__=-1, int onz__=-1, int iind__=0, bool is_add__=false)
      : inz(inz__), onz(onz__), iind(iind__), is_add(is_add__){}
      
      // Sorting by output nonzero
      bool operator<(const Unsorted& other){ return onz<other.onz;}
      
      // Input/output nonzero and dependency index
      int inz, onz, iind;
      
      // Is the operation an addition?
      bool is_add;
    };
    
    /// Not yet sorted operations
    std::vector<Unsorted> unsorted_;
    
    /// Evaluate a block given the data vectors
    void evaluateBlock(int iind, int oind, const std::vector<double>& idata, std::vector<double>& odata, bool fwd) const;
    
    /// Compare output index
    static bool outputSmaller(const std::pair<int,int>& el1, const std::pair<int,int>& el2){return el1.second<el2.second;}
    
    /// Check equality for output index
    static bool outputEqual(const std::pair<int,int>& el1, const std::pair<int,int>& el2){return el1.second==el2.second;}
};

} // namespace CasADi

#endif // MAPPING_HPP
