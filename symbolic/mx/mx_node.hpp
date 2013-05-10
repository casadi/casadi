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
#include "../sx/sx.hpp"
#include "../casadi_math.hpp"
#include "../fx/code_generator.hpp"
#include <vector>
#include <stack>

namespace CasADi{
  //@{
  /** \brief Convenience function, convert vectors to vectors of pointers */
  template<class T>
  std::vector<T*> ptrVec(std::vector<T>& v){
    std::vector<T*> ret(v.size());
    for(int i=0; i<v.size(); ++i) 
      ret[i] = &v[i];
    return ret;
  }

  template<class T>
  const std::vector<T*> ptrVec(const std::vector<T>& v){
    std::vector<T*> ret(v.size());
    for(int i=0; i<v.size(); ++i) 
      ret[i] = const_cast<T*>(&v[i]);
    return ret;
  }
  
  template<class T>
  std::vector<std::vector<T*> > ptrVec(std::vector<std::vector<T> >& v){
    std::vector<std::vector<T*> > ret(v.size());
    for(int i=0; i<v.size(); ++i) 
      ret[i] = ptrVec(v[i]);
    return ret;
  }
  
  template<class T>
  const std::vector<std::vector<T*> > ptrVec(const std::vector<std::vector<T> >& v){
    std::vector<std::vector<T*> > ret(v.size());
    for(int i=0; i<v.size(); ++i) 
      ret[i] = ptrVec(v[i]);
    return ret;
  }
  //@}

  
  /** \brief Node class for MX objects
      \author Joel Andersson 
      \date 2010
      Internal class.
  */
  class MXNode : public SharedObjectNode{
    friend class MX;
  
  public:
    /// Constructor
    MXNode();
  
    /** \brief  Destructor */
    virtual ~MXNode()=0;

    /** \brief  Clone function */
    virtual MXNode* clone() const = 0;

    /** \brief Check the truth value of this node
     */
    virtual bool __nonzero__() const;
    
    /** \brief Check if identically zero */
    virtual bool isZero() const{ return false;}

    /** \brief Check if identically one */
    virtual bool isOne() const{ return false;}

    /** \brief Check if a certain value */
    virtual bool isValue(double val) const{ return false;}

    /** \brief Check if identity matrix */
    virtual bool isIdentity() const{ return false;}

    /** \brief Check if unary operation */
    virtual bool isUnaryOp() const { return false;}

    /** \brief Check if binary operation */
    virtual bool isBinaryOp() const { return false;}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);
    
    /** \brief  Print a representation */
    virtual void repr(std::ostream &stream) const;
    
    /** \brief  Print a description */
    virtual void print(std::ostream &stream) const;
    
    /** \brief  Print expression (make sure number of calls is not exceeded) */
    virtual void print(std::ostream &stream, long& remaining_calls) const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const = 0;

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;
    
    /** \brief  Evaluate the function */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, 
                           const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, 
                           const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens, 
                           std::vector<int>& itmp, std::vector<double>& rtmp){evaluateD(input,output,fwdSeed,fwdSens,adjSeed,adjSens);}

    /** \brief  Evaluate the function, no derivatives*/
    void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp);

    /** \brief  Evaluate symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, 
                            const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, 
                            const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens, 
                            std::vector<int>& itmp, std::vector<SX>& rtmp){ evaluateSX(input,output,fwdSeed,fwdSens,adjSeed,adjSens);}

    /** \brief  Evaluate symbolically (SX), no derivatives */
    void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, 
                    std::vector<int>& itmp, std::vector<SX>& rtmp);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, 
                            const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, 
                            const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given)=0;

    /** \brief  Evaluate symbolically (MX), no derivatives */
    void evaluateMX(const MXPtrV& input, MXPtrV& output);
    
    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp, bool fwd){ propagateSparsity(input,output,fwd);}

    /** \brief  Get the name */
    virtual const std::string& getName() const;
    
    /** \brief  Check if evaluation output */
    virtual bool isOutputNode() const{return false;}

    /** \brief  Check if a multiple output node */
    virtual bool isMultipleOutput() const{return false;}

    /** \brief  Get function reference */
    virtual FX& getFunction();

    /** \brief  Get function reference */
    virtual const FX& getFunction() const{ return const_cast<MXNode*>(this)->getFunction();}

    /** \brief  Get function input */
    virtual int getFunctionInput() const;

    /** \brief  Get function output */
    virtual int getFunctionOutput() const;

    /** \brief Get the operation */
    virtual int getOp() const = 0;

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool isEqual(const MXNode* node, int depth) const{ return false;}
    
    /** \brief Get equality checking depth */
    inline static bool maxDepth(){ return MX::getEqualityCheckingDepth();}

    /** \brief Checks if two nodes have the same operation and have equivalent dependencies up to a given depth */
    bool sameOpAndDeps(const MXNode* node, int depth) const;

    /** \brief  dependencies - functions that have to be evaluated before this one */
    const MX& dep(int ind=0) const;
    MX& dep(int ind=0);
    
    /** \brief  Number of dependencies */
    int ndep() const;
    
    /** \brief  Does the node depend on other nodes*/
    virtual bool hasDep() const{return ndep()>0; }
    
    /** \brief  Number of outputs */
    virtual int getNumOutputs() const{ return 1;}
    
    /** \brief  Get an output */
    virtual MX getOutput(int oind) const;

    /// Get the sparsity
    const CRSSparsity& sparsity() const;

    /// Get the sparsity of output oind
    virtual const CRSSparsity& sparsity(int oind) const;
    
    /** \brief Is the node nonlinear */
    virtual bool isNonLinear(){return false;}
    
    /// Set the sparsity
    void setSparsity(const CRSSparsity& sparsity);
    
    /// Get number of temporary variables needed
    virtual void nTmp(size_t& ni, size_t& nr){ ni=0; nr=0;}

    /// Set unary dependency
    void setDependencies(const MX& dep);
    
    /// Set binary dependencies
    void setDependencies(const MX& dep1, const MX& dep2);
    
    /// Set ternary dependencies
    void setDependencies(const MX& dep1, const MX& dep2, const MX& dep3);
    
    /// Set multiple dependencies
    void setDependencies(const std::vector<MX>& dep);
        
    /// Add a dependency
    int addDependency(const MX& dep);
    
    /// Assign nonzeros (mapping matrix)
    virtual void assign(const MX& d, const std::vector<int>& inz, const std::vector<int>& onz, bool add=false);
    
    /// Assign nonzeros (mapping matrix), output indices sequential
    virtual void assign(const MX& d, const std::vector<int>& inz, bool add=false);

    /// Number of elements
    int numel() const;
    
    /// Get size
    int size() const;
    
    /// Get size
    int size1() const;
    
    /// Get size
    int size2() const;

    /// Convert scalar to matrix
    inline static MX toMatrix(const MX& x, const CRSSparsity& sp){
      if(x.shape()==sp.shape()){
        return x;
      } else {
        return MX(sp,x);
      }
    }

    /// Get the value (only for scalar constant nodes)
    virtual double getValue() const;
    
    /// Get the value (only for constant nodes)
    virtual Matrix<double> getMatrixValue() const;
    
    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual int numInplace() const{ return 0;}

    /// Convert vector of pointers to vector of objects
    template<typename T>
    static std::vector<T> getVector(const std::vector<T*> v);

    /// Convert vector of vectors of pointers to vector of vectors of objects
    template<typename T>
    static std::vector<std::vector<T> > getVector(const std::vector<std::vector<T*> > v);

    /// Simplify the expression (ex is a reference to the node)
    virtual void simplifyMe(MX& ex){}

    /// Get an IMatrix representation of a GetNonzeros or SetNonzeros node
    virtual Matrix<int> mapping() const;

    /// Transpose
    virtual MX getTranspose() const;

    /// Reshape
    virtual MX getReshape(const CRSSparsity& sp) const;
    
    /// Matrix multiplcation
    virtual MX getMultiplication(const MX& y) const;

    /// Solve for square linear system
    virtual MX getSolve(const MX& r, bool tr) const;

    /// Get the nonzeros of matrix
    virtual MX getGetNonzeros(const CRSSparsity& sp, const std::vector<int>& nz) const;

    /// Assign the nonzeros of a matrix to another matrix
    virtual MX getSetNonzeros(const MX& y, const std::vector<int>& nz) const;

    /// Add the nonzeros of a matrix to another matrix
    virtual MX getAddNonzeros(const MX& y, const std::vector<int>& nz) const;

    /// Get submatrix reference
    virtual MX getSubRef(const Slice& i, const Slice& j) const;    

    /// Get submatrix assignment
    virtual MX getSubAssign(const MX& y, const Slice& i, const Slice& j) const;    

    /// Create set sparse
    virtual MX getSetSparse(const CRSSparsity& sp) const;
    
    /// Get a unary operation
    virtual MX getUnary(int op) const;

    /// Get a binary operation operation
    MX getBinarySwitch(int op, const MX& y) const;

    /// Get a binary operation operation (matrix-matrix)
    virtual MX getBinary(int op, const MX& y, bool scX, bool scY) const;

    /// Determinant
    virtual MX getDeterminant() const;

    /// Inverse
    virtual MX getInverse() const;

    /// Inner product
    virtual MX getInnerProd(const MX& y) const;

    /** Temporary variables to be used in user algorithms like sorting, 
        the user is resposible of making sure that use is thread-safe
        The variable is initialized to zero
    */
    int temp;
    
    /** \brief  dependencies - functions that have to be evaluated before this one */
    std::vector<MX> dep_;
    
    /** \brief  The sparsity pattern */
    CRSSparsity sparsity_;
    
  protected:

    /** \brief  Evaluate the function (no work)*/
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, 
                           const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, 
                           const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate symbolically (SX), no work */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, 
                            const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, 
                            const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity, no work */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief Free adjoint memory */
    template<typename T> 
    static void clearVector(const std::vector<std::vector<T*> > v);
  };

  // Implementations

  template<typename T>
  std::vector<T> MXNode::getVector(const std::vector<T*> v){
    std::vector<T> ret(v.size());
    for(int i=0; i<v.size(); i++){
      if(v[i]!=0){
        ret[i] = *v[i];
      }
    }
    return ret;
  }

  template<typename T>
  std::vector<std::vector<T> > MXNode::getVector(const std::vector<std::vector<T*> > v){
    std::vector<std::vector<T> > ret(v.size());
    for(int i=0; i<v.size(); i++){
      ret[i] = getVector(v[i]);
    }
    return ret;
  }

  template<typename T>
  void MXNode::clearVector(const std::vector<std::vector<T*> > v){
    for(int i=0; i<v.size(); ++i){
      for(int j=0; j<v[i].size(); ++j){
        if(v[i][j]!= 0){
          v[i][j]->setZero();
        }
      }
    }
  }


} // namespace CasADi


#endif // MX_NODE_HPP
