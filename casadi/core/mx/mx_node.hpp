/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_MX_NODE_HPP
#define CASADI_MX_NODE_HPP

#include "mx.hpp"
#include "../sx/sx_element.hpp"
#include "../casadi_math.hpp"
#include "../function/code_generator.hpp"
#include "../function/linear_solver.hpp"
#include <vector>
#include <stack>

namespace casadi {
  /// \cond INTERNAL
  ///@{
  /** \brief Convenience function, convert vectors to vectors of pointers */
  template<class T>
  std::vector<T*> ptrVec(std::vector<T>& v) {
    std::vector<T*> ret(v.size());
    for (int i=0; i<v.size(); ++i)
      ret[i] = &v[i];
    return ret;
  }

  template<class T>
  const std::vector<T*> ptrVec(const std::vector<T>& v) {
    std::vector<T*> ret(v.size());
    for (int i=0; i<v.size(); ++i)
      ret[i] = const_cast<T*>(&v[i]);
    return ret;
  }

  template<class T>
  std::vector<std::vector<T*> > ptrVec(std::vector<std::vector<T> >& v) {
    std::vector<std::vector<T*> > ret(v.size());
    for (int i=0; i<v.size(); ++i)
      ret[i] = ptrVec(v[i]);
    return ret;
  }

  template<class T>
  const std::vector<std::vector<T*> > ptrVec(const std::vector<std::vector<T> >& v) {
    std::vector<std::vector<T*> > ret(v.size());
    for (int i=0; i<v.size(); ++i)
      ret[i] = ptrVec(v[i]);
    return ret;
  }
  ///@}


  /** \brief Node class for MX objects
      \author Joel Andersson
      \date 2010
      Internal class.
  */
  class CASADI_CORE_EXPORT MXNode : public SharedObjectNode {
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
    virtual bool isZero() const { return false;}

    /** \brief Check if identically one */
    virtual bool isOne() const { return false;}

    /** \brief Check if a certain value */
    virtual bool isValue(double val) const { return false;}

    /** \brief Check if identity matrix */
    virtual bool isIdentity() const { return false;}

    /** \brief Check if unary operation */
    virtual bool isUnaryOp() const { return false;}

    /** \brief Check if binary operation */
    virtual bool isBinaryOp() const { return false;}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Print a representation */
    virtual void repr(std::ostream &stream) const;

    /** \brief  Print a description */
    virtual void print(std::ostream &stream) const;

    /** \brief  Print expression (make sure number of calls is not exceeded) */
    virtual void print(std::ostream &stream, long& remaining_calls) const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const = 0;

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                   const std::vector<std::string>& res, CodeGenerator& gen) const;

    /** \brief  Evaluate numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                           std::vector<double>& rtmp);

    /** \brief  Evaluate symbolically (SX) */
    virtual void evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                            std::vector<SXElement>& rtmp);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output,
                            const MXPtrVV& fwdSeed, MXPtrVV& fwdSens,
                            const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /** \brief  Evaluate symbolically (MX), no derivatives */
    void evaluateMX(const MXPtrV& input, MXPtrV& output);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                                   std::vector<double>& rtmp, bool fwd)
    { propagateSparsity(input, output, fwd);}

    /** \brief  Get the name */
    virtual const std::string& getName() const;

    /** \brief  Check if evaluation output */
    virtual bool isOutputNode() const {return false;}

    /** \brief  Check if a multiple output node */
    virtual bool isMultipleOutput() const {return false;}

    /** \brief  Get function reference */
    virtual Function& getFunction();

    /** \brief  Get function reference */
    virtual const Function& getFunction() const { return const_cast<MXNode*>(this)->getFunction();}

    /** \brief  Get function input */
    virtual int getFunctionInput() const;

    /** \brief  Get function output */
    virtual int getFunctionOutput() const;

    /** \brief Get the operation */
    virtual int getOp() const = 0;

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool isEqual(const MXNode* node, int depth) const { return false;}

    /** \brief Get equality checking depth */
    inline static bool maxDepth() { return MX::getEqualityCheckingDepth();}

    /** \brief Checks if two nodes have the same operation and have
     * equivalent dependencies up to a given depth */
    bool sameOpAndDeps(const MXNode* node, int depth) const;

    /** \brief  dependencies - functions that have to be evaluated before this one */
    const MX& dep(int ind=0) const;
    MX& dep(int ind=0);

    /** \brief  Number of dependencies */
    int ndep() const;

    /** \brief  Does the node depend on other nodes*/
    virtual bool hasDep() const {return ndep()>0; }

    /** \brief  Number of outputs */
    virtual int getNumOutputs() const { return 1;}

    /** \brief  Get an output */
    virtual MX getOutput(int oind) const;

    /// Get the sparsity
    const Sparsity& sparsity() const { return sparsity_;}

    /// Get the sparsity of output oind
    virtual const Sparsity& sparsity(int oind) const;

    /// Get shape
    int numel() const { return sparsity().numel(); }
    int size() const { return sparsity().size(); }
    int size1() const { return sparsity().size1(); }
    int size2() const { return sparsity().size2(); }
    std::pair<int, int> shape() const { return sparsity().shape();}

    /** \brief Is the node nonlinear */
    virtual bool isNonLinear() {return false;}

    /// Set the sparsity
    void setSparsity(const Sparsity& sparsity);

    /// Get number of temporary variables needed
    virtual void nTmp(size_t& ni, size_t& nr) { ni=0; nr=0;}

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
    virtual void assign(const MX& d, const std::vector<int>& inz,
                        const std::vector<int>& onz, bool add=false);

    /// Assign nonzeros (mapping matrix), output indices sequential
    virtual void assign(const MX& d, const std::vector<int>& inz, bool add=false);

    /// Convert scalar to matrix
    inline static MX toMatrix(const MX& x, const Sparsity& sp) {
      if (x.shape()==sp.shape()) {
        return x;
      } else {
        return MX(sp, x);
      }
    }

    /// Get the value (only for scalar constant nodes)
    virtual double getValue() const;

    /// Get the value (only for constant nodes)
    virtual Matrix<double> getMatrixValue() const;

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual int numInplace() const { return 0;}

    /// Convert vector of pointers to vector of objects
    template<typename T>
    static std::vector<T> getVector(const std::vector<T*> v);

    /// Convert vector of vectors of pointers to vector of vectors of objects
    template<typename T>
    static std::vector<std::vector<T> > getVector(const std::vector<std::vector<T*> > v);

    /// Simplify the expression (ex is a reference to the node)
    virtual void simplifyMe(MX& ex) {}

    /// Get an IMatrix representation of a GetNonzeros or SetNonzeros node
    virtual Matrix<int> mapping() const;

    /// Create a horizontal concatenation node
    virtual MX getHorzcat(const std::vector<MX>& x) const;

    /// Create a horizontal split node
    virtual std::vector<MX> getHorzsplit(const std::vector<int>& output_offset) const;

    /// Create a vertical concatenation node (vectors only)
    virtual MX getVertcat(const std::vector<MX>& x) const;

    /// Create a vertical split node (vectors only)
    virtual std::vector<MX> getVertsplit(const std::vector<int>& output_offset) const;

    /// Create a diagonal concatenation node
    virtual MX getDiagcat(const std::vector<MX>& x) const;

    /// Create a diagonal split node
    virtual std::vector<MX> getDiagsplit(const std::vector<int>& offset1,
                                         const std::vector<int>& offset2) const;

    /// Transpose
    virtual MX getTranspose() const;

    /// Reshape
    virtual MX getReshape(const Sparsity& sp) const;

    /** \brief Matrix multiplication
    *
    *  The optimal argument sp_z will be used as the sparsity pattern of the result
    */
    virtual MX getMultiplication(const MX& y, const Sparsity& sp_z=Sparsity()) const;

    /** \brief Solve a system of linear equations
    *
    *      For system Ax = b:
    *
    *      A->getSolve(b)
    *
    */
    virtual MX getSolve(const MX& r, bool tr, const LinearSolver& linear_solver) const;

    /// Get the nonzeros of matrix
    virtual MX getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const;

    /// Assign the nonzeros of a matrix to another matrix
    virtual MX getSetNonzeros(const MX& y, const std::vector<int>& nz) const;

    /// Add the nonzeros of a matrix to another matrix
    virtual MX getAddNonzeros(const MX& y, const std::vector<int>& nz) const;

    /// Get submatrix reference
    virtual MX getSubRef(const Slice& i, const Slice& j) const;

    /// Get submatrix assignment
    virtual MX getSubAssign(const MX& y, const Slice& i, const Slice& j) const;

    /// Create set sparse
    virtual MX getSetSparse(const Sparsity& sp) const;

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

    /// Frobenius norm
    virtual MX getNormF() const;

    /// Spectral norm
    virtual MX getNorm2() const;

    /// Infinity norm
    virtual MX getNormInf() const;

    /// 1-norm
    virtual MX getNorm1() const;

    /// Assertion
    MX getAssertion(const MX& y, const std::string & fail_message="") const;

    /** Temporary variables to be used in user algorithms like sorting,
        the user is responsible of making sure that use is thread-safe
        The variable is initialized to zero
    */
    int temp;

    /** \brief  dependencies - functions that have to be evaluated before this one */
    std::vector<MX> dep_;

    /** \brief  The sparsity pattern */
    Sparsity sparsity_;

    /** \brief  Propagate sparsity, no work */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief Free adjoint memory (MX) */
    static void clearVector(const std::vector<std::vector<MX*> > v);
  };


  // Implementations

  template<typename T>
  std::vector<T> MXNode::getVector(const std::vector<T*> v) {
    std::vector<T> ret(v.size());
    for (int i=0; i<v.size(); i++) {
      if (v[i]!=0) {
        ret[i] = *v[i];
      }
    }
    return ret;
  }

  template<typename T>
  std::vector<std::vector<T> > MXNode::getVector(const std::vector<std::vector<T*> > v) {
    std::vector<std::vector<T> > ret(v.size());
    for (int i=0; i<v.size(); i++) {
      ret[i] = getVector(v[i]);
    }
    return ret;
  }
  /// \endcond

} // namespace casadi
/// \endcond

#endif // CASADI_MX_NODE_HPP
