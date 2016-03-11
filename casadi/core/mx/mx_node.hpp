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
#include "../sx/sx_elem.hpp"
#include "../calculus.hpp"
#include "../function/code_generator.hpp"
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
  class CASADI_EXPORT MXNode : public SharedObjectNode {
    friend class MX;

  public:
    /// Constructor
    MXNode();

    /** \brief  Destructor */
    virtual ~MXNode()=0;

    /** \brief Check the truth value of this node
     */
    virtual bool __nonzero__() const;

    /** \brief Check if identically zero */
    virtual bool is_zero() const { return false;}

    /** \brief Check if identically one */
    virtual bool is_one() const { return false;}

    /** \brief Check if a certain value */
    virtual bool isValue(double val) const { return false;}

    /** \brief Check if identity matrix */
    virtual bool is_identity() const { return false;}

    /** \brief Check if unary operation */
    virtual bool is_unaryOp() const { return false;}

    /** \brief Check if binary operation */
    virtual bool is_binaryOp() const { return false;}

    /** \brief  Print a representation */
    virtual void repr(std::ostream &stream) const;

    /** \brief  Print a description */
    virtual void print(std::ostream &stream) const;

    /** \brief Find out which nodes can be inlined */
    void can_inline(std::map<const MXNode*, int>& nodeind) const;

    /** \brief Print compact */
    std::string print_compact(std::map<const MXNode*, int>& nodeind,
                             std::vector<std::string>& intermed) const;

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const = 0;

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const {}

    /** \brief Is reference counting needed in codegen? */
    virtual bool has_refcount() const { return false;}

    /** \brief Codegen incref */
    virtual void codegen_incref(CodeGenerator& g, std::set<void*>& added) const {}

    /** \brief Codegen decref */
    virtual void codegen_decref(CodeGenerator& g, std::set<void*>& added) const {}

    /** \brief Generate code for the operation */
    virtual void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const;

    /** \brief  Evaluate numerically */
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem) const;

    /** \brief  Evaluate symbolically (SX) */
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Get the name */
    virtual const std::string& name() const;

    /** \brief  Check if valid function input */
    virtual bool is_valid_input() const { return false;}

    /** \brief Get the number of symbolic primitives */
    virtual int n_primitives() const;

    /** \brief Get symbolic primitives */
    virtual void primitives(std::vector<MX>::iterator& it) const;

    /** \brief Split up an expression along symbolic primitives */
    virtual void split_primitives(const MX& x, std::vector<MX>::iterator& it) const;

    /** \brief Join an expression along symbolic primitives */
    virtual MX join_primitives(std::vector<MX>::const_iterator& it) const;

    /** \brief Detect duplicate symbolic expressions */
    virtual bool has_duplicates();

    /** \brief Reset the marker for an input expression */
    virtual void resetInput();

    /** \brief  Check if evaluation output */
    virtual bool isOutputNode() const {return false;}

    /** \brief  Check if a multiple output node */
    virtual bool isMultipleOutput() const {return false;}

    /** \brief  Number of functions */
    virtual int numFunctions() const {return 0;}

    /** \brief  Get function reference */
    virtual const Function& getFunction(int i) const;

    /** \brief  Get function input */
    virtual int getFunction_input() const;

    /** \brief  Get function output */
    virtual int getFunctionOutput() const;

    /** \brief Get the operation */
    virtual int op() const = 0;

    /** \brief Check if two nodes are equivalent up to a given depth */
    static bool is_equal(const MXNode* x, const MXNode* y, int depth);
    virtual bool is_equal(const MXNode* node, int depth) const { return false;}

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
    virtual int nout() const { return 1;}

    /** \brief  Get an output */
    virtual MX getOutput(int oind) const;

    /// Get the sparsity
    const Sparsity& sparsity() const { return sparsity_;}

    /// Get the sparsity of output oind
    virtual const Sparsity& sparsity(int oind) const;

    /// Get shape
    int numel() const { return sparsity().numel(); }
    int nnz(int i=0) const { return sparsity(i).nnz(); }
    int size1() const { return sparsity().size1(); }
    int size2() const { return sparsity().size2(); }
    std::pair<int, int> size() const { return sparsity().size();}

    /** \brief Is the node nonlinear */
    virtual bool isNonLinear() {return false;}

    /// Set the sparsity
    void setSparsity(const Sparsity& sparsity);

    /** \brief Get required length of arg field */
    virtual size_t sz_arg() const { return ndep();}

    /** \brief Get required length of res field */
    virtual size_t sz_res() const { return nout();}

    /** \brief Get required length of iw field */
    virtual size_t sz_iw() const { return 0;}

    /** \brief Get required length of w field */
    virtual size_t sz_w() const { return 0;}

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
      if (x.size()==sp.size()) {
        return x;
      } else {
        return MX(sp, x);
      }
    }

    /// Get the value (only for scalar constant nodes)
    virtual double to_double() const;

    /// Get the value (only for constant nodes)
    virtual Matrix<double> getMatrixValue() const;

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual int numInplace() const { return 0;}

    /// Simplify the expression (ex is a reference to the node)
    virtual void simplifyMe(MX& ex) {}

    /// Get an IM representation of a GetNonzeros or SetNonzeros node
    virtual Matrix<int> mapping() const;

    /// Create a horizontal concatenation node
    virtual MX getHorzcat(const std::vector<MX>& x) const;

    /// Create a horizontal split node
    virtual std::vector<MX> getHorzsplit(const std::vector<int>& output_offset) const;

    /// Create a repeated matrix node
    virtual MX getRepmat(int m, int n) const;

    /// Create a repeated sum node
    virtual MX getRepsum(int m, int n) const;

    /// Create a vertical concatenation node (vectors only)
    virtual MX getVertcat(const std::vector<MX>& x) const;

    /// Create a vertical split node (vectors only)
    virtual std::vector<MX> getVertsplit(const std::vector<int>& output_offset) const;

    /// Create a diagonal concatenation node
    virtual MX get_diagcat(const std::vector<MX>& x) const;

    /// Create a diagonal split node
    virtual std::vector<MX> get_diagsplit(const std::vector<int>& offset1,
                                         const std::vector<int>& offset2) const;

    /// Transpose
    virtual MX getTranspose() const;

    /// Reshape
    virtual MX getReshape(const Sparsity& sp) const;

    /** \brief Matrix multiplication and addition */
    virtual MX getMultiplication(const MX& y, const MX& z) const;

    /** \brief Bilinear form */
    virtual MX getBilin(const MX& x, const MX& y) const;

    /** \brief Bilinear form */
    virtual MX getRank1(const MX& alpha, const MX& x, const MX& y) const;

    /** \brief Solve a system of linear equations
    *
    *      For system Ax = b:
    *
    *      A->getSolve(b)
    *
    */
    virtual MX getSolve(const MX& r, bool tr, const Function& linear_solver) const;

    /// Get the nonzeros of matrix
    virtual MX getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const;

    /// Assign the nonzeros of a matrix to another matrix
    virtual MX getSetNonzeros(const MX& y, const std::vector<int>& nz) const;

    /// Add the nonzeros of a matrix to another matrix
    virtual MX getAddNonzeros(const MX& y, const std::vector<int>& nz) const;

    /// Get submatrix reference
    virtual MX getRef(const Slice& i, const Slice& j) const;

    /// Get submatrix assignment
    virtual MX getAssign(const MX& y, const Slice& i, const Slice& j) const;

    /// Create set sparse
    virtual MX getProject(const Sparsity& sp) const;

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
    virtual MX getDot(const MX& y) const;

    /// Frobenius norm
    virtual MX getNormF() const;

    /// Spectral norm
    virtual MX getNorm2() const;

    /// Infinity norm
    virtual MX getNormInf() const;

    /// 1-norm
    virtual MX getNorm1() const;

    /// Assertion
    MX getAssertion(const MX& y, const std::string& fail_message) const;

    /// Monitor
    MX getMonitor(const std::string& comment) const;

    /// Find
    MX getFind() const;

    /** Temporary variables to be used in user algorithms like sorting,
        the user is responsible of making sure that use is thread-safe
        The variable is initialized to zero
    */
    int temp;

    /** \brief  dependencies - functions that have to be evaluated before this one */
    std::vector<MX> dep_;

    /** \brief  The sparsity pattern */
    Sparsity sparsity_;

    /** \brief Propagate sparsities forward through a copy operation */
    static void copyFwd(const bvec_t* arg, bvec_t* res, int len);

    /** \brief Propagate sparsities backwards through a copy operation */
    static void copyAdj(bvec_t* arg, bvec_t* res, int len);
  };

  /// \endcond
} // namespace casadi

#endif // CASADI_MX_NODE_HPP
