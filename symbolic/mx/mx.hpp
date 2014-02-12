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

#ifndef MX_HPP
#define MX_HPP

#include "../shared_object.hpp"
#include "../matrix/matrix.hpp"
#include "../matrix/generic_expression.hpp"
#include <vector>
namespace CasADi{
  
  /** \brief  Forward declaration */
  class MXNode;
  class FX;


  /** \brief MX - Matrix expression
      The MX class is used to build up trees made up from MXNodes. It is a more general graph representation than the scalar expression,
      SX, and much less efficient for small objects. On the other hand, the class allows much more general operations than does SX,
      in particular matrix valued operations and calls to arbitrary differentiable functions.
  
      The MX class is designed to have identical syntax with the Matrix<> template class, and uses Matrix<double> as its internal 
      representation of the values at a node. By keeping the syntaxes identical, it is possible to switch from one class to the other, 
      as well as inlining MX functions to SX functions.
  
      Note that an operation is always "lazy", making a matrix multiplication will create a matrix multiplication node, not perform
      the actual multiplication.

      \author Joel Andersson 
      \date 2010-2011
  */
  class MX : public GenericExpression<MX>, public GenericMatrix<MX>, public SharedObject{
  
  public:
  
    /** \brief  Default constructor */
    MX();

    /** \brief  Construct a symbolic matrix (matrix variable) */
    explicit MX(const std::string& name, int n=1, int m=1);

    /** \brief  Construct a symbolic matrix (matrix variable) */
    explicit MX(const std::string& name, const std::pair<int,int> &nm);

    /** \brief  Construct a symbolic matrix (matrix variable) */
    explicit MX(const std::string& name, const CRSSparsity & sp);

    /** \brief  Construct MX with a given sparsity */
    explicit MX(const CRSSparsity& sp, const MX& val=0);
    
    /** \brief  Create scalar constant (also implicit type conversion) */
    MX(double x);

    /** \brief  Copy constructor */
    MX(const MX& x);

    /** \brief  Create vector constant (also implicit type conversion) */
    MX(const std::vector<double> &x);
    
    /** \brief  Create sparse matrix constant (also implicit type conversion) */
    MX(const Matrix<double> &x);

    /** \brief  Matrix with all zeros */
    MX(int nrow, int ncol);
    
    /** \brief  Dense matrix filled with value val */
    MX(int nrow, int ncol, const MX& val);
    
    /** \brief  Destructor */
    virtual ~MX();
    

#ifndef SWIG
    /** \brief  Create from node */
    static MX create(MXNode* node);

    /** \brief  Create from node (multiple-outputs) */
    static std::vector<MX> createMultipleOutput(MXNode* node);

    /// Get a non-zero element, with bounds checking
    const MX at(int k) const;

    /// Access a non-zero element, with bounds checking
    NonZeros<MX,int> at(int k);
    
#endif // SWIG
    
    /// Returns the truth value of an MX expression
    bool __nonzero__() const;
    
    //@{
    /// Indexing for interfaced languages
    
    /// get a non-zero
    const MX indexed_one_based(int k) const{ return at(k-1);}
    const MX indexed_zero_based(int k) const{ return at(k);}
    const MX indexed(const IndexList &k) const{
      return (*this)[k.getAll(size())];
    }
    const MX indexed(const Slice &k) const{ 
      return (*this)[k.getAll(size())];
    }
    
    /// get a matrix element
    const MX indexed_one_based(int i, int j) const{ return (*this)(i-1,j-1);}
    const MX indexed_zero_based(int i, int j) const{ return (*this)(i,j);}
    const MX indexed(const IndexList &i, const IndexList &j) const{ 
      return (*this)(i.getAll(size1()),j.getAll(size2()));
    }
    const MX indexed(const Slice &i, const Slice &j) const{ 
      return (*this)(i.getAll(size1()),j.getAll(size2()));
    }
    const MX indexed(const Matrix<int> &k) const{ 
      return (*this)(k);
    }
    const MX indexed(const CRSSparsity &sp) const{ 
      return (*this)(sp);
    }
    const MX indexed(const Slice &i, const Matrix<int>& k) const{ return (*this)(i,k); }
    const MX indexed(const IndexList &i, const Matrix<int>& k) const{ 
      return (*this)(i.getAll(size1()),k);
    }
    const MX indexed(const Matrix<int>& k, const Slice &j) const{ return (*this)(k,j); }
    const MX indexed(const Matrix<int>& k, const IndexList &j) const{ 
      return (*this)(k,j.getAll(size2()));
    }
    const MX indexed(const Matrix<int>& i, const Matrix<int>& j) const{ 
      return (*this)(i,j);
    }
    
    /// set a non-zero
    void indexed_one_based_assignment(int k, const MX &m){ at(k-1) = m(0,0);}
    void indexed_zero_based_assignment(int k, const MX &m){ at(k) = m(0,0);}
    void indexed_assignment(const IndexList &k, const MX &m){
      (*this)[k.getAll(size())] = m;
    }
    void indexed_assignment(const Slice &k, const MX &m){
      (*this)[k.getAll(size())] = m;
    }
    
    /// set a matrix element
    void indexed_one_based_assignment(int i, int j, const MX &m){ (*this)(i-1,j-1) = m;}
    void indexed_zero_based_assignment(int i, int j, const MX &m){ (*this)(i,j) = m;}
    void indexed_assignment(const IndexList &i, const IndexList &j, const MX &m){
      setSub(m,i.getAll(size1()),j.getAll(size2()));
    }
    
    void indexed_assignment(const Slice &i, const Slice &j, const MX &m){
      (*this)(i.getAll(size1()),j.getAll(size2())) = m;
    }
    
    void indexed_zero_based_assignment(const Matrix<int>& k, const MX &m){
      (*this)[k] = m;
    }
    void indexed_assignment(const CRSSparsity& sp, const MX &m){
      (*this)(sp) = m;
    }
    void indexed_assignment(const Slice &i, const Matrix<int>& j, const MX& m){
      (*this)(i.getAll(size1()),j) = m;
    }
    void indexed_assignment( const Matrix<int>& i, const Slice &j, const MX& m){
      (*this)(i,j.getAll(size2())) = m;
    }
    void indexed_assignment(const IndexList &i, const Matrix<int>& j, const MX& m){
      (*this)(i.getAll(size1()),j) = m;
    }
    void indexed_assignment( const Matrix<int>& i, const IndexList &j, const MX& m){
      (*this)(i,j.getAll(size2())) = m;
    } 
    void indexed_assignment( const Matrix<int>& i, const Matrix<int>& j, const MX& m){
      (*this)(i,j) = m;
    } 
    //@}
    
    /// Scalar type
    typedef MX ScalarType;

    /** \brief Get the sparsity pattern */
    const CRSSparsity& sparsity() const;

    /// Access the sparsity, make a copy if there are multiple references to it
    CRSSparsity& sparsityRef();

    /** \brief Erase a submatrix */
    void erase(const std::vector<int>& ii, const std::vector<int>& jj);

    /** \brief Enlarge matrix
        Make the matrix larger by inserting empty rows and columns, keeping the existing non-zeros */
    void enlarge(int nrow, int ncol, const std::vector<int>& ii, const std::vector<int>& jj);


    MX operator-() const;
  
    //@{
    /** \brief  Access a member of the node */
    MXNode* operator->();

    /** \brief  Const access a member of the node */
    const MXNode* operator->() const;
    //@}
  
    /** \brief Get the nth dependency as MX */
    MX getDep(int ch=0) const;
  
    /** \brief  Number of outputs */
    int getNumOutputs() const;
  
    /** \brief  Get an output */
    MX getOutput(int oind=0) const;

    /** \brief Get the number of dependencies of a binary SX */
    int getNdeps() const;
    
    /// Get the name.
    std::string getName() const;
  
    /// Get the value (only for scalar constant nodes)
    double getValue() const;
    
    /// Get the value (only for constant nodes)
    Matrix<double> getMatrixValue() const;
  
    /// Check if symbolic
    bool isSymbolic () const;
  
    /// Check if constant
    bool isConstant () const;
  
    /// Check if evaluation
    bool isEvaluation () const;
  
    /// Check if evaluation output
    bool isEvaluationOutput () const;
  
    /// Get the index of evaluation output - only valid when isEvaluationoutput() is true
    int getEvaluationOutput() const;
  
    /// Is it a certain operation
    bool isOperation (int op) const;
  
    /// Check if multiplication
    bool isMultiplication() const;
  
    /// Check if commutative operation
    bool isCommutative() const;
    
    /// Check if norm
    bool isNorm () const;
  
    /// Get function
    FX getFunction();

    /// Is binary operation
    bool isBinary() const;
  
    /// Is unary operation
    bool isUnary() const;
  
    /// Get operation type
    int getOp() const;

    /** \brief Check if two nodes are equivalent up to a given depth. 
     *  Depth=0 checks if the expressions are identical, i.e. points to the same node.
     * 
     *  a = x*x
     *  b = x*x
     *
     *  a.isEqual(b,0)  will return false, but a.isEqual(b,1) will return true
     */
    bool isEqual(const MX& y, int depth=0) const;
#ifndef SWIG
    bool isEqual(const MXNode* y, int depth=0) const;
#endif // SWIG
  
    /** \brief Returns a number that is unique for a given MXNode. 
     * If the MX does not point to any node, 0 is returned.
     */
    long __hash__() const;
    
    /// Get the temporary variable
    int getTemp() const;
  
    /// Set the temporary variable
    void setTemp(int t);
  
    //@{
    /** \brief  Create nodes by their ID */
    static MX binary(int op, const MX &x, const MX &y);
    static MX unary(int op, const MX &x);
    //@}

    //@{
    /** \brief  Sparse matrix of all zeros */
    static MX sparse(int nrow, int ncol=1);
    static MX sparse(const std::pair<int, int> &nm);
    //@}
  
    //@{
    /** \brief  Dense matrix of all zeros */
    static MX zeros(const CRSSparsity& sp);
    static MX zeros(int nrow, int ncol=1); 
    static MX zeros(const std::pair<int, int> &nm);
    //@}

    //@{
    /** \brief  Matrix of all ones */  
    static MX ones(const CRSSparsity& sp);
    static MX ones(int nrow, int ncol=1); 
    static MX ones(const std::pair<int, int> &nm);
    //@}

    //@{
    /** \brief  create a matrix with all inf */
    static MX inf(const CRSSparsity& sp);
    static MX inf(int nrow=1, int ncol=1);
    static MX inf(const std::pair<int,int>& nm);
    //@}
  
    //@{
    /** \brief  create a matrix with all nan */
    static MX nan(const CRSSparsity& sp);
    static MX nan(int nrow=1, int ncol=1);
    static MX nan(const std::pair<int,int>& nm);
    //@}
  
    //@{
    /** \brief  create a matrix by repeating an existing matrix */
    static MX repmat(const MX& x, int nrow, int ncol=1);
    static MX repmat(const MX& x, const std::pair<int, int> &nm);
    //@}

    /** \brief  Identity matrix */  
    static MX eye(int nrow);
  
    const MX sub(int i, int j) const;
    const MX sub(int i, const std::vector<int>& j) const;
    const MX sub(const std::vector<int>& i, int j) const;
    const MX sub(const std::vector<int>& i, const std::vector<int>& j) const;
    const MX sub(const Matrix<int>& k, int dummy=0) const;
    const MX sub(const CRSSparsity& sp, int dummy=0) const;
    const MX sub(const std::vector<int>& i, const Matrix<int>& j) const;
    const MX sub(const Matrix<int>& k, const std::vector<int>& j) const;
    const MX sub(const Slice& i, int j) const {return sub(i.getAll(size1()),j);}
    const MX sub(int i, const Slice& j) const {return sub(i,j.getAll(size2()));}
    const MX sub(const Slice& i, const Slice& j) const {return sub(i.getAll(size1()),j.getAll(size2()));}
    const MX sub(const Slice& i, const Matrix<int>& j) const {return sub(i.getAll(size1()),j);}
    const MX sub(const Matrix<int>& i, const Slice& j) const {return sub(i,j.getAll(size2()));}
    const MX sub(const Matrix<int>& i, const Matrix<int>& j) const;

    void setSub(const MX& m, int i, int j);
    void setSub(const MX& m, int i, const std::vector<int>& j);
    void setSub(const MX& m, const std::vector<int>& i, int j);
    void setSub(const MX& m, const std::vector<int>& i, const std::vector<int>& j);
    void setSub(const MX& m, const Matrix<int>& k);
    void setSub(const MX& m, const std::vector<int>& i, const Matrix<int>& k);
    void setSub(const MX& m, const Matrix<int>& k, const std::vector<int>& j);
    void setSub(const MX& m, const Slice& i, const Slice& j);
    //void setSub(const MX& m, const Slice& i, const Matrix<int>& k) {return setSub(m,i.getAll(size1()),k);}
    //void setSub(const MX& m, const Matrix<int>& k, const Slice& j) {return setSub(m,k,j.getAll(size2()));}
    void setSub(const MX& m, const Matrix<int>& i, const Matrix<int>& j);
    void setSub(const MX& m, const CRSSparsity& sp, int dummy);
    
    MX getNZ(int k) const;
    MX getNZ(const std::vector<int>& k) const;
    MX getNZ(const Slice& k) const{ return getNZ(k.getAll(size()));}
    MX getNZ(const Matrix<int>& k) const;
    void setNZ(int k, const MX& el);
    void setNZ(const std::vector<int>& k, const MX& el);
    void setNZ(const Slice& k, const MX& m){ setNZ(k.getAll(size()),m);}
    void setNZ(const Matrix<int>& k, const MX& m);

    /** \brief Append a matrix to the end. */
    void append(const MX& y);
  
    // all binary operations
    MX __add__(const MX& y) const;
    MX __sub__(const MX& y) const;
    MX __mul__(const MX& y) const;
    MX __div__(const MX& y) const;
    MX __lt__(const MX& y) const;
    MX __le__(const MX& y) const;
    MX __eq__(const MX& y) const;
    MX __ne__(const MX& y) const;
    MX __truediv__(const MX& y) const { return __div__(y);};
    MX __pow__(const MX& b) const;
    MX __constpow__(const MX& b) const;
    MX __mrdivide__  (const MX& b) const;
    MX __mpower__(const MX& b) const;
    MX mul(const MX& y, const CRSSparsity &sp_z=CRSSparsity()) const;
    MX mul_full(const MX& y, const CRSSparsity &sp_z=CRSSparsity()) const;
    MX inner_prod(const MX& y) const;
    MX outer_prod(const MX& y) const;
    MX constpow(const MX& y) const;
    MX fmin(const MX& y) const;
    MX fmax(const MX& y) const;
    MX printme(const MX& y) const;
    MX arctan2(const MX& y) const;
    MX logic_and(const MX& y) const;
    MX logic_or(const MX& y) const;
    MX if_else_zero(const MX& y) const;
    MX __copysign__(const MX& y) const;

    // all unary operations
    MX exp() const;
    MX log() const;
    MX log10() const;
    MX sqrt() const;
    MX sin() const;
    MX cos() const;
    MX tan() const;
    MX arcsin() const;
    MX arccos() const;
    MX arctan() const;
    MX floor() const;
    MX ceil() const;
    MX fabs() const;
    MX sign() const;
    MX erfinv() const;
    MX erf() const;
    MX sinh() const;
    MX cosh() const;
    MX tanh() const;
    MX arcsinh() const;
    MX arccosh() const;
    MX arctanh() const;
    MX logic_not() const;
    MX relay() const;
    
    /** \brief returns itself, but with an assertion attached
    *
    *  If y does not evaluate to 1, a runtime error is raised
    */
    MX attachAssert(const MX& y,const std::string &fail_message="") const;

    /** \brief Set sparse */
    MX setSparse(const CRSSparsity& sp, bool intersect=false) const;

    /** \brief Make dense */
    MX makeDense(const MX& val = 0) const;

    /// Lift an expression
    void lift(const MX& x_guess);

    /** \brief Get an IMatrix representation of a GetNonzeros or SetNonzeros node */
    Matrix<int> mapping() const;
    
    /** \brief Set or reset the maximum number of calls to the printing function when printing an expression */
    static void setMaxNumCallsInPrint(long num=10000);

    /** \brief Get the maximum number of calls to the printing function when printing an expression */
    static long getMaxNumCallsInPrint();

    /** \brief Set or reset the depth to which equalities are being checked for simplifications */
    static void setEqualityCheckingDepth(int eq_depth=1);

    /** \brief Get the depth to which equalities are being checked for simplifications */
    static int getEqualityCheckingDepth();

#ifndef SWIG
  private:
    // Maximum number of calls
    static long max_num_calls_in_print_;
  
    // Depth when checking equalities
    static int eq_depth_;
  
#endif // SWIG
  };

  //@{
  /// Some typedefs
  typedef std::vector<MX> MXVector;
  typedef std::vector< std::vector<MX> > MXVectorVector;
  typedef MX* MXPtr;
  typedef std::vector<MXPtr> MXPtrV;
  typedef std::vector<MXPtrV> MXPtrVV;
  //@}

} // namespace CasADi

#endif // MX_HPP
