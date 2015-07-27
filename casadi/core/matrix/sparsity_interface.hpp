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


#ifndef CASADI_SPARSITY_INTERFACE_HPP
#define CASADI_SPARSITY_INTERFACE_HPP

#include "../std_vector_tools.hpp"

namespace casadi {
  /** \brief Sparsity interface class

      This is a common base class for GenericMatrix (i.e. MX and Matrix<>) and Sparsity, introducing a
      uniform syntax and implementing common functionality using the curiously recurring template pattern
      (CRTP) idiom.\n

      \author Joel Andersson
      \date 2014
  */
  template<typename MatType>
  class CASADI_EXPORT SparsityInterface {
#ifndef SWIG
  protected:
    // Helper functions
    inline const MatType& self() const { return static_cast<const MatType&>(*this); }
    inline MatType& self() { return static_cast<MatType&>(*this); }
  public:
    // Create vector with 1 element
    inline friend std::vector<MatType> make_vector(const MatType& x0) {
      return std::vector<MatType>(1, x0);
    }

    // Create vector with 2 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1) {
      MatType x[] = {x0, x1};
      return std::vector<MatType>(x, x+2);
    }

    // Create vector with 3 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1,
                                                   const MatType& x2) {
      MatType x[] = {x0, x1, x2};
      return std::vector<MatType>(x, x+3);
    }

    // Create vector with 4 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1,
                                                   const MatType& x2, const MatType& x3) {
      MatType x[] = {x0, x1, x2, x3};
      return std::vector<MatType>(x, x+4);
    }

    // Create vector with 5 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1,
                                                   const MatType& x2, const MatType& x3,
                                                   const MatType& x4) {
      MatType x[] = {x0, x1, x2, x3, x4};
      return std::vector<MatType>(x, x+5);
    }

    // Create vector with 6 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1,
                                                   const MatType& x2, const MatType& x3,
                                                   const MatType& x4, const MatType& x5) {
      MatType x[] = {x0, x1, x2, x3, x4, x5};
      return std::vector<MatType>(x, x+6);
    }

    // Create vector from map and vector with key order
    inline friend std::vector<MatType> make_vector(const std::map<std::string, MatType>& m,
                                                   const std::vector<std::string>& s) {
      std::vector<MatType> ret(s.size());
      for (size_t i=0; i!=s.size(); ++i) {
        typename std::map<std::string, MatType>::const_iterator it=m.find(s[i]);
        if (it!=m.end()) {
          ret[i]=it->second;
        }
      }
      return ret;
    }

    // Create vector from map and vector with key order
    inline friend std::vector<std::vector<MatType> >
      make_vector(const std::map<std::string, std::vector<MatType> >& m,
                  const std::vector<std::string>& s) {
      std::vector<std::vector<MatType> > ret(s.size());
      for (size_t i=0; i!=s.size(); ++i) {
        typename std::map<std::string, std::vector<MatType> >::const_iterator it=m.find(s[i]);
        if (it!=m.end()) {
          ret[i]=it->second;
        }
      }
      return ret;
    }

    // Create vector from map and vector with key order
    inline friend std::vector<MatType>
      make_vector(const std::pair<std::map<std::string, MatType>, std::vector<std::string> >& ms) {
      return make_vector(ms.first, ms.second);
    }

    // Create vector from map and vector with key order
    inline friend std::vector<std::vector<MatType> >
      make_vector(const std::pair<std::map<std::string, std::vector<MatType> >,
                  std::vector<std::string> >& ms) {
      return make_vector(ms.first, ms.second);
    }

    // Assign 1 element from a vector
    template<typename T0>
      inline friend void assign_vector(T0& x0,
                                       const std::vector<MatType>& x) {
      x.at(0).get(x0);
    }

    // Assign 2 elements from a vector
    template<typename T0, typename T1>
      inline friend void assign_vector(T0& x0, T1& x1,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x);
      x.at(1).get(x1);
    }

    // Assign 3 elements from a vector
    template<typename T0, typename T1, typename T2>
      inline friend void assign_vector(T0& x0, T1& x1, T2& x2,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x1, x);
      x.at(2).get(x2);
    }

    // Assign 4 elements from a vector
    template<typename T0, typename T1, typename T2, typename T3>
      inline friend void assign_vector(T0& x0, T1& x1, T2& x2, T3& x3,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x1, x2, x);
      x.at(3).get(x3);
    }

    // Assign 5 elements from a vector
    template<typename T0, typename T1, typename T2, typename T3, typename T4>
      inline friend void assign_vector(T0& x0, T1& x1, T2& x2, T3& x3, T4& x4,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x1, x2, x3, x);
      x.at(4).get(x4);
    }

    // Assign 6 elements from a vector
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
      inline friend void assign_vector(T0& x0, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x1, x2, x3, x4, x);
      x.at(5).get(x5);
    }

    // Create map with 1 element
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0) {
      std::map<std::string, MatType> ret;
      ret[n0]=x0;
      return ret;
    }

    // Create map with 2 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1) {
      std::map<std::string, MatType> ret=make_map(n0, x0);
      ret[n1]=x1;
      return ret;
    }

    // Create map with 3 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1,
               const std::string& n2, const MatType& x2) {
      std::map<std::string, MatType> ret=make_map(n0, x0, n1, x1);
      ret[n2]=x2;
      return ret;
    }

    // Create map with 4 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1,
               const std::string& n2, const MatType& x2,
               const std::string& n3, const MatType& x3) {
      std::map<std::string, MatType> ret=make_map(n0, x0, n1, x1, n2, x2);
      ret[n3]=x3;
      return ret;
    }

    // Create map with 5 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1,
               const std::string& n2, const MatType& x2,
               const std::string& n3, const MatType& x3,
               const std::string& n4, const MatType& x4) {
      std::map<std::string, MatType> ret=make_map(n0, x0, n1, x1, n2, x2, n3, x3);
      ret[n4]=x4;
      return ret;
    }

    // Create map with 6 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1,
               const std::string& n2, const MatType& x2,
               const std::string& n3, const MatType& x3,
               const std::string& n4, const MatType& x4,
               const std::string& n5, const MatType& x5) {
      std::map<std::string, MatType> ret=make_map(n0, x0, n1, x1, n2, x2, n3, x3, n4, x4);
      ret[n5]=x5;
      return ret;
    }

    // Assign 1 element from a map
    template<typename T0>
    inline friend void assign_map(const std::string& n0, T0& x0,
                                  const std::map<std::string, MatType>& x) {
      x.at(n0).get(x0);
    }

    // Assign 2 elements from a map
    template<typename T0, typename T1>
      inline friend void assign_map(const std::string& n0, T0& x0,
                                    const std::string& n1, T0& x1,
                                    const std::vector<MatType>& x) {
      assign_map(n0, x0, x);
      x.at(n1).get(x1);
    }

    // Assign 3 elements from a map
    template<typename T0, typename T1, typename T2>
    inline friend void assign_map(const std::string& n0, T0& x0,
                                  const std::string& n1, T0& x1,
                                  const std::string& n2, T0& x2,
                                  const std::vector<MatType>& x) {
      assign_map(n0, x0, n1, x1, x);
      x.at(n2).get(x2);
    }

    // Assign 4 elements from a map
    template<typename T0, typename T1, typename T2, typename T3>
      inline friend void assign_map(const std::string& n0, T0& x0,
                                    const std::string& n1, T0& x1,
                                    const std::string& n2, T0& x2,
                                    const std::string& n3, T0& x3,
                                    const std::vector<MatType>& x) {
      assign_map(n0, x0, n1, x1, n2, x2, x);
      x.at(n3).get(x3);
    }

    // Assign 5 elements from a map
    template<typename T0, typename T1, typename T2, typename T3, typename T4>
      inline friend void assign_map(const std::string& n0, T0& x0,
                                    const std::string& n1, T0& x1,
                                    const std::string& n2, T0& x2,
                                    const std::string& n3, T0& x3,
                                    const std::string& n4, T0& x4,
                                    const std::vector<MatType>& x) {
      assign_map(n0, x0, n1, x1, n2, x2, n3, x3, x);
      x.at(n4).get(x4);
    }

    // Assign 6 elements from a map
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
      inline friend void assign_map(const std::string& n0, T0& x0,
                                    const std::string& n1, T0& x1,
                                    const std::string& n2, T0& x2,
                                    const std::string& n3, T0& x3,
                                    const std::string& n4, T0& x4,
                                    const std::string& n5, T0& x5,
                                    const std::vector<MatType>& x) {
      assign_map(n0, x0, n1, x1, n2, x2, n3, x3, n4, x4, x);
      x.at(n5).get(x5);
    }
#endif // SWIG

  public:

    /// \cond CLUTTER
    std::vector< std::vector< MatType > >
      zz_blocksplit(const std::vector<int>& vert_offset,
                    const std::vector<int>& horz_offset) const;
    std::vector< std::vector< MatType > >
      zz_blocksplit(int vert_incr, int horz_incr) const;
    static MatType zz_veccat(const std::vector< MatType >& x);
    MatType zz_vec() const;
    MatType zz_repmat(int n, int m=1) const;
    static std::vector<int> zz_offset(const std::vector< MatType > &v, bool vert=true);
    std::vector< MatType > zz_diagsplit(const std::vector<int>& output_offset) const;
    std::vector< MatType > zz_diagsplit(int incr) const;
    std::vector< MatType > zz_diagsplit(int incr1, int incr2) const;
    static MatType zz_mul(const std::vector<MatType> &args);
    std::vector<MatType > zz_horzsplit(int incr) const;
    std::vector<MatType > zz_vertsplit(int incr) const;
    /// \endcond

    /*! \fn inline friend MatType horzcat(const std::vector<MatType> &v)
     *  \brief Concatenate a list of matrices horizontally
     * Alternative terminology: horizontal stack, hstack, horizontal append, [a b]
     *
     *   horzcat(horzsplit(x, ...)) = x
     */

    /*! \fn friend MatType vertcat(const std::vector<MatType> &v)
     * \brief Concatenate a list of matrices vertically
     * Alternative terminology: vertical stack, vstack, vertical append, [a;b]
     *
     *   vertcat(vertsplit(x, ...)) = x
     */


    /*! \fn inline friend std::vector<MatType > horzsplit(const MatType &v,
      const std::vector<int>& offset)
      \brief  split horizontally, retaining groups of columns
      * \param offset List of all start columns for each group
      *      the last column group will run to the end.
      *
      *   horzcat(horzsplit(x, ...)) = x
      */

    /*! \fn inline friend std::vector<MatType > horzsplit(const MatType &v, int incr=1)
      \brief  split horizontally, retaining fixed-sized groups of columns
      * \param incr Size of each group of columns
      *
      *   horzcat(horzsplit(x, ...)) = x
      */

    /*! \fn friend std::vector<MatType > vertsplit(const MatType &v,
      const std::vector<int>& offset)
      * \brief  split vertically, retaining groups of rows
      * \param output_offset List of all start rows for each group
      *      the last row group will run to the end.
      *
      *   vertcat(vertsplit(x, ...)) = x
      */

    /*! \fn inline friend std::vector<int > offset(const std::vector<MatType> &v, bool vert=true)
      \brief Helper function, get offsets corresponding to a vector of matrices
    */

    /*! \fn inline friend std::vector<MatType > vertsplit(const MatType &v, int incr=1)
      \brief  split vertically, retaining fixed-sized groups of rows
      * \param incr Size of each group of rows
      *
      *   vertcat(vertsplit(x, ...)) = x
     
      \doctest
      print vertsplit(SX.sym("a",4))
      \doctestout
      [SX(a_0), SX(a_1), SX(a_2), SX(a_3)]
      \enddoctest
     
      \doctest
      print vertsplit(SX.sym("a",4),2)
      \doctestout
      [SX([a_0, a_1]), SX([a_2, a_3])]
      \enddoctest
     
      If the number of rows is not a multiple of \p incr,
      the last entry returned will have a size smaller than \p incr.
     
      \doctest
      print vertsplit(DMatrix([0,1,2,3,4]),2)
      \doctestout
      [DMatrix([0, 1]), DMatrix([2, 3]), DMatrix(4)]
      \enddoctest
      *
      */

    /*! \fn inline friend MatType blockcat(const std::vector< std::vector<MatType > > &v)
     * \brief Construct a matrix from a list of list of blocks.
     */

    /*! \fn inline friend MatType blockcat(const MatType &A, const MatType &B,
      const MatType &C, const MatType &D)
      \brief Construct a matrix from 4 blocks
    */

    /*! \fn inline friend std::vector< std::vector< MatType > > 
      blocksplit(const MatType& x,
      const std::vector<int>& vert_offset,
      const std::vector<int>& horz_offset)
  
      \brief  chop up into blocks
      * \param vert_offset Defines the boundaries of the block rows
      * \param horz_offset Defines the boundaries of the block columns
      *
      *   blockcat(blocksplit(x,..., ...)) = x
      */

    /*! \fn inline friend std::vector< std::vector< MatType > >
      blocksplit(const MatType& x, int vert_incr=1, int horz_incr=1)
      \brief  chop up into blocks
      * \param vert_incr Defines the increment for block boundaries in row dimension
      * \param horz_incr Defines the increment for block boundaries in column dimension
      *
      *   blockcat(blocksplit(x,..., ...)) = x
      */

    /*! \fn inline friend MatType diagcat(const std::vector<MatType> &A)
      \brief Construct a matrix with given block on the diagonal
    */


    /*! \fn friend std::vector< MatType > diagsplit(const MatType& x,
      const std::vector<int>& output_offset1,
      const std::vector<int>& output_offset2)
      \brief  split diagonally, retaining square matrices
      * \param output_offset1 List of all start locations (row) for each group
      *      the last matrix will run to the end.
      * \param output_offset2 List of all start locations (row) for each group
      *      the last matrix will run to the end.
      *
      *   diagcat(diagsplit(x, ...)) = x
      */

    /*! \fn inline friend std::vector< MatType > diagsplit(const MatType& x,
      const std::vector<int>& output_offset)
      \brief  split diagonally, retaining square matrices
      * \param output_offset List of all start locations for each group
      *      the last matrix will run to the end.
      *
      *   diagcat(diagsplit(x, ...)) = x
      */

    /*! \fn inline friend std::vector< MatType > diagsplit(const MatType& x, int incr=1)
      \brief  split diagonally, retaining groups of square matrices
      * \param incr Size of each matrix
      *
      *  diagsplit(diagsplit(x, ...)) = x
      */

    /*! \fn inline friend std::vector< MatType > diagsplit(const MatType& x, int incr1, int incr2)
      \brief  split diagonally, retaining fixed-sized matrices
      * \param incr1 Row dimension of each matrix
      * \param incr2 Column dimension of each matrix
      *
      *  diagsplit(diagsplit(x, ...)) = x
      */

    /*! \fn inline friend MatType veccat(const std::vector< MatType >& x)
      \brief  concatenate vertically while vectorizing all arguments with vec
    */

    /*! \fn inline friend MatType mul(const MatType &X, const MatType &Y)
      \brief Matrix product of two matrices
    */

    /*! \fn inline friend MatType mul(const MatType &X, const MatType &Y, const MatType &Z)
      \brief Matrix product and addition
      Matrix product of two matrices (X and Y), adding the result to
      a third matrix Z. The result has the same sparsity pattern as
      C meaning that other entries of (X*Y) are ignored.
      The operation is equivalent to: Z+mul(X,Y).project(Z.sparsity()).
    */

    /*! \fn inline friend MatType mul(const std::vector<MatType> &args)
      \brief Matrix product of n matrices
    */

    /*! \fn inline friend MatType transpose(const MatType& X)
      \brief Transpose
    */

    /*! \fn inline friend MatType vec(const MatType& a)
      \brief  make a vector
      Reshapes/vectorizes the matrix such that the shape becomes (expr.numel(), 1).
      Columns are stacked on top of each other.
      Same as reshape(expr, expr.numel(), 1)
  
      a c \n
      b d \n

      turns into

      a \n
      b \n
      c \n
      d \n
    */

    /*! inline friend MatType vecNZ(const MatType& a)
      \brief Returns a flattened version of the matrix, preserving only nonzeros
    */

    /*! \fn inline friend MatType reshape(const MatType& a, int nrow, int ncol)
      \brief Returns a reshaped version of the matrix
    */

    /*! \fn inline friend MatType reshape(const MatType& a, std::pair<int, int> rc)
      \brief Returns a reshaped version of the matrix, dimensions as a vector
    */

    /*! \fn inline friend MatType reshape(const MatType& a, const Sparsity& sp)
      \brief Reshape the matrix
    */

    /*! inline friend int sprank(const MatType& A)
      \brief Obtain the structural rank of a sparsity-pattern
    */

    /*! inline friend int norm_0_mul(const MatType &x, const MatType &y)
      \brief 0-norm (nonzero count) of a Matrix-matrix product
    */

    /*! inline friend MatType triu(const MatType& a, bool includeDiagonal=true)
      \brief Get the upper triangular part of a matrix
    */

    /*! inline friend MatType tril(const MatType& a, bool includeDiagonal=true)
      \brief Get the lower triangular part of a matrix
    */

    /*! inline friend MatType kron(const MatType& a, const MatType& b)
      \brief Kronecker tensor product
      *
      * Creates a block matrix in which each element (i, j) is a_ij*b
      */

    /*! inline friend MatType repmat(const MatType &A, int n, int m=1)
      \brief Repeat matrix A n times vertically and m times horizontally
    */

    /*! inline friend MatType repmat(const MatType &A, const std::pair<int, int>& rc)
      \brief Repeat matrix A n times vertically and m times horizontally
    */


#define SPARSITY_INTERFACE_FRIENDS(M)                                   \
    inline SWIG_FRIEND M horzcat(const std::vector<M> &v) {             \
      return M::zz_horzcat(v);                                          \
    }                                                                   \
    inline SWIG_FRIEND M vertcat(const std::vector<M> &v) {             \
      return M::zz_vertcat(v);                                          \
    }                                                                   \
    inline SWIG_FRIEND std::vector<M >                                  \
      horzsplit(const M &v, const std::vector<int>& offset) {           \
      return v.zz_horzsplit(offset);                                    \
    }                                                                   \
    inline SWIG_FRIEND std::vector<M >                                  \
      horzsplit(const M &v, int incr=1) {                               \
      return v.zz_horzsplit(incr);                                      \
    }                                                                   \
    inline SWIG_FRIEND std::vector<M >                                  \
      vertsplit(const M &v, const std::vector<int>& offset) {           \
      return v.zz_vertsplit(offset);                                    \
    }                                                                   \
    inline SWIG_FRIEND std::vector<int >                                \
      offset(const std::vector<M> &v, bool vert=true) {                 \
      return M::zz_offset(v, vert);                                     \
    }                                                                   \
    inline SWIG_FRIEND std::vector<M >                                  \
      vertsplit(const M &v, int incr=1) {                               \
      return v.zz_vertsplit(incr);                                      \
    }                                                                   \
    inline SWIG_FRIEND M                                                \
      blockcat(const std::vector< std::vector<M > > &v) {               \
      return M::zz_blockcat(v);                                         \
    }                                                                   \
    inline SWIG_FRIEND M                                                \
      blockcat(const M &A, const M &B, const M &C, const M &D) {        \
      return vertcat(horzcat(A, B), horzcat(C, D));                     \
    }                                                                   \
    inline SWIG_FRIEND std::vector< std::vector< M > >                  \
      blocksplit(const M& x, const std::vector<int>& vert_offset,       \
                 const std::vector<int>& horz_offset) {                 \
      return x.zz_blocksplit(vert_offset, horz_offset);                 \
    }                                                                   \
    inline SWIG_FRIEND std::vector< std::vector< M > >                  \
      blocksplit(const M& x, int vert_incr=1, int horz_incr=1) {        \
      return x.zz_blocksplit(vert_incr, horz_incr);                     \
    }                                                                   \
    inline SWIG_FRIEND M diagcat(const std::vector<M> &A) {             \
      return M::zz_diagcat(A);                                          \
    }                                                                   \
    inline SWIG_FRIEND std::vector< M >                                 \
      diagsplit(const M& x, const std::vector<int>& output_offset1,     \
                const std::vector<int>& output_offset2) {               \
      return x.zz_diagsplit(output_offset1, output_offset2);            \
    }                                                                   \
    inline SWIG_FRIEND std::vector< M >                                 \
      diagsplit(const M& x, const std::vector<int>& output_offset) {    \
      return x.zz_diagsplit(output_offset);                             \
    }                                                                   \
    inline SWIG_FRIEND std::vector< M >                                 \
      diagsplit(const M& x, int incr=1) {                               \
      return x.zz_diagsplit(incr);                                      \
    }                                                                   \
    inline SWIG_FRIEND std::vector< M >                                 \
      diagsplit(const M& x, int incr1, int incr2) {                     \
      return x.zz_diagsplit(incr1, incr2);                              \
    }                                                                   \
    inline SWIG_FRIEND M veccat(const std::vector< M >& x) {            \
      return M::zz_veccat(x);                                           \
    }                                                                   \
    inline SWIG_FRIEND M mul(const M &X, const M &Y) {                  \
      return X.zz_mtimes(Y);                                            \
    }                                                                   \
    inline SWIG_FRIEND M mul(const M &X, const M &Y, const M &Z) {      \
      return X.zz_mtimes(Y, Z);                                         \
    }                                                                   \
    inline SWIG_FRIEND M mul(const std::vector<M> &args) {              \
      return M::zz_mul(args);                                           \
    }                                                                   \
    inline SWIG_FRIEND M transpose(const M& X) {                        \
      return X.T();                                                     \
    }                                                                   \
    inline SWIG_FRIEND M vec(const M& a) {                              \
      return a.zz_vec();                                                \
    }                                                                   \
    inline SWIG_FRIEND M vecNZ(const M& a) {                            \
      return a.zz_vecNZ();                                              \
    }                                                                   \
    inline SWIG_FRIEND M reshape(const M& a, int nrow, int ncol) {      \
      return a.zz_reshape(nrow, ncol);                                  \
    }                                                                   \
    inline SWIG_FRIEND M reshape(const M& a, std::pair<int, int> rc) {  \
      return reshape(a, rc.first, rc.second);                           \
    }                                                                   \
    inline SWIG_FRIEND M reshape(const M& a, const Sparsity& sp) {      \
      return a.zz_reshape(sp);                                          \
    }                                                                   \
    inline SWIG_FRIEND int sprank(const M& A) {                         \
      return A.zz_sprank();                                             \
    }                                                                   \
    inline SWIG_FRIEND int norm_0_mul(const M &x, const M &y) {         \
      return x.zz_norm_0_mul(y);                                        \
    }                                                                   \
    inline SWIG_FRIEND M triu(const M& a, bool includeDiagonal=true) {  \
      return a.zz_triu(includeDiagonal);                                \
    }                                                                   \
    inline SWIG_FRIEND M tril(const M& a, bool includeDiagonal=true) {  \
      return a.zz_tril(includeDiagonal);                                \
    }                                                                   \
    inline SWIG_FRIEND M kron(const M& a, const M& b) {                 \
      return a.zz_kron(b);                                              \
    }                                                                   \
    inline SWIG_FRIEND M repmat(const M &A, int n, int m=1) {           \
      return A.zz_repmat(n, m);                                         \
    }                                                                   \
    inline SWIG_FRIEND M repmat(const M &A, const std::pair<int, int>& rc) { \
      return A.zz_repmat(rc.first, rc.second);                          \
    }                                                                   \

#ifndef SWIG
SPARSITY_INTERFACE_FRIENDS(MatType)

  /** \brief Concatenate horizontally, two matrices */
  inline SWIG_FRIEND MatType horzcat(const MatType &x, const MatType &y) {
  MatType v[] = {x, y};
  return horzcat(std::vector<MatType>(v, v+2));
 }

/** \brief Concatenate horizontally, three matrices */
 inline SWIG_FRIEND MatType horzcat(const MatType &x, const MatType &y, const MatType &z) {
   MatType v[] = {x, y, z};
   return horzcat(std::vector<MatType>(v, v+3));
 }

 /** \brief Concatenate horizontally, four matrices */
 inline SWIG_FRIEND MatType horzcat(const MatType &x, const MatType &y, const MatType &z,
                                    const MatType &w) {
   MatType v[] = {x, y, z, w};
   return horzcat(std::vector<MatType>(v, v+4));
 }

 /** \brief Concatenate vertically, two matrices */
 inline SWIG_FRIEND MatType vertcat(const MatType &x, const MatType &y) {
   MatType v[] = {x, y};
   return vertcat(std::vector<MatType>(v, v+2));
 }

 /** \brief Concatenate vertically, three matrices */
 inline SWIG_FRIEND MatType vertcat(const MatType &x, const MatType &y, const MatType &z) {
   MatType v[] = {x, y, z};
   return vertcat(std::vector<MatType>(v, v+3));
 }

 /** \brief Concatenate vertically, four matrices */
 inline SWIG_FRIEND MatType vertcat(const MatType &x, const MatType &y, const MatType &z,
                                    const MatType &w) {
   MatType v[] = {x, y, z, w};
   return vertcat(std::vector<MatType>(v, v+4));
 }

 /** \brief Concatenate along diagonal, two matrices */
 inline SWIG_FRIEND MatType diagcat(const MatType &x, const MatType &y) {
   MatType v[] = {x, y};
   return diagcat(std::vector<MatType>(v, v+2));
 }

 /** \brief Concatenate along diagonal, three matrices */
 inline SWIG_FRIEND MatType diagcat(const MatType &x, const MatType &y, const MatType &z) {
   MatType v[] = {x, y, z};
   return diagcat(std::vector<MatType>(v, v+3));
 }

 /** \brief Concatenate along diagonal, four matrices */
 inline SWIG_FRIEND MatType diagcat(const MatType &x, const MatType &y, const MatType &z,
                                    const MatType &w) {
   MatType v[] = {x, y, z, w};
   return diagcat(std::vector<MatType>(v, v+4));
 }
#endif // SWIG

  };

#ifndef SWIG
  template<typename MatType>
  MatType SparsityInterface<MatType>::zz_vec() const {
    return reshape(self(), self().numel(), 1);
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::zz_repmat(int n, int m) const {
    MatType allrows = vertcat(std::vector<MatType>(n, self()));
    return horzcat(std::vector<MatType>(m, allrows));
  }

  template<typename MatType>
  std::vector< std::vector< MatType > >
  SparsityInterface<MatType>::zz_blocksplit(const std::vector<int>& vert_offset,
                                            const std::vector<int>& horz_offset) const {
    std::vector< MatType > rows = vertsplit(self(), vert_offset);
    std::vector< std::vector< MatType > > ret;
    for (int i=0;i<rows.size();++i) {
      ret.push_back(horzsplit(rows[i], horz_offset));
    }
    return ret;
  }

  template<typename MatType>
  std::vector< std::vector< MatType > >
  SparsityInterface<MatType>::zz_blocksplit(int vert_incr, int horz_incr) const {
    casadi_assert(horz_incr>=1);
    casadi_assert(vert_incr>=1);
    int sz1 = self().size1();
    std::vector<int> offset1 = range(0, sz1, vert_incr);
    offset1.push_back(sz1);
    int sz2 = self().size2();
    std::vector<int> offset2 = range(0, sz2, horz_incr);
    offset2.push_back(sz2);
    return blocksplit(self(), offset1, offset2);
  }

  template<typename MatType>
  std::vector<int>
  SparsityInterface<MatType>::zz_offset(const std::vector< MatType > &v, bool vert) {
    std::vector<int> ret(v.size()+1);
    ret[0]=0;
    for (int i=0; i<v.size(); ++i) {
      ret[i+1] = ret[i] + (vert ? v[i].size1() : v[i].size2());
    }
    return ret;
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::zz_veccat(const std::vector< MatType >& x) {
    std::vector< MatType > x_vec = x;
    for (typename std::vector< MatType >::iterator it=x_vec.begin();
         it!=x_vec.end(); ++it) {
      *it = vec(*it);
    }
    return vertcat(x_vec);
  }

  template<typename MatType>
  std::vector< MatType >
  SparsityInterface<MatType>::zz_diagsplit(const std::vector<int>& output_offset) const {
    casadi_assert_message(self().issquare(), "diagsplit(x,incr)::input must be square but got "
                          << self().dimString()  << ".");
    return diagsplit(self(), output_offset, output_offset);
  }

  template<typename MatType>
  std::vector< MatType >
  SparsityInterface<MatType>::zz_diagsplit(int incr) const {
    casadi_assert(incr>=1);
    casadi_assert_message(self().issquare(), "diagsplit(x,incr)::input must be square but got "
                          << self().dimString()  << ".");
    std::vector<int> offset2 = range(0, self().size2(), incr);
    offset2.push_back(self().size2());
    return diagsplit(self(), offset2);
  }

  template<typename MatType>
  std::vector< MatType >
  SparsityInterface<MatType>::zz_diagsplit(int incr1, int incr2) const {
    casadi_assert(incr1>=1);
    casadi_assert(incr2>=1);
    std::vector<int> offset1 = range(0, self().size1(), incr1);
    offset1.push_back(self().size1());
    std::vector<int> offset2 = range(0, self().size2(), incr2);
    offset2.push_back(self().size2());
    return diagsplit(self(), offset1, offset2);
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::zz_mul(const std::vector<MatType> &args) {
    casadi_assert_message(args.size()>=1,
                          "mul(std::vector<MatType> &args): "
                          "supplied list must not be empty.");
    MatType ret = args[0];
    for (int i=1; i<args.size(); ++i) ret = mul(ret, args[i]);
    return ret;
  }

  template<typename MatType>
  std::vector<MatType > SparsityInterface<MatType>::zz_horzsplit(int incr) const {
    casadi_assert(incr>=1);
    int sz2 = self().size2();
    std::vector<int> offset2 = range(0, sz2, incr);
    offset2.push_back(sz2);
    return horzsplit(self(), offset2);
  }

  template<typename MatType>
  std::vector<MatType > SparsityInterface<MatType>::zz_vertsplit(int incr) const {
    casadi_assert(incr>=1);
    int sz1 = self().size1();
    std::vector<int> offset1 = range(0, sz1, incr);
    offset1.push_back(sz1);
    return vertsplit(self(), offset1);
  }
#endif // SWIG

} // namespace casadi

#endif // CASADI_SPARSITY_INTERFACE_HPP
