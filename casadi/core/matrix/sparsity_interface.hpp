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
  public:
#ifndef SWIG
    /** \brief Concatenate a list of matrices horizontally
     * Alternative terminology: horizontal stack, hstack, horizontal append, [a b]
     *
     *   horzcat(horzsplit(x, ...)) = x
     */
    inline friend MatType horzcat(const std::vector<MatType> &v) {
      return MatType::zz_horzcat(v);
    }

    /** \brief Concatenate horizontally, two matrices */
    inline friend MatType horzcat(const MatType &x, const MatType &y) {
      std::vector<MatType> v(2);
      v[0] = x;
      v[1] = y;
      return horzcat(v);
    }

    /** \brief Concatenate a list of matrices vertically
     * Alternative terminology: vertical stack, vstack, vertical append, [a;b]
     *
     *   vertcat(vertsplit(x, ...)) = x
     */
    inline friend MatType vertcat(const std::vector<MatType> &v) {
      return MatType::zz_vertcat(v);
    }

    /** \brief Concatenate vertically, two matrices */
    inline friend MatType vertcat(const MatType &x, const MatType &y) {
      std::vector<MatType> v(2);
      v[0] = x;
      v[1] = y;
      return vertcat(v);
    }

    /** \brief Construct a matrix with given block on the diagonal */
    inline friend MatType blkdiag(const std::vector<MatType> &A) {
      return MatType::zz_blkdiag(A);
    }

    /** \brief Concatenate along diagonal, two matrices */
    inline friend MatType blkdiag(const MatType &x, const MatType &y) {
      std::vector<MatType> v(2);
      v[0] = x;
      v[1] = y;
      return blkdiag(v);
    }

    /** \brief Matrix product of two matrices:  */
    inline friend MatType mul(const MatType &X, const MatType &Y) {
      return X.zz_mtimes(Y);
    }

    /** \brief Matrix product and addition
        Matrix product of two matrices (X and Y), adding the result to
        a third matrix Z. The result has the same sparsity pattern as
        C meaning that other entries of (X*Y) are ignored.
        The operation is equivalent to: Z+mul(X,Y).setSparse(Z.sparsity()).
    */
    inline friend MatType mul(const MatType &X, const MatType &Y, const MatType &Z) {
      return X.zz_mtimes(Y, Z);
    }

    /** \brief Matrix product of n matrices */
    inline friend MatType mul(const std::vector<MatType> &args) {
      casadi_assert_message(args.size()>=1,
                            "mul(std::vector<MatType> &args): "
                            "supplied list must not be empty.");
      MatType ret = args[0];
      for (int i=1; i<args.size(); ++i) ret = mul(ret, args[i]);
      return ret;
    }

    /** \brief Transpose */
    inline friend MatType transpose(const MatType& X) {
      return X.T();
    }
#endif // SWIG
  };

} // namespace casadi

#endif // CASADI_SPARSITY_INTERFACE_HPP
