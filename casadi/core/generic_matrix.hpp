/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#ifndef CASADI_GENERIC_MATRIX_HPP
#define CASADI_GENERIC_MATRIX_HPP

#include "slice.hpp"
#include "submatrix.hpp"
#include "nonzeros.hpp"
#include "sparsity.hpp"
#include "calculus.hpp"
#include "sparsity_interface.hpp"
#include "generic_type.hpp"

namespace casadi {
  /** \brief Empty Base

      This class is extended in SWIG.

      \identifier{1al} */
  struct CASADI_EXPORT GenericMatrixCommon {};

  /** \brief Matrix base class

      This is a common base class for MX and Matrix<>, introducing a uniform syntax and implementing
      common functionality using the curiously recurring template pattern (CRTP) idiom.\n

      The class is designed with the idea that "everything is a matrix",
      that is, also scalars and vectors.\n
      This philosophy makes it easy to use and to interface in particularly
      with Python and Matlab/Octave.\n

      The syntax tries to stay as close as possible to the ublas syntax
      when it comes to vector/matrix operations.\n

      Index starts with 0.\n
      Index vec happens as follows: (rr, cc) -> k = rr+cc*size1()\n
      Vectors are column vectors.\n

      The storage format is Compressed Column Storage (CCS),
      similar to that used for sparse matrices in Matlab, \n
      but unlike this format, we do allow for elements to be structurally non-zero
      but numerically zero.\n

      The sparsity pattern, which is reference counted and cached,
      can be accessed with Sparsity& sparsity()\n

      \author Joel Andersson
      \date 2012

      \identifier{1am} */
  template<typename MatType>
  class GenericMatrix
    : public GenericMatrixCommon,
      public SWIG_IF_ELSE(SparsityInterfaceCommon, SparsityInterface<MatType>) {
    using SparsityInterface<MatType>::self;
  public:

    /** \brief Get the number of structural nonzeros

        Counts the entries present in the sparsity pattern. Entries that are stored but
        numerically zero still count; structural zeros (not stored) do not.
        For a dense matrix this equals numel().

        \lang{python}
        \doctest
        >>> DM.eye(3).nnz(), DM.ones(3, 3).nnz(), DM(3, 3).nnz()
        (3, 9, 0)
        >>> DM([[1, 0], [0, 2]]).nnz()
        4
        >>> MX.sym("x", 2, 0).nnz()
        0
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%d %d %d\n', DM.eye(3).nnz(), DM.ones(3, 3).nnz(), DM(3, 3).nnz())
        3 9 0
        >> DM([1 0; 0 2]).nnz()
        4
        >> MX.sym('x', 2, 0).nnz()
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.nnz(ca.DM.eye(3)), ca.nnz(ca.DM.ones(3, 3)), ca.nnz(ca.DM(3, 3))
        (3, 9, 0)
        julia> ca.nnz(ca.DM([1 0; 0 2]))
        4
        julia> ca.nnz(ca.MX.sym("x", 2, 0))
        0
        \enddoctest
        \endlang

        \sa numel, nnz_lower, nnz_upper, nnz_diag, sparsity

        \identifier{1an} */
    casadi_int nnz() const;

    /** \brief Get the number of structural nonzeros on or below the diagonal

        \lang{python}
        \doctest
        >>> DM.ones(3, 3).nnz_lower(), DM.eye(3).nnz_lower(), DM.ones(2, 3).nnz_lower()
        (6, 3, 3)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM.ones(3, 3); I = DM.eye(3); B = DM.ones(2, 3);
        >> fprintf('%d %d %d\n', A.nnz_lower(), I.nnz_lower(), B.nnz_lower())
        6 3 3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM.ones(3, 3); I = ca.DM.eye(3); B = ca.DM.ones(2, 3);
        julia> ca.nnz_lower(A), ca.nnz_lower(I), ca.nnz_lower(B)
        (6, 3, 3)
        \enddoctest
        \endlang

        \sa nnz, nnz_upper, nnz_diag

        \identifier{1ao} */
    casadi_int nnz_lower() const;

    /** \brief Get the number of structural nonzeros on or above the diagonal

        \lang{python}
        \doctest
        >>> DM.ones(3, 3).nnz_upper(), DM.eye(3).nnz_upper(), DM.ones(2, 3).nnz_upper()
        (6, 3, 5)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM.ones(3, 3); I = DM.eye(3); B = DM.ones(2, 3);
        >> fprintf('%d %d %d\n', A.nnz_upper(), I.nnz_upper(), B.nnz_upper())
        6 3 5
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM.ones(3, 3); I = ca.DM.eye(3); B = ca.DM.ones(2, 3);
        julia> ca.nnz_upper(A), ca.nnz_upper(I), ca.nnz_upper(B)
        (6, 3, 5)
        \enddoctest
        \endlang

        \sa nnz, nnz_lower, nnz_diag

        \identifier{1ap} */
    casadi_int nnz_upper() const;

    /** \brief Get the number of structural nonzeros on the diagonal

        Also defined for non-square matrices.

        \lang{python}
        \doctest
        >>> DM.ones(3, 3).nnz_diag(), DM.eye(3).nnz_diag(), DM.ones(2, 3).nnz_diag()
        (3, 3, 2)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM.ones(3, 3); I = DM.eye(3); B = DM.ones(2, 3);
        >> fprintf('%d %d %d\n', A.nnz_diag(), I.nnz_diag(), B.nnz_diag())
        3 3 2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM.ones(3, 3); I = ca.DM.eye(3); B = ca.DM.ones(2, 3);
        julia> ca.nnz_diag(A), ca.nnz_diag(I), ca.nnz_diag(B)
        (3, 3, 2)
        \enddoctest
        \endlang

        \sa nnz, nnz_lower, nnz_upper

        \identifier{1aq} */
    casadi_int nnz_diag() const;

    /** \brief Get the number of elements, size1()*size2()

        Counts all entries, structural zeros included.

        \lang{python}
        \doctest
        >>> DM.eye(3).numel(), MX.sym("x", 2, 0).numel(), SX.sym("x").numel()
        (9, 0, 1)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%d %d %d\n', DM.eye(3).numel(), MX.sym('x', 2, 0).numel(), SX.sym('x').numel())
        9 0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.numel(ca.DM.eye(3)), ca.numel(ca.MX.sym("x", 2, 0)), ca.numel(ca.SX.sym("x"))
        (9, 0, 1)
        julia> length(ca.DM.eye(3))
        9
        \enddoctest
        \endlang

        \sa nnz, size

        \identifier{1ar} */
    casadi_int numel() const;

    /** \brief Get the first dimension (i.e. number of rows)

        \lang{python}
        \doctest
        >>> x = MX.sym("x", 2, 3)
        >>> x.size1(), x.size2()
        (2, 3)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = MX.sym('x', 2, 3);
        >> fprintf('%d %d\n', x.size1(), x.size2())
        2 3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.MX.sym("x", 2, 3);
        julia> ca.size1(x), ca.size2(x)
        (2, 3)
        \enddoctest
        \endlang

        \sa size2, size, rows

        \identifier{1as} */
    casadi_int size1() const;

    /** \brief Get the number of rows, Octave-style syntax (same as size1())

        \identifier{1at} */
    casadi_int rows() const {return size1();}

    /** \brief Get the second dimension (i.e. number of columns)

        \sa size1, size, columns

        \identifier{1au} */
    casadi_int size2() const;

    /** \brief Get the number of columns, Octave-style syntax (same as size2())

        \identifier{1av} */
    casadi_int columns() const {return size2();}

    /** \brief Get string representation of dimensions.

        The representation is e.g. "4x5" or, with \p with_nz true, "4x5,10nz".

        \lang{python}
        \doctest
        >>> DM.eye(3).dim(), DM.eye(3).dim(True), DM.ones(3, 3).dim(True)
        ('3x3', '3x3,3nz', '3x3,9nz')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%s %s %s\n', DM.eye(3).dim(), DM.eye(3).dim(true), DM.ones(3, 3).dim(true))
        3x3 3x3,3nz 3x3,9nz
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.dim(ca.DM.eye(3)), ca.dim(ca.DM.eye(3), true), ca.dim(ca.DM.ones(3, 3), true)
        ("3x3", "3x3,3nz", "3x3,9nz")
        \enddoctest
        \endlang

        \sa size, nnz

        \identifier{1aw} */
    std::string dim(bool with_nz=false) const;

    /** \brief  Get the shape as a pair (number of rows, number of columns)

        \lang{python}
        \doctest
        >>> x = MX.sym("x", 2, 3)
        >>> x.size(), x.size(1), x.size(2)
        ((2, 3), 2, 3)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = MX.sym('x', 2, 3);
        >> size(x)
        2 3
        >> fprintf('%d %d\n', size(x, 1), size(x, 2))
        2 3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.MX.sym("x", 2, 3);
        julia> size(x), size(x, 1), size(x, 2)
        ((2, 3), 2, 3)
        \enddoctest
        \endlang

        \sa size1, size2, numel, dim

        \identifier{1ax} */
    std::pair<casadi_int, casadi_int> size() const;

    /** \brief  Get the size along a particular dimension

        \param axis 1 for the number of rows, 2 for the number of columns; any other value
        raises an error

        \identifier{1ay} */
    casadi_int size(casadi_int axis) const;

    /** \brief Check if the matrix is empty, i.e. if one of the dimensions is zero

        With \p both true, only a 0-by-0 matrix counts as empty.

        \lang{python}
        \doctest
        >>> DM(0, 3).is_empty(), DM(0, 3).is_empty(True), DM(0, 0).is_empty(True)
        (True, False, True)
        >>> DM(1, 1).is_empty()
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> E = DM(0, 3);
        >> fprintf('%d %d %d\n', E.is_empty(), E.is_empty(true), DM(0, 0).is_empty(true))
        1 0 1
        >> DM(1, 1).is_empty()
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> E = ca.DM(0, 3);
        julia> ca.is_empty(E), ca.is_empty(E, true), ca.is_empty(ca.DM(0, 0), true)
        (true, false, true)
        julia> ca.is_empty(ca.DM(1, 1))
        false
        \enddoctest
        \endlang

        \identifier{1az} */
    bool is_empty(bool both=false) const { return sparsity().is_empty(both);}

    /** \brief  Check if the matrix is dense, i.e. nnz() equals numel()

        Based on the sparsity pattern only: stored numerical zeros still count as nonzeros.
        Empty matrices are dense.

        \lang{python}
        \doctest
        >>> DM.zeros(2, 2).is_dense(), DM(2, 2).is_dense(), DM(0, 3).is_dense()
        (True, False, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Z = DM.zeros(2, 2);
        >> fprintf('%d %d %d\n', Z.is_dense(), DM(2, 2).is_dense(), DM(0, 3).is_dense())
        1 0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_dense(ca.DM.zeros(2, 2)), ca.is_dense(ca.DM(2, 2)), ca.is_dense(ca.DM(0, 3))
        (true, false, true)
        \enddoctest
        \endlang

        \sa nnz, numel, densify

        \identifier{1b0} */
    bool is_dense() const { return sparsity().is_dense();}

    /** \brief  Check if the matrix is scalar, i.e. 1-by-1

        \param scalar_and_dense When true, a 1-by-1 matrix with a structural zero does not count

        \lang{python}
        \doctest
        >>> x = SX.sym("x", Sparsity(1, 1))
        >>> x.is_scalar(), x.is_scalar(True), SX.sym("y").is_scalar(True)
        (True, False, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x', Sparsity(1, 1));
        >> fprintf('%d %d %d\n', x.is_scalar(), x.is_scalar(true), SX.sym('y').is_scalar(true))
        1 0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x", ca.Sparsity(1, 1));
        julia> ca.is_scalar(x), ca.is_scalar(x, true), ca.is_scalar(ca.SX.sym("y"), true)
        (true, false, true)
        \enddoctest
        \endlang

        \identifier{1b1} */
    bool is_scalar(bool scalar_and_dense=false) const;

    /** \brief  Check if the matrix is square, i.e. size1() equals size2()

        \identifier{1b2} */
    bool is_square() const { return sparsity().is_square();}

    /** \brief  Check if the matrix is a row or column vector

        True when size1() or size2() is one, so scalars and n-by-1 or 1-by-n empties count.
        A 0-by-0 matrix does not.

        \lang{python}
        \doctest
        >>> DM(3, 1).is_vector(), DM(1, 3).is_vector(), DM(1, 1).is_vector()
        (True, True, True)
        >>> DM(3, 3).is_vector(), DM(0, 1).is_vector(), DM(0, 0).is_vector()
        (False, True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%d %d %d\n', DM(3, 1).is_vector(), DM(1, 3).is_vector(), DM(1, 1).is_vector())
        1 1 1
        >> fprintf('%d %d %d\n', DM(3, 3).is_vector(), DM(0, 1).is_vector(), DM(0, 0).is_vector())
        0 1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_vector(ca.DM(3, 1)), ca.is_vector(ca.DM(1, 3)), ca.is_vector(ca.DM(1, 1))
        (true, true, true)
        julia> ca.is_vector(ca.DM(3, 3)), ca.is_vector(ca.DM(0, 1)), ca.is_vector(ca.DM(0, 0))
        (false, true, false)
        \enddoctest
        \endlang

        \sa is_row, is_column

        \identifier{1b3} */
    bool is_vector() const { return sparsity().is_vector();}

    /** \brief  Check if the matrix is a row vector (i.e. size1()==1)

        \identifier{1b4} */
    bool is_row() const { return sparsity().is_row();}

    /** \brief  Check if the matrix is a column vector (i.e. size2()==1)

        \identifier{1b5} */
    bool is_column() const { return sparsity().is_column();}

    /** \brief Check if the matrix is upper triangular

        Decided on the sparsity pattern alone: all structural nonzeros must lie on or above
        the diagonal. A dense matrix with numerical zeros below the diagonal is not upper
        triangular. Also defined for non-square matrices.

        \lang{python}
        \doctest
        >>> DM(Sparsity.upper(3), 1).is_triu(), DM.zeros(2, 2).is_triu(), DM(2, 2).is_triu()
        (True, False, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> U = DM(Sparsity.upper(3), 1);
        >> fprintf('%d %d %d\n', U.is_triu(), DM.zeros(2, 2).is_triu(), DM(2, 2).is_triu())
        1 0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> U = ca.DM(ca.Sparsity.upper(3), 1);
        julia> ca.is_triu(U), ca.is_triu(ca.DM.zeros(2, 2)), ca.is_triu(ca.DM(2, 2))
        (true, false, true)
        \enddoctest
        \endlang

        \sa is_tril, triu, triu2symm

        \identifier{1b6} */
    bool is_triu() const { return sparsity().is_triu();}

    /** \brief Check if the matrix is lower triangular

        Decided on the sparsity pattern alone, see is_triu.

        \sa is_triu, tril, tril2symm

        \identifier{1b7} */
    bool is_tril() const { return sparsity().is_tril();}

    ///@{
    /** \brief Get the sparsity pattern. See the Sparsity class for details.

        The pattern is stored in compressed column format: colind has size2()+1 entries,
        the nonzeros of column j are the entries colind(j) up to colind(j+1) (exclusive),
        and row(el) is the row index of nonzero \p el. Without argument, the full vectors
        are returned. Out-of-range indices raise an error.

        \lang{python}
        \doctest
        >>> x = DM(Sparsity.diag(2), [5, 6])
        >>> x.row(), x.colind()
        ([0, 1], [0, 1, 2])
        >>> x.row(1), x.colind(2)
        (1, 2)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = DM(Sparsity.diag(2), [5; 6]);
        >> x.row(), x.colind()
        0 1
        0 1 2
        >> fprintf('%d %d\n', x.row(1), x.colind(2))
        1 2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.DM(ca.Sparsity.diag(2), [5, 6]);
        julia> ca.row(x), ca.colind(x)
        ([0, 1], [0, 1, 2])
        julia> ca.row(x, 1), ca.colind(x, 2)
        (1, 2)
        \enddoctest
        \endlang

        \identifier{1b8} */
    std::vector<casadi_int> get_row() const { return sparsity().get_row(); }
    std::vector<casadi_int> get_colind() const { return sparsity().get_colind(); }
#ifndef SWIG
    const casadi_int* row() const { return sparsity().row(); }
    const casadi_int* colind() const { return sparsity().colind(); }
#endif
    casadi_int row(casadi_int el) const { return sparsity().row(el); }
    casadi_int colind(casadi_int col) const { return sparsity().colind(col); }
    ///@}

    /** \brief Get the sparsity pattern

        \lang{python}
        \doctest
        >>> DM.eye(2).sparsity()
        Sparsity(2x2,2nz)
        >>> DM.eye(2).sparsity() == Sparsity.diag(2)
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM.eye(2).sparsity()
        2x2,2nz
        >> DM.eye(2).sparsity() == Sparsity.diag(2)
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.sparsity(ca.DM.eye(2))
        Sparsity(2x2,2nz)
        julia> ca.is_equal(ca.sparsity(ca.DM.eye(2)), ca.Sparsity.diag(2))
        true
        \enddoctest
        \endlang

        \identifier{1b9} */
    SWIG_CONSTREF(Sparsity) sparsity() const;

#ifndef SWIG
    /// \cond CLUTTER
    /**  @{  */
    /// Functions called by friend functions defined here
    static MatType interp1d(const std::vector<double>& x, const MatType &v,
         const std::vector<double>& xq, const std::string& mode, bool equidistant);
    static casadi_int sprank(const MatType &x) { return Sparsity::sprank(x.sparsity());}
    static casadi_int norm_0_mul(const MatType &x, const MatType &y) {
      return Sparsity::norm_0_mul(x.sparsity(), y.sparsity());
    }
    static MatType tril(const MatType &x, bool includeDiagonal=true) {
      return project(x, Sparsity::tril(x.sparsity(), includeDiagonal));
    }
    static MatType triu(const MatType &x, bool includeDiagonal=true) {
      return project(x, Sparsity::triu(x.sparsity(), includeDiagonal));
    }
    static MatType sumsqr(const MatType &x) { return dot(x, x);}
    static MatType linspace(const MatType &a, const MatType &b, casadi_int nsteps);
    static MatType cross(const MatType &a, const MatType &b, casadi_int dim=-1);
    static MatType skew(const MatType &a);
    static MatType inv_skew(const MatType &a);
    static MatType tril2symm(const MatType &x);
    static MatType triu2symm(const MatType &x);
    static MatType repsum(const MatType &x, casadi_int n, casadi_int m=1);
    static MatType diff(const MatType &x, casadi_int n=1, casadi_int axis=-1);

    static bool is_linear(const MatType &expr, const MatType &var);
    static bool is_quadratic(const MatType &expr, const MatType &var);
    static void quadratic_coeff(const MatType &expr, const MatType &var,
        MatType& A, MatType& b, MatType& c, bool check);
    static void linear_coeff(const MatType &expr, const MatType &var,
        MatType& A, MatType& b, bool check);
    /** @}  */
    /// \endcond

    /** \brief  Get vector nonzero or slice of nonzeros

        \identifier{1bb} */
    template<typename K>
    const MatType nz(const K& k) const {
      MatType ret;
      self().get_nz(ret, false, k);
      return ret;
    }

    /** \brief  Access vector nonzero or slice of nonzeros

        \identifier{1bc} */
    template<typename K>
    NonZeros<MatType, K> nz(const K& k) {
      return NonZeros<MatType, K>(self(), k);
    }

    /** \brief  Get vector element or slice

        \identifier{1bd} */
    template<typename RR>
    const MatType operator()(const RR& rr) const {
      MatType ret;
      self().get(ret, false, rr);
      return ret;
    }

    /** \brief  Get Matrix element or slice

        \identifier{1be} */
    template<typename RR, typename CC>
    const MatType operator()(const RR& rr, const CC& cc) const {
      MatType ret;
      self().get(ret, false, rr, cc);
      return ret;
    }

    /** \brief Access Matrix elements (one argument)

        \identifier{1bf} */
    template<typename RR>
    SubIndex<MatType, RR> operator()(const RR& rr) {
      return SubIndex<MatType, RR>(self(), rr);
    }

    /** \brief Access Matrix elements (two arguments)

        \identifier{1bg} */
    template<typename RR, typename CC>
    SubMatrix<MatType, RR, CC> operator()(const RR& rr, const CC& cc) {
      return SubMatrix<MatType, RR, CC>(self(), rr, cc);
    }
#endif // SWIG

#if !defined(SWIG) || defined(DOXYGEN)
/**
\addtogroup expression_tools
@{
*/

     /** \brief Performs 1d linear interpolation

         The data points are (x[i], v[i, :]): row i of \p v holds the values at abscissa x[i],
         so \p v has as many rows as \p x has entries and every column is interpolated
         independently. The result has one row per query point xq[j].

         Queries outside the range of \p x are clamped to the first or last data point
         (no extrapolation). The result is a constant sparse matrix times \p v, so for
         symbolic \p v the expression is linear in \p v.

         \param x Abscissae, strictly increasing, at least two entries
         \param v Data values, x.size()-by-m
         \param xq Query abscissae, any order
         \param mode "linear" (interpolate), "floor" (value at the data point at or below the
         query) or "ceil" (at or above); any other string raises an error
         \param equidistant When true, \p x is assumed uniformly spaced and the interval is
         located arithmetically from x[0] and x[1] instead of by search
         \return Matrix of size xq.size()-by-m

         \lang{python}
         \doctest
         >>> interp1d([0, 1, 2], DM([0, 10, 20]), [0.5, 1.5])
         DM([5, 15])
         >>> interp1d([0, 1, 2], DM([0, 10, 20]), [0.5, 1.5], "floor")
         DM([0, 10])
         >>> interp1d([0, 1, 2], DM([[0, 1], [10, 11], [20, 21]]), [0.5, 2])
         DM(
         [[5, 6],
          [20, 21]])
         >>> interp1d([0, 1, 2], DM([0, 10, 20]), [-1, 3])
         DM([0, 20])
         >>> interp1d([0, 1, 2], SX.sym("v", 3), [0.5])
         SX(@1=0.5, ((@1*v_0)+(@1*v_1)))
         \enddoctest
         \endlang

         \lang{matlab}
         \doctest
         >> interp1d([0, 1, 2], DM([0; 10; 20]), [0.5, 1.5])
         [5, 15]
         >> interp1d([0, 1, 2], DM([0; 10; 20]), [0.5, 1.5], 'floor')
         [0, 10]
         >> interp1d([0, 1, 2], DM([0 1; 10 11; 20 21]), [0.5, 2])
         [[5, 6],
          [20, 21]]
         >> interp1d([0, 1, 2], DM([0; 10; 20]), [-1, 3])
         [0, 20]
         >> interp1d([0, 1, 2], SX.sym('v', 3), 0.5)
         @1=0.5, ((@1*v_0)+(@1*v_1))
         \enddoctest
         \endlang

         \lang{julia}
         \doctest
         julia> ca.interp1d([0, 1, 2], ca.DM([0, 10, 20]), [0.5, 1.5])
         DM([5, 15])
         julia> ca.interp1d([0, 1, 2], ca.DM([0, 10, 20]), [0.5, 1.5], "floor")
         DM([0, 10])
         julia> ca.interp1d([0, 1, 2], ca.DM([0 1; 10 11; 20 21]), [0.5, 2])
         DM(
         [[5, 6],
          [20, 21]])
         julia> ca.interp1d([0, 1, 2], ca.DM([0, 10, 20]), [-1, 3])
         DM([0, 20])
         julia> ca.interp1d([0, 1, 2], ca.SX.sym("v", 3), [0.5])
         SX(@1=0.5, ((@1*v_0)+(@1*v_1)))
         \enddoctest
         \endlang

         \sa interpolant

         \identifier{1bh} */
     inline friend MatType interp1d(const std::vector<double>& x, const MatType&v,
         const std::vector<double>& xq, const std::string& mode, bool equidistant=false) {
       return MatType::interp1d(x, v, xq, mode, equidistant);
     }

    /** \brief Matrix power x^n

        For a square matrix \p x and a constant integer scalar \p n, computes the repeated
        matrix product by squaring. n == 0 gives the (sparse) identity, negative \p n gives
        the inverse of the positive power (via inv). If both arguments are scalar, this is
        the elementwise pow. A non-integer or symbolic \p n for a matrix raises an error.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> mpower(A, 2)
        DM(
        [[7, 10],
         [15, 22]])
        >>> mpower(A, 0)
        DM(
        [[1, 00],
         [00, 1]])
        >>> mpower(A, -1)
        DM(
        [[-2, 1],
         [1.5, -0.5]])
        >>> mpower(DM(2), 0.5)
        DM(1.41421)
        >>> try: mpower(A, 0.5)
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> A^2
        [[7, 10],
         [15, 22]]
        >> A^0
        [[1, 00],
         [00, 1]]
        >> mpower(A, -1)
        [[-2, 1],
         [1.5, -0.5]]
        >> DM(2)^0.5
        1.41421
        >> try, A^0.5; catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> A^2
        DM(
        [[7, 10],
         [15, 22]])
        julia> A^0
        DM(
        [[1, 00],
         [00, 1]])
        julia> A^-1
        DM(
        [[-2, 1],
         [1.5, -0.5]])
        julia> ca.DM(2)^0.5
        DM(1.41421)
        julia> try A^0.5 catch e; println("error") end
        error
        \enddoctest
        \endlang

        \sa mtimes, inv, expm

        \identifier{1bi} */
    inline friend MatType mpower(const MatType& x, const MatType& n) {
      return MatType::mpower(x, n);
    }

    /** \brief Construct second-order-cone matrix

        soc(x, y) computes the (n+1)-by-(n+1) matrix [y*eye(n) x; x' y].
        This matrix is positive semidefinite if and only if norm_2(x) <= y, which makes it
        a way to express a second-order cone constraint as a semidefinite constraint.

        \param x Vector expression of size n (row or column; sparse entries are densified)
        \param y Scalar expression

        \lang{python}
        \doctest
        >>> soc(DM([1, 2]), 3)
        DM(
        [[3, 00, 1],
         [00, 3, 2],
         [1, 2, 3]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> soc(DM([1; 2]), 3)
        [[3, 00, 1],
         [00, 3, 2],
         [1, 2, 3]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.soc(ca.DM([1, 2]), 3)
        DM(
        [[3, 00, 1],
         [00, 3, 2],
         [1, 2, 3]])
        \enddoctest
        \endlang

        \sa norm_2

        \identifier{1bj} */
    inline friend MatType soc(const MatType& x, const MatType& y) {
      return MatType::soc(x, y);
    }

    /** \brief Compute any contraction of two dense tensors, using index/einstein notation

        Given two tensors, A and B, computes a third tensor C such that:

        C_c = A_a * B_b

        With a, b, c representing einstein indices.
        Instead of the classical index labels i,j,k,... we employ -1,-2,-3,...
        An index that appears in a or b but not in c is summed over.

        A, B, C are represented as CasADi column vectors (column-major flattening), with
        dim_a, dim_b, dim_c indicating their tensorial dimensions. The product of dim_a must
        equal the number of elements of A, and likewise for B and C.

        In the overload with an explicit \p C argument, the contraction is added to \p C.

        \lang{python}
        \doctest
        >>> A = DM([1, 2, 3, 4])
        >>> einstein(A, A, [2, 2], [2, 2], [2, 2], [-1, -2], [-2, -3], [-1, -3])
        DM([7, 10, 15, 22])
        >>> mtimes(DM([[1, 3], [2, 4]]), DM([[1, 3], [2, 4]]))
        DM(
        [[7, 15],
         [10, 22]])
        >>> einstein(DM([1, 2, 3]), DM([4, 5, 6]), [3], [3], [], [-1], [-1], [])
        DM(32)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1; 2; 3; 4]);
        >> einstein(A, A, [2, 2], [2, 2], [2, 2], [-1, -2], [-2, -3], [-1, -3])
        [7, 10, 15, 22]
        >> DM([1 3; 2 4])*DM([1 3; 2 4])
        [[7, 15],
         [10, 22]]
        >> einstein(DM([1; 2; 3]), DM([4; 5; 6]), 3, 3, [], -1, -1, [])
        32
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1, 2, 3, 4]);
        julia> ca.einstein(A, A, [2, 2], [2, 2], [2, 2], [-1, -2], [-2, -3], [-1, -3])
        DM([7, 10, 15, 22])
        julia> ca.DM([1 3; 2 4]) * ca.DM([1 3; 2 4])
        DM(
        [[7, 15],
         [10, 22]])
        julia> ca.einstein(ca.DM([1, 2, 3]), ca.DM([4, 5, 6]), [3], [3], Int[], [-1], [-1], Int[])
        DM(32)
        \enddoctest
        \endlang

        \sa mtimes, dot

        \identifier{1bk} */
    ///@{
    inline friend MatType
      einstein(const MatType &A, const MatType &B, const MatType &C,
        const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
        const std::vector<casadi_int>& dim_c,
        const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
        const std::vector<casadi_int>& c) {
      return MatType::einstein(A, B, C, dim_a, dim_b, dim_c, a, b, c);
    }

    inline friend MatType
      einstein(const MatType &A, const MatType &B,
        const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
        const std::vector<casadi_int>& dim_c,
        const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
        const std::vector<casadi_int>& c) {
      return MatType::einstein(A, B, dim_a, dim_b, dim_c, a, b, c);
    }
    ///@}

    /** \brief Matrix right divide, x/A (cf. slash '/' in MATLAB)

        Computes x*inv(A) as solve(A', x')'. If either argument is scalar, the
        elementwise division x/A is returned instead.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> mrdivide(DM([[1, 2]]), A)
        DM([[1, 0]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> DM([1 2])/A
        [[1, 0]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> ca.DM([1 2]) / A
        DM([[1, 0]])
        \enddoctest
        \endlang

        \sa mldivide, solve, inv

        \identifier{1bl} */
    inline friend MatType mrdivide(const MatType& x, const MatType& n) {
      return MatType::mrdivide(x, n);
    }

    /** \brief Matrix left divide, A\\b (cf. backslash '\' in MATLAB)

        Computes inv(A)*b as solve(A, b). If either argument is scalar, the
        elementwise division b/A is returned instead.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> mldivide(A, DM([1, 2]))
        DM([0, 0.5])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> A\DM([1; 2])
        [0, 0.5]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> A \ ca.DM([1, 2])
        DM([0, 0.5])
        \enddoctest
        \endlang

        \sa mrdivide, solve, inv

        \identifier{1bm} */
    inline friend MatType mldivide(const MatType& x, const MatType& n) {
      return MatType::mldivide(x, n);
    }

    /** \brief Get all symbolic primitives present in an expression

        Returns the symbols the expression was built from, in order of first appearance,
        each once. For SX every scalar symbol is a separate entry; for MX the entries are
        the symbolic matrices as created by MX.sym, even if only part of such a matrix
        is used. Presence in the graph is a structural notion; an expression like 0*x
        still lists x. Use depends_on for mathematical dependency. Not defined for DM.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y"); z = SX.sym("z")
        >>> symvar(vertcat(z*x, y))
        [SX(z), SX(x), SX(y)]
        >>> symvar(SX.sym("v", 2)*3)
        [SX(v_0), SX(v_1)]
        >>> v = MX.sym("v", 2)
        >>> symvar(v[0]*3)
        [MX(v)]
        >>> symvar(SX(3))
        []
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z');
        >> v = symvar(vertcat(z*x, y));
        >> numel(v), vertcat(v{:})
        3
        [z, x, y]
        >> v = symvar(SX.sym('v', 2)*3); disp(vertcat(v{:}))
        [v_0, v_1]
        >> v = MX.sym('v', 2);
        >> s = symvar(v(1)*3); disp(numel(s)), disp(s{1})
        1
        v
        >> numel(symvar(SX(3)))
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y"); z = ca.SX.sym("z");
        julia> ca.symvar([z*x; y])
        3-element Vector{SX}:
         SX(z)
         SX(x)
         SX(y)
        julia> ca.symvar(ca.SX.sym("v", 2) * 3)
        2-element Vector{SX}:
         SX(v_0)
         SX(v_1)
        julia> v = ca.MX.sym("v", 2);
        julia> ca.symvar(v[1] * 3)
        1-element Vector{MX}:
         MX(v)
        julia> ca.symvar(ca.SX(3))
        SX[]
        \enddoctest
        \endlang

        \sa depends_on, contains

        \identifier{1bn} */
    inline friend std::vector<MatType> symvar(const MatType& x) {
      return MatType::symvar(x);
    }

    ///@{
    /** \brief Calculate bilinear/quadratic form x^T A y

        \p x and \p y must be vectors (row vectors are transposed, sparse vectors densified)
        with x.numel() == A.size1() and y.numel() == A.size2(); otherwise an error is raised.
        The result is scalar. When \p y is omitted, the quadratic form x^T A x is calculated.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> bilin(A, DM([1, 1]), DM([1, 2]))
        DM(16)
        >>> bilin(A, DM([1, 1]))
        DM(10)
        >>> bilin(MX.sym("A", 2, 2), MX.sym("x", 2), MX.sym("y", 2))
        MX(bilin(A, x, y))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> bilin(A, DM([1; 1]), DM([1; 2]))
        16
        >> bilin(A, DM([1; 1]))
        10
        >> bilin(MX.sym('A', 2, 2), MX.sym('x', 2), MX.sym('y', 2))
        bilin(A, x, y)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> ca.bilin(A, ca.DM([1, 1]), ca.DM([1, 2]))
        DM(16)
        julia> ca.bilin(A, ca.DM([1, 1]))
        DM(10)
        julia> ca.bilin(ca.MX.sym("A", 2, 2), ca.MX.sym("x", 2), ca.MX.sym("y", 2))
        MX(bilin(A, x, y))
        \enddoctest
        \endlang

        \sa rank1, dot, mtimes

        \identifier{1bo} */
    inline friend MatType bilin(const MatType &A, const MatType &x, const MatType &y) {
      return MatType::bilin(A, x, y);
    }
    inline friend MatType bilin(const MatType &A, const MatType &x) {
      return MatType::bilin(A, x, x);
    }
    static MatType bilin(const MatType& A, const MatType& x, const MatType& y);
    ///@}

    ///@{
    /** \brief Make a rank-1 update to a matrix A

        Calculates A + alpha * x*y'

        \param alpha Scalar; if it is a structural zero, \p A is returned unchanged
        \param x Vector with x.numel() == A.size1() (row vectors are transposed)
        \param y Vector with y.numel() == A.size2()

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> rank1(A, 2, DM([1, 0]), DM([0, 1]))
        DM(
        [[1, 4],
         [3, 4]])
        >>> rank1(MX.sym("A", 2, 2), MX.sym("a"), MX.sym("x", 2), MX.sym("y", 2))
        MX(rank1(A, a, x, y))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> rank1(A, 2, DM([1; 0]), DM([0; 1]))
        [[1, 4],
         [3, 4]]
        >> rank1(MX.sym('A', 2, 2), MX.sym('a'), MX.sym('x', 2), MX.sym('y', 2))
        rank1(A, a, x, y)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> ca.rank1(A, 2, ca.DM([1, 0]), ca.DM([0, 1]))
        DM(
        [[1, 4],
         [3, 4]])
        julia> ca.rank1(ca.MX.sym("A", 2, 2), ca.MX.sym("a"), ca.MX.sym("x", 2), ca.MX.sym("y", 2))
        MX(rank1(A, a, x, y))
        \enddoctest
        \endlang

        \sa bilin, mtimes

        \identifier{1bp} */
    inline friend MatType rank1(const MatType &A, const MatType &alpha,
                                const MatType &x, const MatType &y) {
      return MatType::rank1(A, alpha, x, y);
    }
    static MatType rank1(const MatType& A, const MatType& alpha,
                         const MatType& x, const MatType& y);
    ///@}

    /** \brief Calculate sum of squares: sum_ij X_ij^2

        Defined as dot(x, x) for a matrix of any shape; structural zeros contribute nothing.
        Prefer this over norm_2(x)**2 in objectives: it is smooth everywhere.

        \lang{python}
        \doctest
        >>> sumsqr(DM([[1, 2], [3, 4]]))
        DM(30)
        >>> sumsqr(SX.sym("x", 2))
        SX((sq(x_0)+sq(x_1)))
        >>> sumsqr(MX.sym("x", 2))
        MX(dot(x, x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sumsqr(DM([1 2; 3 4]))
        30
        >> sumsqr(SX.sym('x', 2))
        (sq(x_0)+sq(x_1))
        >> sumsqr(MX.sym('x', 2))
        dot(x, x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.sumsqr(ca.DM([1 2; 3 4]))
        DM(30)
        julia> ca.sumsqr(ca.SX.sym("x", 2))
        SX((sq(x_0)+sq(x_1)))
        julia> ca.sumsqr(ca.MX.sym("x", 2))
        MX(dot(x, x))
        \enddoctest
        \endlang

        \sa dot, norm_2, norm_fro

        \identifier{1bq} */
    inline friend MatType sumsqr(const MatType &x) {
      return MatType::sumsqr(x);
    }

    /** \brief x -> log(sum_i exp(x_i))

        Can be used to achieve a smooth max.

        The argument must be a dense column vector; a scalar is returned unchanged.
        Implementations for DM/SX/MX are hardened against overflow by shifting with the
        largest entry. Implementations for DM/MX are more accurate with log1p (*).

        (*) https://nhigham.com/2021/01/05/what-is-the-log-sum-exp-function/

        \lang{python}
        \doctest
        >>> logsumexp(DM([1, 2, 3]))
        DM(3.40761)
        >>> logsumexp(DM([1000, 1000]))
        DM(1000.69)
        >>> logsumexp(MX.sym("x", 3))
        MX(logsumexp(x))
        >>> try: logsumexp(DM([1, 2, 3]).T)
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> logsumexp(DM([1; 2; 3]))
        3.40761
        >> logsumexp(DM([1000; 1000]))
        1000.69
        >> logsumexp(MX.sym('x', 3))
        logsumexp(x)
        >> try, logsumexp(DM([1; 2; 3])'); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.logsumexp(ca.DM([1, 2, 3]))
        DM(3.40761)
        julia> ca.logsumexp(ca.DM([1000, 1000]))
        DM(1000.69)
        julia> ca.logsumexp(ca.MX.sym("x", 3))
        MX(logsumexp(x))
        julia> try ca.logsumexp(ca.DM([1, 2, 3])') catch e; println("error") end
        error
        \enddoctest
        \endlang

        \sa mmax

        \identifier{1br} */
    inline friend MatType logsumexp(const MatType& x) {
      return MatType::logsumexp(x);
    }
    /** \brief Scaled version of logsumexp

        Scaled such that max(x) <= logsumexp(x, margin) <= max(x)+margin,
        computed as logsumexp(alpha*x)/alpha with alpha = log(n)/margin for n = x.size1().
        \p margin must be of the same type as \p x (MX for MX).

        \lang{python}
        \doctest
        >>> logsumexp(DM([1, 2, 3]), 0.1)
        DM(3)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> logsumexp(DM([1; 2; 3]), 0.1)
        3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.logsumexp(ca.DM([1, 2, 3]), 0.1)
        DM(3)
        \enddoctest
        \endlang

        \identifier{1bs} */
    inline friend MatType logsumexp(const MatType& x, const MatType& margin) {
      MatType alpha = log(x.size1()) / margin;
      return MatType::logsumexp(alpha*x)/alpha;
    }
    static MatType logsumexp(const MatType& x);

    /** \brief Matlab's \c linspace command

        Returns the column vector [a, a + (b-a)/(nsteps-1), ..., b] with \p nsteps entries,
        the endpoints included. For nsteps == 1 the result is \p b. When \p a and \p b are
        matrices, the \p nsteps intermediate matrices are stacked vertically (vertcat).
        \p nsteps must be at least 1.

        \lang{python}
        \doctest
        >>> linspace(0, 1, 5)
        DM([0, 0.25, 0.5, 0.75, 1])
        >>> linspace(SX.sym("a"), SX.sym("b"), 3)
        SX([a, (a+((b-a)/2)), b])
        >>> linspace(DM([0, 0]), DM([1, 2]), 3)
        DM([0, 0, 0.5, 1, 1, 2])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> linspace(DM(0), DM(1), 5)
        [0, 0.25, 0.5, 0.75, 1]
        >> linspace(SX.sym('a'), SX.sym('b'), 3)
        [a, (a+((b-a)/2)), b]
        >> linspace(DM([0; 0]), DM([1; 2]), 3)
        [0, 0, 0.5, 1, 1, 2]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.linspace(0, 1, 5)
        DM([0, 0.25, 0.5, 0.75, 1])
        julia> ca.linspace(ca.SX.sym("a"), ca.SX.sym("b"), 3)
        SX([a, (a+((b-a)/2)), b])
        julia> ca.linspace(ca.DM([0, 0]), ca.DM([1, 2]), 3)
        DM([0, 0, 0.5, 1, 1, 2])
        \enddoctest
        \endlang

        \identifier{1bt} */
    inline friend MatType linspace(const MatType &a, const MatType &b, casadi_int nsteps) {
      return MatType::linspace(a, b, nsteps);
    }

    /** \brief Matlab's \c cross command

        Cross product of 3-vectors. \p a and \p b must have the same shape and one of the
        dimensions must be 3; for 3-by-n or n-by-3 matrices the products are taken per column
        or per row respectively.

        \param dim 1 to take the products along the rows (3-by-n input), 2 along the columns,
        -1 (default) to pick the dimension of length 3 (rows first)

        \lang{python}
        \doctest
        >>> cross(DM([1, 0, 0]), DM([0, 1, 0]))
        DM([0, 0, 1])
        >>> cross(DM([[1, 0, 0], [0, 1, 0]]), DM([[0, 1, 0], [0, 0, 1]]))
        DM(
        [[0, 0, 1],
         [1, 0, 0]])
        >>> cross(SX.sym("a", 3), SX.sym("b", 3))
        SX([((a_1*b_2)-(a_2*b_1)), ((a_2*b_0)-(a_0*b_2)), ((a_0*b_1)-(a_1*b_0))])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> cross(DM([1; 0; 0]), DM([0; 1; 0]))
        [0, 0, 1]
        >> cross(DM([1 0 0; 0 1 0]), DM([0 1 0; 0 0 1]))
        [[0, 0, 1],
         [1, 0, 0]]
        >> cross(SX.sym('a', 3), SX.sym('b', 3))
        [((a_1*b_2)-(a_2*b_1)), ((a_2*b_0)-(a_0*b_2)), ((a_0*b_1)-(a_1*b_0))]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.cross(ca.DM([1, 0, 0]), ca.DM([0, 1, 0]))
        DM([0, 0, 1])
        julia> ca.cross(ca.DM([1 0 0; 0 1 0]), ca.DM([0 1 0; 0 0 1]))
        DM(
        [[0, 0, 1],
         [1, 0, 0]])
        julia> ca.cross(ca.SX.sym("a", 3), ca.SX.sym("b", 3))
        SX([((a_1*b_2)-(a_2*b_1)), ((a_2*b_0)-(a_0*b_2)), ((a_0*b_1)-(a_1*b_0))])
        \enddoctest
        \endlang

        \sa skew, dot

        \identifier{1bu} */
    inline friend MatType cross(const MatType &a, const MatType &b, casadi_int dim = -1) {
      return MatType::cross(a, b, dim);
    }

    /** \brief Generate a skew symmetric matrix from a 3-vector

        skew(a) is the 3-by-3 matrix S with S*b == cross(a, b).

        \lang{python}
        \doctest
        >>> skew(DM([1, 2, 3]))
        DM(
        [[0, -3, 2],
         [3, 0, -1],
         [-2, 1, 0]])
        >>> inv_skew(skew(DM([1, 2, 3])))
        DM([1, 2, 3])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> skew(DM([1; 2; 3]))
        [[0, -3, 2],
         [3, 0, -1],
         [-2, 1, 0]]
        >> inv_skew(skew(DM([1; 2; 3])))
        [1, 2, 3]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.skew(ca.DM([1, 2, 3]))
        DM(
        [[0, -3, 2],
         [3, 0, -1],
         [-2, 1, 0]])
        julia> ca.inv_skew(ca.skew(ca.DM([1, 2, 3])))
        DM([1, 2, 3])
        \enddoctest
        \endlang

        \sa inv_skew, cross

        \identifier{1bv} */
    inline friend MatType skew(const MatType &a) {
      return MatType::skew(a);
    }

    /** \brief Generate the 3-vector progenitor of a skew symmetric matrix

        Inverse of skew: for a 3-by-3 matrix S, returns 0.5*[S(2,1)-S(1,2); S(0,2)-S(2,0);
        S(1,0)-S(0,1)], which recovers \p a from skew(a). A non-3-by-3 argument raises an
        error.

        \sa skew

        \identifier{1bw} */
    inline friend MatType inv_skew(const MatType &a) {
      return MatType::inv_skew(a);
    }

    /** \brief Matrix determinant (experimental)

        For DM and SX the determinant is formed by minor expansion along the sparsest
        row or column (a 2-by-2 directly, a scalar is returned as is); the cost grows
        quickly with the matrix size. For MX a determinant node is created that is
        evaluated through the 'qr' linear solver. A non-square argument raises an error.

        \lang{python}
        \doctest
        >>> det(DM([[1, 2], [3, 4]]))
        DM(-2)
        >>> det(SX.sym("A", 2, 2))
        SX(((A_0*A_3)-(A_2*A_1)))
        >>> det(MX.sym("A", 2, 2))
        MX(det(A))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> det(DM([1 2; 3 4]))
        -2
        >> det(SX.sym('A', 2, 2))
        ((A_0*A_3)-(A_2*A_1))
        >> det(MX.sym('A', 2, 2))
        det(A)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.det(ca.DM([1 2; 3 4]))
        DM(-2)
        julia> ca.det(ca.SX.sym("A", 2, 2))
        SX(((A_0*A_3)-(A_2*A_1)))
        julia> ca.det(ca.MX.sym("A", 2, 2))
        MX(det(A))
        \enddoctest
        \endlang

        \sa inv, inv_minor

        \identifier{1bx} */
    inline friend MatType det(const MatType& A) { return MatType::det(A);}

    /** \brief Matrix determinant through a linear solver factorization (experimental)

        \param lsolver Name of a Linsol plugin, e.g. "qr" or "symbolicqr" (for SX only plugins
        with a symbolic determinant, such as "symbolicqr", are accepted)
        \param dict Options passed to the linear solver

        \lang{python}
        \doctest
        >>> det(DM([[1, 2], [3, 4]]), "qr")
        DM(-2)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> det(DM([1 2; 3 4]), 'qr')
        -2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.det(ca.DM([1 2; 3 4]), "qr")
        DM(-2)
        \enddoctest
        \endlang

        \sa Linsol

        \identifier{2hx} */
    inline friend MatType det(const MatType& A, const std::string& lsolver,
                              const Dict& dict=Dict()) {
      return MatType::det(A, lsolver, dict);
    }

    /** \brief Matrix inverse by the adjugate (experimental)

        Computes adj(A)/det(A) by minor expansion. Available for DM and SX; MX raises
        "Not implemented". Use inv for the general case.

        \lang{python}
        \doctest
        >>> inv_minor(DM([[1, 2], [3, 4]]))
        DM(
        [[-2, 1],
         [1.5, -0.5]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> inv_minor(DM([1 2; 3 4]))
        [[-2, 1],
         [1.5, -0.5]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.inv_minor(ca.DM([1 2; 3 4]))
        DM(
        [[-2, 1],
         [1.5, -0.5]])
        \enddoctest
        \endlang

        \sa inv, det

        \identifier{1by} */
    inline friend MatType inv_minor(const MatType& A) { return MatType::inv_minor(A);}

    /** \brief Matrix inverse

        Computed as solve(A, eye(n)), see solve for the algorithm. A non-square argument
        raises an error. A numerically singular DM does not raise; the result contains
        inf or nan. For MX the result is a solve node.

        \lang{python}
        \doctest
        >>> inv(DM([[1, 2], [3, 4]]))
        DM(
        [[-2, 1],
         [1.5, -0.5]])
        >>> inv(SX.sym("A", 2, 2))
        SX(@1=((A_0*A_3)-(A_2*A_1)),
        [[(A_3/@1), (-(A_2/@1))],
         [(-(A_1/@1)), (A_0/@1)]])
        >>> inv(DM(Sparsity.diag(2), [2, 4]))
        DM(
        [[0.5, 00],
         [00, 0.25]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> inv(DM([1 2; 3 4]))
        [[-2, 1],
         [1.5, -0.5]]
        >> inv(SX.sym('A', 2, 2))
        @1=((A_0*A_3)-(A_2*A_1)),
        [[(A_3/@1), (-(A_2/@1))],
         [(-(A_1/@1)), (A_0/@1)]]
        >> inv(DM(Sparsity.diag(2), [2; 4]))
        [[0.5, 00],
         [00, 0.25]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> print(inv(ca.DM([1 2; 3 4])))
        [[-2, 1],
         [1.5, -0.5]]
        julia> inv(ca.SX.sym("A", 2, 2))
        SX(@1=((A_0*A_3)-(A_2*A_1)),
        [[(A_3/@1), (-(A_2/@1))],
         [(-(A_1/@1)), (A_0/@1)]])
        julia> inv(ca.DM(ca.Sparsity.diag(2), [2, 4]))
        DM(
        [[0.5, 00],
         [00, 0.25]])
        \enddoctest
        \endlang

        \sa solve, pinv, det, mldivide

        \identifier{1bz} */
    inline friend MatType inv(const MatType& A) {
        return MatType::inv(A);
    }

    /** \brief Matrix inverse using a linear solver plugin

        Computes solve(A, eye(n), lsolver, options). Available for DM and MX; SX raises.

        \param lsolver Name of a Linsol plugin, e.g. "qr", "csparse", "lapacklu"
        \param options Options passed to the linear solver

        \lang{python}
        \doctest
        >>> inv(DM([[1, 2], [3, 4]]), "qr")
        DM(
        [[-2, 1],
         [1.5, -0.5]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> inv(DM([1 2; 3 4]), 'qr')
        [[-2, 1],
         [1.5, -0.5]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> print(ca.inv(ca.DM([1 2; 3 4]), "qr"))
        [[-2, 1],
         [1.5, -0.5]]
        \enddoctest
        \endlang

        \sa solve, Linsol

        \identifier{1c0} */
    inline friend MatType inv(const MatType& A,
      const std::string& lsolver,
      const Dict& options=Dict()) {
        return MatType::inv(A, lsolver, options);
    }

    /** \brief Matrix trace

        Sum of the diagonal entries of a square matrix (structural zeros contribute nothing).
        A non-square argument raises an error.

        \lang{python}
        \doctest
        >>> trace(DM([[1, 2], [3, 4]]))
        DM(5)
        >>> trace(SX.sym("A", 2, 2))
        SX((A_0+A_3))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> trace(DM([1 2; 3 4]))
        5
        >> trace(SX.sym('A', 2, 2))
        (A_0+A_3)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.tr(ca.DM([1 2; 3 4]))
        DM(5)
        julia> ca.tr(ca.SX.sym("A", 2, 2))
        SX((A_0+A_3))
        \enddoctest
        \endlang

        \sa diag

        \identifier{1c1} */
    inline friend MatType trace(const MatType& x) { return MatType::trace(x);}

    /** \brief Convert a lower triangular matrix to a symmetric one

        Returns a + a' - diag(diag(a)), i.e. mirrors the strictly lower part to the upper
        part. The argument must be square and have no structural nonzeros above the diagonal
        (a dense matrix is rejected even if numerically lower triangular).

        \lang{python}
        \doctest
        >>> L = DM(Sparsity.lower(3), [1, 2, 3, 4, 5, 6])
        >>> tril2symm(L)
        DM(
        [[1, 2, 3],
         [2, 4, 5],
         [3, 5, 6]])
        >>> tril2symm(SX.sym("L", Sparsity.lower(2)))
        SX(
        [[L_0, L_1],
         [L_1, L_2]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> L = DM(Sparsity.lower(3), [1; 2; 3; 4; 5; 6]);
        >> tril2symm(L)
        [[1, 2, 3],
         [2, 4, 5],
         [3, 5, 6]]
        >> tril2symm(SX.sym('L', Sparsity.lower(2)))
        [[L_0, L_1],
         [L_1, L_2]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> L = ca.DM(ca.Sparsity.lower(3), [1, 2, 3, 4, 5, 6]);
        julia> ca.tril2symm(L)
        DM(
        [[1, 2, 3],
         [2, 4, 5],
         [3, 5, 6]])
        julia> ca.tril2symm(ca.SX.sym("L", ca.Sparsity.lower(2)))
        SX(
        [[L_0, L_1],
         [L_1, L_2]])
        \enddoctest
        \endlang

        \sa triu2symm, tril, is_tril

        \identifier{1c2} */
    inline friend MatType tril2symm(const MatType &a) { return MatType::tril2symm(a);}

    /** \brief Convert a upper triangular matrix to a symmetric one

        Mirror image of tril2symm: the argument must be square without structural nonzeros
        below the diagonal.

        \lang{python}
        \doctest
        >>> U = DM(Sparsity.upper(2), [1, 2, 3])
        >>> triu2symm(U)
        DM(
        [[1, 2],
         [2, 3]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> U = DM(Sparsity.upper(2), [1; 2; 3]);
        >> triu2symm(U)
        [[1, 2],
         [2, 3]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> U = ca.DM(ca.Sparsity.upper(2), [1, 2, 3]);
        julia> ca.triu2symm(U)
        DM(
        [[1, 2],
         [2, 3]])
        \enddoctest
        \endlang

        \sa tril2symm, triu, is_triu

        \identifier{1c3} */
    inline friend MatType triu2symm(const MatType &a) { return MatType::triu2symm(a);}

    /** \brief  Frobenius norm

        sqrt(sum_ij x_ij^2) over all entries, for a matrix of any shape.
        Not differentiable at x == 0 (the derivative evaluates to nan there); use sumsqr
        in objectives where possible.

        \lang{python}
        \doctest
        >>> norm_fro(DM([[1, -2], [3, 4]]))
        DM(5.47723)
        >>> norm_fro(SX.sym("x", 2))
        SX(sqrt((sq(x_0)+sq(x_1))))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> norm_fro(DM([1 -2; 3 4]))
        5.47723
        >> norm_fro(SX.sym('x', 2))
        sqrt((sq(x_0)+sq(x_1)))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.norm_fro(ca.DM([1 -2; 3 4]))
        DM(5.47723)
        julia> ca.norm_fro(ca.SX.sym("x", 2))
        SX(sqrt((sq(x_0)+sq(x_1))))
        \enddoctest
        \endlang

        \sa norm_2, norm_1, norm_inf, sumsqr

        \identifier{1c4} */
    inline friend MatType norm_fro(const MatType &x) { return MatType::norm_fro(x);}

    /** \brief  2-norm

        For a vector (row or column) this is the Euclidean norm, identical to norm_fro.
        For a matrix argument DM and SX raise an error (use norm_fro if that is what you
        meant); MX creates a spectral-norm node that cannot be evaluated.
        Not differentiable at x == 0; use sumsqr in objectives where possible.

        \lang{python}
        \doctest
        >>> norm_2(DM([3, 4]))
        DM(5)
        >>> norm_2(MX.sym("x", 2))
        MX(||x||_F)
        >>> try: norm_2(DM([[1, 2], [3, 4]]))
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> norm_2(DM([3; 4]))
        5
        >> norm_2(MX.sym('x', 2))
        ||x||_F
        >> try, norm_2(DM([1 2; 3 4])); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.norm_2(ca.DM([3, 4]))
        DM(5)
        julia> ca.norm_2(ca.MX.sym("x", 2))
        MX(||x||_F)
        julia> try ca.norm_2(ca.DM([1 2; 3 4])) catch e; println("error") end
        error
        \enddoctest
        \endlang

        \sa norm_fro, norm_1, norm_inf, sumsqr

        \identifier{1c5} */
    inline friend MatType norm_2(const MatType &x) { return MatType::norm_2(x);}

    /** \brief 1-norm

        Sum of the absolute values of all entries. For a matrix this is the entrywise norm,
        not the induced (maximum column sum) matrix norm. Nonsmooth at zero entries.

        \lang{python}
        \doctest
        >>> norm_1(DM([[1, -2], [3, 4]]))
        DM(10)
        >>> norm_1(SX.sym("x", 2))
        SX((fabs(x_0)+fabs(x_1)))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> norm_1(DM([1 -2; 3 4]))
        10
        >> norm_1(SX.sym('x', 2))
        (fabs(x_0)+fabs(x_1))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.norm_1(ca.DM([1 -2; 3 4]))
        DM(10)
        julia> ca.norm_1(ca.SX.sym("x", 2))
        SX((fabs(x_0)+fabs(x_1)))
        \enddoctest
        \endlang

        \sa norm_inf, norm_2, norm_fro

        \identifier{1c6} */
    inline friend MatType norm_1(const MatType &x) { return MatType::norm_1(x);}

    /** \brief Infinity-norm

        Largest absolute value among all entries. For a matrix this is the entrywise norm,
        not the induced (maximum row sum) matrix norm. Nonsmooth.

        \lang{python}
        \doctest
        >>> norm_inf(DM([[1, -2], [3, 4]]))
        DM(4)
        >>> norm_inf(SX.sym("x", 2))
        SX(fmax(fmax(0,fabs(x_0)),fabs(x_1)))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> norm_inf(DM([1 -2; 3 4]))
        4
        >> norm_inf(SX.sym('x', 2))
        fmax(fmax(0,fabs(x_0)),fabs(x_1))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.norm_inf(ca.DM([1 -2; 3 4]))
        DM(4)
        julia> ca.norm_inf(ca.SX.sym("x", 2))
        SX(fmax(fmax(0,fabs(x_0)),fabs(x_1)))
        \enddoctest
        \endlang

        \sa norm_1, norm_2, norm_fro, mmax

        \identifier{1c7} */
    inline friend MatType norm_inf(const MatType &x) { return MatType::norm_inf(x);}

    /** \brief Returns difference (n-th order) along given axis (MATLAB convention)

        Each pass subtracts neighbouring rows (axis 0) or columns (axis 1), shrinking that
        dimension by one; \p n passes are applied. With axis -1 the first non-singleton
        dimension is used (columns for a row vector, rows otherwise) and, as in MATLAB,
        a scalar gives an empty 0-by-0 result.

        \param n Order, at least 1
        \param axis 0 (rows), 1 (columns) or -1 (automatic); 1 and 2 in MATLAB

        \lang{python}
        \doctest
        >>> diff(DM([1, 4, 9, 16]))
        DM([3, 5, 7])
        >>> diff(DM([1, 4, 9, 16]), 2)
        DM([2, 2])
        >>> diff(DM([[1, 2], [4, 8], [9, 10]]), 1, 1)
        DM([1, 4, 1])
        >>> diff(SX.sym("x", 3))
        SX([(x_1-x_0), (x_2-x_1)])
        >>> diff(DM(3)).dim()
        '0x0'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> diff(DM([1; 4; 9; 16]))
        [3, 5, 7]
        >> diff(DM([1; 4; 9; 16]), 2)
        [2, 2]
        >> diff(DM([1 2; 4 8; 9 10]), 1, 2)
        [1, 4, 1]
        >> diff(SX.sym('x', 3))
        [(x_1-x_0), (x_2-x_1)]
        >> d = diff(DM(3)); disp(d.dim())
        0x0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> diff(ca.DM([1, 4, 9, 16]))
        DM([3, 5, 7])
        julia> ca.diff(ca.DM([1, 4, 9, 16]), 2)
        DM([2, 2])
        julia> diff(ca.DM([1 2; 4 8; 9 10]); dims=2)
        DM([1, 4, 1])
        julia> diff(ca.SX.sym("x", 3))
        SX([(x_1-x_0), (x_2-x_1)])
        julia> ca.dim(diff(ca.DM(3)))
        "0x1"
        \enddoctest
        \endlang

        \sa cumsum

        \identifier{1c8} */
    inline friend MatType diff(const MatType &x, casadi_int n=1, casadi_int axis=-1) {
      return MatType::diff(x, n, axis);
    }

    /** \brief Returns cumulative sum along given axis (MATLAB convention)

        Running sum down the rows (axis 0) or along the columns (axis 1); the shape is
        preserved. With axis -1, the columns are summed for a row vector and the rows
        otherwise. MATLAB numbers the axes 1 (rows) and 2 (columns).

        \lang{python}
        \doctest
        >>> cumsum(DM([1, 2, 3]))
        DM([1, 3, 6])
        >>> cumsum(DM([[1, 2], [3, 4]]))
        DM(
        [[1, 2],
         [4, 6]])
        >>> cumsum(DM([[1, 2], [3, 4]]), 1)
        DM(
        [[1, 3],
         [3, 7]])
        >>> cumsum(SX.sym("x", 3))
        SX(@1=(x_1+x_0), [x_0, @1, (x_2+@1)])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> cumsum(DM([1; 2; 3]))
        [1, 3, 6]
        >> cumsum(DM([1 2; 3 4]))
        [[1, 2],
         [4, 6]]
        >> cumsum(DM([1 2; 3 4]), 2)
        [[1, 3],
         [3, 7]]
        >> cumsum(SX.sym('x', 3))
        @1=(x_1+x_0), [x_0, @1, (x_2+@1)]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> cumsum(ca.DM([1, 2, 3]))
        DM([1, 3, 6])
        julia> cumsum(ca.DM([1 2; 3 4]); dims=1)
        DM(
        [[1, 2],
         [4, 6]])
        julia> cumsum(ca.DM([1 2; 3 4]); dims=2)
        DM(
        [[1, 3],
         [3, 7]])
        julia> cumsum(ca.SX.sym("x", 3))
        SX(@1=(x_1+x_0), [x_0, @1, (x_2+@1)])
        \enddoctest
        \endlang

        \sa diff, sum1, sum2

        \identifier{1c9} */
    inline friend MatType cumsum(const MatType &x, casadi_int axis=-1) {
      return MatType::cumsum(x, axis);
    }

    /** \brief Inner product of two matrices

        Sum of the elementwise products over all entries, with \p x and \p y of the same
        shape (a matrix is treated as the vector of its entries). Only entries that are
        structurally nonzero in both arguments contribute. A shape mismatch raises an error.

        \lang{python}
        \doctest
        >>> dot(DM([1, 2]), DM([3, 4]))
        DM(11)
        >>> dot(DM([[1, 2], [3, 4]]), DM([[1, 1], [1, 1]]))
        DM(10)
        >>> dot(DM(Sparsity.diag(2), [1, 2]), DM([[1, 1], [1, 1]]))
        DM(3)
        >>> dot(SX.sym("x", 2), SX.sym("y", 2))
        SX(((x_0*y_0)+(x_1*y_1)))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> dot(DM([1; 2]), DM([3; 4]))
        11
        >> dot(DM([1 2; 3 4]), DM([1 1; 1 1]))
        10
        >> dot(DM(Sparsity.diag(2), [1; 2]), DM([1 1; 1 1]))
        3
        >> dot(SX.sym('x', 2), SX.sym('y', 2))
        ((x_0*y_0)+(x_1*y_1))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.dot(ca.DM([1, 2]), ca.DM([3, 4]))
        DM(11)
        julia> ca.dot(ca.DM([1 2; 3 4]), ca.DM([1 1; 1 1]))
        DM(10)
        julia> ca.dot(ca.DM(ca.Sparsity.diag(2), [1, 2]), ca.DM([1 1; 1 1]))
        DM(3)
        julia> ca.dot(ca.SX.sym("x", 2), ca.SX.sym("y", 2))
        SX(((x_0*y_0)+(x_1*y_1)))
        \enddoctest
        \endlang

        \sa mtimes, sumsqr, bilin

        \identifier{1ca} */
    inline friend MatType dot(const MatType &x, const MatType &y) {
      return MatType::dot(x, y);
    }

    /** \brief Computes the nullspace of a matrix A

        Finds Z m-by-(m-n) such that AZ = 0
        with A n-by-m with m >= n (a matrix with more rows than columns raises an error).

        Assumes A is full rank. Computed with Householder reflections, dense; for MX the
        SX algorithm is embedded as a function call node.

        Inspired by Numerical Methods in Scientific Computing by Ake Bjorck

        \lang{python}
        \doctest
        >>> A = DM([[1, 2, 3], [4, 5, 6]])
        >>> Z = nullspace(A)
        >>> Z.dim(), norm_1(mtimes(A, Z))
        ('3x1', DM(0))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2 3; 4 5 6]);
        >> Z = nullspace(A);
        >> Z.dim(), norm_1(A*Z)
        3x1
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2 3; 4 5 6]);
        julia> Z = ca.nullspace(A);
        julia> ca.dim(Z), ca.norm_1(A * Z)
        ("3x1", 0)
        \enddoctest
        \endlang

        \sa qr, pinv

        \identifier{1cb} */
    inline friend MatType nullspace(const MatType& A) {
      return MatType::nullspace(A);
    }

    /** \brief  Evaluate a polynomial with coefficients p in x

        Horner evaluation of p[0]*x^(n-1) + ... + p[n-2]*x + p[n-1]: the first coefficient
        belongs to the highest power, as in MATLAB's polyval. The polynomial is applied to
        every structurally nonzero entry of \p x; the result keeps the shape and sparsity of
        \p x (structural zeros stay zero even when p[n-1] != 0).

        \param p Dense column vector of coefficients, at least one entry
        \param x Point(s) of evaluation, any shape

        \lang{python}
        \doctest
        >>> polyval(DM([1, 2, 3]), 2)
        DM(11)
        >>> polyval(DM([1, 0, -1]), DM([0, 1, 2]))
        DM([-1, 0, 3])
        >>> polyval(SX.sym("p", 3), SX.sym("x"))
        SX(((((p_0*x)+p_1)*x)+p_2))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> polyval(DM([1; 2; 3]), 2)
        11
        >> polyval(DM([1; 0; -1]), DM([0; 1; 2]))
        [-1, 0, 3]
        >> polyval(SX.sym('p', 3), SX.sym('x'))
        ((((p_0*x)+p_1)*x)+p_2)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.polyval(ca.DM([1, 2, 3]), 2)
        DM(11)
        julia> ca.polyval(ca.DM([1, 0, -1]), ca.DM([0, 1, 2]))
        DM([-1, 0, 3])
        julia> ca.polyval(ca.SX.sym("p", 3), ca.SX.sym("x"))
        SX(((((p_0*x)+p_1)*x)+p_2))
        \enddoctest
        \endlang

        \identifier{1cc} */
    inline friend MatType polyval(const MatType& p, const MatType& x) {
      return MatType::polyval(p, x);
    }

    /** \brief   Get the diagonal of a matrix or construct a diagonal

        If the input is a vector (row or column), a square matrix with the vector on its
        diagonal and structural zeros elsewhere is constructed. Otherwise the diagonal
        entries are returned as a column vector of length min(size1(), size2()); only
        structurally nonzero entries are kept.

        \lang{python}
        \doctest
        >>> diag(DM([1, 2]))
        DM(
        [[1, 00],
         [00, 2]])
        >>> diag(DM([[1, 2], [3, 4]]))
        DM([1, 4])
        >>> diag(DM(Sparsity.upper(2), [1, 2, 3]))
        DM([1, 3])
        >>> diag(DM.ones(2, 3))
        DM([1, 1])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> diag(DM([1; 2]))
        [[1, 00],
         [00, 2]]
        >> diag(DM([1 2; 3 4]))
        [1, 4]
        >> diag(DM(Sparsity.upper(2), [1; 2; 3]))
        [1, 3]
        >> diag(DM.ones(2, 3))
        [1, 1]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.diag(ca.DM([1, 2]))
        DM(
        [[1, 00],
         [00, 2]])
        julia> ca.diag(ca.DM([1 2; 3 4]))
        DM([1, 4])
        julia> ca.diag(ca.DM(ca.Sparsity.upper(2), [1, 2, 3]))
        DM([1, 3])
        julia> ca.diag(ca.DM.ones(2, 3))
        DM([1, 1])
        \enddoctest
        \endlang

        \sa trace, diagcat, Sparsity::diag

        \identifier{1cd} */
    inline friend MatType diag(const MatType &A) {
      return MatType::diag(A);
    }

    /** \brief  Unite two matrices with non-overlapping sparsity patterns

        The result has the union of both patterns with the nonzeros of \p A and \p B at their
        positions. Both must have the same shape; a shared structural nonzero raises
        "Pattern intersection not empty".

        \lang{python}
        \doctest
        >>> A = DM(Sparsity.upper(2), [1, 2, 3])
        >>> B = DM(Sparsity.triplet(2, 2, [1], [0]), [4])
        >>> unite(A, B)
        DM(
        [[1, 2],
         [4, 3]])
        >>> try: unite(DM.eye(2), DM.eye(2))
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM(Sparsity.upper(2), [1; 2; 3]);
        >> B = DM(Sparsity.triplet(2, 2, 1, 0), 4);
        >> unite(A, B)
        [[1, 2],
         [4, 3]]
        >> try, unite(DM.eye(2), DM.eye(2)); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM(ca.Sparsity.upper(2), [1, 2, 3]);
        julia> B = ca.DM(ca.Sparsity.triplet(2, 2, [1], [0]), [4]);
        julia> ca.unite(A, B)
        DM(
        [[1, 2],
         [4, 3]])
        julia> try ca.unite(ca.DM.eye(2), ca.DM.eye(2)) catch e; println("error") end
        error
        \enddoctest
        \endlang

        \sa project, Sparsity::unite

        \identifier{1ce} */
    inline friend MatType unite(const MatType& A, const MatType& B) {
      return MatType::unite(A, B);
    }

    /** \brief  Make the matrix dense if not already

        Structural zeros become stored numerical zeros; the shape and the values of the
        existing nonzeros are unchanged. A dense input is returned as is.

        \lang{python}
        \doctest
        >>> densify(DM.eye(2))
        DM(
        [[1, 0],
         [0, 1]])
        >>> densify(SX.sym("x", Sparsity.diag(2)))
        SX(@1=0,
        [[x_0, @1],
         [@1, x_1]])
        >>> densify(MX.sym("x", Sparsity.diag(2)))
        MX(dense(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> densify(DM.eye(2))
        [[1, 0],
         [0, 1]]
        >> densify(SX.sym('x', Sparsity.diag(2)))
        @1=0,
        [[x_0, @1],
         [@1, x_1]]
        >> densify(MX.sym('x', Sparsity.diag(2)))
        dense(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.densify(ca.DM.eye(2))
        DM(
        [[1, 0],
         [0, 1]])
        julia> ca.densify(ca.SX.sym("x", ca.Sparsity.diag(2)))
        SX(@1=0,
        [[x_0, @1],
         [@1, x_1]])
        julia> ca.densify(ca.MX.sym("x", ca.Sparsity.diag(2)))
        MX(dense(x))
        \enddoctest
        \endlang

        \sa project, sparsify, is_dense

        \identifier{1cf} */
    inline friend MatType densify(const MatType& x) {
      return MatType::densify(x);
    }

    /** \brief  Make the matrix dense and assign nonzeros to a value

        Like densify(x) but the structural zeros get the scalar \p val instead of zero.

        \identifier{1cg} */
    inline friend MatType densify(const MatType& x, const MatType& val) {
      return MatType::densify(x, val);
    }

    /** \brief Create a new matrix with a given sparsity pattern but with the

      * nonzeros taken from an existing matrix

        Entries of \p A outside \p sp are dropped; entries of \p sp not present in \p A
        become numerical zeros. The shape of \p sp must equal that of \p A. With \p intersect
        true, the target pattern is the intersection of \p sp and the pattern of \p A, so the
        result only ever loses entries.

        \lang{python}
        \doctest
        >>> project(DM([[1, 2], [3, 4]]), Sparsity.diag(2))
        DM(
        [[1, 00],
         [00, 4]])
        >>> project(DM.eye(2), Sparsity.dense(2, 2))
        DM(
        [[1, 0],
         [0, 1]])
        >>> project(DM(Sparsity.diag(2), [1, 2]), Sparsity.upper(2), True)
        DM(
        [[1, 00],
         [00, 2]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> project(DM([1 2; 3 4]), Sparsity.diag(2))
        [[1, 00],
         [00, 4]]
        >> project(DM.eye(2), Sparsity.dense(2, 2))
        [[1, 0],
         [0, 1]]
        >> project(DM(Sparsity.diag(2), [1; 2]), Sparsity.upper(2), true)
        [[1, 00],
         [00, 2]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.project(ca.DM([1 2; 3 4]), ca.Sparsity.diag(2))
        DM(
        [[1, 00],
         [00, 4]])
        julia> ca.project(ca.DM.eye(2), ca.Sparsity.dense(2, 2))
        DM(
        [[1, 0],
         [0, 1]])
        julia> ca.project(ca.DM(ca.Sparsity.diag(2), [1, 2]), ca.Sparsity.upper(2), true)
        DM(
        [[1, 00],
         [00, 2]])
        \enddoctest
        \endlang

        \sa densify, sparsify, unite

        \identifier{1ch} */
    inline friend MatType project(const MatType& A, const Sparsity& sp,
                                  bool intersect=false) {
      return MatType::project(A, sp, intersect);
    }

    /** \brief Branching on MX nodes

        Ternary operator, "cond ? if_true : if_false"

        Selects \p if_true where \p cond is nonzero and \p if_false elsewhere. \p cond may be
        a scalar or have the shape of the branches (elementwise selection); a scalar branch is
        broadcast. By default the result is if_else_zero(cond, if_true) +
        if_else_zero(!cond, if_false): both branches are part of the expression and are
        evaluated, but the value of the unselected branch is discarded, so a nan or inf there
        does not propagate. With \p short_circuit true (MX only; ignored for SX and DM) the
        branches are wrapped in functions and only the selected one is evaluated.

        \lang{python}
        \doctest
        >>> if_else(DM(1), DM(10), DM(20)), if_else(DM(0), DM(10), DM(20))
        (DM(10), DM(20))
        >>> if_else(DM([1, 0]), DM([1, 2]), DM([3, 4]))
        DM([1, 4])
        >>> x = SX.sym("x")
        >>> if_else(x > 0, x, -x)
        SX(@1=(0<x), ((@1?x:0)+((!@1)?(-x):0)))
        >>> f = Function("f", [x], [if_else(x >= 0, sqrt(x), 0)])
        >>> f(-1)
        DM(0)
        >>> X = MX.sym("X")
        >>> if_else(X > 0, X, -X, True)
        MX(switch((0<X), X){0})
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> if_else(DM(1), DM(10), DM(20)), if_else(DM(0), DM(10), DM(20))
        10
        20
        >> if_else(DM([1; 0]), DM([1; 2]), DM([3; 4]))
        [1, 4]
        >> x = SX.sym('x');
        >> if_else(x > 0, x, -x)
        @1=(0<x), ((@1?x:0)+((!@1)?(-x):0))
        >> f = Function('f', {x}, {if_else(x >= 0, sqrt(x), 0)});
        >> f(-1)
        0
        >> X = MX.sym('X');
        >> if_else(X > 0, X, -X, true)
        switch((0<X), X){0}
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> a = ca.DM(10); b = ca.DM(20);
        julia> ca.if_else(ca.DM(1), a, b), ca.if_else(ca.DM(0), a, b)
        (10, 20)
        julia> ca.if_else(ca.DM([1, 0]), ca.DM([1, 2]), ca.DM([3, 4]))
        DM([1, 4])
        julia> x = ca.SX.sym("x");
        julia> ca.if_else(x > 0, x, -x)
        SX(@1=(0<x), ((@1?x:0)+((!@1)?(-x):0)))
        julia> f = ca.Function("f", [x], [ca.if_else(x >= 0, sqrt(x), 0)]);
        julia> f(-1)
        DM(0)
        julia> X = ca.MX.sym("X");
        julia> ca.if_else(X > 0, X, -X, true)
        MX(switch((0<X), X){0})
        \enddoctest
        \endlang

        \sa conditional, if_else_zero, Function::if_else

        \identifier{1ci} */
    inline friend MatType if_else(const MatType &cond, const MatType &if_true,
                                  const MatType &if_false, bool short_circuit=false) {
      return MatType::if_else(cond, if_true, if_false, short_circuit);
    }

    /** \brief Create a switch

        If the scalar condition \p ind evaluates to the integer k, where 0<=k<x.size(),
        then x[k] is returned, otherwise \p x_default is returned (also for negative or
        non-integer values). By default this is a chain of if_else expressions, so all cases
        are evaluated. With \p short_circuit true (MX only; SX raises an error) a Switch
        function evaluates only the selected case.

        \lang{python}
        \doctest
        >>> conditional(DM(1), [DM(10), DM(20)], DM(-1))
        DM(20)
        >>> conditional(DM(5), [DM(10), DM(20)], DM(-1))
        DM(-1)
        >>> x = SX.sym("x"); y = SX.sym("y"); z = SX.sym("z")
        >>> conditional(x, [y, z], 0)
        SX(@1=(x==1), ((@1?z:0)+((!@1)?((x==0)?y:0):0)))
        >>> X = MX.sym("X"); Y = MX.sym("Y")
        >>> conditional(X, [Y, 2*Y], 0, True)
        MX(switch(X, Y){0})
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> conditional(DM(1), {DM(10), DM(20)}, DM(-1))
        20
        >> conditional(DM(5), {DM(10), DM(20)}, DM(-1))
        -1
        >> x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z');
        >> conditional(x, {y, z}, 0)
        @1=(x==1), ((@1?z:0)+((!@1)?((x==0)?y:0):0))
        >> X = MX.sym('X'); Y = MX.sym('Y');
        >> conditional(X, {Y, 2*Y}, 0, true)
        switch(X, Y){0}
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.conditional(ca.DM(1), [ca.DM(10), ca.DM(20)], ca.DM(-1))
        DM(20)
        julia> ca.conditional(ca.DM(5), [ca.DM(10), ca.DM(20)], ca.DM(-1))
        DM(-1)
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y"); z = ca.SX.sym("z");
        julia> ca.conditional(x, [y, z], 0)
        SX(@1=(x==1), ((@1?z:0)+((!@1)?((x==0)?y:0):0)))
        julia> X = ca.MX.sym("X"); Y = ca.MX.sym("Y");
        julia> ca.conditional(X, [Y, 2 * Y], 0, true)
        MX(switch(X, Y){0})
        \enddoctest
        \endlang

        \sa if_else, Function::conditional

        \identifier{1cj} */
    inline friend MatType conditional(const MatType& ind, const std::vector<MatType> &x,
                                      const MatType &x_default, bool short_circuit=false) {
      return MatType::conditional(ind, x, x_default, short_circuit);
    }

    /** \brief Check if expression depends on the argument

        The argument must be symbolic (for MX: symbolic primitives as created by MX.sym,
        not elements of them).

        The dependency is in a mathematical way:
        Does the value of the argument affect the value of the expression?
        Equivalently, will the corresponding entry of the Jacobian by a non-zero?
        This is decided by sparsity propagation through the graph, so 0*x and x-x do
        not depend on x.

        The symbol could still be present in the expression graph.

        To check that kind of dependency, use contains(symvar(f),arg)

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> depends_on(x*y, x), depends_on(x*y, SX.sym("z"))
        (True, False)
        >>> depends_on(0*x, x), depends_on(x - x, x), depends_on(SX(3), x)
        (False, False, False)
        >>> depends_on(vertcat(x, 3), vertcat(x, y))
        True
        >>> v = MX.sym("v", 2)
        >>> depends_on(v[0], v)
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> fprintf('%d %d\n', depends_on(x*y, x), depends_on(x*y, SX.sym('z')))
        1 0
        >> fprintf('%d %d %d\n', depends_on(0*x, x), depends_on(x - x, x), depends_on(SX(3), x))
        0 0 0
        >> depends_on(vertcat(x, 3), vertcat(x, y))
        1
        >> v = MX.sym('v', 2);
        >> depends_on(v(1), v)
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.depends_on(x * y, x), ca.depends_on(x * y, ca.SX.sym("z"))
        (true, false)
        julia> ca.depends_on(0 * x, x), ca.depends_on(x - x, x), ca.depends_on(ca.SX(3), x)
        (false, false, false)
        julia> ca.depends_on([x; 3], [x; y])
        true
        julia> v = ca.MX.sym("v", 2);
        julia> ca.depends_on(v[1], v)
        true
        \enddoctest
        \endlang

        \sa which_depends, symvar, contains, jacobian_sparsity

        \identifier{1ck} */
    inline friend bool depends_on(const MatType& f, const MatType &arg) {
      return MatType::depends_on(f, arg);
    }

    /** \brief Check if expression n is listed in v
    *
    * This function checks for correspondence between nodes.
    * The detection does not descend into the expression graph.
    *
    * contains_all checks that every entry of n is listed in v, contains_any that at
    * least one is. For SX the entries must be scalar. Not defined for DM.
    *
    * For mathematical dependency, see depends_on
    * \sa depends_on
    *
    * \lang{python}
    * \doctest
    * >>> x = SX.sym("x"); y = SX.sym("y"); z = SX.sym("z")
    * >>> e = x + y
    * >>> contains([x, y], y), contains([e], y), contains([e], e)
    * (True, False, True)
    * >>> contains([SX.sym("x")], x)
    * False
    * >>> contains_all([x, y, z], [x, y]), contains_all([x, y], [x, z])
    * (True, False)
    * >>> contains_any([x, y], [x, z]), contains_any([x, y], [z])
    * (True, False)
    * >>> X = MX.sym("X"); Y = MX.sym("Y")
    * >>> contains([X, Y], Y), contains([X + Y], Y)
    * (True, False)
    * \enddoctest
    * \endlang
    *
    * \lang{matlab}
    * \doctest
    * >> x = SX.sym('x'); y = SX.sym('y');
    * >> e = x + y;
    * >> fprintf('%d %d %d\n', contains({x, y}, y), contains({e}, y), contains({e}, e))
    * 1 0 1
    * >> contains({SX.sym('x')}, x)
    * 0
    * >> X = MX.sym('X'); Y = MX.sym('Y');
    * >> fprintf('%d %d\n', contains({X, Y}, Y), contains({X + Y}, Y))
    * 1 0
    * \enddoctest
    * \endlang
    *
    * \lang{julia}
    * \doctest
    * julia> x = ca.SX.sym("x"); y = ca.SX.sym("y"); z = ca.SX.sym("z");
    * julia> e = x + y;
    * julia> ca.contains([x, y], y), ca.contains([e], y), ca.contains([e], e)
    * (true, false, true)
    * julia> ca.contains([ca.SX.sym("x")], x)
    * false
    * julia> ca.contains_all([x, y, z], [x, y]), ca.contains_all([x, y], [x, z])
    * (true, false)
    * julia> ca.contains_any([x, y], [x, z]), ca.contains_any([x, y], [z])
    * (true, false)
    * julia> X = ca.MX.sym("X"); Y = ca.MX.sym("Y");
    * julia> ca.contains([X, Y], Y), ca.contains([X + Y], Y)
    * (true, false)
    * \enddoctest
    * \endlang

        \identifier{2c8} */
    /// @{
    inline friend bool contains(const std::vector<MatType>& v, const MatType &n) {
      return contains_all(v, std::vector<MatType>{n});
    }

    inline friend bool contains_all(const std::vector<MatType>& v, const std::vector<MatType> &n) {
      return MatType::contains_all(v, n);
    }

    inline friend bool contains_any(const std::vector<MatType>& v, const std::vector<MatType> &n) {
      return MatType::contains_any(v, n);
    }
    /// @}

    /** \brief  Substitute variable v with expression vdef in an expression ex

        \p v must be symbolic and have the same shape as \p vdef (for SX also the same
        sparsity). For MX, \p v must be a symbolic primitive as created by MX.sym; an element
        of such a primitive is rejected. Substituting a symbol that does not occur in \p ex
        is a no-op. Not defined for DM.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y"); z = SX.sym("z")
        >>> substitute(x*y + x, x, 2)
        SX(((2.*y)+2))
        >>> substitute(x*y + x, x, z**2)
        SX(@1=sq(z), ((@1*y)+@1))
        >>> p = SX.sym("p", 2)
        >>> substitute(sumsqr(p), p, DM([1, 2]))
        SX(5)
        >>> X = MX.sym("X"); Y = MX.sym("Y")
        >>> substitute(X*Y + X, X, Y)
        MX((sq(Y)+Y))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z');
        >> substitute(x*y + x, x, 2)
        ((2.*y)+2)
        >> substitute(x*y + x, x, z^2)
        @1=sq(z), ((@1*y)+@1)
        >> p = SX.sym('p', 2);
        >> substitute(sumsqr(p), p, DM([1; 2]))
        5
        >> X = MX.sym('X'); Y = MX.sym('Y');
        >> substitute(X*Y + X, X, Y)
        (sq(Y)+Y)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y"); z = ca.SX.sym("z");
        julia> ca.substitute(x * y + x, x, 2)
        SX(((2.*y)+2))
        julia> ca.substitute(x * y + x, x, z^2)
        SX(@1=sq(z), ((@1*y)+@1))
        julia> p = ca.SX.sym("p", 2);
        julia> ca.substitute(ca.sumsqr(p), p, ca.DM([1, 2]))
        SX(5)
        julia> X = ca.MX.sym("X"); Y = ca.MX.sym("Y");
        julia> ca.substitute(X * Y + X, X, Y)
        MX((sq(Y)+Y))
        \enddoctest
        \endlang

        \sa substitute_inplace, graph_substitute, Function

        \identifier{1cl} */
    friend inline MatType substitute(const MatType& ex, const MatType& v,
                                     const MatType& vdef) {
      return MatType::substitute(ex, v, vdef);
    }

    /** \brief  Substitute variables v with expressions vdef in multiple expressions ex

        All substitutions happen simultaneously: vdef[k] is not itself subject to the
        other substitutions, so swapping two variables works as expected.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> substitute([x*y, x + y], [x, y], [1, 2])
        [SX(2), SX(3)]
        >>> substitute([x + y], [x, y], [y, x])
        [SX((x+y))]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> r = substitute({x*y, x + y}, {x, y}, {1, 2});
        >> r{1}, r{2}
        2
        3
        >> r = substitute({x + y}, {x, y}, {y, x});
        >> r{1}
        (x+y)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.substitute([x * y, x + y], [x, y], [1, 2])
        2-element Vector{SX}:
         SX(2)
         SX(3)
        julia> ca.substitute([x + y], [x, y], [y, x])
        1-element Vector{SX}:
         SX((x+y))
        \enddoctest
        \endlang

        \identifier{1cm} */
    friend inline std::vector<MatType>
      substitute(const std::vector<MatType>& ex, const std::vector<MatType>& v,
                 const std::vector<MatType>& vdef) {
      return MatType::substitute(ex, v, vdef);
    }

    /** \brief Inplace substitution with piggyback expressions

        Substitute variables \p v out of the expressions \p vdef sequentially,
        as well as out of a number of other expressions piggyback

        Interprets vdef[k] as the definition of the symbol v[k] and resolves the definitions
        into each other, so that afterwards neither \p vdef nor \p ex contain any of the
        symbols \p v. Both \p vdef and \p ex are modified in place (and returned, in Python).
        \p v must consist of symbolic primitives with the shapes of \p vdef.
        \p reverse is only accepted for SX; MX raises "Not implemented".

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> substitute_inplace([x, y], [y + 1, 3], [x*y])
        ([SX(4), SX(3)], [SX(12)])
        >>> X = MX.sym("X"); Y = MX.sym("Y")
        >>> substitute_inplace([X, Y], [Y + 1, 3], [X*Y])
        ([MX(1), MX(3)], [MX(00)])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> [vdef, ex] = substitute_inplace({x, y}, {y + 1, 3}, {x*y});
        >> vdef{1}, vdef{2}, ex{1}
        4
        3
        12
        >> X = MX.sym('X'); Y = MX.sym('Y');
        >> [vdef, ex] = substitute_inplace({X, Y}, {Y + 1, 3}, {X*Y});
        >> vdef{1}, vdef{2}, ex{1}
        1
        3
        00
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.substitute_inplace([x, y], [y + 1, ca.SX(3)], [x*y])
        (SX[4, 3], SX[12])
        julia> X = ca.MX.sym("X"); Y = ca.MX.sym("Y");
        julia> ca.substitute_inplace([X, Y], [Y + 1, ca.MX(3)], [X*Y])
        (MX[1, 3], MX[00])
        \enddoctest
        \endlang

        \sa substitute

        \identifier{1cn} */
    inline friend void
      substitute_inplace(const std::vector<MatType>& v,
                        std::vector<MatType>& inout_vdef,
                        std::vector<MatType>& inout_ex, bool reverse=false) {
      return MatType::substitute_inplace(v, inout_vdef, inout_ex, reverse);
    }

    /** \brief Common subexpression elimination

        Returns an equivalent expression in which identical subexpressions that were built
        separately share a single node, reducing the number of nodes and the evaluation
        cost. The shape is preserved. DM is returned unchanged.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> e = sin(x + y)*cos(x + y)
        >>> cse(e)
        SX(@1=(x+y), (sin(@1)*cos(@1)))
        >>> n_nodes(e), n_nodes(cse(e))
        (7, 6)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> e = sin(x + y)*cos(x + y);
        >> cse(e)
        @1=(x+y), (sin(@1)*cos(@1))
        >> fprintf('%d %d\n', n_nodes(e), n_nodes(cse(e)))
        7 6
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> e = sin(x + y) * cos(x + y);
        julia> ca.cse(e)
        SX(@1=(x+y), (sin(@1)*cos(@1)))
        julia> ca.n_nodes(e), ca.n_nodes(ca.cse(e))
        (7, 6)
        \enddoctest
        \endlang

        \sa n_nodes, shared, extract

        \identifier{1co} */
    inline friend MatType cse(const MatType& e) {
      return MatType::cse({e}).at(0);
    }


    /** \brief Common subexpression elimination over a list of expressions

        Subexpressions are shared across all entries of the list.

        \lang{python}
        \doctest
        >>> X = MX.sym("X"); Y = MX.sym("Y")
        >>> vertcat(sin(X + Y), cos(X + Y))
        MX(vertcat(sin((X+Y)), cos((X+Y))))
        >>> vertcat(*cse([sin(X + Y), cos(X + Y)]))
        MX(@1=(X+Y), vertcat(sin(@1), cos(@1)))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> X = MX.sym('X'); Y = MX.sym('Y');
        >> [sin(X + Y); cos(X + Y)]
        vertcat(sin((X+Y)), cos((X+Y)))
        >> r = cse({sin(X + Y), cos(X + Y)});
        >> vertcat(r{:})
        @1=(X+Y), vertcat(sin(@1), cos(@1))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> X = ca.MX.sym("X"); Y = ca.MX.sym("Y");
        julia> [sin(X + Y); cos(X + Y)]
        MX(vertcat(sin((X+Y)), cos((X+Y))))
        julia> vcat(ca.cse([sin(X + Y), cos(X + Y)])...)
        MX(@1=(X+Y), vertcat(sin(@1), cos(@1)))
        \enddoctest
        \endlang

        \identifier{1cp} */
    inline friend std::vector<MatType> cse(const std::vector<MatType>& e) {
      return MatType::cse(e);
    }

    /** \brief  Solve a system of equations: A*x = b

        The solve routine works similar to Matlab's backslash when A is square and nonsingular.
        The algorithm used is the following:
        1. A simple forward or backward substitution if A is upper or lower triangular
        2. If the linear system is at most 3-by-3, form the inverse via minor expansion and multiply
        3. Permute the variables and equations as to get a (structurally) nonzero diagonal,
        then perform a QR factorization without pivoting and solve the factorized system.

        Note 1: If there are entries of the linear system known to be zero, these will be removed.
        Elements that are very small, or will evaluate to be zero, can still cause numerical errors,
        due to the lack of pivoting (which is not possible since cannot compare the size of entries)

        Note 2: When permuting the linear system, a BLT (block lower triangular) transformation is
        formed. Only the permutation part of this is however used. An improvement would be to solve
        block-by-block if there are multiple BLT blocks.

        The description above applies to DM and SX. For MX, triangular systems give dedicated
        solve nodes, an orthonormal pattern is inverted by transposition, and everything else
        falls back to solve(A, b, "qr"). A scalar \p A gives the elementwise division b/A.
        \p b may have several columns (one right-hand side each). A non-square \p A or a row
        count mismatch raises an error; a numerically singular DM returns nan or inf.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> solve(A, DM([1, 2]))
        DM([0, 0.5])
        >>> solve(A, DM.eye(2))
        DM(
        [[-2, 1],
         [1.5, -0.5]])
        >>> solve(SX.sym("A", 2, 2), SX.sym("b", 2))
        SX(@1=((A_0*A_3)-(A_2*A_1)), [(((A_3/@1)*b_0)-((A_2/@1)*b_1)),
        (((A_0/@1)*b_1)-((A_1/@1)*b_0))])
        >>> solve(MX.sym("A", 2, 2), MX.sym("b", 2))
        MX((A\b))
        >>> solve(DM(2), DM([1, 2]))
        DM([0.5, 1])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> solve(A, DM([1; 2]))
        [0, 0.5]
        >> solve(A, DM.eye(2))
        [[-2, 1],
         [1.5, -0.5]]
        >> solve(SX.sym('A', 2, 2), SX.sym('b', 2))
        @1=((A_0*A_3)-(A_2*A_1)), [(((A_3/@1)*b_0)-((A_2/@1)*b_1)), (((A_0/@1)*b_1)-((A_1/@1)*b_0))]
        >> solve(MX.sym('A', 2, 2), MX.sym('b', 2))
        (A\b)
        >> solve(DM(2), DM([1; 2]))
        [0.5, 1]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> ca.solve(A, ca.DM([1, 2]))
        DM([0, 0.5])
        julia> ca.solve(A, ca.DM.eye(2))
        DM(
        [[-2, 1],
         [1.5, -0.5]])
        julia> ca.solve(ca.SX.sym("A", 2, 2), ca.SX.sym("b", 2))
        SX(@1=((A_0*A_3)-(A_2*A_1)), [(((A_3/@1)*b_0)-((A_2/@1)*b_1)), (((A_0/@1)*b_1)-((A_1/@1)*b_0))])
        julia> ca.solve(ca.MX.sym("A", 2, 2), ca.MX.sym("b", 2))
        MX((A\b))
        julia> ca.solve(ca.DM(2), ca.DM([1, 2]))
        DM([0.5, 1])
        \enddoctest
        \endlang

        \sa inv, mldivide, Linsol, qr

        \identifier{1cq} */
    friend inline MatType solve(const MatType& A, const MatType& b) {
      // If A is scalar, just divide
      if (A.is_scalar()) return b/A;
      return MatType::solve(A, b);
    }

    /** \brief Solve a system of equations: A*x = b, using a linear solver plugin

        Available for DM and MX (for MX the solve becomes a call to the solver, evaluated
        numerically); SX raises an error. A scalar \p A gives the elementwise division b/A.

        \param lsolver Name of a Linsol plugin, e.g. "qr", "csparse", "lapacklu", "lapackqr",
        "symbolicqr", "ldl", "ma27"; which ones are available depends on the build
        \param dict Options passed to the linear solver

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> solve(A, DM([1, 2]), "qr")
        DM([0, 0.5])
        >>> solve(MX.sym("A", 2, 2), MX.sym("b", 2), "qr")
        MX((A\b))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> solve(A, DM([1; 2]), 'qr')
        [0, 0.5]
        >> solve(MX.sym('A', 2, 2), MX.sym('b', 2), 'qr')
        (A\b)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> ca.solve(A, ca.DM([1, 2]), "qr")
        DM([0, 0.5])
        julia> ca.solve(ca.MX.sym("A", 2, 2), ca.MX.sym("b", 2), "qr")
        MX((A\b))
        \enddoctest
        \endlang

        \sa Linsol, has_linsol

        \identifier{1cr} */
    friend inline MatType solve(const MatType& A, const MatType& b,
                                const std::string& lsolver,
                                const Dict& dict = Dict()) {
      // If A is scalar, just divide
      if (A.is_scalar()) return b/A;
      return MatType::solve(A, b, lsolver, dict);
    }

#ifdef WITH_DEPRECATED_FEATURES
    /** \brief Linearize an expression
    *
    * [DEPRECATED] confusing behaviour
    *
    * This function linearizes an expression around a point x0,
    * using x as delta around that point:
    *
    *  linearize(sin(x),x,x0) -> sin(x0)+cos(x0)*x
    *
    * \p x0 must have the shape of \p x or be a scalar (broadcast).
    * For a variant that returns sin(x0)+cos(x0)*(x-x0), see taylor.
    *
    * \lang{python}
    * \doctest
    * >>> x = SX.sym("x"); y = SX.sym("y")
    * >>> linearize(sin(x), x, y)
    * SX((sin(y)+(cos(y)*x)))
    * \enddoctest
    * \endlang
    *
    * \lang{matlab}
    * \doctest
    * >> x = SX.sym('x'); y = SX.sym('y');
    * >> linearize(sin(x), x, y)
    * (sin(y)+(cos(y)*x))
    * \enddoctest
    * \endlang
    *
    * \lang{julia}
    * \doctest
    * julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
    * julia> ca.linearize(sin(x), x, y)
    * SX((sin(y)+(cos(y)*x)))
    * \enddoctest
    * \endlang
    *
    * \sa taylor, linear_coeff

        \identifier{1cs} */
    friend inline MatType linearize(const MatType& f, const MatType& x, const MatType& x0,
        const Dict& opts=Dict()) {
      return MatType::linearize(f, x, x0, opts);
    }
#endif // WITH_DEPRECATED_FEATURES

    /** \brief Computes the Moore-Penrose pseudo-inverse

        If the matrix A is fat (size1<size2), mul(A, pinv(A)) is unity.

          pinv(A)' = (AA')^(-1) A

        If the matrix A is slender (size1>size2), mul(pinv(A), A) is unity.

          pinv(A) = (A'A)^(-1) A'

        The normal equations are solved with solve, so \p A is assumed to have full rank.
        For a square matrix this is the inverse.

        \lang{python}
        \doctest
        >>> A = DM([[1, 0, 0], [0, 2, 0]])
        >>> pinv(A)
        DM(
        [[1, 0],
         [0, 0.5],
         [0, 0]])
        >>> mtimes(A, pinv(A))
        DM(
        [[1, 0],
         [0, 1]])
        >>> pinv(SX.sym("A", 1, 2))
        SX(@1=(sq(A_0)+sq(A_1)), [(A_0/@1), (A_1/@1)])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 0 0; 0 2 0]);
        >> pinv(A)
        [[1, 0],
         [0, 0.5],
         [0, 0]]
        >> A*pinv(A)
        [[1, 0],
         [0, 1]]
        >> pinv(SX.sym('A', 1, 2))
        @1=(sq(A_0)+sq(A_1)), [(A_0/@1), (A_1/@1)]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 0 0; 0 2 0]);
        julia> ca.pinv(A)
        DM(
        [[1, 0],
         [0, 0.5],
         [0, 0]])
        julia> A * ca.pinv(A)
        DM(
        [[1, 0],
         [0, 1]])
        julia> ca.pinv(ca.SX.sym("A", 1, 2))
        SX(@1=(sq(A_0)+sq(A_1)), [(A_0/@1), (A_1/@1)])
        \enddoctest
        \endlang

        \sa inv, solve, nullspace

        \identifier{1ct} */
    friend inline MatType pinv(const MatType& A) {
      return MatType::pinv(A);
    }

    /** \brief Computes the Moore-Penrose pseudo-inverse using a linear solver plugin

        Same as pinv(A) but the normal equations are solved with the given Linsol plugin.
        Available for DM and MX; SX raises an error.

        \param lsolver Name of a Linsol plugin, e.g. "qr"
        \param dict Options passed to the linear solver

        \lang{python}
        \doctest
        >>> pinv(DM([[1, 0, 0], [0, 2, 0]]), "qr")
        DM(
        [[1, 0],
         [0, 0.5],
         [0, 0]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> pinv(DM([1 0 0; 0 2 0]), 'qr')
        [[1, 0],
         [0, 0.5],
         [0, 0]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> print(ca.pinv(ca.DM([1 0 0; 0 2 0]), "qr"))
        [[1, 0],
         [0, 0.5],
         [0, 0]]
        \enddoctest
        \endlang

        \identifier{1cu} */
    friend inline MatType pinv(const MatType& A, const std::string& lsolver,
                               const Dict& dict = Dict()) {
      return MatType::pinv(A, lsolver, dict);
    }


    /** \brief Calculate Matrix exponential

        Computes expm(A*t) with A constant, through the 'slicot' expm plugin, which must be
        available in the build. Only defined for DM and MX (SX raises an error).

        \param A Square matrix
        \param t Scalar

        \sa expm, expmsol

        \identifier{23v} */
    friend inline MatType expm_const(const MatType& A, const MatType& t) {
      return MatType::expm_const(A, t);
    }

    /** \brief Calculate Matrix exponential

        Computes expm(A) for a square matrix through the 'slicot' expm plugin, which must be
        available in the build. Only defined for DM and MX (SX raises an error).

        \sa expm_const, expmsol, mpower

        \identifier{23w} */
    friend inline MatType expm(const MatType& A) {
      return MatType::expm(A);
    }

    /** \brief Calculate Jacobian

        The derivative of \p ex with respect to \p arg, with both treated as column-major
        vectors of their elements: the result is ex.numel()-by-arg.numel(), one row per
        element of \p ex and one column per element of \p arg. Entries that are structurally
        zero are left out of the sparsity pattern.

        \p arg must be symbolic (for MX: symbolic primitives as created by MX.sym, or a
        concatenation of them; elements of a primitive are rejected). Not defined for DM.

        \param opts Options: "allow_forward" and "allow_reverse" (bool, default true) restrict
        the AD mode; "symmetric" (bool, default false) hints that the Jacobian is symmetric;
        "verbose"
        \return Sparse matrix

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> J = jacobian(vertcat(x*y, x**2), vertcat(x, y))
        >>> J
        SX(
        [[y, x],
         [(2.*x), 00]])
        >>> J.dim(True)
        '2x2,3nz'
        >>> A = SX.sym("A", 2, 2)
        >>> jacobian(2*A, A).dim()
        '4x4'
        >>> jacobian(SX(0, 1), vertcat(x, y)).dim()
        '0x2'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> J = jacobian(vertcat(x*y, x^2), vertcat(x, y));
        >> J
        [[y, x],
         [(2.*x), 00]]
        >> J.dim(true)
        2x2,3nz
        >> A = SX.sym('A', 2, 2);
        >> J = jacobian(2*A, A); disp(J.dim())
        4x4
        >> J = jacobian(SX(0, 1), vertcat(x, y)); disp(J.dim())
        0x2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> J = ca.jacobian([x * y; x^2], [x; y])
        SX(
        [[y, x],
         [(2.*x), 00]])
        julia> ca.dim(J, true)
        "2x2,3nz"
        julia> A = ca.SX.sym("A", 2, 2);
        julia> ca.dim(ca.jacobian(2 * A, A))
        "4x4"
        julia> ca.dim(ca.jacobian(ca.SX(0, 1), [x; y]))
        "0x2"
        \enddoctest
        \endlang

        \sa gradient, hessian, jtimes, jacobian_sparsity, Function::jacobian

        \identifier{1cv} */
    inline friend MatType jacobian(const MatType &ex, const MatType &arg,
                                   const Dict& opts = Dict()) {
      return MatType::jacobian(ex, arg, opts);
    }
    /** \brief Calculate the gradient of an expression

        Computed with one reverse-mode sweep (jtimes with tr true) and projected onto the
        sparsity of \p arg: the result has the shape and pattern of \p arg, also for a matrix
        \p arg. A non-scalar \p ex raises an error (use jacobian).

        \param ex Scalar expression to take the gradient of
        \param arg Symbolic expression (any shape)
        \param opts Options, see jacobian
        \return Matrix with the sparsity of \p arg

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> gradient(x*y + x**2, vertcat(x, y))
        SX([((2.*x)+y), x])
        >>> gradient(x**2, vertcat(x, y))
        SX([(2.*x), 0])
        >>> A = SX.sym("A", 2, 2)
        >>> gradient(sumsqr(A), A)
        SX(
        [[(2.*A_0), (2.*A_2)],
         [(2.*A_1), (2.*A_3)]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> gradient(x*y + x^2, vertcat(x, y))
        [((2.*x)+y), x]
        >> gradient(x^2, vertcat(x, y))
        [(2.*x), 0]
        >> A = SX.sym('A', 2, 2);
        >> gradient(sumsqr(A), A)
        [[(2.*A_0), (2.*A_2)],
         [(2.*A_1), (2.*A_3)]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.gradient(x * y + x^2, [x; y])
        SX([((2.*x)+y), x])
        julia> ca.gradient(x^2, [x; y])
        SX([(2.*x), 0])
        julia> A = ca.SX.sym("A", 2, 2);
        julia> ca.gradient(ca.sumsqr(A), A)
        SX(
        [[(2.*A_0), (2.*A_2)],
         [(2.*A_1), (2.*A_3)]])
        \enddoctest
        \endlang

        \sa jacobian, hessian, jtimes

        \identifier{23x} */
    inline friend MatType gradient(const MatType &ex, const MatType &arg, const Dict& opts=Dict()) {
      return MatType::gradient(ex, arg, opts);
    }
    /** \brief Calculate the tangent of an expression

        Derivative of \p ex with respect to a scalar \p arg, computed with one forward-mode
        sweep and projected onto the sparsity of \p ex: the result has the shape and pattern
        of \p ex. A non-scalar \p arg raises an error (use jacobian).

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> tangent(vertcat(x**2, 3*x), x)
        SX([(2.*x), 3])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> tangent(vertcat(x^2, 3*x), x)
        [(2.*x), 3]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> ca.tangent([x^2; 3 * x], x)
        SX([(2.*x), 3])
        \enddoctest
        \endlang

        \sa gradient, jacobian, jtimes

        \identifier{23y} */
    inline friend MatType tangent(const MatType &ex, const MatType &arg, const Dict& opts=Dict()) {
      return MatType::tangent(ex, arg, opts);
    }

    /** \brief Calculate the Jacobian and multiply by a vector from the right

        This is equivalent to <tt>mul(jacobian(ex, arg), v)</tt> or
        <tt>mul(jacobian(ex, arg).T, v)</tt> for
        tr set to false and true respectively. If contrast to these
        expressions, it will use directional derivatives which is typically (but
        not necessarily) more efficient if the complete Jacobian is not needed and v has few rows.

        The matrix shapes of \p ex and \p arg are kept: \p v must have the shape of \p arg
        (or of \p ex when \p tr is true), and several directions can be passed at once by
        concatenating them horizontally, i.e. v.size2() a multiple of arg.size2(). The result
        is ex.size1()-by-ex.size2()*n_dir (arg.size1()-by-arg.size2()*n_dir when \p tr is
        true). An inconsistent shape of \p v raises an error.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> ex = vertcat(x*y, x**2); arg = vertcat(x, y)
        >>> jtimes(ex, arg, DM([1, 0]))
        SX([y, (2.*x)])
        >>> jtimes(ex, arg, DM([1, 0]), True)
        SX([y, x])
        >>> jtimes(ex, arg, DM.eye(2))
        SX(
        [[y, x],
         [(2.*x), 0]])
        >>> jtimes(ex, arg, SX.sym("v", 2))
        SX([((y*v_0)+(x*v_1)), ((2.*x)*v_0)])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> ex = vertcat(x*y, x^2); arg = vertcat(x, y);
        >> jtimes(ex, arg, DM([1; 0]))
        [y, (2.*x)]
        >> jtimes(ex, arg, DM([1; 0]), true)
        [y, x]
        >> jtimes(ex, arg, DM.eye(2))
        [[y, x],
         [(2.*x), 0]]
        >> jtimes(ex, arg, SX.sym('v', 2))
        [((y*v_0)+(x*v_1)), ((2.*x)*v_0)]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ex = [x * y; x^2]; arg = [x; y];
        julia> ca.jtimes(ex, arg, ca.DM([1, 0]))
        SX([y, (2.*x)])
        julia> ca.jtimes(ex, arg, ca.DM([1, 0]), true)
        SX([y, x])
        julia> ca.jtimes(ex, arg, ca.DM.eye(2))
        SX(
        [[y, x],
         [(2.*x), 0]])
        julia> ca.jtimes(ex, arg, ca.SX.sym("v", 2))
        SX([((y*v_0)+(x*v_1)), ((2.*x)*v_0)])
        \enddoctest
        \endlang

        \sa jacobian, forward, reverse

        \identifier{1cw} */
    friend inline MatType jtimes(const MatType &ex, const MatType &arg,
                                 const MatType &v, bool tr=false, const Dict& opts=Dict()) {
      return MatType::jtimes(ex, arg, v, tr, opts);
    }

    /** \brief Forward directional derivative

        For each seed direction v[d] (a list with one seed per entry of \p arg, with matching
        shapes), computes the forward sensitivities of every entry of \p ex: the result is a
        list over directions of lists over \p ex, containing jtimes(ex[i], arg, v[d]) summed
        over all arguments. Not defined for DM.

        \param ex List of expressions
        \param arg List of symbolic expressions
        \param v List of directions, each a list of seeds matching \p arg
        \param opts Options, see jacobian

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> forward([x*y, x**2], [x, y], [[1, 0]])
        [[SX(y), SX((2.*x))]]
        >>> forward([x*y, x**2], [x, y], [[1, 0], [0, 1]])
        [[SX(y), SX((2.*x))], [SX(x), SX(0)]]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> r = forward({x*y, x^2}, {x, y}, {{1, 0}});
        >> r{1}{1}, r{1}{2}
        y
        (2.*x)
        >> r = forward({x*y, x^2}, {x, y}, {{1, 0}, {0, 1}});
        >> r{2}{1}, r{2}{2}
        x
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.forward([x*y, x^2], [x, y], [[ca.SX(1), ca.SX(0)]])
        1-element Vector{Vector{SX}}:
         SX[y, (2.*x)]
        julia> ca.forward([x*y, x^2], [x, y], [[ca.SX(1), ca.SX(0)], [ca.SX(0), ca.SX(1)]])
        2-element Vector{Vector{SX}}:
         SX[y, (2.*x)]
         SX[x, 0]
        \enddoctest
        \endlang

        \sa reverse, jtimes, Function::forward

        \identifier{1cx} */
    friend inline std::vector<std::vector<MatType> >
    forward(const std::vector<MatType> &ex, const std::vector<MatType> &arg,
            const std::vector<std::vector<MatType> > &v,
            const Dict& opts = Dict()) {
      return MatType::forward(ex, arg, v, opts);
    }

    /** \brief Reverse directional derivative

        For each adjoint seed v[d] (a list with one seed per entry of \p ex, with matching
        shapes), computes the adjoint sensitivities with respect to every entry of \p arg:
        the result is a list over directions of lists over \p arg. Not defined for DM.

        \param ex List of expressions
        \param arg List of symbolic expressions
        \param v List of directions, each a list of seeds matching \p ex
        \param opts Options, see jacobian

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> reverse([x*y, x**2], [x, y], [[1, 0]])
        [[SX(y), SX(x)]]
        >>> reverse([x*y, x**2], [x, y], [[1, 1]])
        [[SX(((2.*x)+y)), SX(x)]]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> r = reverse({x*y, x^2}, {x, y}, {{1, 0}});
        >> r{1}{1}, r{1}{2}
        y
        x
        >> r = reverse({x*y, x^2}, {x, y}, {{1, 1}});
        >> r{1}{1}, r{1}{2}
        ((2.*x)+y)
        x
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.reverse([x*y, x^2], [x, y], [[ca.SX(1), ca.SX(0)]])
        1-element Vector{Vector{SX}}:
         SX[y, x]
        julia> ca.reverse([x*y, x^2], [x, y], [[ca.SX(1), ca.SX(1)]])
        1-element Vector{Vector{SX}}:
         SX[((2.*x)+y), x]
        \enddoctest
        \endlang

        \sa forward, jtimes, gradient, Function::reverse

        \identifier{1cy} */
    friend inline std::vector<std::vector<MatType> >
    reverse(const std::vector<MatType> &ex, const std::vector<MatType> &arg,
            const std::vector<std::vector<MatType> > &v,
            const Dict& opts = Dict()) {
      return MatType::reverse(ex, arg, v, opts);
    }

    ///@{
    /** \brief Hessian and (optionally) gradient

        Computes the gradient g = gradient(ex, arg) and the Hessian jacobian(g, arg), an
        arg.numel()-by-arg.numel() symmetric matrix; the "symmetric" option is set unless
        given in \p opts. \p ex must be scalar. The overload with \p output_g also returns
        the gradient; in Python, MATLAB and Julia this is the variant that is exposed, so the
        result is always the pair (H, g). Not defined for DM.

        \param opts Options, see jacobian

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> H, g = hessian(x**2*y, vertcat(x, y))
        >>> H
        SX(@1=(2.*x),
        [[(2.*y), @1],
         [@1, 00]])
        >>> g
        SX([((2.*x)*y), sq(x)])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> [H, g] = hessian(x^2*y, vertcat(x, y));
        >> H
        @1=(2.*x),
        [[(2.*y), @1],
         [@1, 00]]
        >> g
        [((2.*x)*y), sq(x)]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> H, g = ca.hessian(x^2 * y, [x; y]);
        julia> H
        SX(@1=(2.*x),
        [[(2.*y), @1],
         [@1, 00]])
        julia> g
        SX([((2.*x)*y), sq(x)])
        \enddoctest
        \endlang

        \sa gradient, jacobian, Function::hessian

        \identifier{23z} */
    inline friend MatType hessian(const MatType &ex, const MatType &arg,
        const Dict& opts = Dict()) {
      return MatType::hessian(ex, arg, opts);
    }
    inline friend MatType hessian(const MatType &ex, const MatType &arg, MatType& output_g,
        const Dict& opts = Dict()) {
      return MatType::hessian(ex, arg, output_g, opts);
    }
    ///@}

    /** \brief Find out which variables enter with some order

        Structural dependency check by sparsity propagation, without forming derivatives.
        With \p order 1, an entry depends if the corresponding Jacobian entry is structurally
        nonzero; with \p order 2, if it enters nonlinearly (second derivative structurally
        nonzero). Operations whose derivative is structurally zero, such as floor or sign,
        count as not depending.

        \param expr Expression
        \param var Symbolic expression
        \param order 1 or 2; other values raise an error
        \param tr When false, the result has one entry per element of \p var: does any
        element of \p expr depend on it. When true, one entry per element of \p expr: does
        it depend on any element of \p var.
        \return List of booleans

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> e = vertcat(x*y, x**2, 3*x)
        >>> which_depends(e, vertcat(x, y), 1, False)
        [True, True]
        >>> which_depends(e, vertcat(x, y), 1, True)
        [True, True, True]
        >>> which_depends(e, vertcat(x, y), 2, True)
        [True, True, False]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> e = vertcat(x*y, x^2, 3*x);
        >> which_depends(e, vertcat(x, y), 1, false)
        1 1
        >> which_depends(e, vertcat(x, y), 1, true)
        1 1 1
        >> which_depends(e, vertcat(x, y), 2, true)
        1 1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> e = [x * y; x^2; 3 * x];
        julia> ca.which_depends(e, [x; y], 1, false)
        2-element Vector{Bool}:
         1
         1
        julia> ca.which_depends(e, [x; y], 1, true)
        3-element Vector{Bool}:
         1
         1
         1
        julia> ca.which_depends(e, [x; y], 2, true)
        3-element Vector{Bool}:
         1
         1
         0
        \enddoctest
        \endlang

        \sa depends_on, is_linear, jacobian_sparsity

        \identifier{1cz} */
    inline friend std::vector<bool> which_depends(const MatType &expr, const MatType &var,
        casadi_int order, bool tr) {
      return MatType::which_depends(expr, var, order, tr);
    }

    /** \brief Get the sparsity pattern of a jacobian

        Equivalent to, but cheaper to compute than, jacobian(f,x).sparsity()

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> jacobian_sparsity(vertcat(x*y, x**2, 3), vertcat(x, y))
        Sparsity(3x2,3nz)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> jacobian_sparsity(vertcat(x*y, x^2, 3), vertcat(x, y))
        3x2,3nz
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.jacobian_sparsity([x * y; x^2; 3], [x; y])
        Sparsity(3x2,3nz)
        \enddoctest
        \endlang

        \sa jacobian, which_depends

        \identifier{259} */
    inline friend Sparsity jacobian_sparsity(const MatType &f, const MatType &x) {
      return MatType::jacobian_sparsity(f, x);
    }

    /** \brief Is expr linear in var?
    *
    * Defined as: no entry of \p expr depends on \p var with order 2 (see which_depends).
    * Only the structure of the graph is inspected, not the values:
    * false negatives are possible
    * (an expression may not be recognised as linear while it really is, e.g. x*x/x).
    * Operations whose derivative is structurally zero, such as floor, sign or comparisons,
    * are treated as linear.
    * A constant expression counts as linear.
    *
    * \lang{python}
    * \doctest
    * >>> x = SX.sym("x"); y = SX.sym("y")
    * >>> is_linear(2*x + y, vertcat(x, y)), is_linear(x*y, vertcat(x, y))
    * (True, False)
    * >>> is_linear(sin(y)*x, x), is_linear(x*x/x, x)
    * (True, False)
    * \enddoctest
    * \endlang
    *
    * \lang{matlab}
    * \doctest
    * >> x = SX.sym('x'); y = SX.sym('y');
    * >> fprintf('%d %d\n', is_linear(2*x + y, vertcat(x, y)), is_linear(x*y, vertcat(x, y)))
    * 1 0
    * >> fprintf('%d %d\n', is_linear(sin(y)*x, x), is_linear(x*x/x, x))
    * 1 0
    * \enddoctest
    * \endlang
    *
    * \lang{julia}
    * \doctest
    * julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
    * julia> ca.is_linear(2 * x + y, [x; y]), ca.is_linear(x * y, [x; y])
    * (true, false)
    * julia> ca.is_linear(sin(y) * x, x), ca.is_linear(x * x / x, x)
    * (true, false)
    * \enddoctest
    * \endlang
    *
    * \sa is_quadratic, linear_coeff, which_depends

        \identifier{1d0} */
    inline friend bool is_linear(const MatType &expr, const MatType &var) {
      return MatType::is_linear(expr, var);
    }

    /** \brief Is expr quadratic in var?
    *
    * Defined as: gradient(expr, var) is linear in \p var; hence \p expr must be scalar.
    * False negatives are possible
    * (an expression may not be recognised as quadratic while it really is),
    * see is_linear for the structural nature of the check.
    * Linear and constant expressions count as quadratic.
    *
    * \lang{python}
    * \doctest
    * >>> x = SX.sym("x"); y = SX.sym("y")
    * >>> is_quadratic(x*y + x**2, vertcat(x, y)), is_quadratic(x**3, x)
    * (True, False)
    * >>> is_quadratic(x, x)
    * True
    * \enddoctest
    * \endlang
    *
    * \lang{matlab}
    * \doctest
    * >> x = SX.sym('x'); y = SX.sym('y');
    * >> fprintf('%d %d\n', is_quadratic(x*y + x^2, vertcat(x, y)), is_quadratic(x^3, x))
    * 1 0
    * >> is_quadratic(x, x)
    * 1
    * \enddoctest
    * \endlang
    *
    * \lang{julia}
    * \doctest
    * julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
    * julia> ca.is_quadratic(x * y + x^2, [x; y]), ca.is_quadratic(x^3, x)
    * (true, false)
    * julia> ca.is_quadratic(x, x)
    * true
    * \enddoctest
    * \endlang
    *
    * \sa is_linear, quadratic_coeff

        \identifier{1d1} */
    inline friend bool is_quadratic(const MatType &expr, const MatType &var) {
      return MatType::is_quadratic(expr, var);
    }

    /** \brief Recognizes quadratic form in scalar expression
    *
    * 1/2*x' A x + b' x + c
    *
    * e = 0.5*bilin(A,x,x)+dot(b,x)+c
    *
    * Returns A = hessian(expr, var), b = jacobian(expr, var) at var = 0 as a column, and
    * c = expr at var = 0. A non-scalar \p expr raises an error.
    *
    * \param check[in] When true (default), A is checked to be independent of x and an
    *                  error is raised otherwise ("called on non-quadratic expression").
    *                  Provided to deal with false positive dependency checks.
    *
    * \lang{python}
    * \doctest
    * >>> x = SX.sym("x"); y = SX.sym("y")
    * >>> A, b, c = quadratic_coeff(x**2 + 2*x*y + 3*x + 4, vertcat(x, y))
    * >>> A
    * SX(@1=2,
    * [[@1, @1],
    *  [@1, 00]])
    * >>> b, c
    * (SX([3, 0]), SX(4))
    * >>> try: quadratic_coeff(x**3, x)
    * ... except Exception as e: print("error")
    * error
    * \enddoctest
    * \endlang
    *
    * \lang{matlab}
    * \doctest
    * >> x = SX.sym('x'); y = SX.sym('y');
    * >> [A, b, c] = quadratic_coeff(x^2 + 2*x*y + 3*x + 4, vertcat(x, y));
    * >> A
    * @1=2,
    * [[@1, @1],
    *  [@1, 00]]
    * >> b, c
    * [3, 0]
    * 4
    * >> try, [A, b, c] = quadratic_coeff(x^3, x); catch e, disp('error'); end
    * error
    * \enddoctest
    * \endlang
    *
    * \lang{julia}
    * \doctest
    * julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
    * julia> A, b, c = ca.quadratic_coeff(x^2 + 2 * x * y + 3 * x + 4, [x; y]);
    * julia> A
    * SX(@1=2,
    * [[@1, @1],
    *  [@1, 00]])
    * julia> b, c
    * ([3, 0], 4)
    * julia> try ca.quadratic_coeff(x^3, x) catch e; println("error") end
    * error
    * \enddoctest
    * \endlang
    *
    * \sa linear_coeff, is_quadratic, hessian

        \identifier{1d2} */
    inline friend void quadratic_coeff(const MatType &expr, const MatType &var,
        MatType& A, MatType& b, MatType& c, bool check=true) {
      MatType::quadratic_coeff(expr, var, A, b, c, check);
    }

    /** \brief Recognizes linear form in vector expression
    *
    * A x + b
    *
    * Returns A = jacobian(expr, var) at var = 0 and b = vec(expr) at var = 0.
    * \p expr must be a vector (a matrix raises an error).
    *
    * \param check[in] When true (default), \p expr is checked to be linear in \p var
    *                  (is_linear) and an error is raised otherwise.
    *                  Provided to deal with false positive dependency checks.
    *
    * \lang{python}
    * \doctest
    * >>> x = SX.sym("x"); y = SX.sym("y")
    * >>> A, b = linear_coeff(vertcat(2*x + y + 1, 3*y), vertcat(x, y))
    * >>> A
    * SX(
    * [[2, 1],
    *  [00, 3]])
    * >>> b
    * SX([1, 0])
    * \enddoctest
    * \endlang
    *
    * \lang{matlab}
    * \doctest
    * >> x = SX.sym('x'); y = SX.sym('y');
    * >> [A, b] = linear_coeff(vertcat(2*x + y + 1, 3*y), vertcat(x, y));
    * >> A
    * [[2, 1],
    *  [00, 3]]
    * >> b
    * [1, 0]
    * \enddoctest
    * \endlang
    *
    * \lang{julia}
    * \doctest
    * julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
    * julia> A, b = ca.linear_coeff([2 * x + y + 1; 3 * y], [x; y]);
    * julia> A
    * SX(
    * [[2, 1],
    *  [00, 3]])
    * julia> b
    * SX([1, 0])
    * \enddoctest
    * \endlang
    *
    * \sa quadratic_coeff, is_linear, jacobian

        \identifier{1d3} */
    inline friend void linear_coeff(const MatType &expr, const MatType &var,
        MatType& A, MatType& b, bool check=true) {
      MatType::linear_coeff(expr, var, A, b, check);
    }

    /** \brief Extract purely parametric parts from an expression graph
     *
     * The purpose of extract_parametric is ultimately to save on evaluation time of an expression,
     * by extracting out the parts that are only solely dependent on parameters.
     *
     * For any:
     * [expr_ret, symbols, parametric] = extract_parametric(expr, par)
     * It holds that:
     * substitute(expr_ret,symbols,parametric) == expr
     *
     * parametric is only dependant on par
     * expr_ret is not dependant on par, but is dependant on symbols
     *
     * Example:
     * [expr_ret, symbols, parametric] = extract_parametric((x-sqrt(p))*y+cos(p)**2, p)
     *
     * expr_ret: (((x-extracted1)*y)+extracted2)
     * symbols: [extracted1, extracted2]
     * parametric: [sqrt(p),cos(p)**2]
     *
     * Options to control the name of new symbols:
     *   - offset (int), default 0: Offset for new symbol counter
     *   - prefix (str), default 'e_': Prefix for new symbols
     *   - suffix (str), default '': Suffix for new symbols
     *
     * Option 'extract_trivial' controls the behaviour for extracted parts
     * that consist of pure symbols. If false (default),
     * these parts are not extracted.
     *
     * \p expr may also be a list of expressions and \p par a list of parameter symbols
     * (concatenated internally). An unknown option raises an error. Not defined for DM.
     *
     * \lang{python}
     * \doctest
     * >>> x = SX.sym("x"); y = SX.sym("y"); p = SX.sym("p")
     * >>> expr_ret, symbols, parametric = extract_parametric((x - sqrt(p))*y + cos(p)**2, p)
     * >>> expr_ret
     * SX((((x-e_0)*y)+e_1))
     * >>> symbols, parametric
     * ([SX(e_0), SX(e_1)], [SX(sqrt(p)), SX(sq(cos(p)))])
     * >>> extract_parametric(x*p, p)
     * (SX((x*p)), [], [])
     * >>> extract_parametric(x*p, p, {"extract_trivial": True})
     * (SX((x*e_0)), [SX(e_0)], [SX(p)])
     * \enddoctest
     * \endlang
     *
     * \lang{matlab}
     * \doctest
     * >> x = SX.sym('x'); y = SX.sym('y'); p = SX.sym('p');
     * >> [expr_ret, symbols, parametric] = extract_parametric((x - sqrt(p))*y + cos(p)^2, p);
     * >> expr_ret
     * (((x-e_0)*y)+e_1)
     * >> vertcat(symbols{:}), vertcat(parametric{:})
     * [e_0, e_1]
     * [sqrt(p), sq(cos(p))]
     * >> [expr_ret, symbols, parametric] = extract_parametric(x*p, p);
     * >> expr_ret, numel(symbols)
     * (x*p)
     * 0
     * >> opts = struct('extract_trivial', true);
     * >> [expr_ret, symbols, parametric] = extract_parametric(x*p, p, opts);
     * >> expr_ret, symbols{1}, parametric{1}
     * (x*e_0)
     * e_0
     * p
     * \enddoctest
     * \endlang
     *
     * \lang{julia}
     * \doctest
     * julia> x = ca.SX.sym("x"); y = ca.SX.sym("y"); p = ca.SX.sym("p");
     * julia> e = (x - sqrt(p)) * y + cos(p)^2;
     * julia> expr_ret, symbols, parametric = ca.extract_parametric(e, p);
     * julia> expr_ret
     * SX((((x-e_0)*y)+e_1))
     * julia> symbols, parametric
     * (SX[e_0, e_1], SX[sqrt(p), sq(cos(p))])
     * julia> ca.extract_parametric(x * p, p)
     * ((x*p), SX[], SX[])
     * julia> ca.extract_parametric(x * p, p, Dict("extract_trivial" => true))
     * ((x*e_0), SX[e_0], SX[p])
     * \enddoctest
     * \endlang
     *
     * \sa separate_linear, substitute

        \identifier{288} */
    inline friend void extract_parametric(const MatType &expr, const MatType& par,
        MatType& SWIG_OUTPUT(expr_ret),
        std::vector<MatType>& SWIG_OUTPUT(symbols),
        std::vector<MatType>& SWIG_OUTPUT(parametric),
        const Dict& opts=Dict()) {
      MatType::extract_parametric(expr, par, expr_ret, symbols, parametric, opts);
    }

    inline friend void extract_parametric(const std::vector<MatType> &expr, const MatType& par,
        std::vector<MatType>& SWIG_OUTPUT(expr_ret),
        std::vector<MatType>& SWIG_OUTPUT(symbols),
        std::vector<MatType>& SWIG_OUTPUT(parametric),
        const Dict& opts=Dict()) {
      // Concatenate all vector elements
      MatType expr_cat = veccat(expr);
      MatType expr_ret_cat;

      // Concatenated extract_parametric
      MatType::extract_parametric(expr_cat, par, expr_ret_cat, symbols, parametric, opts);

      // Compute edges of vertsplit needed to undo concatenate
      std::vector<casadi_int> edges = {0};
      for (const MatType& e : expr) {
        edges.push_back(edges.back() + e.numel());
      }
      // Perform vertsplit
      std::vector<MatType> expr_ret_catv = MatType::vertsplit(expr_ret_cat, edges);

      // Reshape all elements back into original size
      expr_ret.resize(expr_ret_catv.size());
      for (casadi_int i=0; i<expr_ret_catv.size(); ++i) {
        expr_ret[i] = reshape(expr_ret_catv[i], expr[i].size1(), expr[i].size2());
      }
    }

    inline friend void extract_parametric(const std::vector<MatType> &expr,
        const std::vector<MatType>& par,
        std::vector<MatType>& SWIG_OUTPUT(expr_ret),
        std::vector<MatType>& SWIG_OUTPUT(symbols),
        std::vector<MatType>& SWIG_OUTPUT(parametric),
        const Dict& opts=Dict()) {
      extract_parametric(expr, veccat(par), expr_ret, symbols, parametric, opts);
    }

    inline friend void extract_parametric(const MatType &expr, const std::vector<MatType>& par,
        MatType& SWIG_OUTPUT(expr_ret),
        std::vector<MatType>& SWIG_OUTPUT(symbols),
        std::vector<MatType>& SWIG_OUTPUT(parametric),
        const Dict& opts=Dict()) {
      extract_parametric(expr, veccat(par), expr_ret, symbols, parametric, opts);
    }

    /** \brief Separate an expression into subexpressions that are constant, linear and nonlinear
    *
    * \param expr[in] The expression to be separated
    * \param sym_lin[in] The symbolic variables w.r.t. which linearity should be checked
    * \param sym_const[in] The symbolic variables that are deemed constant
    * \param expr_const[out] The constant part of the expression
    * \param expr_lin[out] The linear part of the expression
    * \param expr_nonlin[out] The nonlinear part of the expression
    *
    * Every symbol in \p expr must be listed in \p sym_lin or \p sym_const; otherwise an
    * error about free variables is raised.
    *
    * A post condition is that the following holds:
    * expr = expr_const + expr_lin + expr_nonlin
    *
    * Here, expr_const is only dependant on sym_const,
    *       expr_lin is linear in sym_lin (with coefficients that may depend on sym_const)
    *
    * Example:
    *
    * [expr_const,expr_lin,expr_nonlin] =
    *   separate_linear(cos(p)+7*x+x*y, vertcat(x,y), p)
    *
    * expr_const: cos(p)
    * expr_lin: 7*x
    * expr_nonlin: x*y
    *
    * \p sym_lin and \p sym_const may also be lists of symbols. Not defined for DM.
    *
    * \lang{python}
    * \doctest
    * >>> x = SX.sym("x"); y = SX.sym("y"); p = SX.sym("p")
    * >>> separate_linear(cos(p) + 7*x + x*y, vertcat(x, y), p)
    * (SX(cos(p)), SX((7*x)), SX((x*y)))
    * >>> separate_linear(p*x, x, p)
    * (SX(0), SX((p*x)), SX(0))
    * \enddoctest
    * \endlang
    *
    * \lang{matlab}
    * \doctest
    * >> x = SX.sym('x'); y = SX.sym('y'); p = SX.sym('p');
    * >> [c, l, n] = separate_linear(cos(p) + 7*x + x*y, vertcat(x, y), p);
    * >> c, l, n
    * cos(p)
    * (7*x)
    * (x*y)
    * >> [c, l, n] = separate_linear(p*x, x, p);
    * >> c, l, n
    * 0
    * (p*x)
    * 0
    * \enddoctest
    * \endlang
    *
    * \lang{julia}
    * \doctest
    * julia> x = ca.SX.sym("x"); y = ca.SX.sym("y"); p = ca.SX.sym("p");
    * julia> ca.separate_linear(cos(p) + 7 * x + x * y, [x; y], p)
    * (cos(p), (7*x), (x*y))
    * julia> ca.separate_linear(p * x, x, p)
    * (0, (p*x), 0)
    * \enddoctest
    * \endlang
    *
    * \sa extract_parametric, linear_coeff

        \identifier{2mo} */
    inline friend void separate_linear(const MatType &expr,
      const MatType &sym_lin, const MatType &sym_const,
      MatType& expr_const, MatType& expr_lin, MatType& expr_nonlin) {
        MatType::separate_linear(expr, sym_lin, sym_const, expr_const, expr_lin, expr_nonlin);
    }

    inline friend void separate_linear(const MatType &expr,
      const std::vector<MatType> &sym_lin, const std::vector<MatType> &sym_const,
      MatType& expr_const, MatType& expr_lin, MatType& expr_nonlin) {
      separate_linear(expr, veccat(sym_lin), veccat(sym_const),
        expr_const, expr_lin, expr_nonlin);
    }

    /** \brief Count the number of nodes in an expression graph

        The number of nodes of the Function that evaluates the expression, which for SX
        includes symbols and constants and for MX also output nodes. Shared subexpressions
        are counted once; identical subexpressions built separately are counted separately
        (use cse to merge them). Not defined for DM.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> n_nodes(x), n_nodes(x + y), n_nodes(vertcat(x + y, x + y))
        (1, 3, 4)
        >>> n_nodes(SX(0, 1))
        0
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> fprintf('%d %d %d\n', n_nodes(x), n_nodes(x + y), n_nodes(vertcat(x + y, x + y)))
        1 3 4
        >> n_nodes(SX(0, 1))
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.n_nodes(x), ca.n_nodes(x + y), ca.n_nodes([x + y; x + y])
        (1, 3, 4)
        julia> ca.n_nodes(ca.SX(0, 1))
        0
        \enddoctest
        \endlang

        \sa cse, Function::n_nodes

        \identifier{2mp} */
    inline friend casadi_int n_nodes(const MatType& A) {
      return MatType::n_nodes(A);
    }

    /** \brief Simplify an expression

        For SX, applies simplification passes such as constant folding and removal of
        trivial operations. For MX and DM the expression is returned unchanged.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> simplify(x*2*3), simplify(x*0 + y*1)
        (SX((6*x)), SX(y))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> simplify(x*2*3), simplify(x*0 + y*1)
        (6*x)
        y
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.simplify(x * 2 * 3), ca.simplify(x * 0 + y * 1)
        ((6*x), y)
        \enddoctest
        \endlang

        \sa transform, cse

        \identifier{2mq} */
    friend inline MatType simplify(const MatType &x) {
      return MatType::simplify(x);
    }

    /** \brief Apply transformation passes to an expression

        Wraps the expression in a Function of its free symbols, applies Function::transform
        and re-evaluates it on the original symbols. Each pass is a list starting with a
        verb: "simplify" (followed by task names), "expand" or "external"; an unknown verb
        raises an error. Without \p passes a default set of simplification passes is applied.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> transform(x*2*3)
        SX((6*x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> transform(x*2*3)
        (6*x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> ca.transform(x * 2 * 3)
        SX((6*x))
        \enddoctest
        \endlang

        \sa simplify, Function::transform

        \identifier{2mr} */
    ///@{
    friend inline MatType transform(const MatType &x, const Dict& opts = Dict()) {
      return MatType::transform(x, opts);
    }
    friend inline MatType transform(const MatType &x,
        const std::vector<std::vector<GenericType> >& passes, const Dict& opts = Dict()) {
      return MatType::transform(x, passes, opts);
    }
    friend inline std::vector<MatType> transform(const std::vector<MatType> &x,
        const Dict& opts = Dict()) {
      return MatType::transform(x, opts);
    }
    friend inline std::vector<MatType> transform(const std::vector<MatType> &x,
        const std::vector<std::vector<GenericType> >& passes, const Dict& opts = Dict()) {
      return MatType::transform(x, passes, opts);
    }
    ///@}

    /** \brief Get a string representation for a binary MatType, using custom arguments

        Prints the operation at the root of a scalar expression with the given strings in
        place of its operands. \p xb must be a unary or binary operation (a symbol or
        constant raises an error for SX) and \p args must have one entry per operand.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> print_operator(x + y, ["a", "b"]), print_operator(sin(x), ["a"])
        ('(a+b)', 'sin(a)')
        >>> X = MX.sym("X"); Y = MX.sym("Y")
        >>> print_operator(vertcat(X, Y), ["a", "b"])
        'vertcat(a, b)'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> print_operator(x + y, {'a', 'b'}), print_operator(sin(x), {'a'})
        (a+b)
        sin(a)
        >> X = MX.sym('X'); Y = MX.sym('Y');
        >> print_operator(vertcat(X, Y), {'a', 'b'})
        vertcat(a, b)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.print_operator(x + y, ["a", "b"]), ca.print_operator(sin(x), ["a"])
        ("(a+b)", "sin(a)")
        julia> X = ca.MX.sym("X"); Y = ca.MX.sym("Y");
        julia> ca.print_operator([X; Y], ["a", "b"])
        "vertcat(a, b)"
        \enddoctest
        \endlang

        \identifier{1d4} */
    inline friend std::string
      print_operator(const MatType& xb, const std::vector<std::string>& args) {
      return MatType::print_operator(xb, args);
    }

    /** \brief Introduce intermediate variables for selected nodes in a graph

        Rewrites the expressions \p ex so that selected intermediate results are replaced by
        new symbols \p v with definitions \p vdef (ex, v and vdef are outputs; in Python the
        call is extract(ex, opts) -> (ex, v, vdef)). Substituting \p vdef for \p v in the
        result gives back the original expressions. Only nodes that are physically shared in
        the graph are detected: two separately constructed copies of x+y are distinct nodes
        (apply cse first to merge them). Not defined for DM.

        Options:
          - lift_shared (bool), default true: extract nodes used more than once
          - lift_calls (bool), default false: extract every Function call and its arguments
            (MX only; SX raises for lift_calls true or lift_shared false)
          - prefix (str), default 'v_', suffix (str), default '', offset (int), default 0:
            naming of the new symbols

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> e = x + y
        >>> extract([sin(e), cos(e)])
        ([SX(sin(v_0)), SX(cos(v_0))], [SX(v_0)], [SX((x+y))])
        >>> extract([sin(e), cos(e)], {"prefix": "w_", "offset": 3})
        ([SX(sin(w_3)), SX(cos(w_3))], [SX(w_3)], [SX((x+y))])
        >>> extract([sin(x + y), cos(x + y)])
        ([SX(sin((x+y))), SX(cos((x+y)))], [], [])
        >>> X = MX.sym("X"); Y = MX.sym("Y")
        >>> f = Function("f", [X], [sin(X)]); E = X + Y
        >>> extract([f(E)*f(E)], {"lift_shared": False, "lift_calls": True})
        ([MX((v_1*v_2))], [MX(v_0), MX(v_1), MX(v_2)], [MX((X+Y)), MX(f(v_0){0}), MX(f(v_0){0})])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> e = x + y;
        >> [ex, v, vdef] = extract({sin(e), cos(e)});
        >> ex{1}, ex{2}, v{1}, vdef{1}
        sin(v_0)
        cos(v_0)
        v_0
        (x+y)
        >> [ex, v, vdef] = extract({sin(e), cos(e)}, struct('prefix', 'w_', 'offset', 3));
        >> ex{1}, v{1}
        sin(w_3)
        w_3
        >> [ex, v, vdef] = extract({sin(x + y), cos(x + y)});
        >> ex{1}, numel(v)
        sin((x+y))
        0
        >> X = MX.sym('X'); Y = MX.sym('Y');
        >> f = Function('f', {X}, {sin(X)}); E = X + Y;
        >> opts = struct('lift_shared', false, 'lift_calls', true);
        >> [ex, v, vdef] = extract({f(E)*f(E)}, opts);
        >> ex{1}, numel(v), vdef{2}
        (v_1*v_2)
        3
        f(v_0){0}
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> e = x + y;
        julia> ca.extract([sin(e), cos(e)])
        (SX[sin(v_0), cos(v_0)], SX[v_0], SX[(x+y)])
        julia> ca.extract([sin(e), cos(e)], Dict("prefix" => "w_", "offset" => 3))
        (SX[sin(w_3), cos(w_3)], SX[w_3], SX[(x+y)])
        julia> ca.extract([sin(x + y), cos(x + y)])
        (SX[sin((x+y)), cos((x+y))], SX[], SX[])
        julia> X = ca.MX.sym("X"); Y = ca.MX.sym("Y");
        julia> f = ca.Function("f", [X], [sin(X)]); E = X + Y;
        julia> ca.extract([f(E) * f(E)], Dict("lift_shared" => false, "lift_calls" => true))
        (MX[(v_1*v_2)], MX[v_0, v_1, v_2], MX[(X+Y), f(v_0){0}, f(v_0){0}])
        \enddoctest
        \endlang

        \sa shared, cse, substitute

        \identifier{1d5} */
    inline friend void extract(std::vector<MatType>& ex,
        std::vector<MatType>& v,
        std::vector<MatType>& vdef,
        const Dict& opts = Dict()) {
      MatType::extract(ex, v, vdef, opts);
    }

    /** \brief Extract shared subexpressions from an set of expressions

        Same as extract with the options lift_shared true and lift_calls false: every node
        that is used more than once is replaced by a new symbol named \p v_prefix + counter
        + \p v_suffix. In Python the call is shared(ex, v_prefix, v_suffix) -> (ex, v, vdef).

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> e = x + y
        >>> shared([sin(e)*cos(e) + e])
        ([SX(((sin(v_0)*cos(v_0))+v_0))], [SX(v_0)], [SX((x+y))])
        >>> shared([sin(e), cos(e)], "w_", "_s")
        ([SX(sin(w_0_s)), SX(cos(w_0_s))], [SX(w_0_s)], [SX((x+y))])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> e = x + y;
        >> [ex, v, vdef] = shared({sin(e)*cos(e) + e});
        >> ex{1}, v{1}, vdef{1}
        ((sin(v_0)*cos(v_0))+v_0)
        v_0
        (x+y)
        >> [ex, v, vdef] = shared({sin(e), cos(e)}, 'w_', '_s');
        >> ex{1}, ex{2}, v{1}
        sin(w_0_s)
        cos(w_0_s)
        w_0_s
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> e = x + y;
        julia> ca.shared([sin(e) * cos(e) + e])
        (SX[((sin(v_0)*cos(v_0))+v_0)], SX[v_0], SX[(x+y)])
        julia> ca.shared([sin(e), cos(e)], "w_", "_s")
        (SX[sin(w_0_s), cos(w_0_s)], SX[w_0_s], SX[(x+y)])
        \enddoctest
        \endlang

        \sa extract, cse

        \identifier{1d6} */
    inline friend void shared(std::vector<MatType>& ex,
        std::vector<MatType>& v,
        std::vector<MatType>& vdef,
        const std::string& v_prefix="v_",
        const std::string& v_suffix="") {
      return MatType::shared(ex, v, vdef, v_prefix, v_suffix);
    }

    /** \brief Given a repeated matrix, computes the sum of repeated parts

        Splits \p A into an n-by-m grid of equally sized blocks and returns their sum, a
        matrix of size A.size1()/n by A.size2()/m. The inverse of repmat(x, n, m) in the
        sense that repsum(repmat(x, n, m), n, m) == n*m*x. The dimensions of \p A must be
        multiples of \p n and \p m.

        \lang{python}
        \doctest
        >>> repsum(DM([[1, 2, 3, 4], [5, 6, 7, 8]]), 1, 2)
        DM(
        [[4, 6],
         [12, 14]])
        >>> repsum(SX.sym("x", 4), 2)
        SX([(x_0+x_2), (x_1+x_3)])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> repsum(DM([1 2 3 4; 5 6 7 8]), 1, 2)
        [[4, 6],
         [12, 14]]
        >> repsum(SX.sym('x', 4), 2)
        [(x_0+x_2), (x_1+x_3)]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.repsum(ca.DM([1 2 3 4; 5 6 7 8]), 1, 2)
        DM(
        [[4, 6],
         [12, 14]])
        julia> ca.repsum(ca.SX.sym("x", 4), 2)
        SX([(x_0+x_2), (x_1+x_3)])
        \enddoctest
        \endlang

        \sa repmat, blocksplit

        \identifier{1d7} */
    inline friend MatType repsum(const MatType &A, casadi_int n, casadi_int m=1) {
      return MatType::repsum(A, n, m);
    }

    ///@{
    /** \brief Smallest element in a matrix

        Over all entries, structural zeros included (so a sparse matrix with only positive
        nonzeros has minimum 0). An empty matrix gives an empty 0-by-0 result. For SX a
        chain of fmin operations is built; for MX a single node.

        \lang{python}
        \doctest
        >>> mmin(DM([[3, -1], [2, 5]])), mmax(DM([[3, -1], [2, 5]]))
        (DM(-1), DM(5))
        >>> mmin(DM(Sparsity.diag(2), [3, 4]))
        DM(0)
        >>> mmin(SX.sym("x", 3))
        SX(fmin(fmin(x_0,x_1),x_2))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> mmin(DM([3 -1; 2 5])), mmax(DM([3 -1; 2 5]))
        -1
        5
        >> mmin(DM(Sparsity.diag(2), [3; 4]))
        0
        >> mmin(SX.sym('x', 3))
        fmin(fmin(x_0,x_1),x_2)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.mmin(ca.DM([3 -1; 2 5])), ca.mmax(ca.DM([3 -1; 2 5]))
        (-1, 5)
        julia> minimum(ca.DM([3 -1; 2 5])), maximum(ca.DM([3 -1; 2 5]))
        (-1, 5)
        julia> ca.mmin(ca.DM(ca.Sparsity.diag(2), [3, 4]))
        DM(0)
        julia> ca.mmin(ca.SX.sym("x", 3))
        SX(fmin(fmin(x_0,x_1),x_2))
        \enddoctest
        \endlang

        \sa mmax, fmin, norm_inf

        \identifier{1d8} */
    friend inline MatType mmin(const MatType& x) {
      return MatType::mmin(x);
    }
    ///@}

    ///@{
    /** \brief Largest element in a matrix

        Over all entries, structural zeros included. See mmin.

        \lang{python}
        \doctest
        >>> mmax(DM(Sparsity.diag(2), [-3, -4]))
        DM(0)
        >>> mmax(MX.sym("x", 3))
        MX(max(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> mmax(DM(Sparsity.diag(2), [-3; -4]))
        0
        >> mmax(MX.sym('x', 3))
        max(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.mmax(ca.DM(ca.Sparsity.diag(2), [-3, -4]))
        DM(0)
        julia> ca.mmax(ca.MX.sym("x", 3))
        MX(max(x))
        \enddoctest
        \endlang

        \sa mmin, fmax, logsumexp

        \identifier{1d9} */
    friend inline MatType mmax(const MatType& x) {
      return MatType::mmax(x);
    }
    ///@}

    ///@{
    /// Functions called by friend functions defined here
    static MatType jtimes(const MatType &ex, const MatType &arg,
                          const MatType &v, bool tr=false, const Dict& opts=Dict());
    static MatType gradient(const MatType &ex, const MatType &arg, const Dict& opts=Dict());
    static MatType tangent(const MatType &ex, const MatType &arg, const Dict& opts=Dict());
    static MatType linearize(const MatType& f, const MatType& x, const MatType& x0,
      const Dict& opts=Dict());
    static MatType mpower(const MatType &x, const MatType &y);
    static MatType soc(const MatType &x, const MatType &y);
    ///@}

/** @} */
#endif // SWIG

    /** @name Construct symbolic primitives
        The "sym" function is intended to work in a similar way as "sym" used
        in the Symbolic Toolbox for Matlab but instead creating a
        CasADi symbolic primitive.
    */
    ///@{

    /** \brief Create an nrow-by-ncol symbolic primitive

        Creates a new, dense symbol; two symbols with the same name are still different
        symbols. For SX each entry is a separate scalar symbol named name_k (k the
        column-major index), while an MX symbol is a single node with the given name.
        Only SX and MX have symbols; DM.sym raises an error.

        \param name Name used for printing. For SX, a name of the form "[a,b,c]" assigns
        individual names to the numel() entries.
        \param nrow Number of rows, default 1
        \param ncol Number of columns, default 1

        \lang{python}
        \doctest
        >>> SX.sym("x"), MX.sym("x", 2, 3).dim()
        (SX(x), '2x3')
        >>> SX.sym("x", 2, 2)
        SX(
        [[x_0, x_2],
         [x_1, x_3]])
        >>> SX.sym("[a,b,c]", 3)
        SX([a, b, c])
        >>> is_equal(SX.sym("x"), SX.sym("x"))
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> SX.sym('x'), MX.sym('x', 2, 3).dim()
        x
        2x3
        >> SX.sym('x', 2, 2)
        [[x_0, x_2],
         [x_1, x_3]]
        >> SX.sym('[a,b,c]', 3)
        [a, b, c]
        >> is_equal(SX.sym('x'), SX.sym('x'))
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.SX.sym("x"), ca.dim(ca.MX.sym("x", 2, 3))
        (x, "2x3")
        julia> ca.SX.sym("x", 2, 2)
        SX(
        [[x_0, x_2],
         [x_1, x_3]])
        julia> ca.SX.sym("[a,b,c]", 3)
        SX([a, b, c])
        julia> isequal(ca.SX.sym("x"), ca.SX.sym("x"))
        false
        \enddoctest
        \endlang

        \sa symvar, zeros, ones

        \identifier{1da} */
    static MatType sym(const std::string& name, casadi_int nrow=1, casadi_int ncol=1) {
      return sym(name, Sparsity::dense(nrow, ncol));
    }

    /** \brief  Construct a symbolic primitive with given dimensions, as a pair (nrow, ncol)

        \identifier{1db} */
    static MatType sym(const std::string& name, const std::pair<casadi_int, casadi_int> &rc) {
      return sym(name, rc.first, rc.second);
    }

    /** \brief Create symbolic primitive with a given sparsity pattern

        Only the structural nonzeros of \p sp become symbolic entries; for SX they are
        numbered by nonzero index.

        \lang{python}
        \doctest
        >>> SX.sym("p", Sparsity.diag(2))
        SX(
        [[p_0, 00],
         [00, p_1]])
        >>> MX.sym("p", Sparsity.diag(2)).dim(True)
        '2x2,2nz'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> SX.sym('p', Sparsity.diag(2))
        [[p_0, 00],
         [00, p_1]]
        >> MX.sym('p', Sparsity.diag(2)).dim(true)
        2x2,2nz
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.SX.sym("p", ca.Sparsity.diag(2))
        SX(
        [[p_0, 00],
         [00, p_1]])
        julia> ca.dim(ca.MX.sym("p", ca.Sparsity.diag(2)), true)
        "2x2,2nz"
        \enddoctest
        \endlang

        \identifier{1dc} */
    static MatType sym(const std::string& name, const Sparsity& sp) {
      return MatType::_sym(name, sp);
    }

    /** \brief Create a list of length p of symbolic primitives with a given sparsity

        The k-th symbol is named name followed by k (no separator).

        \lang{python}
        \doctest
        >>> MX.sym("x", Sparsity.dense(2, 1), 3)
        [MX(x0), MX(x1), MX(x2)]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = MX.sym('x', Sparsity.dense(2, 1), 3);
        >> numel(v), v{1}, v{3}
        3
        x0
        x2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.MX.sym("x", ca.Sparsity.dense(2, 1), 3)
        3-element Vector{MX}:
         MX(x0)
         MX(x1)
         MX(x2)
        \enddoctest
        \endlang

        \identifier{1dd} */
    static std::vector<MatType > sym(const std::string& name, const Sparsity& sp, casadi_int p);

    /** \brief Create a list of length p of nrow-by-ncol symbolic primitives

        \lang{python}
        \doctest
        >>> SX.sym("x", 1, 1, 2)
        [SX(x0), SX(x1)]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = SX.sym('x', 1, 1, 2);
        >> v{1}, v{2}
        x0
        x1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.SX.sym("x", 1, 1, 2)
        2-element Vector{SX}:
         SX(x0)
         SX(x1)
        \enddoctest
        \endlang

        \identifier{1de} */
    static std::vector<MatType > sym(const std::string& name, casadi_int nrow,
        casadi_int ncol, casadi_int p) {
      return sym(name, Sparsity::dense(nrow, ncol), p);
    }

    /** \brief Create a list of length r of lists of length p of symbolic primitives

        with given sparsity. The symbol at outer index i and inner index k is named
        name_ik.

        \identifier{1df} */
    static std::vector<std::vector<MatType> >
      sym(const std::string& name, const Sparsity& sp, casadi_int p, casadi_int r);

    /** \brief Create a list of length r of lists of length p of nrow-by-ncol symbolic primitives

        \lang{python}
        \doctest
        >>> MX.sym("x", 1, 1, 2, 2)
        [[MX(x_00), MX(x_01)], [MX(x_10), MX(x_11)]]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = MX.sym('x', 1, 1, 2, 2);
        >> numel(v), v{1}{2}, v{2}{1}
        2
        x_01
        x_10
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.MX.sym("x", 1, 1, 2, 2)
        2-element Vector{Vector{MX}}:
         MX[x_00, x_01]
         MX[x_10, x_11]
        \enddoctest
        \endlang

        \identifier{1dg} */
    static std::vector<std::vector<MatType> >
      sym(const std::string& name, casadi_int nrow, casadi_int ncol, casadi_int p, casadi_int r) {
      return sym(name, Sparsity::dense(nrow, ncol), p, r);
    }
    ///@}

    ///@{
    /** \brief Create a dense matrix or a matrix with specified sparsity with all entries zero

        zeros(nrow, ncol) stores nrow*ncol numerical zeros; zeros(sp) stores a zero at every
        structural nonzero of \p sp. To get structural zeros instead, construct from a shape,
        e.g. DM(2, 2).

        \lang{python}
        \doctest
        >>> DM.zeros(2, 2).nnz(), DM.zeros(Sparsity.diag(2)).nnz(), DM(2, 2).nnz()
        (4, 2, 0)
        >>> MX.zeros(2, 3)
        MX(zeros(2x3))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Z = DM.zeros(Sparsity.diag(2));
        >> fprintf('%d %d %d\n', DM.zeros(2, 2).nnz(), Z.nnz(), DM(2, 2).nnz())
        4 2 0
        >> MX.zeros(2, 3)
        zeros(2x3)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> Z = ca.DM.zeros(ca.Sparsity.diag(2));
        julia> ca.nnz(ca.DM.zeros(2, 2)), ca.nnz(Z), ca.nnz(ca.DM(2, 2))
        (4, 2, 0)
        julia> ca.MX.zeros(2, 3)
        MX(zeros(2x3))
        \enddoctest
        \endlang

        \sa ones, Sparsity

        \identifier{1dh} */
    static MatType zeros(casadi_int nrow=1, casadi_int ncol=1) {
      return zeros(Sparsity::dense(nrow, ncol));
    }
    static MatType zeros(const Sparsity& sp) { return MatType(sp, 0, false);}
    static MatType zeros(const std::pair<casadi_int, casadi_int>& rc) {
      return zeros(rc.first, rc.second);
    }
    ///@}

    ///@{
    /** \brief Create a dense matrix or a matrix with specified sparsity with all entries one

        ones(sp) sets every structural nonzero of \p sp to one; other entries stay
        structurally zero.

        \lang{python}
        \doctest
        >>> DM.ones(2, 3)
        DM(
        [[1, 1, 1],
         [1, 1, 1]])
        >>> DM.ones(Sparsity.diag(2))
        DM(
        [[1, 00],
         [00, 1]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM.ones(2, 3)
        [[1, 1, 1],
         [1, 1, 1]]
        >> DM.ones(Sparsity.diag(2))
        [[1, 00],
         [00, 1]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM.ones(2, 3)
        DM(
        [[1, 1, 1],
         [1, 1, 1]])
        julia> ca.DM.ones(ca.Sparsity.diag(2))
        DM(
        [[1, 00],
         [00, 1]])
        \enddoctest
        \endlang

        \sa zeros, Sparsity

        \identifier{1di} */
    static MatType ones(casadi_int nrow=1, casadi_int ncol=1) {
      return ones(Sparsity::dense(nrow, ncol));
    }
    static MatType ones(const Sparsity& sp) { return MatType(sp, 1, false);}
    static MatType ones(const std::pair<casadi_int, casadi_int>& rc) {
      return ones(rc.first, rc.second);
    }
    ///@}
  };

  // Throw informative error message
  #define CASADI_THROW_ERROR(FNAME, WHAT) \
  throw CasadiException("Error in " + MatType::type_name() \
    + "::" FNAME " at " + CASADI_WHERE + ":\n" + std::string(WHAT));

#ifndef SWIG
  // Implementations
  template<typename MatType>
  const Sparsity& GenericMatrix<MatType>::sparsity() const {
    return self().sparsity();
  }

  template<typename MatType>
  casadi_int GenericMatrix<MatType>::nnz() const {
    return sparsity().nnz();
  }

  template<typename MatType>
  casadi_int GenericMatrix<MatType>::nnz_lower() const {
    return sparsity().nnz_lower();
  }

  template<typename MatType>
  casadi_int GenericMatrix<MatType>::nnz_upper() const {
    return sparsity().nnz_upper();
  }

  template<typename MatType>
  casadi_int GenericMatrix<MatType>::nnz_diag() const {
    return sparsity().nnz_diag();
  }

  template<typename MatType>
  casadi_int GenericMatrix<MatType>::numel() const {
    return sparsity().numel();
  }

  template<typename MatType>
  casadi_int GenericMatrix<MatType>::size1() const {
    return sparsity().size1();
  }

  template<typename MatType>
  casadi_int GenericMatrix<MatType>::size2() const {
    return sparsity().size2();
  }

  template<typename MatType>
  std::pair<casadi_int, casadi_int> GenericMatrix<MatType>::size() const {
    return sparsity().size();
  }

  template<typename MatType>
  casadi_int GenericMatrix<MatType>::size(casadi_int axis) const {
    return sparsity().size(axis);
  }

  template<typename MatType>
  std::string GenericMatrix<MatType>::dim(bool with_nz) const {
    return sparsity().dim(with_nz);
  }

  template<typename MatType>
  bool GenericMatrix<MatType>::is_scalar(bool scalar_and_dense) const {
    return sparsity().is_scalar(scalar_and_dense);
  }

#endif // SWIG

  template<typename MatType>
  std::vector<MatType> GenericMatrix<MatType>::sym(const std::string& name,
                                                   const Sparsity& sp, casadi_int p) {
    std::vector<MatType> ret(p);
    std::stringstream ss;
    for (casadi_int k=0; k<p; ++k) {
      ss.str("");
      ss << name << k;
      ret[k] = sym(ss.str(), sp);
    }
    return ret;
  }

  template<typename MatType>
  std::vector<std::vector<MatType> > GenericMatrix<MatType>::sym(const std::string& name,
                                                                  const Sparsity& sp, casadi_int p,
                                                                  casadi_int r) {
    std::vector<std::vector<MatType> > ret(r);
    for (casadi_int k=0; k<r; ++k) {
      std::stringstream ss;
      ss << name << "_" << k;
      ret[k] = sym(ss.str(), sp, p);
    }
    return ret;
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::linspace(const MatType& a, const MatType& b, casadi_int nsteps) {
    casadi_assert(nsteps>=1, "linspace: nsteps (" + str(nsteps) + ") must be at least 1");
    std::vector<MatType> ret(nsteps);
    ret[0] = a;
    MatType step = (b-a)/static_cast<MatType>(nsteps-1);

    for (casadi_int i=1; i<nsteps-1; ++i)
      ret[i] = a + i * step;

    ret[nsteps-1] = b;
    return vertcat(ret);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::cross(const MatType& a, const MatType& b, casadi_int dim) {
    casadi_assert(a.size1()==b.size1() && a.size2()==b.size2(),
      "cross(a, b): Inconsistent dimensions. Dimension of a ("
      + a.dim() + " ) must equal that of b (" + b.dim() + ").");

    casadi_assert(a.size1()==3 || a.size2()==3,
      "cross(a, b): One of the dimensions of a should have length 3, but got "
        + a.dim() + ".");
    casadi_assert(dim==-1 || dim==1 || dim==2,
      "cross(a, b, dim): Dim must be 1, 2 or -1 (automatic).");

    std::vector<MatType> ret(3);

    bool t = a.size1()==3;

    if (dim==1) t = true;
    if (dim==2) t = false;

    MatType a1 = t ? a(0, Slice()) : a(Slice(), 0);
    MatType a2 = t ? a(1, Slice()) : a(Slice(), 1);
    MatType a3 = t ? a(2, Slice()) : a(Slice(), 2);

    MatType b1 = t ? b(0, Slice()) : b(Slice(), 0);
    MatType b2 = t ? b(1, Slice()) : b(Slice(), 1);
    MatType b3 = t ? b(2, Slice()) : b(Slice(), 2);

    ret[0] = a2*b3-a3*b2;
    ret[1] = a3*b1-a1*b3;
    ret[2] = a1*b2-a2*b1;

    return t ? vertcat(ret) : horzcat(ret);
  }

  double CASADI_EXPORT index_interp1d(const std::vector<double>& x, double xq,
    bool equidistant=false);

  template<typename MatType>
  MatType GenericMatrix<MatType>::interp1d(const std::vector<double>& x, const MatType& v,
      const std::vector<double>& xq, const std::string& mode, bool equidistant) {

    bool mode_floor = false;
    bool mode_ceil = false;
    if (mode=="floor") {
      mode_floor = true;
    } else if (mode=="ceil") {
      mode_ceil = true;
    } else if (mode=="linear") {
      //
    } else {
      casadi_error("interp1d(x, v, xq, mode): "
        "Mode must be 'floor', 'ceil' or 'linear'. Got '" + mode + "' instead.");
    }

    casadi_assert(is_increasing(x), "interp1d(x, v, xq): x must be increasing.");

    casadi_assert(x.size()==v.size1(),
      "interp1d(x, v, xq): dimensions mismatch. v expected to have " + str(x.size()) + " rows,"
      " but got " + str(v.size1()) + " instead.");

    // Need at least two elements
    casadi_assert(x.size()>=2, "interp1d(x, v, xq): x must be at least length 2.");

    // Vectors to compose a sparse matrix
    std::vector<double> val;
    std::vector<casadi_int> colind(1, 0);
    std::vector<casadi_int> row;

    // Number of nonzeros in to-be composed matrix
    casadi_int nnz = 0;
    for (casadi_int i=0;i<xq.size();++i) {
      // Obtain index corresponding to xq[i]
      double ind = index_interp1d(x, xq[i], equidistant);

      if (mode_floor) ind = floor(ind);
      if (mode_ceil) ind = ceil(ind);

      // Split into integer and fractional part
      double int_partd;
      double frac_part = modf(ind, &int_partd);
      casadi_int int_part = static_cast<casadi_int>(int_partd);

      if (frac_part==0) {
         // Create a single entry
         val.push_back(1);
         row.push_back(int_part);
         nnz+=1;
         colind.push_back(nnz);
       } else {
         // Create a double entry
         val.push_back(1-frac_part);
         val.push_back(frac_part);
         row.push_back(int_part);
         row.push_back(int_part+1);
         nnz+=2;
         colind.push_back(nnz);
       }
     }

     // Construct sparsity for composed matrix
     Sparsity sp(x.size(), xq.size() , colind, row);

     return MatType::mtimes(MatType(sp, val).T(), v);

  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::skew(const MatType& a) {
    casadi_assert(a.is_vector() && (a.size1()==3 || a.size2()==3),
      "skew(a): Expecting 3-vector, got " + a.dim() + ".");

    MatType x = a(0);
    MatType y = a(1);
    MatType z = a(2);
    return blockcat(std::vector< std::vector<MatType> >({{0, -z, y}, {z, 0, -x}, {-y, x, 0}}));
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::inv_skew(const MatType& a) {
    casadi_assert(a.size1()==3 && a.size2()==3,
      "inv_skew(a): Expecting 3-by-3 matrix, got " + a.dim() + ".");

    return 0.5*vertcat(std::vector<MatType>({a(2, 1)-a(1, 2), a(0, 2)-a(2, 0), a(1, 0)-a(0, 1)}));
  }


  template<typename MatType>
  MatType GenericMatrix<MatType>::tril2symm(const MatType& x) {
    casadi_assert(x.is_square(),
      "Shape error in tril2symm. Expecting square shape but got " + x.dim());
    casadi_assert(x.nnz_upper()-x.nnz_diag()==0,
      "Sparsity error in tril2symm. Found above-diagonal entries in argument: " + x.dim());
    return x +  x.T() - diag(diag(x));
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::repsum(const MatType& x, casadi_int n, casadi_int m) {
    casadi_assert(n>0 && x.size1() % n==0,
      "repsum: x.size1() (" + str(x.size1()) + ") must be a multiple of n (" + str(n) + ")");
    casadi_assert(m>0 && x.size2() % m==0,
      "repsum: x.size2() (" + str(x.size2()) + ") must be a multiple of m (" + str(m) + ")");
    std::vector< std::vector< MatType> > s =
      blocksplit(x, x.size1()/n, x.size2()/m);
    MatType sum = 0;
    for (casadi_int i=0;i<s.size();++i) {
      for (casadi_int j=0;j<s[i].size();++j) {
        sum = sum + s[i][j];
      }
    }
    return sum;
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::triu2symm(const MatType& x) {
    casadi_assert(x.is_square(),
      "Shape error in triu2symm. Expecting square shape but got " + x.dim());
    casadi_assert(x.nnz_lower()-x.nnz_diag()==0,
      "Sparsity error in triu2symm. Found below-diagonal entries in argument: " + x.dim());
    return x + x.T() - diag(diag(x));
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::bilin(const MatType& A, const MatType& x,
                                        const MatType& y) {
    // Check/correct x
    casadi_assert_dev(x.is_vector());
    if (!x.is_column()) return bilin(A, x.T(), y);
    if (!x.is_dense()) return bilin(A, densify(x), y);

    // Check/correct y
    casadi_assert_dev(y.is_vector());
    if (!y.is_column()) return bilin(A, x, y.T());
    if (!y.is_dense()) return bilin(A, x, densify(y));

    // Assert dimensions
    casadi_assert(x.size1()==A.size1() && y.size1()==A.size2(),
      "Dimension mismatch. Got x.size1() = " + str(x.size1())
      + " and y.size1() = " + str(y.size1()) + " but A.size() = " + str(A.size()));

    // Call the class specific method
    return MatType::_bilin(A, x, y);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::rank1(const MatType& A, const MatType& alpha,
                                        const MatType& x, const MatType& y) {
    // Check/correct x
    casadi_assert_dev(x.is_vector());
    if (!x.is_column()) return rank1(A, alpha, x.T(), y);
    if (!x.is_dense()) return rank1(A, alpha, densify(x), y);

    // Check/correct y
    casadi_assert_dev(y.is_vector());
    if (!y.is_column()) return rank1(A, alpha, x, y.T());
    if (!y.is_dense()) return rank1(A, alpha, x, densify(y));

    // Check alpha, quick return
    casadi_assert_dev(alpha.is_scalar());
    if (!alpha.is_dense()) return A;

    // Assert dimensions
    casadi_assert(x.size1()==A.size1() && y.size1()==A.size2(),
      "Dimension mismatch. Got x.size1() = " + str(x.size1())
        + " and y.size1() = " + str(y.size1())
        + " but A.size() = " + str(A.size()));

    // Call the class specific method
    return MatType::_rank1(A, alpha, x, y);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::logsumexp(const MatType& x) {
    casadi_assert(x.is_dense(), "Argument must be dense");
    casadi_assert(x.is_column(), "Argument must be column vector");
    casadi_assert(!x.is_empty(), "Argument must be nonempty");
    // Call the class specific method
    return MatType::_logsumexp(x);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::jtimes(const MatType &ex, const MatType &arg,
                                         const MatType &v, bool tr, const Dict& opts) {
    try {
      // Assert consistent input dimensions
      if (tr) {
        if (ex.size2()==0 && v.size2()>0) {
          casadi_error("Ambiguous dimensions.");
        }
        casadi_assert(v.size1() == ex.size1() &&
                      (v.size2()==0 || ex.size2()==0 || v.size2() % ex.size2() == 0),
                      "'v' has inconsistent dimensions: "
                      " v " + v.dim(false) + ", ex " + ex.dim(false) + ".");
      } else {
        if (arg.size2()==0 && v.size2()>0) {
          casadi_error("Ambiguous dimensions.");
        }
        casadi_assert(v.size1() == arg.size1() &&
                      (v.size2()==0 || arg.size2()==0 || v.size2() % arg.size2() == 0),
                      "'v' has inconsistent dimensions: "
                      " v " + v.dim(false) + ", arg " + arg.dim(false) + ".");
      }

      casadi_int n_seeds = 1;
      if (tr) {
        if (ex.size2()>0) n_seeds = v.size2() / ex.size2();
      } else {
        if (arg.size2()>0) n_seeds = v.size2() / arg.size2();
      }

      // Quick return if no seeds
      if (v.is_empty() || ex.is_empty()) {
        return MatType(tr ? arg.size1() : ex.size1(),
                       tr ? arg.size2()*n_seeds : ex.size2()*n_seeds);
      }

      // Split up the seed into its components
      std::vector<MatType> w = horzsplit(v, tr ? ex.size2() : arg.size2());

      // Seeds as a vector of vectors
      std::vector<std::vector<MatType> > ww(w.size());
      for (casadi_int i=0; i<w.size(); ++i) ww[i] = {w[i]};

      // Calculate directional derivatives
      if (tr) {
        ww = reverse({ex}, {arg}, ww, opts);
      } else {
        ww = forward({ex}, {arg}, ww, opts);
      }

      // Get results
      for (casadi_int i=0; i<w.size(); ++i) w[i] = ww[i][0];
      return horzcat(w);
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("jtimes", e.what());
    }
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::gradient(const MatType &ex, const MatType &arg,
      const Dict& opts) {
    try {
      casadi_assert(ex.is_scalar(),
                    "'gradient' only defined for scalar outputs: Use 'jacobian' instead.");
      return project(jtimes(ex, arg, MatType::ones(ex.sparsity()), true, opts), arg.sparsity());
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("gradient", e.what());
    }
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::tangent(const MatType &ex, const MatType &arg,
      const Dict& opts) {
    try {
      casadi_assert(arg.is_scalar(),
                    "'tangent' only defined for scalar inputs: Use 'jacobian' instead.");
      return project(jtimes(ex, arg, MatType::ones(arg.sparsity()), false, opts), ex.sparsity());
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("tangent", e.what());
    }
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::
  linearize(const MatType& f, const MatType& x, const MatType& x0, const Dict& opts) {
    MatType x_lin = MatType::sym("x_lin", x.sparsity());
    // mismatching dimensions
    if (x0.size() != x.size()) {
      // Scalar x0 is ok
      if (x0.sparsity().is_scalar()) {
        return linearize(f, x, MatType(x.sparsity(), x0));
      }
      casadi_error("Dimension mismatch in 'linearize'");
    }
    return substitute(f + jtimes(f, x, x_lin, false, opts),
      MatType::vertcat({x_lin, x}), MatType::vertcat({x, x0}));
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::mpower(const MatType& a,
                                            const MatType& b) {
    if (a.is_scalar() && b.is_scalar()) return pow(a, b);
    casadi_assert(a.is_square() && b.is_constant() && b.is_scalar(),
      "Not Implemented");
    double bv = static_cast<double>(b);
    casadi_int N = static_cast<casadi_int>(bv);
    casadi_assert(bv-static_cast<double>(N)==0, "mpower only defined for integer powers.");
    casadi_assert(bv==N, "Not Implemented");
    if (N<0) return inv(mpower(a, -N));
    if (N==0) return MatType::eye(a.size1());
    if (N==1) return a;
    if (N % 2 == 0) {
      MatType h = mpower(a, N/2); // NOLINT
      return MatType::mtimes(h, h);
    } else {
      return MatType::mtimes(mpower(a, N-1), a);
    }
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::soc(const MatType& x,
                                            const MatType& y) {
    casadi_assert(y.is_scalar(), "y needs to be scalar. Got " + y.dim() + ".");
    casadi_assert(x.is_vector(), "x needs to be a vector. Got " + x.dim() + ".");

    MatType x_col = x.is_column() ? x : x.T();

    x_col = x_col.nz(Slice()); // NOLINT

    casadi_int n = x_col.numel();
    return blockcat(y*MatType::eye(n), x_col, x_col.T(), y);
  }

  template<typename MatType>
  bool GenericMatrix<MatType>::is_linear(const MatType &expr, const MatType &var) {
    return !any(MatType::which_depends(expr, var, 2, true));
  }

  template<typename MatType>
  bool GenericMatrix<MatType>::is_quadratic(const MatType &expr, const MatType &var) {
    return is_linear(gradient(expr, var), var);
  }

  template<typename MatType>
  void GenericMatrix<MatType>::quadratic_coeff(const MatType &expr, const MatType &var,
          MatType& A, MatType& b, MatType& c, bool check) {
    casadi_assert(expr.is_scalar(), "'quadratic_coeff' only defined for scalar expressions.");
    A = hessian(expr, var);
    b = substitute(jacobian(expr, var), var, 0).T();
    if (check)
      casadi_assert(!depends_on(A, var), "'quadratic_coeff' called on non-quadratic expression.");
    c = substitute(expr, var, 0);
  }

  template<typename MatType>
  void GenericMatrix<MatType>::linear_coeff(const MatType &expr, const MatType &var,
          MatType& A, MatType& b, bool check) {
    casadi_assert(expr.is_vector(), "'linear_coeff' only defined for vector expressions.");
    if (check)
      casadi_assert(is_linear(expr, var), "'linear_coeff' called on non-linear expression.");
    A = substitute(jacobian(expr, var), var, 0);
    b = vec(substitute(expr, var, 0));
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::diff(const MatType& x, casadi_int n, casadi_int axis) {
    casadi_assert(axis==-1 || axis==0 || axis==1, "Axis argument invalid");
    casadi_assert(n>=1, "n argument invalid");

    MatType ret = x;
    for (casadi_int i=0;i<n;++i) {
      // Matlab's special case
      if (axis==-1 && ret.is_scalar()) return MatType();

      casadi_int local_axis = (axis==-1) ? ret.is_row() : axis;
      if (local_axis==0) {
        if (ret.size1()<=1) {
          ret = MatType::zeros(0, ret.size2());
        } else {
          ret = ret(Slice(1, ret.size1()), Slice())-ret(Slice(0, ret.size1()-1), Slice());
        }
      } else {
        if (ret.size2()<=1) {
          ret = MatType::zeros(ret.size1(), 0);
        } else {
          ret = ret(Slice(), Slice(1, ret.size2()))-ret(Slice(), Slice(0, ret.size2()-1));
        }
      }
    }
    return ret;
  }

#undef CASADI_THROW_ERROR

} // namespace casadi

#endif // CASADI_GENERIC_MATRIX_HPP
