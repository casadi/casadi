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


#ifndef CASADI_MATRIX_DECL_HPP
#define CASADI_MATRIX_DECL_HPP

#include "matrix_fwd.hpp"
#include "exception.hpp"
#include "casadi_limits.hpp"
#include "runtime/casadi_runtime.hpp"
#include "generic_matrix.hpp"
#include "generic_expression.hpp"
#include "printable.hpp"

#include <random>
#include <typeinfo>
#include <vector>
#include <initializer_list>
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

namespace casadi {

  /** \brief Empty Base

      This class is extended in SWIG.

      \identifier{18c} */
  struct CASADI_EXPORT MatrixCommon {};

/// \cond CLUTTER
  ///@{
  /** \brief Get typename

      \identifier{18d} */
  template <typename Scalar> inline std::string matrixName()
  { return std::string("Matrix<") + typeid(Scalar).name() + std::string(">");}
  template<> inline std::string matrixName<double>() { return "DM"; }
  template<> inline std::string matrixName<casadi_int>() { return "IM"; }
  ///@}
/// \endcond

  /** \brief Sparse matrix class. SX and DM are specializations.

      General sparse matrix class that is designed with the idea that "everything is a matrix",
      that is, also scalars and vectors.\n
      This philosophy makes it easy to use and to interface in particularly
      with Python and Matlab/Octave.\n

      Index starts with 0.\n
      Index vec happens as follows: (rr, cc) -> k = rr+cc*size1()\n
      Vectors are column vectors.\n

      The storage format is Compressed Column Storage (CCS), similar to that used for
      sparse matrices in Matlab, \n
      but unlike this format, we do allow for elements to be structurally non-zero
      but numerically zero.\n

      Matrix<Scalar> is polymorphic with a std::vector<Scalar> that contain
      all non-identical-zero elements.\n
      The sparsity can be accessed with Sparsity& sparsity()\n

      SX holds symbolic scalar entries (SXElem), DM holds doubles. A structural zero prints
      as 00, a numerical zero as 0. The structural nonzeros are available as a flat
      column-major vector with nonzeros().

      \lang{python}
      \doctest
      >>> A = DM([[1, 2], [3, 4]])
      >>> A.sparsity().dim(), A.nonzeros()
      ('2x2', [1.0, 3.0, 2.0, 4.0])
      >>> S = DM(Sparsity.diag(2))
      >>> S
      DM(
      [[1, 00],
       [00, 1]])
      >>> S.nnz(), S.numel()
      (2, 4)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> A = DM([1 2; 3 4]);
      >> A.sparsity().dim(), A.nonzeros()
      2x2
      1 3 2 4
      >> S = DM(Sparsity.diag(2));
      >> S
      [[1, 00],
       [00, 1]]
      >> fprintf('%d %d\n', S.nnz(), S.numel())
      2 4
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> A = ca.DM([1 2; 3 4]);
      julia> ca.dim(ca.sparsity(A)), ca.nonzeros(A)
      ("2x2", [1.0, 3.0, 2.0, 4.0])
      julia> S = ca.DM(ca.Sparsity.diag(2))
      DM(
      [[1, 00],
       [00, 1]])
      julia> ca.nnz(S), ca.numel(S)
      (2, 4)
      \enddoctest
      \endlang

      \author Joel Andersson
      \date 2010-2014

      \identifier{18e} */
  template<typename Scalar>
  class CASADI_EXPORT Matrix :
    public MatrixCommon,
    public SWIG_IF_ELSE(GenericExpressionCommon, GenericExpression<Matrix<Scalar> >),
    public GenericMatrix<Matrix<Scalar> >,
    public SWIG_IF_ELSE(PrintableCommon, Printable<Matrix<Scalar> >) {
  public:

    /** \brief Empty 0-by-0 matrix

        \lang{python}
        \doctest
        >>> DM().dim(), SX().is_empty()
        ('0x0', True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%s %d\n', DM().dim(), SX().is_empty())
        0x0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.dim(ca.DM()), ca.is_empty(ca.SX())
        ("0x0", true)
        \enddoctest
        \endlang

        \identifier{18f} */
    Matrix();

    /// Copy constructor: the copy has its own sparsity pattern and nonzeros
    Matrix(const Matrix<Scalar>& m);

#ifndef SWIG
    /// Assignment (normal)
    Matrix<Scalar>& operator=(const Matrix<Scalar>& m);
#endif // SWIG

    /** \brief Create a matrix of the given shape with all structural zeros (no nonzeros)

        The result has nnz()==0. Use zeros(nrow, ncol) for a dense matrix of numerical zeros.
        Structural zeros print as 00.

        \param nrow Number of rows, nonnegative
        \param ncol Number of columns, nonnegative

        \lang{python}
        \doctest
        >>> DM(2, 3)
        DM(
        [[00, 00, 00],
         [00, 00, 00]])
        >>> DM(2, 3).nnz(), DM.zeros(2, 3).nnz()
        (0, 6)
        >>> SX(2, 0).dim()
        '2x0'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(2, 3)
        [[00, 00, 00],
         [00, 00, 00]]
        >> fprintf('%d %d\n', DM(2, 3).nnz(), DM.zeros(2, 3).nnz())
        0 6
        >> SX(2, 0).dim()
        2x0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(2, 3)
        DM(
        [[00, 00, 00],
         [00, 00, 00]])
        julia> ca.nnz(ca.DM(2, 3)), ca.nnz(ca.DM.zeros(2, 3))
        (0, 6)
        julia> ca.dim(ca.SX(2, 0))
        "2x0"
        \enddoctest
        \endlang

        \sa zeros, ones

        \identifier{18g} */
    Matrix(casadi_int nrow, casadi_int ncol);

#ifndef SWIG
    /** \brief Create a sparse matrix with all structural zeros

        \identifier{18h} */
    explicit Matrix(const std::pair<casadi_int, casadi_int>& rc);

    /** \brief  Access functions of the node

        \identifier{18i} */
    std::vector<Scalar>* operator->();

    /** \brief  Const access functions of the node

        \identifier{18j} */
    const std::vector<Scalar>* operator->() const;
#endif // SWIG

    /** \brief Create a matrix with the given sparsity pattern and all nonzeros equal to one

        Same as ones(sp).

        \lang{python}
        \doctest
        >>> DM(Sparsity.diag(2))
        DM(
        [[1, 00],
         [00, 1]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.diag(2))
        [[1, 00],
         [00, 1]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.diag(2))
        DM(
        [[1, 00],
         [00, 1]])
        \enddoctest
        \endlang

        \identifier{18k} */
    explicit Matrix(const Sparsity& sp);

    /** \brief Create a matrix with a given sparsity pattern and given nonzeros

        \param sp Sparsity pattern of the result
        \param d Either a scalar, which is used for every structural nonzero, or a dense
                 vector (row or column) with exactly sp.nnz() entries, in column-major order
                 of the pattern. Anything else raises an error.

        \lang{python}
        \doctest
        >>> DM(Sparsity.diag(2), [5, 6])
        DM(
        [[5, 00],
         [00, 6]])
        >>> DM(Sparsity.diag(2), 5)
        DM(
        [[5, 00],
         [00, 5]])
        >>> SX(Sparsity.lower(2), SX.sym("v", 3))
        SX(
        [[v_0, 00],
         [v_1, v_2]])
        >>> try: DM(Sparsity.diag(3), [1, 2])
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.diag(2), [5; 6])
        [[5, 00],
         [00, 6]]
        >> DM(Sparsity.diag(2), 5)
        [[5, 00],
         [00, 5]]
        >> SX(Sparsity.lower(2), SX.sym('v', 3))
        [[v_0, 00],
         [v_1, v_2]]
        >> try, DM(Sparsity.diag(3), [1; 2]); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.diag(2), [5, 6])
        DM(
        [[5, 00],
         [00, 6]])
        julia> ca.DM(ca.Sparsity.diag(2), 5)
        DM(
        [[5, 00],
         [00, 5]])
        julia> ca.SX(ca.Sparsity.lower(2), ca.SX.sym("v", 3))
        SX(
        [[v_0, 00],
         [v_1, v_2]])
        julia> try ca.DM(ca.Sparsity.diag(3), [1, 2]); catch; println("error") end
        error
        \enddoctest
        \endlang

        \identifier{18l} */
    Matrix(const Sparsity& sp, const Matrix<Scalar>& d);

    /// Dense 1-by-1 matrix holding a number; enables implicit conversion from a numeric type
    Matrix(double val);

#if !(defined(SWIG) && defined(SWIGMATLAB))
    /** \brief Dense matrix from a list of rows

        Every row must have the same length; a ragged input raises an error.
        An empty outer list gives a 0-by-1 matrix.
        At the Python level, 2-D NumPy arrays and SciPy sparse matrices are accepted as well;
        a sparse input keeps its sparsity pattern.

        \lang{python}
        \doctest
        >>> DM([[1, 2], [3, 4]])
        DM(
        [[1, 2],
         [3, 4]])
        >>> DM(np.array([[1, 2], [3, 4]])).dim()
        '2x2'
        >>> try: DM([[1, 2], [3]])
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM([1 2; 3 4])
        DM(
        [[1, 2],
         [3, 4]])
        julia> try ca.DM([[1, 2], [3]]); catch; println("error") end
        error
        \enddoctest
        \endlang

        \identifier{2oo} */
    explicit Matrix(const std::vector< std::vector<double> >& m);

    /** \brief Dense column vector from a flat list of numbers

        A 1-D NumPy array also gives a column vector.

        \lang{python}
        \doctest
        >>> DM([1, 2, 3]).dim(), DM(np.array([1, 2, 3])).dim(), DM([]).dim()
        ('3x1', '3x1', '0x1')
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.dim(ca.DM([1, 2, 3])), ca.dim(ca.DM(Float64[]))
        ("3x1", "0x1")
        \enddoctest
        \endlang

        \identifier{18m} */
    template<typename A>
    Matrix(const std::vector<A>& x) : sparsity_(Sparsity::dense(x.size(), 1)),
      nonzeros_(std::vector<Scalar>(x.size())) {
        auto x_it = x.begin();
        for (auto&& d : nonzeros_) d = static_cast<Scalar>(*x_it++);
    }
#endif

    /** \brief Create a matrix from another matrix with a different entry type

        The sparsity pattern is kept and every nonzero is converted.
        SX(DM) always works. DM(SX) requires constant entries: a symbolic entry is
        converted to nan without an error. Use evalf to get an error instead.

        \lang{python}
        \doctest
        >>> SX(DM(Sparsity.diag(2))).dim()
        '2x2'
        >>> DM(SX([1, 2]))
        DM([1, 2])
        >>> DM(SX.sym("x"))
        DM(nan)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> SX(DM(Sparsity.diag(2))).dim()
        2x2
        >> DM(SX([1; 2]))
        [1, 2]
        >> DM(SX.sym('x'))
        nan
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.dim(ca.SX(ca.DM(ca.Sparsity.diag(2))))
        "2x2"
        julia> ca.DM(ca.SX([1, 2]))
        DM([1, 2])
        julia> ca.DM(ca.SX.sym("x"))
        DM(nan)
        \enddoctest
        \endlang

        \sa evalf

        \identifier{18n} */
    template<typename A>
    Matrix(const Matrix<A>& x) : sparsity_(x.sparsity()), nonzeros_(std::vector<Scalar>(x.nnz())) {
      auto x_it = x->begin();
      for (auto&& d : nonzeros_) d = static_cast<Scalar>(*x_it++);
    }

#ifndef SWIG
    /// Construct from a vector
    Matrix(const std::vector<Scalar>& x);

    /// Construct from initializer list
    Matrix(std::initializer_list<Scalar> x) : Matrix<Scalar>(std::vector<Scalar>(x)) {}

    /// Convert to scalar type
    const Scalar scalar() const;

    /// Scalar type
    typedef Scalar ScalarType;

    /// Base class
    typedef GenericMatrix<Matrix<Scalar> > B;

    /// Expose base class functions
    using B::nnz;
    using B::nnz_lower;
    using B::nnz_upper;
    using B::numel;
    using B::size1;
    using B::size2;
    using B::size;
    using B::is_empty;
    using B::is_scalar;
    using B::is_dense;
    using B::is_vector;
    using B::is_row;
    using B::is_column;
    using B::is_tril;
    using B::is_triu;
    using B::colind;
    using B::row;
    using B::dim;
    using B::sym;
    using B::zeros;
    using B::ones;
    using B::nz;
    using B::operator();
    using B::horzsplit;
    using B::vertsplit;
    using B::diagsplit;
    using B::mtimes;
#endif // SWIG

    /** \brief Check whether element (rr, cc) is a structural nonzero

        Looks only at the sparsity pattern; a stored numerical zero still counts as a nonzero.
        Negative indices count from the end. An index past the end raises an error.

        \lang{python}
        \doctest
        >>> S = DM(Sparsity.diag(2))
        >>> S.has_nz(0, 0), S.has_nz(0, 1), DM([0, 1]).has_nz(0, 0)
        (True, False, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> S = DM(Sparsity.diag(2));
        >> fprintf('%d %d %d\n', S.has_nz(0, 0), S.has_nz(0, 1), DM([0; 1]).has_nz(0, 0))
        1 0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> S = ca.DM(ca.Sparsity.diag(2));
        julia> ca.has_nz(S, 0, 0), ca.has_nz(S, 0, 1), ca.has_nz(ca.DM([0, 1]), 0, 0)
        (true, false, true)
        \enddoctest
        \endlang

        \sa has_zeros, Sparsity::has_nz

        \identifier{2op} */
    bool has_nz(casadi_int rr, casadi_int cc) const;

    /** \brief Truth value of a 1-by-1 matrix (Python: bool(x))

        Only defined for 1-by-1 matrices; other shapes raise an error.
        A structurally empty 1-by-1 matrix is false. A symbolic SX entry raises an error
        unless the on-the-fly simplifications already decided its value (x==x is true).

        \lang{python}
        \doctest
        >>> bool(DM(2)), bool(DM(0)), bool(DM(1, 1))
        (True, False, False)
        >>> x = SX.sym("x")
        >>> bool(x == x)
        True
        >>> try: bool(x > 0)
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> Bool(ca.DM(1)), Bool(ca.DM(0)), Bool(ca.DM(1, 1))
        (true, false, false)
        julia> x = ca.SX.sym("x");
        julia> Bool(x == x)
        true
        julia> try Bool(x > 0); catch; println("error") end
        error
        \enddoctest
        \endlang

        \identifier{2oq} */
    bool __nonzero__() const;

    ///@{
    /** \brief Get elements by a single (linear) index

        Backs the Python expressions x[k] and x[sp]. The index addresses elements, not
        nonzeros, in column-major order: element (i, j) has linear index i + j*size1().
        Negative indices count from the end; out-of-range indices raise an error.

        - A scalar index returns a 1-by-1 matrix. If that element is a structural zero,
          the result is structurally empty (prints 00).
        - A Slice (Python slice) or a list of indices returns a column vector,
          except that indexing a row vector gives a row vector.
        - A matrix of indices (IM) returns a matrix of the same shape.
        - A Sparsity pattern must have the same shape as the matrix; the result is the
          projection of the matrix onto that pattern (entries outside it are dropped).

        Note that the result of a comparison such as A > 2 is a DM of zeros and ones, not a
        boolean mask: indexing with it uses those values as linear indices. Use a Sparsity
        pattern or the nz accessor for masking.

        \param ind1 True for one-based indices (MATLAB); Python always passes false
        \param rr Linear index, Slice or list/matrix of indices
        \param sp Sparsity pattern of the same shape as the matrix

        \lang{python}
        \doctest
        >>> A = DM([[1, 2, 3], [4, 5, 6]])
        >>> A[1], A[2], A[-1]
        (DM(4), DM(2), DM(6))
        >>> A[[0, 5]]
        DM([1, 6])
        >>> A[Sparsity.diag(2, 3)]
        DM(
        [[1, 00, 00],
         [00, 5, 00]])
        >>> S = DM(Sparsity.diag(2))
        >>> S[1]
        DM(00)
        >>> A.get(False, Slice(0, 2))
        DM([1, 4])
        >>> try: A[6]
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2 3; 4 5 6]);
        >> A(2), A(3), A(end)
        4
        2
        6
        >> A([1 6])
        [[1, 6]]
        >> A(Sparsity.diag(2, 3))
        [[1, 00, 00],
         [00, 5, 00]]
        >> S = DM(Sparsity.diag(2));
        >> S(2)
        00
        >> A.get(true, [1; 2])
        [1, 4]
        >> try, A(7); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2 3; 4 5 6]);
        julia> A[2], A[3], A[end]
        (4, 2, 6)
        julia> A[[1, 6]]
        DM([1, 6])
        julia> S = ca.DM(ca.Sparsity.diag(2));
        julia> S[2]
        DM(00)
        julia> try A[7]; catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa set, get_nz, Sparsity::sub

        \identifier{2or} */
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Slice& rr) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Sparsity& sp) const;
    ///@}

    ///@{
    /** \brief Get a submatrix by row and column indices

        Backs the Python expression x[i, j]. Each of \p rr and \p cc is a scalar, a Slice or
        a list/vector of indices (any order, repetitions allowed, negative counts from the end).
        The result has len(rr) rows and len(cc) columns and keeps the sparsity of the
        selected entries. Two scalar indices give a 1-by-1 matrix (structurally empty when the
        entry is a structural zero). Out-of-range indices raise an error.

        \param ind1 True for one-based indices (MATLAB); Python always passes false
        \param rr Row index, Slice or list/matrix of row indices
        \param cc Column index, Slice or list/matrix of column indices

        \lang{python}
        \doctest
        >>> A = DM([[1, 2, 3], [4, 5, 6]])
        >>> A[0, 1], A[1, -1]
        (DM(2), DM(6))
        >>> A[:, 1]
        DM([2, 5])
        >>> A[[1, 0], [0, 2]]
        DM(
        [[4, 6],
         [1, 3]])
        >>> A[1:1, :].dim()
        '0x3'
        >>> x = SX.sym("x", 2, 2)
        >>> x[0, :]
        SX([[x_0, x_2]])
        >>> DM(Sparsity.diag(2))[1, 0]
        DM(00)
        >>> try: A[2, 0]
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2 3; 4 5 6]);
        >> A(1, 2), A(2, end)
        2
        6
        >> A(:, 2)
        [2, 5]
        >> A([2 1], [1 3])
        [[4, 6],
         [1, 3]]
        >> B = A(2:1, :); disp(B.dim())
        0x3
        >> x = SX.sym('x', 2, 2);
        >> x(1, :)
        [[x_0, x_2]]
        >> S = DM(Sparsity.diag(2));
        >> S(2, 1)
        00
        >> try, A(3, 1); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2 3; 4 5 6]);
        julia> A[1, 2], A[2, end]
        (2, 6)
        julia> A[:, 2]
        DM([2, 5])
        julia> A[[2, 1], [1, 3]]
        DM(
        [[4, 6],
         [1, 3]])
        julia> ca.dim(A[2:1, :])
        "0x3"
        julia> x = ca.SX.sym("x", 2, 2);
        julia> x[1, :]
        SX([[x_0, x_2]])
        julia> ca.DM(ca.Sparsity.diag(2))[2, 1]
        DM(00)
        julia> try A[3, 1]; catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa set

        \identifier{2os} */
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Slice& rr, const Slice& cc) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Slice& rr, const Matrix<casadi_int>& cc) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Matrix<casadi_int>& rr, const Slice& cc) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) const;
    ///@}

    ///@{
    /** \brief Assign elements by a single (linear) index

        Backs the Python statements x[k] = v and x[sp] = v. Indices are linear, column-major
        and zero-based (see get). The matrix is modified in place.

        - \p m must be a scalar (assigned to every selected element) or have as many
          entries as there are indices (a vector of the other orientation is transposed).
          Any other shape raises an error.
        - Assigning to a structural zero inserts a new nonzero; the sparsity pattern grows.
          Assigning a structurally empty value (e.g. DM(1, 1)) removes the entry instead.
        - With a Sparsity index of the same shape, the entries at the pattern's positions
          receive \p m (a scalar) or the corresponding entries of \p m (a matrix).

        \param m Value(s) to assign
        \param ind1 True for one-based indices (MATLAB); Python always passes false
        \param rr Linear index, Slice or list/matrix of indices
        \param sp Sparsity pattern of the same shape as the matrix

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> A[1] = 10; A[[0, 3]] = [7, 8]
        >>> A
        DM(
        [[7, 2],
         [10, 8]])
        >>> A[Sparsity.diag(2)] = 0
        >>> A
        DM(
        [[0, 2],
         [10, 0]])
        >>> S = DM(Sparsity.diag(2))
        >>> S[1] = 5
        >>> S.nnz(), S[0, 1]
        (3, DM(00))
        >>> try: A[[0, 1]] = [1, 2, 3]
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> A(2) = 10; A([1 4]) = [7 8];
        >> A
        [[7, 2],
         [10, 8]]
        >> A(Sparsity.diag(2)) = 0;
        >> A
        [[0, 2],
         [10, 0]]
        >> S = DM(Sparsity.diag(2));
        >> S(2) = 5;
        >> S.nnz(), S(1, 2)
        3
        00
        >> try, A([1 2]) = [1 2 3]; catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> A[2] = 10; A[[1, 4]] = [7, 8];
        julia> A
        DM(
        [[7, 2],
         [10, 8]])
        julia> S = ca.DM(ca.Sparsity.diag(2));
        julia> S[2] = 5;
        julia> ca.nnz(S), S[1, 2]
        (3, 00)
        julia> try A[[1, 2]] = [1, 2, 3]; catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa get, set_nz, erase

        \identifier{2ot} */
    void set(const Matrix<Scalar>& m, bool ind1, const Slice& rr);
    void set(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr);
    void set(const Matrix<Scalar>& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /** \brief Assign a submatrix by row and column indices

        Backs the Python statement x[i, j] = v. The matrix is modified in place.
        \p m must be len(rr)-by-len(cc), a scalar (broadcast to all selected entries) or a
        vector of the transposed orientation (transposed automatically); other shapes raise
        an error. Assigning into structural zeros inserts nonzeros; assigning a sparse \p m
        first erases the selected block, so its structural zeros become structural zeros.

        \param m Value(s) to assign
        \param ind1 True for one-based indices (MATLAB); Python always passes false
        \param rr Row index, Slice or list/matrix of row indices
        \param cc Column index, Slice or list/matrix of column indices

        \lang{python}
        \doctest
        >>> A = DM([[1, 2, 3], [4, 5, 6]])
        >>> A[0, 1] = 9; A[1, :] = 0
        >>> A
        DM(
        [[1, 9, 3],
         [0, 0, 0]])
        >>> A[0, :] = DM([7, 8, 9])
        >>> A[0, :]
        DM([[7, 8, 9]])
        >>> A[:, 0] = DM(2, 1)
        >>> A.nnz()
        4
        >>> x = SX.sym("x", 2, 2)
        >>> x[0, 0] = 0
        >>> x
        SX(
        [[0, x_2],
         [x_1, x_3]])
        >>> try: A[0, :] = [1, 2]
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2 3; 4 5 6]);
        >> A(1, 2) = 9; A(2, :) = 0;
        >> A
        [[1, 9, 3],
         [0, 0, 0]]
        >> A(1, :) = DM([7; 8; 9]);
        >> A(1, :)
        [[7, 8, 9]]
        >> A(:, 1) = DM(2, 1);
        >> A.nnz()
        4
        >> x = SX.sym('x', 2, 2);
        >> x(1, 1) = 0;
        >> x
        [[0, x_2],
         [x_1, x_3]]
        >> try, A(1, :) = [1 2]; catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2 3; 4 5 6]);
        julia> A[1, 2] = 9; A[2, :] = 0;
        julia> A
        DM(
        [[1, 9, 3],
         [0, 0, 0]])
        julia> A[1, :] = ca.DM([7, 8, 9]);
        julia> A[1, :]
        DM([[7, 8, 9]])
        julia> A[:, 1] = ca.DM(2, 1);
        julia> ca.nnz(A)
        4
        julia> x = ca.SX.sym("x", 2, 2);
        julia> x[1, 1] = 0;
        julia> x
        SX(
        [[0, x_2],
         [x_1, x_3]])
        julia> try A[1, :] = [1, 2]; catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa get

        \identifier{2ou} */
    void set(const Matrix<Scalar>& m, bool ind1, const Slice& rr, const Slice& cc);
    void set(const Matrix<Scalar>& m, bool ind1, const Slice& rr, const Matrix<casadi_int>& cc);
    void set(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr, const Slice& cc);
    void set(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr,
                                                 const Matrix<casadi_int>& cc);
    ///@}

    ///@{
    /** \brief Get structural nonzeros by their position in the nonzero storage

        Backs the Python expression x.nz[k]. Index k runs over the nonzeros in column-major
        order, 0 <= k < nnz(); negative values count from the end and out-of-range values
        raise an error. Structural zeros are skipped, so for a dense matrix this coincides
        with linear indexing. The result is a column vector shaped like the index
        (a row vector when indexing a row vector; a matrix for a matrix index).
        Iterating over x.nz in Python yields the nonzeros one by one.

        \param ind1 True for one-based indices (MATLAB); Python always passes false
        \param k Nonzero index, Slice or list/matrix of indices

        \lang{python}
        \doctest
        >>> S = DM(Sparsity.diag(3)); S[0, 2] = 7
        >>> S.nz[:]
        DM([1, 1, 7, 1])
        >>> S.nz[2], S.nz[-1], len(S.nz)
        (DM(7), DM(1), 4)
        >>> x = SX.sym("x", Sparsity.lower(2))
        >>> x.nz[[0, 2]]
        SX([x_0, x_2])
        >>> SX.sym("e", 2, 0).nz[:].dim()
        '0x1'
        >>> try: S.nz[4]
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> S = DM(Sparsity.diag(3)); S(1, 3) = 7;
        >> S.get_nz(true, Slice())
        [1, 1, 7, 1]
        >> S.get_nz(true, 3), S.get_nz(true, [1; 3])
        7
        [1, 7]
        >> x = SX.sym('x', Sparsity.lower(2));
        >> x.get_nz(true, [1; 3])
        [x_0, x_2]
        >> SX.sym('e', 2, 0).get_nz(true, Slice()).dim()
        0x1
        >> try, S.get_nz(true, 5); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> S = ca.DM(ca.Sparsity.diag(3)); S[1, 3] = 7;
        julia> S.nz[:]
        DM([1, 1, 7, 1])
        julia> S.nz[3], S.nz[end], length(S.nz)
        (7, 1, 4)
        julia> x = ca.SX.sym("x", ca.Sparsity.lower(2));
        julia> x.nz[[1, 3]]
        SX([x_0, x_2])
        julia> ca.dim(ca.SX.sym("e", 2, 0).nz[:])
        "0x1"
        julia> try S.nz[5]; catch e; println(typeof(e)) end
        BoundsError
        \enddoctest
        \endlang

        \sa set_nz, get, nonzeros

        \identifier{2ov} */
    void get_nz(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Slice& k) const;
    void get_nz(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& k) const;
    ///@}

    ///@{
    /** \brief Assign structural nonzeros by their position in the nonzero storage

        Backs the Python statement x.nz[k] = v. The sparsity pattern never changes: only
        stored values are overwritten (writing 0 keeps the entry as a structural nonzero).
        \p m must be a scalar (assigned to every selected nonzero) or have as many entries
        as there are indices; a mismatch raises an error.

        \param m Value(s) to assign
        \param ind1 True for one-based indices (MATLAB); Python always passes false
        \param k Nonzero index, Slice or list/matrix of indices

        \lang{python}
        \doctest
        >>> S = DM(Sparsity.diag(3))
        >>> S.nz[:] = [4, 5, 6]; S.nz[0] = 0
        >>> S
        DM(
        [[0, 00, 00],
         [00, 5, 00],
         [00, 00, 6]])
        >>> S.nnz()
        3
        >>> try: S.nz[[0, 1]] = [1, 2, 3]
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> S = DM(Sparsity.diag(3));
        >> S.set_nz([4 5 6], true, 1:3); S.set_nz(0, true, 1);
        >> S
        [[0, 00, 00],
         [00, 5, 00],
         [00, 00, 6]]
        >> S.nnz()
        3
        >> try, S.set_nz([1 2 3], true, [1 2]); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> S = ca.DM(ca.Sparsity.diag(3));
        julia> S.nz[:] = [4, 5, 6]; S.nz[1] = 0;
        julia> S
        DM(
        [[0, 00, 00],
         [00, 5, 00],
         [00, 00, 6]])
        julia> ca.nnz(S)
        3
        julia> try S.nz[[1, 2]] = [1, 2, 3]; catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa get_nz, set

        \identifier{2ow} */
    void set_nz(const Matrix<Scalar>& m, bool ind1, const Slice& k);
    void set_nz(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& k);
    ///@}

    Matrix<Scalar> operator+() const;
    Matrix<Scalar> operator-() const;

    /// \cond INTERNAL
    ///@{
    /** \brief  Create nodes by their ID

        \identifier{18o} */
    static Matrix<Scalar> binary(casadi_int op, const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> unary(casadi_int op, const Matrix<Scalar> &x);
    static Matrix<Scalar> scalar_matrix(casadi_int op,
                                          const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> matrix_scalar(casadi_int op,
                                          const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> matrix_matrix(casadi_int op,
                                          const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static std::vector< Matrix<Scalar> > call(const Function& f,
        const std::vector< Matrix<Scalar> > &x);
    static std::vector< Scalar > call(const Function& f, const std::vector< Scalar > &x);
    ///@}
    /// \endcond

#ifndef SWIG
    /// \cond CLUTTER
    ///@{
    /// Functions called by friend functions defined for GenericExpression
    static bool is_equal(const Matrix<Scalar> &x, const Matrix<Scalar> &y, casadi_int depth=0);
    static Matrix<Scalar> mmin(const Matrix<Scalar> &x);
    static Matrix<Scalar> mmax(const Matrix<Scalar> &x);
    ///@}

    ///@{
    /// Functions called by friend functions defined for GenericMatrix
    static Matrix<Scalar> simplify(const Matrix<Scalar> &x);
    static Matrix<Scalar> transform(const Matrix<Scalar> &x, const Dict& opts = Dict());
    static Matrix<Scalar> transform(const Matrix<Scalar> &x,
        const std::vector<std::vector<GenericType> >& passes, const Dict& opts = Dict());
    static std::vector<Matrix<Scalar> > transform(const std::vector<Matrix<Scalar> >& x,
        const Dict& opts = Dict());
    static std::vector<Matrix<Scalar> > transform(const std::vector<Matrix<Scalar> >& x,
        const std::vector<std::vector<GenericType> >& passes, const Dict& opts = Dict());
    static Matrix<Scalar> jacobian(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                                   const Dict& opts = Dict());
    static Sparsity jacobian_sparsity(const Matrix<Scalar> &f, const Matrix<Scalar> &x);
    static Matrix<Scalar> hessian(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
      const Dict& opts = Dict());
    static Matrix<Scalar> hessian(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                                    Matrix<Scalar>& g, const Dict& opts = Dict());
    static Matrix<Scalar>
      substitute(const Matrix<Scalar>& ex,
                 const Matrix<Scalar>& v,
                 const Matrix<Scalar>& vdef);
    static std::vector<Matrix<Scalar> >
      substitute(const std::vector<Matrix<Scalar> >& ex,
                 const std::vector<Matrix<Scalar> >& v,
                 const std::vector<Matrix<Scalar> >& vdef);
    static void substitute_inplace(const std::vector<Matrix<Scalar> >& v,
                                  std::vector<Matrix<Scalar> >& vdef,
                                  std::vector<Matrix<Scalar> >& ex,
                                  bool revers);
    static Matrix<Scalar> pinv(const Matrix<Scalar> &x);
    static Matrix<Scalar> pinv(const Matrix<Scalar> &A,
                                 const std::string& lsolver, const Dict& opts);
    static Matrix<Scalar> expm_const(const Matrix<Scalar> &A, const Matrix<Scalar> &t);
    static Matrix<Scalar> expm(const Matrix<Scalar> &A);
    static Matrix<Scalar> solve(const Matrix<Scalar> &A, const Matrix<Scalar>& b);
    static Matrix<Scalar> solve(const Matrix<Scalar> &A, const Matrix<Scalar>& b,
                                  const std::string& lsolver, const Dict& opts);
    static Matrix<Scalar> inv(const Matrix<Scalar> &A);
    static Matrix<Scalar> inv(const Matrix<Scalar> &A,
                                  const std::string& lsolver, const Dict& opts);

    static casadi_int n_nodes(const Matrix<Scalar> &x);
    static std::string print_operator(const Matrix<Scalar> &x,
                                      const std::vector<std::string>& args);
    static void extract(std::vector<Matrix<Scalar>>& ex, std::vector<Matrix<Scalar>>& v,
      std::vector<Matrix<Scalar>>& vdef, const Dict& opts = Dict());
    static void shared(std::vector<Matrix<Scalar> >& ex,
                              std::vector<Matrix<Scalar> >& v,
                              std::vector<Matrix<Scalar> >& vdef,
                              const std::string& v_prefix,
                              const std::string& v_suffix);
    static Matrix<Scalar> _bilin(const Matrix<Scalar>& A,
                                   const Matrix<Scalar>& x,
                                   const Matrix<Scalar>& y);
    static Matrix<Scalar> _rank1(const Matrix<Scalar>& A,
                                   const Matrix<Scalar>& alpha,
                                   const Matrix<Scalar>& x,
                                   const Matrix<Scalar>& y);
    static Matrix<Scalar> if_else(const Matrix<Scalar> &x,
                                    const Matrix<Scalar> &if_true,
                                    const Matrix<Scalar> &if_false,
                                    bool short_circuit=false);
    static Matrix<Scalar> conditional(const Matrix<Scalar> &ind,
                                        const std::vector<Matrix<Scalar> > &x,
                                        const Matrix<Scalar> &x_default,
                                        bool short_circuit=false);
    static bool depends_on(const Matrix<Scalar> &x, const Matrix<Scalar> &arg);
    static bool contains(const std::vector<Matrix<Scalar> >& v, const Matrix<Scalar> &n);
    static bool contains_all(const std::vector<Matrix<Scalar> >& v,
        const std::vector<Matrix<Scalar> > &n);
    static bool contains_any(const std::vector<Matrix<Scalar> >& v,
        const std::vector<Matrix<Scalar> > &n);

    static Matrix<Scalar> mrdivide(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> mldivide(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static std::vector<Matrix<Scalar> > symvar(const Matrix<Scalar> &x);
    static Matrix<Scalar> det(const Matrix<Scalar> &x);
    static Matrix<Scalar> det(const Matrix<Scalar> &x, const std::string& lsolver,
                              const Dict& dict=Dict());
    static Matrix<Scalar> inv_minor(const Matrix<Scalar> &x);
    static Matrix<Scalar> trace(const Matrix<Scalar> &x);
    static Matrix<Scalar> norm_1(const Matrix<Scalar> &x);
    static Matrix<Scalar> norm_2(const Matrix<Scalar> &x);
    static Matrix<Scalar> norm_fro(const Matrix<Scalar> &x);
    static Matrix<Scalar> norm_inf(const Matrix<Scalar> &x);
    static Matrix<Scalar> sum2(const Matrix<Scalar> &x);
    static Matrix<Scalar> sum1(const Matrix<Scalar> &x);
    static Matrix<Scalar> dot(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> nullspace(const Matrix<Scalar> &x);
    static Matrix<Scalar> diag(const Matrix<Scalar> &x);
    static Matrix<Scalar> unite(const Matrix<Scalar> &A, const Matrix<Scalar>& B);
    static Matrix<Scalar> project(const Matrix<Scalar> &x,
                                    const Sparsity& sp, bool intersect=false);
    static Matrix<Scalar> polyval(const Matrix<Scalar> &p, const Matrix<Scalar>& x);
    static Matrix<Scalar> densify(const Matrix<Scalar> &x, const Matrix<Scalar>& val);
    static Matrix<Scalar> densify(const Matrix<Scalar> &x);
    static Matrix<Scalar> einstein(const Matrix<Scalar>& A, const Matrix<Scalar>& B,
      const Matrix<Scalar>& C,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);

    static Matrix<Scalar> einstein(const Matrix<Scalar>& A, const Matrix<Scalar>& B,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);
    static Matrix<Scalar> cumsum(const Matrix<Scalar> &x, casadi_int axis=-1);
    static Matrix<Scalar> _logsumexp(const Matrix<Scalar>& x);
    static std::vector< Matrix<Scalar> > cse(const std::vector< Matrix<Scalar> >& e);
    ///@}

    ///@{
    /// Functions called by friend functions defined for SparsityInterface
    static Matrix<Scalar> blockcat(const std::vector< std::vector<Matrix<Scalar> > > &v);
    static Matrix<Scalar> horzcat(const std::vector<Matrix<Scalar> > &v);
    static std::vector<Matrix<Scalar> >
      horzsplit(const Matrix<Scalar>& x,
                const std::vector<casadi_int>& offset);
    static Matrix<Scalar> vertcat(const std::vector<Matrix<Scalar> > &v);
    static std::vector< Matrix<Scalar> >
      vertsplit(const Matrix<Scalar>& x,
                const std::vector<casadi_int>& offset);
    static std::vector< Matrix<Scalar> >
      diagsplit(const Matrix<Scalar>& x,
                const std::vector<casadi_int>& offset1,
                const std::vector<casadi_int>& offset2);
    static Matrix<Scalar> reshape(const Matrix<Scalar> &x, casadi_int nrow, casadi_int ncol);
    static Matrix<Scalar> reshape(const Matrix<Scalar> &x, const Sparsity& sp);
    static Matrix<Scalar> sparsity_cast(const Matrix<Scalar> &x, const Sparsity& sp);
    static Matrix<Scalar> kron(const Matrix<Scalar> &x, const Matrix<Scalar>& y);
    static Matrix<Scalar> mtimes(const Matrix<Scalar> &x, const Matrix<Scalar> &y,
                                 const std::string& blas = "reference");
    static Matrix<Scalar> mac(const Matrix<Scalar> &x,
                                const Matrix<Scalar> &y,
                                const Matrix<Scalar> &z,
                                const std::string& blas = "reference");
    static void extract_parametric(const Matrix<Scalar> &expr,
        const Matrix<Scalar>& par,
        Matrix<Scalar>& expr_ret,
        std::vector<Matrix<Scalar> >& symbols,
        std::vector<Matrix<Scalar>>& parametric,
        const Dict& opts);
    static void separate_linear(const Matrix<Scalar> &expr,
      const Matrix<Scalar> &sym_lin, const Matrix<Scalar> &sym_const,
      Matrix<Scalar>& expr_const, Matrix<Scalar>& expr_lin, Matrix<Scalar>& expr_nonlin);
    ///@}

    ///@{
    /// Functions called by friend functions defined here
    static Matrix<Scalar> sparsify(const Matrix<Scalar> &x, double tol=0);
    static void expand(const Matrix<Scalar>& x,
                       Matrix<Scalar>& weights,
                       Matrix<Scalar>& terms);
    static Matrix<Scalar> pw_const(const Matrix<Scalar> &t,
                                     const Matrix<Scalar> &tval, const Matrix<Scalar> &val);
    static Matrix<Scalar> pw_lin(const Matrix<Scalar> &t,
                                   const Matrix<Scalar> &tval, const Matrix<Scalar> &val);
    static Matrix<Scalar> heaviside(const Matrix<Scalar> &x);
    static Matrix<Scalar> rectangle(const Matrix<Scalar> &x);
    static Matrix<Scalar> triangle(const Matrix<Scalar> &x);
    static Matrix<Scalar> ramp(const Matrix<Scalar> &x);
    static Matrix<Scalar> gauss_quadrature(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x, const Matrix<Scalar> &a,
                                             const Matrix<Scalar> &b, casadi_int order=5);
    static Matrix<Scalar> gauss_quadrature(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x, const Matrix<Scalar> &a,
                                             const Matrix<Scalar> &b, casadi_int order,
                                             const Matrix<Scalar>& w);
    static std::vector<std::vector<Matrix<Scalar> > >
    forward(const std::vector<Matrix<Scalar> > &ex,
            const std::vector<Matrix<Scalar> > &arg,
            const std::vector<std::vector<Matrix<Scalar> > > &v,
            const Dict& opts = Dict());
    static std::vector<std::vector<Matrix<Scalar> > >
    reverse(const std::vector<Matrix<Scalar> > &ex,
            const std::vector<Matrix<Scalar> > &arg,
            const std::vector<std::vector<Matrix<Scalar> > > &v,
            const Dict& opts = Dict());
    static std::vector<bool> which_depends(const Matrix<Scalar> &expr, const Matrix<Scalar> &var,
        casadi_int order=1, bool tr=false);
    static Matrix<Scalar> taylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                   const Matrix<Scalar>& a, casadi_int order);
    static Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                    const Matrix<Scalar>& a, casadi_int order);
    static Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex,
                                    const Matrix<Scalar>& x, const Matrix<Scalar>& a,
                                    casadi_int order,
                                    const std::vector<casadi_int>& order_contributions);
    static Matrix<Scalar> poly_coeff(const Matrix<Scalar>& ex, const Matrix<Scalar>&x);
    static Matrix<Scalar> poly_roots(const Matrix<Scalar>& p);
    static Matrix<Scalar> eig_symbolic(const Matrix<Scalar>& m);
    static Matrix<double> evalf(const Matrix<Scalar>& m);
    static void qr_sparse(const Matrix<Scalar>& A, Matrix<Scalar>& V, Matrix<Scalar>& R,
                          Matrix<Scalar>& beta, std::vector<casadi_int>& prinv,
                          std::vector<casadi_int>& pc, bool amd=true);
    static Matrix<Scalar> qr_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& v,
                                   const Matrix<Scalar>& r, const Matrix<Scalar>& beta,
                                   const std::vector<casadi_int>& prinv,
                                   const std::vector<casadi_int>& pc, bool tr=false);
    static void qr(const Matrix<Scalar>& A, Matrix<Scalar>& Q, Matrix<Scalar>& R);
    static void ldl(const Matrix<Scalar>& A, Matrix<Scalar>& D, Matrix<Scalar>& LT,
                    std::vector<casadi_int>& p, bool amd=true);
    static Matrix<Scalar> ldl_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& D,
                                    const Matrix<Scalar>& LT, const std::vector<casadi_int>& p);
    static Matrix<Scalar> all(const Matrix<Scalar>& x);
    static Matrix<Scalar> any(const Matrix<Scalar>& x);
    static Matrix<Scalar> adj(const Matrix<Scalar>& x);
    static Matrix<Scalar> minor(const Matrix<Scalar>& x, casadi_int i, casadi_int j);
    static Matrix<Scalar> cofactor(const Matrix<Scalar>& A, casadi_int i, casadi_int j);
    static Matrix<Scalar> chol(const Matrix<Scalar>& A);
    static Matrix<Scalar> norm_inf_mul(const Matrix<Scalar>& x, const Matrix<Scalar> &y);
    static Matrix<Scalar> diagcat(const std::vector< Matrix<Scalar> > &A);
    ///@}

    // Simplification with reference counting awareness
    static bool simplify_ref_count(std::vector< Matrix<Scalar> >& arg,
                                   std::vector< Matrix<Scalar> >& res,
                                   const Dict& opts = Dict());

    // Simplification with constant folding
    static bool simplify_const_folding(std::vector< Matrix<Scalar> >& arg,
                                   std::vector< Matrix<Scalar> >& res,
                                   const Dict& opts = Dict());

    // Simplification by combining like terms in linear combinations
    static bool simplify_combine_terms(std::vector< Matrix<Scalar> >& arg,
                                   std::vector< Matrix<Scalar> >& res,
                                   const Dict& opts = Dict());

    /// \endcond
#endif // SWIG

    Matrix<Scalar> printme(const Matrix<Scalar>& y) const;

    /// Transpose the matrix
    Matrix<Scalar> T() const;

#if !defined(SWIG) || defined(DOXYGEN)
/**
\addtogroup expression_tools
@{
*/
    /** \brief Adjugate (classical adjoint) of a square matrix

        Transpose of the cofactor matrix, computed by Laplace expansion. For an invertible
        matrix inv(A) == adj(A)/det(A). Cost grows factorially with the size; intended for
        small (symbolic) matrices. A structural zero cofactor stays structurally zero.
        A non-square input raises an error.

        \lang{python}
        \doctest
        >>> adj(DM([[1, 2], [3, 4]]))
        DM(
        [[4, -2],
         [-3, 1]])
        >>> X = SX.sym("X", 2, 2)
        >>> adj(X)
        SX(
        [[X_3, (-X_2)],
         [(-X_1), X_0]])
        >>> adj(DM(5)), adj(DM(Sparsity.diag(3))).nnz()
        (DM(1), 3)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> adj(DM([1 2; 3 4]))
        [[4, -2],
         [-3, 1]]
        >> X = SX.sym('X', 2, 2);
        >> adj(X)
        [[X_3, (-X_2)],
         [(-X_1), X_0]]
        >> adj(DM(5)), adj(DM(Sparsity.diag(3))).nnz()
        1
        3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.adj(ca.DM([1 2; 3 4]))
        DM(
        [[4, -2],
         [-3, 1]])
        julia> X = ca.SX.sym("X", 2, 2);
        julia> ca.adj(X)
        SX(
        [[X_3, (-X_2)],
         [(-X_1), X_0]])
        julia> ca.adj(ca.DM(5)), ca.nnz(ca.adj(ca.DM(ca.Sparsity.diag(3))))
        (1, 3)
        \enddoctest
        \endlang

        \sa inv_minor, cofactor, det

        \identifier{18p} */
    friend inline Matrix<Scalar> adj(const Matrix<Scalar>& A) {
      return Matrix<Scalar>::adj(A);
    }

    /** \brief Minor of a square matrix: determinant with one column and one row removed

        Note the argument order: \p i is the column that is removed, \p j the row.
        The input must be square; a 1-by-1 input gives 1.

        \param x Square matrix
        \param i Column index to remove (zero-based)
        \param j Row index to remove (zero-based)

        \lang{python}
        \doctest
        >>> A = DM([[1, 2, 3], [0, 1, 4], [5, 6, 0]])
        >>> minor(A, 0, 1)
        DM(-18)
        >>> X = SX.sym("X", 2, 2)
        >>> minor(X, 0, 1)
        SX(X_2)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2 3; 0 1 4; 5 6 0]);
        >> minor(A, 0, 1)
        -18
        >> X = SX.sym('X', 2, 2);
        >> minor(X, 0, 1)
        X_2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2 3; 0 1 4; 5 6 0]);
        julia> ca.minor(A, 0, 1)
        DM(-18)
        julia> X = ca.SX.sym("X", 2, 2);
        julia> ca.minor(X, 0, 1)
        SX(X_2)
        \enddoctest
        \endlang

        \sa cofactor, adj, det

        \identifier{18q} */
    friend inline Matrix<Scalar> minor(const Matrix<Scalar> &x, casadi_int i, casadi_int j) {
      return Matrix<Scalar>::minor(x, i, j);
    }

    /** \brief Cofactor of a square matrix: (-1)^(i+j) times minor(x, i, j)

        Same argument convention as minor: \p i is a column index, \p j a row index.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2, 3], [0, 1, 4], [5, 6, 0]])
        >>> cofactor(A, 0, 1), cofactor(A, 0, 0)
        (DM(18), DM(-24))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2 3; 0 1 4; 5 6 0]);
        >> cofactor(A, 0, 1), cofactor(A, 0, 0)
        18
        -24
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2 3; 0 1 4; 5 6 0]);
        julia> ca.cofactor(A, 0, 1), ca.cofactor(A, 0, 0)
        (18, -24)
        \enddoctest
        \endlang

        \sa minor, adj

        \identifier{18r} */
    friend inline Matrix<Scalar> cofactor(const Matrix<Scalar> &x, casadi_int i, casadi_int j) {
      return Matrix<Scalar>::cofactor(x, i, j);
    }

    /** \brief Dense QR factorization using the modified Gram-Schmidt algorithm

        Factorizes A = Q*R with Q of the same shape as A (orthonormal columns) and R square
        upper triangular. More stable than classical Gram-Schmidt, but breaks down (nan
        entries, no error) when the columns of A are (nearly) linearly dependent.
        A must have at least as many rows as columns. Works symbolically on SX.
        See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1).
        In Python and MATLAB, Q and R are returned as a tuple.

        \lang{python}
        \doctest
        >>> Q, R = qr(DM([[3, 1], [4, 2]]))
        >>> Q
        DM(
        [[0.6, -0.8],
         [0.8, 0.6]])
        >>> R
        DM(
        [[5, 2.2],
         [00, 0.4]])
        >>> mtimes(Q, R)
        DM(
        [[3, 1],
         [4, 2]])
        >>> try: qr(DM([[1, 2, 3]]))
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Q, R] = qr(DM([3 1; 4 2]));
        >> Q
        [[0.6, -0.8],
         [0.8, 0.6]]
        >> R
        [[5, 2.2],
         [00, 0.4]]
        >> Q*R
        [[3, 1],
         [4, 2]]
        >> try, [Q, R] = qr(DM([1 2 3])); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> Q, R = ca.qr(ca.DM([3 1; 4 2]));
        julia> Q
        DM(
        [[0.6, -0.8],
         [0.8, 0.6]])
        julia> R
        DM(
        [[5, 2.2],
         [00, 0.4]])
        julia> Q * R
        DM(
        [[3, 1],
         [4, 2]])
        julia> try ca.qr(ca.DM([1 2 3])); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa qr_sparse, chol, ldl

        \identifier{18s} */
    friend inline void qr(const Matrix<Scalar>& A, Matrix<Scalar>& Q, Matrix<Scalar>& R) {
      return Matrix<Scalar>::qr(A, Q, R);
    }

    /** \brief Sparse direct QR factorization with Householder reflections

        Computes a compact representation of Q (Householder vectors \p V and scaling factors
        \p beta), the upper triangular factor \p R and the row and column permutations
        \p prinv and \p pc that were applied, so that the factors can be reused by qr_solve.
        See T. Davis: Direct Methods for Sparse Linear Systems. Works symbolically on SX.

        \param A Matrix to factorize
        \param amd Use an approximate minimum degree column ordering to reduce fill-in
        \return V, R, beta, prinv, pc (a tuple in Python)

        \lang{python}
        \doctest
        >>> A = DM([[3, 1], [4, 2]])
        >>> V, R, beta, prinv, pc = qr_sparse(A)
        >>> R
        DM(
        [[5, 2.2],
         [00, 0.4]])
        >>> qr_solve(DM([1, 2]), V, R, beta, prinv, pc)
        DM([0, 1])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([3 1; 4 2]);
        >> [V, R, beta, prinv, pc] = qr_sparse(A);
        >> R
        [[5, 2.2],
         [00, 0.4]]
        >> qr_solve(DM([1; 2]), V, R, beta, prinv, pc)
        [0, 1]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([3 1; 4 2]);
        julia> V, R, beta, prinv, pc = ca.qr_sparse(A);
        julia> R
        DM(
        [[5, 2.2],
         [00, 0.4]])
        julia> ca.qr_solve(ca.DM([1, 2]), V, R, beta, prinv, pc)
        DM([0, 1])
        \enddoctest
        \endlang

        \sa qr_solve, Sparsity::qr_sparse, qr

        \identifier{18t} */
    friend inline void qr_sparse(const Matrix<Scalar>& A, Matrix<Scalar>& V, Matrix<Scalar>& R,
                                 Matrix<Scalar>& beta, std::vector<casadi_int>& prinv,
                                 std::vector<casadi_int>& pc, bool amd=true) {
      return Matrix<Scalar>::qr_sparse(A, V, R, beta, prinv, pc, amd);
    }

    /** \brief Solve A*x = b (or A^T*x = b) with the factors returned by qr_sparse

        \param b Right-hand side, one column per system; a sparse \p b is densified
        \param v, r, beta, prinv, pc Output of qr_sparse(A)
        \param tr Solve with the transpose of A instead
        \return Dense solution with the shape of \p b

        \lang{python}
        \doctest
        >>> A = DM([[3, 1], [4, 2]])
        >>> V, R, beta, prinv, pc = qr_sparse(A)
        >>> qr_solve(DM([1, 2]), V, R, beta, prinv, pc, True)
        DM([-3, 2.5])
        >>> solve(A.T, DM([1, 2]))
        DM([-3, 2.5])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([3 1; 4 2]);
        >> [V, R, beta, prinv, pc] = qr_sparse(A);
        >> qr_solve(DM([1; 2]), V, R, beta, prinv, pc, true)
        [-3, 2.5]
        >> solve(A', DM([1; 2]))
        [-3, 2.5]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([3 1; 4 2]);
        julia> V, R, beta, prinv, pc = ca.qr_sparse(A);
        julia> ca.qr_solve(ca.DM([1, 2]), V, R, beta, prinv, pc, true)
        DM([-3, 2.5])
        julia> A' \ [1, 2]
        DM([-3, 2.5])
        \enddoctest
        \endlang

        \sa qr_sparse, solve

        \identifier{18u} */
    friend inline Matrix<Scalar>
    qr_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& v,
             const Matrix<Scalar>& r, const Matrix<Scalar>& beta,
             const std::vector<casadi_int>& prinv,
             const std::vector<casadi_int>& pc, bool tr=false) {
        return Matrix<Scalar>::qr_solve(b, v, r, beta, prinv, pc, tr);
    }

    /** \brief Cholesky factorization of a symmetric positive definite matrix

        Returns the upper triangular factor R with R^T*R == A, obtained from the sparse
        LDL^T factorization (without fill-reducing permutation) as diag(sqrt(D))*L^T.
        The sparsity pattern of A must be symmetric; otherwise an error is raised. An
        indefinite A gives nan entries without an error. Works symbolically on SX.

        \lang{python}
        \doctest
        >>> R = chol(DM([[4, 2], [2, 3]]))
        >>> R
        DM(
        [[2, 1],
         [00, 1.41421]])
        >>> mtimes(R.T, R)
        DM(
        [[4, 2],
         [2, 3]])
        >>> chol(DM(Sparsity.diag(3)) * 4).nnz()
        3
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> R = chol(DM([4 2; 2 3]));
        >> R
        [[2, 1],
         [00, 1.41421]]
        >> R'*R
        [[4, 2],
         [2, 3]]
        >> chol(DM(Sparsity.diag(3)) * 4).nnz()
        3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> R = ca.chol(ca.DM([4 2; 2 3]))
        DM(
        [[2, 1],
         [00, 1.41421]])
        julia> R' * R
        DM(
        [[4, 2],
         [2, 3]])
        julia> ca.nnz(ca.chol(ca.DM(ca.Sparsity.diag(3)) * 4))
        3
        \enddoctest
        \endlang

        \sa ldl, qr

        \identifier{18v} */
    friend inline Matrix<Scalar> chol(const Matrix<Scalar>& A) {
      return Matrix<Scalar>::chol(A);
    }

    /** \brief Sparse LDL^T factorization of a symmetric matrix

        Computes A(p, p) = L*diag(D)*L^T with L unit lower triangular. Returns \p D as a
        vector, \p LT holding the strictly upper triangular entries of L^T (the unit diagonal
        is omitted), and the permutation \p p. Only the upper triangle of A is referenced.
        The sparsity pattern of A must be symmetric. Only guaranteed to work for positive
        definite matrices (no pivoting). Works symbolically on SX.

        \param A Square matrix with symmetric sparsity pattern
        \param amd Use an approximate minimum degree ordering to reduce fill-in
        \return D, LT, p (a tuple in Python)

        \lang{python}
        \doctest
        >>> D, LT, p = ldl(DM([[4, 2], [2, 3]]))
        >>> D, p
        (DM([4, 2]), [0, 1])
        >>> LT
        DM(
        [[00, 0.5],
         [00, 00]])
        >>> L = DM.eye(2) + LT.T
        >>> mtimes(mtimes(L, diag(D)), L.T)
        DM(
        [[4, 2],
         [2, 3]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [D, LT, p] = ldl(DM([4 2; 2 3]));
        >> D, p
        [4, 2]
        0 1
        >> LT
        [[00, 0.5],
         [00, 00]]
        >> L = DM.eye(2) + LT';
        >> L*diag(D)*L'
        [[4, 2],
         [2, 3]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> D, LT, p = ca.ldl(ca.DM([4 2; 2 3]));
        julia> D
        DM([4, 2])
        julia> p
        2-element Vector{Int64}:
         0
         1
        julia> LT
        DM(
        [[00, 0.5],
         [00, 00]])
        julia> L = ca.DM.eye(2) + LT';
        julia> L * ca.diag(D) * L'
        DM(
        [[4, 2],
         [2, 3]])
        \enddoctest
        \endlang

        \sa ldl_solve, chol, Sparsity::ldl

        \identifier{18w} */
    friend inline void ldl(const Matrix<Scalar>& A, Matrix<Scalar>& D, Matrix<Scalar>& LT,
                           std::vector<casadi_int>& p, bool amd=true) {
      return Matrix<Scalar>::ldl(A, D, LT, p, amd);
    }

    /** \brief Solve A*x = b with the factors returned by ldl

        \param b Right-hand side, one column per system; a sparse \p b is densified
        \param D, LT, p Output of ldl(A)
        \return Dense solution with the shape of \p b

        \lang{python}
        \doctest
        >>> A = DM([[4, 2], [2, 3]])
        >>> D, LT, p = ldl(A)
        >>> ldl_solve(DM([1, 2]), D, LT, p)
        DM([-0.125, 0.75])
        >>> solve(A, DM([1, 2]))
        DM([-0.125, 0.75])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([4 2; 2 3]);
        >> [D, LT, p] = ldl(A);
        >> ldl_solve(DM([1; 2]), D, LT, p)
        [-0.125, 0.75]
        >> solve(A, DM([1; 2]))
        [-0.125, 0.75]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([4 2; 2 3]);
        julia> D, LT, p = ca.ldl(A);
        julia> ca.ldl_solve(ca.DM([1, 2]), D, LT, p)
        DM([-0.125, 0.75])
        julia> A \ [1, 2]
        DM([-0.125, 0.75])
        \enddoctest
        \endlang

        \sa ldl, solve

        \identifier{18x} */
    friend inline Matrix<Scalar>
    ldl_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& D, const Matrix<Scalar>& LT,
              const std::vector<casadi_int>& p) {
      return Matrix<Scalar>::ldl_solve(b, D, LT, p);
    }

    /** \brief True if any entry equals one

        Returns a 1-by-1 matrix. A matrix with structural zeros gives false.
        Not available in Python (shadowed by the builtin); use logic_any or sum1/sum2.

        \lang{matlab}
        \doctest
        >> any(DM([0; 1])), any(DM([0; 0]))
        1
        0
        >> any(DM(Sparsity.diag(2)))
        0
        \enddoctest
        \endlang

        \sa all, logic_any

        \identifier{18y} */
    friend inline Matrix<Scalar> any(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::any(x);
    }

    /** \brief True if every entry equals one

        Returns a 1-by-1 matrix. A matrix with structural zeros gives false.
        Not available in Python (shadowed by the builtin); use logic_all.

        \lang{matlab}
        \doctest
        >> all(DM([1; 1])), all(DM([0; 1]))
        1
        0
        >> all(DM(Sparsity.diag(2)))
        0
        \enddoctest
        \endlang

        \sa any, logic_all

        \identifier{18z} */
    friend inline Matrix<Scalar> all(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::all(x);
    }

    /** \brief Infinity norm of a matrix product, without forming the product

        Equals norm_inf(mtimes(x, y)): the largest absolute entry of x*y. Exploits sparsity;
        a shape mismatch raises an error. Works symbolically on SX.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> norm_inf_mul(A, A), norm_inf(mtimes(A, A))
        (DM(22), DM(22))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> norm_inf_mul(A, A), norm_inf(A*A)
        22
        22
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> ca.norm_inf_mul(A, A), ca.norm_inf(A * A)
        (22, 22)
        \enddoctest
        \endlang

        \sa norm_inf, mtimes

        \identifier{190} */
    friend inline Matrix<Scalar>
      norm_inf_mul(const Matrix<Scalar> &x, const Matrix<Scalar> &y) {
      return Matrix<Scalar>::norm_inf_mul(x, y);
    }

    /** \brief Make a matrix sparse by turning numerical zeros into structural zeros

        Every stored entry with absolute value at most \p tol is dropped from the sparsity
        pattern. With the default tol=0 only exact zeros are removed; SX entries that the
        symbolic simplifications already reduced to the constant 0 count as zeros.
        Returns a new matrix; the input is unchanged.

        \param A Matrix to compact
        \param tol Absolute threshold below which an entry counts as zero

        \lang{python}
        \doctest
        >>> sparsify(DM([[1, 0], [0, 2]]))
        DM(
        [[1, 00],
         [00, 2]])
        >>> A = DM([[1, 0], [1e-9, 2]])
        >>> sparsify(A).nnz(), sparsify(A, 1e-6).nnz()
        (3, 2)
        >>> x = SX.sym("x")
        >>> sparsify(vertcat(x, 0, x * 0))
        SX([x, 00, 00])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sparsify(DM([1 0; 0 2]))
        [[1, 00],
         [00, 2]]
        >> A = DM([1 0; 1e-9 2]);
        >> fprintf('%d %d\n', sparsify(A).nnz(), sparsify(A, 1e-6).nnz())
        3 2
        >> x = SX.sym('x');
        >> sparsify([x; 0; x * 0])
        [x, 00, 00]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.sparsify(ca.DM([1 0; 0 2]))
        DM(
        [[1, 00],
         [00, 2]])
        julia> A = ca.DM([1 0; 1e-9 2]);
        julia> ca.nnz(ca.sparsify(A)), ca.nnz(ca.sparsify(A, 1e-6))
        (3, 2)
        julia> x = ca.SX.sym("x");
        julia> ca.sparsify([x; 0; x * 0])
        SX([x, 00, 00])
        \enddoctest
        \endlang

        \sa densify, has_zeros

        \identifier{191} */
    friend inline Matrix<Scalar>
      sparsify(const Matrix<Scalar>& A, double tol=0) {
      return Matrix<Scalar>::sparsify(A, tol);
    }

    /** \brief Expand a scalar SX expression as a weighted sum of terms with constant weights

        Collects additions, subtractions and multiplications by constants so that
        ex == dot(weights, terms), merging repeated terms. Any other operation (a product of
        two non-constant factors, a function call, ...) is kept as a single term with
        weight 1. The input must be a scalar SX; DM is not supported.

        \return weights and terms, both column vectors (a tuple in Python)

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> w, t = expand(3*x - x + 5*y)
        >>> w, t
        (SX([2, 5]), SX([x, y]))
        >>> w, t = expand(3*(x + y) + 1)
        >>> t
        SX([x, y, 1])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> [w, t] = expand(3*x - x + 5*y);
        >> w, t
        [2, 5]
        [x, y]
        >> [w, t] = expand(3*(x + y) + 1);
        >> t
        [x, y, 1]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> w, t = ca.expand(3*x - x + 5*y);
        julia> w
        SX([2, 5])
        julia> t
        SX([x, y])
        julia> w, t = ca.expand(3*(x + y) + 1);
        julia> t
        SX([x, y, 1])
        \enddoctest
        \endlang

        \sa linear_coeff

        \identifier{192} */
    friend inline void expand(const Matrix<Scalar>& ex, Matrix<Scalar> &weights,
                              Matrix<Scalar>& terms) {
      Matrix<Scalar>::expand(ex, weights, terms);
    }

    /** \brief Piecewise constant function of a scalar symbolic variable

        Builds an SX expression that equals val[0] for t < tval[0], val[k] for
        tval[k-1] <= t < tval[k] and val[n-1] for t >= tval[n-2], with n = val.numel().
        \p tval must be increasing. The result is a sum of step functions, so it is not
        smooth; its derivative with respect to \p t is zero almost everywhere.
        Only defined for SX.

        \param t Scalar SX variable (e.g. time)
        \param tval Breakpoints, length n-1; a length mismatch raises an error
        \param val Function value on each interval, length n

        \lang{python}
        \doctest
        >>> t = SX.sym("t")
        >>> e = pw_const(t, [0, 2], [7, 1, 3])
        >>> f = Function("f", [t], [e])
        >>> [float(f(v)) for v in [-1, 0, 1, 2, 5]]
        [7.0, 1.0, 1.0, 3.0, 3.0]
        >>> try: pw_const(t, [0], [7, 1, 3])
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> t = SX.sym('t');
        >> e = pw_const(t, [0 2], [7 1 3]);
        >> f = Function('f', {t}, {e});
        >> f([-1; 0; 1; 2; 5])
        [7, 1, 1, 3, 3]
        >> try, pw_const(t, 0, [7 1 3]); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> t = ca.SX.sym("t");
        julia> e = ca.pw_const(t, [0, 2], [7, 1, 3]);
        julia> f = ca.Function("f", [t], [e]);
        julia> [Float64(f(v)) for v in [-1, 0, 1, 2, 5]]
        5-element Vector{Float64}:
         7.0
         1.0
         1.0
         3.0
         3.0
        julia> try ca.pw_const(t, [0], [7, 1, 3]); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa pw_lin, if_else, interpolant

        \identifier{193} */
    friend inline Matrix<Scalar> pw_const(const Matrix<Scalar> &t,
                                            const Matrix<Scalar> &tval,
                                            const Matrix<Scalar> &val) {
      return Matrix<Scalar>::pw_const(t, tval, val);
    }

    /** \brief Piecewise linear function of a scalar symbolic variable

        Builds an SX expression that interpolates linearly between the points
        (tval[k], val[k]). Outside [tval[0], tval[N-1]] the first and last segments are
        extrapolated. Only defined for SX.

        \param t Scalar SX variable (e.g. time)
        \param tval Breakpoints, increasing, at least two
        \param val Function values at the breakpoints, same length as \p tval

        \lang{python}
        \doctest
        >>> t = SX.sym("t")
        >>> f = Function("f", [t], [pw_lin(t, [0, 1, 3], [0, 2, 1])])
        >>> [float(f(v)) for v in [-1, 0, 0.5, 1, 2, 3, 4]]
        [-2.0, 0.0, 1.0, 2.0, 1.5, 1.0, 0.5]
        >>> try: pw_lin(t, [0], [7])
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> t = SX.sym('t');
        >> f = Function('f', {t}, {pw_lin(t, [0 1 3], [0 2 1])});
        >> f([-1; 0; 0.5; 1; 2; 3; 4])
        [-2, 0, 1, 2, 1.5, 1, 0.5]
        >> try, pw_lin(t, 0, 7); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> t = ca.SX.sym("t");
        julia> f = ca.Function("f", [t], [ca.pw_lin(t, [0, 1, 3], [0, 2, 1])]);
        julia> [Float64(f(v)) for v in [-1, 0, 0.5, 1, 2, 3, 4]]
        7-element Vector{Float64}:
         -2.0
          0.0
          1.0
          2.0
          1.5
          1.0
          0.5
        julia> try ca.pw_lin(t, [0], [7]); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa pw_const, interpolant

        \identifier{194} */
    friend inline Matrix<Scalar>
      pw_lin(const Matrix<Scalar> &t, const Matrix<Scalar> &tval,
             const Matrix<Scalar> &val) {
      return Matrix<Scalar>::pw_lin(t, tval, val);
    }

    /** \brief Heaviside step function, elementwise

        H(x) = 0 for x < 0, 1/2 for x == 0 and 1 for x > 0, computed as (1 + sign(x))/2.
        Structural zeros of a sparse input become 1/2. Not defined for MX.

        \lang{python}
        \doctest
        >>> heaviside(DM([-1, 0, 1]))
        DM([0, 0.5, 1])
        >>> heaviside(SX.sym("x"))
        SX(((1+sign(x))/2))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> heaviside(DM([-1; 0; 1]))
        [0, 0.5, 1]
        >> heaviside(SX.sym('x'))
        ((1+sign(x))/2)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.heaviside(ca.DM([-1, 0, 1]))
        DM([0, 0.5, 1])
        julia> ca.heaviside(ca.SX.sym("x"))
        SX(((1+sign(x))/2))
        \enddoctest
        \endlang

        \sa ramp, rectangle, sign

        \identifier{195} */
    friend inline Matrix<Scalar> heaviside(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::heaviside(x);
    }

    /** \brief Rectangle (gate, pulse, window) function, elementwise

        1 for |x| < 1/2, 1/2 for |x| == 1/2 and 0 for |x| > 1/2, built from sign.
        Not defined for MX.

        \lang{python}
        \doctest
        >>> rectangle(DM([-1, -0.5, 0, 0.5, 1]))
        DM([0, 0.5, 1, 0.5, 0])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> rectangle(DM([-1; -0.5; 0; 0.5; 1]))
        [0, 0.5, 1, 0.5, 0]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.rectangle(ca.DM([-1, -0.5, 0, 0.5, 1]))
        DM([0, 0.5, 1, 0.5, 0])
        \enddoctest
        \endlang

        \sa triangle, heaviside

        \identifier{23n} */
    friend inline Matrix<Scalar> rectangle(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::rectangle(x);
    }

    /** \brief Triangle (hat) function, elementwise

        1 - |x| for |x| < 1 and 0 otherwise, built from rectangle and fabs. Not defined for MX.

        \lang{python}
        \doctest
        >>> triangle(DM([-1, -0.5, 0, 0.5, 1]))
        DM([0, 0.5, 1, 0.5, 0])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> triangle(DM([-1; -0.5; 0; 0.5; 1]))
        [0, 0.5, 1, 0.5, 0]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.triangle(ca.DM([-1, -0.5, 0, 0.5, 1]))
        DM([0, 0.5, 1, 0.5, 0])
        \enddoctest
        \endlang

        \sa rectangle, ramp

        \identifier{23o} */
    friend inline Matrix<Scalar> triangle(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::triangle(x);
    }

    /** \brief Ramp function, elementwise

        x for x > 0 and 0 for x <= 0, computed as x*heaviside(x). Also called the slope
        function or relu. Not defined for MX; use fmax(x, 0) for an MX-compatible alternative.

        \lang{python}
        \doctest
        >>> ramp(DM([0, 0.5, 3]))
        DM([0, 0.5, 3])
        >>> ramp(SX.sym("x"))
        SX((x*((1+sign(x))/2)))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> ramp(DM([0; 0.5; 3]))
        [0, 0.5, 3]
        >> ramp(SX.sym('x'))
        (x*((1+sign(x))/2))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.ramp(ca.DM([0, 0.5, 3]))
        DM([0, 0.5, 3])
        julia> ca.ramp(ca.SX.sym("x"))
        SX((x*((1+sign(x))/2)))
        \enddoctest
        \endlang

        \sa heaviside, fmax

        \identifier{23p} */
    friend inline Matrix<Scalar> ramp(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::ramp(x);
    }

    ///@{
    /** \brief Integrate a scalar SX expression over [a, b] with 5-point Gauss-Legendre quadrature

        Evaluates \p f at the five Gauss points of [a, b] and returns the weighted sum,
        which is exact for polynomials in \p x up to degree 9. Only defined for SX.

        \param f Scalar expression to integrate; it may depend only on \p x
        \param x Scalar symbolic integration variable
        \param a Lower bound, a numeric constant
        \param b Upper bound, a numeric constant
        \param order Must be 5; any other value raises an error
        \param w Reserved; must be empty

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> try: gauss_quadrature(x**2, x, 0, 1, 3)
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> try, gauss_quadrature(x^2, x, 0, 1, 3); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> try ca.gauss_quadrature(x^2, x, 0, 1, 3); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa integrator

        \identifier{196} */
    friend inline Matrix<Scalar>
      gauss_quadrature(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                       const Matrix<Scalar> &a, const Matrix<Scalar> &b,
                       casadi_int order=5) {
      return Matrix<Scalar>::gauss_quadrature(f, x, a, b, order);
    }
    friend inline Matrix<Scalar>
      gauss_quadrature(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                       const Matrix<Scalar> &a, const Matrix<Scalar> &b,
                       casadi_int order, const Matrix<Scalar>& w) {
      return Matrix<Scalar>::gauss_quadrature(f, x, a, b, order, w);
    }
    ///@}

    ///@{
    /** \brief Univariate Taylor series expansion

        Expands \p ex around \p a with respect to the scalar symbol \p x up to the given
        order:
        ex(a) + ex'(a)*(x-a) + ex''(a)*(x-a)^2/2! + ... + ex^(order)(a)*(x-a)^order/order!
        The expansion is done entrywise for a dense matrix \p ex; sparse \p ex raises an
        error. \p a may itself be symbolic. The two-argument form expands around 0 to first
        order. Only defined for SX.

        \param ex Dense SX expression to expand
        \param x Scalar SX symbol
        \param a Scalar expansion point, default 0
        \param order Highest power of (x-a) kept, default 1; 0 gives ex(a)

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); a = SX.sym("a")
        >>> taylor(sin(x), x)
        SX(x)
        >>> taylor(sin(x), x, a, 1)
        SX((sin(a)+(cos(a)*(x-a))))
        >>> taylor(exp(x), x, 0, 3)
        SX(@1=sq(x), (((1+x)+(0.5*@1))+(0.166667*(@1*x))))
        >>> taylor(vertcat(sin(x), cos(x)), x, 0, 2)
        SX([x, (1+(-0.5*sq(x)))])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); a = SX.sym('a');
        >> taylor(sin(x), x)
        x
        >> taylor(sin(x), x, a, 1)
        (sin(a)+(cos(a)*(x-a)))
        >> taylor(exp(x), x, 0, 3)
        @1=sq(x), (((1+x)+(0.5*@1))+(0.166667*(@1*x)))
        >> taylor([sin(x); cos(x)], x, 0, 2)
        [x, (1+(-0.5*sq(x)))]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); a = ca.SX.sym("a");
        julia> ca.taylor(sin(x), x)
        SX(x)
        julia> ca.taylor(sin(x), x, a, 1)
        SX((sin(a)+(cos(a)*(x-a))))
        julia> ca.taylor(exp(x), x, 0, 3)
        SX(@1=sq(x), (((1+x)+(0.5*@1))+(0.166667*(@1*x))))
        julia> ca.taylor([sin(x); cos(x)], x, 0, 2)
        SX([x, (1+(-0.5*sq(x)))])
        \enddoctest
        \endlang

        \sa mtaylor, linearize, linear_coeff

        \identifier{23q} */
    friend inline Matrix<Scalar> taylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                          const Matrix<Scalar>& a, casadi_int order=1) {
      return Matrix<Scalar>::taylor(ex, x, a, order);
    }
    friend inline Matrix<Scalar> taylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x) {
      return Matrix<Scalar>::taylor(ex, x, 0, 1);
    }
    ///@}

    /** \brief Multivariate Taylor series expansion

        Expands \p ex around the point \p a with respect to the symbols in \p x, keeping all
        terms whose total degree is at most \p order. The total degree of x^n*y^m is n+m.
        \p ex and \p x must be dense; \p a must have the same number of entries as \p x.
        Only defined for SX.

        \param ex Dense SX expression to expand
        \param x Vector of SX symbols
        \param a Expansion point, a vector of the same length as \p x
        \param order Maximum total degree, default 1 (linearization)

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y"); a = SX.sym("a"); b = SX.sym("b")
        >>> mtaylor(sin(x + y), vertcat(x, y), vertcat(a, b), 1)
        SX(((sin((a+b))+(cos((a+b))*(x-a)))+(cos((a+b))*(y-b))))
        >>> mtaylor(exp(x * y), vertcat(x, y), DM([0, 0]), 2)
        SX(@1=2, ((1+((x*y)/@1))+((y*x)/@1)))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y'); a = SX.sym('a'); b = SX.sym('b');
        >> mtaylor(sin(x + y), [x; y], [a; b], 1)
        ((sin((a+b))+(cos((a+b))*(x-a)))+(cos((a+b))*(y-b)))
        >> mtaylor(exp(x * y), [x; y], DM([0; 0]), 2)
        @1=2, ((1+((x*y)/@1))+((y*x)/@1))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y"); a = ca.SX.sym("a"); b = ca.SX.sym("b");
        julia> ca.mtaylor(sin(x + y), [x; y], [a; b], 1)
        SX(((sin((a+b))+(cos((a+b))*(x-a)))+(cos((a+b))*(y-b))))
        julia> ca.mtaylor(exp(x * y), [x; y], ca.DM([0, 0]), 2)
        SX(@1=2, ((1+((x*y)/@1))+((y*x)/@1)))
        \enddoctest
        \endlang

        \sa taylor

        \identifier{23r} */
    friend inline Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                           const Matrix<Scalar>& a, casadi_int order=1) {
      return Matrix<Scalar>::mtaylor(ex, x, a, order);
    }

    /** \brief Multivariate Taylor series expansion with weighted degrees

        Same as mtaylor, but each variable contributes a given weight to the total degree:
        with x = [x, y] and order_contributions = [1, 2] the degree of x^n*y^m is n + 2*m.
        The length of \p order_contributions must equal the number of entries of \p x.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> mtaylor(sin(x + y), vertcat(x, y), vertcat(0, 0), 3, [1, 2])
        SX(((x-((sq(x)*x)/6))+y))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> mtaylor(sin(x + y), [x; y], [0; 0], 3, [1 2])
        ((x-((sq(x)*x)/6))+y)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.mtaylor(sin(x + y), [x; y], [0; 0], 3, [1, 2])
        SX(((x-((sq(x)*x)/6))+y))
        \enddoctest
        \endlang

        \sa taylor

        \identifier{23s} */
    friend inline Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                           const Matrix<Scalar>& a, casadi_int order,
                                           const std::vector<casadi_int>& order_contributions) {
      return Matrix<Scalar>::mtaylor(ex, x, a, order, order_contributions);
    }

    /** \brief Extract the polynomial coefficients of a scalar expression in one symbol

        Returns the coefficients in decreasing power order (highest degree first), as a
        column vector; the constant term is last. Coefficients may be symbolic in other
        variables. Found by repeated differentiation, so an expression that is not polynomial
        in \p x raises an error. Only defined for SX.

        \param ex Scalar SX expression that is a polynomial in \p x
        \param x Scalar SX symbol

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); a = SX.sym("a")
        >>> poly_coeff(3*x**2 + 2*x + 1, x)
        SX([3, 2, 1])
        >>> poly_coeff((x - a)*(x + a), x)
        SX([1, 0, (-sq(a))])
        >>> try: poly_coeff(sin(x), x)
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); a = SX.sym('a');
        >> poly_coeff(3*x^2 + 2*x + 1, x)
        [3, 2, 1]
        >> poly_coeff((x - a)*(x + a), x)
        [1, 0, (-sq(a))]
        >> try, poly_coeff(sin(x), x); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); a = ca.SX.sym("a");
        julia> ca.poly_coeff(3*x^2 + 2*x + 1, x)
        SX([3, 2, 1])
        julia> ca.poly_coeff((x - a)*(x + a), x)
        SX([1, 0, (-sq(a))])
        julia> try ca.poly_coeff(sin(x), x); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa poly_roots, polyval

        \identifier{197} */
    friend inline Matrix<Scalar> poly_coeff(const Matrix<Scalar>& f,
                                              const Matrix<Scalar>& x) {
      return Matrix<Scalar>::poly_coeff(f, x);
    }

    /** \brief Roots of a polynomial given by its coefficients, in closed form

        \p p is a dense column vector of coefficients in decreasing power order, as returned
        by poly_coeff. Closed-form formulas are used for degrees 1 to 4; a zero constant term
        is split off as a root 0 first. Higher degrees raise an error. All roots are assumed
        real: complex roots give nan. The result stays symbolic for symbolic coefficients.
        Only defined for SX; pass SX(p) for numeric coefficients.

        \lang{python}
        \doctest
        >>> poly_roots(SX([1, -3, 2]))
        SX([1, 2])
        >>> evalf(poly_roots(SX([1, -6, 11, -6])))
        DM([3, 1, 2])
        >>> a = SX.sym("a")
        >>> poly_roots(vertcat(a, 2))
        SX((-2/a))
        >>> try: poly_roots(SX([1, 1, 1, 1, 1, 1]))
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> poly_roots(SX([1; -3; 2]))
        [1, 2]
        >> evalf(poly_roots(SX([1; -6; 11; -6])))
        [3, 1, 2]
        >> a = SX.sym('a');
        >> poly_roots([a; 2])
        (-2/a)
        >> try, poly_roots(SX([1; 1; 1; 1; 1; 1])); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.poly_roots(ca.SX([1, -3, 2]))
        SX([1, 2])
        julia> ca.evalf(ca.poly_roots(ca.SX([1, -6, 11, -6])))
        DM([3, 1, 2])
        julia> a = ca.SX.sym("a");
        julia> ca.poly_roots([a; 2])
        SX((-2/a))
        julia> try ca.poly_roots(ca.SX([1, 1, 1, 1, 1, 1])); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa poly_coeff, eig_symbolic

        \identifier{198} */
    friend inline Matrix<Scalar> poly_roots(const Matrix<Scalar>& p) {
      return Matrix<Scalar>::poly_roots(p);
    }

    /** \brief Eigenvalues of a square SX matrix in closed form

        Splits the matrix into its strongly connected (block diagonal) components and solves
        the characteristic polynomial of each block with poly_roots, so every block may be at
        most 4-by-4. Larger blocks raise an error. Eigenvalues are assumed real.
        Returns a column vector; the order follows the block structure, not the magnitude.
        Only defined for SX; pass SX(m) for a numeric matrix.

        \lang{python}
        \doctest
        >>> evalf(eig_symbolic(SX([[2, 1], [1, 2]])))
        DM([1, 3])
        >>> M = SX.sym("M", 2, 2)
        >>> eig_symbolic(M)
        SX(@1=(M_3+M_0), @2=sqrt((sq(@1)-(4*((M_0*M_3)-(M_2*M_1))))), @3=2, [((@1-@2)/@3),
        ((@1+@2)/@3)])
        >>> evalf(eig_symbolic(SX.eye(5) * 3)).T
        DM([[3, 3, 3, 3, 3]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> evalf(eig_symbolic(SX([2 1; 1 2])))
        [1, 3]
        >> a = SX.sym('a'); b = SX.sym('b'); c = SX.sym('c'); d = SX.sym('d');
        >> eig_symbolic([a b; c d])
        @1=(d+a), @2=sqrt((sq(@1)-(4*((a*d)-(b*c))))), @3=2, [((@1-@2)/@3), ((@1+@2)/@3)]
        >> evalf(eig_symbolic(SX.eye(5) * 3))'
        [[3, 3, 3, 3, 3]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.evalf(ca.eig_symbolic(ca.SX([2 1; 1 2])))
        DM([1, 3])
        julia> M = ca.SX.sym("M", 2, 2);
        julia> l = ca.eig_symbolic(M); size(l)
        (2, 1)
        julia> l[1]
        SX(@1=(M_3+M_0), ((@1-sqrt((sq(@1)-(4*((M_0*M_3)-(M_2*M_1))))))/2))
        julia> l[2]
        SX(@1=(M_3+M_0), ((@1+sqrt((sq(@1)-(4*((M_0*M_3)-(M_2*M_1))))))/2))
        julia> ca.evalf(ca.eig_symbolic(ca.SX.eye(5) * 3))'
        DM([[3, 3, 3, 3, 3]])
        \enddoctest
        \endlang

        \sa poly_roots, Sparsity::scc

        \identifier{199} */
    friend inline Matrix<Scalar> eig_symbolic(const Matrix<Scalar>& m) {
      return Matrix<Scalar>::eig_symbolic(m);
    }


    /** \brief Evaluate a constant expression numerically, returning a DM

        Works on SX and MX expressions without free symbols (including constant results of
        function calls). The sparsity pattern is preserved. An expression that still
        depends on a symbol raises an error, unlike the DM(SX) constructor which yields nan.

        \lang{python}
        \doctest
        >>> evalf(SX([[1, 2], [3, 4]]) * 2)
        DM(
        [[2, 4],
         [6, 8]])
        >>> evalf(MX(DM([1, 2])) + 1)
        DM([2, 3])
        >>> try: evalf(SX.sym("x"))
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> evalf(SX([1 2; 3 4]) * 2)
        [[2, 4],
         [6, 8]]
        >> evalf(MX(DM([1; 2])) + 1)
        [2, 3]
        >> try, evalf(SX.sym('x')); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.evalf(ca.SX([1 2; 3 4]) * 2)
        DM(
        [[2, 4],
         [6, 8]])
        julia> ca.evalf(ca.MX(ca.DM([1, 2])) + 1)
        DM([2, 3])
        julia> try ca.evalf(ca.SX.sym("x")); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa substitute, Function

        \identifier{19a} */
    inline friend Matrix<double> evalf(const Matrix<Scalar>& expr) {
      return Matrix<Scalar>::evalf(expr);
    }
/** @} */
#endif

    /** \brief Set the depth to which SX subexpressions are compared during simplification

        The on-the-fly simplifications (x - x -> 0, x/x -> 1, x == x -> 1, ...) recognise two
        operands as equal when they are the same node, or, up to this depth, when they are
        the same operation applied to equal operands. Depth 1 (the default) recognises
        (x+y) - (x+y); depth 2 also sin(x+y) - sin(x+y). Depth 0 compares nodes only.
        Global setting, SX only; calling without an argument restores the default.

        \param eq_depth Comparison depth, default 1

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> sin(x + y) - sin(x + y)
        SX((sin((x+y))-sin((x+y))))
        >>> SX.set_max_depth(2)
        >>> sin(x + y) - sin(x + y)
        SX(0)
        >>> SX.set_max_depth()
        >>> SX.get_max_depth()
        1
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> sin(x + y) - sin(x + y)
        (sin((x+y))-sin((x+y)))
        >> SX.set_max_depth(2);
        >> sin(x + y) - sin(x + y)
        0
        >> SX.set_max_depth();
        >> SX.get_max_depth()
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> sin(x + y) - sin(x + y)
        SX((sin((x+y))-sin((x+y))))
        julia> ca.SX.set_max_depth(2)
        julia> sin(x + y) - sin(x + y)
        SX(0)
        julia> ca.SX.set_max_depth()
        julia> ca.SX.get_max_depth()
        1
        \enddoctest
        \endlang

        \sa get_max_depth, is_equal

        \identifier{19b} */
    static void set_max_depth(casadi_int eq_depth=1);

    /** \brief Get the depth to which SX subexpressions are compared during simplification

        \sa set_max_depth

        \identifier{19c} */
    static casadi_int get_max_depth();

    /** \brief Get the symbolic input expressions of an SX Function

        Same as Function::sx_in. Only defined for SX.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> SX.get_input(Function("f", [x], [x + 1]))
        [SX(x)]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> v = SX.get_input(Function('f', {x}, {x + 1}));
        >> numel(v), v{1}
        1
        x
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> ca.SX.get_input(ca.Function("f", [x], [x + 1]))
        1-element Vector{SX}:
         SX(x)
        \enddoctest
        \endlang

        \sa Function::sx_in, get_free

        \identifier{19d} */
    static std::vector<Matrix<Scalar> > get_input(const Function& f);

    /** \brief Get the free symbols of an SX Function

        Same as Function::free_sx. Only defined for SX.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> SX.get_free(Function("f", [x], [x + y], {"allow_free": True}))
        [SX(y)]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> v = SX.get_free(Function('f', {x}, {x + y}, struct('allow_free', true)));
        >> v{1}
        y
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.SX.get_free(ca.Function("f", [x], [x + y], Dict("allow_free" => true)))
        1-element Vector{SX}:
         SX(y)
        \enddoctest
        \endlang

        \sa Function::free_sx, get_input

        \identifier{19e} */
    static std::vector<Matrix<Scalar> > get_free(const Function& f);

    /** \brief Name of the class: "DM" or "SX"

        \lang{python}
        \doctest
        >>> DM.type_name(), SX.type_name()
        ('DM', 'SX')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%s %s\n', DM.type_name(), SX.type_name())
        DM SX
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM.type_name(), ca.SX.type_name()
        ("DM", "SX")
        \enddoctest
        \endlang

        \identifier{2ox} */
    static std::string type_name();

    /** \brief String representation of each nonzero and of the shared subexpressions

        Returns two lists: one string per nonzero (in storage order) and the definitions of
        the intermediate expressions "@1", "@2", ... that the nonzero strings refer to
        (always empty for DM).

        \lang{python}
        \doctest
        >>> DM([[1, 2], [3, 4]]).print_split()
        (['1', '3', '2', '4'], [])
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> e = x * y
        >>> vertcat(e + 1, e - 1).print_split()
        (['(@1+@2)', '(@1-@2)'], ['(x*y)', '1'])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [nz, inter] = DM([1 2; 3 4]).print_split();
        >> disp(strjoin(cellstr(nz)', ' ')), disp(isempty(inter))
        1 3 2 4
        1
        >> x = SX.sym('x'); y = SX.sym('y');
        >> e = x * y;
        >> [nz, inter] = vertcat(e + 1, e - 1).print_split();
        >> disp(strjoin(cellstr(nz)', ' ')), disp(strjoin(cellstr(inter)', ' '))
        (@1+@2) (@1-@2)
        (x*y) 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.print_split(ca.DM([1 2; 3 4]))
        (["1", "3", "2", "4"], String[])
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> e = x * y;
        julia> ca.print_split([e + 1; e - 1])
        (["(@1+@2)", "(@1-@2)"], ["(x*y)", "1"])
        \enddoctest
        \endlang

        \identifier{2oy} */
    void print_split(std::vector<std::string>& SWIG_OUTPUT(nz),
                    std::vector<std::string>& SWIG_OUTPUT(inter)) const;

    /** \brief Print a representation of the object to a stream

        Picks the layout automatically: scalar, vector, dense matrix (small or at least
        half full) or sparse coordinate list. This is what str() and print() use.
        \p more has no effect for matrices.

        \identifier{2oz} */
    void disp(std::ostream& stream, bool more=false) const;

    /// Get the string representation that disp prints
    std::string get_str(bool more=false) const;

    /** \brief Print a 1-by-1 matrix as a bare number or expression

        A structurally empty 1-by-1 matrix prints as 00; other shapes raise an error.

        \lang{python}
        \doctest
        >>> DM(3).print_scalar()
        3
        >>> DM(1, 1).print_scalar()
        00
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.print_scalar(ca.DM(3))
        3
        julia> ca.print_scalar(ca.DM(1, 1))
        00
        \enddoctest
        \endlang

        \identifier{2p0} */
    void print_scalar(std::ostream &stream) const;

    /** \brief Print a column vector as [a, b, ...], structural zeros as 00

        Only column vectors are accepted; a row vector raises an error.

        \param truncate Print only the first and last 3 entries of vectors longer than 1000

        \lang{python}
        \doctest
        >>> DM(Sparsity.unit(3, 1)).print_vector()
        [00, 1, 00]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.print_vector(ca.DM(ca.Sparsity.unit(3, 1)))
        [00, 1, 00]
        \enddoctest
        \endlang

        \identifier{2p1} */
    void print_vector(std::ostream &stream, bool truncate=true) const;

    /** \brief Print as a dense matrix, one row per line, structural zeros as 00

        Shared SX subexpressions are printed first as @1=..., @2=....

        \param truncate Print only 3 leading and trailing rows and columns when the matrix
                        has more than 1000 elements and at least 7 rows or columns

        \lang{python}
        \doctest
        >>> DM(Sparsity.diag(2)).print_dense()
        [[1, 00],
         [00, 1]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.print_dense(ca.DM(ca.Sparsity.diag(2)))
        [[1, 00],
         [00, 1]]
        \enddoctest
        \endlang

        \sa print_sparse

        \identifier{2p2} */
    void print_dense(std::ostream &stream, bool truncate=true) const;

    /** \brief Print as a list of (row, col) -> value triplets, column-major order

        \param truncate Print only the first and last 3 nonzeros when there are more than 1000

        \lang{python}
        \doctest
        >>> DM(Sparsity.diag(2)).print_sparse()
        sparse: 2-by-2, 2 nnz
         (0, 0) -> 1
         (1, 1) -> 1
        >>> DM(2, 2).print_sparse()
        all zero sparse: 2-by-2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.print_sparse(ca.DM(ca.Sparsity.diag(2)))
        sparse: 2-by-2, 2 nnz
         (0, 0) -> 1
         (1, 1) -> 1
        julia> ca.print_sparse(ca.DM(2, 2))
        all zero sparse: 2-by-2
        \enddoctest
        \endlang

        \sa print_dense

        \identifier{2p3} */
    void print_sparse(std::ostream &stream, bool truncate=true) const;

#ifndef SWIG
    /// Print default style
    static void print_default(std::ostream &stream, const Sparsity& sp, const Scalar* nonzeros,
      bool truncate=true);

    /// Print scalar
    static void print_scalar(std::ostream &stream, const Scalar& e);

    /// Print vector-style
    static void print_vector(std::ostream &stream, const Sparsity& sp, const Scalar* nonzeros,
      bool truncate=true);

    /// Print canonical style
    static void print_canonical(std::ostream &stream, const Sparsity& sp, const Scalar* nonzeros,
        bool truncate=true);

    /// Print scalar
    static void print_sparse(std::ostream &stream, const Sparsity& sp, const Scalar* nonzeros,
      bool truncate=true);

    /// Get strings corresponding to the nonzeros and the interdependencies
    static void print_split(casadi_int nnz, const Scalar* nonzeros, std::vector<std::string>& nz,
                    std::vector<std::string>& inter);

    /// Print dense matrix-stype
    static void print_dense(std::ostream &stream,  const Sparsity& sp, const Scalar* nonzeros,
      bool truncate=true);
#endif

    /// Reset to an empty 0-by-0 matrix, in place
    void clear();

    /** \brief Change the shape in place, keeping nonzeros that remain inside the new bounds

        Growing appends structural zeros; shrinking drops the entries outside the new shape.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> A.resize(3, 3)
        >>> A
        DM(
        [[1, 2, 00],
         [3, 4, 00],
         [00, 00, 00]])
        >>> A.resize(1, 1)
        >>> A
        DM(1)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> A.resize(3, 3);
        >> A
        [[1, 2, 00],
         [3, 4, 00],
         [00, 00, 00]]
        >> A.resize(1, 1);
        >> A
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> ca.resize(A, 3, 3)
        julia> A
        DM(
        [[1, 2, 00],
         [3, 4, 00],
         [00, 00, 00]])
        julia> ca.resize(A, 1, 1)
        julia> A
        DM(1)
        \enddoctest
        \endlang

        \sa enlarge, reshape

        \identifier{2p4} */
    void resize(casadi_int nrow, casadi_int ncol);

    /// Reserve storage for the given number of nonzeros (performance hint, no visible effect)
    void reserve(casadi_int nnz);
    void reserve(casadi_int nnz, casadi_int ncol);

    /** \brief Erase the block at the given rows and columns, leaving structural zeros

        The shape is unchanged; the entries of the submatrix rr x cc are removed from the
        sparsity pattern. Negative indices count from the end; out-of-range indices raise
        an error.

        \param rr Row indices
        \param cc Column indices
        \param ind1 True for one-based indices

        \lang{python}
        \doctest
        >>> A = DM([[1, 2, 3], [4, 5, 6]])
        >>> A.erase([0], [1])
        >>> A
        DM(
        [[1, 00, 3],
         [4, 5, 6]])
        >>> A.dim(), A.nnz()
        ('2x3', 5)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2 3; 4 5 6]);
        >> A.erase(1, 2, true);
        >> A
        [[1, 00, 3],
         [4, 5, 6]]
        >> fprintf('%s %d\n', A.dim(), A.nnz())
        2x3 5
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2 3; 4 5 6]);
        julia> ca.erase(A, [0], [1])
        julia> A
        DM(
        [[1, 00, 3],
         [4, 5, 6]])
        julia> ca.dim(A), ca.nnz(A)
        ("2x3", 5)
        \enddoctest
        \endlang

        \sa remove, sparsify

        \identifier{19f} */
    void erase(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
               bool ind1=false);

    /** \brief Erase the elements at the given linear indices, leaving structural zeros

        Column-major linear indices, as in get. The shape is unchanged.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2, 3], [4, 5, 6]])
        >>> A.erase([0, -1])
        >>> A
        DM(
        [[00, 2, 3],
         [4, 5, 00]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2 3; 4 5 6]);
        >> A.erase([1 6], true);
        >> A
        [[00, 2, 3],
         [4, 5, 00]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2 3; 4 5 6]);
        julia> ca.erase(A, [0, 5])
        julia> A
        DM(
        [[00, 2, 3],
         [4, 5, 00]])
        \enddoctest
        \endlang

        \identifier{19g} */
    void erase(const std::vector<casadi_int>& rr, bool ind1=false);

    /** \brief Delete rows and columns, shrinking the matrix in place

        Both lists may be empty; duplicates are ignored, out-of-range indices raise an error.

        \lang{python}
        \doctest
        >>> A = DM([[1, 2, 3], [4, 5, 6]])
        >>> A.remove([0], [1])
        >>> A
        DM([[4, 6]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2 3; 4 5 6]);
        >> A.remove(0, 1);
        >> A
        [[4, 6]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2 3; 4 5 6]);
        julia> ca.remove(A, [0], [1])
        julia> A
        DM([[4, 6]])
        \enddoctest
        \endlang

        \sa erase

        \identifier{19h} */
    void remove(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc);

    /** \brief Enlarge the matrix in place, moving the existing rows and columns to new positions

        The new shape is nrow-by-ncol. \p rr lists, for every existing row, its index in the
        enlarged matrix; \p cc does the same for the columns. Both must have exactly size1()
        and size2() entries. All other rows and columns are structurally zero.

        \param nrow New number of rows
        \param ncol New number of columns
        \param rr New position of each existing row, length size1()
        \param cc New position of each existing column, length size2()
        \param ind1 True for one-based indices

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> A.enlarge(3, 4, [0, 2], [0, 1])
        >>> A
        DM(
        [[1, 2, 00, 00],
         [00, 00, 00, 00],
         [3, 4, 00, 00]])
        >>> try: A.enlarge(4, 4, [0, 1], [0, 1])
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> A.enlarge(3, 4, [1 3], [1 2], true);
        >> A
        [[1, 2, 00, 00],
         [00, 00, 00, 00],
         [3, 4, 00, 00]]
        >> try, A.enlarge(4, 4, [1 2], [1 2], true); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> ca.enlarge(A, 3, 4, [0, 2], [0, 1])
        julia> A
        DM(
        [[1, 2, 00, 00],
         [00, 00, 00, 00],
         [3, 4, 00, 00]])
        julia> try ca.enlarge(A, 4, 4, [0, 1], [0, 1]); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa resize, Sparsity::enlarge

        \identifier{19i} */
    void enlarge(casadi_int nrow, casadi_int ncol,
                  const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                  bool ind1=false);

#ifndef SWIG
    ///@{
    /// Access the non-zero elements
    std::vector<Scalar>& nonzeros();
    const std::vector<Scalar>& nonzeros() const;
    ///@}

    ///@{
    /// Get a pointer to the data
    Scalar* ptr();
    const Scalar* ptr() const;
    friend inline Scalar* get_ptr(Matrix<Scalar>& v) { return v.ptr(); }
    friend inline const Scalar* get_ptr(const Matrix<Scalar>& v) { return v.ptr(); }
    ///@}

    /// Const access the sparsity - reference to data member
    const Sparsity& sparsity() const;

#endif // SWIG

    /** \brief Get the sparsity pattern (Python: sparsity())

        Returns a copy of the Sparsity object; patterns are reference counted, so this is cheap.

        \lang{python}
        \doctest
        >>> DM(Sparsity.diag(2)).sparsity()
        Sparsity(2x2,2nz)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.diag(2)).sparsity()
        2x2,2nz
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.sparsity(ca.DM(ca.Sparsity.diag(2)))
        Sparsity(2x2,2nz)
        \enddoctest
        \endlang

        \sa nonzeros

        \identifier{19j} */
    Sparsity get_sparsity() const;

    /** \brief  Get an output

        \identifier{28i} */
    Matrix<Scalar> get_output(casadi_int oind) const;

    /** \brief Construct a sparse matrix from triplet form (row indices, column indices, values)

        Entry k of the result is at (row[k], col[k]) with value d[k]. The three lists must have
        the same length; indices are zero-based and must be within the matrix bounds.
        The nonzeros are reordered into column-major order. When the same (row, col) pair
        appears more than once, only the first value is kept; duplicates are not summed.

        \param row Row index of each entry
        \param col Column index of each entry
        \param d Dense vector of values, one per entry
        \param nrow Number of rows of the result [default: largest row index plus one]
        \param ncol Number of columns of the result [default: largest column index plus one]

        \lang{python}
        \doctest
        >>> DM.triplet([0, 1], [1, 0], [5, 6], 3, 3)
        DM(
        [[00, 5, 00],
         [6, 00, 00],
         [00, 00, 00]])
        >>> DM.triplet([0, 1], [1, 0], [5, 6], (3, 3)).nonzeros()
        [6.0, 5.0]
        >>> SX.triplet([0, 1], [0, 1], SX.sym("v", 2), 2, 2)
        SX(
        [[v_0, 00],
         [00, v_1]])
        >>> DM.triplet([], [], [], 2, 2).nnz()
        0
        >>> try: DM.triplet([0, 1], [1, 0], [5], 3, 3)
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM.triplet([0 1], [1 0], [5 6], 3, 3)
        [[00, 5, 00],
         [6, 00, 00],
         [00, 00, 00]]
        >> DM.triplet([0 1], [1 0], [5 6], [3 3]).nonzeros()
        6 5
        >> SX.triplet([0 1], [0 1], SX.sym('v', 2), 2, 2)
        [[v_0, 00],
         [00, v_1]]
        >> DM.triplet([], [], [], 2, 2).nnz()
        0
        >> try, DM.triplet([0 1], [1 0], 5, 3, 3); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM.triplet([0, 1], [1, 0], [5, 6], 3, 3)
        DM(
        [[00, 5, 00],
         [6, 00, 00],
         [00, 00, 00]])
        julia> ca.SX.triplet([0, 1], [0, 1], ca.SX.sym("v", 2), 2, 2)
        SX(
        [[v_0, 00],
         [00, v_1]])
        julia> ca.nnz(ca.DM.triplet([], [], [], 2, 2))
        0
        julia> try ca.DM.triplet([0, 1], [1, 0], [5], 3, 3); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa Sparsity::triplet

        \identifier{23t} */
    ///@{
    static Matrix<Scalar> triplet(const std::vector<casadi_int>& row,
                                  const std::vector<casadi_int>& col,
                                  const Matrix<Scalar>& d);
    static Matrix<Scalar> triplet(const std::vector<casadi_int>& row,
                                  const std::vector<casadi_int>& col,
                                  const Matrix<Scalar>& d, casadi_int nrow, casadi_int ncol);
    static Matrix<Scalar> triplet(const std::vector<casadi_int>& row,
                                  const std::vector<casadi_int>& col,
                                  const Matrix<Scalar>& d,
                                  const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    ///@{
    /** \brief Create a matrix with every nonzero equal to +inf

        Takes a shape (dense result, default 1-by-1) or a sparsity pattern.

        \lang{python}
        \doctest
        >>> DM.inf(), DM.inf(2, 1).is_inf()
        (DM(inf), True)
        >>> DM.inf(Sparsity.diag(2))
        DM(
        [[inf, 00],
         [00, inf]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM.inf(), DM.inf(2, 1).is_inf()
        inf
        1
        >> DM.inf(Sparsity.diag(2))
        [[inf, 00],
         [00, inf]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM.inf(), ca.is_inf(ca.DM.inf(2, 1))
        (inf, true)
        julia> ca.DM.inf(ca.Sparsity.diag(2))
        DM(
        [[inf, 00],
         [00, inf]])
        \enddoctest
        \endlang

        \sa nan, is_inf

        \identifier{19k} */
    static Matrix<Scalar> inf(const Sparsity& sp);
    static Matrix<Scalar> inf(casadi_int nrow=1, casadi_int ncol=1);
    static Matrix<Scalar> inf(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    ///@{
    /** \brief Create a matrix with every nonzero equal to nan

        Takes a shape (dense result, default 1-by-1) or a sparsity pattern.
        Useful as a placeholder that is easy to detect with is_regular.

        \lang{python}
        \doctest
        >>> DM.nan(), DM.nan((2, 1)).dim()
        (DM(nan), '2x1')
        >>> DM.nan(2, 2).is_regular()
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM.nan(), DM.nan([2 1]).dim()
        nan
        2x1
        >> DM.nan(2, 2).is_regular()
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM.nan(), ca.dim(ca.DM.nan(2, 1))
        (nan, "2x1")
        julia> ca.is_regular(ca.DM.nan(2, 2))
        false
        \enddoctest
        \endlang

        \sa inf, is_regular

        \identifier{19l} */
    static Matrix<Scalar> nan(const Sparsity& sp);
    static Matrix<Scalar> nan(casadi_int nrow=1, casadi_int ncol=1);
    static Matrix<Scalar> nan(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    /** \brief Create an n-by-n identity matrix

        The result is sparse: only the n diagonal entries are structural nonzeros.

        \lang{python}
        \doctest
        >>> DM.eye(2)
        DM(
        [[1, 00],
         [00, 1]])
        >>> DM.eye(3).nnz(), DM.eye(0).dim()
        (3, '0x0')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM.eye(2)
        [[1, 00],
         [00, 1]]
        >> fprintf('%d %s\n', DM.eye(3).nnz(), DM.eye(0).dim())
        3 0x0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM.eye(2)
        DM(
        [[1, 00],
         [00, 1]])
        julia> ca.nnz(ca.DM.eye(3)), ca.dim(ca.DM.eye(0))
        (3, "0x0")
        \enddoctest
        \endlang

        \sa is_eye, Sparsity::diag

        \identifier{19m} */
    static Matrix<Scalar> eye(casadi_int n);

    /** \brief Hash of the node of a 1-by-1 SX expression

        Equal expressions (same node) have the same hash. Only defined for a scalar SX;
        a vector or a DM raises an error.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> x.element_hash() == x.element_hash(), x.element_hash() == y.element_hash()
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> h = x.element_hash();
        >> fprintf('%d %d\n', h == x.element_hash(), h == y.element_hash())
        1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.element_hash(x) == ca.element_hash(x), ca.element_hash(x) == ca.element_hash(y)
        (true, false)
        \enddoctest
        \endlang

        \identifier{19n} */
    casadi_int element_hash() const;

    /** \brief Check that no entry is nan, inf or -inf

        Structural zeros are ignored. For SX, every entry must be a constant: the check
        raises an error when it reaches a symbolic entry (an infinite constant found before
        it still gives false).

        \lang{python}
        \doctest
        >>> DM([1, 2]).is_regular(), DM([1, inf]).is_regular(), DM.nan(2, 2).is_regular()
        (True, False, False)
        >>> DM(2, 2).is_regular(), SX([1, 2]).is_regular()
        (True, True)
        >>> try: SX.sym("x").is_regular()
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%d %d\n', DM([1; 2]).is_regular(), DM([1; inf]).is_regular())
        1 0
        >> fprintf('%d %d\n', DM.nan(2, 2).is_regular(), DM(2, 2).is_regular())
        0 1
        >> SX([1; 2]).is_regular()
        1
        >> try, SX.sym('x').is_regular(); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_regular(ca.DM([1, 2])), ca.is_regular(ca.DM([1, Inf]))
        (true, false)
        julia> ca.is_regular(ca.DM.nan(2, 2))
        false
        julia> ca.is_regular(ca.DM(2, 2)), ca.is_regular(ca.SX([1, 2]))
        (true, true)
        julia> try ca.is_regular(ca.SX.sym("x")); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa nan, inf

        \identifier{2p5} */
    bool is_regular() const;

    /** \brief Check whether an SX expression contains only smooth operations

        False when any node is a non-smooth operation such as floor, ceil, sign, a comparison
        or if_else. Operations such as fabs, fmin and fmax count as smooth. Always true for DM.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> sin(x).is_smooth(), fabs(x).is_smooth(), floor(x).is_smooth()
        (True, True, False)
        >>> if_else(x > 0, x, 0).is_smooth(), (x > 0).is_smooth(), DM(1).is_smooth()
        (False, False, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> fprintf('%d %d %d\n', sin(x).is_smooth(), abs(x).is_smooth(), floor(x).is_smooth())
        1 1 0
        >> c = x > 0; b = if_else(c, x, 0);
        >> fprintf('%d %d %d\n', b.is_smooth(), c.is_smooth(), DM(1).is_smooth())
        0 0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> ca.is_smooth(sin(x)), ca.is_smooth(abs(x)), ca.is_smooth(floor(x))
        (true, true, false)
        julia> ca.is_smooth(ifelse(x > 0, x, 0)), ca.is_smooth(x > 0), ca.is_smooth(ca.DM(1))
        (false, false, true)
        \enddoctest
        \endlang

        \identifier{19o} */
    bool is_smooth() const;

    /** \brief Check if SX is a leaf of the SX graph

        Only defined if symbolic scalar.

        \identifier{19p} */
    bool is_leaf() const;

    /** \brief Check whether a binary SX is commutative

        Only defined if symbolic scalar.

        \identifier{19q} */
    bool is_commutative() const;

    /** \brief Check if symbolic (Dense)

        Check if an expression is a pure symbol.
        Pure means that no operations should happen on the symbol
        (not even vec, transpose, index, concatenation, ...)
        By consequence, a slice of a vector-shaped MX symbol is not a symbol.
        However, the SX type is really a container format with scalar entries.
        Therefore, a slice of a vector-shaped SX symbol is still a symbol.

        Sparse matrices invariable return false

        \seealso is_valid_input

        \identifier{19r} */
    bool is_symbolic() const;

    /** \brief Check if matrix can be used to define function inputs.


        is_valid_input is more forgiving than is_symbolic.
        Some compositions are allowed: vec, vertcat.

        Sparse matrices can return true if all non-zero elements are symbolic

        \seealso is_symbolic

        \identifier{19s} */
    bool is_valid_input() const;

    /// \cond INTERNAL
    /** \brief Detect symbolic primitives that appear more than once in an SX matrix

        Returns true and emits a warning naming each duplicate. Used when checking Function
        inputs ("input expressions are not independent"). Marks the visited nodes; a second
        call on the same symbols reports them as duplicates until reset_input is called.
        Only defined for SX.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> v = vertcat(x, y)
        >>> v.has_duplicates()
        False
        >>> v.reset_input()
        >>> vertcat(x, x).has_duplicates()
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> v = vertcat(x, y);
        >> v.has_duplicates()
        0
        >> v.reset_input();
        >> vertcat(x, x).has_duplicates()
        ...
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> v = [x; y];
        julia> ca.has_duplicates(v)
        false
        julia> ca.reset_input(v)
        julia> ca.has_duplicates([x; x])
        true
        \enddoctest
        \endlang

        \sa reset_input, is_valid_input

        \identifier{19t} */
    bool has_duplicates() const;

    /** \brief Clear the markers set by has_duplicates

        \identifier{19u} */
    void reset_input() const;
  /// \endcond

    /** \brief Check if every entry is a numerical constant (no symbols)

        Always true for DM. For SX, true when no entry depends on a symbol; the
        simplifications make x*0 and x-x constant. A false negative is possible for
        expressions that are constant only in value, e.g. x - x after other operations.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> SX([1, 2]).is_constant(), x.is_constant(), (x * 0).is_constant()
        (True, False, True)
        >>> vertcat(x, 1).is_constant(), SX.sym("e", 2, 0).is_constant()
        (False, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); z = x * 0;
        >> fprintf('%d %d %d\n', SX([1; 2]).is_constant(), x.is_constant(), z.is_constant())
        1 0 1
        >> fprintf('%d %d\n', vertcat(x, 1).is_constant(), SX.sym('e', 2, 0).is_constant())
        0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> ca.is_constant(ca.SX([1, 2])), ca.is_constant(x), ca.is_constant(x * 0)
        (true, false, true)
        julia> ca.is_constant([x; 1]), ca.is_constant(ca.SX.sym("e", 2, 0))
        (false, true)
        \enddoctest
        \endlang

        \sa is_symbolic, evalf

        \identifier{19v} */
    bool is_constant() const;

    /** \brief Check if function call

        \identifier{28j} */
    bool is_call() const;

    /** \brief  Check if evaluation output

        \identifier{28k} */
    bool is_output() const;

    /** \brief  Check if a multiple output node

        \identifier{28l} */
    bool has_output() const;

    /** \brief Get the index of evaluation output - only valid when is_output() is true

        \identifier{28m} */
    casadi_int which_output() const;

    /** \brief Get function - only valid when is_call() is true

        \identifier{28n} */
    Function which_function() const;

    /** \brief Check if the matrix is integer-valued

     * (note that false negative answers are possible)

        \identifier{19w} */
    bool is_integer() const;

    /** \brief  check if the matrix is 0 (note that false negative answers are possible)

        \identifier{19x} */
    bool is_zero() const;

    /** \brief  check if the matrix is 1 (note that false negative answers are possible)

        \identifier{19y} */
    bool is_one() const;

    /** \brief  check if the matrix is -1 (note that false negative answers are possible)

        \identifier{19z} */
    bool is_minus_one() const;

    /** \brief  check if the matrix is 0.5 (note that false negative answers are possible)

        \identifier{2en} */
    bool is_half() const;

    /** \brief  check if the matrix is value

        \identifier{2eo} */
    bool is_value(double val) const;

    /** \brief  check if the matrix is inf (note that false negative answers are possible)

        \identifier{2e9} */
    bool is_inf() const;

    /** \brief  check if the matrix is -inf (note that false negative answers are possible)

        \identifier{2ea} */
    bool is_minus_inf() const;

    /** \brief  check if the matrix is not negative 

     (note that false negative answers are possible)

        \identifier{2eb} */
    bool is_nonnegative() const;

    /** \brief  check if the matrix is an identity matrix (note that false negative answers

     * are possible)

        \identifier{1a0} */
    bool is_eye() const;

    /// Get operation type
    casadi_int op() const;

    /// Is it a certain operation
    bool is_op(casadi_int op) const;

    /** \brief Check if any structural nonzero is known to be numerically zero

        For SX only entries that are the constant 0 count. Use sparsify to remove them.

        \lang{python}
        \doctest
        >>> DM([1, 0]).has_zeros(), DM([1, 2]).has_zeros(), DM(2, 2).has_zeros()
        (True, False, False)
        >>> x = SX.sym("x")
        >>> vertcat(x, 0).has_zeros(), vertcat(x, x * 0).has_zeros()
        (True, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%d %d\n', DM([1; 0]).has_zeros(), DM([1; 2]).has_zeros())
        1 0
        >> DM(2, 2).has_zeros()
        0
        >> x = SX.sym('x');
        >> fprintf('%d %d\n', vertcat(x, 0).has_zeros(), vertcat(x, x * 0).has_zeros())
        1 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.has_zeros(ca.DM([1, 0])), ca.has_zeros(ca.DM([1, 2])), ca.has_zeros(ca.DM(2, 2))
        (true, false, false)
        julia> x = ca.SX.sym("x");
        julia> ca.has_zeros([x; 0]), ca.has_zeros([x; x * 0])
        (true, true)
        \enddoctest
        \endlang

        \sa sparsify, has_nz

        \identifier{1a1} */
    bool has_zeros() const;

    /** \brief Get all structural nonzeros as a flat list in column-major order (Python: nonzeros())

        Returns floats for DM and 1-by-1 SX objects for SX. Structural zeros are skipped.

        \lang{python}
        \doctest
        >>> DM([[1, 2], [3, 4]]).nonzeros()
        [1.0, 3.0, 2.0, 4.0]
        >>> DM(Sparsity.diag(2)).nonzeros(), DM(2, 0).nonzeros()
        ([1.0, 1.0], [])
        >>> SX.sym("x", 2).nonzeros()
        [SX(x_0), SX(x_1)]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM([1 2; 3 4]).nonzeros()
        1 3 2 4
        >> DM(Sparsity.diag(2)).nonzeros(), size(DM(2, 0).nonzeros())
        1 1
        1 0
        >> v = SX.sym('x', 2).nonzeros();
        >> class(v), v{2}
        cell
        x_1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.nonzeros(ca.DM([1 2; 3 4]))
        4-element Vector{Float64}:
         1.0
         3.0
         2.0
         4.0
        julia> ca.nonzeros(ca.DM(ca.Sparsity.diag(2))), ca.nonzeros(ca.DM(2, 0))
        ([1.0, 1.0], Float64[])
        \enddoctest
        \endlang

        \sa elements, get_nz

        \identifier{1a2} */
    std::vector<Scalar> get_nonzeros() const;

    /** \brief Get all elements as a flat list in column-major order (Python: elements())

        Like nonzeros, but structural zeros are included as 0, so the list has numel() entries.

        \lang{python}
        \doctest
        >>> DM(Sparsity.diag(2)).elements()
        [1.0, 0.0, 0.0, 1.0]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.diag(2)).elements()
        1 0 0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.elements(ca.DM(ca.Sparsity.diag(2)))
        4-element Vector{Float64}:
         1.0
         0.0
         0.0
         1.0
        \enddoctest
        \endlang

        \sa nonzeros, densify

        \identifier{1a3} */
    std::vector<Scalar> get_elements() const;

#ifndef SWIG
    /** \brief Get all nonzeros

        \identifier{1a4} */
    template<typename A>
    std::vector<A> get_nonzeros() const;
#endif // SWIG

    /** \brief Conversion of a 1-by-1 matrix to a double (Python: float(x))

        Other shapes raise an error. A symbolic SX entry converts to nan without an error.

        \lang{python}
        \doctest
        >>> float(DM(2.5)), float(SX(3))
        (2.5, 3.0)
        >>> try: float(DM([1, 2]))
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%g %g\n', full(DM(2.5)), SX(3).to_double())
        2.5 3
        >> try, DM([1; 2]).to_double(); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> Float64(ca.DM(2.5)), Float64(ca.SX(3))
        (2.5, 3.0)
        julia> try Float64(ca.DM([1, 2])); catch; println("error") end
        error
        \enddoctest
        \endlang

        \identifier{1a5} */
    explicit operator double() const;

    /** \brief Conversion of a 1-by-1 matrix to an integer (Python: int(x)), truncating

        \lang{python}
        \doctest
        >>> int(DM(3.7))
        3
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(3.7).to_int()
        3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.to_int(ca.DM(3.7))
        3
        \enddoctest
        \endlang

        \identifier{1a6} */
    explicit operator casadi_int() const;

#ifndef SWIG
    /** \brief Type conversion to a vector

        \identifier{1a7} */
    template<typename A>
    explicit operator std::vector<A>() const;
#endif // SWIG

    /** \brief Name of a 1-by-1 symbolic SX primitive

        Raises an error for a non-scalar, for an expression that is not a pure symbol, and for DM.

        \lang{python}
        \doctest
        >>> SX.sym("x").name(), SX.sym("v", 2)[1].name()
        ('x', 'v_1')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = SX.sym('v', 2); v2 = v(2);
        >> fprintf('%s %s\n', SX.sym('x').name(), v2.name())
        x v_1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.name(ca.SX.sym("x")), ca.name(ca.SX.sym("v", 2)[2])
        ("x", "v_1")
        \enddoctest
        \endlang

        \sa is_symbolic

        \identifier{1a8} */
    std::string name() const;

    /** \brief Get operand ch of a 1-by-1 SX expression

        Only defined for a scalar SX that is the result of an operation; a symbol, a
        constant or a DM raises an error.

        \param ch Operand index, 0 or 1

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> e = x * y
        >>> e.dep(0), e.dep(1), e.n_dep()
        (SX(x), SX(y), 2)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> e = x * y;
        >> e.dep(0), e.dep(1), e.n_dep()
        x
        y
        2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> e = x * y;
        julia> ca.dep(e, 0), ca.dep(e, 1), ca.n_dep(e)
        (x, y, 2)
        \enddoctest
        \endlang

        \sa n_dep, op

        \identifier{1a9} */
    Matrix<Scalar> dep(casadi_int ch=0) const;

    /** \brief Number of operands of a 1-by-1 SX expression (0 for symbols and constants)

        \sa dep

        \identifier{1aa} */
    casadi_int n_dep() const;

    // @{
    /** \brief Set the precision, field width and scientific flag used when printing nonzeros

        Global settings per class (DM and SX each keep their own), applied by str(), print()
        and the print_* methods. Defaults: precision 6, width 0, scientific false.

        \lang{python}
        \doctest
        >>> DM.set_precision(3)
        >>> DM(1/3)
        DM(0.333)
        >>> DM.set_precision(6)
        >>> DM(1/3)
        DM(0.333333)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM.set_precision(3);
        >> DM(1/3)
        0.333
        >> DM.set_precision(6);
        >> DM(1/3)
        0.333333
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM.set_precision(3)
        julia> ca.DM(1/3)
        DM(0.333)
        julia> ca.DM.set_precision(6)
        julia> ca.DM(1/3)
        DM(0.333333)
        \enddoctest
        \endlang

        \sa get_precision

        \identifier{2p6} */
    static void set_precision(casadi_int precision);
    static void set_width(casadi_int width);
    static void set_scientific(bool scientific);
    // @}

    // @{
    /// Get the precision, field width and scientific flag used when printing nonzeros
    static casadi_int get_precision();
    static casadi_int get_width();
    static bool get_scientific();
    // @}

    /** \brief Seed the random number generator used by rand

        The generator is global to the class and is seeded from the clock at startup;
        seeding it makes subsequent rand calls reproducible.

        \lang{python}
        \doctest
        >>> DM.rng(1); a = DM.rand(2, 2)
        >>> DM.rng(1); b = DM.rand(2, 2)
        >>> norm_inf(a - b)
        DM(0)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM.rng(1); a = DM.rand(2, 2);
        >> DM.rng(1); b = DM.rand(2, 2);
        >> norm_inf(a - b)
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM.rng(1); a = ca.DM.rand(2, 2);
        julia> ca.DM.rng(1); b = ca.DM.rand(2, 2);
        julia> ca.norm_inf(a - b)
        DM(0)
        \enddoctest
        \endlang

        \sa rand

        \identifier{2p7} */
    static void rng(casadi_int seed);

    ///@{
    /** \brief Create a matrix with uniformly distributed random numbers in [0, 1)

        Takes a shape (dense result, default 1-by-1) or a sparsity pattern, in which case
        only the structural nonzeros get a random value. Only available for DM.

        \lang{python}
        \doctest
        >>> DM.rand(2, 3).dim(), DM.rand().dim(), DM.rand(Sparsity.diag(3)).nnz()
        ('2x3', '1x1', 3)
        >>> try: SX.rand(2)
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%s %s\n', DM.rand(2, 3).dim(), DM.rand().dim())
        2x3 1x1
        >> DM.rand(Sparsity.diag(3)).nnz()
        3
        >> try, SX.rand(2); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.dim(ca.DM.rand(2, 3)), ca.dim(ca.DM.rand())
        ("2x3", "1x1")
        julia> ca.nnz(ca.DM.rand(ca.Sparsity.diag(3)))
        3
        julia> try ca.SX.rand(2); catch; println("error") end
        error
        \enddoctest
        \endlang

        \sa rng

        \identifier{1ab} */
    static Matrix<Scalar> rand(  // NOLINT(runtime/threadsafe_fn)
        casadi_int nrow=1,
        casadi_int ncol=1);
    static Matrix<Scalar> rand(const Sparsity& sp);  // NOLINT(runtime/threadsafe_fn)
    static Matrix<Scalar> rand(  // NOLINT(runtime/threadsafe_fn)
        const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    /** \brief Export matrix in specific language
     *
     * lang: only 'matlab' supported for now
     * \verbatim
     * options:
     *   inline: Indicates if you want everything on a single line (default: False)
     *   name: Name of exported variable (default: 'm')
     *   indent_level: Level of indentation (default: 0)
     *   spoof_zero: Replace numerical zero by a 1e-200 (default: false)
     *               might be needed for matlab sparse construct,
     *               which doesn't allow numerical zero
     * \endverbatim

        \identifier{1ac} */
    void export_code(const std::string& lang,
        std::ostream &stream=casadi::uout(), const Dict& options=Dict()) const;

    /** \brief Dictionary with the sparsity pattern and the nonzeros (DM only)

        \lang{python}
        \doctest
        >>> sorted(DM([1, 2]).info().keys())
        ['data', 'sparsity']
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> s = DM([1; 2]).info();
        >> strjoin(sort(fieldnames(s))', ' ')
        data sparsity
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sort(collect(keys(ca.info(ca.DM([1, 2])))))
        2-element Vector{String}:
         "data"
         "sparsity"
        \enddoctest
        \endlang

        \identifier{2p8} */
    Dict info() const;
    #ifndef SWIG
        /** \brief Serialize an object

            \identifier{1ad} */
        void serialize(std::ostream &stream) const;
    #endif

    /** \brief Serialize to a string

        The string can be restored with deserialize of the same class (DM or SX).
        For SX the whole expression graph is stored, symbols included.

        \lang{python}
        \doctest
        >>> s = DM([[1, 2], [3, 4]]).serialize()
        >>> DM.deserialize(s)
        DM(
        [[1, 2],
         [3, 4]])
        >>> x = SX.sym("x")
        >>> SX.deserialize((sin(x) + x).serialize())
        SX((sin(x)+x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> s = DM([1 2; 3 4]).serialize();
        >> DM.deserialize(s)
        [[1, 2],
         [3, 4]]
        >> x = SX.sym('x'); e = sin(x) + x;
        >> SX.deserialize(e.serialize())
        (sin(x)+x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> s = ca.serialize(ca.DM([1 2; 3 4]));
        julia> ca.DM.deserialize(s)
        DM(
        [[1, 2],
         [3, 4]])
        julia> x = ca.SX.sym("x");
        julia> ca.SX.deserialize(ca.serialize(sin(x) + x))
        SX((sin(x)+x))
        \enddoctest
        \endlang

        \sa deserialize, StringSerializer

        \identifier{1ae} */
    std::string serialize() const;

    /** \brief Build a matrix from a serialization stream

        \identifier{1af} */
    static Matrix<Scalar> deserialize(std::istream& stream);

    /** \brief Build a matrix from a string produced by serialize

        A string that is not a CasADi serialization raises an error.

        \sa serialize

        \identifier{1ag} */
    static Matrix<Scalar> deserialize(const std::string& s);

    /** \brief Serialize an object

        \identifier{1ah} */
    void serialize(SerializingStream& s) const;

    static Matrix<Scalar> deserialize(DeserializingStream& s);

    // @{
    /** \brief Write a numerical matrix to a text file (DM only)

        Supported formats, selected by the file extension unless \p format is given:
        - mtx: Matrix Market coordinate format (sparse, one-based row/col/value triplets)
        - txt: one line per row, whitespace separated, aligned columns, full precision,
               structural zeros written as 00. Does not scale to large sparse matrices.
        Values are written with the C locale; nan and inf round-trip. An unknown extension
        or format raises an error. SX is not supported.

        \param filename Output path
        \param format "mtx" or "txt"; empty means deduce from the extension

        \lang{python}
        \doctest
        >>> import os, tempfile
        >>> fn = os.path.join(tempfile.mkdtemp(), "a.txt")
        >>> A = DM(Sparsity.diag(2)) * 2.5
        >>> A.to_file(fn)
        >>> DM.from_file(fn)
        DM(
        [[2.5, 00],
         [00, 2.5]])
        >>> os.remove(fn); os.rmdir(os.path.dirname(fn))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fn = [tempname '.txt'];
        >> A = DM(Sparsity.diag(2)) * 2.5;
        >> A.to_file(fn);
        >> DM.from_file(fn)
        [[2.5, 00],
         [00, 2.5]]
        >> delete(fn);
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> fn = joinpath(mktempdir(), "a.txt");
        julia> A = ca.DM(ca.Sparsity.diag(2)) * 2.5;
        julia> ca.to_file(A, fn)
        julia> ca.DM.from_file(fn)
        DM(
        [[2.5, 00],
         [00, 2.5]])
        julia> rm(dirname(fn), recursive=true)
        \enddoctest
        \endlang

        \sa from_file, Sparsity::to_file

        \identifier{2p9} */
    void to_file(const std::string& filename, const std::string& format="") const;
#ifndef SWIG
    static void to_file(const std::string& filename, const Sparsity& sp,
      const Scalar* nonzeros, const std::string& format="");
#endif

    /** \brief Read a numerical matrix from a text file written by to_file (DM only)

        Same formats as to_file. In txt format, lines starting with #, % or / and empty lines
        are skipped, 00 marks a structural zero and a plain 0 a stored zero; every row must
        have the same number of columns, otherwise an error is raised. An empty file gives a
        0-by-0 matrix. A missing file raises an error.

        \param filename Input path
        \param format_hint "mtx" or "txt"; empty means deduce from the extension

        \lang{python}
        \doctest
        >>> import os, tempfile
        >>> fn = os.path.join(tempfile.mkdtemp(), "a.txt")
        >>> with open(fn, "w") as f: n = f.write("# comment\n1 00\n0 2\n")
        >>> A = DM.from_file(fn)
        >>> A
        DM(
        [[1, 00],
         [0, 2]])
        >>> A.nnz()
        3
        >>> os.remove(fn); os.rmdir(os.path.dirname(fn))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fn = [tempname '.txt'];
        >> fid = fopen(fn, 'w'); fprintf(fid, '# comment\n1 00\n0 2\n'); fclose(fid);
        >> A = DM.from_file(fn);
        >> A
        [[1, 00],
         [0, 2]]
        >> A.nnz()
        3
        >> delete(fn);
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> fn = joinpath(mktempdir(), "a.txt");
        julia> write(fn, "# comment\n1 00\n0 2\n");
        julia> A = ca.DM.from_file(fn)
        DM(
        [[1, 00],
         [0, 2]])
        julia> ca.nnz(A)
        3
        julia> rm(dirname(fn), recursive=true)
        \enddoctest
        \endlang

        \sa to_file

        \identifier{2pa} */
    static Matrix<double> from_file(const std::string& filename, const std::string& format_hint="");
    //@}

#ifndef SWIG
    /// Sparse matrix with a given sparsity with all values same
    Matrix(const Sparsity& sp, const Scalar& val, bool dummy);

    /// Sparse matrix with a given sparsity and non-zero elements.
    Matrix(const Sparsity& sp, const std::vector<Scalar>& d, bool dummy);

    /// \cond INTERNAL
    static Matrix<Scalar> _sym(const std::string& name, const Sparsity& sp);
    /// \endcond

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex& get_mutex_temp();
#endif //CASADI_WITH_THREADSAFE_SYMBOLICS

  private:
    /// Sparsity of the matrix in a compressed column storage (CCS) format
    Sparsity sparsity_;

    /// Nonzero elements
    std::vector<Scalar> nonzeros_;

    /// Precision used in streams
    static casadi_int stream_precision_;
    static casadi_int stream_width_;
    static bool stream_scientific_;

    /// Random number generator
    static std::default_random_engine rng_;

#endif // SWIG
  };

  /// Implementation of Matrix::get_nonzeros (in public API)
  template<typename Scalar>
  template<typename A>
  std::vector<A> Matrix<Scalar>::get_nonzeros() const {
    std::vector<A> ret(nnz());
    auto r = ret.begin();
    for (auto&& e : nonzeros()) *r++ = static_cast<A>(e);
    return ret;
  }

  /// Implementation of std::vector(Matrix) (in public API)
  template<typename Scalar>
  template<typename A>
  Matrix<Scalar>::operator std::vector<A>() const {
    // Get sparsity pattern
    casadi_int size1 = this->size1(), size2 = this->size2();
    const casadi_int *colind = this->colind(), *row = this->row();
    // Copy the nonzeros
    auto it = nonzeros().begin();
    std::vector<A> ret(numel(), 0);
    for (casadi_int cc=0; cc<size2; ++cc) {
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        ret[row[el] + cc*size1] = static_cast<A>(*it++);
      }
    }
    return ret;
  }

} // namespace casadi

#endif // CASADI_MATRIX_DECL_HPP
