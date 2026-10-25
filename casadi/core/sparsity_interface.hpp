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


#ifndef CASADI_SPARSITY_INTERFACE_HPP
#define CASADI_SPARSITY_INTERFACE_HPP

#include "casadi_misc.hpp"

namespace casadi {
  /** \brief Empty Base

      This class is extended in SWIG.

      \identifier{3c} */
  struct CASADI_EXPORT SparsityInterfaceCommon {};

#ifndef SWIG
  /** \brief Sparsity interface class

      This is a common base class for GenericMatrix (i.e. MX and Matrix<>) and Sparsity, introducing a
      uniform syntax and implementing common functionality using the curiously recurring template pattern
      (CRTP) idiom.\n

      \author Joel Andersson
      \date 2014

      \identifier{3d} */
  template<typename MatType>
  class SparsityInterface : public SparsityInterfaceCommon {
  protected:
    // Helper functions
    inline const MatType& self() const { return static_cast<const MatType&>(*this); }
    inline MatType& self() { return static_cast<MatType&>(*this); }
  public:

    /// \cond CLUTTER
    static std::vector< std::vector< MatType > >
      blocksplit(const MatType& x, const std::vector<casadi_int>& vert_offset,
                 const std::vector<casadi_int>& horz_offset);
    static std::vector< std::vector< MatType > >
      blocksplit(const MatType& x, casadi_int vert_incr, casadi_int horz_incr);
    static MatType veccat(const std::vector< MatType >& x);
    static MatType vec(const MatType& x);
    static MatType repmat(const MatType& x, casadi_int n, casadi_int m=1);
    static std::vector<casadi_int> offset(const std::vector< MatType > &v, bool vert=true);
    static std::vector< MatType > diagsplit(const MatType& x,
                                            const std::vector<casadi_int>& output_offset);
    static std::vector< MatType > diagsplit(const MatType& x, casadi_int incr);
    static std::vector< MatType > diagsplit(const MatType& x, casadi_int incr1, casadi_int incr2);
    static MatType mtimes(const std::vector<MatType> &args,
                          const std::string& blas = "reference");
    static std::vector<MatType > horzsplit(const MatType& x, casadi_int incr);
    static std::vector<MatType > vertsplit(const MatType& x, casadi_int incr);
    static std::vector<MatType > horzsplit_n(const MatType& x, casadi_int n);
    static std::vector<MatType > vertsplit_n(const MatType& x, casadi_int n);
    static MatType repmat(const MatType &A, const std::pair<casadi_int, casadi_int>& rc) {
      return MatType::repmat(A, rc.first, rc.second);
    }
    /// \endcond

/** \defgroup expression_tools Expression tools
* Functions for manipulating DM, SX, MX or Sparsity
*
*/
/**
\addtogroup expression_tools
@{
*/
    /** \brief Concatenate a list of matrices horizontally

        Places the matrices side by side, in order: the result has the row count of the inputs
        and a column count equal to the sum of the column counts. All inputs must have the same
        number of rows; a 0-by-0 matrix is ignored. The result is of the type of the arguments,
        numeric arguments are converted to the symbolic type (DM mixes with SX or MX, but SX and MX
        cannot be mixed). The sparsity pattern of each block is preserved. With no arguments the
        result is a 1-by-0 matrix.

        Alternative terminology: horizontal stack, hstack, horizontal append, [a b].
        Inverse of horzsplit: horzcat(horzsplit(x, ...)) == x.

        \param v List of matrices with matching row count
        \return Matrix with the blocks of \p v side by side

        \lang{python}
        \doctest
        >>> horzcat(DM([1, 2]), DM([3, 4]), DM([5, 6]))
        DM(
        [[1, 3, 5],
         [2, 4, 6]])
        >>> horzcat(SX.sym("x", 2), DM.ones(2, 1))
        SX(@1=1,
        [[x_0, @1],
         [x_1, @1]])
        >>> horzcat(MX.sym("x", 2, 2), MX.sym("y", 2, 0))
        MX(x)
        >>> horzcat().shape, horzcat(DM(2, 3), DM(0, 0)).shape
        ((1, 0), (2, 3))
        >>> horzcat(DM.ones(2, 2), DM.ones(3, 1))
        Traceback (most recent call last):
          ...
        RuntimeError: ...Mismatching number of rows...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [DM([1; 2]), DM([3; 4]), DM([5; 6])]
        [[1, 3, 5],
         [2, 4, 6]]
        >> horzcat(SX.sym('x', 2), DM.ones(2, 1))
        @1=1,
        [[x_0, @1],
         [x_1, @1]]
        >> [MX.sym('x', 2, 2), MX.sym('y', 2, 0)]
        x
        >> size([DM(2, 3), DM(0, 0)])
        2 3
        >> try, [DM.ones(2, 2), DM.ones(3, 1)]; catch e, disp(e.message); end
        ...Mismatching number of rows
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> [ca.DM([1, 2]) ca.DM([3, 4]) ca.DM([5, 6])]
        DM(
        [[1, 3, 5],
         [2, 4, 6]])
        julia> [ca.SX.sym("x", 2) ca.DM.ones(2, 1)]
        SX(@1=1,
        [[x_0, @1],
         [x_1, @1]])
        julia> [ca.MX.sym("x", 2, 2) ca.MX.sym("y", 2, 0)]
        MX(x)
        julia> size([ca.DM(2, 3) ca.DM(0, 0)])
        (2, 3)
        julia> try
                   [ca.DM.ones(2, 2) ca.DM.ones(3, 1)]
               catch err
                   println(occursin("dimension mismatch", sprint(showerror, err)))
               end
        false
        \enddoctest
        \endlang

        \sa vertcat, diagcat, veccat, horzsplit

        \identifier{3e} */
    inline friend MatType horzcat(const std::vector<MatType> &v) {
      return MatType::horzcat(v);
    }

    /** \brief Concatenate a list of matrices vertically

        Stacks the matrices on top of each other, in order: the result has the column count of the
        inputs and a row count equal to the sum of the row counts. All inputs must have the same
        number of columns; a 0-by-0 matrix is ignored. Type and sparsity rules are those of
        horzcat. With no arguments the result is a 0-by-1 matrix.

        Column vectors stack cheaply; for MX, stacking matrices with more than one column is
        rewritten as a horizontal concatenation of the transposes.

        Alternative terminology: vertical stack, vstack, vertical append, [a;b].
        Inverse of vertsplit: vertcat(vertsplit(x, ...)) == x.

        \param v List of matrices with matching column count
        \return Matrix with the blocks of \p v stacked vertically

        \lang{python}
        \doctest
        >>> vertcat(SX.sym("x", 2), 3)
        SX([x_0, x_1, 3])
        >>> vertcat(DM([[1, 2]]), DM.eye(2))
        DM(
        [[1, 2],
         [1, 00],
         [00, 1]])
        >>> vertcat(MX.sym("a", 2, 3), MX.sym("b", 1, 3))
        MX(horzcat(a', b')')
        >>> vertcat().shape
        (0, 1)
        >>> vertcat(MX.sym("a", 2, 2), MX.sym("b", 2, 3))
        Traceback (most recent call last):
          ...
        RuntimeError: ...vertcat dimension mismatch...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [SX.sym('x', 2); 3]
        [x_0, x_1, 3]
        >> vertcat(DM([1 2]), DM.eye(2))
        [[1, 2],
         [1, 00],
         [00, 1]]
        >> [MX.sym('a', 2, 3); MX.sym('b', 1, 3)]
        horzcat(a', b')'
        >> try, [MX.sym('a', 2, 2); MX.sym('b', 2, 3)]; catch e, disp(e.message); end
        ...vertcat dimension mismatch...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> [ca.SX.sym("x", 2); 3]
        SX([x_0, x_1, 3])
        julia> [ca.DM([1 2]); ca.DM.eye(2)]
        DM(
        [[1, 2],
         [1, 00],
         [00, 1]])
        julia> print([ca.MX.sym("a", 2, 3); ca.MX.sym("b", 1, 3)])
        horzcat(a', b')'
        julia> try
                   [ca.MX.sym("a", 2, 2); ca.MX.sym("b", 2, 3)]
               catch err
                   println(occursin("vertcat dimension mismatch", sprint(showerror, err)))
               end
        false
        \enddoctest
        \endlang

        \sa horzcat, diagcat, veccat, vertsplit

        \identifier{3f} */
    inline friend MatType vertcat(const std::vector<MatType> &v) {
      return MatType::vertcat(v);
    }

    /** \brief Split a matrix horizontally into groups of columns

        Returns the column blocks \p x[:, offset[i]:offset[i+1]], in order. The offsets must
        start at 0, be non-decreasing and end at the number of columns of \p x, otherwise an
        error is raised. Two equal consecutive offsets give an empty block. Each block carries the
        corresponding part of the sparsity pattern.
        Inverse of horzcat: horzcat(horzsplit(x, offset)) == x.

        \param x Matrix (or Sparsity) to split; any shape
        \param offset Start column of every block followed by the column count of \p x
        \return List of matrices of the same type as \p x, one per block

        \lang{python}
        \doctest
        >>> x = DM([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]])
        >>> horzsplit(x, [0, 2, 5])
        [DM(
        [[1, 2],
         [6, 7]]), DM(
        [[3, 4, 5],
         [8, 9, 10]])]
        >>> horzsplit(x, [0, 0, 5])[0].shape
        (2, 0)
        >>> horzsplit(x, [0, 2])
        Traceback (most recent call last):
          ...
        RuntimeError: ...must start at 0, be nondecreasing and end at 5...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = DM([1 2 3 4 5; 6 7 8 9 10]);
        >> v = horzsplit(x, [0, 2, 5]); disp(numel(v)), disp(v{1}), disp(v{2})
        2
        [[1, 2],
         [6, 7]]
        [[3, 4, 5],
         [8, 9, 10]]
        >> v = horzsplit(x, [0, 0, 5]); disp(size(v{1}))
        2 0
        >> try, horzsplit(x, [0, 2]); catch e, disp(e.message); end
        ...must start at 0, be nondecreasing and end at 5...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.DM([1 2 3 4 5; 6 7 8 9 10]);
        julia> v = ca.horzsplit(x, [0, 2, 5]);
        julia> length(v)
        2
        julia> v[1]
        DM(
        [[1, 2],
         [6, 7]])
        julia> v[2]
        DM(
        [[3, 4, 5],
         [8, 9, 10]])
        julia> size(ca.horzsplit(x, [0, 0, 5])[1])
        (2, 0)
        julia> try
                   ca.horzsplit(x, [0, 2])
               catch err
                   println(occursin("end at 5", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa horzsplit_n, vertsplit, offset

        \identifier{3g} */
    inline friend std::vector<MatType >
      horzsplit(const MatType &x, const std::vector<casadi_int>& offset) {
      return MatType::horzsplit(x, offset);
    }

    /** \brief Split a matrix horizontally into fixed-size groups of columns

        Returns the consecutive column blocks of \p x, each with \p incr columns, in order.
        The last block is smaller when the column count is not a multiple of \p incr.
        A matrix without columns gives an empty list. Each block carries the corresponding part
        of the sparsity pattern.
        Inverse of horzcat: horzcat(horzsplit(x, incr)) == x.

        \param x Matrix (or Sparsity) to split; any shape
        \param incr Number of columns per block, must be positive
        \return List of matrices of the same type as \p x

        \lang{python}
        \doctest
        >>> horzsplit(SX.sym("x", 1, 3))
        [SX(x_0), SX(x_1), SX(x_2)]
        >>> horzsplit(DM([[1, 2, 3, 4, 5]]), 2)
        [DM([[1, 2]]), DM([[3, 4]]), DM(5)]
        >>> horzsplit(MX.sym("x", 2, 4), 2)
        [MX(horzsplit(x){0}), MX(horzsplit(x){1})]
        >>> horzsplit(Sparsity.diag(3), 2)
        [Sparsity(3x2,2nz), Sparsity(3x1,1nz)]
        >>> horzsplit(MX.sym("x", 2, 0))
        []
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = horzsplit(SX.sym('x', 1, 3)); disp(numel(v)), disp(v{1}), disp(v{3})
        3
        x_0
        x_2
        >> v = horzsplit(DM([1 2 3 4 5]), 2); disp(numel(v)), disp(v{1}), disp(v{3})
        3
        [[1, 2]]
        5
        >> v = horzsplit(MX.sym('x', 2, 4), 2); disp(v{1}), disp(v{2})
        horzsplit(x){0}
        horzsplit(x){1}
        >> v = horzsplit(Sparsity.diag(3), 2); disp(v{1}), disp(v{2})
        3x2,2nz
        3x1,1nz
        >> v = horzsplit(MX.sym('x', 2, 0)); disp(numel(v))
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.horzsplit(ca.SX.sym("x", 1, 3))
        3-element Vector{SX}:
         SX(x_0)
         SX(x_1)
         SX(x_2)
        julia> ca.horzsplit(ca.DM([1 2 3 4 5]), 2)
        3-element Vector{DM}:
         DM([[1, 2]])
         DM([[3, 4]])
         DM(5)
        julia> ca.horzsplit(ca.MX.sym("x", 2, 4), 2)
        2-element Vector{MX}:
         MX(horzsplit(x){0})
         MX(horzsplit(x){1})
        julia> ca.horzsplit(ca.Sparsity.diag(3), 2)
        2-element Vector{Sparsity}:
         Sparsity(3x2,2nz)
         Sparsity(3x1,1nz)
        julia> ca.horzsplit(ca.MX.sym("x", 2, 0))
        MX[]
        \enddoctest
        \endlang

        \sa horzsplit_n, vertsplit

        \identifier{3h} */
    inline friend std::vector<MatType > horzsplit(const MatType& x, casadi_int incr=1) {
      return MatType::horzsplit(x, incr);
    }

    /** \brief Split a matrix horizontally into a given number of equal groups of columns

        Returns \p n column blocks of equal width, in order. An error is raised when the column
        count of \p x is not a multiple of \p n or \p n is negative. A matrix without columns gives
        \p n copies of itself.
        Inverse of horzcat: horzcat(horzsplit_n(x, n)) == x.

        \param x Matrix (or Sparsity) to split; any shape
        \param n Number of blocks
        \return List of \p n matrices of the same type as \p x

        \lang{python}
        \doctest
        >>> horzsplit_n(DM([[1, 2, 3, 4]]), 2)
        [DM([[1, 2]]), DM([[3, 4]])]
        >>> horzsplit_n(DM(2, 0), 3)
        [DM(2x0), DM(2x0), DM(2x0)]
        >>> horzsplit_n(DM([[1, 2, 3]]), 2)
        Traceback (most recent call last):
          ...
        RuntimeError: ...must be a multiple of n...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = horzsplit_n(DM([1 2 3 4]), 2); disp(v{1}), disp(v{2})
        [[1, 2]]
        [[3, 4]]
        >> v = horzsplit_n(DM(2, 0), 3); disp(numel(v)), disp(size(v{1}))
        3
        2 0
        >> try, horzsplit_n(DM([1 2 3]), 2); catch e, disp(e.message); end
        ...must be a multiple of n...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.horzsplit_n(ca.DM([1 2 3 4]), 2)
        2-element Vector{DM}:
         DM([[1, 2]])
         DM([[3, 4]])
        julia> ca.horzsplit_n(ca.DM(2, 0), 3)
        3-element Vector{DM}:
         DM(2x0)
         DM(2x0)
         DM(2x0)
        julia> try
                   ca.horzsplit_n(ca.DM([1 2 3]), 2)
               catch err
                   println(occursin("must be a multiple of n", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa horzsplit, vertsplit_n

        \identifier{277} */
    inline friend std::vector<MatType > horzsplit_n(const MatType& x, casadi_int n) {
      return MatType::horzsplit_n(x, n);
    }

    /** \brief Split a matrix vertically into groups of rows

        Returns the row blocks \p x[offset[i]:offset[i+1], :], in order. The offsets must start
        at 0, be non-decreasing and end at the number of rows of \p x, otherwise an error is
        raised. Each block carries the corresponding part of the sparsity pattern. The offsets of
        a list of matrices that were stacked with vertcat can be recovered with offset.
        Inverse of vertcat: vertcat(vertsplit(x, offset)) == x.

        \param x Matrix (or Sparsity) to split; any shape
        \param offset Start row of every block followed by the row count of \p x
        \return List of matrices of the same type as \p x, one per block

        \lang{python}
        \doctest
        >>> x = SX.sym("x", 5)
        >>> vertsplit(x, [0, 2, 5])
        [SX([x_0, x_1]), SX([x_2, x_3, x_4])]
        >>> parts = [DM.ones(2, 1), DM.ones(3, 1)]
        >>> vertsplit(x, offset(parts))
        [SX([x_0, x_1]), SX([x_2, x_3, x_4])]
        >>> vertsplit(DM([[1, 2], [3, 4]]), [0, 1, 2])
        [DM([[1, 2]]), DM([[3, 4]])]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x', 5);
        >> v = vertsplit(x, [0, 2, 5]); disp(numel(v)), disp(v{1}), disp(v{2})
        2
        [x_0, x_1]
        [x_2, x_3, x_4]
        >> v = vertsplit(DM([1 2; 3 4]), [0, 1, 2]); disp(v{1}), disp(v{2})
        [[1, 2]]
        [[3, 4]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x", 5);
        julia> ca.vertsplit(x, [0, 2, 5])
        2-element Vector{SX}:
         SX([x_0, x_1])
         SX([x_2, x_3, x_4])
        julia> parts = [ca.DM.ones(2, 1), ca.DM.ones(3, 1)];
        julia> ca.vertsplit(x, ca.offset(parts))
        2-element Vector{SX}:
         SX([x_0, x_1])
         SX([x_2, x_3, x_4])
        julia> ca.vertsplit(ca.DM([1 2; 3 4]), [0, 1, 2])
        2-element Vector{DM}:
         DM([[1, 2]])
         DM([[3, 4]])
        \enddoctest
        \endlang

        \sa vertsplit_n, horzsplit, offset

        \identifier{3i} */
    inline friend std::vector<MatType >
      vertsplit(const MatType& x, const std::vector<casadi_int>& offset) {
      return MatType::vertsplit(x, offset);
    }

    /** \brief Get the block offsets of a list of matrices

        Returns the cumulative row counts (or column counts when \p vert is false) of the
        matrices in \p v, starting with 0. The result has one entry more than \p v and is the
        offset list that vertsplit (or horzsplit) needs to undo a vertcat (or horzcat) of \p v.

        \param v List of matrices (or Sparsity patterns)
        \param vert True for row offsets, false for column offsets
        \return List of offsets of length len(v)+1

        \lang{python}
        \doctest
        >>> offset([DM(2, 3), DM(4, 3)])
        [0, 2, 6]
        >>> offset([DM(2, 3), DM(2, 4)], False)
        [0, 3, 7]
        >>> offset([])
        [0]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.offset([ca.DM(2, 3), ca.DM(4, 3)])
        3-element Vector{Int64}:
         0
         2
         6
        julia> ca.offset([ca.DM(2, 3), ca.DM(2, 4)], false)
        3-element Vector{Int64}:
         0
         3
         7
        julia> ca.offset(ca.DM[])
        1-element Vector{Int64}:
         0
        \enddoctest
        \endlang

        \sa vertsplit, horzsplit

        \identifier{3j} */
    inline friend std::vector<casadi_int > offset(const std::vector<MatType> &v, bool vert=true) {
      return MatType::offset(v, vert);
    }

    /** \brief Split a matrix vertically into fixed-size groups of rows

        Returns the consecutive row blocks of \p x, each with \p incr rows, in order.
        The last block is smaller when the row count is not a multiple of \p incr.
        A matrix without rows gives an empty list. Each block carries the corresponding part of
        the sparsity pattern.
        Inverse of vertcat: vertcat(vertsplit(x, incr)) == x.

        \param x Matrix (or Sparsity) to split; any shape
        \param incr Number of rows per block, must be positive
        \return List of matrices of the same type as \p x

        \lang{python}
        \doctest
        >>> vertsplit(SX.sym("a", 4))
        [SX(a_0), SX(a_1), SX(a_2), SX(a_3)]
        >>> vertsplit(SX.sym("a", 4), 2)
        [SX([a_0, a_1]), SX([a_2, a_3])]
        >>> vertsplit(DM([0, 1, 2, 3, 4]), 2)
        [DM([0, 1]), DM([2, 3]), DM(4)]
        >>> vertsplit(MX.sym("x", 4), 2)
        [MX(vertsplit(x){0}), MX(vertsplit(x){1})]
        >>> vertsplit(MX.sym("x", 0, 2))
        []
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = vertsplit(SX.sym('a', 4)); disp(numel(v)), disp(v{1}), disp(v{4})
        4
        a_0
        a_3
        >> v = vertsplit(SX.sym('a', 4), 2); disp(v{1}), disp(v{2})
        [a_0, a_1]
        [a_2, a_3]
        >> v = vertsplit(DM([0; 1; 2; 3; 4]), 2); disp(numel(v)), disp(v{3})
        3
        4
        >> v = vertsplit(MX.sym('x', 4), 2); disp(v{1}), disp(v{2})
        vertsplit(x){0}
        vertsplit(x){1}
        >> v = vertsplit(MX.sym('x', 0, 2)); disp(numel(v))
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.vertsplit(ca.SX.sym("a", 4))
        4-element Vector{SX}:
         SX(a_0)
         SX(a_1)
         SX(a_2)
         SX(a_3)
        julia> ca.vertsplit(ca.SX.sym("a", 4), 2)
        2-element Vector{SX}:
         SX([a_0, a_1])
         SX([a_2, a_3])
        julia> ca.vertsplit(ca.DM([0, 1, 2, 3, 4]), 2)
        3-element Vector{DM}:
         DM([0, 1])
         DM([2, 3])
         DM(4)
        julia> ca.vertsplit(ca.MX.sym("x", 4), 2)
        2-element Vector{MX}:
         MX(vertsplit(x){0})
         MX(vertsplit(x){1})
        julia> ca.vertsplit(ca.MX.sym("x", 0, 2))
        MX[]
        \enddoctest
        \endlang

        \sa vertsplit_n, horzsplit

        \identifier{3k} */
    inline friend std::vector<MatType > vertsplit(const MatType &x, casadi_int incr=1) {
      return MatType::vertsplit(x, incr);
    }

    /** \brief Split a matrix vertically into a given number of equal groups of rows

        Returns \p n row blocks of equal height, in order. An error is raised when the row count
        of \p x is not a multiple of \p n or \p n is negative. A matrix without rows gives \p n
        copies of itself.
        Inverse of vertcat: vertcat(vertsplit_n(x, n)) == x.

        \param x Matrix (or Sparsity) to split; any shape
        \param n Number of blocks
        \return List of \p n matrices of the same type as \p x

        \lang{python}
        \doctest
        >>> vertsplit_n(DM([1, 2, 3, 4]), 2)
        [DM([1, 2]), DM([3, 4])]
        >>> vertsplit_n(DM([1, 2, 3]), 2)
        Traceback (most recent call last):
          ...
        RuntimeError: ...must be a multiple of n...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = vertsplit_n(DM([1; 2; 3; 4]), 2); disp(v{1}), disp(v{2})
        [1, 2]
        [3, 4]
        >> try, vertsplit_n(DM([1; 2; 3]), 2); catch e, disp(e.message); end
        ...must be a multiple of n...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.vertsplit_n(ca.DM([1, 2, 3, 4]), 2)
        2-element Vector{DM}:
         DM([1, 2])
         DM([3, 4])
        julia> try
                   ca.vertsplit_n(ca.DM([1, 2, 3]), 2)
               catch err
                   println(occursin("must be a multiple of n", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa vertsplit, horzsplit_n

        \identifier{278} */
    inline friend std::vector<MatType > vertsplit_n(const MatType& x, casadi_int n) {
      return MatType::vertsplit_n(x, n);
    }

    /** \brief Construct a matrix from a list of lists of blocks

        Each inner list is a block row: its entries are concatenated horizontally, and the
        resulting block rows are concatenated vertically. The usual concatenation rules apply:
        blocks in one block row need the same number of rows, block rows need the same total
        number of columns, otherwise an error is raised.
        Inverse of blocksplit: blockcat(blocksplit(x, ...)) == x.

        \param v List of block rows, each a list of matrices
        \return Matrix of the same type as the blocks

        \lang{python}
        \doctest
        >>> blockcat([[DM(1), DM([[2, 3]])], [DM([4, 5]), DM.eye(2)]])
        DM(
        [[1, 2, 3],
         [4, 1, 00],
         [5, 00, 1]])
        >>> a, b, c, d = MX.sym("a", 2, 2), MX.sym("b", 2, 1), MX.sym("c", 1, 2), MX.sym("d")
        >>> blockcat([[a, b], [c, d]])
        MX(horzcat(horzcat(a, b)', horzcat(c, d)')')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> blockcat({{DM(1), DM([2 3])}, {DM([4; 5]), DM.eye(2)}})
        [[1, 2, 3],
         [4, 1, 00],
         [5, 00, 1]]
        >> a = MX.sym('a', 2, 2); b = MX.sym('b', 2, 1); c = MX.sym('c', 1, 2); d = MX.sym('d');
        >> blockcat({{a, b}, {c, d}})
        horzcat(horzcat(a, b)', horzcat(c, d)')'
        >> [a, b; c, d]
        horzcat(horzcat(a, b)', horzcat(c, d)')'
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.blockcat([[ca.DM(1), ca.DM([2 3])], [ca.DM([4, 5]), ca.DM.eye(2)]])
        DM(
        [[1, 2, 3],
         [4, 1, 00],
         [5, 00, 1]])
        julia> a, b = ca.MX.sym("a", 2, 2), ca.MX.sym("b", 2, 1);
        julia> c, d = ca.MX.sym("c", 1, 2), ca.MX.sym("d");
        julia> ca.blockcat([[a, b], [c, d]])
        MX(horzcat(horzcat(a, b)', horzcat(c, d)')')
        julia> print([a b; c d])
        horzcat(horzcat(a, b)', horzcat(c, d)')'
        \enddoctest
        \endlang

        \sa horzcat, vertcat, blocksplit

        \identifier{3l} */
    inline friend MatType blockcat(const std::vector< std::vector<MatType > > &v) {
      return MatType::blockcat(v);
    }

    /** \brief Construct a matrix from 4 blocks

        Shorthand for the two-by-two block layout [[A, B], [C, D]].

        \lang{python}
        \doctest
        >>> blockcat(1, 2, 3, 4)
        DM(
        [[1, 2],
         [3, 4]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [DM(1), 2; 3, 4]
        [[1, 2],
         [3, 4]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.blockcat(ca.DM(1), ca.DM(2), ca.DM(3), ca.DM(4))
        DM(
        [[1, 2],
         [3, 4]])
        \enddoctest
        \endlang

        \identifier{3m} */
    inline friend MatType
      blockcat(const MatType &A, const MatType &B, const MatType &C, const MatType &D) {
      return vertcat(horzcat(A, B), horzcat(C, D));
    }

    /** \brief Split a matrix into a grid of blocks

        The rows are split at \p vert_offset and every row block is then split at \p horz_offset.
        Both offset lists follow the rules of vertsplit and horzsplit: start at 0, non-decreasing,
        end at the row (column) count.
        Inverse of blockcat: blockcat(blocksplit(x, vert_offset, horz_offset)) == x.

        \param x Matrix (or Sparsity) to split; any shape
        \param vert_offset Row boundaries of the block rows
        \param horz_offset Column boundaries of the block columns
        \return List of block rows, each a list of matrices of the same type as \p x

        \lang{python}
        \doctest
        >>> x = DM([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> blocksplit(x, [0, 1, 3], [0, 2, 3])
        [[DM([[1, 2]]), DM(3)], [DM(
        [[4, 5],
         [7, 8]]), DM([6, 9])]]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = DM([1 2 3; 4 5 6; 7 8 9]);
        >> b = blocksplit(x, [0, 1, 3], [0, 2, 3]); disp(size(b)), disp(b{1}{1}), disp(b{2}{1})
        1 2
        [[1, 2]]
        [[4, 5],
         [7, 8]]
        >> b{1}{2}, b{2}{2}
        3
        [6, 9]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.DM([1 2 3; 4 5 6; 7 8 9]);
        julia> ca.blocksplit(x, [0, 1, 3], [0, 2, 3])
        2-element Vector{Vector{DM}}:
         DM[[[1, 2]], 3]
         DM[
        [[4, 5],
         [7, 8]], [6, 9]]
        \enddoctest
        \endlang

        \sa blockcat, vertsplit, horzsplit

        \identifier{3n} */
    inline friend std::vector< std::vector< MatType > >
      blocksplit(const MatType& x,
                 const std::vector<casadi_int>& vert_offset,
                 const std::vector<casadi_int>& horz_offset) {
      return MatType::blocksplit(x, vert_offset, horz_offset);
    }

    /** \brief Split a matrix into a grid of fixed-size blocks

        Every block has \p vert_incr rows and \p horz_incr columns; the blocks in the last block
        row or column are smaller when the dimensions are not multiples of the increments.
        Inverse of blockcat: blockcat(blocksplit(x, vert_incr, horz_incr)) == x.

        \param x Matrix (or Sparsity) to split; any shape
        \param vert_incr Rows per block, must be positive
        \param horz_incr Columns per block, must be positive
        \return List of block rows, each a list of matrices of the same type as \p x

        \lang{python}
        \doctest
        >>> x = DM([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> blocksplit(x, 2, 2)
        [[DM(
        [[1, 2],
         [4, 5]]), DM([3, 6])], [DM([[7, 8]]), DM(9)]]
        >>> blockcat(blocksplit(x, 2, 2))
        DM(
        [[1, 2, 3],
         [4, 5, 6],
         [7, 8, 9]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = DM([1 2 3; 4 5 6; 7 8 9]);
        >> b = blocksplit(x, 2, 2); disp(b{1}{1}), disp(b{1}{2}), disp(b{2}{1}), disp(b{2}{2})
        [[1, 2],
         [4, 5]]
        [3, 6]
        [[7, 8]]
        9
        >> blockcat(blocksplit(x, 2, 2))
        [[1, 2, 3],
         [4, 5, 6],
         [7, 8, 9]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.DM([1 2 3; 4 5 6; 7 8 9]);
        julia> ca.blocksplit(x, 2, 2)
        2-element Vector{Vector{DM}}:
         DM[
        [[1, 2],
         [4, 5]], [3, 6]]
         DM[[[7, 8]], 9]
        julia> ca.blockcat(ca.blocksplit(x, 2, 2))
        DM(
        [[1, 2, 3],
         [4, 5, 6],
         [7, 8, 9]])
        \enddoctest
        \endlang

        \sa blockcat

        \identifier{3o} */
    inline friend std::vector< std::vector< MatType > >
      blocksplit(const MatType& x, casadi_int vert_incr=1, casadi_int horz_incr=1) {
      return MatType::blocksplit(x, vert_incr, horz_incr);
    }

    /** \brief Construct a block-diagonal matrix from a list of matrices

        Places the matrices along the diagonal, in order; all other entries are structural
        zeros. The blocks need not be square: the result has the sum of the row counts by the sum
        of the column counts. Type rules are those of horzcat. With no arguments the result is a
        0-by-0 matrix.
        Inverse of diagsplit: diagcat(diagsplit(x, ...)) == x.

        \param A List of matrices of any shape
        \return Block-diagonal matrix of the same type as the blocks

        \lang{python}
        \doctest
        >>> diagcat(DM.ones(1, 2), 3)
        DM(
        [[1, 1, 00],
         [00, 00, 3]])
        >>> diagcat(SX.sym("x"), SX.sym("y", 2))
        SX(
        [[x, 00],
         [00, y_0],
         [00, y_1]])
        >>> diagcat(MX.sym("x"), MX.sym("y", 2, 2))
        MX(diagcat(x, y))
        >>> diagcat().shape
        (0, 0)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> diagcat(DM.ones(1, 2), 3)
        [[1, 1, 00],
         [00, 00, 3]]
        >> diagcat(SX.sym('x'), SX.sym('y', 2))
        [[x, 00],
         [00, y_0],
         [00, y_1]]
        >> diagcat(MX.sym('x'), MX.sym('y', 2, 2))
        diagcat(x, y)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.diagcat(ca.DM.ones(1, 2), 3)
        DM(
        [[1, 1, 00],
         [00, 00, 3]])
        julia> ca.diagcat(ca.SX.sym("x"), ca.SX.sym("y", 2))
        SX(
        [[x, 00],
         [00, y_0],
         [00, y_1]])
        julia> ca.diagcat(ca.MX.sym("x"), ca.MX.sym("y", 2, 2))
        MX(diagcat(x, y))
        \enddoctest
        \endlang

        \sa horzcat, vertcat, diagsplit

        \identifier{3p} */
    inline friend MatType diagcat(const std::vector<MatType> &A) {
      return MatType::diagcat(A);
    }

    /** \brief Split a matrix into diagonal blocks of given row and column boundaries

        Block i is \p x[output_offset1[i]:output_offset1[i+1],
        output_offset2[i]:output_offset2[i+1]]. Both lists must have the same length, start at 0,
        be non-decreasing and end at the row (column) count of \p x. Entries of \p x outside the
        diagonal blocks must be zero (structurally zero for MX and Sparsity), otherwise an
        error is raised.
        Inverse of diagcat: diagcat(diagsplit(x, ...)) == x for block-diagonal \p x.

        \param x Matrix (or Sparsity) to split; any shape
        \param output_offset1 Row boundaries of the blocks
        \param output_offset2 Column boundaries of the blocks
        \return List of matrices of the same type as \p x, one per block

        \lang{python}
        \doctest
        >>> x = DM([[1, 2, 0], [0, 0, 6]])
        >>> diagsplit(x, [0, 1, 2], [0, 2, 3])
        [DM([[1, 2]]), DM(6)]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = DM([1 2 0; 0 0 6]);
        >> v = diagsplit(x, [0, 1, 2], [0, 2, 3]); disp(v{1}), disp(v{2})
        [[1, 2]]
        6
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.DM([1 2 0; 0 0 6]);
        julia> ca.diagsplit(x, [0, 1, 2], [0, 2, 3])
        2-element Vector{DM}:
         DM([[1, 2]])
         DM(6)
        \enddoctest
        \endlang

        \sa diagcat

        \identifier{3q} */
    friend std::vector< MatType >
      diagsplit(const MatType& x,
                const std::vector<casadi_int>& output_offset1,
                const std::vector<casadi_int>& output_offset2) {
      return MatType::diagsplit(x, output_offset1, output_offset2);
    }

    /** \brief Split a square matrix into square diagonal blocks of given boundaries

        Same as the two-offset variant with \p output_offset used for both rows and columns.
        An error is raised when \p x is not square.

        \param x Square matrix (or Sparsity) to split
        \param output_offset Boundaries of the blocks

        \lang{python}
        \doctest
        >>> x = DM([[1, 2, 0], [3, 4, 0], [0, 0, 5]])
        >>> diagsplit(x, [0, 2, 3])
        [DM(
        [[1, 2],
         [3, 4]]), DM(5)]
        >>> diagsplit(DM.ones(2, 3), [0, 2])
        Traceback (most recent call last):
          ...
        RuntimeError: ...must be square...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = DM([1 2 0; 3 4 0; 0 0 5]);
        >> v = diagsplit(x, [0, 2, 3]); disp(v{1}), disp(v{2})
        [[1, 2],
         [3, 4]]
        5
        >> try, diagsplit(DM.ones(2, 3), [0, 2]); catch e, disp(e.message); end
        ...must be square...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.DM([1 2 0; 3 4 0; 0 0 5]);
        julia> v = ca.diagsplit(x, [0, 2, 3]);
        julia> v[1]
        DM(
        [[1, 2],
         [3, 4]])
        julia> v[2]
        DM(5)
        julia> try
                   ca.diagsplit(ca.DM.ones(2, 3), [0, 2])
               catch err
                   println(occursin("must be square", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa diagcat

        \identifier{3r} */
    inline friend std::vector< MatType >
      diagsplit(const MatType& x, const std::vector<casadi_int>& output_offset) {
      return MatType::diagsplit(x, output_offset);
    }

    /** \brief Split a square matrix into square diagonal blocks of fixed size

        Returns the consecutive \p incr-by-\p incr diagonal blocks of \p x; the last block is
        smaller when the dimension is not a multiple of \p incr. An error is raised when \p x is
        not square. Entries outside the blocks must be zero (structurally zero for MX and
        Sparsity).
        Inverse of diagcat: diagcat(diagsplit(x, incr)) == x for block-diagonal \p x.

        \param x Square matrix (or Sparsity) to split
        \param incr Size of each block, must be positive

        \lang{python}
        \doctest
        >>> diagsplit(DM([[1, 0], [0, 4]]))
        [DM(1), DM(4)]
        >>> diagsplit(diagcat(SX.sym("a", 2, 2), SX.sym("b", 2, 2)), 2)
        [SX(
        [[a_0, a_2],
         [a_1, a_3]]), SX(
        [[b_0, b_2],
         [b_1, b_3]])]
        >>> diagsplit(MX.sym("x", Sparsity.diag(4)), 2)
        [MX(diagsplit(x){0}), MX(diagsplit(x){1})]
        >>> diagsplit(MX.sym("x", 4, 4), 2)
        Traceback (most recent call last):
          ...
        RuntimeError: ...nonzeros outside the diagonal blocks...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = diagsplit(DM([1 0; 0 4])); disp(v{1}), disp(v{2})
        1
        4
        >> v = diagsplit(diagcat(SX.sym('a', 2, 2), SX.sym('b', 2, 2)), 2); disp(v{1}), disp(v{2})
        [[a_0, a_2],
         [a_1, a_3]]
        [[b_0, b_2],
         [b_1, b_3]]
        >> v = diagsplit(MX.sym('x', Sparsity.diag(4)), 2); disp(v{1}), disp(v{2})
        diagsplit(x){0}
        diagsplit(x){1}
        >> try, diagsplit(MX.sym('x', 4, 4), 2); catch e, disp(e.message); end
        ...nonzeros outside the diagonal blocks...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.diagsplit(ca.DM([1 0; 0 4]))
        2-element Vector{DM}:
         DM(1)
         DM(4)
        julia> v = ca.diagsplit(ca.diagcat(ca.SX.sym("a", 2, 2), ca.SX.sym("b", 2, 2)), 2);
        julia> v[1]
        SX(
        [[a_0, a_2],
         [a_1, a_3]])
        julia> v[2]
        SX(
        [[b_0, b_2],
         [b_1, b_3]])
        julia> ca.diagsplit(ca.MX.sym("x", ca.Sparsity.diag(4)), 2)
        2-element Vector{MX}:
         MX(diagsplit(x){0})
         MX(diagsplit(x){1})
        julia> try
                   ca.diagsplit(ca.MX.sym("x", 4, 4), 2)
               catch err
                   println(occursin("nonzeros outside the diagonal blocks", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa diagcat

        \identifier{3s} */
    inline friend std::vector< MatType >
      diagsplit(const MatType& x, casadi_int incr=1) {
      return MatType::diagsplit(x, incr);
    }

    /** \brief Split a matrix into diagonal blocks of fixed size

        Block i covers rows i*incr1 onwards and columns i*incr2 onwards. The number of row
        blocks must equal the number of column blocks.

        \param x Matrix (or Sparsity) to split
        \param incr1 Rows per block, must be positive
        \param incr2 Columns per block, must be positive

        \lang{python}
        \doctest
        >>> diagsplit(DM([[1, 2, 0, 0], [0, 0, 7, 8]]), 1, 2)
        [DM([[1, 2]]), DM([[7, 8]])]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> v = diagsplit(DM([1 2 0 0; 0 0 7 8]), 1, 2); disp(v{1}), disp(v{2})
        [[1, 2]]
        [[7, 8]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.diagsplit(ca.DM([1 2 0 0; 0 0 7 8]), 1, 2)
        2-element Vector{DM}:
         DM([[1, 2]])
         DM([[7, 8]])
        \enddoctest
        \endlang

        \sa diagcat

        \identifier{3t} */
    inline friend std::vector< MatType >
      diagsplit(const MatType& x, casadi_int incr1, casadi_int incr2) {
      return MatType::diagsplit(x, incr1, incr2);
    }

    /** \brief Concatenate vertically after vectorizing every argument

        Every argument is first turned into a column with vec (column-major), then the columns
        are stacked. Arguments can have any shape, so this is the common way to collect decision
        variables into one vector. With no arguments the result is a 0-by-1 matrix.

        \param x List of matrices of any shape
        \return Column vector of the same type as the arguments

        \lang{python}
        \doctest
        >>> veccat(DM([[1, 2], [3, 4]]), 5)
        DM([1, 3, 2, 4, 5])
        >>> veccat(MX.sym("x", 2, 2), MX.sym("y", 1, 3))
        MX(vertcat(vec(x), y'))
        >>> veccat().shape
        (0, 1)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> veccat(DM([1 2; 3 4]), 5)
        [1, 3, 2, 4, 5]
        >> veccat(MX.sym('x', 2, 2), MX.sym('y', 1, 3))
        vertcat(vec(x), y')
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.veccat([ca.MX.sym("x", 2, 2), ca.MX.sym("y", 1, 3)])
        MX(vertcat(vec(x), y'))
        julia> [vec(ca.DM([1 2; 3 4])); 5]
        DM([1, 3, 2, 4, 5])
        \enddoctest
        \endlang

        \sa vec, vertcat

        \identifier{3u} */
    inline friend MatType veccat(const std::vector< MatType >& x) {
      return MatType::veccat(x);
    }

    /** \brief Matrix product of two matrices

        Computes x*y in the linear-algebra sense; the column count of \p x must equal the row
        count of \p y, otherwise an error is raised. When either argument is a scalar the result
        is the elementwise product instead. The sparsity pattern of the result is the structural
        product of the patterns (entries that are structurally zero in every term stay structural
        zeros). For Sparsity arguments only the pattern of the product is computed.

        Optional \p blas selects the BLAS plugin used for dense numeric kernels. The choice only
        takes effect on numeric DM / MX dense * dense; symbolic types (Sparsity, SX) and the
        sparse fallback path ignore it.

        Values of \p blas:
        - `reference` (default): built-in reference kernel (triple loop)
        - `classic`: whatever BLAS was built into casadi via WITH_LAPACK (typically OpenBLAS)
        - `blasfeo`: BLASFEO

        Plugins (`classic`, `blasfeo`) are loaded lazily on first use; if the corresponding
        plugin DLL isn't present in the build, calling them raises an error.

        \param x Left factor
        \param y Right factor
        \param blas Name of the BLAS plugin for dense numeric products
        \return Matrix of size x.size1() by y.size2()

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> mtimes(A, A)
        DM(
        [[7, 10],
         [15, 22]])
        >>> mtimes(A, 2)
        DM(
        [[2, 4],
         [6, 8]])
        >>> mtimes(Sparsity.diag(2), Sparsity.dense(2, 1))
        Sparsity(2x1)
        >>> mtimes(MX.sym("x", 2, 3), MX.sym("y", 3, 1))
        MX(mac(x,y,zeros(2x1)))
        >>> mtimes(A, DM.ones(3, 1))
        Traceback (most recent call last):
          ...
        RuntimeError: ...incompatible dimensions...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> A * A
        [[7, 10],
         [15, 22]]
        >> mtimes(A, 2)
        [[2, 4],
         [6, 8]]
        >> MX.sym('x', 2, 3) * MX.sym('y', 3, 1)
        mac(x,y,zeros(2x1))
        >> try, A * DM.ones(3, 1); catch e, disp(e.message); end
        ...incompatible dimensions...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> A * A
        DM(
        [[7, 10],
         [15, 22]])
        julia> A * 2
        DM(
        [[2, 4],
         [6, 8]])
        julia> ca.mtimes(ca.Sparsity.diag(2), ca.Sparsity.dense(2, 1))
        Sparsity(2x1)
        julia> ca.MX.sym("x", 2, 3) * ca.MX.sym("y", 3, 1)
        MX(mac(x,y,zeros(2x1)))
        julia> try
                   A * ca.DM.ones(3, 1)
               catch err
                   println(occursin("incompatible dimensions", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa mac, times, norm_0_mul

        \identifier{3v} */
    inline friend MatType mtimes(const MatType &x, const MatType &y,
                                 const std::string& blas = "reference") {
      return MatType::mtimes(x, y, blas);
    }

    /** \brief Matrix product of n matrices

        Multiplies the matrices left to right: mtimes([a, b, c]) == mtimes(mtimes(a, b), c).
        An empty list raises an error; a single matrix is returned as is.
        See the 2-arg `mtimes` for the meaning of `blas`.

        \lang{python}
        \doctest
        >>> mtimes([DM([[1, 2]]), DM([[1, 0], [0, 2]]), DM([1, 1])])
        DM(5)
        >>> mtimes([])
        Traceback (most recent call last):
          ...
        RuntimeError: ...must not be empty...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM([1 2]) * ca.DM([1 0; 0 2]) * ca.DM([1, 1])
        DM(5)
        julia> try
                   ca.mtimes(ca.DM[])
               catch err
                   println(occursin("must not be empty", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \identifier{3w} */
    inline friend MatType mtimes(const std::vector<MatType> &args,
                                 const std::string& blas = "reference") {
      return MatType::mtimes(args, blas);
    }

    /** \brief Multiply-accumulate operation

        Matrix product of two matrices (x and y), adding the result to a third matrix z.
        The result has the same sparsity pattern as \p z, meaning that other entries of the
        product are ignored. The operation is equivalent to
        z + project(mtimes(x, y), z.sparsity()). The dimensions of \p z must match the product,
        otherwise an error is raised. When \p x or \p y is a scalar the result is z + x*y.
        See `mtimes` for the meaning of the `blas` argument.

        mtimes is implemented as mac with a structurally zero \p z of the product pattern.

        \param x Left factor
        \param y Right factor
        \param z Matrix to add to; its sparsity pattern selects the computed entries
        \return Matrix with the sparsity pattern of \p z

        \lang{python}
        \doctest
        >>> A = DM([[1, 2], [3, 4]])
        >>> mac(A, A, DM.eye(2))
        DM(
        [[8, 00],
         [00, 23]])
        >>> mac(A, A, DM.ones(2, 2))
        DM(
        [[8, 11],
         [16, 23]])
        >>> mac(MX.sym("x", 2, 3), MX.sym("y", 3, 2), MX.sym("z", Sparsity.diag(2)))
        MX(mac(x,y,z))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([1 2; 3 4]);
        >> mac(A, A, DM.eye(2))
        [[8, 00],
         [00, 23]]
        >> mac(A, A, DM.ones(2, 2))
        [[8, 11],
         [16, 23]]
        >> mac(MX.sym('x', 2, 3), MX.sym('y', 3, 2), MX.sym('z', Sparsity.diag(2)))
        mac(x,y,z)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([1 2; 3 4]);
        julia> ca.mac(A, A, ca.DM.eye(2))
        DM(
        [[8, 00],
         [00, 23]])
        julia> ca.mac(A, A, ca.DM.ones(2, 2))
        DM(
        [[8, 11],
         [16, 23]])
        julia> x, y = ca.MX.sym("x", 2, 3), ca.MX.sym("y", 3, 2);
        julia> ca.mac(x, y, ca.MX.sym("z", ca.Sparsity.diag(2)))
        MX(mac(x,y,z))
        \enddoctest
        \endlang

        \sa mtimes

        \identifier{3x} */
    inline friend MatType
      mac(const MatType &x, const MatType &y, const MatType &z,
          const std::string& blas = "reference") {
      return MatType::mac(x, y, z, blas);
    }

    /** \brief Transpose

        Returns the transpose of \p X; the sparsity pattern is transposed along. Same as X.T().

        \lang{python}
        \doctest
        >>> transpose(DM([[1, 2, 3]]))
        DM([1, 2, 3])
        >>> transpose(MX.sym("x", 2, 3))
        MX(x')
        >>> transpose(Sparsity.lower(2)) == Sparsity.upper(2)
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM([1 2 3])'
        [1, 2, 3]
        >> MX.sym('x', 2, 3)'
        x'
        >> transpose(Sparsity.lower(2)) == Sparsity.upper(2)
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM([1 2 3])'
        DM([1, 2, 3])
        julia> ca.MX.sym("x", 2, 3)'
        MX(x')
        julia> ca.is_equal(ca.transpose(ca.Sparsity.lower(2)), ca.Sparsity.upper(2))
        true
        \enddoctest
        \endlang

        \identifier{3y} */
    inline friend MatType transpose(const MatType& X) {
      return X.T();
    }

    /** \brief Reshape a matrix into a column vector

        Stacks the columns of \p x on top of each other (column-major order), giving a column
        of length x.numel(). Same as reshape(x, x.numel(), 1); a column is returned unchanged.
        The sparsity pattern is carried along: no entries are densified.

        a c \n
        b d \n

        turns into

        a \n
        b \n
        c \n
        d \n

        \lang{python}
        \doctest
        >>> vec(DM([[1, 2], [3, 4]]))
        DM([1, 3, 2, 4])
        >>> vec(DM.eye(2))
        DM([1, 00, 00, 1])
        >>> vec(MX.sym("x", 1, 3))
        MX(x')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> vec(DM([1 2; 3 4]))
        [1, 3, 2, 4]
        >> vec(DM.eye(2))
        [1, 00, 00, 1]
        >> vec(MX.sym('x', 1, 3))
        x'
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> vec(ca.DM([1 2; 3 4]))
        DM([1, 3, 2, 4])
        julia> vec(ca.DM.eye(2))
        DM([1, 00, 00, 1])
        julia> vec(ca.MX.sym("x", 1, 3))
        MX(x')
        \enddoctest
        \endlang

        \sa veccat, reshape

        \identifier{3z} */
    inline friend MatType vec(const MatType& x) {
      return MatType::vec(x);
    }

    /** \brief Returns a reshaped version of the matrix

        Reinterprets the entries of \p x in column-major order as a matrix of \p nrow by \p ncol,
        as reshape does in MATLAB and numpy with order='F'. The number of elements must not
        change, otherwise an error is raised. One of the dimensions can be -1, it is then inferred
        from the other. The sparsity pattern is reshaped along; no entries are densified. For MX
        this is a cheap reordering node.

        \param x Matrix (or Sparsity) to reshape; any shape
        \param nrow Number of rows of the result, or -1
        \param ncol Number of columns of the result, or -1

        \lang{python}
        \doctest
        >>> x = DM([[1, 2, 3], [4, 5, 6]])
        >>> reshape(x, 3, 2)
        DM(
        [[1, 5],
         [4, 3],
         [2, 6]])
        >>> reshape(x, -1, 1)
        DM([1, 4, 2, 5, 3, 6])
        >>> reshape(DM.eye(2), 1, 4)
        DM([[1, 00, 00, 1]])
        >>> reshape(MX.sym("x", 2, 3), 3, 2)
        MX(reshape(x))
        >>> reshape(x, 4, 2)
        Traceback (most recent call last):
          ...
        RuntimeError: ...number of elements must remain the same...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = DM([1 2 3; 4 5 6]);
        >> reshape(x, 3, 2)
        [[1, 5],
         [4, 3],
         [2, 6]]
        >> reshape(x, -1, 1)
        [1, 4, 2, 5, 3, 6]
        >> reshape(DM.eye(2), 1, 4)
        [[1, 00, 00, 1]]
        >> reshape(MX.sym('x', 2, 3), 3, 2)
        reshape(x)
        >> try, reshape(x, 4, 2); catch e, disp(e.message); end
        ...number of elements must remain the same...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.DM([1 2 3; 4 5 6]);
        julia> reshape(x, 3, 2)
        DM(
        [[1, 5],
         [4, 3],
         [2, 6]])
        julia> reshape(x, :, 1)
        DM([1, 4, 2, 5, 3, 6])
        julia> reshape(ca.DM.eye(2), 1, 4)
        DM([[1, 00, 00, 1]])
        julia> reshape(ca.MX.sym("x", 2, 3), 3, 2)
        MX(reshape(x))
        julia> try
                   reshape(x, 4, 2)
               catch err
                   println(occursin("number of elements must remain", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa vec, sparsity_cast, transpose

        \identifier{40} */
    inline friend MatType reshape(const MatType& x, casadi_int nrow, casadi_int ncol) {
      return MatType::reshape(x, nrow, ncol);
    }

    /** \brief Returns a reshaped version of the matrix, dimensions as a vector

        Same as reshape(x, rc[0], rc[1]).

        \lang{python}
        \doctest
        >>> reshape(DM([1, 2, 3, 4]), (2, 2))
        DM(
        [[1, 3],
         [2, 4]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> reshape(DM([1; 2; 3; 4]), [2, 2])
        [[1, 3],
         [2, 4]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> reshape(ca.DM([1, 2, 3, 4]), (2, 2))
        DM(
        [[1, 3],
         [2, 4]])
        \enddoctest
        \endlang

        \identifier{41} */
    inline friend MatType reshape(const MatType& x, std::pair<casadi_int, casadi_int> rc) {
      return MatType::reshape(x, rc.first, rc.second);
    }

    /** \brief Reshape the matrix to a given sparsity pattern

        \p sp must be the reshaped pattern of x.sparsity() (same number of elements and the same
        nonzeros in column-major order), otherwise an error is raised. Useful when the target
        pattern is already available, e.g. to undo a vec.

        \lang{python}
        \doctest
        >>> reshape(SX.sym("x", 2, 2), Sparsity.dense(4, 1))
        SX([x_0, x_1, x_2, x_3])
        >>> reshape(MX.sym("x", Sparsity.diag(2)), Sparsity.dense(4, 1))
        Traceback (most recent call last):
          ...
        RuntimeError: ...Reshape mismatch...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> reshape(SX.sym('x', 2, 2), Sparsity.dense(4, 1))
        [x_0, x_1, x_2, x_3]
        >> x = MX.sym('x', Sparsity.diag(2));
        >> try, reshape(x, Sparsity.dense(4, 1)); catch e, disp(e.message); end
        ...Reshape mismatch...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.reshape(ca.SX.sym("x", 2, 2), ca.Sparsity.dense(4, 1))
        SX([x_0, x_1, x_2, x_3])
        julia> try
                   ca.reshape(ca.MX.sym("x", ca.Sparsity.diag(2)), ca.Sparsity.dense(4, 1))
               catch err
                   println(occursin("Reshape mismatch", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa sparsity_cast

        \identifier{42} */
    inline friend MatType reshape(const MatType& x, const Sparsity& sp) {
      return MatType::reshape(x, sp);
    }

    /** \brief Cast matrix nonzeros to different Sparsity

        Keeps the nonzero values of \p x in storage order and gives them the pattern \p sp, which
        must have the same number of nonzeros (an error is raised otherwise). Unlike reshape, the
        dimensions and the positions of the nonzeros can change freely; only the nonzero count
        is preserved. This is how the nonzeros of a sparse matrix are extracted as a dense vector.

        \param x Matrix (or Sparsity) whose nonzeros are reused
        \param sp Target pattern with x.nnz() nonzeros
        \return Matrix of the same type as \p x with pattern \p sp

        \lang{python}
        \doctest
        >>> sparsity_cast(DM([[1, 2], [3, 4]]), Sparsity.dense(1, 4))
        DM([[1, 3, 2, 4]])
        >>> sparsity_cast(DM.eye(3), Sparsity.dense(3, 1))
        DM([1, 1, 1])
        >>> sparsity_cast(MX.sym("x", Sparsity.diag(3)), Sparsity.dense(3, 1))
        MX(nonzeros(x))
        >>> sparsity_cast(MX.sym("x", 2, 2), Sparsity.dense(3, 1))
        Traceback (most recent call last):
          ...
        RuntimeError: ...Mismatching nonzero count...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sparsity_cast(DM([1 2; 3 4]), Sparsity.dense(1, 4))
        [[1, 3, 2, 4]]
        >> sparsity_cast(DM.eye(3), Sparsity.dense(3, 1))
        [1, 1, 1]
        >> sparsity_cast(MX.sym('x', Sparsity.diag(3)), Sparsity.dense(3, 1))
        nonzeros(x)
        >> x = MX.sym('x', 2, 2);
        >> try, sparsity_cast(x, Sparsity.dense(3, 1)); catch e, disp(e.message); end
        ...Mismatching nonzero count...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.sparsity_cast(ca.DM([1 2; 3 4]), ca.Sparsity.dense(1, 4))
        DM([[1, 3, 2, 4]])
        julia> ca.sparsity_cast(ca.DM.eye(3), ca.Sparsity.dense(3, 1))
        DM([1, 1, 1])
        julia> ca.sparsity_cast(ca.MX.sym("x", ca.Sparsity.diag(3)), ca.Sparsity.dense(3, 1))
        MX(nonzeros(x))
        julia> try
                   ca.sparsity_cast(ca.MX.sym("x", 2, 2), ca.Sparsity.dense(3, 1))
               catch err
                   println(occursin("Mismatching nonzero count", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa reshape

        \identifier{24z} */
    inline friend MatType sparsity_cast(const MatType& x, const Sparsity& sp) {
      return MatType::sparsity_cast(x, sp);
    }

    /** \brief Obtain the structural rank of a sparsity-pattern

        The structural rank is the size of a maximum matching between rows and columns of the
        pattern (via the Dulmage-Mendelsohn decomposition): an upper bound on the numerical rank
        that depends only on which entries are structurally nonzero, not on their values. A
        dense matrix always has full structural rank, even when it is numerically singular.

        \param x Matrix (or Sparsity); only x.sparsity() is used
        \return Structural rank, at most min(size1, size2)

        \lang{python}
        \doctest
        >>> sprank(Sparsity.diag(3))
        3
        >>> sprank(Sparsity.triplet(3, 3, [0, 1, 2], [0, 0, 0]))
        1
        >>> sprank(DM([[1, 2], [2, 4]]))
        2
        >>> sprank(DM(2, 3))
        0
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sprank(Sparsity.diag(3))
        3
        >> sprank(Sparsity.triplet(3, 3, [0, 1, 2], [0, 0, 0]))
        1
        >> sprank(DM([1 2; 2 4]))
        2
        >> sprank(DM(2, 3))
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.sprank(ca.Sparsity.diag(3))
        3
        julia> ca.sprank(ca.Sparsity.triplet(3, 3, [0, 1, 2], [0, 0, 0]))
        1
        julia> ca.sprank(ca.DM([1 2; 2 4]))
        2
        julia> ca.sprank(ca.DM(2, 3))
        0
        \enddoctest
        \endlang

        \identifier{43} */
    inline friend casadi_int sprank(const MatType& x) {
      return MatType::sprank(x);
    }

    /** \brief 0-norm (nonzero count) of a Matrix-matrix product

        Number of structural nonzeros that mtimes(x, y) would have, computed from the sparsity
        patterns without forming the product. The inner dimensions must match, otherwise an
        error is raised.

        \lang{python}
        \doctest
        >>> norm_0_mul(DM.eye(3), DM.ones(3, 2))
        6
        >>> norm_0_mul(Sparsity.diag(3), Sparsity.diag(3))
        3
        >>> norm_0_mul(DM.ones(2, 3), DM.ones(2, 3))
        Traceback (most recent call last):
          ...
        RuntimeError: ...Dimension error...
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> norm_0_mul(DM.eye(3), DM.ones(3, 2))
        6
        >> norm_0_mul(Sparsity.diag(3), Sparsity.diag(3))
        3
        >> try, norm_0_mul(DM.ones(2, 3), DM.ones(2, 3)); catch e, disp(e.message); end
        ...Dimension error...
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.norm_0_mul(ca.DM.eye(3), ca.DM.ones(3, 2))
        6
        julia> ca.norm_0_mul(ca.Sparsity.diag(3), ca.Sparsity.diag(3))
        3
        julia> try
                   ca.norm_0_mul(ca.DM.ones(2, 3), ca.DM.ones(2, 3))
               catch err
                   println(occursin("Dimension error", sprint(showerror, err)))
               end
        true
        \enddoctest
        \endlang

        \sa mtimes

        \identifier{44} */
    inline friend casadi_int norm_0_mul(const MatType &x, const MatType &y) {
      return MatType::norm_0_mul(x, y);
    }

    /** \brief Get the upper triangular part of a matrix

        Entries below the diagonal become structural zeros; the matrix need not be square.
        For MX this is a projection node.

        \param x Matrix (or Sparsity); any shape
        \param includeDiagonal False to drop the diagonal as well

        \lang{python}
        \doctest
        >>> x = DM([[1, 2, 3], [4, 5, 6]])
        >>> triu(x)
        DM(
        [[1, 2, 3],
         [00, 5, 6]])
        >>> triu(x, False)
        DM(
        [[00, 2, 3],
         [00, 00, 6]])
        >>> triu(MX.sym("x", 2, 2))
        MX(project(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = DM([1 2 3; 4 5 6]);
        >> triu(x)
        [[1, 2, 3],
         [00, 5, 6]]
        >> triu(x, false)
        [[00, 2, 3],
         [00, 00, 6]]
        >> triu(MX.sym('x', 2, 2))
        project(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.DM([1 2 3; 4 5 6]);
        julia> ca.triu(x)
        DM(
        [[1, 2, 3],
         [00, 5, 6]])
        julia> ca.triu(x, false)
        DM(
        [[00, 2, 3],
         [00, 00, 6]])
        julia> ca.triu(ca.MX.sym("x", 2, 2))
        MX(project(x))
        \enddoctest
        \endlang

        \sa tril

        \identifier{45} */
    inline friend MatType triu(const MatType& x, bool includeDiagonal=true) {
      return MatType::triu(x, includeDiagonal);
    }

    /** \brief Get the lower triangular part of a matrix

        Entries above the diagonal become structural zeros; the matrix need not be square.

        \param x Matrix (or Sparsity); any shape
        \param includeDiagonal False to drop the diagonal as well

        \lang{python}
        \doctest
        >>> x = DM([[1, 2, 3], [4, 5, 6]])
        >>> tril(x)
        DM(
        [[1, 00, 00],
         [4, 5, 00]])
        >>> tril(x, False)
        DM(
        [[00, 00, 00],
         [4, 00, 00]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = DM([1 2 3; 4 5 6]);
        >> tril(x)
        [[1, 00, 00],
         [4, 5, 00]]
        >> tril(x, false)
        [[00, 00, 00],
         [4, 00, 00]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.DM([1 2 3; 4 5 6]);
        julia> ca.tril(x)
        DM(
        [[1, 00, 00],
         [4, 5, 00]])
        julia> ca.tril(x, false)
        DM(
        [[00, 00, 00],
         [4, 00, 00]])
        \enddoctest
        \endlang

        \sa triu

        \identifier{46} */
    inline friend MatType tril(const MatType& x, bool includeDiagonal=true) {
      return MatType::tril(x, includeDiagonal);
    }

    /** \brief Kronecker tensor product

        Creates a block matrix in which block (i, j) is a_ij*b, for any shapes of \p a and \p b.
        The result is a.size1()*b.size1() by a.size2()*b.size2() with the Kronecker product of
        the sparsity patterns.

        \lang{python}
        \doctest
        >>> kron(DM([[1, 2]]), DM([[1], [10]]))
        DM(
        [[1, 2],
         [10, 20]])
        >>> kron(DM.eye(2), DM([[1, 2], [3, 4]]))
        DM(
        [[1, 2, 00, 00],
         [3, 4, 00, 00],
         [00, 00, 1, 2],
         [00, 00, 3, 4]])
        >>> kron(MX.sym("a", 2, 2), MX.sym("b", 1, 2)).shape
        (2, 4)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> kron(DM([1 2]), DM([1; 10]))
        [[1, 2],
         [10, 20]]
        >> kron(DM.eye(2), DM([1 2; 3 4]))
        [[1, 2, 00, 00],
         [3, 4, 00, 00],
         [00, 00, 1, 2],
         [00, 00, 3, 4]]
        >> size(kron(MX.sym('a', 2, 2), MX.sym('b', 1, 2)))
        2 4
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> kron(ca.DM([1 2]), ca.DM([1; 10]))
        DM(
        [[1, 2],
         [10, 20]])
        julia> kron(ca.DM.eye(2), ca.DM([1 2; 3 4]))
        DM(
        [[1, 2, 00, 00],
         [3, 4, 00, 00],
         [00, 00, 1, 2],
         [00, 00, 3, 4]])
        julia> size(kron(ca.MX.sym("a", 2, 2), ca.MX.sym("b", 1, 2)))
        (2, 4)
        \enddoctest
        \endlang

        \identifier{47} */
    inline friend MatType kron(const MatType& a, const MatType& b) {
      return MatType::kron(a, b);
    }

    /** \brief Repeat matrix A n times vertically and m times horizontally

        Tiles \p A into an n-by-m grid of copies; the sparsity pattern is tiled along. A zero
        count gives an empty matrix with the remaining dimension kept. For MX the horizontal
        repetition is a single node, the vertical one a vertcat of \p n references.

        \param A Matrix (or Sparsity) to tile; any shape
        \param n Number of vertical repetitions, non-negative
        \param m Number of horizontal repetitions, non-negative

        \lang{python}
        \doctest
        >>> repmat(DM([[1, 2]]), 2, 2)
        DM(
        [[1, 2, 1, 2],
         [1, 2, 1, 2]])
        >>> repmat(DM.eye(2), 1, 2)
        DM(
        [[1, 00, 1, 00],
         [00, 1, 00, 1]])
        >>> repmat(MX.sym("x", 2), 1, 3)
        MX(repmat(x, 3))
        >>> repmat(MX.sym("x", 2), 3, 1)
        MX(vertcat(x, x, x))
        >>> repmat(DM.ones(2, 3), 0, 2).shape
        (0, 6)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> repmat(DM([1 2]), 2, 2)
        [[1, 2, 1, 2],
         [1, 2, 1, 2]]
        >> repmat(DM.eye(2), 1, 2)
        [[1, 00, 1, 00],
         [00, 1, 00, 1]]
        >> repmat(MX.sym('x', 2), 1, 3)
        repmat(x, 3)
        >> repmat(MX.sym('x', 2), 3, 1)
        vertcat(x, x, x)
        >> size(repmat(DM.ones(2, 3), 0, 2))
        0 6
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> repeat(ca.DM([1 2]), 2, 2)
        DM(
        [[1, 2, 1, 2],
         [1, 2, 1, 2]])
        julia> repeat(ca.DM.eye(2), 1, 2)
        DM(
        [[1, 00, 1, 00],
         [00, 1, 00, 1]])
        julia> repeat(ca.MX.sym("x", 2), 1, 3)
        MX(repmat(x, 3))
        julia> repeat(ca.MX.sym("x", 2), 3, 1)
        MX(vertcat(x, x, x))
        julia> size(repeat(ca.DM.ones(2, 3), 0, 2))
        (0, 6)
        \enddoctest
        \endlang

        \identifier{48} */
    inline friend MatType repmat(const MatType &A, casadi_int n, casadi_int m=1) {
      return MatType::repmat(A, n, m);
    }

    /** \brief Repeat matrix A n times vertically and m times horizontally

        Same as repmat(A, rc[0], rc[1]).

        \lang{python}
        \doctest
        >>> repmat(DM(7), (1, 3))
        DM([[7, 7, 7]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> repmat(DM(7), [1, 3])
        [[7, 7, 7]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.repmat(ca.DM(7), 1, 3)
        DM([[7, 7, 7]])
        \enddoctest
        \endlang

        \identifier{49} */
    inline friend MatType repmat(const MatType &A, const std::pair<casadi_int, casadi_int>& rc) {
      return MatType::repmat(A, rc);
    }

    /** \brief Concatenate horizontally, two matrices

        Same as the list form of horzcat.

        \identifier{4a} */
    inline friend MatType horzcat(const MatType &x, const MatType &y) {
      return horzcat(std::vector<MatType>{x, y});
    }

    /** \brief Concatenate horizontally, three matrices

        Same as the list form of horzcat.

        \identifier{4b} */
    inline friend MatType horzcat(const MatType &x, const MatType &y, const MatType &z) {
      return horzcat(std::vector<MatType>{x, y, z});
    }

    /** \brief Concatenate horizontally, four matrices

        Same as the list form of horzcat.

        \identifier{4c} */
    inline friend MatType horzcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w) {
      return horzcat(std::vector<MatType>{x, y, z, w});
    }

    /** \brief Concatenate horizontally, five matrices

        Same as the list form of horzcat.

        \identifier{4d} */
    inline friend MatType horzcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w, const MatType &v) {
      return horzcat(std::vector<MatType>{x, y, z, w, v});
    }

    /** \brief Concatenate horizontally, six matrices

        Same as the list form of horzcat.

        \identifier{4e} */
    inline friend MatType horzcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w, const MatType &v, const MatType &u) {
      return horzcat(std::vector<MatType>{x, y, z, w, v, u});
    }

    /** \brief Concatenate vertically, two matrices

        Same as the list form of vertcat.

        \identifier{4f} */
    inline friend MatType vertcat(const MatType &x, const MatType &y) {
      return vertcat(std::vector<MatType>{x, y});
    }

    /** \brief Concatenate vertically, three matrices

        Same as the list form of vertcat.

        \identifier{4g} */
    inline friend MatType vertcat(const MatType &x, const MatType &y, const MatType &z) {
      return vertcat(std::vector<MatType>{x, y, z});
    }

    /** \brief Concatenate vertically, four matrices

        Same as the list form of vertcat.

        \identifier{4h} */
    inline friend MatType vertcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w) {
      return vertcat(std::vector<MatType>{x, y, z, w});
    }

    /** \brief Concatenate vertically, five matrices

        Same as the list form of vertcat.

        \identifier{4i} */
    inline friend MatType vertcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w, const MatType &v) {
      return vertcat(std::vector<MatType>{x, y, z, w, v});
    }

    /** \brief Concatenate vertically, six matrices

        Same as the list form of vertcat.

        \identifier{4j} */
    inline friend MatType vertcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w, const MatType &v, const MatType &u) {
      return vertcat(std::vector<MatType>{x, y, z, w, v, u});
    }

    /** \brief Concatenate along diagonal, two matrices

        Same as the list form of diagcat.

        \identifier{4k} */
    inline friend MatType diagcat(const MatType &x, const MatType &y) {
      return diagcat(std::vector<MatType>{x, y});
    }

    /** \brief Concatenate along diagonal, three matrices

        Same as the list form of diagcat.

        \identifier{4l} */
    inline friend MatType diagcat(const MatType &x, const MatType &y, const MatType &z) {
      return diagcat(std::vector<MatType>{x, y, z});
    }

    /** \brief Concatenate along diagonal, four matrices

        Same as the list form of diagcat.

        \identifier{4m} */
    inline friend MatType diagcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w) {
      return diagcat(std::vector<MatType>{x, y, z, w});
    }

    /** \brief Concatenate along diagonal, five matrices

        Same as the list form of diagcat.

        \identifier{4n} */
    inline friend MatType diagcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w, const MatType &v) {
      return diagcat(std::vector<MatType>{x, y, z, w, v});
    }

    /** \brief Concatenate along diagonal, six matrices

        Same as the list form of diagcat.

        \identifier{4o} */
    inline friend MatType diagcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w, const MatType &v, const MatType &u) {
      return diagcat(std::vector<MatType>{x, y, z, w, v, u});
    }

    /** \brief Sum the elements of every column

        Returns a row vector with the sum over the rows, as MATLAB sum(x, 1); implemented as the
        product ones(1, size1) * x, so it is a mtimes node for MX. A column of structural zeros
        sums to a structural zero. For Sparsity only the pattern is computed.

        \lang{python}
        \doctest
        >>> sum1(DM([[1, 2, 3], [4, 5, 6]]))
        DM([[5, 7, 9]])
        >>> sum1(SX.sym("x", 2, 2))
        SX([[(x_0+x_1), (x_2+x_3)]])
        >>> sum1(MX.sym("x", 2, 3))
        MX(mac(ones(1x2),x,zeros(1x3)))
        >>> sum1(DM(0, 3))
        DM([[00, 00, 00]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sum1(DM([1 2 3; 4 5 6]))
        [[5, 7, 9]]
        >> sum(DM([1 2 3; 4 5 6]), 1)
        [[5, 7, 9]]
        >> sum1(SX.sym('x', 2, 2))
        [[(x_0+x_1), (x_2+x_3)]]
        >> sum1(MX.sym('x', 2, 3))
        mac(ones(1x2),x,zeros(1x3))
        >> sum1(DM(0, 3))
        [[00, 00, 00]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.sum1(ca.DM([1 2 3; 4 5 6]))
        DM([[5, 7, 9]])
        julia> ca.sum1(ca.SX.sym("x", 2, 2))
        SX([[(x_0+x_1), (x_2+x_3)]])
        julia> ca.sum1(ca.MX.sym("x", 2, 3))
        MX(mac(ones(1x2),x,zeros(1x3)))
        julia> ca.sum1(ca.DM(0, 3))
        DM([[00, 00, 00]])
        \enddoctest
        \endlang

        \sa sum2, sum

        \identifier{4p} */
    inline friend MatType sum1(const MatType &x) { return MatType::sum1(x);}

    /** \brief Sum the elements of every row

        Returns a column vector with the sum over the columns, as MATLAB sum(x, 2); implemented
        as the product x * ones(size2, 1). A row of structural zeros sums to a structural zero.

        \lang{python}
        \doctest
        >>> sum2(DM([[1, 2, 3], [4, 5, 6]]))
        DM([6, 15])
        >>> sum2(SX.sym("x", 2, 2))
        SX([(x_0+x_2), (x_1+x_3)])
        >>> sum2(DM(3, 3))
        DM([00, 00, 00])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sum2(DM([1 2 3; 4 5 6]))
        [6, 15]
        >> sum(DM([1 2 3; 4 5 6]), 2)
        [6, 15]
        >> sum2(SX.sym('x', 2, 2))
        [(x_0+x_2), (x_1+x_3)]
        >> sum2(DM(3, 3))
        [00, 00, 00]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.sum2(ca.DM([1 2 3; 4 5 6]))
        DM([6, 15])
        julia> ca.sum2(ca.SX.sym("x", 2, 2))
        SX([(x_0+x_2), (x_1+x_3)])
        julia> ca.sum2(ca.DM(3, 3))
        DM([00, 00, 00])
        \enddoctest
        \endlang

        \sa sum1, sum

        \identifier{4q} */
    inline friend MatType sum2(const MatType &x) { return MatType::sum2(x);}

    /** \brief Sum all elements

        Returns a 1-by-1 matrix with the sum of all elements, combining sum1 and sum2. An empty
        matrix sums to a structural zero.

        \lang{python}
        sum also accepts an axis argument: sum(x, 0) is sum1(x) and sum(x, 1) is sum2(x).

        \doctest
        >>> sum(DM([[1, 2, 3], [4, 5, 6]]))
        DM(21)
        >>> sum(DM([[1, 2, 3], [4, 5, 6]]), 0)
        DM([[5, 7, 9]])
        >>> sum(SX.sym("x", 3))
        SX(((x_0+x_1)+x_2))
        >>> sum(DM(3, 0))
        DM(00)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sum(DM([1 2 3; 4 5 6]), 'all')
        21
        >> sum(DM([1 2 3; 4 5 6]))
        [[5, 7, 9]]
        >> sum(SX.sym('x', 3))
        ((x_0+x_1)+x_2)
        >> sum(DM(3, 0), 'all')
        00
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sum(ca.DM([1 2 3; 4 5 6]))
        DM(21)
        julia> sum(ca.DM([1 2 3; 4 5 6]), dims=1)
        DM([[5, 7, 9]])
        julia> sum(ca.SX.sym("x", 3))
        SX(((x_0+x_1)+x_2))
        julia> sum(ca.DM(3, 0))
        DM(00)
        \enddoctest
        \endlang

        \sa sum1, sum2

        \identifier{2d6} */
    inline friend MatType sum(const MatType &x) {
      if (x.is_vector()) {
        if (x.is_column()) {
          return MatType::sum1(x);
        } else {
          return MatType::sum2(x);
        }
      }
      if (x.size1()>x.size2()) {
        return MatType::sum2(MatType::sum1(x));
      } else {
        return MatType::sum1(MatType::sum2(x));
      }
    }

    /** \@} */
  };
#endif // SWIG

#ifndef SWIG
  template<typename MatType>
  MatType SparsityInterface<MatType>::vec(const MatType& x) {
    if (x.size2()==1) return x;
    return reshape(x, x.numel(), 1);
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::repmat(const MatType& x, casadi_int n, casadi_int m) {
    casadi_assert(n>=0 && m>=0,
      "repmat(x,n,m): n (" + str(n) + ") and m (" + str(m) + ") must be non-negative");
    if (n==1 && m==1) return x;
    MatType allrows = vertcat(std::vector<MatType>(n, x));
    if (n==0) allrows = MatType(0, x.size2());
    MatType ret = horzcat(std::vector<MatType>(m, allrows));
    if (m==0) ret = MatType(allrows.size1(), 0);
    return ret;
  }

  template<typename MatType>
  std::vector< std::vector< MatType > >
  SparsityInterface<MatType>::blocksplit(const MatType& x,
                                         const std::vector<casadi_int>& vert_offset,
                                         const std::vector<casadi_int>& horz_offset) {
    std::vector<MatType> rows = MatType::vertsplit(x, vert_offset);
    std::vector< std::vector< MatType > > ret;
    for (auto&& r : rows) ret.push_back(MatType::horzsplit(r, horz_offset));
    return ret;
  }

  template<typename MatType>
  std::vector< std::vector< MatType > >
  SparsityInterface<MatType>::blocksplit(const MatType& x,
      casadi_int vert_incr, casadi_int horz_incr) {
    casadi_assert(vert_incr>=1, "blocksplit: vert_incr (" + str(vert_incr) + ") must be positive");
    casadi_assert(horz_incr>=1, "blocksplit: horz_incr (" + str(horz_incr) + ") must be positive");
    casadi_int sz1 = x.size1();
    std::vector<casadi_int> offset1 = range(0, sz1, vert_incr);
    offset1.push_back(sz1);
    casadi_int sz2 = x.size2();
    std::vector<casadi_int> offset2 = range(0, sz2, horz_incr);
    offset2.push_back(sz2);
    return blocksplit(x, offset1, offset2);
  }

  template<typename MatType>
  std::vector<casadi_int>
  SparsityInterface<MatType>::offset(const std::vector< MatType > &v, bool vert) {
    std::vector<casadi_int> ret(v.size()+1);
    ret[0]=0;
    for (casadi_int i=0; i<v.size(); ++i) {
      ret[i+1] = ret[i] + (vert ? v[i].size1() : v[i].size2());
    }
    return ret;
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::veccat(const std::vector< MatType >& x) {
    std::vector< MatType > x_vec = x;
    for (typename std::vector< MatType >::iterator it=x_vec.begin();
         it!=x_vec.end(); ++it) {
      *it = vec(*it);
    }
    if (x_vec.empty()) {
      return MatType(0, 1);
    } else {
      return vertcat(x_vec);
    }
  }

  template<typename MatType>
  std::vector< MatType >
  SparsityInterface<MatType>::diagsplit(const MatType& x,
      const std::vector<casadi_int>& output_offset) {
    casadi_assert(x.is_square(), "diagsplit(x,incr)::input must be square but got "
                          + x.dim()  + ".");
    return MatType::diagsplit(x, output_offset, output_offset);
  }

  template<typename MatType>
  std::vector< MatType >
  SparsityInterface<MatType>::diagsplit(const MatType& x, casadi_int incr) {
    casadi_assert(incr>=1, "diagsplit(x,incr): incr (" + str(incr) + ") must be positive");
    casadi_assert(x.is_square(), "diagsplit(x,incr)::input must be square but got "
                          + x.dim()  + ".");
    std::vector<casadi_int> offset2 = range(0, x.size2(), incr);
    offset2.push_back(x.size2());
    return MatType::diagsplit(x, offset2);
  }

  template<typename MatType>
  std::vector< MatType >
  SparsityInterface<MatType>::diagsplit(const MatType& x, casadi_int incr1, casadi_int incr2) {
    casadi_assert(incr1>=1, "diagsplit: incr1 (" + str(incr1) + ") must be positive");
    casadi_assert(incr2>=1, "diagsplit: incr2 (" + str(incr2) + ") must be positive");
    std::vector<casadi_int> offset1 = range(0, x.size1(), incr1);
    offset1.push_back(x.size1());
    std::vector<casadi_int> offset2 = range(0, x.size2(), incr2);
    offset2.push_back(x.size2());
    return MatType::diagsplit(x, offset1, offset2);
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::mtimes(const std::vector<MatType> &args,
                                             const std::string& blas) {
    casadi_assert(!args.empty(),
                          "mul(std::vector<MatType> &args): "
                          "supplied list must not be empty.");
    MatType ret = args[0];
    for (casadi_int i=1; i<args.size(); ++i) ret = MatType::mtimes(ret, args[i], blas);
    return ret;
  }

  template<typename MatType>
  std::vector<MatType > SparsityInterface<MatType>::horzsplit(const MatType& x, casadi_int incr) {
    casadi_assert(incr>=1, "horzsplit(x,incr): incr (" + str(incr) + ") must be positive");
    casadi_int sz2 = x.size2();
    std::vector<casadi_int> offset2 = range(0, sz2, incr);
    offset2.push_back(sz2);
    return MatType::horzsplit(x, offset2);
  }

  template<typename MatType>
  std::vector<MatType > SparsityInterface<MatType>::vertsplit(const MatType& x, casadi_int incr) {
    casadi_assert(incr>=1, "vertsplit(x,incr): incr (" + str(incr) + ") must be positive");
    casadi_int sz1 = x.size1();
    std::vector<casadi_int> offset1 = range(0, sz1, incr);
    offset1.push_back(sz1);
    return MatType::vertsplit(x, offset1);
  }
  template<typename MatType>
  std::vector<MatType > SparsityInterface<MatType>::horzsplit_n(const MatType& x, casadi_int n) {
    casadi_assert(n>=0, "horzsplit_n(x,n): n (" + str(n) + ") must be non-negative");
    if (x.size2()==0) return std::vector<MatType>(n, x);
    casadi_assert(n>0 && x.size2() % n==0, "horzsplit_n(x,n): x.size2() (" + str(x.size2()) +
                     ") must be a multiple of n (" + str(n) + ")");
    return horzsplit(x, x.size2()/n);
  }
  template<typename MatType>
  std::vector<MatType > SparsityInterface<MatType>::vertsplit_n(const MatType& x, casadi_int n) {
    casadi_assert(n>=0, "vertsplit_n(x,n): n (" + str(n) + ") must be non-negative");
    if (x.size1()==0) return std::vector<MatType>(n, x);
    casadi_assert(n>0 && x.size1() % n==0, "vertsplit_n(x,n): x.size1() (" + str(x.size1()) +
                     ") must be a multiple of n (" + str(n) + ")");
    return vertsplit(x, x.size1()/n);
  }

#endif // SWIG

} // namespace casadi

#endif // CASADI_SPARSITY_INTERFACE_HPP
