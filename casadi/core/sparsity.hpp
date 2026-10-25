/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *    Copyright (C) 2005-2013 Timothy A. Davis
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


#ifndef CASADI_SPARSITY_HPP
#define CASADI_SPARSITY_HPP

#include "shared_object.hpp"
#include "printable.hpp"
#include "casadi_common.hpp"
#include "sparsity_interface.hpp"
#include "generic_type.hpp"
#include <vector>
#include <list>
#include <limits>
#include <unordered_map>
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

namespace casadi {
  // Forward declaration
  class SparsityInternal;
  class SerializingStream;
  class DeserializingStream;

  #ifndef SWIG
    /** \brief Compact representation of a sparsity pattern

        \identifier{b8} */
    struct CASADI_EXPORT SparsityStruct {
      casadi_int nrow;
      casadi_int ncol;
      const casadi_int* colind;
      const casadi_int* row;
    };
  #endif // SWIG

  /** \brief General sparsity class
   *
   * The storage format is a compressed column storage (CCS) format.\n
   *
   In this format, the structural non-zero elements are stored in column-major order, starting from
   the upper left corner of the matrix and ending in the lower right corner.

   In addition to the dimension (size1(), size2()), (i.e. the number of rows and the number of
   columns respectively), there are also two vectors of integers:

   1. "colind" [length size2()+1], which contains the index to the first non-zero element on or
   after the corresponding column. All the non-zero elements of a particular column i are thus
   the elements with index el that fulfills: colind[i] <= el < colind[i+1].

   2. "row" [same length as the number of non-zero elements, nnz()] The rows for each of the
   structural non-zeros.

   Note that with this format, it is cheap to loop over all the non-zero elements of a particular
   column, at constant time per element, but expensive to jump to access a location (i, j).

   If the matrix is dense, i.e. length(row) == size1()*size2(), the format reduces to standard dense
   column major format, which allows access to an arbitrary element in constant time.

   Since the object is reference counted (it inherits from SharedObject), several matrices are
   allowed to share the same sparsity pattern.

   The implementations of methods marked as such in this class has been taken from the
   CSparse package and modified to fit CasADi data structures and separation of
   sparsity pattern calculation and numerical evaluation.
   These functions are Copyright(c) Timothy A. Davis, 2006-2009
   and licensed as a derivative work under the GNU LGPL

   Every DM, SX and MX carries a pattern, obtained with sparsity(), and a symbolic or numeric
   matrix can be created from one, e.g. SX.sym("x", sp) or DM(sp, 1). Entries outside the
   pattern are structural zeros: they are never stored and never computed with. The row indices
   within a column are strictly increasing.

   Patterns are immutable and cached: constructing a pattern that already exists returns a
   reference to the cached object, so equal patterns are cheap to compare. The in-place methods
   (resize, enlarge, erase, append, add_nz, ...) rebind the object to a new pattern rather than
   modify the shared one. Use is_equal or the == operator to compare patterns; object identity
   is not meaningful in the target languages.

   A default-constructed Sparsity is a null object (is_null() is true): it carries no pattern
   and most methods raise an error on it. It is used as "no pattern" sentinel, e.g. as the
   result of a coloring that exceeds its cutoff. Empty patterns are allowed: a 0-by-n or
   n-by-0 pattern has no nonzeros, is dense, and is_empty() is true for it.

      \lang{python}
      \doctest
      >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
      >>> DM(sp, 1)
      DM(
      [[1, 00, 00, 1],
       [00, 00, 1, 00],
       [1, 00, 00, 00]])
      >>> sp
      Sparsity(3x4,4nz)
      >>> sp.size1(), sp.size2(), sp.nnz(), sp.numel()
      (3, 4, 4, 12)
      >>> sp.colind(), sp.row()
      ([0, 2, 2, 3, 4], [0, 2, 1, 0])
      >>> sp == Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
      True
      >>> SX.sym("x", sp)
      SX(
      [[x_0, 00, 00, x_3],
       [00, 00, x_2, 00],
       [x_1, 00, 00, 00]])
      >>> Sparsity.dense(0, 3), Sparsity().is_null()
      (Sparsity(0x3), True)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
      >> DM(sp, 1)
      [[1, 00, 00, 1],
       [00, 00, 1, 00],
       [1, 00, 00, 00]]
      >> sp
      3x4,4nz
      >> fprintf('%d %d %d %d\n', sp.size1(), sp.size2(), sp.nnz(), sp.numel())
      3 4 4 12
      >> sp.colind(), sp.row()
      0 2 2 3 4
      0 2 1 0
      >> sp == Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
      1
      >> SX.sym('x', sp)
      [[x_0, 00, 00, x_3],
       [00, 00, x_2, 00],
       [x_1, 00, 00, 00]]
      >> Sparsity.dense(0, 3), Sparsity().is_null()
      0x3
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
      julia> ca.DM(sp, 1)
      DM(
      [[1, 00, 00, 1],
       [00, 00, 1, 00],
       [1, 00, 00, 00]])
      julia> sp
      Sparsity(3x4,4nz)
      julia> ca.size1(sp), ca.size2(sp), ca.nnz(sp), ca.numel(sp)
      (3, 4, 4, 12)
      julia> ca.colind(sp), ca.row(sp)
      ([0, 2, 2, 3, 4], [0, 2, 1, 0])
      julia> ca.is_equal(sp, ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]))
      true
      julia> ca.SX.sym("x", sp)
      SX(
      [[x_0, 00, 00, x_3],
       [00, 00, x_2, 00],
       [x_1, 00, 00, 00]])
      julia> ca.Sparsity.dense(0, 3), ca.is_null(ca.Sparsity())
      (0x3, true)
      \enddoctest
      \endlang

   * \see Matrix
   *
   * \author Joel Andersson
   * \date 2010-2015

      \identifier{b9} */
  class CASADI_EXPORT Sparsity
    : public SharedObject,
      public SWIG_IF_ELSE(SparsityInterfaceCommon, SparsityInterface<Sparsity>),
      public SWIG_IF_ELSE(PrintableCommon, Printable<Sparsity>) {
  public:

    /** \brief Default constructor: a null object that holds no pattern

        is_null() is true for the result and most other methods raise an error on it.
        A null Sparsity is the "no pattern" sentinel used by e.g. uni_coloring (as default
        for \p AT) and returned by the colorings when the cutoff is exceeded.

        \param dummy Must be 0 (any other value raises an error)

        \lang{python}
        \doctest
        >>> sp = Sparsity()
        >>> sp.is_null(), sp
        (True, Sparsity(NULL))
        >>> try: sp.nnz()
        ... except Exception as e: print("error")
        error
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity();
        >> sp.is_null(), sp
        1
        NULL
        >> try, sp.nnz(); catch e, disp('error'); end
        error
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity();
        julia> ca.is_null(sp), sp
        (true, NULL)
        julia> try
                   ca.nnz(sp)
               catch e
                   println("error")
               end
        error
        \enddoctest
        \endlang

        \identifier{2nj} */
    explicit Sparsity(casadi_int dummy=0);

    /** \brief Pattern of the given shape with all structural zeros

        \param nrow Number of rows, nonnegative
        \param ncol Number of columns, nonnegative

        \lang{python}
        \doctest
        >>> sp = Sparsity(2, 3)
        >>> sp, sp.nnz(), sp.is_dense()
        (Sparsity(2x3,0nz), 0, False)
        >>> Sparsity(0, 0), Sparsity(0, 0).is_dense()
        (Sparsity(0x0), True)
        >>> DM(Sparsity(2, 2), 1)
        DM(
        [[00, 00],
         [00, 00]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity(2, 3);
        >> sp, sp.nnz(), sp.is_dense()
        2x3,0nz
        0
        0
        >> Sparsity(0, 0), Sparsity(0, 0).is_dense()
        0x0
        1
        >> DM(Sparsity(2, 2), 1)
        [[00, 00],
         [00, 00]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity(2, 3)
        Sparsity(2x3,0nz)
        julia> ca.nnz(sp), ca.is_dense(sp)
        (0, false)
        julia> ca.Sparsity(0, 0), ca.is_dense(ca.Sparsity(0, 0))
        (0x0, true)
        julia> ca.DM(ca.Sparsity(2, 2), 1)
        DM(
        [[00, 00],
         [00, 00]])
        \enddoctest
        \endlang

        \sa dense, is_empty

        \identifier{ba} */
    Sparsity(casadi_int nrow, casadi_int ncol);

    /** \brief Construct from compressed column storage vectors

        \param nrow Number of rows
        \param ncol Number of columns
        \param colind Column offsets, length \p ncol+1, colind[0] == 0, nondecreasing
        \param row Row index of each nonzero, length colind[ncol], each in [0, nrow),
               strictly increasing within a column
        \param order_rows If true, row indices that are not sorted within a column are
               sorted (duplicates merged) instead of raising an error

        An invalid colind or an out-of-range row index raises an error. A pattern
        with colind[ncol] == nrow*ncol is the dense pattern.

        \lang{python}
        \doctest
        >>> sp = Sparsity(3, 2, [0, 2, 3], [0, 2, 1])
        >>> DM(sp, 1)
        DM(
        [[1, 00],
         [00, 1],
         [1, 00]])
        >>> Sparsity(2, 2, [0, 2, 4], [0, 1, 0, 1]).is_dense()
        True
        >>> try: Sparsity(3, 2, [0, 2, 3], [2, 0, 1])
        ... except Exception as e: print("not sorted")
        not sorted
        >>> Sparsity(3, 2, [0, 2, 3], [2, 0, 1], True).row()
        [0, 2, 1]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity(3, 2, [0, 2, 3], [0, 2, 1]);
        >> DM(sp, 1)
        [[1, 00],
         [00, 1],
         [1, 00]]
        >> Sparsity(2, 2, [0, 2, 4], [0, 1, 0, 1]).is_dense()
        1
        >> try, Sparsity(3, 2, [0, 2, 3], [2, 0, 1]); catch e, disp('not sorted'); end
        not sorted
        >> Sparsity(3, 2, [0, 2, 3], [2, 0, 1], true).row()
        0 2 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity(3, 2, [0, 2, 3], [0, 2, 1]);
        julia> ca.DM(sp, 1)
        DM(
        [[1, 00],
         [00, 1],
         [1, 00]])
        julia> ca.is_dense(ca.Sparsity(2, 2, [0, 2, 4], [0, 1, 0, 1]))
        true
        julia> try
                   ca.Sparsity(3, 2, [0, 2, 3], [2, 0, 1])
               catch e
                   println("not sorted")
               end
        not sorted
        julia> ca.row(ca.Sparsity(3, 2, [0, 2, 3], [2, 0, 1], true))
        3-element Vector{Int64}:
         0
         2
         1
        \enddoctest
        \endlang

        \sa compressed, triplet, get_ccs

        \identifier{2nk} */
    Sparsity(casadi_int nrow, casadi_int ncol,
             const std::vector<casadi_int>& colind, const std::vector<casadi_int>& row,
             bool order_rows=false);

    /** \brief Pattern with all structural zeros, shape given as a pair (nrow, ncol)

        \identifier{bb} */
    explicit Sparsity(const std::pair<casadi_int, casadi_int>& rc);

#ifndef SWIG
    /// Construct from sparsity pattern vectors given in compressed column storage format
    Sparsity(casadi_int nrow, casadi_int ncol, const casadi_int* colind, const casadi_int* row,
            bool order_rows=false);

    /** \brief  Create from node

        \identifier{bc} */
    static Sparsity create(SparsityInternal *node);

    /// Base class
    typedef SparsityInterface<Sparsity> B;

    /// Expose base class functions
    using B::horzsplit;
    using B::diagsplit;
    using B::vertsplit;
    using B::mtimes;

    SparsityInternal* get() const;
#endif

    /** \brief Create a 1-by-1 sparsity pattern

        \param dense_scalar If true (default), the single entry is a nonzero; if false,
               the pattern is 1-by-1 with zero nonzeros (a structurally zero scalar)

        \lang{python}
        \doctest
        >>> Sparsity.scalar(), Sparsity.scalar().nnz()
        (Sparsity(1x1), 1)
        >>> Sparsity.scalar(False), Sparsity.scalar(False).nnz()
        (Sparsity(1x1,0nz), 0)
        >>> Sparsity.scalar() == Sparsity.dense(1, 1)
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.scalar(), Sparsity.scalar().nnz()
        1x1
        1
        >> Sparsity.scalar(false), Sparsity.scalar(false).nnz()
        1x1,0nz
        0
        >> Sparsity.scalar() == Sparsity.dense(1, 1)
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.Sparsity.scalar(), ca.nnz(ca.Sparsity.scalar())
        (1x1, 1)
        julia> ca.Sparsity.scalar(false), ca.nnz(ca.Sparsity.scalar(false))
        (1x1,0nz, 0)
        julia> ca.is_equal(ca.Sparsity.scalar(), ca.Sparsity.dense(1, 1))
        true
        \enddoctest
        \endlang

        \sa dense, is_scalar

        \identifier{bd} */
    ///@{
    static Sparsity scalar(bool dense_scalar=true)
    { return dense_scalar ? dense(1, 1) : Sparsity(1, 1); }
    ///@}

    /** \brief Create a dense rectangular sparsity pattern

        All nrow*ncol entries are nonzeros, stored column by column. This is the pattern
        of every matrix created without an explicit sparsity, e.g. SX.sym("x", 2, 3).

        \param nrow Number of rows, nonnegative
        \param ncol Number of columns, nonnegative; default 1 (a dense column vector)

        The shape can also be given as a single pair (nrow, ncol).

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(2, 3)
        >>> sp, sp.nnz(), sp.is_dense()
        (Sparsity(2x3), 6, True)
        >>> sp.colind(), sp.row()
        ([0, 2, 4, 6], [0, 1, 0, 1, 0, 1])
        >>> Sparsity.dense(3), Sparsity.dense((2, 3))
        (Sparsity(3x1), Sparsity(2x3))
        >>> Sparsity.dense(0, 3).nnz(), Sparsity.dense(0, 3).is_dense()
        (0, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(2, 3);
        >> sp, sp.nnz(), sp.is_dense()
        2x3
        6
        1
        >> sp.colind(), sp.row()
        0 2 4 6
        0 1 0 1 0 1
        >> Sparsity.dense(3), Sparsity.dense([2, 3])
        3x1
        2x3
        >> Sparsity.dense(0, 3).nnz(), Sparsity.dense(0, 3).is_dense()
        0
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(2, 3)
        Sparsity(2x3)
        julia> ca.nnz(sp), ca.is_dense(sp)
        (6, true)
        julia> ca.colind(sp), ca.row(sp)
        ([0, 2, 4, 6], [0, 1, 0, 1, 0, 1])
        julia> ca.Sparsity.dense(3), ca.Sparsity.dense(2, 3)
        (3x1, 2x3)
        julia> ca.nnz(ca.Sparsity.dense(0, 3)), ca.is_dense(ca.Sparsity.dense(0, 3))
        (0, true)
        \enddoctest
        \endlang

        \sa Sparsity(nrow, ncol), scalar, is_dense

        \identifier{be} */
    ///@{
    static Sparsity dense(casadi_int nrow, casadi_int ncol=1);
    static Sparsity dense(const std::pair<casadi_int, casadi_int> &rc) {
      return dense(rc.first, rc.second);
    }
    ///@}

    /** \brief Create the pattern of a unit vector: a column of length n with one nonzero

        \param n Length of the column vector
        \param el Row of the nonzero, in [0, n); an index outside this range raises an error

        \lang{python}
        \doctest
        >>> DM(Sparsity.unit(4, 1), 1)
        DM([00, 1, 00, 00])
        >>> Sparsity.unit(4, 1).row()
        [1]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.unit(4, 1), 1)
        [00, 1, 00, 00]
        >> Sparsity.unit(4, 1).row()
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.unit(4, 1), 1)
        DM([00, 1, 00, 00])
        julia> ca.row(ca.Sparsity.unit(4, 1))
        1-element Vector{Int64}:
         1
        \enddoctest
        \endlang

        \sa nonzeros, triplet

        \identifier{bf} */
    ///@{
    static Sparsity unit(casadi_int n, casadi_int el);
    ///@}

    /** \brief Create an upper triangular square sparsity pattern

        Column j contains the rows 0, ..., j, so the diagonal is included and the pattern
        has n*(n+1)/2 nonzeros.

        \param n Size of the square pattern, nonnegative

        \lang{python}
        \doctest
        >>> DM(Sparsity.upper(3), 1)
        DM(
        [[1, 1, 1],
         [00, 1, 1],
         [00, 00, 1]])
        >>> Sparsity.upper(3).nnz(), Sparsity.upper(3).is_triu()
        (6, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.upper(3), 1)
        [[1, 1, 1],
         [00, 1, 1],
         [00, 00, 1]]
        >> Sparsity.upper(3).nnz(), Sparsity.upper(3).is_triu()
        6
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.upper(3), 1)
        DM(
        [[1, 1, 1],
         [00, 1, 1],
         [00, 00, 1]])
        julia> ca.nnz(ca.Sparsity.upper(3)), ca.is_triu(ca.Sparsity.upper(3))
        (6, true)
        \enddoctest
        \endlang

        \sa lower, triu, is_triu

        \identifier{bg} */
    static Sparsity upper(casadi_int n);

    /** \brief Create a lower triangular square sparsity pattern

        Column j contains the rows j, ..., n-1, so the diagonal is included and the pattern
        has n*(n+1)/2 nonzeros.

        \param n Size of the square pattern, nonnegative

        \lang{python}
        \doctest
        >>> DM(Sparsity.lower(3), 1)
        DM(
        [[1, 00, 00],
         [1, 1, 00],
         [1, 1, 1]])
        >>> SX.sym("x", Sparsity.lower(2))
        SX(
        [[x_0, 00],
         [x_1, x_2]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.lower(3), 1)
        [[1, 00, 00],
         [1, 1, 00],
         [1, 1, 1]]
        >> SX.sym('x', Sparsity.lower(2))
        [[x_0, 00],
         [x_1, x_2]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.lower(3), 1)
        DM(
        [[1, 00, 00],
         [1, 1, 00],
         [1, 1, 1]])
        julia> ca.SX.sym("x", ca.Sparsity.lower(2))
        SX(
        [[x_0, 00],
         [x_1, x_2]])
        \enddoctest
        \endlang

        \sa upper, tril, is_tril

        \identifier{bh} */
    static Sparsity lower(casadi_int n);

    /** \brief Create a diagonal sparsity pattern

        The nonzeros are the entries (i, i) for i < min(nrow, ncol); a rectangular
        diagonal pattern is allowed. The shape can be given as one integer (square),
        as two integers or as a pair.

        \lang{python}
        \doctest
        >>> DM(Sparsity.diag(3), 1)
        DM(
        [[1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]])
        >>> DM(Sparsity.diag(2, 4), 1)
        DM(
        [[1, 00, 00, 00],
         [00, 1, 00, 00]])
        >>> Sparsity.diag((4, 2)).nnz(), Sparsity.diag(0).nnz()
        (2, 0)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.diag(3), 1)
        [[1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]]
        >> DM(Sparsity.diag(2, 4), 1)
        [[1, 00, 00, 00],
         [00, 1, 00, 00]]
        >> Sparsity.diag([4, 2]).nnz(), Sparsity.diag(0).nnz()
        2
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.diag(3), 1)
        DM(
        [[1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]])
        julia> ca.DM(ca.Sparsity.diag(2, 4), 1)
        DM(
        [[1, 00, 00, 00],
         [00, 1, 00, 00]])
        julia> ca.nnz(ca.Sparsity.diag(4, 2)), ca.nnz(ca.Sparsity.diag(0))
        (2, 0)
        \enddoctest
        \endlang

        \sa band, is_diag, get_diag

        \identifier{bi} */
    ///@{
    static Sparsity diag(casadi_int nrow) { return diag(nrow, nrow);}
    static Sparsity diag(casadi_int nrow, casadi_int ncol);
    static Sparsity diag(const std::pair<casadi_int, casadi_int> &rc) {
      return diag(rc.first, rc.second);
    }
    ///@}

    /** \brief Create a single band in a square sparsity pattern

        The nonzeros are the entries (i+p, i): band(n, 0) equals diag(n), band(n, 1) is the
        first subdiagonal and band(n, -1) the first superdiagonal.

        \param n Size of the square pattern
        \param p Offset of the band from the diagonal, positive below the diagonal and
               negative above it; |p| must be smaller than \p n, otherwise an error is raised

        \lang{python}
        \doctest
        >>> DM(Sparsity.band(4, 1), 1)
        DM(
        [[00, 00, 00, 00],
         [1, 00, 00, 00],
         [00, 1, 00, 00],
         [00, 00, 1, 00]])
        >>> DM(Sparsity.band(4, -2), 1)
        DM(
        [[00, 00, 1, 00],
         [00, 00, 00, 1],
         [00, 00, 00, 00],
         [00, 00, 00, 00]])
        >>> Sparsity.band(4, 0) == Sparsity.diag(4)
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.band(4, 1), 1)
        [[00, 00, 00, 00],
         [1, 00, 00, 00],
         [00, 1, 00, 00],
         [00, 00, 1, 00]]
        >> DM(Sparsity.band(4, -2), 1)
        [[00, 00, 1, 00],
         [00, 00, 00, 1],
         [00, 00, 00, 00],
         [00, 00, 00, 00]]
        >> Sparsity.band(4, 0) == Sparsity.diag(4)
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.band(4, 1), 1)
        DM(
        [[00, 00, 00, 00],
         [1, 00, 00, 00],
         [00, 1, 00, 00],
         [00, 00, 1, 00]])
        julia> ca.DM(ca.Sparsity.band(4, -2), 1)
        DM(
        [[00, 00, 1, 00],
         [00, 00, 00, 1],
         [00, 00, 00, 00],
         [00, 00, 00, 00]])
        julia> ca.is_equal(ca.Sparsity.band(4, 0), ca.Sparsity.diag(4))
        true
        \enddoctest
        \endlang

        \sa banded, diag

        \identifier{bj} */
    static Sparsity band(casadi_int n, casadi_int p);

    /** \brief Create a banded square sparsity pattern

        Union of band(n, k) for all k in [-p, p]: banded(n, 0) equals diag(n) and banded(n, 1)
        is tridiagonal.

        \param n Size of the square pattern
        \param p Half-bandwidth; \p p >= n gives a dense pattern and a negative \p p an
               all-zero pattern

        \lang{python}
        \doctest
        >>> DM(Sparsity.banded(4, 1), 1)
        DM(
        [[1, 1, 00, 00],
         [1, 1, 1, 00],
         [00, 1, 1, 1],
         [00, 00, 1, 1]])
        >>> Sparsity.banded(4, 1).bw_upper(), Sparsity.banded(4, 1).bw_lower()
        (1, 1)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.banded(4, 1), 1)
        [[1, 1, 00, 00],
         [1, 1, 1, 00],
         [00, 1, 1, 1],
         [00, 00, 1, 1]]
        >> [Sparsity.banded(4, 1).bw_upper(), Sparsity.banded(4, 1).bw_lower()]
        1 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.banded(4, 1), 1)
        DM(
        [[1, 1, 00, 00],
         [1, 1, 1, 00],
         [00, 1, 1, 1],
         [00, 00, 1, 1]])
        julia> ca.bw_upper(ca.Sparsity.banded(4, 1)), ca.bw_lower(ca.Sparsity.banded(4, 1))
        (1, 1)
        \enddoctest
        \endlang

        \sa band, bw_upper, bw_lower

        \identifier{bk} */
    static Sparsity banded(casadi_int n, casadi_int p);

    /** \brief Create the pattern of a block: all entries (i, j) with i in \p row and j in \p col

        Every combination of a row index and a column index becomes a nonzero; the block
        need not be contiguous. Duplicate indices are merged.

        \param row Row indices, each in [0, nrow)
        \param col Column indices, each in [0, ncol)
        \param nrow Number of rows of the pattern
        \param ncol Number of columns of the pattern

        \lang{python}
        \doctest
        >>> DM(Sparsity.rowcol([0, 2], [1, 2], 3, 4), 1)
        DM(
        [[00, 1, 1, 00],
         [00, 00, 00, 00],
         [00, 1, 1, 00]])
        >>> Sparsity.rowcol([], [1], 3, 4).nnz()
        0
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.rowcol([0, 2], [1, 2], 3, 4), 1)
        [[00, 1, 1, 00],
         [00, 00, 00, 00],
         [00, 1, 1, 00]]
        >> Sparsity.rowcol([], [1], 3, 4).nnz()
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.rowcol([0, 2], [1, 2], 3, 4), 1)
        DM(
        [[00, 1, 1, 00],
         [00, 00, 00, 00],
         [00, 1, 1, 00]])
        julia> ca.nnz(ca.Sparsity.rowcol([], [1], 3, 4))
        0
        \enddoctest
        \endlang

        \sa triplet

        \identifier{bl} */
    static Sparsity rowcol(const std::vector<casadi_int>& row,
                           const std::vector<casadi_int>& col,
                           casadi_int nrow, casadi_int ncol);

    ///@{
    /** \brief Create a sparsity pattern from the (row, col) index of every nonzero

        The entries may be given in any order; duplicates are merged. An index out of
        range or vectors of different length raise an error.

        \param nrow Number of rows
        \param ncol Number of columns
        \param row Row index of each entry, in [0, nrow)
        \param col Column index of each entry, in [0, ncol), same length as \p row
        \param invert_mapping Selects the direction of the returned mapping (see below)

        The overload with \p invert_mapping additionally returns a mapping between the
        given entries and the nonzeros of the result (which are in column-major order).
        With \p invert_mapping false, mapping[k] is the position in the input of nonzero k
        (length nnz(); for duplicates the first occurrence). With \p invert_mapping true,
        mapping[i] is the nonzero index of input entry i (length of the input).

        \lang{python}
        \doctest
        >>> DM(Sparsity.triplet(3, 3, [0, 2, 1], [0, 0, 2]), 1)
        DM(
        [[1, 00, 00],
         [00, 00, 1],
         [1, 00, 00]])
        >>> Sparsity.triplet(2, 2, [0, 0, 1], [0, 0, 1]).nnz()
        2
        >>> sp, mapping = Sparsity.triplet(3, 3, [2, 0, 1], [1, 0, 0], False)
        >>> mapping
        [1, 2, 0]
        >>> sp, mapping = Sparsity.triplet(3, 3, [2, 0, 1], [1, 0, 0], True)
        >>> mapping
        [2, 0, 1]
        >>> try: Sparsity.triplet(2, 2, [0, 2], [0, 0])
        ... except Exception as e: print("out of bounds")
        out of bounds
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.triplet(3, 3, [0, 2, 1], [0, 0, 2]), 1)
        [[1, 00, 00],
         [00, 00, 1],
         [1, 00, 00]]
        >> Sparsity.triplet(2, 2, [0, 0, 1], [0, 0, 1]).nnz()
        2
        >> [sp, mapping] = Sparsity.triplet(3, 3, [2, 0, 1], [1, 0, 0], false); disp(mapping)
        1 2 0
        >> [sp, mapping] = Sparsity.triplet(3, 3, [2, 0, 1], [1, 0, 0], true); disp(mapping)
        2 0 1
        >> try, Sparsity.triplet(2, 2, [0, 2], [0, 0]); catch e, disp('out of bounds'); end
        out of bounds
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.triplet(3, 3, [0, 2, 1], [0, 0, 2]), 1)
        DM(
        [[1, 00, 00],
         [00, 00, 1],
         [1, 00, 00]])
        julia> ca.nnz(ca.Sparsity.triplet(2, 2, [0, 0, 1], [0, 0, 1]))
        2
        julia> sp, mapping = ca.Sparsity.triplet(3, 3, [2, 0, 1], [1, 0, 0], false);
        julia> mapping
        3-element Vector{Int64}:
         1
         2
         0
        julia> sp, mapping = ca.Sparsity.triplet(3, 3, [2, 0, 1], [1, 0, 0], true);
        julia> mapping
        3-element Vector{Int64}:
         2
         0
         1
        julia> try
                   ca.Sparsity.triplet(2, 2, [0, 2], [0, 0])
               catch e
                   println("out of bounds")
               end
        out of bounds
        \enddoctest
        \endlang

        \sa get_triplet, rowcol, nonzeros, Sparsity(nrow, ncol, colind, row)

        \identifier{bm} */
    static Sparsity triplet(casadi_int nrow, casadi_int ncol,
                            const std::vector<casadi_int>& row, const std::vector<casadi_int>& col,
                            std::vector<casadi_int>& SWIG_OUTPUT(mapping), bool invert_mapping);
    static Sparsity triplet(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& row,
                            const std::vector<casadi_int>& col);
    ///@}

    /** \brief Create a sparsity pattern from the flat (column-major) index of every nonzero

        Inverse of find(): Sparsity.nonzeros(sp.size1(), sp.size2(), sp.find()) equals sp.
        Duplicates are merged.

        \param nrow Number of rows, must be positive
        \param ncol Number of columns
        \param nz Flat indices, each in [0, nrow*ncol); an index out of range raises an error
        \param ind1 If true, the indices are one-based

        \lang{python}
        \doctest
        >>> DM(Sparsity.nonzeros(3, 2, [0, 4]), 1)
        DM(
        [[1, 00],
         [00, 1],
         [00, 00]])
        >>> sp = Sparsity.banded(3, 1)
        >>> Sparsity.nonzeros(3, 3, sp.find()) == sp
        True
        >>> Sparsity.nonzeros(3, 2, [1, 5], True) == Sparsity.nonzeros(3, 2, [0, 4])
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.nonzeros(3, 2, [1, 5]), 1)
        [[1, 00],
         [00, 1],
         [00, 00]]
        >> sp = Sparsity.banded(3, 1);
        >> Sparsity.nonzeros(3, 3, sp.find()) == sp
        1
        >> Sparsity.nonzeros(3, 2, [0, 4], false) == Sparsity.nonzeros(3, 2, [1, 5])
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.nonzeros(3, 2, [0, 4]), 1)
        DM(
        [[1, 00],
         [00, 1],
         [00, 00]])
        julia> sp = ca.Sparsity.banded(3, 1);
        julia> ca.is_equal(ca.Sparsity.nonzeros(3, 3, ca.find(sp)), sp)
        true
        julia> sp1 = ca.Sparsity.nonzeros(3, 2, [1, 5], true);
        julia> ca.is_equal(sp1, ca.Sparsity.nonzeros(3, 2, [0, 4]))
        true
        \enddoctest
        \endlang

        \sa find, triplet

        \identifier{23g} */
    static Sparsity nonzeros(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& nz,
                              bool ind1=SWIG_IND1);

    /** \brief Create from the compact single-vector representation of a pattern

        The vector is laid out as [nrow, ncol, colind (ncol+1 entries), row (colind[ncol]
        entries)]. This is the format returned by compress() and used by generated C code
        (casadi_int arrays). When colind[ncol] equals nrow*ncol the row entries may be omitted
        (a dense pattern). An inconsistent vector raises an error.

        \param v Compact representation
        \param order_rows If true, unsorted row indices within a column are sorted instead
               of raising an error

        \lang{python}
        \doctest
        >>> sp = Sparsity.compressed([3, 2, 0, 1, 2, 0, 2])
        >>> DM(sp, 1)
        DM(
        [[1, 00],
         [00, 00],
         [00, 1]])
        >>> Sparsity.compressed(sp.compress()) == sp
        True
        >>> Sparsity.compressed([2, 2, 0, 2, 4])
        Sparsity(2x2)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.compressed([3, 2, 0, 1, 2, 0, 2]);
        >> DM(sp, 1)
        [[1, 00],
         [00, 00],
         [00, 1]]
        >> Sparsity.compressed(sp.compress()) == sp
        1
        >> Sparsity.compressed([2, 2, 0, 2, 4])
        2x2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.compressed([3, 2, 0, 1, 2, 0, 2]);
        julia> ca.DM(sp, 1)
        DM(
        [[1, 00],
         [00, 00],
         [00, 1]])
        julia> ca.is_equal(ca.Sparsity.compressed(ca.compress(sp)), sp)
        true
        julia> ca.Sparsity.compressed([2, 2, 0, 2, 4])
        Sparsity(2x2)
        \enddoctest
        \endlang

        \sa compress, Sparsity(nrow, ncol, colind, row)

        \identifier{2nl} */
    ///@{
    static Sparsity compressed(const std::vector<casadi_int>& v, bool order_rows=false);
#ifndef SWIG
    static Sparsity compressed(const casadi_int* v, bool order_rows=false);
#endif // SWIG
    ///@}

    /** \brief Construct the pattern of a permutation matrix P from a permutation vector p

        P is square with one nonzero per row and column, placed such that multiplying a
        column vector v from the left gives the same result as indexing it: P @ v == v[p].
        The inverse of a permutation matrix is its transpose, so \p invert true returns
        the transpose pattern, for which P @ v puts v[i] at position p[i].

        \param p Permutation vector: every integer in [0, len(p)) exactly once;
               anything else raises an error
        \param invert If true, build the inverse (transposed) permutation

        \lang{python}
        \doctest
        >>> P = Sparsity.permutation([2, 0, 1])
        >>> DM(P, 1)
        DM(
        [[00, 00, 1],
         [1, 00, 00],
         [00, 1, 00]])
        >>> DM(P, 1) @ DM([10, 20, 30])
        DM([30, 10, 20])
        >>> DM(Sparsity.permutation([2, 0, 1], True), 1) @ DM([10, 20, 30])
        DM([20, 30, 10])
        >>> Sparsity.permutation([2, 0, 1], True) == P.T
        True
        >>> try: Sparsity.permutation([0, 0, 1])
        ... except Exception as e: print("not a permutation")
        not a permutation
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> P = Sparsity.permutation([2, 0, 1]);
        >> DM(P, 1)
        [[00, 00, 1],
         [1, 00, 00],
         [00, 1, 00]]
        >> DM(P, 1) * DM([10; 20; 30])
        [30, 10, 20]
        >> DM(Sparsity.permutation([2, 0, 1], true), 1) * DM([10; 20; 30])
        [20, 30, 10]
        >> Sparsity.permutation([2, 0, 1], true) == P.T
        1
        >> try, Sparsity.permutation([0, 0, 1]); catch e, disp('not a permutation'); end
        not a permutation
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> P = ca.Sparsity.permutation([2, 0, 1]);
        julia> ca.DM(P, 1)
        DM(
        [[00, 00, 1],
         [1, 00, 00],
         [00, 1, 00]])
        julia> ca.DM(P, 1) * ca.DM([10, 20, 30])
        DM([30, 10, 20])
        julia> ca.DM(ca.Sparsity.permutation([2, 0, 1], true), 1) * ca.DM([10, 20, 30])
        DM([20, 30, 10])
        julia> ca.is_equal(ca.Sparsity.permutation([2, 0, 1], true), ca.T(P))
        true
        julia> try
                   ca.Sparsity.permutation([0, 0, 1])
               catch e
                   println("not a permutation")
               end
        not a permutation
        \enddoctest
        \endlang

        \sa permutation_vector, is_permutation, pmult

        \identifier{244} */
    static Sparsity permutation(const std::vector<casadi_int>& p, bool invert=false);

    /** \brief Recover the permutation vector p from the pattern of a permutation matrix

        Inverse of permutation(): Sparsity.permutation(p).permutation_vector() returns p.
        Raises an error if is_permutation() is false.

        \param invert If true, return the inverse permutation instead

        \lang{python}
        \doctest
        >>> Sparsity.permutation([2, 0, 1]).permutation_vector()
        [2, 0, 1]
        >>> Sparsity.permutation([2, 0, 1]).permutation_vector(True)
        [1, 2, 0]
        >>> try: Sparsity.dense(2, 2).permutation_vector()
        ... except Exception as e: print("not a permutation matrix")
        not a permutation matrix
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.permutation([2, 0, 1]).permutation_vector()
        2 0 1
        >> Sparsity.permutation([2, 0, 1]).permutation_vector(true)
        1 2 0
        >> try, Sparsity.dense(2, 2).permutation_vector(); catch e, disp('not a permutation'); end
        not a permutation
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.permutation_vector(ca.Sparsity.permutation([2, 0, 1]))
        3-element Vector{Int64}:
         2
         0
         1
        julia> ca.permutation_vector(ca.Sparsity.permutation([2, 0, 1]), true)
        3-element Vector{Int64}:
         1
         2
         0
        julia> try
                   ca.permutation_vector(ca.Sparsity.dense(2, 2))
               catch e
                   println("not a permutation matrix")
               end
        not a permutation matrix
        \enddoctest
        \endlang

        \sa permutation, is_permutation

        \identifier{245} */
    const std::vector<casadi_int> permutation_vector(bool invert=false) const;

    /** \brief Pattern of the diagonal of a matrix, or of the diagonal matrix built from a vector

        If the pattern is not a vector, the result is the min(nrow, ncol)-by-1 pattern of the
        main diagonal: it has a nonzero for every diagonal entry (i, i) present in the pattern.
        If the pattern is a row or column vector of length n, the result is the n-by-n pattern
        with the vector entries on the diagonal.

        \return The pattern and a mapping: mapping[k] is the nonzero index in the original
                pattern of nonzero k of the result

        \lang{python}
        \doctest
        >>> sp, mapping = Sparsity.banded(3, 1).get_diag()
        >>> sp, mapping
        (Sparsity(3x1), [0, 3, 6])
        >>> sp, mapping = Sparsity.unit(3, 1).get_diag()
        >>> DM(sp, 1), mapping
        (DM(
        [[00, 00, 00],
         [00, 1, 00],
         [00, 00, 00]]), [0])
        >>> sp, mapping = Sparsity(2, 3).get_diag()
        >>> sp, mapping
        (Sparsity(2x1,0nz), [])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [sp, mapping] = Sparsity.banded(3, 1).get_diag(); disp(sp), disp(mapping)
        3x1
        0 3 6
        >> [sp, mapping] = Sparsity.unit(3, 1).get_diag(); disp(DM(sp, 1)), disp(mapping)
        [[00, 00, 00],
         [00, 1, 00],
         [00, 00, 00]]
        0
        >> [sp, mapping] = Sparsity(2, 3).get_diag(); disp(sp), disp(isempty(mapping))
        2x1,0nz
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.get_diag(ca.Sparsity.banded(3, 1))
        (3x1, [0, 3, 6])
        julia> ca.get_diag(ca.Sparsity.unit(3, 1))
        (3x3,1nz, [0])
        julia> ca.get_diag(ca.Sparsity(2, 3))
        (2x1,0nz, Int64[])
        \enddoctest
        \endlang

        \sa diag, nnz_diag

        \identifier{2nm} */
    Sparsity get_diag(std::vector<casadi_int>& SWIG_OUTPUT(mapping)) const;

    /** \brief Compact single-vector representation of the pattern

        Returns [nrow, ncol, colind..., row...], the layout accepted by compressed() and used
        for sparsity patterns in generated C code. With \p canonical false, a dense pattern
        is abbreviated to [nrow, ncol, 1] (an abbreviation that compressed() given a
        vector does not accept).

        \lang{python}
        \doctest
        >>> Sparsity.diag(2).compress()
        [2, 2, 0, 1, 2, 0, 1]
        >>> Sparsity.dense(2, 2).compress()
        [2, 2, 0, 2, 4, 0, 1, 0, 1]
        >>> Sparsity.dense(2, 2).compress(False)
        [2, 2, 1]
        >>> Sparsity(0, 2).compress()
        [0, 2, 0, 0, 0]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.diag(2).compress()
        2 2 0 1 2 0 1
        >> Sparsity.dense(2, 2).compress()
        2 2 0 2 4 0 1 0 1
        >> Sparsity.dense(2, 2).compress(false)
        2 2 1
        >> Sparsity(0, 2).compress()
        0 2 0 0 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> println(ca.compress(ca.Sparsity.diag(2)))
        [2, 2, 0, 1, 2, 0, 1]
        julia> println(ca.compress(ca.Sparsity.dense(2, 2)))
        [2, 2, 0, 2, 4, 0, 1, 0, 1]
        julia> println(ca.compress(ca.Sparsity.dense(2, 2), false))
        [2, 2, 1]
        julia> println(ca.compress(ca.Sparsity(0, 2)))
        [0, 2, 0, 0, 0]
        \enddoctest
        \endlang

        \sa compressed, get_ccs

        \identifier{2nn} */
    std::vector<casadi_int> compress(bool canonical=true) const;

#ifndef SWIG
    /// Access a member function or object
    const SparsityInternal* operator->() const;

    /// Reference to internal structure
    const SparsityInternal& operator*() const;
#endif // SWIG
    /** \brief Check whether two sparsity patterns are identical

        True when the dimensions and the nonzero positions coincide. This is the comparison
        behind the == operator. Shared (cached) patterns are detected in constant time;
        otherwise colind and row are compared. The second overload compares against a
        pattern given in compressed column storage vectors.

        \lang{python}
        \doctest
        >>> a = Sparsity.diag(2)
        >>> a.is_equal(Sparsity.diag(2)), a == Sparsity.diag(2), a != Sparsity.dense(2, 2)
        (True, True, True)
        >>> a.is_equal(2, 2, [0, 1, 2], [0, 1])
        True
        >>> Sparsity.dense(2, 2) == Sparsity.dense(2, 3)
        False
        >>> Sparsity(2, 2) == Sparsity.dense(2, 2)
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> a = Sparsity.diag(2);
        >> a.is_equal(Sparsity.diag(2)), a == Sparsity.diag(2)
        1
        1
        >> a ~= Sparsity.dense(2, 2)
        1
        >> a.is_equal(2, 2, [0, 1, 2], [0, 1])
        1
        >> Sparsity.dense(2, 2) == Sparsity.dense(2, 3)
        0
        >> Sparsity(2, 2) == Sparsity.dense(2, 2)
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> a = ca.Sparsity.diag(2);
        julia> ca.is_equal(a, ca.Sparsity.diag(2)), !ca.is_equal(a, ca.Sparsity.dense(2, 2))
        (true, true)
        julia> ca.is_equal(a, 2, 2, [0, 1, 2], [0, 1])
        true
        julia> ca.is_equal(ca.Sparsity.dense(2, 2), ca.Sparsity.dense(2, 3))
        false
        julia> ca.is_equal(ca.Sparsity(2, 2), ca.Sparsity.dense(2, 2))
        false
        \enddoctest
        \endlang

        \sa is_subset, is_transpose, is_reshape, hash

        \identifier{2no} */
    /// @{
    bool is_equal(const Sparsity& y) const;
    bool is_equal(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& colind,
                 const std::vector<casadi_int>& row) const;
#ifndef SWIG
    bool is_equal(casadi_int nrow, casadi_int ncol,
                  const casadi_int* colind, const casadi_int* row) const;
#endif // SWIG

    bool operator==(const Sparsity& y) const { return is_equal(y);}
    /// @}

    /// Check whether two sparsity patterns differ (negation of is_equal)
    bool operator!=(const Sparsity& y) const {return !is_equal(y);}

    /** \brief Check whether the pattern is \p n copies of \p y placed side by side

        True when this pattern equals horzcat(y, ..., y) with \p n repetitions, i.e. the
        pattern of repmat(y, 1, n). The patterns are compared structurally, so a dense
        pattern of the right shape does not count as a stack of a sparse one.

        \lang{python}
        \doctest
        >>> sp = horzcat(Sparsity.diag(2), Sparsity.diag(2), Sparsity.diag(2))
        >>> sp.is_stacked(Sparsity.diag(2), 3), sp.is_stacked(Sparsity.diag(2), 2)
        (True, False)
        >>> Sparsity.dense(2, 6).is_stacked(Sparsity.diag(2), 3)
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = [Sparsity.diag(2), Sparsity.diag(2), Sparsity.diag(2)];
        >> sp.is_stacked(Sparsity.diag(2), 3), sp.is_stacked(Sparsity.diag(2), 2)
        1
        0
        >> Sparsity.dense(2, 6).is_stacked(Sparsity.diag(2), 3)
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.repmat(ca.Sparsity.diag(2), 1, 3)
        Sparsity(2x6,6nz)
        julia> ca.is_stacked(sp, ca.Sparsity.diag(2), 3), ca.is_stacked(sp, ca.Sparsity.diag(2), 2)
        (true, false)
        julia> ca.is_stacked(ca.Sparsity.dense(2, 6), ca.Sparsity.diag(2), 3)
        false
        \enddoctest
        \endlang

        \sa is_equal, horzcat

        \identifier{2np} */
    bool is_stacked(const Sparsity& y, casadi_int n) const;

#ifndef SWIG
    /** \brief Implicit or explicit type conversion to C representation

        In the C runtime, sparsity patterns are represented as a "const casadi_int*".
        This enables using the C runtime functions using a natural syntax.

        \identifier{bn} */
    operator const casadi_int*() const;

    /** \brief Implicit or explicit type conversion to compact representation

        \identifier{bo} */
    operator const std::vector<casadi_int>&() const;

    /** \brief Implicit or explicit type conversion to C representation

        \identifier{bp} */
    operator SparsityStruct() const;
#endif // SWIG

    /// \name Size and element counting
    /// @{

    /** \brief Number of rows

        \lang{python}
        \doctest
        >>> sp = Sparsity.diag(2, 3)
        >>> sp.size1(), sp.rows(), sp.size2(), sp.columns()
        (2, 2, 3, 3)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.diag(2, 3);
        >> fprintf('%d %d %d %d\n', sp.size1(), sp.rows(), sp.size2(), sp.columns())
        2 2 3 3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.diag(2, 3);
        julia> ca.size1(sp), ca.rows(sp), ca.size2(sp), ca.columns(sp)
        (2, 2, 3, 3)
        \enddoctest
        \endlang

        \sa size2, size, numel

        \identifier{2nq} */
    casadi_int size1() const;

    /// Number of rows, Octave-style alias of size1
    casadi_int rows() const {return size1();}

    /// Number of columns
    casadi_int size2() const;

    /// Number of columns, Octave-style alias of size2
    casadi_int columns() const {return size2();}

    /** \brief Total number of entries, structural zeros included, i.e. size1()*size2()

        The product may overflow for very large shapes; nnz() never exceeds it.

        \lang{python}
        \doctest
        >>> sp = Sparsity.diag(2, 3)
        >>> sp.numel(), sp.nnz(), Sparsity(0, 3).numel()
        (6, 2, 0)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.diag(2, 3);
        >> fprintf('%d %d %d\n', sp.numel(), sp.nnz(), Sparsity(0, 3).numel())
        6 2 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.diag(2, 3);
        julia> ca.numel(sp), ca.nnz(sp), ca.numel(ca.Sparsity(0, 3))
        (6, 2, 0)
        \enddoctest
        \endlang

        \sa nnz, size

        \identifier{bq} */
    casadi_int numel() const;

    /** \brief Percentage of nonzeros: 100*nnz()/numel()

        Computed in floating point, so it does not overflow for large shapes.
        An empty pattern (numel() == 0) is dense and gives 100.

        \lang{python}
        \doctest
        >>> Sparsity.diag(4).density(), Sparsity.dense(2, 2).density()
        (25.0, 100.0)
        >>> Sparsity(0, 3).density()
        100.0
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%g %g\n', Sparsity.diag(4).density(), Sparsity.dense(2, 2).density())
        25 100
        >> Sparsity(0, 3).density()
        100
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.density(ca.Sparsity.diag(4)), ca.density(ca.Sparsity.dense(2, 2))
        (25.0, 100.0)
        julia> ca.density(ca.Sparsity(0, 3))
        100.0
        \enddoctest
        \endlang

        \sa nnz, numel, is_dense

        \identifier{br} */
    double density() const;

    /** \brief Check whether the pattern has a zero dimension

        \param both If false (default), true when the number of rows or the number of
               columns is zero; if true, only when both are zero

        A pattern without nonzeros but with positive dimensions is not empty.

        \lang{python}
        \doctest
        >>> Sparsity(0, 3).is_empty(), Sparsity(0, 3).is_empty(True)
        (True, False)
        >>> Sparsity(0, 0).is_empty(True)
        True
        >>> Sparsity(3, 3).is_empty(), Sparsity(3, 3).nnz()
        (False, 0)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity(0, 3).is_empty(), Sparsity(0, 3).is_empty(true)
        1
        0
        >> Sparsity(0, 0).is_empty(true)
        1
        >> Sparsity(3, 3).is_empty(), Sparsity(3, 3).nnz()
        0
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_empty(ca.Sparsity(0, 3)), ca.is_empty(ca.Sparsity(0, 3), true)
        (true, false)
        julia> ca.is_empty(ca.Sparsity(0, 0), true)
        true
        julia> ca.is_empty(ca.Sparsity(3, 3)), ca.nnz(ca.Sparsity(3, 3))
        (false, 0)
        \enddoctest
        \endlang

        \sa numel, nnz

        \identifier{bs} */
    bool is_empty(bool both=false) const;

    /** \brief Number of structural nonzeros

        Equals colind[size2()] and the length of the row vector. Zero for an empty pattern
        and for a pattern created with Sparsity(nrow, ncol).

        \sa numel, density, nnz_upper, nnz_lower, nnz_diag

        \identifier{bt} */
    casadi_int nnz() const;

    /** \brief Number of nonzeros in the upper triangular part, i.e. entries (i, j) with i <= j

        \param strictly If true, exclude the diagonal (count i < j only)

        Defined for any shape, not only square patterns.

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(3, 3)
        >>> sp.nnz_upper(), sp.nnz_upper(True), sp.nnz_lower(), sp.nnz_lower(True)
        (6, 3, 6, 3)
        >>> sp.nnz_diag()
        3
        >>> sp = Sparsity.dense(2, 4)
        >>> sp.nnz_upper(), sp.nnz_lower(), sp.nnz_diag()
        (7, 3, 2)
        >>> Sparsity.lower(3).nnz_upper(), Sparsity.lower(3).nnz_upper(True)
        (3, 0)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(3, 3);
        >> [sp.nnz_upper(), sp.nnz_upper(true), sp.nnz_lower(), sp.nnz_lower(true)]
        6 3 6 3
        >> sp.nnz_diag()
        3
        >> sp = Sparsity.dense(2, 4);
        >> [sp.nnz_upper(), sp.nnz_lower(), sp.nnz_diag()]
        7 3 2
        >> [Sparsity.lower(3).nnz_upper(), Sparsity.lower(3).nnz_upper(true)]
        3 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(3, 3);
        julia> ca.nnz_upper(sp), ca.nnz_upper(sp, true), ca.nnz_lower(sp), ca.nnz_lower(sp, true)
        (6, 3, 6, 3)
        julia> ca.nnz_diag(sp)
        3
        julia> sp = ca.Sparsity.dense(2, 4);
        julia> ca.nnz_upper(sp), ca.nnz_lower(sp), ca.nnz_diag(sp)
        (7, 3, 2)
        julia> ca.nnz_upper(ca.Sparsity.lower(3)), ca.nnz_upper(ca.Sparsity.lower(3), true)
        (3, 0)
        \enddoctest
        \endlang

        \sa nnz_lower, nnz_diag, get_upper, is_triu

        \identifier{bu} */
    casadi_int nnz_upper(bool strictly=false) const;

    /** \brief Number of nonzeros in the lower triangular part, i.e. entries (i, j) with i >= j

        \param strictly If true, exclude the diagonal (count i > j only)

        \sa nnz_upper, nnz_diag, get_lower, is_tril

        \identifier{bv} */
    casadi_int nnz_lower(bool strictly=false) const;

    /** \brief Number of nonzeros on the main diagonal, i.e. entries (i, i)

        \sa nnz_upper, nnz_lower, get_diag

        \identifier{bw} */
    casadi_int nnz_diag() const;

    /** \brief Upper half-bandwidth: the largest j - i over all nonzeros (i, j) with j > i

        Zero for a diagonal, lower triangular or all-zero pattern.

        \lang{python}
        \doctest
        >>> Sparsity.banded(5, 2).bw_upper(), Sparsity.band(5, -2).bw_upper()
        (2, 2)
        >>> Sparsity.band(5, 2).bw_upper(), Sparsity.diag(3).bw_upper()
        (0, 0)
        >>> Sparsity(3, 3).bw_upper(), Sparsity.dense(2, 5).bw_upper()
        (0, 4)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.banded(5, 2).bw_upper(), Sparsity.band(5, -2).bw_upper()]
        2 2
        >> [Sparsity.band(5, 2).bw_upper(), Sparsity.diag(3).bw_upper()]
        0 0
        >> [Sparsity(3, 3).bw_upper(), Sparsity.dense(2, 5).bw_upper()]
        0 4
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.bw_upper(ca.Sparsity.banded(5, 2)), ca.bw_upper(ca.Sparsity.band(5, -2))
        (2, 2)
        julia> ca.bw_upper(ca.Sparsity.band(5, 2)), ca.bw_upper(ca.Sparsity.diag(3))
        (0, 0)
        julia> ca.bw_upper(ca.Sparsity(3, 3)), ca.bw_upper(ca.Sparsity.dense(2, 5))
        (0, 4)
        \enddoctest
        \endlang

        \sa bw_lower, banded

        \identifier{bx} */
    casadi_int bw_upper() const;

    /** \brief Lower half-bandwidth: the largest i - j over all nonzeros (i, j) with i > j

        Zero for a diagonal, upper triangular or all-zero pattern.

        \lang{python}
        \doctest
        >>> Sparsity.banded(5, 2).bw_lower(), Sparsity.band(5, 2).bw_lower()
        (2, 2)
        >>> Sparsity.band(5, -2).bw_lower()
        0
        >>> Sparsity.dense(5, 2).bw_lower()
        4
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.banded(5, 2).bw_lower(), Sparsity.band(5, 2).bw_lower()]
        2 2
        >> Sparsity.band(5, -2).bw_lower()
        0
        >> Sparsity.dense(5, 2).bw_lower()
        4
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.bw_lower(ca.Sparsity.banded(5, 2)), ca.bw_lower(ca.Sparsity.band(5, 2))
        (2, 2)
        julia> ca.bw_lower(ca.Sparsity.band(5, -2))
        0
        julia> ca.bw_lower(ca.Sparsity.dense(5, 2))
        4
        \enddoctest
        \endlang

        \sa bw_upper, banded

        \identifier{by} */
    casadi_int bw_lower() const;

    /** \brief Shape as a pair (number of rows, number of columns)

        \lang{python}
        \doctest
        >>> sp = Sparsity.diag(2, 3)
        >>> sp.size(), sp.size(1), sp.size(2), sp.shape
        ((2, 3), 2, 3, (2, 3))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.diag(2, 3);
        >> size(sp), sp.size(1), sp.size(2)
        2 3
        2
        3
        >> [m, n] = size(sp); fprintf('%d %d\n', m, n)
        2 3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.diag(2, 3);
        julia> ca.size(sp), ca.size(sp, 1), ca.size(sp, 2)
        ((2, 3), 2, 3)
        \enddoctest
        \endlang

        \sa size1, size2, numel

        \identifier{bz} */
    std::pair<casadi_int, casadi_int> size() const;

    /** \brief Size along one dimension

        \param axis 1 for the number of rows, 2 for the number of columns; any other value
               raises an error

        \identifier{c0} */
    casadi_int size(casadi_int axis) const;
    /// @}

    /** \brief The pattern as a dictionary with entries nrow, ncol, colind and row

        An empty dictionary for a null object.

        \lang{python}
        \doctest
        >>> Sparsity.diag(2, 3).info()
        {'colind': [0, 1, 2, 2], 'ncol': 3, 'nrow': 2, 'row': [0, 1]}
        >>> Sparsity().info()
        {}
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> s = Sparsity.diag(2, 3).info(); disp(s)
        colind: [0 1 2 2]
        ncol: 3
        nrow: 2
        row: [0 1]
        >> isempty(fieldnames(Sparsity().info()))
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.info(ca.Sparsity.diag(2, 3))
        Dict{String, Any} with 4 entries:
          "nrow" => 2
          "colind" => [0, 1, 2, 2]
          "ncol" => 3
          "row" => [0, 1]
        julia> ca.info(ca.Sparsity())
        Dict{String, Any}()
        \enddoctest
        \endlang

        \sa get_ccs, compress

        \identifier{2nr} */
    Dict info() const;

    /** \brief Write the pattern to a file

        The only supported format is Matrix Market ("mtx", coordinate pattern general),
        with one-based (row, col) pairs. The format is taken from the file extension unless
        \p format_hint is given; an unknown format raises an error. The file is overwritten.

        \param filename Path of the file to write
        \param format_hint Format name, e.g. "mtx", when the extension does not tell

        \lang{python}
        \doctest
        >>> import os, tempfile
        >>> d = tempfile.mkdtemp()
        >>> f = os.path.join(d, "pattern.mtx")
        >>> sp = Sparsity.triplet(3, 2, [0, 2, 1], [0, 0, 1])
        >>> sp.to_file(f)
        >>> print(open(f).read())
        %%MatrixMarket matrix coordinate pattern general
        3 2 3
        1 1
        3 1
        2 2
        >>> print(Sparsity.from_file(f) == sp)
        True
        >>> os.remove(f); os.rmdir(d)
        >>> try: sp.to_file("pattern.txt")
        ... except Exception as e: print("unknown extension")
        unknown extension
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = [tempname '.mtx'];
        >> sp = Sparsity.triplet(3, 2, [0, 2, 1], [0, 0, 1]);
        >> sp.to_file(f);
        >> fileread(f)
        %%MatrixMarket matrix coordinate pattern general
        3 2 3
        1 1
        3 1
        2 2
        >> Sparsity.from_file(f) == sp
        1
        >> delete(f);
        >> try, sp.to_file('pattern.txt'); catch e, disp('unknown extension'); end
        unknown extension
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> d = mktempdir();
        julia> f = joinpath(d, "pattern.mtx");
        julia> sp = ca.Sparsity.triplet(3, 2, [0, 2, 1], [0, 0, 1]);
        julia> ca.to_file(sp, f)
        julia> print(read(f, String))
        %%MatrixMarket matrix coordinate pattern general
        3 2 3
        1 1
        3 1
        2 2
        julia> ca.is_equal(ca.Sparsity.from_file(f), sp)
        true
        julia> rm(d, recursive=true)
        julia> try
                   ca.to_file(sp, "pattern.txt")
               catch e
                   println("unknown extension")
               end
        unknown extension
        \enddoctest
        \endlang

        \sa from_file, DM::to_file

        \identifier{2ns} */
    void to_file(const std::string& filename, const std::string& format_hint="") const;

    /** \brief Read a pattern from a file written by to_file (Matrix Market pattern format)

        \param filename Path of the file to read; a missing file raises an error
        \param format_hint Format name, e.g. "mtx", when the extension does not tell

        \sa to_file

        \identifier{2nt} */
    static Sparsity from_file(const std::string& filename, const std::string& format_hint="");

#ifndef SWIG
    /** \brief Serialize

        \identifier{c1} */
    void serialize(std::ostream &stream) const;
#endif

    /** \brief Serialize the pattern to a string

        The string can be stored and turned back into an equal pattern with deserialize().
        A null object serializes and deserializes as a null object. Pickling in Python uses
        this mechanism.

        \lang{python}
        \doctest
        >>> s = Sparsity.diag(2, 3).serialize()
        >>> type(s).__name__, Sparsity.deserialize(s) == Sparsity.diag(2, 3)
        ('str', True)
        >>> Sparsity.deserialize(Sparsity().serialize()).is_null()
        True
        >>> import pickle
        >>> pickle.loads(pickle.dumps(Sparsity.diag(2, 3)))
        Sparsity(2x3,2nz)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> s = Sparsity.diag(2, 3).serialize();
        >> class(s), Sparsity.deserialize(s) == Sparsity.diag(2, 3)
        char
        1
        >> Sparsity.deserialize(Sparsity().serialize()).is_null()
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> s = ca.serialize(ca.Sparsity.diag(2, 3));
        julia> typeof(s), ca.is_equal(ca.Sparsity.deserialize(s), ca.Sparsity.diag(2, 3))
        (String, true)
        julia> ca.is_null(ca.Sparsity.deserialize(ca.serialize(ca.Sparsity())))
        true
        \enddoctest
        \endlang

        \sa deserialize, compress

        \identifier{c2} */
    std::string serialize() const;

    /** \brief Build a Sparsity from a serialized stream

        \identifier{c3} */
    static Sparsity deserialize(std::istream& stream);

    /** \brief Build a Sparsity from a string produced by serialize()

        \sa serialize

        \identifier{c4} */
    static Sparsity deserialize(const std::string& s);

    /** \brief Serialize an object

        \identifier{c5} */
    void serialize(SerializingStream& s) const;

    /** \brief Deserialize

        \identifier{c6} */
    static Sparsity deserialize(DeserializingStream& s);

#ifndef SWIG
    /** \brief Get a reference to row-vector,

     * containing rows for all non-zero elements (see class description)

        \identifier{c7} */
    const casadi_int* row() const;

    /** \brief Get a reference to the colindex of all column element (see class description)

        \identifier{c8} */
    const casadi_int* colind() const;
#endif

    /** \brief Row index of every nonzero, in nonzero (column-major) order

        Length nnz(). Together with get_colind() this is the compressed column storage of
        the pattern; together with get_col() it is the triplet form. In the target languages
        this method is called row(), the same name as the single-element accessor.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> DM(sp, 1)
        DM(
        [[1, 00, 00, 1],
         [00, 00, 1, 00],
         [1, 00, 00, 00]])
        >>> sp.row(), sp.colind(), sp.get_col()
        ([0, 2, 1, 0], [0, 2, 2, 3, 4], [0, 0, 2, 3])
        >>> Sparsity.dense(2, 2).row(), Sparsity(2, 2).colind()
        ([0, 1, 0, 1], [0, 0, 0])
        >>> Sparsity(0, 2).colind(), Sparsity(2, 0).colind()
        ([0, 0, 0], [0])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> DM(sp, 1)
        [[1, 00, 00, 1],
         [00, 00, 1, 00],
         [1, 00, 00, 00]]
        >> sp.row(), sp.colind(), sp.get_col()
        0 2 1 0
        0 2 2 3 4
        0 0 2 3
        >> Sparsity.dense(2, 2).row(), Sparsity(2, 2).colind()
        0 1 0 1
        0 0 0
        >> Sparsity(0, 2).colind(), Sparsity(2, 0).colind()
        0 0 0
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.row(sp), ca.colind(sp), ca.get_col(sp)
        ([0, 2, 1, 0], [0, 2, 2, 3, 4], [0, 0, 2, 3])
        julia> ca.row(ca.Sparsity.dense(2, 2)), ca.colind(ca.Sparsity(2, 2))
        ([0, 1, 0, 1], [0, 0, 0])
        julia> ca.colind(ca.Sparsity(0, 2)), ca.colind(ca.Sparsity(2, 0))
        ([0, 0, 0], [0])
        \enddoctest
        \endlang

        \sa get_colind, get_col, get_ccs, get_triplet

        \identifier{c9} */
    std::vector<casadi_int> get_row() const;

    /** \brief Column offsets of the compressed column storage

        Length size2()+1: entry j is the index of the first nonzero of column j and the last
        entry equals nnz(). In the target languages this method is called colind(), the
        same name as the single-column accessor.

        \sa get_row, get_ccs

        \identifier{ca} */
    std::vector<casadi_int> get_colind() const;

    /** \brief Index of the first nonzero of column \p cc

        \param cc Column, in [0, size2()]; colind(size2()) equals nnz(). An index out of
               range raises an error.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sp.row(1), sp.colind(1), sp.colind(4)
        (2, 2, 4)
        >>> try: sp.row(4)
        ... except Exception as e: print("out of range")
        out of range
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> [sp.row(1), sp.colind(1), sp.colind(4)]
        2 2 4
        >> try, sp.row(4); catch e, disp('out of range'); end
        out of range
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.row(sp, 1), ca.colind(sp, 1), ca.colind(sp, 4)
        (2, 2, 4)
        julia> try
                   ca.row(sp, 4)
               catch e
                   println("out of range")
               end
        out of range
        \enddoctest
        \endlang

        \sa row, get_colind

        \identifier{cb} */
    casadi_int colind(casadi_int cc) const;

    /** \brief Row index of nonzero \p el

        \param el Nonzero index, in [0, nnz()); out of range raises an error

        \sa colind, get_row

        \identifier{cc} */
    casadi_int row(casadi_int el) const;

    /** \brief Column index of every nonzero, in nonzero (column-major) order

        Length nnz(). Together with get_row() this gives the pattern in sparse triplet
        format.

        \sa get_row, get_triplet

        \identifier{cd} */
    std::vector<casadi_int> get_col() const;

    /** \brief Change the shape in place, keeping the nonzeros that still fit

        Nonzeros with row >= \p nrow or column >= \p ncol are dropped; growing adds
        structural zeros. Unlike reshape, the (row, col) position of the kept nonzeros is
        unchanged. The object is rebound to the new pattern; other references to the old
        pattern are unaffected.

        \param nrow New number of rows, nonnegative
        \param ncol New number of columns, nonnegative

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(3, 3); sp.resize(2, 2); print(sp)
        2x2
        >>> sp.resize(4, 4); print(DM(sp, 1))
        [[1, 1, 00, 00], 
         [1, 1, 00, 00], 
         [00, 00, 00, 00], 
         [00, 00, 00, 00]]
        >>> sp.resize(0, 2); print(sp)
        0x2
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(3, 3); sp.resize(2, 2); disp(sp)
        2x2
        >> sp.resize(4, 4); disp(DM(sp, 1))
        [[1, 1, 00, 00],
         [1, 1, 00, 00],
         [00, 00, 00, 00],
         [00, 00, 00, 00]]
        >> sp.resize(0, 2); disp(sp)
        0x2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(3, 3); ca.resize(sp, 2, 2); sp
        Sparsity(2x2)
        julia> ca.resize(sp, 4, 4); ca.DM(sp, 1)
        DM(
        [[1, 1, 00, 00],
         [1, 1, 00, 00],
         [00, 00, 00, 00],
         [00, 00, 00, 00]])
        julia> ca.resize(sp, 0, 2); sp
        Sparsity(0x2)
        \enddoctest
        \endlang

        \sa enlarge, reshape

        \identifier{2nu} */
    void resize(casadi_int nrow, casadi_int ncol);

    /** \brief Nonzero index of entry (rr, cc), inserting the entry into the pattern if needed

        Returns the position of the entry among the nonzeros. If the entry is a structural
        zero, the object is rebound to a pattern that includes it (nonzero indices after it
        shift by one). Negative indices count from the end. An index out of bounds raises
        an error.

        \lang{python}
        \doctest
        >>> sp = Sparsity(3, 3)
        >>> sp.add_nz(1, 1), sp.add_nz(0, 0), sp.add_nz(1, 1), sp
        (0, 0, 1, Sparsity(3x3,2nz))
        >>> DM(sp, 1)
        DM(
        [[1, 00, 00],
         [00, 1, 00],
         [00, 00, 00]])
        >>> try: sp.add_nz(3, 0)
        ... except Exception as e: print("out of bounds")
        out of bounds
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity(3, 3);
        >> [sp.add_nz(1, 1), sp.add_nz(0, 0), sp.add_nz(1, 1)], sp
        0 0 1
        3x3,2nz
        >> DM(sp, 1)
        [[1, 00, 00],
         [00, 1, 00],
         [00, 00, 00]]
        >> try, sp.add_nz(3, 0); catch e, disp('out of bounds'); end
        out of bounds
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity(3, 3);
        julia> ca.add_nz(sp, 1, 1), ca.add_nz(sp, 0, 0), ca.add_nz(sp, 1, 1)
        (0, 0, 1)
        julia> ca.DM(sp, 1)
        DM(
        [[1, 00, 00],
         [00, 1, 00],
         [00, 00, 00]])
        julia> try
                   ca.add_nz(sp, 3, 0)
               catch e
                   println("out of bounds")
               end
        out of bounds
        \enddoctest
        \endlang

        \sa get_nz, has_nz

        \identifier{ce} */
    casadi_int add_nz(casadi_int rr, casadi_int cc);

    /** \brief Nonzero index of entry (rr, cc), or -1 if it is a structural zero

        Negative indices count from the end. An index out of bounds raises an error.
        For a dense pattern the result is rr + cc*size1(). The returned index addresses
        the nonzeros() vector of a matrix with this pattern.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sp.get_nz(0, 0), sp.get_nz(2, 0), sp.get_nz(1, 0), sp.get_nz(-1, 0)
        (0, 1, -1, 1)
        >>> try: sp.get_nz(3, 0)
        ... except Exception as e: print("out of bounds")
        out of bounds
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> [sp.get_nz(0, 0), sp.get_nz(2, 0), sp.get_nz(1, 0), sp.get_nz(-1, 0)]
        0 1 -1 1
        >> try, sp.get_nz(3, 0); catch e, disp('out of bounds'); end
        out of bounds
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.get_nz(sp, 0, 0), ca.get_nz(sp, 2, 0), ca.get_nz(sp, 1, 0), ca.get_nz(sp, -1, 0)
        (0, 1, -1, 1)
        julia> try
                   ca.get_nz(sp, 3, 0)
               catch e
                   println("out of bounds")
               end
        out of bounds
        \enddoctest
        \endlang

        \sa has_nz, add_nz, find

        \identifier{cf} */
    casadi_int get_nz(casadi_int rr, casadi_int cc) const;

    /** \brief Check whether entry (rr, cc) is a structural nonzero

        Same index rules as get_nz.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sp.has_nz(2, 0), sp.has_nz(1, 0)
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> [sp.has_nz(2, 0), sp.has_nz(1, 0)]
        1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.has_nz(sp, 2, 0), ca.has_nz(sp, 1, 0)
        (true, false)
        \enddoctest
        \endlang

        \sa get_nz

        \identifier{2nv} */
    bool has_nz(casadi_int rr, casadi_int cc) const;

    /** \brief Nonzero indices of the block of entries with rows \p rr and columns \p cc

        Returns one index per (row, column) combination, in column-major order over the
        block (all rows of the first column first); -1 marks a structural zero. Indices must
        be in range (no negative indices); otherwise an error is raised.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sp.get_nz([0, 2], [0, 3])
        [0, 1, 3, -1]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> sp.get_nz([0, 2], [0, 3])
        0 1 3 -1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> println(ca.get_nz(sp, [0, 2], [0, 3]))
        [0, 1, 3, -1]
        \enddoctest
        \endlang

        \sa get_nz(rr, cc), sub

        \identifier{cg} */
    std::vector<casadi_int> get_nz(const std::vector<casadi_int>& rr,
                                    const std::vector<casadi_int>& cc) const;

    /** \brief Map flat (column-major) element indices to nonzero indices

        Each entry of \p indices is replaced by the nonzero index of that element, or by -1
        if the element is a structural zero. Entries that are already negative are left
        alone. The indices need not be sorted. An index of numel() or more raises an error.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sp.get_nz([0, 2, 4, 9, 1])
        [0, 1, -1, 3, -1]
        >>> Sparsity.dense(2, 2).get_nz([3, 0])
        [3, 0]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> sp.get_nz([0, 2, 4, 9, 1])
        0 1 -1 3 -1
        >> Sparsity.dense(2, 2).get_nz([3, 0])
        3 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.get_nz(sp, [0, 2, 4, 9, 1])
        5-element Vector{Int64}:
          0
          1
         -1
          3
         -1
        julia> ca.get_nz(ca.Sparsity.dense(2, 2), [3, 0])
        2-element Vector{Int64}:
         3
         0
        \enddoctest
        \endlang

        \sa find, nonzeros

        \identifier{ch} */
    void get_nz(std::vector<casadi_int>& SWIG_INOUT(indices)) const;

    /** \brief Nonzero indices of the entries in the lower triangular part (row >= col)

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sp.get_lower(), sp.get_upper()
        ([0, 1], [0, 2, 3])
        >>> Sparsity.dense(2, 2).get_lower(), Sparsity.dense(2, 2).get_upper()
        ([0, 1, 3], [0, 2, 3])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> sp.get_lower(), sp.get_upper()
        0 1
        0 2 3
        >> Sparsity.dense(2, 2).get_lower(), Sparsity.dense(2, 2).get_upper()
        0 1 3
        0 2 3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.get_lower(sp), ca.get_upper(sp)
        ([0, 1], [0, 2, 3])
        julia> ca.get_lower(ca.Sparsity.dense(2, 2)), ca.get_upper(ca.Sparsity.dense(2, 2))
        ([0, 1, 3], [0, 2, 3])
        \enddoctest
        \endlang

        \sa get_upper, nnz_lower, tril

        \identifier{2nw} */
    std::vector<casadi_int> get_lower() const;

    /** \brief Nonzero indices of the entries in the upper triangular part (row <= col)

        \sa get_lower, nnz_upper, triu

        \identifier{2nx} */
    std::vector<casadi_int> get_upper() const;

    /** \brief The pattern in compressed column storage: the colind and row vectors

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> colind, row = sp.get_ccs()
        >>> colind, row
        ([0, 2, 2, 3, 4], [0, 2, 1, 0])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> [colind, row] = sp.get_ccs(); disp(colind), disp(row)
        0 2 2 3 4
        0 2 1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> colind, row = ca.get_ccs(sp);
        julia> colind, row
        ([0, 2, 2, 3, 4], [0, 2, 1, 0])
        \enddoctest
        \endlang

        \sa get_crs, get_triplet, get_colind, get_row

        \identifier{2ny} */
    void get_ccs(std::vector<casadi_int>& SWIG_OUTPUT(colind),
                std::vector<casadi_int>& SWIG_OUTPUT(row)) const;

    /** \brief The pattern in compressed row storage: the rowind and col vectors

        Equivalent to the compressed column storage of the transpose. Note that the
        resulting order of the nonzeros (row-major) differs from the storage order of the
        pattern.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> rowind, col = sp.get_crs()
        >>> rowind, col
        ([0, 2, 3, 4], [0, 3, 2, 0])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> [rowind, col] = sp.get_crs(); disp(rowind), disp(col)
        0 2 3 4
        0 3 2 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> rowind, col = ca.get_crs(sp);
        julia> rowind, col
        ([0, 2, 3, 4], [0, 3, 2, 0])
        \enddoctest
        \endlang

        \sa get_ccs, T

        \identifier{2nz} */
    void get_crs(std::vector<casadi_int>& SWIG_OUTPUT(rowind),
                std::vector<casadi_int>& SWIG_OUTPUT(col)) const;

    /** \brief The pattern in sparse triplet format: the row and column of every nonzero

        Inverse of triplet(): Sparsity.triplet(nrow, ncol, row, col) rebuilds the pattern.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> row, col = sp.get_triplet()
        >>> row, col
        ([0, 2, 1, 0], [0, 0, 2, 3])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> [row, col] = sp.get_triplet(); disp(row), disp(col)
        0 2 1 0
        0 0 2 3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> row, col = ca.get_triplet(sp);
        julia> row, col
        ([0, 2, 1, 0], [0, 0, 2, 3])
        \enddoctest
        \endlang

        \sa triplet, get_row, get_col

        \identifier{2o0} */
    void get_triplet(std::vector<casadi_int>& SWIG_OUTPUT(row),
                    std::vector<casadi_int>& SWIG_OUTPUT(col)) const;

    /** \brief Pattern of the submatrix with rows \p rr and columns \p cc

        This is the sparsity counterpart of slicing a matrix as A[rr, cc]. Indices may
        repeat and may be negative (counting from the end); an index out of range raises an
        error.

        \param rr Row indices of the submatrix
        \param cc Column indices of the submatrix
        \param ind1 If true, indices are one-based
        \return The len(rr)-by-len(cc) pattern and a mapping with one entry per nonzero
                of the result: nonzero k of the submatrix is nonzero mapping[k] of this
                pattern

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sub, mapping = sp.sub([0, 2], [0, 3])
        >>> DM(sub, 1), mapping
        (DM(
        [[1, 1],
         [1, 00]]), [0, 1, 3])
        >>> sub, mapping = sp.sub([-1], [0])
        >>> sub, mapping
        (Sparsity(1x1), [1])
        >>> sub, mapping = sp.sub([0, 0], [0])
        >>> sub, mapping
        (Sparsity(2x1), [0, 0])
        >>> try: sp.sub([3], [0])
        ... except Exception as e: print("out of range")
        out of range
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> [sub, mapping] = sp.sub([0, 2], [0, 3]); disp(DM(sub, 1)), disp(mapping)
        [[1, 1],
         [1, 00]]
        0 1 3
        >> [sub, mapping] = sp.sub(-1, 0); disp(sub), disp(mapping)
        1x1
        1
        >> [sub, mapping] = sp.sub([0, 0], 0); disp(sub), disp(mapping)
        2x1
        0 0
        >> try, sp.sub(3, 0); catch e, disp('out of range'); end
        out of range
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.sub(sp, [0, 2], [0, 3])
        (2x2,3nz, [0, 1, 3])
        julia> ca.sub(sp, [-1], [0])
        (1x1, [1])
        julia> ca.sub(sp, [0, 0], [0])
        (2x1, [0, 0])
        julia> try
                   ca.sub(sp, [3], [0])
               catch e
                   println("out of range")
               end
        out of range
        \enddoctest
        \endlang

        \sa get_nz, erase, pmult

        \identifier{ci} */
    Sparsity sub(const std::vector<casadi_int>& rr,
                 const std::vector<casadi_int>& cc,
                 std::vector<casadi_int>& SWIG_OUTPUT(mapping), bool ind1=false) const;

    /** \brief Pattern of a set of elements, arranged according to a given pattern

        \p rr lists one flat (column-major) element index per nonzero of \p sp. The result
        has the shape of \p sp and keeps the nonzeros of \p sp whose element in this
        pattern is a structural nonzero. This is the sparsity counterpart of indexing a
        matrix with a Sparsity or with a list of element indices.

        \param rr Flat element indices, length sp.nnz(); negative indices count from the
               end, out of range raises an error
        \param sp Pattern that lays out the result
        \param ind1 If true, indices are one-based
        \return The pattern and a mapping: nonzero k of the result is nonzero mapping[k]
                of this pattern

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sub, mapping = sp.sub([0, 2, 9], Sparsity.dense(3, 1))
        >>> sub, mapping
        (Sparsity(3x1), [0, 1, 3])
        >>> sub, mapping = sp.sub([0, 1, 2, 3], Sparsity.dense(2, 2))
        >>> DM(sub, 1), mapping
        (DM(
        [[1, 1],
         [00, 00]]), [0, 1])
        >>> try: sp.sub([0, 2], Sparsity.dense(3, 1))
        ... except Exception as e: print("length mismatch")
        length mismatch
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> [sub, mapping] = sp.sub([0, 2, 9], Sparsity.dense(3, 1)); disp(sub), disp(mapping)
        3x1
        0 1 3
        >> [sub, mapping] = sp.sub([0, 1, 2, 3], Sparsity.dense(2, 2));
        >> DM(sub, 1), mapping
        [[1, 1],
         [00, 00]]
        0 1
        >> try, sp.sub([0, 2], Sparsity.dense(3, 1)); catch e, disp('length mismatch'); end
        length mismatch
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.sub(sp, [0, 2, 9], ca.Sparsity.dense(3, 1))
        (3x1, [0, 1, 3])
        julia> ca.sub(sp, [0, 1, 2, 3], ca.Sparsity.dense(2, 2))
        (2x2,2nz, [0, 1])
        julia> try
                   ca.sub(sp, [0, 2], ca.Sparsity.dense(3, 1))
               catch e
                   println("length mismatch")
               end
        length mismatch
        \enddoctest
        \endlang

        \sa get_nz, find

        \identifier{cj} */
    Sparsity sub(const std::vector<casadi_int>& rr, const Sparsity& sp,
                 std::vector<casadi_int>& SWIG_OUTPUT(mapping), bool ind1=false) const;

    /** \brief Transposed pattern

        In Python, T is a property (sp.T); the method transpose() also returns the nonzero
        reordering.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> DM(sp.T, 1)
        DM(
        [[1, 00, 1],
         [00, 00, 00],
         [00, 1, 00],
         [1, 00, 00]])
        >>> sp.T.colind(), sp.T.row()
        ([0, 2, 3, 4], [0, 3, 2, 0])
        >>> Sparsity.dense(2, 3).T, Sparsity(0, 3).T
        (Sparsity(3x2), Sparsity(3x0))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> DM(sp.T, 1)
        [[1, 00, 1],
         [00, 00, 00],
         [00, 1, 00],
         [1, 00, 00]]
        >> sp.T.colind(), sp.T.row()
        0 2 3 4
        0 3 2 0
        >> Sparsity.dense(2, 3).T, Sparsity(0, 3).T
        3x2
        3x0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.DM(ca.T(sp), 1)
        DM(
        [[1, 00, 1],
         [00, 00, 00],
         [00, 1, 00],
         [1, 00, 00]])
        julia> ca.colind(ca.T(sp)), ca.row(ca.T(sp))
        ([0, 2, 3, 4], [0, 3, 2, 0])
        julia> ca.T(ca.Sparsity.dense(2, 3)), ca.T(ca.Sparsity(0, 3))
        (3x2, 3x0)
        \enddoctest
        \endlang

        \sa transpose, is_transpose, get_crs

        \identifier{2o1} */
    Sparsity T() const;

    /** \brief Transposed pattern together with the reordering of the nonzeros

        \param invert_mapping Direction of the mapping
        \return The transposed pattern and a mapping of length nnz(). By default
                mapping[k] is the nonzero of the original pattern that becomes nonzero k
                of the transpose; with \p invert_mapping true, mapping[k] is the nonzero of
                the transpose that original nonzero k becomes.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> spT, mapping = sp.transpose()
        >>> spT, mapping
        (Sparsity(4x3,4nz), [0, 3, 2, 1])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> [spT, mapping] = sp.transpose(); disp(spT), disp(mapping)
        4x3,4nz
        0 3 2 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.transpose(sp, false)
        (4x3,4nz, [0, 3, 2, 1])
        \enddoctest
        \endlang

        \sa T

        \identifier{ck} */
    Sparsity transpose(std::vector<casadi_int>& SWIG_OUTPUT(mapping),
                        bool invert_mapping=false) const;

    /** \brief Check whether \p y is the transpose of this pattern

        \lang{python}
        \doctest
        >>> Sparsity.band(3, 1).is_transpose(Sparsity.band(3, -1))
        True
        >>> Sparsity.band(3, 1).is_transpose(Sparsity.band(3, 1))
        False
        >>> Sparsity.dense(2, 3).is_transpose(Sparsity.dense(3, 2))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.band(3, 1).is_transpose(Sparsity.band(3, -1))
        1
        >> Sparsity.band(3, 1).is_transpose(Sparsity.band(3, 1))
        0
        >> Sparsity.dense(2, 3).is_transpose(Sparsity.dense(3, 2))
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_transpose(ca.Sparsity.band(3, 1), ca.Sparsity.band(3, -1))
        true
        julia> ca.is_transpose(ca.Sparsity.band(3, 1), ca.Sparsity.band(3, 1))
        false
        julia> ca.is_transpose(ca.Sparsity.dense(2, 3), ca.Sparsity.dense(3, 2))
        true
        \enddoctest
        \endlang

        \sa T, is_symmetric

        \identifier{2o2} */
    bool is_transpose(const Sparsity& y) const;

    /** \brief Check whether \p y is a reshape of this pattern

        True when both have the same number of elements and the same flat (column-major)
        nonzero positions, i.e. reshape(A, y.size()) has pattern \p y.

        \lang{python}
        \doctest
        >>> Sparsity.dense(2, 3).is_reshape(Sparsity.dense(6, 1))
        True
        >>> Sparsity.diag(2).is_reshape(Sparsity.nonzeros(4, 1, [0, 3]))
        True
        >>> Sparsity.diag(2).is_reshape(Sparsity.nonzeros(4, 1, [0, 2]))
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.dense(2, 3).is_reshape(Sparsity.dense(6, 1))
        1
        >> Sparsity.diag(2).is_reshape(Sparsity.nonzeros(4, 1, [1, 4]))
        1
        >> Sparsity.diag(2).is_reshape(Sparsity.nonzeros(4, 1, [1, 3]))
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_reshape(ca.Sparsity.dense(2, 3), ca.Sparsity.dense(6, 1))
        true
        julia> ca.is_reshape(ca.Sparsity.diag(2), ca.Sparsity.nonzeros(4, 1, [0, 3]))
        true
        julia> ca.is_reshape(ca.Sparsity.diag(2), ca.Sparsity.nonzeros(4, 1, [0, 2]))
        false
        \enddoctest
        \endlang

        \sa reshape, find

        \identifier{2o3} */
    bool is_reshape(const Sparsity& y) const;

    /** \brief Check if the nonzero pattern is the Cartesian product of a row and a column subset

        If true, returns the (sorted) row and column index sets such that
        nz(i, j) <-> i in row and j in col. The nonzero buffer of a compactible
        matrix in CCS is laid out exactly as a column-major dense
        row.size() x col.size() matrix, so it can be consumed by dense kernels
        without gather/scatter. The index sets are returned in either case
        (the nonempty rows and columns); a pattern without nonzeros is compactible.

        \lang{python}
        \doctest
        >>> ok, row, col = Sparsity.rowcol([0, 2], [0, 2], 3, 3).is_compactible()
        >>> ok, row, col
        (True, [0, 2], [0, 2])
        >>> ok, row, col = Sparsity.diag(3).is_compactible()
        >>> ok, row, col
        (False, [0, 1, 2], [0, 1, 2])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [ok, row, col] = Sparsity.rowcol([0, 2], [0, 2], 3, 3).is_compactible();
        >> ok, row, col
        1
        0 2
        0 2
        >> [ok, row, col] = Sparsity.diag(3).is_compactible(); disp(ok), disp(row), disp(col)
        0
        0 1 2
        0 1 2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_compactible(ca.Sparsity.rowcol([0, 2], [0, 2], 3, 3))
        (true, [0, 2], [0, 2])
        julia> ca.is_compactible(ca.Sparsity.diag(3))
        (false, [0, 1, 2], [0, 1, 2])
        \enddoctest
        \endlang

        \sa rowcol

        \identifier{2gf} */
    bool is_compactible(std::vector<casadi_int>& SWIG_OUTPUT(row),
                        std::vector<casadi_int>& SWIG_OUTPUT(col)) const;

    /// @{
    /** \brief Combine two patterns of the same shape, as the result of a binary operation

        Models the pattern of f(x, y) applied elementwise, where \p y has this pattern's
        shape. An entry is kept in the result when both arguments are nonzero; an entry
        where only y is nonzero is kept unless \p f0x_is_zero (f(0, y) == 0) and an entry
        where only x is nonzero is kept unless \p function0_is_zero (f(x, 0) == 0).
        So combine(y, False, False) is the union (unite) and combine(y, True, True) is the
        intersection (intersect). A shape mismatch raises an error.

        The C++ overload with \p mapping returns, per nonzero of the union, bit 1 if the
        first argument is nonzero, bit 2 if the second is and bit 4 if the entry was
        dropped from the result.

        \param y Pattern with the same shape
        \param f0x_is_zero True if entries where only \p y is nonzero are structural zeros
        \param function0_is_zero True if entries where only this pattern is nonzero are
               structural zeros

        \lang{python}
        \doctest
        >>> a = Sparsity.diag(2)
        >>> b = Sparsity.triplet(2, 2, [0, 1], [1, 1])
        >>> a.combine(b, False, False) == a + b, a.combine(b, True, True) == a * b
        (True, True)
        >>> DM(a.combine(b, True, False), 1)
        DM(
        [[1, 00],
         [00, 1]])
        >>> DM(a.combine(b, False, True), 1)
        DM(
        [[00, 1],
         [00, 1]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> a = Sparsity.diag(2);
        >> b = Sparsity.triplet(2, 2, [0, 1], [1, 1]);
        >> a.combine(b, false, false) == a + b, a.combine(b, true, true) == a .* b
        1
        1
        >> DM(a.combine(b, true, false), 1)
        [[1, 00],
         [00, 1]]
        >> DM(a.combine(b, false, true), 1)
        [[00, 1],
         [00, 1]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> a = ca.Sparsity.diag(2);
        julia> b = ca.Sparsity.triplet(2, 2, [0, 1], [1, 1]);
        julia> ca.is_equal(ca.combine(a, b, false, false), ca.unite(a, b))
        true
        julia> ca.DM(ca.combine(a, b, true, false), 1)
        DM(
        [[1, 00],
         [00, 1]])
        julia> ca.DM(ca.combine(a, b, false, true), 1)
        DM(
        [[00, 1],
         [00, 1]])
        \enddoctest
        \endlang

        \sa unite, intersect

        \identifier{cl} */
#ifndef SWIG
    Sparsity combine(const Sparsity& y, bool f0x_is_zero, bool function0_is_zero,
                            std::vector<unsigned char>& mapping) const;
#endif // SWIG
    Sparsity combine(const Sparsity& y, bool f0x_is_zero, bool function0_is_zero) const;
    /// @}

    /// @{
    /** \brief Union of two sparsity patterns of the same shape (also available as the + operator)

        The pattern of x + y for matrices with these patterns. A shape mismatch raises an
        error.

        \lang{python}
        \doctest
        >>> a = Sparsity.diag(2)
        >>> b = Sparsity.triplet(2, 2, [0, 1], [1, 1])
        >>> DM(a.unite(b), 1)
        DM(
        [[1, 1],
         [00, 1]])
        >>> a + b == a.unite(b), a * b == a.intersect(b)
        (True, True)
        >>> try: a + Sparsity.dense(3, 3)
        ... except Exception as e: print("dimension mismatch")
        dimension mismatch
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> a = Sparsity.diag(2);
        >> b = Sparsity.triplet(2, 2, [0, 1], [1, 1]);
        >> DM(a.unite(b), 1)
        [[1, 1],
         [00, 1]]
        >> a + b == a.unite(b), a .* b == a.intersect(b)
        1
        1
        >> try, a + Sparsity.dense(3, 3); catch e, disp('dimension mismatch'); end
        dimension mismatch
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> a = ca.Sparsity.diag(2);
        julia> b = ca.Sparsity.triplet(2, 2, [0, 1], [1, 1]);
        julia> ca.DM(ca.unite(a, b), 1)
        DM(
        [[1, 1],
         [00, 1]])
        julia> try
                   ca.unite(a, ca.Sparsity.dense(3, 3))
               catch e
                   println("dimension mismatch")
               end
        dimension mismatch
        \enddoctest
        \endlang

        \sa intersect, combine, is_subset

        \identifier{cm} */
#ifndef SWIG
    Sparsity unite(const Sparsity& y, std::vector<unsigned char>& mapping) const;
#endif // SWIG
    Sparsity unite(const Sparsity& y) const;
    Sparsity operator+(const Sparsity& b) const;
    /// @}

    /// @{
    /** \brief Intersection of two sparsity patterns of the same shape (also the * operator)

        The pattern of the elementwise product x * y for matrices with these patterns.
        A shape mismatch raises an error. In MATLAB this is the .* operator; * and mtimes
        give the pattern of the matrix product.

        \lang{python}
        \doctest
        >>> a = Sparsity.diag(2)
        >>> b = Sparsity.triplet(2, 2, [0, 1], [1, 1])
        >>> DM(a.intersect(b), 1)
        DM(
        [[00, 00],
         [00, 1]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> a = Sparsity.diag(2);
        >> b = Sparsity.triplet(2, 2, [0, 1], [1, 1]);
        >> DM(a.intersect(b), 1)
        [[00, 00],
         [00, 1]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> a = ca.Sparsity.diag(2);
        julia> b = ca.Sparsity.triplet(2, 2, [0, 1], [1, 1]);
        julia> ca.DM(ca.intersect(a, b), 1)
        DM(
        [[00, 00],
         [00, 1]])
        \enddoctest
        \endlang

        \sa unite, combine

        \identifier{cn} */
#ifndef SWIG
    Sparsity intersect(const Sparsity& y,
                       std::vector<unsigned char>& mapping) const;
#endif // SWIG
    Sparsity intersect(const Sparsity& y) const;
    Sparsity operator*(const Sparsity& b) const;
    /// @}

    /** \brief Check whether every nonzero of this pattern is also a nonzero of \p rhs

        Both patterns must have the same shape, otherwise an error is raised. A pattern is
        a subset of itself and an all-zero pattern is a subset of any pattern of its shape.

        \lang{python}
        \doctest
        >>> Sparsity.diag(3).is_subset(Sparsity.dense(3, 3))
        True
        >>> Sparsity.dense(3, 3).is_subset(Sparsity.diag(3))
        False
        >>> Sparsity(3, 3).is_subset(Sparsity.diag(3))
        True
        >>> Sparsity.diag(3).is_subset(Sparsity.diag(3))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.diag(3).is_subset(Sparsity.dense(3, 3))
        1
        >> Sparsity.dense(3, 3).is_subset(Sparsity.diag(3))
        0
        >> Sparsity(3, 3).is_subset(Sparsity.diag(3))
        1
        >> Sparsity.diag(3).is_subset(Sparsity.diag(3))
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_subset(ca.Sparsity.diag(3), ca.Sparsity.dense(3, 3))
        true
        julia> ca.is_subset(ca.Sparsity.dense(3, 3), ca.Sparsity.diag(3))
        false
        julia> ca.is_subset(ca.Sparsity(3, 3), ca.Sparsity.diag(3))
        true
        julia> ca.is_subset(ca.Sparsity.diag(3), ca.Sparsity.diag(3))
        true
        \enddoctest
        \endlang

        \sa unite, intersect, is_equal

        \identifier{2o4} */
    bool is_subset(const Sparsity& rhs) const;

    /** \brief Propagate a subset of a pattern through a sparsity cast

        \p X and \p Y have the same number of nonzeros and sparsity_cast(X, Y) maps nonzero k
        of X to nonzero k of Y. Given that this pattern is a subset of \p X, the result is
        the subset of \p Y made of the nonzeros of \p Y whose counterpart in \p X is in this
        pattern.

        \lang{python}
        \doctest
        >>> X = Sparsity.diag(4)
        >>> Y = Sparsity.dense(2, 2)
        >>> x = Sparsity.nonzeros(4, 4, [0, 10, 15])
        >>> DM(x.sparsity_cast_mod(X, Y), 1)
        DM(
        [[1, 1],
         [00, 1]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> X = Sparsity.diag(4);
        >> Y = Sparsity.dense(2, 2);
        >> x = Sparsity.nonzeros(4, 4, [1, 11, 16]);
        >> DM(x.sparsity_cast_mod(X, Y), 1)
        [[1, 1],
         [00, 1]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> X = ca.Sparsity.diag(4);
        julia> Y = ca.Sparsity.dense(2, 2);
        julia> x = ca.Sparsity.nonzeros(4, 4, [0, 10, 15]);
        julia> ca.DM(ca.sparsity_cast_mod(x, X, Y), 1)
        DM(
        [[1, 1],
         [00, 1]])
        \enddoctest
        \endlang

        \sa sparsity_cast, is_subset

        \identifier{246} */
    Sparsity sparsity_cast_mod(const Sparsity& X, const Sparsity& Y) const;

    /** \brief Complement of the pattern: structural zeros become nonzeros and vice versa

        The result has the same shape; a dense pattern becomes all-zero and an all-zero
        pattern becomes dense.

        \lang{python}
        \doctest
        >>> DM(Sparsity.diag(2).pattern_inverse(), 1)
        DM(
        [[00, 1],
         [1, 00]])
        >>> Sparsity.dense(2, 2).pattern_inverse(), Sparsity(2, 2).pattern_inverse()
        (Sparsity(2x2,0nz), Sparsity(2x2))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.diag(2).pattern_inverse(), 1)
        [[00, 1],
         [1, 00]]
        >> Sparsity.dense(2, 2).pattern_inverse(), Sparsity(2, 2).pattern_inverse()
        2x2,0nz
        2x2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.pattern_inverse(ca.Sparsity.diag(2)), 1)
        DM(
        [[00, 1],
         [1, 00]])
        julia> ca.pattern_inverse(ca.Sparsity.dense(2, 2)), ca.pattern_inverse(ca.Sparsity(2, 2))
        (2x2,0nz, 2x2)
        \enddoctest
        \endlang

        \sa unite, intersect

        \identifier{2o5} */
    Sparsity pattern_inverse() const;

#ifndef SWIG
    /** \brief Propagate sparsity using 0-1 logic through a matrix product,

     * no memory allocation: <tt>z = mul(x, y)</tt> with work vector
     * Forward mode.

        \identifier{co} */
    static void mul_sparsityF(const bvec_t* x, const Sparsity& x_sp,
                              const bvec_t* y, const Sparsity& y_sp,
                              bvec_t* z, const Sparsity& z_sp,
                              bvec_t* w);

    /** \brief Propagate signal activity through a matrix product, forward mode
     *
     * Like mul_sparsityF, but combines the two factors of each term with AND rather
     * than OR (an inactive factor annihilates that term), so an output is active only
     * where some term has both factors active. No memory allocation; work vector w.

        \identifier{2in} */
    static void mul_activityF(const bvec_t* x, const Sparsity& x_sp,
                              const bvec_t* y, const Sparsity& y_sp,
                              bvec_t* z, const Sparsity& z_sp,
                              bvec_t* w);

    /** \brief Propagate sparsity using 0-1 logic through a matrix product,

     * no memory allocation: <tt>z = mul(x, y)</tt> with work vector
     * Reverse mode.

        \identifier{cp} */
    static void mul_sparsityR(bvec_t* x, const Sparsity& x_sp,
                              bvec_t* y, const Sparsity& y_sp,
                              bvec_t* z, const Sparsity& z_sp,
                              bvec_t* w);

    /// \cond INTERNAL
    /// @{
    /** \brief Accessed by SparsityInterface

        \identifier{cq} */
    static Sparsity horzcat(const std::vector<Sparsity> & sp);
    static Sparsity vertcat(const std::vector<Sparsity> & sp);
    static Sparsity blockcat(const std::vector< std::vector< Sparsity > > &v);
    static Sparsity diagcat(const std::vector< Sparsity > &v);
    static std::vector<Sparsity>
      horzsplit(const Sparsity& x, const std::vector<casadi_int>& offset);
    static std::vector<Sparsity>
      vertsplit(const Sparsity& x, const std::vector<casadi_int>& offset);
    static std::vector<Sparsity>
      diagsplit(const Sparsity& x,
                const std::vector<casadi_int>& offset1,
                const std::vector<casadi_int>& offset2);
    static Sparsity mtimes(const Sparsity& x, const Sparsity& y,
                           const std::string& blas = "reference");
    static Sparsity mac(const Sparsity& x, const Sparsity& y, const Sparsity& z,
                        const std::string& /*blas*/ = "reference") { return z;}
    static Sparsity reshape(const Sparsity& x, casadi_int nrow, casadi_int ncol);
    static Sparsity reshape(const Sparsity& x, const Sparsity& sp);
    static Sparsity sparsity_cast(const Sparsity& x, const Sparsity& sp);
    static casadi_int sprank(const Sparsity& x);
    static casadi_int norm_0_mul(const Sparsity& x, const Sparsity& A);
    static Sparsity kron(const Sparsity& a, const Sparsity& b);

    /** \brief Output sparsity of casadi::KronContract
     *
     * Given an (mA*mB, nA*nB) pattern `sp_m` and an (xrow, xcol) pattern
     * `sp_x`, computes the output sparsity. Both modes share the same
     * `(yrow, ycol) = sp_m.size() / sp_x.size()` shape; what differs is
     * which axis pair indexes the output:
     *
     * - `inner = true`:  (i, j) ∈ y iff ∃ (r, s) ∈ sp_x with
     *                    (i*xrow + r, j*xcol + s) ∈ sp_m;
     *                    output indexed by the outer block of M.
     * - `inner = false`: (r, s) ∈ y iff ∃ (i, j) ∈ sp_x with
     *                    (i*yrow + r, j*ycol + s) ∈ sp_m;
     *                    output indexed by the inner block of M.

        \identifier{2hn} */
    static Sparsity kron_contract(const Sparsity& sp_m, const Sparsity& sp_x, bool inner);
    static Sparsity triu(const Sparsity& x, bool includeDiagonal=true);
    static Sparsity tril(const Sparsity& x, bool includeDiagonal=true);

    static Sparsity sum2(const Sparsity &x);
    static Sparsity sum1(const Sparsity &x);
    /// @}
    /// \endcond
#endif //SWIG

    /** \brief Enlarge the pattern in place by inserting empty rows and columns

        The existing nonzeros are kept and scattered: row i of the old pattern becomes row
        rr[i] and column j becomes column cc[j] of the new nrow-by-ncol pattern, so that
        B[rr, cc] has the old pattern. The new rows and columns are structural zeros.
        Nonzero indices are unchanged. The object is rebound to the new pattern.

        \param nrow New number of rows
        \param ncol New number of columns
        \param rr New row of every old row: length size1(), strictly increasing, entries
               in [0, nrow); negative entries count from the end
        \param cc New column of every old column: length size2(), strictly increasing,
               entries in [0, ncol)
        \param ind1 If true, indices are one-based

        A wrong length or an index out of range raises an error.

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(2, 2)
        >>> sp.enlarge(4, 3, [0, 3], [1, 2])
        >>> DM(sp, 1)
        DM(
        [[00, 1, 1],
         [00, 00, 00],
         [00, 00, 00],
         [00, 1, 1]])
        >>> try: Sparsity.dense(2, 2).enlarge(4, 3, [0], [1, 2])
        ... except Exception as e: print("wrong length")
        wrong length
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(2, 2);
        >> sp.enlarge(4, 3, [0, 3], [1, 2]);
        >> DM(sp, 1)
        [[00, 1, 1],
         [00, 00, 00],
         [00, 00, 00],
         [00, 1, 1]]
        >> try, Sparsity.dense(2, 2).enlarge(4, 3, 0, [1, 2]); catch e, disp('wrong length'); end
        wrong length
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(2, 2);
        julia> ca.enlarge(sp, 4, 3, [0, 3], [1, 2])
        julia> ca.DM(sp, 1)
        DM(
        [[00, 1, 1],
         [00, 00, 00],
         [00, 00, 00],
         [00, 1, 1]])
        julia> try
                   ca.enlarge(ca.Sparsity.dense(2, 2), 4, 3, [0], [1, 2])
               catch e
                   println("wrong length")
               end
        wrong length
        \enddoctest
        \endlang

        \sa enlargeRows, enlargeColumns, resize, erase

        \identifier{cr} */
    void enlarge(casadi_int nrow, casadi_int ncol,
                  const std::vector<casadi_int>& rr,
                  const std::vector<casadi_int>& cc, bool ind1=false);

    /** \brief Enlarge the pattern in place by inserting empty rows

        Same as enlarge with the columns left unchanged.

        \param nrow New number of rows
        \param rr New row of every old row, length size1(), strictly increasing
        \param ind1 If true, indices are one-based

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(2, 2); sp.enlargeRows(3, [0, 2]); print(DM(sp, 1))
        [[1, 1], 
         [00, 00], 
         [1, 1]]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(2, 2); sp.enlargeRows(3, [0, 2]); disp(DM(sp, 1))
        [[1, 1],
         [00, 00],
         [1, 1]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(2, 2); ca.enlargeRows(sp, 3, [0, 2]); ca.DM(sp, 1)
        DM(
        [[1, 1],
         [00, 00],
         [1, 1]])
        \enddoctest
        \endlang

        \sa enlarge, enlargeColumns

        \identifier{cs} */
    void enlargeRows(casadi_int nrow, const std::vector<casadi_int>& rr, bool ind1=false);

    /** \brief Enlarge the pattern in place by inserting empty columns

        Same as enlarge with the rows left unchanged.

        \param ncol New number of columns
        \param cc New column of every old column, length size2(), strictly increasing
        \param ind1 If true, indices are one-based

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(2, 2); sp.enlargeColumns(3, [0, 2]); print(DM(sp, 1))
        [[1, 00, 1], 
         [1, 00, 1]]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(2, 2); sp.enlargeColumns(3, [0, 2]); disp(DM(sp, 1))
        [[1, 00, 1],
         [1, 00, 1]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(2, 2); ca.enlargeColumns(sp, 3, [0, 2]); ca.DM(sp, 1)
        DM(
        [[1, 00, 1],
         [1, 00, 1]])
        \enddoctest
        \endlang

        \sa enlarge, enlargeRows

        \identifier{ct} */
    void enlargeColumns(casadi_int ncol, const std::vector<casadi_int>& cc, bool ind1=false);

    /** \brief Dense pattern of the same shape, with the flat index of every nonzero

        \return The dense pattern and a mapping of length nnz(): mapping[k] is the flat
                (column-major) element index of nonzero k, as returned by find()

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 2, [0, 2], [0, 1])
        >>> dense, mapping = sp.makeDense()
        >>> dense, mapping
        (Sparsity(3x2), [0, 5])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 2, [0, 2], [0, 1]);
        >> [dense, mapping] = sp.makeDense(); disp(dense), disp(mapping)
        3x2
        0 5
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 2, [0, 2], [0, 1]);
        julia> ca.makeDense(sp)
        (3x2, [0, 5])
        \enddoctest
        \endlang

        \sa find, dense

        \identifier{cu} */
    Sparsity makeDense(std::vector<casadi_int>& SWIG_OUTPUT(mapping)) const;

    /** \brief Turn the entries with rows \p rr and columns \p cc into structural zeros, in place

        The shape is unchanged: only the nonzeros in the block rr x cc are removed. An
        entry is removed when its row is in \p rr and its column is in \p cc, so an empty
        \p rr or \p cc removes nothing. Indices may be unsorted and negative (counting from
        the end); out of range raises an error. The object is rebound to the new pattern.

        \param rr Rows of the block to clear
        \param cc Columns of the block to clear
        \param ind1 If true, indices are one-based
        \return Mapping of length nnz() of the result: nonzero k of the result was nonzero
                mapping[k] of the original pattern

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(3, 3)
        >>> mapping = sp.erase([0, 2], [1])
        >>> DM(sp, 1), mapping
        (DM(
        [[1, 00, 1],
         [1, 1, 1],
         [1, 00, 1]]), [0, 1, 2, 4, 6, 7, 8])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(3, 3);
        >> mapping = sp.erase([0, 2], 1);
        >> DM(sp, 1), mapping
        [[1, 00, 1],
         [1, 1, 1],
         [1, 00, 1]]
        0 1 2 4 6 7 8
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(3, 3);
        julia> mapping = ca.erase(sp, [0, 2], [1]);
        julia> println(mapping)
        [0, 1, 2, 4, 6, 7, 8]
        julia> ca.DM(sp, 1)
        DM(
        [[1, 00, 1],
         [1, 1, 1],
         [1, 00, 1]])
        \enddoctest
        \endlang

        \sa erase(rr), sub, enlarge

        \identifier{cv} */
    std::vector<casadi_int> erase(const std::vector<casadi_int>& rr,
                                  const std::vector<casadi_int>& cc, bool ind1=false);

    /** \brief Turn the entries with flat (column-major) indices \p rr into structural zeros

        The shape is unchanged. Indices may be unsorted and negative; out of range raises an
        error. The object is rebound to the new pattern.

        \param rr Flat element indices, in [0, numel())
        \param ind1 If true, indices are one-based
        \return Mapping of length nnz() of the result: nonzero k of the result was nonzero
                mapping[k] of the original pattern

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(3, 3)
        >>> sp.erase([0, 4, 8]), sp
        ([1, 2, 3, 5, 6, 7], Sparsity(3x3,6nz))
        >>> DM(sp, 1)
        DM(
        [[00, 1, 1],
         [1, 00, 1],
         [1, 1, 00]])
        >>> try: Sparsity.dense(3, 3).erase([9])
        ... except Exception as e: print("out of range")
        out of range
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(3, 3);
        >> sp.erase([0, 4, 8]), sp
        1 2 3 5 6 7
        3x3,6nz
        >> DM(sp, 1)
        [[00, 1, 1],
         [1, 00, 1],
         [1, 1, 00]]
        >> try, Sparsity.dense(3, 3).erase(9); catch e, disp('out of range'); end
        out of range
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(3, 3);
        julia> println(ca.erase(sp, [0, 4, 8]))
        [1, 2, 3, 5, 6, 7]
        julia> ca.DM(sp, 1)
        DM(
        [[00, 1, 1],
         [1, 00, 1],
         [1, 1, 00]])
        julia> try
                   ca.erase(ca.Sparsity.dense(3, 3), [9])
               catch e
                   println("out of range")
               end
        out of range
        \enddoctest
        \endlang

        \sa erase(rr, cc), find

        \identifier{cw} */
    std::vector<casadi_int> erase(const std::vector<casadi_int>& rr, bool ind1=false);

    /** \brief Append a pattern below this one, in place (vertical concatenation)

        The number of columns must match, otherwise an error is raised; a 0-by-0 pattern on
        either side is ignored. Equivalent to rebinding to vertcat(self, sp); only cheap
        when both are column vectors.

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(2, 2); sp.append(Sparsity.diag(2)); print(DM(sp, 1))
        [[1, 1], 
         [1, 1], 
         [1, 00], 
         [00, 1]]
        >>> sp = Sparsity(0, 0); sp.append(Sparsity.diag(2)); print(sp)
        2x2,2nz
        >>> try: Sparsity.dense(2, 2).append(Sparsity.dense(2, 3))
        ... except Exception as e: print("column mismatch")
        column mismatch
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(2, 2); sp.append(Sparsity.diag(2)); disp(DM(sp, 1))
        [[1, 1],
         [1, 1],
         [1, 00],
         [00, 1]]
        >> sp = Sparsity(0, 0); sp.append(Sparsity.diag(2)); disp(sp)
        2x2,2nz
        >> try, sp.append(Sparsity.dense(2, 3)); catch e, disp('column mismatch'); end
        column mismatch
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(2, 2); ca.append(sp, ca.Sparsity.diag(2)); ca.DM(sp, 1)
        DM(
        [[1, 1],
         [1, 1],
         [1, 00],
         [00, 1]])
        julia> sp = ca.Sparsity(0, 0); ca.append(sp, ca.Sparsity.diag(2)); print(sp)
        2x2,2nz
        julia> try
                   ca.append(ca.Sparsity.dense(2, 2), ca.Sparsity.dense(2, 3))
               catch e
                   println("column mismatch")
               end
        column mismatch
        \enddoctest
        \endlang

        \sa appendColumns, vertcat

        \identifier{2o6} */
    void append(const Sparsity& sp);

    /** \brief Append a pattern to the right of this one, in place (horizontal concatenation)

        The number of rows must match, otherwise an error is raised; a 0-by-0 pattern on
        either side is ignored. Equivalent to rebinding to horzcat(self, sp).

        \lang{python}
        \doctest
        >>> sp = Sparsity.dense(2, 2); sp.appendColumns(Sparsity.diag(2)); print(DM(sp, 1))
        [[1, 1, 1, 00], 
         [1, 1, 00, 1]]
        >>> try: Sparsity.dense(2, 2).appendColumns(Sparsity.dense(3, 2))
        ... except Exception as e: print("row mismatch")
        row mismatch
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.dense(2, 2); sp.appendColumns(Sparsity.diag(2)); disp(DM(sp, 1))
        [[1, 1, 1, 00],
         [1, 1, 00, 1]]
        >> try, sp.appendColumns(Sparsity.dense(3, 2)); catch e, disp('row mismatch'); end
        row mismatch
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.dense(2, 2); ca.appendColumns(sp, ca.Sparsity.diag(2)); ca.DM(sp, 1)
        DM(
        [[1, 1, 1, 00],
         [1, 1, 00, 1]])
        julia> try
                   ca.appendColumns(ca.Sparsity.dense(2, 2), ca.Sparsity.dense(3, 2))
               catch e
                   println("row mismatch")
               end
        row mismatch
        \enddoctest
        \endlang

        \sa append, horzcat

        \identifier{2o7} */
    void appendColumns(const Sparsity& sp);

    /** \brief Check whether the pattern is 1-by-1

        \param scalar_and_dense If true, additionally require the single entry to be a
               nonzero (a 1-by-1 pattern with zero nonzeros is a structural zero)

        \lang{python}
        \doctest
        >>> Sparsity.dense(1, 1).is_scalar(), Sparsity(1, 1).is_scalar()
        (True, True)
        >>> Sparsity(1, 1).is_scalar(True), Sparsity.dense(2, 1).is_scalar()
        (False, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.dense(1, 1).is_scalar(), Sparsity(1, 1).is_scalar()]
        1 1
        >> [Sparsity(1, 1).is_scalar(true), Sparsity.dense(2, 1).is_scalar()]
        0 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_scalar(ca.Sparsity.dense(1, 1)), ca.is_scalar(ca.Sparsity(1, 1))
        (true, true)
        julia> ca.is_scalar(ca.Sparsity(1, 1), true), ca.is_scalar(ca.Sparsity.dense(2, 1))
        (false, false)
        \enddoctest
        \endlang

        \sa scalar, is_dense

        \identifier{2o8} */
    bool is_scalar(bool scalar_and_dense=false) const;

    /** \brief Check whether every entry is a structural nonzero, i.e. nnz() == numel()

        An empty pattern (a zero dimension) is dense. Note that a DM built from a nested
        list is dense even if it holds numerical zeros.

        \lang{python}
        \doctest
        >>> Sparsity.dense(2, 3).is_dense(), Sparsity.diag(3).is_dense()
        (True, False)
        >>> Sparsity(0, 3).is_dense(), Sparsity(1, 1).is_dense()
        (True, False)
        >>> DM([[1, 0], [0, 1]]).sparsity().is_dense()
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.dense(2, 3).is_dense(), Sparsity.diag(3).is_dense()]
        1 0
        >> [Sparsity(0, 3).is_dense(), Sparsity(1, 1).is_dense()]
        1 0
        >> sparsity(DM([1, 0; 0, 1])).is_dense()
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_dense(ca.Sparsity.dense(2, 3)), ca.is_dense(ca.Sparsity.diag(3))
        (true, false)
        julia> ca.is_dense(ca.Sparsity(0, 3)), ca.is_dense(ca.Sparsity(1, 1))
        (true, false)
        julia> ca.is_dense(ca.sparsity(ca.DM([1 0; 0 1])))
        true
        \enddoctest
        \endlang

        \sa dense, density, makeDense

        \identifier{2o9} */
    bool is_dense() const;

    /** \brief Check whether the pattern is a row vector, i.e. size1() == 1

        The number of columns (including zero) and the nonzeros play no role.

        \lang{python}
        \doctest
        >>> Sparsity.dense(1, 3).is_row(), Sparsity.dense(3, 1).is_row()
        (True, False)
        >>> Sparsity(1, 0).is_row()
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.dense(1, 3).is_row(), Sparsity.dense(3, 1).is_row()]
        1 0
        >> Sparsity(1, 0).is_row()
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_row(ca.Sparsity.dense(1, 3)), ca.is_row(ca.Sparsity.dense(3, 1))
        (true, false)
        julia> ca.is_row(ca.Sparsity(1, 0))
        true
        \enddoctest
        \endlang

        \sa is_column, is_vector

        \identifier{cx} */
    bool is_row() const;

    /** \brief Check whether the pattern is a column vector, i.e. size2() == 1

        \lang{python}
        \doctest
        >>> Sparsity.dense(3, 1).is_column(), Sparsity.dense(1, 3).is_column()
        (True, False)
        >>> Sparsity(1, 1).is_column()
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.dense(3, 1).is_column(), Sparsity.dense(1, 3).is_column()]
        1 0
        >> Sparsity(1, 1).is_column()
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_column(ca.Sparsity.dense(3, 1)), ca.is_column(ca.Sparsity.dense(1, 3))
        (true, false)
        julia> ca.is_column(ca.Sparsity(1, 1))
        true
        \enddoctest
        \endlang

        \sa is_row, is_vector

        \identifier{cy} */
    bool is_column() const;

    /** \brief Check whether the pattern is a row or a column vector

        A 1-by-1 pattern is a vector; a 0-by-3 pattern is not.

        \lang{python}
        \doctest
        >>> Sparsity.dense(3, 1).is_vector(), Sparsity.dense(1, 3).is_vector()
        (True, True)
        >>> Sparsity.dense(1, 1).is_vector(), Sparsity.dense(2, 2).is_vector()
        (True, False)
        >>> Sparsity(0, 3).is_vector()
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.dense(3, 1).is_vector(), Sparsity.dense(1, 3).is_vector()]
        1 1
        >> [Sparsity.dense(1, 1).is_vector(), Sparsity.dense(2, 2).is_vector()]
        1 0
        >> Sparsity(0, 3).is_vector()
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_vector(ca.Sparsity.dense(3, 1)), ca.is_vector(ca.Sparsity.dense(1, 3))
        (true, true)
        julia> ca.is_vector(ca.Sparsity.dense(1, 1)), ca.is_vector(ca.Sparsity.dense(2, 2))
        (true, false)
        julia> ca.is_vector(ca.Sparsity(0, 3))
        false
        \enddoctest
        \endlang

        \sa is_row, is_column

        \identifier{cz} */
    bool is_vector() const;

    /** \brief Check whether the pattern is square with exactly the full diagonal as nonzeros

        True only for Sparsity.diag(n) (including n == 0); a rectangular diagonal pattern,
        an all-zero square pattern or a pattern with off-diagonal entries is not diagonal.

        \lang{python}
        \doctest
        >>> Sparsity.diag(3).is_diag(), Sparsity.diag(2, 3).is_diag()
        (True, False)
        >>> Sparsity.dense(1, 1).is_diag(), Sparsity(2, 2).is_diag()
        (True, False)
        >>> Sparsity(0, 0).is_diag(), Sparsity.dense(2, 2).is_diag()
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.diag(3).is_diag(), Sparsity.diag(2, 3).is_diag()]
        1 0
        >> [Sparsity.dense(1, 1).is_diag(), Sparsity(2, 2).is_diag()]
        1 0
        >> [Sparsity(0, 0).is_diag(), Sparsity.dense(2, 2).is_diag()]
        1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_diag(ca.Sparsity.diag(3)), ca.is_diag(ca.Sparsity.diag(2, 3))
        (true, false)
        julia> ca.is_diag(ca.Sparsity.dense(1, 1)), ca.is_diag(ca.Sparsity(2, 2))
        (true, false)
        julia> ca.is_diag(ca.Sparsity(0, 0)), ca.is_diag(ca.Sparsity.dense(2, 2))
        (true, false)
        \enddoctest
        \endlang

        \sa diag, nnz_diag

        \identifier{247} */
    bool is_diag() const;

    /** \brief Check whether the number of rows equals the number of columns

        \lang{python}
        \doctest
        >>> Sparsity.dense(2, 2).is_square(), Sparsity.dense(2, 3).is_square()
        (True, False)
        >>> Sparsity(0, 0).is_square()
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.dense(2, 2).is_square(), Sparsity.dense(2, 3).is_square()]
        1 0
        >> Sparsity(0, 0).is_square()
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_square(ca.Sparsity.dense(2, 2)), ca.is_square(ca.Sparsity.dense(2, 3))
        (true, false)
        julia> ca.is_square(ca.Sparsity(0, 0))
        true
        \enddoctest
        \endlang

        \identifier{248} */
    bool is_square() const;

    /** \brief Check whether the pattern equals its transpose

        Structural symmetry only: (i, j) is a nonzero if and only if (j, i) is. Required by
        amd, ldl and the star colorings. A rectangular pattern is never symmetric.

        \lang{python}
        \doctest
        >>> Sparsity.banded(3, 1).is_symmetric(), Sparsity.upper(3).is_symmetric()
        (True, False)
        >>> Sparsity(3, 3).is_symmetric(), Sparsity.dense(2, 3).is_symmetric()
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.banded(3, 1).is_symmetric(), Sparsity.upper(3).is_symmetric()]
        1 0
        >> [Sparsity(3, 3).is_symmetric(), Sparsity.dense(2, 3).is_symmetric()]
        1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_symmetric(ca.Sparsity.banded(3, 1)), ca.is_symmetric(ca.Sparsity.upper(3))
        (true, false)
        julia> ca.is_symmetric(ca.Sparsity(3, 3)), ca.is_symmetric(ca.Sparsity.dense(2, 3))
        (true, false)
        \enddoctest
        \endlang

        \sa is_transpose, T

        \identifier{249} */
    bool is_symmetric() const;

    /** \brief Check whether all nonzeros are on or above the diagonal (row <= col)

        \param strictly If true, require row < col (no diagonal entries)

        Defined for any shape; a pattern without nonzeros is upper triangular.

        \lang{python}
        \doctest
        >>> Sparsity.upper(3).is_triu(), Sparsity.upper(3).is_triu(True)
        (True, False)
        >>> Sparsity.band(3, -1).is_triu(True), Sparsity(3, 3).is_triu()
        (True, True)
        >>> Sparsity.dense(1, 3).is_triu(), Sparsity.dense(3, 1).is_triu()
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.upper(3).is_triu(), Sparsity.upper(3).is_triu(true)]
        1 0
        >> [Sparsity.band(3, -1).is_triu(true), Sparsity(3, 3).is_triu()]
        1 1
        >> [Sparsity.dense(1, 3).is_triu(), Sparsity.dense(3, 1).is_triu()]
        1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_triu(ca.Sparsity.upper(3)), ca.is_triu(ca.Sparsity.upper(3), true)
        (true, false)
        julia> ca.is_triu(ca.Sparsity.band(3, -1), true), ca.is_triu(ca.Sparsity(3, 3))
        (true, true)
        julia> ca.is_triu(ca.Sparsity.dense(1, 3)), ca.is_triu(ca.Sparsity.dense(3, 1))
        (true, false)
        \enddoctest
        \endlang

        \sa is_tril, upper, triu, nnz_upper

        \identifier{24a} */
    bool is_triu(bool strictly = false) const;

    /** \brief Check whether all nonzeros are on or below the diagonal (row >= col)

        \param strictly If true, require row > col (no diagonal entries)

        \lang{python}
        \doctest
        >>> Sparsity.lower(3).is_tril(), Sparsity.lower(3).is_tril(True)
        (True, False)
        >>> Sparsity.band(3, 1).is_tril(True), Sparsity.diag(3).is_tril()
        (True, True)
        >>> Sparsity.dense(3, 1).is_tril(), Sparsity.dense(1, 3).is_tril()
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.lower(3).is_tril(), Sparsity.lower(3).is_tril(true)]
        1 0
        >> [Sparsity.band(3, 1).is_tril(true), Sparsity.diag(3).is_tril()]
        1 1
        >> [Sparsity.dense(3, 1).is_tril(), Sparsity.dense(1, 3).is_tril()]
        1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_tril(ca.Sparsity.lower(3)), ca.is_tril(ca.Sparsity.lower(3), true)
        (true, false)
        julia> ca.is_tril(ca.Sparsity.band(3, 1), true), ca.is_tril(ca.Sparsity.diag(3))
        (true, true)
        julia> ca.is_tril(ca.Sparsity.dense(3, 1)), ca.is_tril(ca.Sparsity.dense(1, 3))
        (true, false)
        \enddoctest
        \endlang

        \sa is_triu, lower, tril, nnz_lower

        \identifier{24b} */
    bool is_tril(bool strictly = false) const;

    /** \brief Check whether a square pattern is structurally singular

        True when the structural rank (sprank) is smaller than the size, i.e. no choice
        of nonzero values can make a matrix with this pattern invertible. A non-square
        pattern raises an error.

        \lang{python}
        \doctest
        >>> Sparsity.diag(3).is_singular(), Sparsity(3, 3).is_singular()
        (False, True)
        >>> Sparsity.nonzeros(3, 3, [0, 1, 2]).is_singular(), Sparsity(0, 0).is_singular()
        (True, False)
        >>> try: Sparsity.dense(2, 3).is_singular()
        ... except Exception as e: print("not square")
        not square
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Sparsity.diag(3).is_singular(), Sparsity(3, 3).is_singular()]
        0 1
        >> [Sparsity.nonzeros(3, 3, [1, 2, 3]).is_singular(), Sparsity(0, 0).is_singular()]
        1 0
        >> try, Sparsity.dense(2, 3).is_singular(); catch e, disp('not square'); end
        not square
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_singular(ca.Sparsity.diag(3)), ca.is_singular(ca.Sparsity(3, 3))
        (false, true)
        julia> ca.is_singular(ca.Sparsity.nonzeros(3, 3, [0, 1, 2]))
        true
        julia> ca.is_singular(ca.Sparsity(0, 0))
        false
        julia> try
                   ca.is_singular(ca.Sparsity.dense(2, 3))
               catch e
                   println("not square")
               end
        not square
        \enddoctest
        \endlang

        \sa sprank, btf

        \identifier{24c} */
    bool is_singular() const;

    /** \brief Check whether the pattern is that of a permutation matrix

        A square pattern with exactly one nonzero in every row and every column, so that
        multiplying a dense vector permutes its entries. Equivalent to is_orthonormal().

        \lang{python}
        \doctest
        >>> Sparsity.permutation([1, 2, 0]).is_permutation()
        True
        >>> Sparsity.diag(3).is_permutation()
        True
        >>> Sparsity.dense(2, 2).is_permutation(), Sparsity.diag(2, 3).is_permutation()
        (False, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.permutation([1, 2, 0]).is_permutation()
        1
        >> Sparsity.diag(3).is_permutation()
        1
        >> [Sparsity.dense(2, 2).is_permutation(), Sparsity.diag(2, 3).is_permutation()]
        0 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_permutation(ca.Sparsity.permutation([1, 2, 0]))
        true
        julia> ca.is_permutation(ca.Sparsity.diag(3))
        true
        julia> ca.is_permutation(ca.Sparsity.dense(2, 2)), ca.is_permutation(ca.Sparsity.diag(2, 3))
        (false, false)
        \enddoctest
        \endlang

        \sa permutation, permutation_vector, is_selection, is_orthonormal

        \identifier{24d} */
    bool is_permutation() const;

    /** \brief Check whether the pattern is that of a selection matrix

        A selection matrix S picks a subset of the entries of a dense vector: S @ v is a
        vector of distinct entries of v. The pattern has at most one nonzero per column and
        exactly one per row. Equivalent to is_orthonormal_rows(allow_empty).

        \param allow_empty If true, rows without nonzeros are allowed (the result may
               contain structural zeros)

        \lang{python}
        \doctest
        >>> S = Sparsity.triplet(2, 3, [0, 1], [2, 0])
        >>> DM(S, 1) @ DM([10, 20, 30])
        DM([30, 10])
        >>> S.is_selection(), Sparsity.triplet(2, 3, [0], [2]).is_selection()
        (True, False)
        >>> Sparsity.triplet(2, 3, [0], [2]).is_selection(True)
        True
        >>> Sparsity.diag(3, 2).is_selection()
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> S = Sparsity.triplet(2, 3, [0, 1], [2, 0]);
        >> DM(S, 1) * DM([10; 20; 30])
        [30, 10]
        >> [S.is_selection(), Sparsity.triplet(2, 3, 0, 2).is_selection()]
        1 0
        >> Sparsity.triplet(2, 3, 0, 2).is_selection(true)
        1
        >> Sparsity.diag(3, 2).is_selection()
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> S = ca.Sparsity.triplet(2, 3, [0, 1], [2, 0]);
        julia> ca.DM(S, 1) * ca.DM([10, 20, 30])
        DM([30, 10])
        julia> ca.is_selection(S), ca.is_selection(ca.Sparsity.triplet(2, 3, [0], [2]))
        (true, false)
        julia> ca.is_selection(ca.Sparsity.triplet(2, 3, [0], [2]), true)
        true
        julia> ca.is_selection(ca.Sparsity.diag(3, 2))
        false
        \enddoctest
        \endlang

        \sa is_permutation, is_orthonormal_rows

        \identifier{24e} */
    bool is_selection(bool allow_empty=false) const;

    /** \brief Check whether every row and every column has at most one nonzero

        Without \p allow_empty the pattern must also be square with exactly one nonzero
        per row and column (a permutation pattern). With \p allow_empty, rows and columns
        without nonzeros are disregarded and the shape may be rectangular.

        \param allow_empty Disregard empty rows and columns in the analysis

        \lang{python}
        \doctest
        >>> Sparsity.diag(3).is_orthonormal()
        True
        >>> Sparsity.triplet(3, 3, [0, 1], [0, 1]).is_orthonormal()
        False
        >>> Sparsity.triplet(3, 3, [0, 1], [0, 1]).is_orthonormal(True)
        True
        >>> Sparsity.triplet(3, 3, [0, 0], [0, 1]).is_orthonormal(True)
        False
        >>> Sparsity.diag(2, 3).is_orthonormal(), Sparsity.diag(2, 3).is_orthonormal(True)
        (False, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.diag(3).is_orthonormal()
        1
        >> Sparsity.triplet(3, 3, [0, 1], [0, 1]).is_orthonormal()
        0
        >> Sparsity.triplet(3, 3, [0, 1], [0, 1]).is_orthonormal(true)
        1
        >> Sparsity.triplet(3, 3, [0, 0], [0, 1]).is_orthonormal(true)
        0
        >> [Sparsity.diag(2, 3).is_orthonormal(), Sparsity.diag(2, 3).is_orthonormal(true)]
        0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.is_orthonormal(ca.Sparsity.diag(3))
        true
        julia> ca.is_orthonormal(ca.Sparsity.triplet(3, 3, [0, 1], [0, 1]))
        false
        julia> ca.is_orthonormal(ca.Sparsity.triplet(3, 3, [0, 1], [0, 1]), true)
        true
        julia> ca.is_orthonormal(ca.Sparsity.triplet(3, 3, [0, 0], [0, 1]), true)
        false
        julia> sp = ca.Sparsity.diag(2, 3);
        julia> ca.is_orthonormal(sp), ca.is_orthonormal(sp, true)
        (false, true)
        \enddoctest
        \endlang

        \sa is_orthonormal_rows, is_orthonormal_columns, is_permutation

        \identifier{24f} */
    bool is_orthonormal(bool allow_empty=false) const;

    /** \brief Check whether the rows are orthonormal: one nonzero per row, at most one per column

        Without \p allow_empty the pattern must have at most as many rows as columns and
        exactly one nonzero per row. With \p allow_empty, rows without nonzeros are
        disregarded.

        \param allow_empty Disregard empty rows in the analysis

        \lang{python}
        \doctest
        >>> sp = Sparsity.diag(2, 3)
        >>> sp.is_orthonormal_rows(), sp.is_orthonormal_columns()
        (True, False)
        >>> sp.T.is_orthonormal_rows(), sp.T.is_orthonormal_columns()
        (False, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.diag(2, 3);
        >> [sp.is_orthonormal_rows(), sp.is_orthonormal_columns()]
        1 0
        >> [sp.T.is_orthonormal_rows(), sp.T.is_orthonormal_columns()]
        0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.diag(2, 3);
        julia> ca.is_orthonormal_rows(sp), ca.is_orthonormal_columns(sp)
        (true, false)
        julia> ca.is_orthonormal_rows(ca.T(sp)), ca.is_orthonormal_columns(ca.T(sp))
        (false, true)
        \enddoctest
        \endlang

        \sa is_orthonormal_columns, is_selection

        \identifier{24g} */
    bool is_orthonormal_rows(bool allow_empty=false) const;

    /** \brief Check whether the columns are orthonormal (see is_orthonormal_rows)

        Without \p allow_empty the pattern must have at most as many columns as rows and
        exactly one nonzero per column. With \p allow_empty, columns without nonzeros are
        disregarded.

        \param allow_empty Disregard empty columns in the analysis

        \sa is_orthonormal_rows, is_orthonormal

        \identifier{24h} */
    bool is_orthonormal_columns(bool allow_empty=false) const;

    /** \brief Check whether the row indices are increasing within every column

        A pattern constructed through the public interface always satisfies this: row
        indices are sorted and unique within a column.

        \param strictly If true (default), duplicate rows within a column are not allowed

        \lang{python}
        \doctest
        >>> Sparsity.banded(3, 1).rowsSequential()
        True
        >>> Sparsity.dense(2, 2).rowsSequential(False)
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.banded(3, 1).rowsSequential()
        1
        >> Sparsity.dense(2, 2).rowsSequential(false)
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.rowsSequential(ca.Sparsity.banded(3, 1))
        true
        julia> ca.rowsSequential(ca.Sparsity.dense(2, 2), false)
        true
        \enddoctest
        \endlang

        \sa removeDuplicates

        \identifier{d0} */
    bool rowsSequential(bool strictly=true) const;

    /** \brief Remove duplicate entries, in place

        Patterns created through the public interface have no duplicates, so this is a
        no-op for them. The entries of \p mapping (one per nonzero) corresponding to the
        removed duplicates are removed as well.

        \param mapping Vector of length nnz(), filtered in the same way as the nonzeros

        \lang{python}
        \doctest
        >>> sp = Sparsity.banded(3, 1)
        >>> m = sp.removeDuplicates([0, 1, 2, 3, 4, 5, 6])
        >>> m, sp
        ([0, 1, 2, 3, 4, 5, 6], Sparsity(3x3,7nz))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.banded(3, 1);
        >> m = sp.removeDuplicates([0, 1, 2, 3, 4, 5, 6]);
        >> m, sp
        0 1 2 3 4 5 6
        3x3,7nz
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.banded(3, 1);
        julia> m = ca.removeDuplicates(sp, [0, 1, 2, 3, 4, 5, 6]);
        julia> m, sp
        ([0, 1, 2, 3, 4, 5, 6], 3x3,7nz)
        \enddoctest
        \endlang

        \sa rowsSequential, triplet

        \identifier{d1} */
    void removeDuplicates(std::vector<casadi_int>& SWIG_INOUT(mapping));

#ifndef SWIG
    typedef std::unordered_multimap<std::size_t, WeakRef> CachingMap;

    /// Cached sparsity patterns
    static CachingMap& getCache();

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    // Safe access to CachingMap
    static std::mutex cachingmap_mtx;
#endif //CASADI_WITH_THREADSAFE_SYMBOLICS

    /// (Dense) scalar
    static const Sparsity& getScalar();

    /// (Sparse) scalar
    static const Sparsity& getScalarSparse();

    /// Empty zero-by-zero
    static const Sparsity& getEmpty();

#endif //SWIG

    /** \brief Elimination tree of the pattern

        Returns, for every column, the index of its parent in the elimination tree, or -1
        for a root. See Direct Methods for Sparse Linear Systems by Davis (2006).
        With \p ata false the result matches MATLAB's etree(A) with zero-based indices
        (the pattern is treated as symmetric, only the upper triangular part is used);
        with \p ata true it matches etree(A, 'col'), the tree of A'*A.

        The implementation is a modified version of cs_etree in CSparse
        Copyright(c) Timothy A. Davis, 2006-2009
        Licensed as a derivative work under the GNU LGPL

        \lang{python}
        \doctest
        >>> Sparsity.banded(4, 1).etree(), Sparsity.diag(4).etree()
        ([1, 2, 3, -1], [-1, -1, -1, -1])
        >>> Sparsity.dense(2, 3).etree(True)
        [1, 2, -1]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.banded(4, 1).etree(), Sparsity.diag(4).etree()
        1 2 3 -1
        -1 -1 -1 -1
        >> Sparsity.dense(2, 3).etree(true)
        1 2 -1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.etree(ca.Sparsity.banded(4, 1)), ca.etree(ca.Sparsity.diag(4))
        ([1, 2, 3, -1], [-1, -1, -1, -1])
        julia> println(ca.etree(ca.Sparsity.dense(2, 3), true))
        [1, 2, -1]
        \enddoctest
        \endlang

        \sa ldl, qr_sparse, amd

        \identifier{d2} */
    std::vector<casadi_int> etree(bool ata=false) const;

    /** \brief Symbolic LDL factorization of a symmetric pattern

        Computes the pattern of the factor in A[p, p] = L*D*L^T, where L is unit lower
        triangular. Returns the pattern of the strictly upper triangular L^T (the unit
        diagonal is not stored) and the ordering p. A non-symmetric pattern raises an error.
        The numerical counterpart is ldl() applied to a matrix.

        \param amd If true (default), a fill-reducing AMD ordering is computed and returned
               in p; if false, p is the identity
        \return Pattern of L^T and the permutation vector p

        \lang{python}
        \doctest
        >>> Lt, p = Sparsity.banded(4, 1).ldl(False)
        >>> DM(Lt, 1), p
        (DM(
        [[00, 1, 00, 00],
         [00, 00, 1, 00],
         [00, 00, 00, 1],
         [00, 00, 00, 00]]), [0, 1, 2, 3])
        >>> A = Sparsity.nonzeros(4, 4, [0, 3, 5, 10, 12, 15])
        >>> Lt, p = A.ldl()
        >>> Lt.nnz(), p
        (1, [1, 2, 0, 3])
        >>> try: Sparsity.upper(3).ldl()
        ... except Exception as e: print("not symmetric")
        not symmetric
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [Lt, p] = Sparsity.banded(4, 1).ldl(false);
        >> DM(Lt, 1), p
        [[00, 1, 00, 00],
         [00, 00, 1, 00],
         [00, 00, 00, 1],
         [00, 00, 00, 00]]
        0 1 2 3
        >> A = Sparsity.nonzeros(4, 4, [1, 4, 6, 11, 13, 16]);
        >> [Lt, p] = A.ldl(); disp(Lt.nnz()), disp(p)
        1
        1 2 0 3
        >> try, Sparsity.upper(3).ldl(); catch e, disp('not symmetric'); end
        not symmetric
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.ldl(ca.Sparsity.banded(4, 1), false)
        (4x4,3nz, [0, 1, 2, 3])
        julia> A = ca.Sparsity.nonzeros(4, 4, [0, 3, 5, 10, 12, 15]);
        julia> ca.nnz(first(ca.ldl(A)))
        1
        julia> try
                   ca.ldl(ca.Sparsity.upper(3))
               catch e
                   println("not symmetric")
               end
        not symmetric
        \enddoctest
        \endlang

        The implementation is a modified version of LDL
        Copyright(c) Timothy A. Davis, 2005-2013
        Licensed as a derivative work under the GNU LGPL

        \sa amd, etree, qr_sparse

        \identifier{d3} */
    Sparsity ldl(std::vector<casadi_int>& SWIG_OUTPUT(p), bool amd=true) const;

    /** \brief Symbolic QR factorization

        Computes the patterns of V (the Householder vectors, a compact representation of Q)
        and R for the factorization of this pattern with its columns permuted by pc, plus
        the row permutation prinv used by the numerical factorization (qr_sparse on a
        matrix) and solution. V and R may have more rows than the pattern (the matrix is
        extended with zero rows when it has fewer rows than columns in the structural sense).

        \param amd If true (default), the columns are reordered with the AMD ordering of
               A^T*A, returned in pc; if false, pc is the identity
        \return V, R, prinv, pc

        \lang{python}
        \doctest
        >>> V, R, prinv, pc = Sparsity.dense(3, 2).qr_sparse()
        >>> V, R, prinv, pc
        (Sparsity(3x2,5nz), Sparsity(3x2,3nz), [0, 1, 2], [0, 1])
        >>> V, R, prinv, pc = Sparsity.banded(4, 1).qr_sparse(False)
        >>> V, R, prinv, pc
        (Sparsity(4x4,7nz), Sparsity(4x4,9nz), [0, 1, 2, 3], [0, 1, 2, 3])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> [V, R, prinv, pc] = Sparsity.dense(3, 2).qr_sparse();
        >> V, R, prinv, pc
        3x2,5nz
        3x2,3nz
        0 1 2
        0 1
        >> [V, R, prinv, pc] = Sparsity.banded(4, 1).qr_sparse(false);
        >> V, R, prinv, pc
        4x4,7nz
        4x4,9nz
        0 1 2 3
        0 1 2 3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> V, R, prinv, pc = ca.qr_sparse(ca.Sparsity.dense(3, 2));
        julia> V, R, prinv, pc
        (3x2,5nz, 3x2,3nz, [0, 1, 2], [0, 1])
        julia> V, R, prinv, pc = ca.qr_sparse(ca.Sparsity.banded(4, 1), false);
        julia> V, R, prinv, pc
        (4x4,7nz, 4x4,9nz, [0, 1, 2, 3], [0, 1, 2, 3])
        \enddoctest
        \endlang

        The implementation is a modified version of CSparse
        Copyright(c) Timothy A. Davis, 2006-2009
        Licensed as a derivative work under the GNU LGPL

        \sa ldl, amd, etree

        \identifier{d4} */
    void qr_sparse(Sparsity& SWIG_OUTPUT(V), Sparsity& SWIG_OUTPUT(R),
                   std::vector<casadi_int>& SWIG_OUTPUT(prinv),
                   std::vector<casadi_int>& SWIG_OUTPUT(pc), bool amd=true) const;

    /** \brief Depth-first search on the adjacency graph of the sparsity

        Low-level building block of scc (cs_dfs in CSparse). Starting from column \p j, it
        visits the columns reachable through the nonzeros and writes them, in reverse
        finishing order, into \p xi below position \p top. See Direct Methods for Sparse
        Linear Systems by Davis (2006).

        \param j Start node (column)
        \param top Current top of the output stack in \p xi
        \param xi Output stack, length at least size2(); the visited nodes are stored in
               xi[new_top:top]
        \param pstack Work vector, length at least size2()
        \param pinv Inverse row permutation, or an empty vector for none
        \param marked Visited flags per node, updated in place
        \return The new top of the stack, with the updated xi, pstack and marked vectors

        \lang{python}
        \doctest
        >>> Sparsity.banded(3, 1).dfs(0, 3, [0, 0, 0], [0, 0, 0], [], [False] * 3)
        (0, [0, 1, 2], [1, 4, 5], [True, True, True])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.banded(3, 1); xi = [0, 0, 0]; pstack = [0, 0, 0];
        >> [top, xi, pstack, marked] = sp.dfs(0, 3, xi, pstack, [], [false, false, false]);
        >> top, xi, pstack, marked
        0
        0 1 2
        1 4 5
        1 1 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.dfs(ca.Sparsity.banded(3, 1), 0, 3, [0, 0, 0], [0, 0, 0], [], [false, false, false])
        (0, [0, 1, 2], [1, 4, 5], Bool[1, 1, 1])
        \enddoctest
        \endlang

        \sa scc

        \identifier{d5} */
    casadi_int dfs(casadi_int j, casadi_int top, std::vector<casadi_int>& SWIG_INOUT1(xi),
            std::vector<casadi_int>& SWIG_INOUT2(pstack),
            const std::vector<casadi_int>& pinv, std::vector<bool>& SWIG_INOUT3(marked)) const;

    /** \brief Strongly connected components of the directed graph of a square pattern

        Node j is the j-th row/column and nonzero (i, j) is an edge from j to i.
        See Direct Methods for Sparse Linear Systems by Davis (2006).

        \return The number of components n, the node indices grouped by component, and the
                offsets (length n+1): component i consists of the nodes
                index[offset[i]], ..., index[offset[i+1]-1], in increasing order.

        The components are ordered so that the symmetric permutation A[index, index] is
        block upper triangular with n diagonal blocks; for a symmetric pattern it is block
        diagonal, with the block boundaries given by offset.

        The implementation is a modified version of cs_scc in CSparse
        Copyright(c) Timothy A. Davis, 2006-2009
        Licensed as a derivative work under the GNU LGPL

        \lang{python}
        \doctest
        >>> A = Sparsity.nonzeros(4, 4, [0, 2, 5, 7, 10, 13, 15])
        >>> DM(A, 1)
        DM(
        [[1, 00, 00, 00],
         [00, 1, 00, 1],
         [1, 00, 1, 00],
         [00, 1, 00, 1]])
        >>> n, index, offset = A.scc()
        >>> n, index, offset
        (3, [2, 0, 1, 3], [0, 1, 2, 4])
        >>> DM(A, 1)[index, index]
        DM(
        [[1, 1, 00, 00],
         [00, 1, 00, 00],
         [00, 00, 1, 1],
         [00, 00, 1, 1]])
        >>> Sparsity.diag(3).scc(), Sparsity.dense(3, 3).scc(), Sparsity(0, 0).scc()
        ((3, [0, 1, 2], [0, 1, 2, 3]), (1, [0, 1, 2], [0, 3]), (0, [], [0]))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = Sparsity.nonzeros(4, 4, [1, 3, 6, 8, 11, 14, 16]);
        >> DM(A, 1)
        [[1, 00, 00, 00],
         [00, 1, 00, 1],
         [1, 00, 1, 00],
         [00, 1, 00, 1]]
        >> [n, index, offset] = A.scc(); disp(n), disp(index), disp(offset)
        3
        2 0 1 3
        0 1 2 4
        >> Af = DM(A, 1); disp(Af(index + 1, index + 1))
        [[1, 1, 00, 00],
         [00, 1, 00, 00],
         [00, 00, 1, 1],
         [00, 00, 1, 1]]
        >> [n, index, offset] = Sparsity.diag(3).scc(); disp(n), disp(index), disp(offset)
        3
        0 1 2
        0 1 2 3
        >> [n, index, offset] = Sparsity.dense(3, 3).scc(); disp(n), disp(index), disp(offset)
        1
        0 1 2
        0 3
        >> [n, index, offset] = Sparsity(0, 0).scc(); disp(n), disp(isempty(index)), disp(offset)
        0
        1
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.Sparsity.nonzeros(4, 4, [0, 2, 5, 7, 10, 13, 15]);
        julia> ca.DM(A, 1)
        DM(
        [[1, 00, 00, 00],
         [00, 1, 00, 1],
         [1, 00, 1, 00],
         [00, 1, 00, 1]])
        julia> ca.scc(A)
        (3, [2, 0, 1, 3], [0, 1, 2, 4])
        julia> ca.scc(ca.Sparsity.diag(3)), ca.scc(ca.Sparsity.dense(3, 3))
        ((3, [0, 1, 2], [0, 1, 2, 3]), (1, [0, 1, 2], [0, 3]))
        julia> ca.scc(ca.Sparsity(0, 0))
        (0, Int64[], [0])
        \enddoctest
        \endlang

        \sa btf, dfs

        \identifier{d6} */
    casadi_int scc(std::vector<casadi_int>& SWIG_OUTPUT(index),
            std::vector<casadi_int>& SWIG_OUTPUT(offset)) const;

    /** \brief Block triangular form (BTF) through the Dulmage-Mendelsohn decomposition

        Computes a row permutation and a column permutation that bring the pattern into
        block upper triangular form, see Direct Methods for Sparse Linear Systems by
        Davis (2006). The pattern may be rectangular.

        \return
          - the number of fine blocks nb
          - rowperm, colperm: row and column permutation vectors, so that
            A[rowperm, colperm] is in block upper triangular form
          - rowblock, colblock (length nb+1): block k spans rows rowblock[k] to
            rowblock[k+1]-1 and columns colblock[k] to colblock[k+1]-1 of the permuted
            matrix
          - coarse_rowblock, coarse_colblock (length 5): boundaries of the coarse
            decomposition into the underdetermined, square (well-determined) and
            overdetermined parts; coarse_colblock[3] is the structural rank (sprank)

        Only structure is used: the distance of entries to the diagonal plays no role, so
        a randomly permuted block-diagonal matrix is not necessarily brought back to block
        diagonal form, and the row and column permutations differ in general, which
        destroys symmetry. The result is cached in the pattern.

        The implementation is a modified version of cs_dmperm in CSparse
        Copyright(c) Timothy A. Davis, 2006-2009
        Licensed as a derivative work under the GNU LGPL

        \lang{python}
        \doctest
        >>> B = Sparsity.triplet(3, 3, [0, 0, 1, 2], [0, 2, 1, 2])
        >>> DM(B, 1)
        DM(
        [[1, 00, 1],
         [00, 1, 00],
         [00, 00, 1]])
        >>> nb, rowperm, colperm, rowblock, colblock, crb, ccb = B.btf()
        >>> nb, rowperm, colperm, rowblock, colblock
        (3, [2, 0, 1], [2, 0, 1], [0, 1, 2, 3], [0, 1, 2, 3])
        >>> crb, ccb
        ([0, 0, 0, 3, 3], [0, 0, 3, 3, 3])
        >>> DM(B, 1)[rowperm, colperm]
        DM(
        [[1, 00, 00],
         [1, 1, 00],
         [00, 00, 1]])
        >>> sprank(Sparsity.nonzeros(3, 3, [0, 1, 2])), sprank(Sparsity.diag(3))
        (1, 3)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> B = Sparsity.triplet(3, 3, [0, 0, 1, 2], [0, 2, 1, 2]);
        >> DM(B, 1)
        [[1, 00, 1],
         [00, 1, 00],
         [00, 00, 1]]
        >> [nb, rowperm, colperm, rowblock, colblock, crb, ccb] = B.btf();
        >> nb, rowperm, colperm, rowblock, colblock
        3
        2 0 1
        2 0 1
        0 1 2 3
        0 1 2 3
        >> crb, ccb
        0 0 0 3 3
        0 0 3 3 3
        >> Bf = DM(B, 1); disp(Bf(rowperm + 1, colperm + 1))
        [[1, 00, 00],
         [1, 1, 00],
         [00, 00, 1]]
        >> [sprank(Sparsity.nonzeros(3, 3, [1, 2, 3])), sprank(Sparsity.diag(3))]
        1 3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> B = ca.Sparsity.triplet(3, 3, [0, 0, 1, 2], [0, 2, 1, 2]);
        julia> ca.DM(B, 1)
        DM(
        [[1, 00, 1],
         [00, 1, 00],
         [00, 00, 1]])
        julia> ca.btf(B)
        (3, [2, 0, 1], [2, 0, 1], [0, 1, 2, 3], [0, 1, 2, 3], [0, 0, 0, 3, 3], [0, 0, 3, 3, 3])
        julia> ca.sprank(ca.Sparsity.nonzeros(3, 3, [0, 1, 2])), ca.sprank(ca.Sparsity.diag(3))
        (1, 3)
        \enddoctest
        \endlang

        \sa scc, sprank, is_singular

        \identifier{d7} */
    casadi_int btf(std::vector<casadi_int>& SWIG_OUTPUT(rowperm),
            std::vector<casadi_int>& SWIG_OUTPUT(colperm),
            std::vector<casadi_int>& SWIG_OUTPUT(rowblock),
            std::vector<casadi_int>& SWIG_OUTPUT(colblock),
            std::vector<casadi_int>& SWIG_OUTPUT(coarse_rowblock),
            std::vector<casadi_int>& SWIG_OUTPUT(coarse_colblock)) const;

    /** \brief Approximate minimum degree (AMD) ordering of a symmetric pattern

        Fill-reducing ordering p to apply to a symmetric linear system before factorization:
        factorize A[p, p] instead of A. The pattern must be symmetric, otherwise an error is
        raised; for an unsymmetric A, order the pattern of A^T*A (as qr_sparse does).
        Used by ldl and qr_sparse when their \p amd argument is true.

        The implementation is a modified version of cs_amd in CSparse
        Copyright(c) Timothy A. Davis, 2006-2009
        Licensed as a derivative work under the GNU LGPL

        \lang{python}
        \doctest
        >>> Sparsity.banded(4, 1).amd(), Sparsity.dense(3, 3).amd(), Sparsity(0, 0).amd()
        ([3, 2, 0, 1], [0, 1, 2], [])
        >>> try: Sparsity.upper(3).amd()
        ... except Exception as e: print("not symmetric")
        not symmetric
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.banded(4, 1).amd(), Sparsity.dense(3, 3).amd()
        3 2 0 1
        0 1 2
        >> isempty(Sparsity(0, 0).amd())
        1
        >> try, Sparsity.upper(3).amd(); catch e, disp('not symmetric'); end
        not symmetric
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.amd(ca.Sparsity.banded(4, 1)), ca.amd(ca.Sparsity.dense(3, 3))
        ([3, 2, 0, 1], [0, 1, 2])
        julia> ca.amd(ca.Sparsity(0, 0))
        Int64[]
        julia> try
                   ca.amd(ca.Sparsity.upper(3))
               catch e
                   println("not symmetric")
               end
        not symmetric
        \enddoctest
        \endlang

        \sa ldl, qr_sparse, etree, is_symmetric

        \identifier{d8} */
    std::vector<casadi_int> amd() const;

#ifndef SWIG
    /** \brief Propagate sparsity through a linear solve

        \identifier{d9} */
    void spsolve(bvec_t* X, bvec_t* B, bool tr) const;
#endif // SWIG

    /** \brief Flat (column-major) element index of every nonzero

        Returns, in nonzero order, the position each nonzero would have in a dense
        column-major matrix of the same shape: row + col*size1(). Indexing a dense matrix
        A of this shape with the result picks the entries that are nonzeros of the pattern.
        Inverse of nonzeros(). For a dense pattern this is 0, ..., numel()-1.

        \param ind1 If true, the indices are one-based

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sp.find(), sp.find(True)
        ([0, 2, 7, 9], [1, 3, 8, 10])
        >>> Sparsity.dense(2, 2).find(), Sparsity(2, 2).find(), Sparsity(0, 3).find()
        ([0, 1, 2, 3], [], [])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> sp.find(), sp.find(false)
        1 3 8 10
        0 2 7 9
        >> Sparsity.dense(2, 2).find()
        1 2 3 4
        >> [isempty(Sparsity(2, 2).find()), isempty(Sparsity(0, 3).find())]
        1 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.find(sp), ca.find(sp, true)
        ([0, 2, 7, 9], [1, 3, 8, 10])
        julia> println(ca.find(ca.Sparsity.dense(2, 2)))
        [0, 1, 2, 3]
        julia> ca.find(ca.Sparsity(2, 2)), ca.find(ca.Sparsity(0, 3))
        (Int64[], Int64[])
        \enddoctest
        \endlang

        \sa nonzeros, get_nz, makeDense

        \identifier{da} */
    std::vector<casadi_int> find(bool ind1=SWIG_IND1) const;

#ifndef SWIG
    /// Get the location of all nonzero elements (inplace version)
    void find(std::vector<casadi_int>& loc, bool ind1=false) const;
#endif // SWIG

    /** \brief Unidirectional (column) coloring for compressed Jacobian evaluation

        Greedy distance-2 coloring (Algorithm 3.1 in Gebremedhin, Manne and Pothen,
        What Color Is Your Jacobian?, SIAM Rev. 47(4), 2006): columns that have no row in
        common may share a color, so that one directional derivative per color recovers
        all nonzeros. This is the coloring used for forward-mode Jacobians (and, applied
        to the transpose, for reverse mode).

        \param AT The transpose of this pattern, if already available; a null Sparsity
               (default) computes it
        \param cutoff If more than this many colors are needed, a null Sparsity is
               returned instead
        \return A size2()-by-ncolors pattern: column c lists the columns of this pattern
                that received color c

        \lang{python}
        \doctest
        >>> A = Sparsity.triplet(3, 4, [0, 1, 2, 0], [0, 1, 2, 3])
        >>> DM(A, 1)
        DM(
        [[1, 00, 00, 1],
         [00, 1, 00, 00],
         [00, 00, 1, 00]])
        >>> c = A.uni_coloring()
        >>> DM(c, 1)
        DM(
        [[1, 00],
         [1, 00],
         [1, 00],
         [00, 1]])
        >>> c.size2(), [c.row()[c.colind()[j]:c.colind()[j + 1]] for j in range(c.size2())]
        (2, [[0, 1, 2], [3]])
        >>> A.uni_coloring(Sparsity(), 1).is_null()
        True
        >>> Sparsity.dense(3, 3).uni_coloring(), Sparsity.diag(3).uni_coloring()
        (Sparsity(3x3,3nz), Sparsity(3x1))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = Sparsity.triplet(3, 4, [0, 1, 2, 0], [0, 1, 2, 3]);
        >> DM(A, 1)
        [[1, 00, 00, 1],
         [00, 1, 00, 00],
         [00, 00, 1, 00]]
        >> c = A.uni_coloring();
        >> DM(c, 1)
        [[1, 00],
         [1, 00],
         [1, 00],
         [00, 1]]
        >> disp(c.size2()), ci = c.colind(); r = c.row();
        2
        >> for j = 1:c.size2(), disp(r(ci(j) + 1:ci(j + 1))), end
        0 1 2
        3
        >> A.uni_coloring(Sparsity(), 1).is_null()
        1
        >> Sparsity.dense(3, 3).uni_coloring(), Sparsity.diag(3).uni_coloring()
        3x3,3nz
        3x1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.Sparsity.triplet(3, 4, [0, 1, 2, 0], [0, 1, 2, 3]);
        julia> ca.DM(A, 1)
        DM(
        [[1, 00, 00, 1],
         [00, 1, 00, 00],
         [00, 00, 1, 00]])
        julia> c = ca.uni_coloring(A, ca.Sparsity(), 1000);
        julia> ca.DM(c, 1)
        DM(
        [[1, 00],
         [1, 00],
         [1, 00],
         [00, 1]])
        julia> ca.size2(c), [ca.row(c)[ca.colind(c)[j]+1:ca.colind(c)[j+1]] for j in 1:ca.size2(c)]
        (2, [[0, 1, 2], [3]])
        julia> ca.is_null(ca.uni_coloring(A, ca.Sparsity(), 1))
        true
        julia> ca.uni_coloring(ca.Sparsity.dense(3, 3), ca.Sparsity(), 1000)
        Sparsity(3x3,3nz)
        julia> ca.uni_coloring(ca.Sparsity.diag(3), ca.Sparsity(), 1000)
        Sparsity(3x1)
        \enddoctest
        \endlang

        \sa star_coloring, star_coloring_new

        \identifier{db} */
    Sparsity uni_coloring(const Sparsity& AT=Sparsity(),
                          casadi_int cutoff = std::numeric_limits<casadi_int>::max()) const;

    /** \brief Star coloring of a symmetric pattern, with the color that recovers each nonzero

        Front end to star_coloring and star_coloring2 for compressed Hessian evaluation.
        A non-symmetric pattern raises an error.

        \param opts Options: "new_algo" (bool, default false) selects star_coloring2
               instead of star_coloring; "largest_first" (bool, default true) selects the
               largest-first ordering; "cutoff" (int) returns a null pattern when more
               colors are needed
        \return The coloring (a size2()-by-ncolors pattern as in uni_coloring) and
                which_color, one entry per nonzero: the color whose directional derivative
                yields that nonzero, or -1 when the nonzero is obtained from its symmetric
                counterpart instead (the lowest color index is preferred)

        \lang{python}
        \doctest
        >>> H = Sparsity.banded(4, 1)
        >>> coloring, which_color = H.star_coloring_new()
        >>> DM(coloring, 1)
        DM(
        [[00, 1, 00],
         [1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]])
        >>> for k in range(H.nnz()): print(H.row()[k], H.get_col()[k], which_color[k])
        0 0 1
        1 0 -1
        0 1 0
        1 1 0
        2 1 0
        1 2 -1
        2 2 1
        3 2 1
        2 3 -1
        3 3 2
        >>> try: Sparsity.upper(3).star_coloring_new()
        ... except Exception as e: print("not symmetric")
        not symmetric
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> H = Sparsity.banded(4, 1);
        >> [coloring, which_color] = H.star_coloring_new();
        >> DM(coloring, 1)
        [[00, 1, 00],
         [1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]]
        >> [H.row(); H.get_col(); which_color]'
        0 0 1
        1 0 -1
        0 1 0
        1 1 0
        2 1 0
        1 2 -1
        2 2 1
        3 2 1
        2 3 -1
        3 3 2
        >> try, Sparsity.upper(3).star_coloring_new(); catch e, disp('not symmetric'); end
        not symmetric
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> H = ca.Sparsity.banded(4, 1);
        julia> ca.star_coloring_new(H)
        (4x3,4nz, [1, -1, 0, 0, 0, -1, 1, 1, -1, 2])
        julia> try
                   ca.star_coloring_new(ca.Sparsity.upper(3))
               catch e
                   println("not symmetric")
               end
        not symmetric
        \enddoctest
        \endlang

        \sa star_coloring, star_coloring2, uni_coloring

        \identifier{2jz} */
    Sparsity star_coloring_new(std::vector<casadi_int>& SWIG_OUTPUT(which_color),
        const Dict& opts=Dict()) const;

    /** \brief Star coloring of a symmetric pattern (greedy distance-2 coloring)

        Algorithm 4.1 in Gebremedhin, Manne and Pothen, What Color Is Your Jacobian?,
        SIAM Rev. 47(4), 2006. Exploits symmetry so that fewer colors than uni_coloring
        are needed to recover all nonzeros of a Hessian. A non-symmetric pattern raises an
        error.

        \param ordering 0 for the natural column order, 1 (default) to color in
               largest_first() order
        \param cutoff If more than this many colors are needed, a null Sparsity is
               returned
        \return A size2()-by-ncolors pattern: column c lists the columns with color c

        \lang{python}
        \doctest
        >>> H = Sparsity.banded(4, 1)
        >>> DM(H.star_coloring(), 1)
        DM(
        [[00, 1, 00],
         [1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]])
        >>> Sparsity.diag(5).star_coloring()
        Sparsity(5x1)
        >>> Sparsity.dense(4, 4).star_coloring(1, 2).is_null()
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> H = Sparsity.banded(4, 1);
        >> DM(H.star_coloring(), 1)
        [[00, 1, 00],
         [1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]]
        >> Sparsity.diag(5).star_coloring()
        5x1
        >> Sparsity.dense(4, 4).star_coloring(1, 2).is_null()
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> H = ca.Sparsity.banded(4, 1);
        julia> ca.DM(ca.star_coloring(H, 1, 1000), 1)
        DM(
        [[00, 1, 00],
         [1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]])
        julia> ca.star_coloring(ca.Sparsity.diag(5), 1, 1000)
        Sparsity(5x1)
        julia> ca.is_null(ca.star_coloring(ca.Sparsity.dense(4, 4), 0, 2))
        true
        \enddoctest
        \endlang

        \sa star_coloring_new, star_coloring2, uni_coloring, largest_first

        \identifier{dc} */
    Sparsity star_coloring(casadi_int ordering = 1,
                            casadi_int cutoff = std::numeric_limits<casadi_int>::max()) const;

    /** \brief Star coloring of a symmetric pattern with the newer algorithm

        Algorithm 4.1 in Gebremedhin, Tarafdar, Manne and Pothen, New Acyclic and Star
        Coloring Algorithms with Application to Computing Hessians, SIAM J. Sci. Comput.
        29(3), 2007. Same arguments and return value as star_coloring.

        \lang{python}
        \doctest
        >>> H = Sparsity.banded(4, 1)
        >>> DM(H.star_coloring2(), 1)
        DM(
        [[00, 1, 00],
         [1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> H = Sparsity.banded(4, 1);
        >> DM(H.star_coloring2(), 1)
        [[00, 1, 00],
         [1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> H = ca.Sparsity.banded(4, 1);
        julia> ca.DM(ca.star_coloring2(H, 1, 1000), 1)
        DM(
        [[00, 1, 00],
         [1, 00, 00],
         [00, 1, 00],
         [00, 00, 1]])
        \enddoctest
        \endlang

        \sa star_coloring, star_coloring_new

        \identifier{dd} */
    Sparsity star_coloring2(casadi_int ordering = 1,
                            casadi_int cutoff = std::numeric_limits<casadi_int>::max()) const;

    /** \brief Column indices ordered by decreasing number of nonzeros

        Columns with the same number of nonzeros keep their relative order. Used as the
        ordering of the star colorings.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sp.largest_first(), Sparsity.diag(3).largest_first()
        ([0, 2, 3, 1], [0, 1, 2])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> sp.largest_first(), Sparsity.diag(3).largest_first()
        0 2 3 1
        0 1 2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.largest_first(sp), ca.largest_first(ca.Sparsity.diag(3))
        ([0, 2, 3, 1], [0, 1, 2])
        \enddoctest
        \endlang

        \sa star_coloring

        \identifier{de} */
    std::vector<casadi_int> largest_first() const;

    /** \brief Permute the rows and/or columns of the pattern

        Nonzero (i, j) moves to (p[i], p[j]) when both rows and columns are permuted, to
        (p[i], j) for rows only and to (i, p[j]) for columns only. In matrix terms this is
        P^T*A*P, P^T*A or A*P with P the matrix of Sparsity.permutation(p). With
        \p invert_permutation the inverse of p is used instead, so that
        pmult(p, True, True, True) is the pattern of A[p, p].

        \param p Permutation vector; its length must equal the number of rows
               (if \p permute_rows) and the number of columns (if \p permute_columns)
        \param permute_rows Permute the rows
        \param permute_columns Permute the columns
        \param invert_permutation Use the inverse permutation

        \lang{python}
        \doctest
        >>> A = Sparsity.triplet(3, 3, [0, 1, 2, 0], [1, 2, 0, 0])
        >>> DM(A, 1)
        DM(
        [[1, 1, 00],
         [00, 00, 1],
         [1, 00, 00]])
        >>> p = [2, 0, 1]
        >>> DM(A.pmult(p), 1)
        DM(
        [[00, 1, 00],
         [00, 00, 1],
         [1, 00, 1]])
        >>> P = DM(Sparsity.permutation(p), 1)
        >>> A.pmult(p) == (P.T @ DM(A, 1) @ P).sparsity()
        True
        >>> A.pmult(p, True, True, True) == A.sub(p, p)[0]
        True
        >>> DM(A.pmult(p, True, False), 1)
        DM(
        [[00, 00, 1],
         [1, 00, 00],
         [1, 1, 00]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = Sparsity.triplet(3, 3, [0, 1, 2, 0], [1, 2, 0, 0]);
        >> DM(A, 1)
        [[1, 1, 00],
         [00, 00, 1],
         [1, 00, 00]]
        >> p = [2, 0, 1];
        >> DM(A.pmult(p), 1)
        [[00, 1, 00],
         [00, 00, 1],
         [1, 00, 1]]
        >> P = DM(Sparsity.permutation(p), 1);
        >> A.pmult(p) == sparsity(P' * DM(A, 1) * P)
        1
        >> [sub, mapping] = A.sub(p, p); disp(A.pmult(p, true, true, true) == sub)
        1
        >> DM(A.pmult(p, true, false), 1)
        [[00, 00, 1],
         [1, 00, 00],
         [1, 1, 00]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.Sparsity.triplet(3, 3, [0, 1, 2, 0], [1, 2, 0, 0]);
        julia> ca.DM(A, 1)
        DM(
        [[1, 1, 00],
         [00, 00, 1],
         [1, 00, 00]])
        julia> p = [2, 0, 1];
        julia> ca.DM(ca.pmult(A, p), 1)
        DM(
        [[00, 1, 00],
         [00, 00, 1],
         [1, 00, 1]])
        julia> P = ca.DM(ca.Sparsity.permutation(p), 1);
        julia> ca.is_equal(ca.pmult(A, p), ca.sparsity(P' * ca.DM(A, 1) * P))
        true
        julia> ca.is_equal(ca.pmult(A, p, true, true, true), first(ca.sub(A, p, p)))
        true
        julia> ca.DM(ca.pmult(A, p, true, false), 1)
        DM(
        [[00, 00, 1],
         [1, 00, 00],
         [1, 1, 00]])
        \enddoctest
        \endlang

        \sa permutation, sub

        \identifier{df} */
    Sparsity pmult(const std::vector<casadi_int>& p,
                    bool permute_rows=true, bool permute_columns=true,
                    bool invert_permutation=false) const;

    /** \brief Dimensions as a string "nrow x ncol", optionally with the number of nonzeros

        \param with_nz If true, append ",<nnz>nz"

        This is the string printed for a Sparsity: with the nonzero count when the pattern
        is not dense.

        \lang{python}
        \doctest
        >>> sp = Sparsity.diag(2, 3)
        >>> sp.dim(), sp.dim(True)
        ('2x3', '2x3,2nz')
        >>> Sparsity.dense(2, 2).dim(), Sparsity.dense(2, 2).dim(True)
        ('2x2', '2x2,4nz')
        >>> sp
        Sparsity(2x3,2nz)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.diag(2, 3);
        >> sp.dim(), sp.dim(true)
        2x3
        2x3,2nz
        >> Sparsity.dense(2, 2).dim(), Sparsity.dense(2, 2).dim(true)
        2x2
        2x2,4nz
        >> sp
        2x3,2nz
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.diag(2, 3);
        julia> ca.dim(sp), ca.dim(sp, true)
        ("2x3", "2x3,2nz")
        julia> ca.dim(ca.Sparsity.dense(2, 2)), ca.dim(ca.Sparsity.dense(2, 2), true)
        ("2x2", "2x2,4nz")
        julia> sp
        Sparsity(2x3,2nz)
        \enddoctest
        \endlang

        \sa postfix_dim, repr_el

        \identifier{2oa} */
    std::string dim(bool with_nz=false) const;

    /** \brief Dimension string to append to a name, as used when printing symbolic primitives

        Rules:
        1. Dense scalar: ""
        2. 0-by-0: "[]"
        3. Dense column vector: "[5]"
        4. Other dense pattern: "[5x10]"
        5. Otherwise: "[5x10,3nz]"

        \lang{python}
        \doctest
        >>> repr(Sparsity.dense(1, 1).postfix_dim()), repr(Sparsity(0, 0).postfix_dim())
        ("''", "'[]'")
        >>> Sparsity.dense(5).postfix_dim(), Sparsity.dense(5, 10).postfix_dim()
        ('[5]', '[5x10]')
        >>> Sparsity.diag(5, 10).postfix_dim(), Sparsity(1, 1).postfix_dim()
        ('[5x10,5nz]', '[1x1,0nz]')
        >>> Sparsity.dense(1, 5).postfix_dim(), Sparsity(0, 3).postfix_dim()
        ('[1x5]', '[0x3]')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('"%s" "%s"\n', Sparsity.dense(1, 1).postfix_dim(), Sparsity(0, 0).postfix_dim())
        "" "[]"
        >> Sparsity.dense(5).postfix_dim(), Sparsity.dense(5, 10).postfix_dim()
        [5]
        [5x10]
        >> Sparsity.diag(5, 10).postfix_dim(), Sparsity(1, 1).postfix_dim()
        [5x10,5nz]
        [1x1,0nz]
        >> Sparsity.dense(1, 5).postfix_dim(), Sparsity(0, 3).postfix_dim()
        [1x5]
        [0x3]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.postfix_dim(ca.Sparsity.dense(1, 1)), ca.postfix_dim(ca.Sparsity(0, 0))
        ("", "[]")
        julia> ca.postfix_dim(ca.Sparsity.dense(5)), ca.postfix_dim(ca.Sparsity.dense(5, 10))
        ("[5]", "[5x10]")
        julia> ca.postfix_dim(ca.Sparsity.diag(5, 10)), ca.postfix_dim(ca.Sparsity(1, 1))
        ("[5x10,5nz]", "[1x1,0nz]")
        julia> ca.postfix_dim(ca.Sparsity.dense(1, 5)), ca.postfix_dim(ca.Sparsity(0, 3))
        ("[1x5]", "[0x3]")
        \enddoctest
        \endlang

        \sa dim

        \identifier{dg} */
    std::string postfix_dim() const;

    /** \brief Describe the position of nonzero \p k as a string

        Gives "(row i, col j)", prefixed with "nonzero index k " when the pattern is not
        dense. Indices follow the global start_index option (zero-based by default).
        \p k must be in [0, nnz()); it is not checked.

        \lang{python}
        \doctest
        >>> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3])
        >>> sp.repr_el(1)
        'nonzero index 1 (row 2, col 0)'
        >>> Sparsity.dense(2, 2).repr_el(3)
        '(row 1, col 1)'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        >> sp.repr_el(1)
        nonzero index 2 (row 3, col 1)
        >> Sparsity.dense(2, 2).repr_el(3)
        (row 2, col 2)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.triplet(3, 4, [0, 2, 1, 0], [0, 0, 2, 3]);
        julia> ca.repr_el(sp, 1)
        "nonzero index 1 (row 2, col 0)"
        julia> ca.repr_el(ca.Sparsity.dense(2, 2), 3)
        "(row 1, col 1)"
        \enddoctest
        \endlang

        \sa dim, get_nz

        \identifier{2ob} */
    std::string repr_el(casadi_int k) const;

    /** \brief Print a textual picture of the pattern

        Writes one line per row, with "*" for a structural nonzero and "." for a structural
        zero, to the standard output (or the given stream). Large patterns produce
        nrow lines of ncol characters.

        \sa spy_matlab, export_code

        \identifier{dh} */
    void spy(std::ostream &stream=casadi::uout()) const;

    /** \brief Write a MATLAB/Octave script that builds the pattern as a sparse matrix and calls spy

        The file is overwritten. The script defines a sparse matrix A of ones with this
        pattern (see export_code) and ends with spy(A).

        \param mfile Path of the .m file to write

        \lang{python}
        \doctest
        >>> import os, tempfile
        >>> d = tempfile.mkdtemp()
        >>> f = os.path.join(d, "show_pattern.m")
        >>> Sparsity.diag(2, 3).spy_matlab(f)
        >>> print(open(f).read())
        % This function was automatically generated by CasADi
        A_m = 2;
        A_n = 3;
        A_j = [1, 2];
        A_i = [1, 2];
        A_v = ones(size(A_i));
        A = sparse(A_i, A_j, A_v, A_m, A_n);
        spy(A);
        >>> os.remove(f); os.rmdir(d)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = [tempname '.m'];
        >> Sparsity.diag(2, 3).spy_matlab(f);
        >> fileread(f)
        % This function was automatically generated by CasADi
        A_m = 2;
        A_n = 3;
        A_j = [1, 2];
        A_i = [1, 2];
        A_v = ones(size(A_i));
        A = sparse(A_i, A_j, A_v, A_m, A_n);
        spy(A);
        >> delete(f);
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> d = mktempdir();
        julia> f = joinpath(d, "show_pattern.m");
        julia> ca.spy_matlab(ca.Sparsity.diag(2, 3), f)
        julia> print(read(f, String))
        % This function was automatically generated by CasADi
        A_m = 2;
        A_n = 3;
        A_j = [1, 2];
        A_i = [1, 2];
        A_v = ones(size(A_i));
        A = sparse(A_i, A_j, A_v, A_m, A_n);
        spy(A);
        julia> rm(d, recursive=true)
        \enddoctest
        \endlang

        \sa spy, export_code

        \identifier{di} */
    void spy_matlab(const std::string& mfile) const;

    /** \brief Write code that reconstructs the pattern in another language

        Only "matlab" is supported: the code defines the dimensions, one-based column and
        row indices and the nonzero values of a sparse matrix. The code is written to the
        standard output or to the given stream.

        \param lang Target language, must be "matlab"
        \param options Dictionary with
          - "name": name of the generated variable (default "sp")
          - "as_matrix": also emit "<name> = sparse(...)" (default true)
          - "nonzeros": list of strings used as the nonzero values (default: ones)
          - "indent_level": number of two-space indents (default 0)
          - "inline": not supported, true raises an error

        An unknown option raises an error.

        \lang{python}
        \doctest
        >>> Sparsity.diag(2, 3).export_code("matlab")
        sp_m = 2;
        sp_n = 3;
        sp_j = [1, 2];
        sp_i = [1, 2];
        sp_v = ones(size(sp_i));
        sp = sparse(sp_i, sp_j, sp_v, sp_m, sp_n);
        >>> opts = {"name": "A", "as_matrix": False, "nonzeros": ["1.5", "2"]}
        >>> Sparsity.diag(2).export_code("matlab", opts)
        A_m = 2;
        A_n = 2;
        A_j = [1, 2];
        A_i = [1, 2];
        A_v = [1.5, 2];
        >>> try: Sparsity.diag(2).export_code("python")
        ... except Exception as e: print("unsupported language")
        unsupported language
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Sparsity.diag(2, 3).export_code('matlab');
        sp_m = 2;
        sp_n = 3;
        sp_j = [1, 2];
        sp_i = [1, 2];
        sp_v = ones(size(sp_i));
        sp = sparse(sp_i, sp_j, sp_v, sp_m, sp_n);
        >> opts = struct('name', 'A', 'as_matrix', false, 'nonzeros', {{'1.5', '2'}});
        >> Sparsity.diag(2).export_code('matlab', opts);
        A_m = 2;
        A_n = 2;
        A_j = [1, 2];
        A_i = [1, 2];
        A_v = [1.5, 2];
        >> try, Sparsity.diag(2).export_code('python'); catch e, disp('unsupported language'); end
        unsupported language
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.export_code(ca.Sparsity.diag(2, 3), "matlab")
        sp_m = 2;
        sp_n = 3;
        sp_j = [1, 2];
        sp_i = [1, 2];
        sp_v = ones(size(sp_i));
        sp = sparse(sp_i, sp_j, sp_v, sp_m, sp_n);
        julia> opts = Dict("name" => "A", "as_matrix" => false, "nonzeros" => ["1.5", "2"]);
        julia> ca.export_code(ca.Sparsity.diag(2), "matlab", opts)
        A_m = 2;
        A_n = 2;
        A_j = [1, 2];
        A_i = [1, 2];
        A_v = [1.5, 2];
        julia> try
                   ca.export_code(ca.Sparsity.diag(2), "python")
               catch e
                   println("unsupported language")
               end
        unsupported language
        \enddoctest
        \endlang

        \sa spy_matlab, compress

        \identifier{dj} */
    void export_code(const std::string& lang, std::ostream &stream=casadi::uout(),
       const Dict& options=Dict()) const;

    /// Readable name of the public class
    static std::string type_name() {return "Sparsity";}

    /** \brief Hash of the pattern (dimensions, colind and row)

        Equal patterns have equal hashes; the hash is used for the pattern cache.
        In Python, Sparsity objects are not hashable with the built-in hash(); use this
        method.

        \lang{python}
        \doctest
        >>> sp = Sparsity.diag(2, 3)
        >>> sp.hash() == Sparsity.diag(2, 3).hash(), sp.hash() == Sparsity.dense(2, 3).hash()
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sp = Sparsity.diag(2, 3);
        >> class(sp.hash())
        uint64
        >> [sp.hash() == Sparsity.diag(2, 3).hash(), sp.hash() == Sparsity.dense(2, 3).hash()]
        1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sp = ca.Sparsity.diag(2, 3);
        julia> ca.hash(sp) == ca.hash(ca.Sparsity.diag(2, 3)), ca.hash(sp) == ca.hash(ca.Sparsity.dense(2, 3))
        (true, false)
        julia> sp == ca.Sparsity.diag(2, 3), hash(sp) == hash(ca.Sparsity.diag(2, 3))
        (true, true)
        \enddoctest
        \endlang

        \sa is_equal

        \identifier{2oc} */
    std::size_t hash() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    /** \brief Pattern of a KKT system [H + I1, J'; J, I2]

        Block pattern of the KKT matrix of an NLP with Hessian \p H (square, n-by-n) and
        constraint Jacobian \p J (m-by-n), where I1 and I2 are the n-by-n and m-by-m
        diagonal patterns that are included on request. Shape (n+m)-by-(n+m).

        \param H Hessian pattern, must be square
        \param J Jacobian pattern, its number of columns must equal the size of \p H
        \param with_x_diag Add the diagonal to the H block
        \param with_lam_g_diag Make the lower right block diagonal instead of all-zero

        \lang{python}
        \doctest
        >>> DM(Sparsity.kkt(Sparsity.diag(2), Sparsity.dense(1, 2)), 1)
        DM(
        [[1, 00, 1],
         [00, 1, 1],
         [1, 1, 1]])
        >>> DM(Sparsity.kkt(Sparsity(2, 2), Sparsity.dense(1, 2), False, False), 1)
        DM(
        [[00, 00, 1],
         [00, 00, 1],
         [1, 1, 00]])
        >>> try: Sparsity.kkt(Sparsity.diag(2), Sparsity.dense(1, 3))
        ... except Exception as e: print("dimension mismatch")
        dimension mismatch
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> DM(Sparsity.kkt(Sparsity.diag(2), Sparsity.dense(1, 2)), 1)
        [[1, 00, 1],
         [00, 1, 1],
         [1, 1, 1]]
        >> DM(Sparsity.kkt(Sparsity(2, 2), Sparsity.dense(1, 2), false, false), 1)
        [[00, 00, 1],
         [00, 00, 1],
         [1, 1, 00]]
        >> J = Sparsity.dense(1, 3);
        >> try, Sparsity.kkt(Sparsity.diag(2), J); catch e, disp('dimension mismatch'); end
        dimension mismatch
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.DM(ca.Sparsity.kkt(ca.Sparsity.diag(2), ca.Sparsity.dense(1, 2)), 1)
        DM(
        [[1, 00, 1],
         [00, 1, 1],
         [1, 1, 1]])
        julia> ca.DM(ca.Sparsity.kkt(ca.Sparsity(2, 2), ca.Sparsity.dense(1, 2), false, false), 1)
        DM(
        [[00, 00, 1],
         [00, 00, 1],
         [1, 1, 00]])
        julia> try
                   ca.Sparsity.kkt(ca.Sparsity.diag(2), ca.Sparsity.dense(1, 3))
               catch e
                   println("dimension mismatch")
               end
        dimension mismatch
        \enddoctest
        \endlang

        \sa blockcat

        \identifier{dk} */
    static Sparsity kkt(const Sparsity& H, const Sparsity& J,
                        bool with_x_diag=true, bool with_lam_g_diag=true);

#ifndef SWIG
    /** \brief Assign the nonzero entries of one sparsity pattern to the nonzero

     * entries of another sparsity pattern

        \identifier{dl} */
    template<typename T>
    void set(T* data, const T* val_data, const Sparsity& val_sp) const;

    /** \brief Add the nonzero entries of one sparsity pattern to the nonzero entries

     * of another sparsity pattern

        \identifier{dm} */
    template<typename T>
    void add(T* data, const T* val_data, const Sparsity& val_sp) const;

    /** \brief Bitwise or of the nonzero entries of one sparsity pattern and the nonzero

     * entries of another sparsity pattern

        \identifier{dn} */
    template<typename T>
    void bor(T* data, const T* val_data, const Sparsity& val_sp) const;

    static std::string file_format(const std::string& filename,
      const std::string& format_hint, const std::set<std::string>& file_formats);
    static std::set<std::string> file_formats;
  private:

    /// Construct a sparsity pattern from vectors, reuse cached pattern if possible
    void assign_cached(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& colind,
                      const std::vector<casadi_int>& row, bool order_rows=false);

    /// Construct a sparsity pattern from vectors, reuse cached pattern if possible
    void assign_cached(casadi_int nrow, casadi_int ncol,
                        const casadi_int* colind, const casadi_int* row, bool order_rows=false);

#endif //SWIG
  };

  /** \brief Hash a sparsity pattern

      \identifier{ds} */
  CASADI_EXPORT std::size_t hash_sparsity(casadi_int nrow, casadi_int ncol,
                                          const std::vector<casadi_int>& colind,
                                          const std::vector<casadi_int>& row);

  CASADI_EXPORT std::size_t hash_sparsity(casadi_int nrow, casadi_int ncol,
                                          const casadi_int* colind,
                                          const casadi_int* row);

#ifndef SWIG
  // Template instantiations
  template<typename DataType>
  void Sparsity::set(DataType* data, const DataType* val_data, const Sparsity& val_sp) const {
    // Get dimensions of this
    const casadi_int sz = nnz();
    const casadi_int sz1 = size1();
    const casadi_int sz2 = size2();

    // Get dimensions of assigning matrix
    const casadi_int val_sz = val_sp.nnz();
    const casadi_int val_sz1 = val_sp.size1();
    const casadi_int val_sz2 = val_sp.size2();
    const casadi_int val_nel = val_sz1*val_sz2;

    // Check if sparsity matches
    if (val_sp==*this) {
      std::copy(val_data, val_data+sz, data);
    } else if (this->is_empty()) {
      // Quick return
      return;
    } else if (val_sp.is_empty()) {
      // Quick return
      return;
    } else if (val_nel==1) { // if scalar
      std::fill(data, data+sz, val_sz==0 ? DataType(0) : val_data[0]);
    } else if (sz2==val_sz2 && sz1==val_sz1) {
      // Matching dimensions
      // Sparsity
      const casadi_int* c = row();
      const casadi_int* rind = colind();
      const casadi_int* v_c = val_sp.row();
      const casadi_int* v_rind = val_sp.colind();

      // For all columns
      for (casadi_int i=0; i<sz2; ++i) {

        // Nonzero of the assigning matrix
        casadi_int v_el = v_rind[i];

        // First nonzero of the following column
        casadi_int v_el_end = v_rind[i+1];

        // Next row of the assigning matrix
        casadi_int v_j = v_el<v_el_end ? v_c[v_el] : sz1;

        // Assign all nonzeros
        for (casadi_int el=rind[i]; el!=rind[i+1]; ++el) {

          //  Get row
          casadi_int j=c[el];

          // Forward the assigning nonzero
          while (v_j<j) {
            v_el++;
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }

          // Assign nonzero
          if (v_j==j) {
            data[el] = val_data[v_el++];
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          } else {
            data[el] = 0;
          }
        }
      }
    } else if (sz1==val_sz2 && sz2==val_sz1 && sz2 == 1) {
      // Assign transposed (this is column)
      const casadi_int* v_cind = val_sp.colind();
      const casadi_int* r = row();
      for (casadi_int el=0; el<sz; ++el) {
        casadi_int rr=r[el];
        data[el] = v_cind[rr]==v_cind[rr+1] ? 0 : val_data[v_cind[rr]];
      }
    } else if (sz1==val_sz2 && sz2==val_sz1 && sz1 == 1) {
      // Assign transposed (this is row)
      for (casadi_int el=0; el<sz; ++el) data[el] = 0;
      const casadi_int* cind = colind();
      const casadi_int* v_r = val_sp.row();
      for (casadi_int el=0; el<val_sz; ++el) {
        casadi_int rr=v_r[el];
        if (cind[rr]!=cind[rr+1]) {
          data[cind[rr]] = val_data[el];
        }
      }
    } else {
      // Make sure that dimension matches
      casadi_error("Sparsity::set<DataType>: shape mismatch. lhs is "
                   + dim() + ", while rhs is " + val_sp.dim() + ".");
    }
  }

  template<typename DataType>
  void Sparsity::add(DataType* data, const DataType* val_data, const Sparsity& val_sp) const {
    // Get dimensions of this
    const casadi_int sz = nnz();
    const casadi_int sz1 = size1();
    const casadi_int sz2 = size2();
    const casadi_int nel = sz1*sz2;

    // Get dimensions of assigning matrix
    const casadi_int val_sz = val_sp.nnz();
    const casadi_int val_sz1 = val_sp.size1();
    const casadi_int val_sz2 = val_sp.size2();
    const casadi_int val_nel = val_sz1*val_sz2;

    // Check if sparsity matches
    if (val_sp==*this) {
      for (casadi_int k=0; k<sz; ++k) {
        data[k] += val_data[k];
      }
    } else if (this->is_empty()) {
      // Quick return
      return;
    } else if (val_sp.is_empty()) {
      // Quick return
      return;
    }  else if (val_nel==1) { // if scalar
      if (val_sz!=0) {
        for (casadi_int k=0; k<sz; ++k) {
          data[k] += val_data[0];
        }
      }
    } else {
      // Quick return if empty
      if (nel==0 && val_nel==0) return;

      // Make sure that dimension matches
      casadi_assert(sz2==val_sz2 && sz1==val_sz1,
                            "Sparsity::add<DataType>: shape mismatch. lhs is "
                            + dim() + ", while rhs is " + val_sp.dim() + ".");

      // Sparsity
      const casadi_int* c = row();
      const casadi_int* rind = colind();
      const casadi_int* v_c = val_sp.row();
      const casadi_int* v_rind = val_sp.colind();

      // For all columns
      for (casadi_int i=0; i<sz2; ++i) {

        // Nonzero of the assigning matrix
        casadi_int v_el = v_rind[i];

        // First nonzero of the following column
        casadi_int v_el_end = v_rind[i+1];

        // Next row of the assigning matrix
        casadi_int v_j = v_el<v_el_end ? v_c[v_el] : sz1;

        // Assign all nonzeros
        for (casadi_int el=rind[i]; el!=rind[i+1]; ++el) {

          //  Get row
          casadi_int j=c[el];

          // Forward the assigning nonzero
          while (v_j<j) {
            v_el++;
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }

          // Assign nonzero
          if (v_j==j) {
            data[el] += val_data[v_el++];
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }
        }
      }
    }
  }

  template<typename DataType>
  void Sparsity::bor(DataType* data, const DataType* val_data, const Sparsity& val_sp) const {
    // Get dimensions of this
    const casadi_int sz = nnz();
    const casadi_int sz1 = size1();
    const casadi_int sz2 = size2();
    const casadi_int nel = sz1*sz2;

    // Get dimensions of assigning matrix
    const casadi_int val_sz = val_sp.nnz();
    const casadi_int val_sz1 = val_sp.size1();
    const casadi_int val_sz2 = val_sp.size2();
    const casadi_int val_nel = val_sz1*val_sz2;

    // Check if sparsity matches
    if (val_sp==*this) {
      for (casadi_int k=0; k<sz; ++k) {
        data[k] |= val_data[k];
      }
    } else if (this->is_empty()) {
      // Quick return
      return;
    } else if (val_sp.is_empty()) {
      // Quick return
      return;
    }  else if (val_nel==1) { // if scalar
      if (val_sz!=0) {
        for (casadi_int k=0; k<sz; ++k) {
          data[k] |= val_data[0];
        }
      }
    } else {
      // Quick return if empty
      if (nel==0 && val_nel==0) return;

      // Make sure that dimension matches
      casadi_assert(sz2==val_sz2 && sz1==val_sz1,
                            "Sparsity::add<DataType>: shape mismatch. lhs is "
                            + dim() + ", while rhs is " + val_sp.dim() + ".");

      // Sparsity
      const casadi_int* c = row();
      const casadi_int* rind = colind();
      const casadi_int* v_c = val_sp.row();
      const casadi_int* v_rind = val_sp.colind();

      // For all columns
      for (casadi_int i=0; i<sz2; ++i) {

        // Nonzero of the assigning matrix
        casadi_int v_el = v_rind[i];

        // First nonzero of the following column
        casadi_int v_el_end = v_rind[i+1];

        // Next row of the assigning matrix
        casadi_int v_j = v_el<v_el_end ? v_c[v_el] : sz1;

        // Assign all nonzeros
        for (casadi_int el=rind[i]; el!=rind[i+1]; ++el) {

          //  Get row
          casadi_int j=c[el];

          // Forward the assigning nonzero
          while (v_j<j) {
            v_el++;
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }

          // Assign nonzero
          if (v_j==j) {
            data[el] |= val_data[v_el++];
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }
        }
      }
    }
  }

#endif //SWIG

  ///@{
  /// Readability typedefs
  typedef std::map<std::string, Sparsity> SpDict;
  ///@}

} // namespace casadi

#endif // CASADI_SPARSITY_HPP
