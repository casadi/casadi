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


#ifndef CASADI_MX_HPP
#define CASADI_MX_HPP
#include "shared_object.hpp"
#include "matrix_fwd.hpp"
#include "sx_fwd.hpp"
#include "dm.hpp"
#include "generic_matrix.hpp"
#include "generic_expression.hpp"
#include "generic_type.hpp"
#include "printable.hpp"
#include <vector>
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

namespace casadi {

  /** \brief  Forward declaration

      \identifier{px} */
  class MXNode;
  class Function;
  class SerializingStream;
  class DeserializingStream;

#ifndef SWIG
  struct ConvexifyData {
    std::vector<casadi_int> scc_offset, scc_mapping;
    Sparsity scc_sp;
    Sparsity Hrsp;
    Sparsity Hsp;
    casadi_convexify_config<double> config;
    casadi_int sz_iw;
    casadi_int sz_w;
    casadi_int *iw;
    double* w;
  };
#endif

  /** \brief MX - Matrix expression

      The MX class is used to build up trees made up from MXNodes. It is a more general
      graph representation than the scalar expression, SX, and much less efficient for small
      objects. On the other hand, the class allows much more general operations than does SX,
      in particular matrix valued operations and calls to arbitrary differentiable functions.

      The MX class is designed to have identical syntax with the Matrix<> template class,
      and uses DM (i.e. Matrix<double>) as its internal representation of the values at a node. By keeping
      the syntaxes identical, it is possible to switch from one class to the other,
      as well as inlining MX functions to SXElem functions.

      Note that an operation is always "lazy", making a matrix multiplication will create a
      matrix multiplication node, not perform the actual multiplication.

      Each node is a matrix-valued operation (symbolic primitive, constant, binary operation,
      concatenation, matrix product, call to a Function, ...), so an expression such as
      A@b+1 is a small number of nodes regardless of the matrix sizes. The graph is immutable:
      every operation creates new nodes. Symbolic primitives are created with MX.sym;
      MX(...) creates constants.

      \doctest
      >>> x = MX.sym("x", 2)
      >>> A = MX.sym("A", 2, 2)
      >>> e = A @ x + 1
      >>> e
      MX((ones(2x1)+mac(A,x,zeros(2x1))))
      >>> e.shape, e.n_dep(), e.is_symbolic()
      ((2, 1), 2, False)
      >>> f = Function("f", [A, x], [e])
      >>> f(DM([[1, 2], [3, 4]]), DM([1, 1]))
      DM([4, 8])
      \enddoctest

      \doctest{matlab}
      >> x = MX.sym('x', 2); A = MX.sym('A', 2, 2);
      >> e = A*x + 1;
      >> e
      (ones(2x1)+mac(A,x,zeros(2x1)))
      >> size(e), e.n_dep(), e.is_symbolic()
      2 1
      2
      0
      >> f = Function('f', {A, x}, {e});
      >> f(DM([1 2; 3 4]), DM([1; 1]))
      [4, 8]
      \enddoctest

      \sa SX, DM, Function

      \author Joel Andersson
      \date 2010-2011

      \identifier{py} */
  class CASADI_EXPORT MX :
    public SWIG_IF_ELSE(GenericExpressionCommon, GenericExpression<MX>),
    public SWIG_IF_ELSE(PrintableCommon, Printable<MX>),
    public GenericMatrix<MX>,
    public SharedObject {
  public:
    /** \brief Get type name

        \identifier{pz} */
    static std::string type_name() {return "MX";}

    /** \brief  Default constructor: a 0-by-0 constant

        The result is an empty constant, not a symbolic primitive (use MX.sym for that).

        \doctest
        >>> MX(), MX().shape, MX().is_constant()
        (MX(0x0), (0, 0), True)
        \enddoctest

        \doctest{matlab}
        >> MX(), size(MX()), MX().is_constant()
        0x0
        0 0
        1
        \enddoctest

        \identifier{q0} */
    MX();

    /** \brief Create a sparse matrix with all structural zeros

        The result is a constant of shape \p nrow-by-\p ncol without any nonzeros.
        It is the usual starting point for building a matrix by assignment: assigning to
        an element grows the sparsity pattern.

        \doctest
        >>> Z = MX(2, 3)
        >>> Z, Z.nnz(), Z.is_zero()
        (MX(zeros(2x3,0nz)), 0, True)
        >>> Z[1, 2] = MX.sym("x")
        >>> Z.sparsity()
        Sparsity(2x3,1nz)
        \enddoctest

        \doctest{matlab}
        >> Z = MX(2, 3);
        >> Z, Z.nnz(), Z.is_zero()
        zeros(2x3,0nz)
        0
        1
        >> Z(2, 3) = MX.sym('x');
        >> Z.sparsity()
        2x3,1nz
        \enddoctest

        \sa MX.zeros, MX.sym

        \identifier{q1} */
    MX(casadi_int nrow, casadi_int ncol);

#ifndef SWIG
    /** \brief Create a sparse matrix with all structural zeros

        \identifier{q2} */
    explicit MX(const std::pair<casadi_int, casadi_int>& rc);
#endif // SWIG

    /** \brief Create a constant with a given sparsity pattern and all nonzeros equal to one

        Same as MX::ones(sp).

        \doctest
        >>> MX(Sparsity.lower(2))
        MX(ones(2x2,3nz))
        >>> evalf(MX(Sparsity.lower(2)))
        DM(
        [[1, 00],
         [1, 1]])
        \enddoctest

        \doctest{matlab}
        >> MX(Sparsity.lower(2))
        ones(2x2,3nz)
        >> evalf(MX(Sparsity.lower(2)))
        [[1, 00],
         [1, 1]]
        \enddoctest

        \sa MX.ones

        \identifier{q3} */
    explicit MX(const Sparsity& sp);

    /** \brief Construct a matrix with a given sparsity pattern and given nonzeros

        \param sp Sparsity pattern of the result
        \param val Either a scalar, which is repeated over all nonzeros of \p sp, or a
               column vector with exactly sp.nnz() entries, which become the nonzeros of the
               result in column-major order. A matrix whose shape is a reshape of \p sp is
               reshaped. Anything else is an error.

        \p val may be symbolic; the result is then an expression, not a constant.

        \doctest
        >>> v = MX.sym("v", 3)
        >>> L = MX(Sparsity.lower(2), v)
        >>> L.sparsity(), L.is_constant()
        (Sparsity(2x2,3nz), False)
        >>> Function("f", [v], [L])(DM([1, 2, 3]))
        DM(
        [[1, 00],
         [2, 3]])
        >>> evalf(MX(Sparsity.lower(2), MX(7)))
        DM(
        [[7, 00],
         [7, 7]])
        \enddoctest

        \doctest{matlab}
        >> v = MX.sym('v', 3);
        >> L = MX(Sparsity.lower(2), v);
        >> L.sparsity(), L.is_constant()
        2x2,3nz
        0
        >> f = Function('f', {v}, {L});
        >> f(DM([1; 2; 3]))
        [[1, 00],
         [2, 3]]
        >> evalf(MX(Sparsity.lower(2), MX(7)))
        [[7, 00],
         [7, 7]]
        \enddoctest

        \sa MX.sym, sparsity_cast

        \identifier{q4} */
    MX(const Sparsity& sp, const MX& val);

    /** \brief Construct a constant whose nonzeros are read from a text file

        The file must contain exactly sp.nnz() whitespace-separated numbers, giving the
        nonzeros in column-major order; otherwise an error is raised. In generated C code
        the file is read at initialization time rather than embedded in the source.
        The node cannot be converted with to_DM.

        \param fname Path of the file, read when the constructor runs

        \doctest
        >>> try:
        ...   MX(Sparsity.lower(2), "no_such_file.txt")
        ... except Exception as e:
        ...   print("Cannot open file" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> try, MX(Sparsity.lower(2), 'no_such_file.txt'); catch e, disp('cannot open file'); end
        cannot open file
        \enddoctest

        \identifier{q5} */
    MX(const Sparsity& sp, const std::string& fname);


    /** \brief Construct a named constant that generated code exposes as a modifiable pool

        Behaves like the constant \p val during evaluation. In generated C code the
        nonzeros are placed in a named double pool so that they can be changed after
        compilation. The node cannot be converted with to_DM.

        \doctest
        >>> p = MX(DM([1, 2]), "p")
        >>> p, p.is_constant()
        (MX(constant_pool('p'): [1, 2]), True)
        >>> evalf(2*p)
        DM([2, 4])
        \enddoctest

        \doctest{matlab}
        >> p = MX(DM([1; 2]), 'p');
        >> p, p.is_constant()
        constant_pool('p'): [1, 2]
        1
        >> evalf(2*p)
        [2, 4]
        \enddoctest

        \sa CodeGenerator

        \identifier{2aa} */
    MX(const Matrix<double>& val, const std::string& name);

    /** \brief  Create a 1-by-1 constant (also used for implicit conversion from a number)

        \doctest
        >>> c = MX(3.5)
        >>> c, c.is_constant(), c.shape, float(c)
        (MX(3.5), True, (1, 1), 3.5)
        \enddoctest

        \doctest{matlab}
        >> c = MX(3.5);
        >> c, c.is_constant(), size(c), c.to_double()
        3.5
        1
        1 1
        3.5000
        \enddoctest

        \identifier{q6} */
    MX(double x);

#ifndef SWIG
    /** \brief  Create vector constant (also implicit type conversion)

        \identifier{q7} */
    MX(const std::vector<double> &x);
#endif

    /** \brief  Create a constant from a DM (also used for implicit conversion)

        The sparsity pattern of \p x is kept. Lists and numpy arrays are accepted through
        the DM conversion. SX cannot be converted to MX this way.

        \doctest
        >>> c = MX(DM(Sparsity.diag(2), [1, 2]))
        >>> c.sparsity(), c.is_constant()
        (Sparsity(2x2,2nz), True)
        >>> MX([1, 2, 3]).shape
        (3, 1)
        >>> MX(np.array([[1, 2], [3, 4]]))
        MX(
        [[1, 2],
         [3, 4]])
        \enddoctest

        \doctest{matlab}
        >> c = MX(DM(Sparsity.diag(2), [1; 2]));
        >> c.sparsity(), c.is_constant()
        2x2,2nz
        1
        >> size(MX([1; 2; 3]))
        3 1
        >> MX([1 2; 3 4])
        [[1, 2],
         [3, 4]]
        \enddoctest

        \sa to_DM, evalf

        \identifier{q8} */
    MX(const Matrix<double> &x);

/// \cond INTERNAL
    /** \brief  Destructor

        \identifier{q9} */
    ~MX();
/// \endcond

#ifndef SWIG
/// \cond INTERNAL
    /** \brief  Create from node

        \identifier{qa} */
    static MX create(MXNode* node);

    /** \brief  Create from node (multiple-outputs)

        \identifier{qb} */
    static std::vector<MX> createMultipleOutput(MXNode* node);
/// \endcond

    /** \brief Get the sparsity pattern

        Every node has a fixed sparsity pattern, known at construction time.

        \doctest
        >>> MX.sym("S", Sparsity.lower(2)).sparsity()
        Sparsity(2x2,3nz)
        >>> (MX.sym("S", Sparsity.lower(2)) + 1).sparsity()
        Sparsity(2x2)
        \enddoctest

        \doctest{matlab}
        >> MX.sym('S', Sparsity.lower(2)).sparsity()
        2x2,3nz
        >> e = MX.sym('S', Sparsity.lower(2)) + 1;
        >> e.sparsity()
        2x2
        \enddoctest

        \identifier{qc} */
    const Sparsity& sparsity() const;

    /// \cond INTERNAL
    /// Scalar type
    typedef MX ScalarType;
    /// \endcond

    /// Base class
    typedef GenericMatrix<MX> B;

    /// Expose base class functions
    using B::horzsplit;
    using B::diagsplit;
    using B::vertsplit;
    using B::mtimes;
    using B::repmat;
#endif // SWIG

    /** \brief Truth value of a 1-by-1 constant

        Nonzero constants are true, zero is false. Symbolic expressions and non-scalar
        constants raise an error: use a Function or evalf to evaluate first.

        \doctest
        >>> bool(MX(2)), bool(MX(0))
        (True, False)
        >>> try:
        ...   bool(MX.sym("x"))
        ... except Exception as e:
        ...   print("Can only determine truth value of a numeric MX" in str(e))
        True
        \enddoctest

        \identifier{2pp} */
    bool __nonzero__() const;

    /** \brief Get an owning reference to the sparsity pattern

        \identifier{qd} */
    Sparsity get_sparsity() const { return sparsity();}

    /** \brief Get the nonzeros as a list of 1-by-1 expressions

        Since an MX is not a container, each entry is an expression that selects one nonzero
        of this expression. When the expression satisfies is_valid_input, the entries are
        the scalar primitives themselves: for vertcat(x, y) with scalar symbols x and y the
        result is [x, y]. The list has nnz() entries, in column-major order.

        \sa x.nz[:], primitives

        \identifier{2bh} */
    std::vector<MX> get_nonzeros() const;

    /** \brief Erase a submatrix (leaving structural zeros in its place)

        Removes the nonzeros at the intersection of rows \p rr and columns \p cc from the
        sparsity pattern. The shape is unchanged. Modifies the expression in place by
        replacing it with a new node that references the remaining nonzeros.

        \param rr Row indices to erase; out-of-range values raise an error
        \param cc Column indices to erase
        \param ind1 True for 1-based indices (MATLAB), default 0-based

        \doctest
        >>> B = MX.sym("B", 3, 3)
        >>> B.erase([0], [1])
        >>> B.sparsity(), B.shape
        (Sparsity(3x3,8nz), (3, 3))
        \enddoctest

        \doctest{matlab}
        >> B = MX.sym('B', 3, 3);
        >> B.erase(1, 2, true);
        >> B.sparsity(), size(B)
        3x3,8nz
        3 3
        \enddoctest

        \sa enlarge, Sparsity.erase

        \identifier{qe} */
    void erase(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                bool ind1=false);

    /** \brief Erase elements (leaving structural zeros in their place)

        Same as the two-index overload, but \p rr holds linear (column-major) element
        indices.

        \doctest
        >>> B = MX.sym("B", 3, 3)
        >>> B.erase([0, 4])
        >>> B.sparsity()
        Sparsity(3x3,7nz)
        \enddoctest

        \doctest{matlab}
        >> B = MX.sym('B', 3, 3);
        >> B.erase([1, 5], true);
        >> B.sparsity()
        3x3,7nz
        \enddoctest

        \identifier{qf} */
    void erase(const std::vector<casadi_int>& rr, bool ind1=false);

    /** \brief Enlarge matrix

        Makes the matrix \p nrow-by-\p ncol by inserting empty rows and columns, keeping
        the existing nonzeros. Modifies the expression in place.

        \param nrow New number of rows
        \param ncol New number of columns
        \param rr Rows of the new matrix that receive the existing rows; one entry per
               existing row
        \param cc Columns of the new matrix that receive the existing columns; one entry
               per existing column
        \param ind1 True for 1-based indices (MATLAB), default 0-based

        \doctest
        >>> B = MX.sym("B", 2, 2)
        >>> B.enlarge(3, 3, [0, 2], [0, 2])
        >>> B.shape, B.nnz(), B.sparsity().row()
        ((3, 3), 4, [0, 2, 0, 2])
        \enddoctest

        \doctest{matlab}
        >> B = MX.sym('B', 2, 2);
        >> B.enlarge(3, 3, [1, 3], [1, 3], true);
        >> size(B), B.nnz(), B.sparsity().row()
        3 3
        4
        0 2 0 2
        \enddoctest

        \sa erase, Sparsity.enlarge

        \identifier{qg} */
    void enlarge(casadi_int nrow, casadi_int ncol,
                  const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                  bool ind1=false);

    /** \brief Unary minus

        Negating a negation returns the original node.

        \doctest
        >>> x = MX.sym("x")
        >>> -x, -(-x)
        (MX((-x)), MX(x))
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> -x, -(-x)
        (-x)
        x
        \enddoctest

        \identifier{2pq} */
    MX operator-() const;

    /// Unary plus: returns the expression itself
    MX operator+() const { return *this; }

    /** \brief Element-wise inverse, 1/x

        Not the matrix inverse: use inv or inv_node for that. Applied twice it returns
        the original node. Not available in MATLAB, where inv is the matrix inverse;
        use 1./x there.

        \doctest
        >>> A = MX.sym("A", 2, 2)
        >>> A.inv(), A.inv().inv()
        (MX((1./A)), MX(A))
        >>> evalf(MX(4).inv())
        DM(0.25)
        \enddoctest

        \sa inv, inv_node

        \identifier{2dz} */
    MX inv() const;

#ifndef SWIG
    /// \cond INTERNAL
    ///@{
    /** \brief  Access a member of the node

        \identifier{qh} */
    MXNode* operator->();

    /** \brief  Const access a member of the node

        \identifier{qi} */
    const MXNode* operator->() const;
    ///@}
    /// \endcond
#endif // SWIG

    /** \brief Get the \p ch-th dependency (child node) of this node

        Together with n_dep, op and is_op this allows walking the expression graph.
        Symbolic primitives and constants have no dependencies.

        \param ch Index in [0, n_dep()); out of range raises an error

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> e = x*y + 1
        >>> e.n_dep(), e.dep(0), e.dep(1)
        (2, MX(1), MX((x*y)))
        >>> e.dep(1).dep(0), x.n_dep()
        (MX(x), 0)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> e = x*y + 1;
        >> e.n_dep(), e.dep(0), e.dep(1)
        2
        1
        (x*y)
        >> e.dep(1).dep(0), x.n_dep()
        x
        0
        \enddoctest

        \sa n_dep, op, is_op

        \identifier{qj} */
    MX dep(casadi_int ch=0) const;

    /** \brief  Number of outputs of this node

        Greater than one only for multiple-output nodes, i.e. a call to a Function with
        several outputs. Every other node, including the OutputNode that selects one
        output of such a call, reports 1.

        \doctest
        >>> a = SX.sym("a")
        >>> f = Function("f", [a], [a + 1, 2*a])
        >>> x = MX.sym("x")
        >>> r0, r1 = f(x, )
        >>> call = r0.dep(0)
        >>> call, call.n_out(), r0.n_out()
        (MX(f(x)), 2, 1)
        \enddoctest

        \doctest{matlab}
        >> a = SX.sym('a');
        >> f = Function('f', {a}, {a + 1, 2*a});
        >> x = MX.sym('x');
        >> [r0, r1] = f(x);
        >> call = r0.dep(0);
        >> call, call.n_out(), r0.n_out()
        f(x)
        2
        1
        \enddoctest

        \sa get_output, has_output, is_output

        \identifier{qk} */
    casadi_int n_out() const;

    /** \brief  Get an output of a multiple-output node

        For a multiple-output node (a Function call), returns the OutputNode selecting
        output \p oind. For any other node only \p oind == 0 is valid and the node itself is
        returned.

        \doctest
        >>> a = SX.sym("a")
        >>> f = Function("f", [a], [a + 1, 2*a])
        >>> x = MX.sym("x")
        >>> r0, r1 = f(x, )
        >>> call = r0.dep(0)
        >>> call.get_output(1), is_equal(call.get_output(0), r0)
        (MX(f(x){1}), True)
        >>> x.get_output(0)
        MX(x)
        \enddoctest

        \doctest{matlab}
        >> a = SX.sym('a');
        >> f = Function('f', {a}, {a + 1, 2*a});
        >> x = MX.sym('x');
        >> [r0, r1] = f(x);
        >> call = r0.dep(0);
        >> call.get_output(1), is_equal(call.get_output(0), r0)
        f(x){1}
        1
        >> x.get_output(0)
        x
        \enddoctest

        \sa n_out, which_output

        \identifier{ql} */
    MX get_output(casadi_int oind) const;

    /** \brief Get the number of dependencies (child nodes) of this node

        Zero for symbolic primitives and constants, one for unary operations and
        projections, two for binary operations, the number of arguments for
        concatenations and Function calls.

        \doctest
        >>> x = MX.sym("x")
        >>> A = MX.sym("A", 2, 2)
        >>> x.n_dep(), sin(x).n_dep(), (x + 1).n_dep(), vertcat(x, x, x).n_dep()
        (0, 1, 2, 3)
        >>> (A @ A).n_dep(), A[0, 0].n_dep(), MX(1).n_dep()
        (3, 1, 0)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); A = MX.sym('A', 2, 2);
        >> e = x + 1; v = vertcat(x, x, x);
        >> x.n_dep(), sin(x).n_dep(), e.n_dep(), v.n_dep()
        0
        1
        2
        3
        >> P = A*A; a = A(1, 1);
        >> P.n_dep(), a.n_dep(), MX(1).n_dep()
        3
        1
        0
        \enddoctest

        \sa dep

        \identifier{qm} */
    casadi_int n_dep() const;

    /** \brief Get the name of a symbolic primitive

        Only defined for nodes created by MX.sym; any other node raises an error.
        Check with is_symbolic first.

        \doctest
        >>> x = MX.sym("x", 2)
        >>> x.name()
        'x'
        >>> try:
        ...   (x + 1).name()
        ... except Exception as e:
        ...   print("'name' not defined" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x', 2);
        >> x.name()
        x
        >> e = x + 1;
        >> try, e.name(); catch err, disp('error'); end
        error
        \enddoctest

        \sa is_symbolic, symvar

        \identifier{2pr} */
    std::string name() const;

    /** \brief Get the value of a 1-by-1 constant node (float conversion)

        Raises an error for symbolic nodes and for constants that are not 1-by-1.

        \doctest
        >>> float(MX(2))
        2.0
        >>> try:
        ...   float(MX.sym("x"))
        ... except Exception as e:
        ...   print("'to_double' not defined" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> MX(2).to_double()
        2
        >> try, MX.sym('x').to_double(); catch e, disp('error'); end
        error
        \enddoctest

        \sa to_DM, evalf

        \identifier{2ps} */
    explicit operator double() const;

    /** \brief Get the value of a constant node as a DM

        Only defined for constant nodes (is_constant). For a general expression without
        free symbols use evalf instead. Not defined for constants created from a file or
        with a pool name.

        \doctest
        >>> MX(DM([[1, 2], [3, 4]])).to_DM()
        DM(
        [[1, 2],
         [3, 4]])
        >>> try:
        ...   (MX.sym("x") + 1).to_DM()
        ... except Exception as e:
        ...   print("'get_DM' not defined" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> MX(DM([1 2; 3 4])).to_DM()
        [[1, 2],
         [3, 4]]
        >> e = MX.sym('x') + 1;
        >> try, e.to_DM(); catch err, disp('error'); end
        error
        \enddoctest

        \sa evalf, is_constant

        \identifier{2pt} */
    explicit operator Matrix<double>() const;

    /** \brief Check if the node is a symbolic primitive (created by MX.sym)

        False for any expression built from primitives, including concatenations of
        primitives (use is_valid_input for those). A symbol without nonzeros, e.g.
        MX.sym("x", 2, 0), is stored as an empty constant and is not symbolic.

        \doctest
        >>> x = MX.sym("x")
        >>> x.is_symbolic(), (x + 1).is_symbolic(), vertcat(x, x).is_symbolic()
        (True, False, False)
        >>> MX.sym("e", 2, 0).is_symbolic()
        False
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); e = x + 1; v = vertcat(x, x);
        >> x.is_symbolic(), e.is_symbolic(), v.is_symbolic()
        1
        0
        0
        >> MX.sym('e', 2, 0).is_symbolic()
        0
        \enddoctest

        \sa is_constant, is_valid_input, name

        \identifier{2pu} */
    bool is_symbolic() const;

    /** \brief Check if the node is a constant

        True only for constant nodes. Operations between constants are folded into a new
        constant where possible, but an expression without symbols is not constant in
        general (e.g. a Function call with constant arguments); use evalf to evaluate it.

        \doctest
        >>> MX(1).is_constant(), MX(2, 2).is_constant(), MX.sym("x").is_constant()
        (True, True, False)
        >>> (MX(1) + MX(2)).is_constant()
        True
        \enddoctest

        \doctest{matlab}
        >> MX(1).is_constant(), MX(2, 2).is_constant(), MX.sym('x').is_constant()
        1
        1
        0
        >> e = MX(1) + MX(2);
        >> e.is_constant()
        1
        \enddoctest

        \sa is_symbolic, to_DM, evalf

        \identifier{2pv} */
    bool is_constant() const;

    /** \brief Check if the node is a constant with integer-valued nonzeros

        False for non-constant nodes.

        \doctest
        >>> MX(3).is_integer(), MX(3.5).is_integer(), MX.sym("x").is_integer()
        (True, False, False)
        \enddoctest

        \doctest{matlab}
        >> MX(3).is_integer(), MX(3.5).is_integer(), MX.sym('x').is_integer()
        1
        0
        0
        \enddoctest

        \identifier{2pw} */
    bool is_integer() const;

    /** \brief Check if the node is the sum of two equal expressions, or 2*x

        True for x + x (which is stored as a doubling node) and for 2*x, and for a sum
        whose two operands compare equal up to get_max_depth.

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> (x + x).is_doubled(), (2*x).is_doubled(), (x + y).is_doubled()
        (True, True, False)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> a = x + x; b = 2*x; c = x + y;
        >> a.is_doubled(), b.is_doubled(), c.is_doubled()
        1
        1
        0
        \enddoctest

        \identifier{2e0} */
    bool is_doubled() const;

    /** \brief Check if the node is a call to a Function

        A call to a Function with several outputs is a single call node; the values you get
        back from calling the Function are OutputNodes that reference it (see is_output).
        MX Functions are inlined when called with MX arguments unless never_inline is set,
        so no call node appears in that case.

        \doctest
        >>> a = SX.sym("a")
        >>> f = Function("f", [a], [a + 1, 2*a])
        >>> x = MX.sym("x")
        >>> r0, r1 = f(x, )
        >>> r0.is_call(), r0.is_output(), r0.dep(0).is_call()
        (False, True, True)
        >>> r0.dep(0).which_function().name()
        'f'
        \enddoctest

        \doctest{matlab}
        >> a = SX.sym('a');
        >> f = Function('f', {a}, {a + 1, 2*a});
        >> x = MX.sym('x');
        >> [r0, r1] = f(x);
        >> call = r0.dep(0);
        >> r0.is_call(), r0.is_output(), call.is_call()
        0
        1
        1
        >> call.which_function().name()
        f
        \enddoctest

        \sa which_function, is_output, which_output, n_out

        \identifier{2px} */
    bool is_call() const;

    /** \brief Get the Function of a call node

        Only valid when is_call() is true; raises an error otherwise, also for the
        OutputNode of a call (use dep(0) to reach the call node).

        \doctest
        >>> a = SX.sym("a")
        >>> f = Function("f", [a], [2*a])
        >>> x = MX.sym("x")
        >>> r = f(x)
        >>> r.dep(0).which_function().name()
        'f'
        >>> try:
        ...   r.which_function()
        ... except Exception as e:
        ...   print("'which_function' not defined" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> a = SX.sym('a');
        >> f = Function('f', {a}, {2*a});
        >> x = MX.sym('x');
        >> r = f(x);
        >> call = r.dep(0);
        >> call.which_function().name()
        f
        >> try, r.which_function(); catch e, disp('error'); end
        error
        \enddoctest

        \sa is_call

        \identifier{2py} */
    Function which_function() const;

    /** \brief Check if the node selects one output of a multiple-output node

        The expressions obtained by calling a Function with several outputs are OutputNodes:
        they have one dependency, the call node, and which_output gives the index.

        \doctest
        >>> a = SX.sym("a")
        >>> f = Function("f", [a], [a + 1, 2*a])
        >>> x = MX.sym("x")
        >>> r0, r1 = f(x, )
        >>> r1, r1.is_output(), r1.which_output(), r1.dep(0)
        (MX(f(x){1}), True, 1, MX(f(x)))
        >>> x.is_output(), (x + 1).is_output()
        (False, False)
        \enddoctest

        \doctest{matlab}
        >> a = SX.sym('a');
        >> f = Function('f', {a}, {a + 1, 2*a});
        >> x = MX.sym('x');
        >> [r0, r1] = f(x);
        >> r1, r1.is_output(), r1.which_output(), r1.dep(0)
        f(x){1}
        1
        1
        f(x)
        >> e = x + 1;
        >> x.is_output(), e.is_output()
        0
        0
        \enddoctest

        \sa which_output, is_call, has_output

        \identifier{2pz} */
    bool is_output() const;

    /** \brief  Check if this is a multiple-output node

        True for the call node of a Function with several outputs. Such a node cannot be
        used directly in expressions; select an output with get_output.

        \doctest
        >>> a = SX.sym("a")
        >>> f = Function("f", [a], [a + 1, 2*a])
        >>> x = MX.sym("x")
        >>> r0, r1 = f(x, )
        >>> r0.dep(0).has_output(), r0.has_output(), x.has_output()
        (True, False, False)
        \enddoctest

        \doctest{matlab}
        >> a = SX.sym('a');
        >> f = Function('f', {a}, {a + 1, 2*a});
        >> x = MX.sym('x');
        >> [r0, r1] = f(x);
        >> call = r0.dep(0);
        >> call.has_output(), r0.has_output(), x.has_output()
        1
        0
        0
        \enddoctest

        \sa n_out, get_output, is_output

        \identifier{284} */
    bool has_output() const;

    /** \brief Get the output index of an OutputNode

        Only valid when is_output() is true; raises an error otherwise.

        \doctest
        >>> a = SX.sym("a")
        >>> f = Function("f", [a], [a + 1, 2*a])
        >>> x = MX.sym("x")
        >>> r0, r1 = f(x, )
        >>> r0.which_output(), r1.which_output()
        (0, 1)
        \enddoctest

        \doctest{matlab}
        >> a = SX.sym('a');
        >> f = Function('f', {a}, {a + 1, 2*a});
        >> x = MX.sym('x');
        >> [r0, r1] = f(x);
        >> r0.which_output(), r1.which_output()
        0
        1
        \enddoctest

        \sa is_output

        \identifier{2q0} */
    casadi_int which_output() const;

    /** \brief Check if the node is a certain operation

        \param op One of the OP_* constants (OP_ADD, OP_SIN, OP_CALL, OP_GETNONZEROS, ...)

        \doctest
        >>> x = MX.sym("x", 2)
        >>> (x + 1).is_op(OP_ADD), sin(x).is_op(OP_SIN), x.is_op(OP_PARAMETER)
        (True, True, True)
        >>> x[0].is_op(OP_GETNONZEROS), MX(1).is_op(OP_CONST)
        (True, True)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x', 2);
        >> a = x + 1; b = sin(x);
        >> a.is_op(OP_ADD), b.is_op(OP_SIN), x.is_op(OP_PARAMETER)
        1
        1
        1
        >> c = x(1);
        >> c.is_op(OP_GETNONZEROS), MX(1).is_op(OP_CONST)
        1
        1
        \enddoctest

        \sa op

        \identifier{2q1} */
    bool is_op(casadi_int op) const;

    /** \brief Check if the node is a matrix product (mtimes)

        Elementwise multiplication does not count.

        \doctest
        >>> A = MX.sym("A", 2, 2)
        >>> (A @ A).is_multiplication(), (A * A).is_multiplication()
        (True, False)
        \enddoctest

        \doctest{matlab}
        >> A = MX.sym('A', 2, 2);
        >> a = A*A; b = A.*A;
        >> a.is_multiplication(), b.is_multiplication()
        1
        0
        \enddoctest

        \identifier{2q2} */
    bool is_multiplication() const;

    /** \brief Check if a binary operation is commutative

        Unary operations count as commutative. Any other node raises an error.

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> (x + y).is_commutative(), (x - y).is_commutative()
        (True, False)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> a = x + y; b = x - y;
        >> a.is_commutative(), b.is_commutative()
        1
        0
        \enddoctest

        \sa is_binary, is_unary

        \identifier{2q3} */
    bool is_commutative() const;

    /** \brief Check if the node is a norm (norm_1, norm_2, norm_inf or norm_fro)

        \doctest
        >>> A = MX.sym("A", 2, 2)
        >>> norm_2(A).is_norm(), norm_fro(A).is_norm(), sqrt(dot(A, A)).is_norm()
        (True, True, False)
        \enddoctest

        \doctest{matlab}
        >> A = MX.sym('A', 2, 2);
        >> a = norm_2(A); b = norm_fro(A); c = sqrt(dot(A, A));
        >> a.is_norm(), b.is_norm(), c.is_norm()
        1
        1
        0
        \enddoctest

        \identifier{2q4} */
    bool is_norm() const;

    /** \brief Check if the expression can be used as a Function input

        Valid inputs are symbolic primitives, reshapes and transposes of them, and
        (nested) concatenations of valid inputs. Constants without nonzeros are also
        accepted. An expression that is not a valid input raises an error when passed as a
        Function input.

        \doctest
        >>> x = MX.sym("x", 2)
        >>> y = MX.sym("y", 3)
        >>> vertcat(x, y).is_valid_input(), vertcat(x, y).T.is_valid_input()
        (True, True)
        >>> (x + 1).is_valid_input(), MX(1).is_valid_input(), MX(2, 2).is_valid_input()
        (False, False, True)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x', 2); y = MX.sym('y', 3);
        >> v = [x; y]; vt = v';
        >> v.is_valid_input(), vt.is_valid_input()
        1
        1
        >> e = x + 1;
        >> e.is_valid_input(), MX(1).is_valid_input(), MX(2, 2).is_valid_input()
        0
        0
        1
        \enddoctest

        \sa primitives, n_primitives, is_symbolic

        \identifier{qn} */
    bool is_valid_input() const;

    /** \brief Get the number of symbolic primitives in a valid Function input

        For an expression that is not a valid input the count is 1 (the expression
        itself); check is_valid_input first.

        \doctest
        >>> x = MX.sym("x", 2)
        >>> y = MX.sym("y", 3)
        >>> vertcat(x, y).n_primitives(), x.n_primitives()
        (2, 1)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x', 2); y = MX.sym('y', 3);
        >> v = [x; y];
        >> v.n_primitives(), x.n_primitives()
        2
        1
        \enddoctest

        \sa primitives, is_valid_input

        \identifier{qo} */
    casadi_int n_primitives() const;

    /** \brief Get the symbolic primitives of a valid Function input

        Returns the symbolic primitives in the order in which they appear in the
        concatenation. For an expression that is not a valid input, the expression itself
        is returned.

        \doctest
        >>> x = MX.sym("x", 2)
        >>> y = MX.sym("y", 3)
        >>> vertcat(x, y).primitives()
        [MX(x), MX(y)]
        >>> reshape(x, 1, 2).primitives()
        [MX(x)]
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x', 2); y = MX.sym('y', 3);
        >> v = [x; y];
        >> p = v.primitives();
        >> numel(p), p{1}, p{2}
        2
        x
        y
        >> p = reshape(x, 1, 2).primitives();
        >> numel(p), p{1}
        1
        x
        \enddoctest

        \sa n_primitives, split_primitives, join_primitives, is_valid_input

        \identifier{qp} */
    std::vector<MX> primitives() const;

    /// @{
    /** \brief Split a value along the symbolic primitives of this expression

        Uses this expression (a valid Function input, see is_valid_input) as a template:
        \p x must have the same shape and is cut into one piece per primitive, each with the
        shape of that primitive. The typical use is to split a solver result back into the
        variables it was concatenated from. Works with numeric (DM), SX and MX values.
        A shape mismatch raises an error.

        \doctest
        >>> x = MX.sym("x", 2)
        >>> y = MX.sym("y", 3)
        >>> X = vertcat(x, y)
        >>> X.split_primitives(DM([1, 2, 3, 4, 5]))
        [DM([1, 2]), DM([3, 4, 5])]
        >>> X.split_primitives(SX.sym("z", 5))
        [SX([z_0, z_1]), SX([z_2, z_3, z_4])]
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x', 2); y = MX.sym('y', 3);
        >> X = [x; y];
        >> parts = X.split_primitives(DM([1; 2; 3; 4; 5]));
        >> parts{1}, parts{2}
        [1, 2]
        [3, 4, 5]
        >> parts = X.split_primitives(SX.sym('z', 5));
        >> parts{1}, parts{2}
        [z_0, z_1]
        [z_2, z_3, z_4]
        \enddoctest

        \sa join_primitives, primitives

        \identifier{qq} */
    std::vector<MX> split_primitives(const MX& x) const;
    std::vector<SX> split_primitives(const SX& x) const;
    std::vector<DM> split_primitives(const DM& x) const;
    /// @}

    /// @{
    /** \brief Join values along the symbolic primitives of this expression

        Inverse of split_primitives: \p v holds one value per primitive, with the shape of
        that primitive, and the result has the shape of this expression. A wrong number of
        pieces raises an error.

        \doctest
        >>> x = MX.sym("x", 2)
        >>> y = MX.sym("y", 3)
        >>> X = vertcat(x, y)
        >>> X.join_primitives([DM([1, 2]), DM([3, 4, 5])])
        DM([1, 2, 3, 4, 5])
        >>> try:
        ...   X.join_primitives([DM([1, 2])])
        ... except Exception as e:
        ...   print("Wrong number of primitives" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x', 2); y = MX.sym('y', 3);
        >> X = [x; y];
        >> X.join_primitives({DM([1; 2]), DM([3; 4; 5])})
        [1, 2, 3, 4, 5]
        >> try, X.join_primitives({DM([1; 2])}); catch e, disp('wrong number of primitives'); end
        wrong number of primitives
        \enddoctest

        \sa split_primitives, primitives

        \identifier{qr} */
    MX join_primitives(const std::vector<MX>& v) const;
    SX join_primitives(const std::vector<SX>& v) const;
    DM join_primitives(const std::vector<DM>& v) const;
    /// @}

    /// \cond INTERNAL
    /** \brief Detect duplicate symbolic primitives in a Function input

        If a symbolic primitive appears more than once in this expression (a valid Function
        input), returns true and warns with the name of the duplicate. Only defined for
        valid inputs; other nodes raise an error.
        Note: marks the visited primitives using MX::set_temp; call reset_input()
        afterwards.

        \doctest
        >>> x = MX.sym("x", 2)
        >>> vertcat(x, x).has_duplicates()
        True
        >>> vertcat(x, x).reset_input()
        >>> vertcat(x, MX.sym("y")).has_duplicates()
        False
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x', 2);
        >> vertcat(x, x).has_duplicates()
        ...
        1
        >> vertcat(x, x).reset_input();
        >> vertcat(x, MX.sym('y')).has_duplicates()
        0
        \enddoctest

        \sa reset_input, is_valid_input

        \identifier{qs} */
    bool has_duplicates() const;

    /** \brief Reset the marker set by has_duplicates on an input expression

        \identifier{qt} */
    void reset_input() const;
  /// \endcond

    /** \brief  Check if the node is a constant identity matrix

        Only recognizes constants created as an identity (MX.eye, DM.eye); a dense
        constant that happens to equal the identity is not detected. A 1-by-1 constant one
        counts.

        \doctest
        >>> MX.eye(2).is_eye(), MX(DM([[1, 0], [0, 1]])).is_eye(), MX(1).is_eye()
        (True, False, True)
        \enddoctest

        \doctest{matlab}
        >> MX.eye(2).is_eye(), MX(DM([1 0; 0 1])).is_eye(), MX(1).is_eye()
        1
        0
        1
        \enddoctest

        \identifier{qu} */
    bool is_eye() const;

    /** \brief  Check if the node is a constant zero (false negatives are possible)

        True for constants with all nonzeros zero and for any matrix without nonzeros,
        whether symbolic or not. Expressions that simplify to zero (x - x, 0*x) are folded
        to a constant at construction and therefore report true.

        \doctest
        >>> x = MX.sym("x")
        >>> MX(0).is_zero(), MX(2, 2).is_zero(), (x - x).is_zero(), x.is_zero()
        (True, True, True, False)
        >>> MX.sym("e", 2, 0).is_zero()
        True
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); e = x - x;
        >> MX(0).is_zero(), MX(2, 2).is_zero(), e.is_zero(), x.is_zero()
        1
        1
        1
        0
        >> MX.sym('e', 2, 0).is_zero()
        1
        \enddoctest

        \identifier{qv} */
    bool is_zero() const;

    /** \brief  Check if the node is a constant with all nonzeros one (false negatives possible)

        Structural zeros are ignored: MX(Sparsity.lower(2)) is not "one" because it has
        structural zeros that would be zero when densified.

        \doctest
        >>> MX(1).is_one(), MX.ones(2, 2).is_one(), MX(Sparsity.lower(2)).is_one()
        (True, True, False)
        \enddoctest

        \doctest{matlab}
        >> MX(1).is_one(), MX.ones(2, 2).is_one(), MX(Sparsity.lower(2)).is_one()
        1
        1
        0
        \enddoctest

        \identifier{qw} */
    bool is_one() const;

    /** \brief  Check if the node is a constant with all nonzeros -1 (false negatives possible)

        \doctest
        >>> MX(-1).is_minus_one(), (-MX(1)).is_minus_one(), MX(1).is_minus_one()
        (True, True, False)
        \enddoctest

        \doctest{matlab}
        >> m = -MX(1);
        >> MX(-1).is_minus_one(), m.is_minus_one(), MX(1).is_minus_one()
        1
        1
        0
        \enddoctest

        \identifier{qx} */
    bool is_minus_one() const;

    /** \brief  Check if the node is a constant with all nonzeros 0.5 (false negatives possible)

        \identifier{2ep} */
    bool is_half() const;

    /** \brief  Check if the node is a constant with all nonzeros equal to \p val

        False negatives are possible: only constant nodes are inspected; symbolic
        expressions always give false.

        \doctest
        >>> MX(3).is_value(3), MX(DM([3, 3])).is_value(3), MX(DM([3, 4])).is_value(3)
        (True, True, False)
        >>> MX.sym("x").is_value(3)
        False
        \enddoctest

        \doctest{matlab}
        >> MX(3).is_value(3), MX(DM([3; 3])).is_value(3)
        1
        1
        >> MX(DM([3; 4])).is_value(3)
        0
        >> MX.sym('x').is_value(3)
        0
        \enddoctest

        \identifier{2eq} */
    bool is_value(double val) const;

    /** \brief  Check if the node is a constant with all nonzeros +inf (false negatives possible)

        \doctest
        >>> MX.inf().is_inf(), MX.inf(2, 2).is_inf(), MX(1).is_inf()
        (True, True, False)
        \enddoctest

        \doctest{matlab}
        >> MX.inf().is_inf(), MX.inf(2, 2).is_inf(), MX(1).is_inf()
        1
        1
        0
        \enddoctest

        \sa is_minus_inf, MX.inf

        \identifier{2e1} */
    bool is_inf() const;

    /** \brief  Check if the node is a constant with all nonzeros -inf (false negatives possible)

        \doctest
        >>> (-MX.inf()).is_minus_inf(), MX.inf().is_minus_inf()
        (True, False)
        \enddoctest

        \doctest{matlab}
        >> m = -MX.inf();
        >> m.is_minus_inf(), MX.inf().is_minus_inf()
        1
        0
        \enddoctest

        \identifier{2e2} */
    bool is_minus_inf() const;

    /** \brief  Check if the node is a transpose operation

        The transpose of a scalar is the scalar itself and that of a vector is a reshape,
        so only transposes of genuine matrices create a transpose node.

        \doctest
        >>> A = MX.sym("A", 2, 2)
        >>> A.T.is_transpose(), A.is_transpose(), MX.sym("v", 2).T.is_transpose()
        (True, False, False)
        \enddoctest

        \doctest{matlab}
        >> A = MX.sym('A', 2, 2); At = A'; vt = MX.sym('v', 2)';
        >> At.is_transpose(), A.is_transpose(), vt.is_transpose()
        1
        0
        0
        \enddoctest

        \sa T

        \identifier{qy} */
    bool is_transpose() const;

    /** \brief Check that a constant node contains neither NaN nor Inf

        Only defined for constant nodes; a symbolic expression raises an error.

        \doctest
        >>> MX(1).is_regular(), MX.inf().is_regular(), MX(DM([1, float("nan")])).is_regular()
        (True, False, False)
        >>> try:
        ...   MX.sym("x").is_regular()
        ... except Exception as e:
        ...   print("Cannot check regularity for symbolic MX" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> MX(1).is_regular(), MX.inf().is_regular()
        1
        0
        >> MX(DM([1; nan])).is_regular()
        0
        >> try, MX.sym('x').is_regular(); catch e, disp('error'); end
        error
        \enddoctest

        \identifier{2q5} */
    bool is_regular() const;

    /** \brief Check if the node is an elementwise binary operation

        Elementwise operations with two operands (OP_ADD, OP_MUL, OP_POW, comparisons, ...).
        Matrix products, concatenations and calls are not binary nodes. Note that x + x
        and x*x are stored as unary nodes (doubling, square).

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> A = MX.sym("A", 2, 2)
        >>> (x + y).is_binary(), (x**y).is_binary(), (A @ A).is_binary(), (x + x).is_binary()
        (True, True, False, False)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y'); A = MX.sym('A', 2, 2);
        >> a = x + y; b = x^y; c = A*A; d = x + x;
        >> a.is_binary(), b.is_binary(), c.is_binary(), d.is_binary()
        1
        1
        0
        0
        \enddoctest

        \sa is_unary, op

        \identifier{2q6} */
    bool is_binary() const;

    /** \brief Check if a value is always nonnegative (false negatives are allowed)

        Recognizes nonnegative constants and operations with a nonnegative range such as
        squares, absolute values, exponentials and square roots.

        \doctest
        >>> x = MX.sym("x")
        >>> (x**2).is_nonnegative(), fabs(x).is_nonnegative(), exp(x).is_nonnegative()
        (True, True, True)
        >>> x.is_nonnegative(), MX(-3).is_nonnegative(), MX(0).is_nonnegative()
        (False, False, True)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> a = x^2; b = abs(x); c = exp(x);
        >> a.is_nonnegative(), b.is_nonnegative(), c.is_nonnegative()
        1
        1
        1
        >> x.is_nonnegative(), MX(-3).is_nonnegative(), MX(0).is_nonnegative()
        0
        0
        1
        \enddoctest

        \identifier{2e3} */
    bool is_nonnegative() const;

    /** \brief Check if the node is an elementwise unary operation

        OP_SIN, OP_NEG, OP_SQ, ... Transposes, reshapes and projections are not unary
        nodes even though they have one dependency.

        \doctest
        >>> x = MX.sym("x")
        >>> A = MX.sym("A", 2, 2)
        >>> sin(x).is_unary(), (-x).is_unary(), x.is_unary(), A.T.is_unary()
        (True, True, False, False)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); A = MX.sym('A', 2, 2);
        >> a = sin(x); b = -x; c = A';
        >> a.is_unary(), b.is_unary(), x.is_unary(), c.is_unary()
        1
        1
        0
        0
        \enddoctest

        \sa is_binary, op

        \identifier{2q7} */
    bool is_unary() const;

    /** \brief Get the operation type of the node

        Returns one of the OP_* constants. Every node, not only elementwise operations,
        has an operation code (OP_PARAMETER, OP_CONST, OP_CALL, OP_MTIMES,
        OP_GETNONZEROS, ...). Function print_operator shows the symbol of an operation.

        \doctest
        >>> x = MX.sym("x")
        >>> x.op() == OP_PARAMETER, (x + 1).op() == OP_ADD, MX(1).op() == OP_CONST
        (True, True, True)
        >>> print_operator(x + 1, ["a", "b"])
        '(a+b)'
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); e = x + 1;
        >> x.op() == OP_PARAMETER, e.op() == OP_ADD, MX(1).op() == OP_CONST
        1
        1
        1
        >> print_operator(x + 1, {'a', 'b'})
        (a+b)
        \enddoctest

        \sa is_op, print_operator

        \identifier{2q8} */
    casadi_int op() const;

    /** \brief Obtain information about the node as a dictionary

        The contents depend on the node type; most nodes return an empty dictionary.
        For example, a slicing node reports the slice it takes.

        \doctest
        >>> X = MX.sym("X", 2, 2)
        >>> X[0, 1].info()
        {'slice': {'start': 2, 'step': 1, 'stop': 3}}
        >>> X.info()
        {}
        \enddoctest

        \doctest{matlab}
        >> X = MX.sym('X', 2, 2);
        >> e = X(1, 2); s = e.info();
        >> s.slice
        start: 2
        step: 1
        stop: 3
        >> isempty(fieldnames(X.info()))
        1
        \enddoctest

        \identifier{2q9} */
    Dict info() const;

    /** \brief Serialize the expression to a stream

        From Python, serialize through StringSerializer, FileSerializer or pickle
        within a pickle context; the stream object is not constructed directly.

        \doctest
        >>> x = MX.sym("x")
        >>> s = StringSerializer()
        >>> s.pack(x + 1)
        >>> d = StringDeserializer(s.encode())
        >>> d.unpack()
        MX((1+x))
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> s = StringSerializer();
        >> s.pack(x + 1);
        >> d = StringDeserializer(s.encode());
        >> d.unpack()
        (1+x)
        \enddoctest

        \sa deserialize, StringSerializer

        \identifier{qz} */
    void serialize(SerializingStream& s) const;

    /** \brief Deserialize an expression from a stream, with type disambiguation

        Counterpart of serialize; use StringDeserializer.unpack from Python.

        \identifier{r0} */
    static MX deserialize(DeserializingStream& s);

    /// \cond INTERNAL
    /// Get the temporary variable
    casadi_int get_temp() const;

    /// Set the temporary variable
    void set_temp(casadi_int t) const;
    /// \endcond

    ///@{
    /** \brief  Create an elementwise binary or unary node from an operation code

        Low-level constructors behind the operators: MX.binary(OP_ADD, x, y) is x + y.
        Operands must have the same shape, or one of them is a scalar, or one is a
        horizontal repetition of the other (which is repmat-ed). The usual
        simplifications are applied. An invalid operation code is not checked.

        \param op One of the OP_* constants
        \param unique_x,unique_y Internal hints for in-place simplification; leave false

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> MX.binary(OP_ADD, x, y), MX.unary(OP_SIN, x)
        (MX((x+y)), MX(sin(x)))
        >>> MX.binary(OP_ADD, MX.sym("V", 2, 4), MX.sym("W", 2, 2))
        MX((V+repmat(W, 2)))
        >>> try:
        ...   MX.binary(OP_ADD, MX.sym("v", 2), MX.sym("w", 3))
        ... except Exception as e:
        ...   print("Dimension mismatch" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> MX.binary(OP_ADD, x, y), MX.unary(OP_SIN, x)
        (x+y)
        sin(x)
        >> MX.binary(OP_ADD, MX.sym('V', 2, 4), MX.sym('W', 2, 2))
        (V+repmat(W, 2))
        >> v = MX.sym('v', 2); w = MX.sym('w', 3);
        >> try, MX.binary(OP_ADD, v, w); catch e, disp('dimension mismatch'); end
        dimension mismatch
        \enddoctest

        \sa op

        \identifier{r1} */
    static MX binary(casadi_int op, const MX &x, const MX &y,
        bool unique_x=false, bool unique_y=false);
    static MX unary(casadi_int op, const MX &x, bool unique=false);
    ///@}

    ///@{
    /** \brief  Create a constant with all nonzeros +inf

        Dense \p nrow-by-\p ncol (default 1-by-1), or with the given sparsity pattern.

        \doctest
        >>> MX.inf(), MX.inf(2, 1), MX.inf(Sparsity.lower(2))
        (MX(inf), MX(inf(2x1)), MX(inf(2x2,3nz)))
        \enddoctest

        \doctest{matlab}
        >> MX.inf(), MX.inf(2, 1), MX.inf(Sparsity.lower(2))
        inf
        inf(2x1)
        inf(2x2,3nz)
        \enddoctest

        \sa is_inf, nan

        \identifier{r2} */
    static MX inf(const Sparsity& sp);
    static MX inf(casadi_int nrow=1, casadi_int ncol=1);
    static MX inf(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    ///@{
    /** \brief  Create a constant with all nonzeros NaN

        Dense \p nrow-by-\p ncol (default 1-by-1), or with the given sparsity pattern.

        \doctest
        >>> MX.nan(), MX.nan(Sparsity.lower(2))
        (MX(nan), MX(nan(2x2,3nz)))
        \enddoctest

        \doctest{matlab}
        >> MX.nan(), MX.nan(Sparsity.lower(2))
        nan
        nan(2x2,3nz)
        \enddoctest

        \sa inf

        \identifier{r3} */
    static MX nan(const Sparsity& sp);
    static MX nan(casadi_int nrow=1, casadi_int ncol=1);
    static MX nan(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    /** \brief  Identity matrix constant of size \p n-by-\p n, with a diagonal sparsity pattern

        \doctest
        >>> MX.eye(2), MX.eye(2).is_eye()
        (MX(ones(2x2,2nz)), True)
        \enddoctest

        \doctest{matlab}
        >> MX.eye(2), MX.eye(2).is_eye()
        ones(2x2,2nz)
        1
        \enddoctest

        \identifier{r4} */
    static MX eye(casadi_int n);

#ifndef SWIG
    /// Get a const pointer to the node
    MXNode* get() const;
#endif // SWIG

    /// \cond INTERNAL
    /** \brief Low-level access to get_nzref
    *
    * Intended for writing unittests

        \identifier{2er} */
    MX nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const;
    /// \endcond

    ///@{
    /** \brief Get a submatrix by linear index (what x[k] does)

        Backs single-index subscripting. Elements are addressed in column-major order over
        all numel() entries, including structural zeros. Indices may be negative (counted
        from the end); out-of-range indices raise an error. The result is a new
        GetNonzeros node referencing the selected nonzeros of this expression; a structural
        zero is returned as an empty-sparsity constant. Accepted index types: integer,
        slice, list of integers, list of booleans (mask), integer DM, Sparsity pattern of
        the same shape (returns the projection on that pattern) and MX (a symbolic index
        into the nonzeros of a dense expression, evaluated at run time with truncation and
        NaN for out-of-range values).

        \param ind1 True for 1-based indices (MATLAB), false for 0-based
        \param rr Element indices

        \doctest
        >>> X = MX.sym("X", 3, 3)
        >>> X[4], X[-1], X[[0, 4]], X[1:3].shape
        (MX(X[4]), MX(X[8]), MX(X[:8:4]), (2, 1))
        >>> try:
        ...   X[9]
        ... except Exception as e:
        ...   print("Out of bounds" in str(e))
        True
        >>> S = MX.sym("S", Sparsity.lower(2))
        >>> S[2], S[2].nnz()
        (MX(00), 0)
        >>> X[Sparsity.lower(3)]
        MX(project(X))
        >>> k = MX.sym("k")
        >>> v = MX.sym("v", 4)
        >>> v[k], Function("f", [v, k], [v[k]])(DM([10, 20, 30, 40]), 1.7)
        (MX(v[k]), DM(20))
        \enddoctest

        \doctest{matlab}
        >> X = MX.sym('X', 3, 3);
        >> a = X(5); b = X(end); c = X([1, 5]); d = X([2; 3]);
        >> a, b, c, size(d)
        X[4]
        X[8]
        X[:8:4]
        2 1
        >> try, X(10); catch e, disp('out of bounds'); end
        out of bounds
        >> S = MX.sym('S', Sparsity.lower(2));
        >> e = S(3);
        >> e, e.nnz()
        00
        0
        >> X(Sparsity.lower(3))
        project(X)
        >> k = MX.sym('k'); v = MX.sym('v', 4);
        >> f = Function('f', {v, k}, {v(k)});
        >> v(k), f(DM([10; 20; 30; 40]), 2.7)
        v[(k-1)]
        20
        \enddoctest

        \sa set, get_nz, nz

        \identifier{2qa} */
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Sparsity& sp) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const MX& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const casadi_int rr) const {
      get(m, ind1, Matrix<casadi_int>(rr));
    }
    ///@}

    ///@{
    /** \brief Get a submatrix by row and column index (what x[i, j] does)

        Backs two-index subscripting. \p rr and \p cc select rows and columns (integer,
        slice, list, integer DM); the result has shape len(rr)-by-len(cc) and the
        corresponding part of the sparsity pattern. Out-of-range indices raise an error.
        An MX index is only allowed for dense expressions and selects at run time.

        \doctest
        >>> X = MX.sym("X", 3, 3)
        >>> X[0, 1], X[1, :].shape, X[[0, 2], [0, 2]].shape
        (MX(X[3]), (1, 3), (2, 2))
        >>> S = MX.sym("S", Sparsity.lower(3))
        >>> S[:, 0].sparsity(), S[0, 2].nnz()
        (Sparsity(3x1), 0)
        >>> k = MX.sym("k")
        >>> X[k, 1]
        MX(X[(k;3:6:3)])
        \enddoctest

        \doctest{matlab}
        >> X = MX.sym('X', 3, 3);
        >> a = X(1, 2); b = X(2, :); c = X([1 3], [1 3]);
        >> a, size(b), size(c)
        X[3]
        1 3
        2 2
        >> S = MX.sym('S', Sparsity.lower(3));
        >> a = S(:, 1); b = S(1, 3);
        >> a.sparsity(), b.nnz()
        3x1
        0
        >> k = MX.sym('k');
        >> X(k, 2)
        X[((k-1);3)]
        \enddoctest

        \sa set

        \identifier{2qb} */
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const Matrix<casadi_int>& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, casadi_int cc) const {
      get(m, ind1, rr, Matrix<casadi_int>(cc));
    }
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, casadi_int rr, const Slice& cc) const {
      get(m, ind1, Matrix<casadi_int>(rr), cc);
    }
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr,
                                            const Matrix<casadi_int>& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, casadi_int rr,
                                            casadi_int cc) const {
      get(m, ind1, Matrix<casadi_int>(rr), Matrix<casadi_int>(cc));
    }
    void get(MX& SWIG_OUTPUT(m), bool ind1, const MX& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const MX& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const MX& rr, const MX& cc) const;
    ///@}

    ///@{
    /** \brief Assign to a submatrix by linear index (what x[k] = m does)

        Expressions are immutable: assignment rebinds this object to a new SetNonzeros
        node whose dependencies are the old expression and \p m. Expressions built earlier
        from the old value are unaffected, but a second name bound to the same object
        (Y0 = Y) sees the new node; copy the handle with MX(Y) to keep the old value.
        The sparsity pattern grows to include the assigned
        elements; assigning a sparse \p m first erases the targeted elements. \p m must have
        the shape of the index set, be a scalar (broadcast), or be its transpose when the
        index set is a vector; anything else raises a dimension mismatch error. A Sparsity
        index assigns to the nonzero positions of that pattern. An MX index is only
        supported for the nonzero (x.nz[k]) form.

        \doctest
        >>> Y = MX.sym("Y", 3, 3)
        >>> alias, copy, twice = Y, MX(Y), 2*Y
        >>> Y[[4]] = 5
        >>> Y, Y.n_dep(), Y.dep(0), Y.dep(1)
        (MX((Y[4] = 5)), 2, MX(Y), MX(5))
        >>> alias, copy, twice
        (MX((Y[4] = 5)), MX(Y), MX((2.*Y)))
        >>> Z = MX(3, 3)
        >>> Z[[0, 8]] = MX.sym("a", 2)
        >>> Z.sparsity()
        Sparsity(3x3,2nz)
        >>> Y = MX.sym("Y", 3, 3)
        >>> Y[Sparsity.lower(3)] = 1
        >>> Y
        MX((Y[0, 1, 2, 4, 5, 8] = ones(6x1)))
        \enddoctest

        \doctest{matlab}
        >> Y = MX.sym('Y', 3, 3);
        >> alias = Y; copy = MX(Y); twice = 2*Y;
        >> Y(5) = 5;
        >> Y, Y.n_dep(), Y.dep(0), Y.dep(1)
        (Y[4] = 5)
        2
        Y
        5
        >> alias, copy, twice
        (Y[4] = 5)
        Y
        (2.*Y)
        >> Z = MX(3, 3);
        >> Z([1, 9]) = MX.sym('a', 2);
        >> Z.sparsity()
        3x3,2nz
        >> Y = MX.sym('Y', 3, 3);
        >> Y(Sparsity.lower(3)) = 1;
        >> Y
        (Y[0, 1, 2, 4, 5, 8] = ones(6x1))
        \enddoctest

        \sa get, set_nz

        \identifier{2qc} */
    void set(const MX& m, bool ind1, const Slice& rr);
    void set(const MX& m, bool ind1, const Matrix<casadi_int>& rr);
    void set(const MX& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /** \brief Assign to a submatrix by row and column index (what x[i, j] = m does)

        Same semantics as the single-index form. \p m must be len(rr)-by-len(cc), a
        scalar, or the transpose of a vector block.

        \doctest
        >>> Y = MX.sym("Y", 3, 3)
        >>> Y[:, 0] = MX.sym("c", 1, 3)
        >>> Y
        MX((Y[:3] = c'))
        >>> Y[0:2, 0:2] = 0
        >>> Y.sparsity()
        Sparsity(3x3)
        >>> try:
        ...   Y[:, 0] = MX.sym("d", 2)
        ... except Exception as e:
        ...   print("Dimension mismatch" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> Y = MX.sym('Y', 3, 3);
        >> Y(:, 1) = MX.sym('c', 1, 3);
        >> Y
        (Y[:3] = c')
        >> Y(1:2, 1:2) = 0;
        >> Y.sparsity()
        3x3
        >> try, Y(:, 1) = MX.sym('d', 2); catch e, disp('dimension mismatch'); end
        dimension mismatch
        \enddoctest

        \sa get

        \identifier{2qd} */
    void set(const MX& m, bool ind1, const Slice& rr, const Slice& cc);
    void set(const MX& m, bool ind1, const Slice& rr, const Matrix<casadi_int>& cc);
    void set(const MX& m, bool ind1, const Matrix<casadi_int>& rr, const Slice& cc);
    void set(const MX& m, bool ind1, const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc);
    ///@}

    ///@{
    /** \brief Get a set of nonzeros (what x.nz[k] does)

        Indexes the nnz() nonzeros directly, skipping structural zeros, in column-major
        order. Negative indices count from the end; out-of-range indices raise an error.
        An MX index (or an inner/outer pair of indices) creates a parametric lookup that
        is resolved at evaluation time; such a node has no derivative with respect to the
        index. For dense expressions x[k] and x.nz[k] coincide.

        \doctest
        >>> S = MX.sym("S", Sparsity.lower(3))
        >>> S.nz[1], S.nz[:].shape, len(S.nz)
        (MX(S[1]), (6, 1), 6)
        >>> k = MX.sym("k")
        >>> S.nz[k]
        MX(S[k])
        >>> MX.sym("x", 2, 2).nz[[]].shape
        (0, 1)
        \enddoctest

        \doctest{matlab}
        >> S = MX.sym('S', Sparsity.lower(3));
        >> a = S.get_nz(true, 2); b = S.get_nz(true, Slice());
        >> a, size(b), S.nnz()
        S[1]
        6 1
        6
        >> k = MX.sym('k');
        >> S.get_nz(true, k)
        S[(k-1)]
        >> size(MX.sym('x', 2, 2).get_nz(true, []))
        0 0
        \enddoctest

        \sa set_nz, get

        \identifier{2qe} */
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const Slice& kk) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& kk) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const MX& kk) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, casadi_int kk) const {
      get_nz(m, ind1, Matrix<casadi_int>(kk));
    }
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const MX& inner, const Slice& outer) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const Slice& inner, const MX& outer) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const MX& inner, const MX& outer) const;
    ///@}

    ///@{
    /** \brief Assign to a set of nonzeros (what x.nz[k] = m does)

        Rebinds this object to a SetNonzeros node; the sparsity pattern is unchanged.
        \p m must have as many nonzeros as the index set or be a scalar. With an MX index
        the assignment is parametric and resolved at evaluation time.

        \doctest
        >>> Y = MX.sym("Y", 3, 3)
        >>> Y.nz[0] = 7
        >>> Y
        MX((Y[0] = 7))
        >>> k = MX.sym("k")
        >>> Y.nz[k] = 8
        >>> Y, Y.is_op(OP_SETNONZEROS_PARAM)
        (MX(((Y[0] = 7)[k] = 8)), True)
        \enddoctest

        \doctest{matlab}
        >> Y = MX.sym('Y', 3, 3);
        >> Y.set_nz(7, true, 1);
        >> Y
        (Y[0] = 7)
        >> k = MX.sym('k');
        >> Y.set_nz(8, true, k);
        >> Y, Y.is_op(OP_SETNONZEROS_PARAM)
        ((Y[0] = 7)[(k-1)] = 8)
        1
        \enddoctest

        \sa get_nz, set

        \identifier{2qf} */
    void set_nz(const MX& m, bool ind1, const Slice& kk);
    void set_nz(const MX& m, bool ind1, const Matrix<casadi_int>& kk);
    void set_nz(const MX& m, bool ind1, const MX& kk);
    void set_nz(const MX& m, bool ind1, casadi_int kk) { set_nz(m, ind1, Matrix<casadi_int>(kk)); }
    ///@}

    ///@{
    /** \brief Compute a dense tensor contraction in Einstein notation

        Computes C_c += A_a * B_b, where a, b and c are index labels in Einstein notation:
        a label appearing in both \p a and \p b is summed over, the remaining labels index
        the result. Tensors are passed as flat column vectors in column-major (Fortran)
        order together with their dimensions. Labels are encoded as strictly negative
        integers. For the matrix product C_ij = A_ik B_kj one writes
        a = [-1, -3], b = [-3, -2], c = [-1, -2].

        \param A,B Dense column vectors with product(dim_a), product(dim_b) entries
        \param C Optional accumulator with product(dim_c) entries, added to the result;
               when omitted the result starts from zero
        \param dim_a,dim_b,dim_c Tensor dimensions
        \param a,b,c Negative labels, one per dimension
        \return Column vector with product(dim_c) entries

        Inconsistent dimensions or labels are not validated with a user-facing message.
        Also available for SX and DM.

        \doctest
        >>> A = MX.sym("A", 6)
        >>> B = MX.sym("B", 6)
        >>> C = einstein(A, B, [2, 3], [3, 2], [2, 2], [-1, -3], [-3, -2], [-1, -2])
        >>> C.shape
        (4, 1)
        >>> f = Function("f", [A, B], [C])
        >>> f(DM(range(1, 7)), DM(range(1, 7)))
        DM([22, 28, 49, 64])
        >>> reshape(DM(range(1, 7)), 2, 3) @ reshape(DM(range(1, 7)), 3, 2)
        DM(
        [[22, 49],
         [28, 64]])
        \enddoctest

        \doctest{matlab}
        >> A = MX.sym('A', 6); B = MX.sym('B', 6);
        >> C = einstein(A, B, [2 3], [3 2], [2 2], [-1 -3], [-3 -2], [-1 -2]);
        >> size(C)
        4 1
        >> f = Function('f', {A, B}, {C});
        >> f(DM((1:6)'), DM((1:6)'))
        [22, 28, 49, 64]
        >> reshape(DM((1:6)'), 2, 3) * reshape(DM((1:6)'), 3, 2)
        [[22, 49],
         [28, 64]]
        \enddoctest

        \identifier{r5} */
    static MX einstein(const MX& A, const MX& B, const MX& C,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);

    static MX einstein(const MX& A, const MX& B,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);
    ///@}

#ifndef SWIG
    /// \cond CLUTTER
    ///@{
    /// Functions called by friend functions defined for GenericExpression
    static bool is_equal(const MX& x, const MX& y, casadi_int depth=0);
    static MX mmin(const MX &x);
    static MX mmax(const MX &x);
    ///@}

    ///@{
    /// Functions called by friend functions defined for SparsityInterface
    static MX horzcat(const std::vector<MX>& x);
    static MX diagcat(const std::vector<MX>& x);
    static MX vertcat(const std::vector<MX>& x);
    static std::vector<MX> horzsplit(const MX& x, const std::vector<casadi_int>& offset);
    static std::vector<MX> diagsplit(const MX& x, const std::vector<casadi_int>& offset1,
                                     const std::vector<casadi_int>& offset2);
    static std::vector<MX> vertsplit(const MX& x, const std::vector<casadi_int>& offset);
    static MX blockcat(const std::vector< std::vector<MX > > &v);
    static MX mtimes(const MX& x, const MX& y,
                     const std::string& blas = "reference");
    static MX mac(const MX& x, const MX& y, const MX& z,
                  const std::string& blas = "reference");
    static MX reshape(const MX& x, casadi_int nrow, casadi_int ncol);
    static MX reshape(const MX& x, const Sparsity& sp);
    static MX sparsity_cast(const MX& x, const Sparsity& sp);
    static MX kron(const MX& x, const MX& b);
    static MX kron_contract(const MX& m, const MX& x, bool inner);
    static MX repmat(const MX& x, casadi_int n, casadi_int m=1);
    ///@}

    ///@{
    /// Functions called by friend functions defined for GenericMatrix
    static MX linspace(const MX& a, const MX& b, casadi_int nsteps);
    static MX jacobian(const MX& f, const MX& x, const Dict& opts = Dict());
    static MX hessian(const MX& f, const MX& x, const Dict& opts = Dict());
    static MX hessian(const MX& f, const MX& x, MX& g, const Dict& opts = Dict());
    static std::vector<std::vector<MX> >
    forward(const std::vector<MX> &ex,
            const std::vector<MX> &arg,
            const std::vector<std::vector<MX> > &v,
            const Dict& opts = Dict());
    static std::vector<std::vector<MX> >
    reverse(const std::vector<MX> &ex,
            const std::vector<MX> &arg,
            const std::vector<std::vector<MX> > &v,
            const Dict& opts = Dict());
    static std::vector<bool> which_depends(const MX &expr, const MX &var,
        casadi_int order=1, bool tr=false);
    static std::vector<bool> which_jump(const MX &expr, const MX &var, bool tr=false);
    static Sparsity jacobian_sparsity(const MX& f, const MX& x);
    static MX substitute(const MX& ex, const MX& v, const MX& vdef);
    static std::vector<MX> substitute(const std::vector<MX> &ex,
                                         const std::vector<MX> &v,
                                         const std::vector<MX> &vdef);
    static void substitute_inplace(const std::vector<MX>& v,
                                  std::vector<MX>& vdef,
                                  std::vector<MX>& ex, bool reverse);
    static MX solve(const MX& a, const MX& b);
    static MX solve(const MX& a, const MX& b, const std::string& lsolver,
                    const Dict& dict = Dict());
    static MX inv_minor(const MX& A);
    static MX inv_node(const MX& A);
    static MX inv(const MX& A, const std::string& lsolver="qr", const Dict& dict = Dict());
    static MX pinv(const MX& A, const std::string& lsolver="qr",
               const Dict& dict = Dict());
    static MX expm_const(const MX& A, const MX& t);
    static MX expm(const MX& A);
    static casadi_int n_nodes(const MX& x);
    static std::string print_operator(const MX& x, const std::vector<std::string>& args);
    static void extract(std::vector<MX>& ex, std::vector<MX>& v,
      std::vector<MX>& vdef, const Dict& opts = Dict());
    static void shared(std::vector<MX>& ex, std::vector<MX>& v,
      std::vector<MX>& vdef, const std::string& v_prefix, const std::string& v_suffix);
    static MX if_else(const MX& cond, const MX& if_true,
                      const MX& if_false, bool short_circuit=false);
    static MX conditional(const MX& ind, const std::vector<MX> &x, const MX& x_default,
                          bool short_circuit=false);
    static bool depends_on(const MX& x, const MX& arg);
    static bool contains_all(const std::vector<MX>& v, const std::vector<MX> &n);
    static bool contains_any(const std::vector<MX>& v, const std::vector<MX> &n);
    static MX simplify(const MX& x);
    static MX transform(const MX& x, const Dict& opts = Dict());
    static MX transform(const MX& x,
        const std::vector<std::vector<GenericType> >& passes, const Dict& opts = Dict());
    static std::vector<MX> transform(const std::vector<MX>& x, const Dict& opts = Dict());
    static std::vector<MX> transform(const std::vector<MX>& x,
        const std::vector<std::vector<GenericType> >& passes, const Dict& opts = Dict());
    static MX dot(const MX& x, const MX& y);
    static MX mrdivide(const MX& a, const MX& b);
    static MX mldivide(const MX& a, const MX& b);
    static MX norm_2(const MX& x);
    static MX norm_fro(const MX& x);
    static MX norm_1(const MX& x);
    static MX norm_inf(const MX& x);
    static MX unite(const MX& A, const MX& B);
    static MX trace(const MX& x);
    static MX diag(const MX& x);
    static MX sum2(const MX& x);
    static MX sum1(const MX& x);
    static MX polyval(const MX& p, const MX& x);
    static MX det(const MX& x);
    static MX det(const MX& x, const std::string& lsolver, const Dict& opts=Dict());
    static std::vector<MX> symvar(const MX& x);
    static MX nullspace(const MX& A);
    static MX repsum(const MX& x, casadi_int n, casadi_int m=1);
    static MX densify(const MX& x, const MX& val=0);
    static MX _bilin(const MX& A, const MX& x, const MX& y);
    static MX _rank1(const MX& A, const MX& alpha, const MX& x, const MX& y);
    static MX project(const MX& x, const Sparsity& sp, bool intersect=false);
    static MX cumsum(const MX &x, casadi_int axis=-1);
    static MX _logsumexp(const MX& x);
    static std::vector<MX> cse(const std::vector<MX>& e);
    static void extract_parametric(const MX &expr, const MX& par,
            MX& expr_ret, std::vector<MX>& symbols, std::vector<MX>& parametric, const Dict& opts);
    static void separate_linear(const MX &expr,
      const MX &sym_lin, const MX &sym_const,
      MX& expr_const, MX& expr_lin, MX& expr_nonlin);
    ///@}

    ///@{
    /// Functions called by friend functions defined for this class
    static MX find(const MX& x);
    static MX low(const MX& v, const MX& p, const Dict& options = Dict());
    static MX graph_substitute(const MX& x, const std::vector<MX> &v,
                               const std::vector<MX> &vdef);
    static MX graph_substitute(const MX& x, const std::vector<MX> &v,
                               const std::vector<MX> &vdef, bool& updated);
    static std::vector<MX> graph_substitute(const std::vector<MX> &ex,
                                            const std::vector<MX> &v,
                                            const std::vector<MX> &vdef);
    static std::vector<MX> graph_substitute(const std::vector<MX> &ex,
                                            const std::vector<MX> &v,
                                            const std::vector<MX> &vdef,
                                            bool& updated);
    static MX matrix_expand(const MX& e, const std::vector<MX> &boundary,
                            const Dict& options);
    static std::vector<MX> matrix_expand(const std::vector<MX>& e,
                                         const std::vector<MX>& boundary,
                                         const Dict& options);
    static MX lift(const MX& x, const MX& x_guess);
    static DM evalf(const MX& m);
    static MX bspline(const MX& x,
            const DM& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict());
    static MX bspline(const MX& x, const MX& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict());
    static MX bspline(const MX& x, const MX& coeffs,
            const std::vector<MX>& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict());
    static MX convexify(const MX& H, const Dict& opts = Dict());
    static MX stop_diff(const MX& expr, casadi_int order);
    static MX stop_diff(const MX& expr, const MX& var, casadi_int order);
    static std::vector<MX> difference(const std::vector<MX>& a, const std::vector<MX>& b);
    ///@}
    /// \endcond

#endif // SWIG

    // Simplification with reference counting awareness
    static bool simplify_ref_count(std::vector<MX>& arg,
                                   std::vector<MX>& res,
                                   const Dict& opts = Dict());

    // Simplification with constant folding
    static bool simplify_const_folding(std::vector<MX>& arg,
                                   std::vector<MX>& res,
                                   const Dict& opts = Dict());

    // Simplification by combining like terms in linear combinations
    static bool simplify_combine_terms(std::vector<MX>& arg,
                                   std::vector<MX>& res,
                                   const Dict& opts = Dict());

    /** \brief Evaluate the B-spline basis functions at given points

        Returns the sparse matrix of basis function values such that
        bspline_dual(x, knots, degree) @ coeffs equals bspline(x, coeffs, knots, degree, 1)
        at the same points. Each row corresponds to one evaluation point and each row sums
        to one inside the knot range. Used to fit coefficients to data.

        \param x Evaluation points, flattened: n_points * n_dims numbers, point by point
        \param knots One knot vector per dimension
        \param degree Degree per dimension
        \param opts "lookup_mode": list with one of "linear", "exact", "binary" per dimension
        \return n_points-by-n_coeffs DM

        \doctest
        >>> knots = [[0, 0, 0, 0, 1, 2, 3, 3, 3, 3]]
        >>> B = MX.bspline_dual([0.5, 1.5], knots, [3])
        >>> B.shape, sum2(B)
        ((2, 6), DM([1, 1]))
        >>> B @ DM(range(6))
        DM([1.17708, 2.5])
        \enddoctest

        \doctest{matlab}
        >> knots = {[0 0 0 0 1 2 3 3 3 3]};
        >> B = MX.bspline_dual([0.5, 1.5], knots, 3);
        >> size(B), sum2(B)
        2 6
        [1, 1]
        >> B * DM((0:5)')
        [1.17708, 2.5]
        \enddoctest

        \sa bspline

        \identifier{2qg} */
    static DM bspline_dual(const std::vector<double>& x,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            const Dict& opts = Dict());

    /** \brief Low-level inlined multilinear interpolation

        Builds the interpolation as explicit MX operations (low, nonzero lookups and
        arithmetic) instead of an interpolant Function, so that the grid and values can be
        symbolic. Usually you want interpolant('name', 'linear', ...) instead.
        Extrapolates linearly outside the grid.

        \param x Grid points, one vector per dimension
        \param v Values on the grid, flattened in column-major order over the grid; its
               length must be an integer multiple m of the number of grid points, giving
               m outputs per query point
        \param xq Query coordinates, one vector per dimension, all of the same length nq
        \param opts "lookup_mode": list with one of "linear", "exact", "binary", "auto"
               (default) per dimension, see low
        \return nq-by-m matrix

        \doctest
        >>> g = MX.sym("g", 3)
        >>> v = MX.sym("v", 3)
        >>> xq = MX.sym("xq", 2)
        >>> e = MX.interpn_linear([g], v, [xq])
        >>> f = Function("f", [g, v, xq], [e])
        >>> f(DM([0, 1, 2]), DM([0, 10, 40]), DM([0.5, 3]))
        DM([5, 70])
        >>> try:
        ...   MX.interpn_linear([g], MX.sym("w", 4), [xq])
        ... except Exception as e:
        ...   print("Dimension mismatch" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> g = MX.sym('g', 3); v = MX.sym('v', 3); xq = MX.sym('xq', 2);
        >> e = MX.interpn_linear({g}, v, {xq});
        >> f = Function('f', {g, v, xq}, {e});
        >> f(DM([0; 1; 2]), DM([0; 10; 40]), DM([0.5; 3]))
        [5, 70]
        >> w = MX.sym('w', 4);
        >> try, MX.interpn_linear({g}, w, {xq}); catch e, disp('dimension mismatch'); end
        dimension mismatch
        \enddoctest

        \sa interpolant, low

        \identifier{r6} */
    static MX interpn_linear(const std::vector<MX>& x, const MX& v, const std::vector<MX>& xq,
      const Dict& opts=Dict());

    /** \brief Return the expression with the side effect of printing its value with a tag

        Returns an expression equal to this one; when it is evaluated, every nonzero is
        printed to standard output as "|> b : value", where \p b is a scalar expression
        evaluated alongside it (an identifier for the print). Useful to trace values
        inside a Function. The printing is not part of the returned value and not seen by
        the doctest runner's output capture.

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> e = x.printme(y)
        >>> e, e.is_op(OP_PRINTME)
        (MX(printme(x,y)), True)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> e = x.printme(y);
        >> e, e.is_op(OP_PRINTME)
        printme(x,y)
        1
        \enddoctest

        \sa monitor, attachAssert

        \identifier{2qh} */
    MX printme(const MX& b) const;

#if !defined(SWIG) || defined(DOXYGEN)
/**
\addtogroup expression_tools
@{
*/
    /** \brief Find the index of the first nonzero entry of a vector

        Returns a 1-by-1 expression that evaluates to the index of the first entry of \p x
        whose value is not zero, or to the length of \p x when all entries are zero.
        Structural zeros are skipped. The result is piecewise constant, so its derivative
        is zero. \p x must be a column or row vector; matrices raise an error.

        \doctest
        >>> v = MX.sym("v", 3)
        >>> f = Function("f", [v], [find(v)])
        >>> find(v).shape, f(DM([0, 3, 1])), f(DM([0, 0, 0]))
        ((1, 1), DM(1), DM(3))
        >>> evalf(find(MX(DM(Sparsity.unit(3, 2), [1]))))
        DM(2)
        >>> try:
        ...   find(MX.sym("A", 2, 2))
        ... except Exception as e:
        ...   print("Argument must be vector" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> v = MX.sym('v', 3);
        >> f = Function('f', {v}, {find(v)});
        >> size(find(v)), f(DM([0; 3; 1])), f(DM([0; 0; 0]))
        1 1
        1
        3
        >> evalf(find(MX(DM(Sparsity.unit(3, 2), 1))))
        2
        >> try, find(MX.sym('A', 2, 2)); catch e, disp('argument must be vector'); end
        argument must be vector
        \enddoctest

        \sa low

        \identifier{r7} */
    inline friend MX find(const MX& x) {
      return MX::find(x);
    }

    /** \brief Find the interval of a sorted grid that contains each query value

        For each entry p_i of \p p, returns the index k such that v[k] <= p_i < v[k+1]
        for the ascending grid \p v, clipped to [0, len(v)-2] so that values outside the
        grid map to the first or last interval. This is the lookup step of linear
        interpolation. The result has the shape of \p p and is piecewise constant, so its
        derivative is zero.

        \param v Dense ascending grid, vector
        \param p Query values, dense vector
        \param options "lookup_mode": "linear" (scan), "binary" (bisection), "exact"
               (assumes an equidistant grid and computes the index arithmetically) or
               "auto" (default: "binary" for grids longer than 100 entries, else "linear").
               Any other key or mode raises an error.

        \doctest
        >>> g = MX.sym("g", 4)
        >>> p = MX.sym("p", 5)
        >>> f = Function("f", [g, p], [low(g, p)])
        >>> f(DM([0, 1, 2, 3]), DM([-1, 0, 1.5, 3, 4]))
        DM([0, 0, 1, 2, 2])
        >>> jacobian(low(g, p), p).nnz()
        0
        >>> try:
        ...   low(g, p, {"lookup_mode": "foo"})
        ... except Exception as e:
        ...   print("Invalid lookup mode" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> g = MX.sym('g', 4); p = MX.sym('p', 5);
        >> f = Function('f', {g, p}, {low(g, p)});
        >> f(DM([0; 1; 2; 3]), DM([-1; 0; 1.5; 3; 4]))
        [0, 0, 1, 2, 2]
        >> J = jacobian(low(g, p), p);
        >> J.nnz()
        0
        >> try, low(g, p, struct('lookup_mode', 'foo')); catch e, disp('invalid lookup mode'); end
        invalid lookup mode
        \enddoctest

        \sa find, interpn_linear, interpolant

        \identifier{r8} */
    inline friend MX low(const MX& v, const MX& p, const Dict& options=Dict()) {
      return MX::low(v, p, options);
    }

    /** \brief Substitute nodes in the expression graph, by node identity

        Replaces every occurrence of the node v[i] in \p ex by vdef[i]. Unlike substitute,
        which only accepts symbolic primitives for \p v, graph_substitute matches any
        node, so intermediate results and Function call results can be replaced. Matching
        is by identity: an expression that was rebuilt (x + y created twice) is a different
        node and is not matched. Nodes that do not depend on any substituted node are
        reused, not copied. Shapes of v[i] and vdef[i] must agree.

        \param ex Expression to substitute in
        \param v Nodes to replace
        \param vdef Replacements, same length as \p v
        \param updated Set to true if any substitution took place (C++ only)

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> s = x + y
        >>> e = sin(s) + s
        >>> graph_substitute(e, [s], [MX(3)])
        MX(3.14112)
        >>> graph_substitute(e, [x + y], [MX(3)])
        MX(@1=(x+y), (sin(@1)+@1))
        >>> try:
        ...   substitute(e, s, MX(3))
        ... except Exception as err:
        ...   print("not symbolic" in str(err))
        True
        >>> try:
        ...   graph_substitute(e, [s], [MX.sym("z", 2)])
        ... except Exception as err:
        ...   print("Inconsistent shapes" in str(err))
        True
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> s = x + y;
        >> e = sin(s) + s;
        >> graph_substitute(e, {s}, {MX(3)})
        3.14112
        >> graph_substitute(e, {x + y}, {MX(3)})
        @1=(x+y), (sin(@1)+@1)
        >> try, substitute(e, s, MX(3)); catch err, disp('not symbolic'); end
        not symbolic
        >> z = MX.sym('z', 2);
        >> try, graph_substitute(e, {s}, {z}); catch err, disp('inconsistent shapes'); end
        inconsistent shapes
        \enddoctest

        \sa substitute

        \identifier{r9} */
    inline friend MX graph_substitute(const MX& ex, const std::vector<MX> &v,
                                      const std::vector<MX> &vdef) {
      return MX::graph_substitute(ex, v, vdef);
    }

    inline friend MX graph_substitute(const MX& ex, const std::vector<MX> &v,
                                      const std::vector<MX> &vdef, bool& updated) {
      return MX::graph_substitute(ex, v, vdef, updated);
    }

    /** \brief Substitute nodes in several expressions at once, by node identity

        Same as the single-expression overload, applied to a list of expressions sharing
        one traversal.

        \doctest
        >>> x = MX.sym("x")
        >>> s = x + 1
        >>> graph_substitute([sin(s), 2*s], [s], [MX(3)])
        [MX(0.14112), MX(6)]
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> s = x + 1;
        >> r = graph_substitute({sin(s), 2*s}, {s}, {MX(3)});
        >> r{1}, r{2}
        0.14112
        6
        \enddoctest

        \identifier{ra} */
    inline friend std::vector<MX>
      graph_substitute(const std::vector<MX> &ex,
                       const std::vector<MX> &v,
                       const std::vector<MX> &vdef) {
      return MX::graph_substitute(ex, v, vdef);
    }

    inline friend std::vector<MX>
      graph_substitute(const std::vector<MX> &ex,
                       const std::vector<MX> &v,
                       const std::vector<MX> &vdef,
                       bool& updated) {
      return MX::graph_substitute(ex, v, vdef, updated);
    }

    /** \brief Expand an MX expression into a call to an SX Function

        Converts the whole matrix graph of \p e to scalar SX operations (including inlined
        Function calls, as Function.expand does) and returns a call to the resulting
        SXFunction, with the free symbols of \p e as arguments. Expressions listed in
        \p boundary are not expanded: they stay as MX inputs of the call.

        \param e Expression to expand
        \param boundary Sub-expressions at which the expansion stops
        \param options Options for the SXFunction, forwarded to Function.expand
        \return An OutputNode of the SXFunction call

        \doctest
        >>> x = MX.sym("x")
        >>> g = Function("g", [x], [x**2 + 1])
        >>> sx = sin(x)
        >>> e = g(x) + sx
        >>> r = matrix_expand(e)
        >>> r, r.is_output(), r.dep(0).which_function()
        (MX(expand_tmp_matrix_expand(x){0}), True, Function(expand_tmp_matrix_expand:(i0)->(o0)
        SXFunction))
        >>> Function("f", [x], [r])(2)
        DM(5.9093)
        >>> matrix_expand(e, [sx]).dep(0).which_function()
        Function(expand_tmp_matrix_expand:(i0,i1)->(o0) SXFunction)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> g = Function('g', {x}, {x^2 + 1});
        >> sx = sin(x);
        >> e = g(x) + sx;
        >> r = matrix_expand(e);
        >> r, r.is_output(), r.dep(0).which_function()
        expand_tmp_matrix_expand(x){0}
        1
        expand_tmp_matrix_expand:(i0)->(o0) SXFunction
        >> f = Function('f', {x}, {r});
        >> f(2)
        5.9093
        >> r = matrix_expand(e, {sx});
        >> r.dep(0).which_function()
        expand_tmp_matrix_expand:(i0,i1)->(o0) SXFunction
        \enddoctest

        \sa Function.expand

        \identifier{rb} */
    inline friend MX
      matrix_expand(const MX& e, const std::vector<MX> &boundary = std::vector<MX>(),
        const Dict& options = Dict()) {
      return MX::matrix_expand(e, boundary, options);
    }

    /** \brief Expand several MX expressions into one call to an SX Function

        Same as the single-expression overload; the result is one output per expression.

        \identifier{rc} */
    inline friend std::vector<MX>
      matrix_expand(const std::vector<MX>& e,
                    const std::vector<MX> &boundary = std::vector<MX>(),
                    const Dict& options = Dict()) {
      return MX::matrix_expand(e, boundary, options);
    }


    /** \brief Evaluate a tensor-product B-spline with fixed coefficients

        Returns an \p m-by-1 expression with the spline values at \p x. The spline is
        defined by one knot vector and one degree per dimension; the number of basis
        functions per dimension is len(knots[i]) - degree[i] - 1 and the coefficient
        vector must hold m times the product of these counts, in column-major order.
        Outside the knot range the spline is zero. Derivatives are available.

        \param x Column vector with one entry per dimension
        \param coeffs Coefficient vector; its length is not validated
        \param knots One knot vector per dimension, numeric
        \param degree Degree per dimension; must have as many entries as \p knots
        \param m Number of outputs per evaluation
        \param opts "lookup_mode" (see low) and "inline" (bool, default false: build an
               explicit expression instead of a BSpline node). Unknown keys are ignored.

        \doctest
        >>> x = MX.sym("x")
        >>> knots = [[0, 0, 0, 0, 1, 2, 3, 3, 3, 3]]
        >>> c = DM(range(6))
        >>> y = bspline(x, c, knots, [3], 1)
        >>> f = Function("f", [x], [y])
        >>> y, f(1.5), f(4)
        (MX(BSpline(x)), DM(2.5), DM(0))
        >>> Function("df", [x], [jacobian(y, x)])(1.5)
        DM(1.125)
        >>> try:
        ...   bspline(MX.sym("v", 2), c, knots, [3], 1)
        ... except Exception as e:
        ...   print("must match" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> knots = {[0 0 0 0 1 2 3 3 3 3]};
        >> c = DM((0:5)');
        >> y = bspline(x, c, knots, 3, 1);
        >> f = Function('f', {x}, {y});
        >> y, f(1.5), f(4)
        BSpline(x)
        2.5
        0
        >> df = Function('df', {x}, {jacobian(y, x)});
        >> df(1.5)
        1.125
        >> try, bspline(MX.sym('v', 2), c, knots, 3, 1); catch e, disp('must match'); end
        must match
        \enddoctest

        \sa bspline_dual, interpolant

        \identifier{2qi} */
    inline friend MX bspline(const MX& x,
            const DM& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict()) {
      return MX::bspline(x, coeffs, knots, degree, m, opts);
    }

    /** \brief Evaluate a B-spline whose coefficients are an expression

        Same as the fixed-coefficient overload, but \p coeffs may be symbolic, so the
        spline is differentiable with respect to its coefficients as well.

        \doctest
        >>> x = MX.sym("x")
        >>> C = MX.sym("C", 6)
        >>> knots = [[0, 0, 0, 0, 1, 2, 3, 3, 3, 3]]
        >>> y = bspline(x, C, knots, [3], 1)
        >>> y
        MX(BSplineParametric(x, C))
        >>> Function("f", [x, C], [y])(1.5, DM(range(6)))
        DM(2.5)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> C = MX.sym('C', 6);
        >> knots = {[0 0 0 0 1 2 3 3 3 3]};
        >> y = bspline(x, C, knots, 3, 1);
        >> y
        BSplineParametric(x, C)
        >> f = Function('f', {x, C}, {y});
        >> f(1.5, DM((0:5)'))
        2.5
        \enddoctest

        \identifier{2qj} */
    inline friend MX bspline(const MX& x, const MX& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict()) {
      return MX::bspline(x, coeffs, knots, degree, m, opts);
    }

    /** \brief Evaluate a B-spline whose knots and coefficients are expressions

        With symbolic knots the evaluation is always inlined into explicit MX operations.

        \doctest
        >>> x = MX.sym("x")
        >>> C = MX.sym("C", 6)
        >>> K = MX.sym("K", 10)
        >>> y = bspline(x, C, [K], [3], 1)
        >>> f = Function("f", [x, C, K], [y])
        >>> f(1.5, DM(range(6)), DM([0, 0, 0, 0, 1, 2, 3, 3, 3, 3]))
        DM(2.5)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> C = MX.sym('C', 6);
        >> K = MX.sym('K', 10);
        >> y = bspline(x, C, {K}, 3, 1);
        >> f = Function('f', {x, C, K}, {y});
        >> f(1.5, DM((0:5)'), DM([0; 0; 0; 0; 1; 2; 3; 3; 3; 3]))
        2.5
        \enddoctest

        \identifier{2qk} */
    inline friend MX bspline(const MX& x, const MX& coeffs,
            const std::vector<MX>& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict()) {
      return MX::bspline(x, coeffs, knots, degree, m, opts);
    }

    /// Evaluate the B-spline basis functions at given points, see MX::bspline_dual
    inline friend DM bspline_dual(const std::vector<double>& x,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            const Dict& opts = Dict()) {
      return MX::bspline_dual(x, knots, degree, opts);
    }

    /** \brief Make a symmetric matrix positive definite at evaluation time

        Returns an expression with the sparsity pattern of \p H that, when evaluated,
        modifies the value of \p H so that it is positive definite. Intended for Hessians
        in SQP-type methods. \p H must be square and symmetric or triangular (a triangular
        pattern is interpreted as one half of a symmetric matrix); other patterns raise an
        error. The node has no derivatives.

        \param H Square expression with a symmetric or triangular sparsity pattern
        \param opts Options:
               "strategy": "eigen-clip" (default, eigenvalues below the margin are set to
               the margin), "eigen-reflect" (their absolute value is taken) or
               "regularize" (a multiple of the identity is added; symmetric patterns only);
               "margin": smallest allowed eigenvalue, default 1e-7, must be nonnegative;
               "max_iter_eig": iteration limit of the eigenvalue computation, default 200;
               "verbose": print the block structure found.
               Unknown keys raise an error.

        \doctest
        >>> H = MX.sym("H", 2, 2)
        >>> f = Function("f", [H], [convexify(H)])
        >>> f(DM([[1, 0], [0, -2]]))
        DM(
        [[1, 0],
         [0, 1e-07]])
        >>> f = Function("f", [H], [convexify(H, {"strategy": "eigen-reflect"})])
        >>> f(DM([[1, 0], [0, -2]]))
        DM(
        [[1, 0],
         [0, 2]])
        >>> f = Function("f", [H], [convexify(H, {"strategy": "regularize", "margin": 0.5})])
        >>> f(DM([[1, 0], [0, -2]]))
        DM(
        [[3.5, 0],
         [0, 0.5]])
        >>> try:
        ...   convexify(H, {"strategy": "foo"})
        ... except Exception as e:
        ...   print("Invalid convexify strategy" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> H = MX.sym('H', 2, 2);
        >> f = Function('f', {H}, {convexify(H)});
        >> f(DM([1 0; 0 -2]))
        [[1, 0],
         [0, 1e-07]]
        >> f = Function('f', {H}, {convexify(H, struct('strategy', 'eigen-reflect'))});
        >> f(DM([1 0; 0 -2]))
        [[1, 0],
         [0, 2]]
        >> f = Function('f', {H}, {convexify(H, struct('strategy', 'regularize', 'margin', 0.5))});
        >> f(DM([1 0; 0 -2]))
        [[3.5, 0],
         [0, 0.5]]
        >> try, convexify(H, struct('strategy', 'foo')); catch e, disp('invalid strategy'); end
        invalid strategy
        \enddoctest

        \sa nlpsol

        \identifier{2ql} */
    inline friend MX convexify(const MX& H,
            const Dict& opts = Dict()) {
      return MX::convexify(H, opts);
    }

    /** \brief Lift the expression

        Experimental feature. Marks \p x as a lifted intermediate variable with initial
        guess \p x_guess, for use by lifted Newton-type methods. Evaluates to \p x.
        Both arguments must have the same sparsity pattern.

        \identifier{rd} */
    inline friend MX lift(const MX& x, const MX& x_guess) {
      return MX::lift(x, x_guess);
    }

    /** \brief Symbolic matrix inverse node

        Creates a single node representing the inverse of the square matrix \p x. The node
        defines its derivatives (via -inv(x) dx inv(x)) but cannot be evaluated
        numerically or expanded to SX; it is meant as a placeholder in symbolic
        manipulations. For an inverse that can be evaluated use inv or solve.

        \doctest
        >>> A = MX.sym("A", 2, 2)
        >>> inv_node(A), inv_node(A).is_op(OP_INVERSE), jacobian(inv_node(A), A).shape
        (MX(inv(A)), True, (4, 4))
        >>> try:
        ...   inv_node(MX.sym("B", 2, 3))
        ... except Exception as e:
        ...   print("must be square" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> A = MX.sym('A', 2, 2);
        >> e = inv_node(A); J = jacobian(e, A);
        >> e, e.is_op(OP_INVERSE), size(J)
        inv(A)
        1
        4 4
        >> try, inv_node(MX.sym('B', 2, 3)); catch e, disp('must be square'); end
        must be square
        \enddoctest

        \sa inv, solve

        \identifier{re} */
    inline friend MX inv_node(const MX& x) {
      return MX::inv_node(x);
    }

    /** \brief Evaluate an expression without free symbols numerically

        Builds a Function without inputs around \p expr and evaluates it. Function
        calls, solver calls and sparse constants are handled. An error is raised when the
        expression depends on a symbolic primitive.

        \doctest
        >>> evalf(MX(DM([[1, 2], [3, 4]])) @ MX([1, 1]))
        DM([3, 7])
        >>> evalf(MX(Sparsity.lower(2)) * 2)
        DM(
        [[2, 00],
         [2, 2]])
        >>> try:
        ...   evalf(MX.sym("x") + 1)
        ... except Exception as e:
        ...   print("are free" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> evalf(MX(DM([1 2; 3 4])) * MX([1; 1]))
        [3, 7]
        >> evalf(MX(Sparsity.lower(2)) * 2)
        [[2, 00],
         [2, 2]]
        >> try, evalf(MX.sym('x') + 1); catch e, disp('free variables'); end
        free variables
        \enddoctest

        \sa to_DM, Function

        \identifier{rf} */
    inline friend DM evalf(const MX& expr) {
      return MX::evalf(expr);
    }

    /** \brief Stop derivatives of an expression with respect to all its symbolic variables

        Wraps \p expr in a non-inlined Function call whose derivatives of the given order
        are declared zero. The value is unchanged. With \p order 1 the expression is
        treated as a constant by all differentiation; with \p order 2 first derivatives
        are propagated normally but second derivatives are zero. Other orders raise an
        error.

        \doctest
        >>> x = MX.sym("x")
        >>> e = x * stop_diff(x**2, 1)
        >>> Function("f", [x], [e, jacobian(e, x)])(3)
        (DM(27), DM(9))
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> e = x * stop_diff(x^2, 1);
        >> f = Function('f', {x}, {e, jacobian(e, x)});
        >> [v, J] = f(3);
        >> v, J
        27
        9
        \enddoctest

        \sa no_grad, no_hess

        \identifier{25l} */
    inline friend MX stop_diff(const MX& expr, casadi_int order) {
      return MX::stop_diff(expr, order);
    }

    /** \brief Treat an expression as constant for differentiation (stop first derivatives)

        Same as stop_diff(expr, 1): the value is unchanged but all derivatives with
        respect to the symbols in \p expr are zero.

        \doctest
        >>> x = MX.sym("x")
        >>> y = no_grad(x**2)
        >>> Function("f", [x], [y])(3), jacobian(y, x).nnz()
        (DM(9), 0)
        >>> Function("f", [x], [jacobian(x * y, x)])(3)
        DM(9)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> y = no_grad(x^2); J = jacobian(y, x);
        >> f = Function('f', {x}, {y});
        >> f(3), J.nnz()
        9
        0
        >> f = Function('f', {x}, {jacobian(x * y, x)});
        >> f(3)
        9
        \enddoctest

        \sa stop_diff, no_hess

        \identifier{25m} */
    inline friend MX no_grad(const MX& expr) {
      return MX::stop_diff(expr, 1);
    }

    /** \brief Stop second derivatives of an expression wrt all its symbolic variables

        Same as stop_diff(expr, 2): first derivatives are exact, the Hessian contribution
        is zero. Useful for Gauss-Newton-type Hessians.

        \doctest
        >>> x = MX.sym("x")
        >>> y = no_hess(x**3)
        >>> Function("f", [x], [jacobian(y, x)])(2), hessian(y, x)[0].nnz()
        (DM(12), 0)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> y = no_hess(x^3);
        >> f = Function('f', {x}, {jacobian(y, x)});
        >> [H, g] = hessian(y, x);
        >> f(2), H.nnz()
        12
        0
        \enddoctest

        \sa stop_diff, no_grad

        \identifier{25n} */
    inline friend MX no_hess(const MX& expr) {
      return MX::stop_diff(expr, 2);
    }


    /** \brief Stop derivatives of an expression wrt a selected set of symbolic variables

        Like stop_diff(expr, order), but only the derivatives with respect to the symbols
        occurring in \p var are stopped; derivatives with respect to the other symbols are
        propagated. Emits a warning that this variant is not well tested.

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> e = stop_diff(x * y, x, 1)
        >>> jacobian(e, x).nnz(), Function("f", [x, y], [jacobian(e, y)])(2, 3)
        (0, DM(2))
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> e = stop_diff(x * y, x, 1);
        ...
        >> J = jacobian(e, x);
        >> f = Function('f', {x, y}, {jacobian(e, y)});
        >> J.nnz(), f(2, 3)
        0
        2
        \enddoctest

        \identifier{25o} */
    inline friend MX stop_diff(const MX& expr, const MX& var, casadi_int order) {
      return MX::stop_diff(expr, var, order);
    }

    /** \brief Return all elements of \p a that do not occur in \p b, preserving order

        Comparison is by node identity, not by value: two separately built expressions
        with the same formula are different nodes. Duplicates in \p a are kept.

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> difference([x, y, x], [y]), difference([x + 1], [x + 1])
        ([MX(x), MX(x)], [MX((1+x))])
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> d = difference({x, y, x}, {y});
        >> numel(d), d{1}, d{2}
        2
        x
        x
        >> d = difference({x + 1}, {x + 1});
        >> numel(d), d{1}
        1
        (1+x)
        \enddoctest

        \sa symvar

        \identifier{2qm} */
    inline friend std::vector<MX> difference(const std::vector<MX>& a, const std::vector<MX>& b) {
      return MX::difference(a, b);
    }

/** @} */
#endif // SWIG

    /** \brief Return the expression with a run-time assertion attached

        The result is equal to this expression, but when it is evaluated the scalar
        expression \p y is evaluated too and a runtime error with \p fail_message is
        raised unless \p y equals 1. Useful to check preconditions inside a Function.
        A non-scalar \p y raises an error at construction.

        \doctest
        >>> x = MX.sym("x")
        >>> e = x.attachAssert(x > 0, "x must be positive")
        >>> e
        MX(assertion(x, (0<x)))
        >>> f = Function("f", [x], [e])
        >>> f(1)
        DM(1)
        >>> try:
        ...   f(-1)
        ... except Exception as e:
        ...   print("Assertion error: x must be positive" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> e = x.attachAssert(x > 0, 'x must be positive');
        >> e
        assertion(x, (0<x))
        >> f = Function('f', {x}, {e});
        >> f(1)
        1
        >> try, f(-1); catch err, disp('assertion error'); end
        assertion error
        \enddoctest

        \sa monitor, printme

        \identifier{rg} */
    MX attachAssert(const MX& y, const std::string& fail_message="") const;

    /** \brief Return the expression with the side effect of printing its value

        When the returned expression is evaluated, \p comment followed by the nonzeros is
        printed to standard output. The value is unchanged. A matrix without nonzeros is
        returned as is.

        \doctest
        >>> x = MX.sym("x")
        >>> e = x.monitor("value of x")
        >>> e, e.is_op(OP_MONITOR)
        (MX(monitor(x, value of x)), True)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> e = x.monitor('value of x');
        >> e, e.is_op(OP_MONITOR)
        monitor(x, value of x)
        1
        \enddoctest

        \sa printme, dump, attachAssert

        \identifier{rh} */
    MX monitor(const std::string& comment) const;

    /** \brief Return the expression with the side effect of writing its value to a file

        Each evaluation writes the value to a new file named after \p base_filename and a
        running counter. The value is unchanged. A matrix without nonzeros is returned as
        is.

        \param base_filename Prefix of the file names
        \param opts "dir" (directory, default "."), "format" (file format, default
               "mtx", Matrix Market), "verbose" (print the file name on each dump, default
               false). Unknown keys raise an error.

        \sa monitor, reset_dump_count

        \identifier{2f9} */
    MX dump(const std::string& base_filename, const Dict& opts=Dict()) const;

    /// Reset the file counter of a dump node; raises an error on any other node
    void reset_dump_count();

    /** \brief Transpose the matrix

        In Python use the property x.T. Scalars are returned as is, vectors become a
        reshape node and matrices a transpose node. Transposing twice returns the
        original node.

        \doctest
        >>> A = MX.sym("A", 2, 3)
        >>> A.T, A.T.shape, A.T.T
        (MX(A'), (3, 2), MX(A))
        \enddoctest

        \doctest{matlab}
        >> A = MX.sym('A', 2, 3);
        >> At = A';
        >> At, size(At), At'
        A'
        3 2
        A
        \enddoctest

        \sa is_transpose, transpose

        \identifier{2qn} */
    MX T() const;

    /** \brief Get the nonzero mapping of a GetNonzeros or SetNonzeros node

        For a slicing node (x[...], x.nz[...]) returns the indices of the nonzeros of the
        dependency that make up the result; for an assignment node, the target indices.
        Any other node raises an error.

        \doctest
        >>> X = MX.sym("X", 3, 3)
        >>> X[0, 1].mapping(), X[:, 1].mapping()
        (DM(3), DM([3, 4, 5]))
        >>> try:
        ...   X.mapping()
        ... except Exception as e:
        ...   print("'mapping' not defined" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> X = MX.sym('X', 3, 3);
        >> a = X(1, 2); b = X(:, 2);
        >> a.mapping(), b.mapping()
        3
        [3, 4, 5]
        >> try, X.mapping(); catch e, disp('error'); end
        error
        \enddoctest

        \sa info

        \identifier{ri} */
    Matrix<casadi_int> mapping() const;

    /** \brief Set the depth to which equalities are checked for simplifications

        Global setting (default 1) used when building expressions: operations such as
        x - y or x/y test whether their operands are equal up to this depth with is_equal
        and simplify to 0 or 1 when they are. Depth 0 compares node identity only; each
        extra level compares the operation and recurses into the dependencies. Larger
        depths find more simplifications at a cost per created node. Calling without an
        argument restores the default.

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> MX.get_max_depth(), (x + y) - (x + y)
        (1, MX(0))
        >>> MX.set_max_depth(0)
        >>> (x + y) - (x + y)
        MX(((x+y)-(x+y)))
        >>> MX.set_max_depth()
        >>> MX.get_max_depth()
        1
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> MX.get_max_depth(), (x + y) - (x + y)
        1
        0
        >> MX.set_max_depth(0);
        >> (x + y) - (x + y)
        ((x+y)-(x+y))
        >> MX.set_max_depth();
        >> MX.get_max_depth()
        1
        \enddoctest

        \sa get_max_depth, is_equal

        \identifier{rj} */
    static void set_max_depth(casadi_int eq_depth=1);

    /** \brief Get the depth to which equalities are checked for simplifications

        \sa set_max_depth

        \identifier{rk} */
    static casadi_int get_max_depth();

    /// Check if a SharedObjectInternal pointer refers to an MX node (C++ only)
    static bool test_cast(const SharedObjectInternal* ptr);

    /** \brief Get the symbolic inputs of a Function, as MX

        Same as f.mx_in(). For an MX Function these are the original input expressions;
        for other Functions fresh symbols with the input names and sparsities.

        \doctest
        >>> x = MX.sym("x")
        >>> f = Function("f", [x], [2*x])
        >>> MX.get_input(f), is_equal(MX.get_input(f)[0], x)
        ([MX(x)], True)
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x');
        >> f = Function('f', {x}, {2*x});
        >> v = MX.get_input(f);
        >> numel(v), v{1}, is_equal(v{1}, x)
        1
        x
        1
        \enddoctest

        \sa Function.mx_in, get_free

        \identifier{rl} */
    static std::vector<MX> get_input(const Function& f);

    /** \brief Get the free variables of a Function

        Same as f.free_mx(): the symbolic primitives that occur in the outputs of an MX
        Function without being an input. Such a Function can only be created with the
        option "allow_free" and cannot be evaluated until the free variables are bound.

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> f = Function("f", [x], [x + y], {"allow_free": True})
        >>> MX.get_free(f), MX.get_free(Function("g", [x], [x]))
        ([MX(y)], [])
        >>> try:
        ...   Function("f", [x], [x + y])
        ... except Exception as e:
        ...   print("variables [y] are free" in str(e))
        True
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> f = Function('f', {x}, {x + y}, struct('allow_free', true));
        >> v = MX.get_free(f);
        >> numel(v), v{1}
        1
        y
        >> numel(MX.get_free(Function('g', {x}, {x})))
        0
        >> try, Function('f', {x}, {x + y}); catch e, disp('free variables'); end
        free variables
        \enddoctest

        \sa Function.free_mx, symvar

        \identifier{rm} */
    static std::vector<MX> get_free(const Function& f);

    /// Readability typedef
    typedef std::map<std::string, MX> MXDict;

    /** \brief Re-apply the operation of this node to new dependencies

        Creates the same operation as this node, applied to \p arg instead of the node's
        own dependencies, with the usual simplifications. \p arg must have exactly
        n_dep() entries of matching shapes; a mismatch is a developer error. Symbolic
        primitives and constants return the node itself. Call nodes are re-called with
        the new arguments; OutputNodes cannot be evaluated this way. Mainly useful for
        writing graph transformations.

        \param arg New dependencies
        \param unique Internal hints that the arguments are not referenced elsewhere;
               leave empty
        \return List with n_out() results

        \doctest
        >>> x = MX.sym("x")
        >>> y = MX.sym("y")
        >>> (x * y).eval_mx([MX(2), MX(3)]), (x + y).eval_mx([y, x])
        ([MX(6)], [MX((y+x))])
        >>> MX.sym("v", 3)[1].eval_mx([MX([4, 5, 6])])
        [MX([4, 5, 6][1])]
        \enddoctest

        \doctest{matlab}
        >> x = MX.sym('x'); y = MX.sym('y');
        >> e = x * y; r = e.eval_mx({MX(2), MX(3)});
        >> numel(r), r{1}
        1
        6
        >> e = x + y; r = e.eval_mx({y, x});
        >> r{1}
        (y+x)
        >> v = MX.sym('v', 3); e = v(2); r = e.eval_mx({MX([4; 5; 6])});
        >> r{1}
        [4, 5, 6][1]
        \enddoctest

        \sa dep, n_dep, graph_substitute

        \identifier{rn} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& SWIG_OUTPUT(res),
        const std::vector<bool>& unique=std::vector<bool>()) const;

#ifndef SWIG
    ///@{
    /** \brief Called from MXFunction

        \identifier{ro} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                    std::vector<std::vector<MX> >& fsens) const;
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                    std::vector<std::vector<MX> >& asens) const;
    ///@}

    /// Construct constant matrix with a given sparsity and values
    MX(const Sparsity& sp, double val, bool dummy);

    // Create matrix symbolic primitive
    static MX _sym(const std::string& name, const Sparsity& sp);

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex& get_mutex_temp() { return mutex_temp; }
    static std::mutex mutex_temp;
#endif //CASADI_WITH_THREADSAFE_SYMBOLICS
  private:

    /// Create an expression from a node: extra dummy arguments to avoid ambiguity for 0/NULL
    MX(MXNode* node, bool dummy1, bool dummy2, bool dummy3, bool dummy4);

    // Depth when checking equalities
    static casadi_int eq_depth_;

#endif // SWIG
  };


  ///@{
  /// Readability typedefs
  typedef std::vector<MX> MXVector;
  typedef std::initializer_list<MX> MXIList;
  typedef std::vector<MXVector> MXVectorVector;
  typedef std::map<std::string, MX> MXDict;
  ///@}

  /** \brief Kronecker contraction, the adjoint of kron

      Contracts a matrix \p m whose dimensions are multiples of those of \p x
      (as produced by kron) against \p x. With mB-by-nB the shape of \p x and blocks
      M[i*mB+r, j*nB+s] of \p m:

      - \p inner true: Y[i, j] = sum over (r, s) of M[i*mB+r, j*nB+s] * X[r, s]
      - \p inner false: Y[r, s] = sum over (i, j) of X[i, j] * M[i*mB+r, j*nB+s]

      Closes the AD algebra of kron under arbitrary-order differentiation. A shape that
      is not a multiple raises an error.

      \doctest
      >>> M = MX.sym("M", 4, 4)
      >>> X = MX.sym("X", 2, 2)
      >>> f = Function("f", [M, X], [kron_contract(M, X, True), kron_contract(M, X, False)])
      >>> f(DM(np.arange(16).reshape(4, 4)), DM.eye(2))
      (DM(
      [[5, 9],
       [21, 25]]), DM(
      [[10, 12],
       [18, 20]]))
      \enddoctest

      \doctest{matlab}
      >> M = MX.sym('M', 4, 4); X = MX.sym('X', 2, 2);
      >> f = Function('f', {M, X}, {kron_contract(M, X, true), kron_contract(M, X, false)});
      >> [a, b] = f(DM(reshape(0:15, 4, 4)'), DM.eye(2));
      >> a, b
      [[5, 9],
       [21, 25]]
      [[10, 12],
       [18, 20]]
      \enddoctest

      \sa kron

      \identifier{2ho} */
  inline MX kron_contract(const MX& m, const MX& x, bool inner) {
    return MX::kron_contract(m, x, inner);
  }

#ifndef SWIG
  /** \brief Export expressions as .html, .dot, or .casadi_viz, defaulting to expression view.

   * Options are forwarded to Function::export_graph.

      \identifier{2k6} */
  ///@{
  /// Without a filename, return the .casadi_viz JSON bundle.
  CASADI_EXPORT std::string export_graph(const MX& expression, const Dict& opts=Dict());
  CASADI_EXPORT std::string export_graph(const std::vector<MX>& expressions,
    const Dict& opts=Dict());
  CASADI_EXPORT void export_graph(const MX& expression, const std::string& fname,
    const Dict& opts=Dict());
  CASADI_EXPORT void export_graph(const std::vector<MX>& expressions, const std::string& fname,
    const Dict& opts=Dict());
  ///@}
#endif // SWIG

} // namespace casadi

#endif // CASADI_MX_HPP
