/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
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


#ifndef CASADI_LINSOL_HPP
#define CASADI_LINSOL_HPP

#include "function.hpp"
#include "printable.hpp"

namespace casadi {

  // Forward declaration of internal class
  class LinsolInternal;

  /** \brief Linear solver

      Solves the square linear system A*X = B or A'*X = B for X with a sparse direct solver
      plugin. A Linsol is created for a fixed sparsity pattern of A; the numerical values are
      given at every call. A matrix with another sparsity pattern is projected onto the fixed
      one (entries outside it are dropped). The solve is split in a symbolic factorization
      (sfact, pivot selection from the pattern), a numeric factorization (nfact) and the
      triangular solves (solve); solve(A, B) performs all three. A numerically singular matrix
      makes nfact (or solve) raise an error for most plugins; for a symbolic (MX) A, solve
      returns an expression node that is differentiable with respect to A and B.

      Available plugins include "qr" and "ldl" (sparse direct QR and LDL, codegen support),
      "symbolicqr" (QR as a symbolic expression), "csparse", "csparsecholesky", "lapacklu",
      "lapackqr", "ma27", "mumps". neig and rank are only implemented by "ldl", det by "qr"
      and "csparse". The constructor raises an error for a non-square pattern only when the
      plugin checks it; otherwise the solve fails.

      \lang{python}
      \doctest
      >>> A = DM([[3, 1], [1, 2]]); b = DM([4, 3])
      >>> L = Linsol("L", "qr", A.sparsity())
      >>> L.plugin_name(), L.sparsity()
      ('qr', Sparsity(2x2))
      >>> L.solve(A, b), L.solve(DM([[1, 2], [0, 1]]), DM([1, 1]), True)
      (DM([1, 1]), DM([1, -1]))
      >>> L.solve(A, DM([[4, 1], [3, 0]]))
      DM(
      [[1, 0.4],
       [1, -0.2]])
      >>> L.det(A)
      5.0
      >>> try:
      ...     L.solve(DM([[1, 1], [1, 1]]), b)
      ... except RuntimeError as e:
      ...     print("'nfact' failed" in str(e))
      True
      >>> Am = MX.sym("A", A.sparsity()); bm = MX.sym("b", 2)
      >>> f = Function("f", [Am, bm], [L.solve(Am, bm), L.det(Am)])
      >>> f(A, b)
      (DM([1, 1]), DM(5))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> A = DM([3 1; 1 2]); b = DM([4; 3]);
      >> L = Linsol('L', 'qr', A.sparsity());
      >> L.plugin_name(), L.sparsity()
      qr
      2x2
      >> L.solve(A, b), L.solve(DM([1 2; 0 1]), DM([1; 1]), true)
      [1, 1]
      [1, -1]
      >> L.solve(A, DM([4 1; 3 0]))
      [[1, 0.4],
       [1, -0.2]]
      >> L.det(A)
      5
      >> try, L.solve(DM([1 1; 1 1]), b); catch e, msg = e.message; end
      >> contains(msg, '''nfact'' failed')
      1
      >> Am = MX.sym('A', A.sparsity()); bm = MX.sym('b', 2);
      >> f = Function('f', {Am, bm}, {L.solve(Am, bm), L.det(Am)});
      >> [xs, d] = f(A, b); disp(xs), disp(d)
      [1, 1]
      5
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> A = ca.DM([3 1; 1 2]); b = ca.DM([4, 3]);
      julia> L = ca.Linsol("L", "qr", ca.sparsity(A));
      julia> ca.plugin_name(L), ca.sparsity(L)
      ("qr", 2x2)
      julia> ca.solve(L, A, b)
      DM([1, 1])
      julia> ca.solve(L, ca.DM([1 2; 0 1]), ca.DM([1, 1]), true)
      DM([1, -1])
      julia> ca.solve(L, A, ca.DM([4 1; 3 0]))
      DM(
      [[1, 0.4],
       [1, -0.2]])
      julia> L.det(A)
      5.0
      julia> try
                 ca.solve(L, ca.DM([1 1; 1 1]), b, false)
             catch e
                 println(occursin("'nfact' failed", sprint(showerror, e)))
             end
      true
      julia> Am = ca.MX.sym("A", ca.sparsity(A)); bm = ca.MX.sym("b", 2);
      julia> f = ca.Function("f", [Am, bm], [ca.solve(L, Am, bm), L.det(Am)]);
      julia> f(A, b)
      ([1, 1], 5)
      \enddoctest
      \endlang

      \generalsection{Linsol}
      \pluginssection{Linsol}

      \author Joel Andersson
      \date 2011-2016

      \identifier{1kh} */
  class CASADI_EXPORT Linsol
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<Linsol>) {
  public:
    /** \brief Get type name

        \identifier{1ki} */
    static std::string type_name() {return "Linsol";}

    /// Default constructor (a null object)
    Linsol();

    /** \brief Create a linear solver

        \param name Name of the solver object
        \param solver Plugin name, e.g. "qr", "ldl", "csparse"; see has_linsol
        \param sp Sparsity pattern of A (square); the numerical values are given at every call
        \param opts Options; an unknown option raises an error

        \lang{python}
        \doctest
        >>> L = Linsol("L", "ldl", Sparsity.dense(2, 2))
        >>> L
        Linsol(Linear solver LinsolLdl)
        >>> try:
        ...     Linsol("L", "qr", Sparsity.dense(2, 2), {"no_such_option": 1})
        ... except RuntimeError as e:
        ...     print("Unknown option" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> L = Linsol('L', 'ldl', Sparsity.dense(2, 2)); disp(L)
        Linear solver LinsolLdl
        >> opts = struct('no_such_option', 1);
        >> try, Linsol('L', 'qr', Sparsity.dense(2, 2), opts); catch e, msg = e.message; end
        >> contains(msg, 'Unknown option')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> L = ca.Linsol("L", "ldl", ca.Sparsity.dense(2, 2));
        julia> ca.str(L)
        "Linear solver LinsolLdl"
        julia> try
                   ca.Linsol("L", "qr", ca.Sparsity.dense(2, 2), Dict("no_such_option" => 1))
               catch e
                   println(occursin("Unknown option", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \identifier{2pe} */
    explicit Linsol(const std::string& name, const std::string& solver,
                    const Sparsity& sp, const Dict& opts=Dict());

    /// Access functions of the node
    LinsolInternal* operator->();
    const LinsolInternal* operator->() const;

#ifndef SWIG
    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);
#endif //SWIG

    /// \cond INTERNAL
#ifndef SWIG
    /** \brief  Create from node

        \identifier{1kj} */
    static Linsol create(LinsolInternal* node);
#endif // SWIG
    /// \endcond

    /** \brief Check if a plugin is available

        Same as has_linsol.

        \lang{python}
        \doctest
        >>> Linsol.has_plugin("qr"), Linsol.has_plugin("no_such_plugin")
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Linsol.has_plugin('qr'), Linsol.has_plugin('no_such_plugin')
        1
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.Linsol.has_plugin("qr"), ca.Linsol.has_plugin("no_such_plugin")
        (true, false)
        \enddoctest
        \endlang

        \identifier{2pf} */
    static bool has_plugin(const std::string& name);

    /** \brief Explicitly load a plugin dynamically

        Same as load_linsol; raises an error when the plugin cannot be found.

        \lang{python}
        \doctest
        >>> try:
        ...     Linsol.load_plugin("no_such_plugin")
        ... except RuntimeError as e:
        ...     print("no_such_plugin" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> try, Linsol.load_plugin('no_such_plugin'); catch e, msg = e.message; end
        >> contains(msg, 'no_such_plugin')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> try
                   ca.Linsol.load_plugin("no_such_plugin")
               catch e
                   println(occursin("no_such_plugin", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \identifier{2pg} */
    static void load_plugin(const std::string& name);

    /** \brief Get solver specific documentation

        Same as doc_linsol.

        \lang{python}
        \doctest
        >>> "LDL" in Linsol.doc("ldl")
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> contains(Linsol.doc('ldl'), 'LDL')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> occursin("LDL", ca.Linsol.doc("ldl"))
        true
        \enddoctest
        \endlang

        \identifier{2ph} */
    static std::string doc(const std::string& name);

    /// Query plugin name
    std::string plugin_name() const;

    /// Get the sparsity pattern of A fixed at construction
    const Sparsity& sparsity() const;

    /** \brief Symbolic factorization of the linear system, e.g. selecting pivots

        Only uses the sparsity pattern of \p A. Called automatically by nfact and solve when
        needed; raises an error when the plugin fails.

        \identifier{2pi} */
    void sfact(const DM& A) const;

    /** \brief Numeric factorization of the linear system

        Factorizes the values of \p A (after sfact if not done yet), so that neig, rank and the
        low-level solve can use the factorization. Raises an error when the factorization fails
        with option "regularity_check" set; otherwise a failure is only reported by solve.

        \lang{python}
        \doctest
        >>> A = DM([[3, 1], [1, 2]])
        >>> L = Linsol("L", "ldl", A.sparsity())
        >>> L.nfact(A)
        >>> L.neig(A), L.rank(A)
        (0, 2)
        >>> L.nfact(DM([[-1, 0], [0, 1]]))
        >>> L.neig(DM([[-1, 0], [0, 1]]))
        1
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([3 1; 1 2]); L = Linsol('L', 'ldl', A.sparsity());
        >> L.nfact(A); fprintf('%d %d\n', L.neig(A), L.rank(A))
        0 2
        >> B = DM([-1 0; 0 1]); L.nfact(B); disp(L.neig(B))
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([3 1; 1 2]);
        julia> L = ca.Linsol("L", "ldl", ca.sparsity(A));
        julia> ca.nfact(L, A)
        julia> ca.neig(L, A), ca.rank(L, A)
        (0, 2)
        julia> ca.nfact(L, ca.DM([-1 0; 0 1]))
        julia> ca.neig(L, ca.DM([-1 0; 0 1]))
        1
        \enddoctest
        \endlang

        \identifier{2pj} */
    void nfact(const DM& A) const;

    ///@{
    /** \brief Solve the linear system A*X = B (or A'*X = B with \p tr true)

        For DM arguments, factorizes \p A and solves for all columns of \p B at once; the result
        is dense. For MX arguments, returns an expression node that is differentiable. The row
        count of \p A and \p B must match.

        \param A Square matrix with the sparsity pattern given at construction (projected
               onto it otherwise)
        \param B Right-hand side, one column per system to solve
        \param tr Solve with the transpose of \p A instead
        \return X with the shape of \p B

        \lang{python}
        \doctest
        >>> L = Linsol("L", "qr", Sparsity.dense(2, 2))
        >>> L.solve(DM([[3, 1], [1, 2]]), DM([4, 3]))
        DM([1, 1])
        >>> try:
        ...     L.solve(DM([[3, 1], [1, 2]]), DM([1, 2, 3]))
        ... except RuntimeError as e:
        ...     print("matching row count" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> L = Linsol('L', 'qr', Sparsity.dense(2, 2));
        >> L.solve(DM([3 1; 1 2]), DM([4; 3]))
        [1, 1]
        >> try, L.solve(DM([3 1; 1 2]), DM([1; 2; 3])); catch e, msg = e.message; end
        >> contains(msg, 'matching row count')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> L = ca.Linsol("L", "qr", ca.Sparsity.dense(2, 2));
        julia> ca.solve(L, ca.DM([3 1; 1 2]), ca.DM([4, 3]))
        DM([1, 1])
        julia> try
                   ca.solve(L, ca.DM([3 1; 1 2]), ca.DM([1, 2, 3]))
               catch e
                   println(occursin("matching row count", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \identifier{2pk} */
    DM solve(const DM& A, const DM& B, bool tr=false) const;
    MX solve(const MX& A, const MX& B, bool tr=false) const;
    ///@}

    /** \brief Number of negative eigenvalues

      * Not available for all solvers (only "ldl"); uses the factorization of the most recent
      * nfact call, which must have been made with the same matrix.

        \lang{python}
        \doctest
        >>> A = DM([[-2, 0], [0, 1]])
        >>> L = Linsol("L", "ldl", A.sparsity())
        >>> L.nfact(A)
        >>> L.neig(A)
        1
        >>> try:
        ...     Linsol("L", "qr", A.sparsity()).neig(A)
        ... except RuntimeError as e:
        ...     print("'neig' not defined" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([-2 0; 0 1]); L = Linsol('L', 'ldl', A.sparsity());
        >> L.nfact(A); disp(L.neig(A))
        1
        >> L = Linsol('L', 'qr', A.sparsity());
        >> try, L.neig(A); catch e, disp(contains(e.message, '''neig'' not defined')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([-2 0; 0 1]);
        julia> L = ca.Linsol("L", "ldl", ca.sparsity(A));
        julia> ca.nfact(L, A)
        julia> ca.neig(L, A)
        1
        julia> try
                   ca.neig(ca.Linsol("L", "qr", ca.sparsity(A)), A)
               catch e
                   println(occursin("'neig' not defined", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \identifier{1kk} */
    casadi_int neig(const DM& A) const;

    /** \brief Matrix rank

      * Not available for all solvers (only "ldl"); uses the factorization of the most recent
      * nfact call, which must have been made with the same matrix.

        A rank-deficient matrix also triggers a warning in nfact.

        \lang{python}
        \doctest
        >>> A = DM([[2, 1], [1, 3]])
        >>> L = Linsol("L", "ldl", A.sparsity())
        >>> L.nfact(A)
        >>> L.rank(A)
        2
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([2 1; 1 3]); L = Linsol('L', 'ldl', A.sparsity());
        >> L.nfact(A); disp(L.rank(A))
        2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([2 1; 1 3]);
        julia> L = ca.Linsol("L", "ldl", ca.sparsity(A));
        julia> ca.nfact(L, A)
        julia> ca.rank(L, A)
        2
        \enddoctest
        \endlang

        \identifier{1kl} */
    casadi_int rank(const DM& A) const;

    /** \brief Matrix determinant

      * Not available for all solvers ("qr" and "csparse"); factorizes \p A first.

        \lang{python}
        \doctest
        >>> A = DM([[3, 1], [1, 2]])
        >>> Linsol("L", "qr", A.sparsity()).det(A)
        5.0
        >>> try:
        ...     Linsol("L", "ldl", A.sparsity()).det(A)
        ... except RuntimeError as e:
        ...     print("'det' not defined" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> A = DM([3 1; 1 2]); L = Linsol('L', 'qr', A.sparsity()); disp(L.det(A))
        5
        >> L = Linsol('L', 'ldl', A.sparsity());
        >> try, L.det(A); catch e, disp(contains(e.message, '''det'' not defined')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.DM([3 1; 1 2]);
        julia> ca.Linsol("L", "qr", ca.sparsity(A)).det(A)
        5.0
        julia> try
                   ca.Linsol("L", "ldl", ca.sparsity(A)).det(A)
               catch e
                   println(occursin("'det' not defined", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \identifier{2hw} */
    double det(const DM& A) const;

    /// Matrix determinant as an MX graph node
    MX det(const MX& A) const;

    /** \brief Get all statistics obtained at the end of the last evaluate call

        Plugin-specific; empty for most plugins. Raises an error before the first solve.
        With the option "record_time" set, it holds call counts and timings (n_call_*,
        t_proc_*, t_wall_*) of sfact, nfact and solve.

        \lang{python}
        \doctest
        >>> L = Linsol("L", "qr", Sparsity.dense(1, 1))
        >>> try:
        ...     L.stats()
        ... except RuntimeError as e:
        ...     print("did not solve a problem yet" in str(e))
        True
        >>> L.solve(DM(2), DM(4)), L.stats()
        (DM(2), {})
        >>> L = Linsol("L", "qr", Sparsity.dense(1, 1), {"record_time": True})
        >>> L.solve(DM(2), DM(4)); L.stats()["n_call_solve"]
        DM(2)
        1
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> L = Linsol('L', 'qr', Sparsity.dense(1, 1));
        >> try, L.stats(); catch e, disp(contains(e.message, 'did not solve a problem yet')); end
        1
        >> L.solve(DM(2), DM(4)), isempty(fieldnames(L.stats()))
        2
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> L = ca.Linsol("L", "qr", ca.Sparsity.dense(1, 1));
        julia> try
                   ca.stats(L)
               catch e
                   println(occursin("did not solve a problem yet", sprint(showerror, e)))
               end
        true
        julia> ca.solve(L, ca.DM(2), ca.DM(4)), ca.stats(L)
        (2, Dict{String, Any}())
        julia> L = ca.Linsol("L", "qr", ca.Sparsity.dense(1, 1), Dict("record_time" => true));
        julia> ca.solve(L, ca.DM(2), ca.DM(4)); ca.stats(L)["n_call_solve"]
        1
        \enddoctest
        \endlang

        \identifier{2pl} */
    Dict stats(int mem=1) const;

    #ifndef SWIG
    ///@{
    /// Low-level API
    int sfact(const double* A, int mem=0) const;
    int nfact(const double* A, int mem=0) const;
    int solve(const double* A, double* x, casadi_int nrhs=1, bool tr=false, int mem=0) const;
    double det(const double* A, int mem = 0) const;
    casadi_int neig(const double* A, int mem=0) const;
    casadi_int rank(const double* A, int mem=0) const;
    ///@}

    /// Checkout a memory object
    casadi_int checkout() const;

    /// Release a memory object
    void release(int mem) const;

    /** \brief Serialize an object

        \identifier{1km} */
    void serialize(SerializingStream &s) const;

    /** \brief Deserialize with type disambiguation

        \identifier{1kn} */
    static Linsol deserialize(DeserializingStream& s);

    #endif // SWIG
  };

  /** \brief Check if a linear solver plugin is available

      Tries to load the plugin \p name if it is not loaded yet; false when it is not found.

      \lang{python}
      \doctest
      >>> has_linsol("qr"), has_linsol("no_such_plugin")
      (True, False)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> has_linsol('qr'), has_linsol('no_such_plugin')
      1
      0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.has_linsol("qr"), ca.has_linsol("no_such_plugin")
      (true, false)
      \enddoctest
      \endlang

      \identifier{2pm} */
  CASADI_EXPORT bool has_linsol(const std::string& name);

  /** \brief Explicitly load a linear solver plugin

      Normally plugins are loaded on first use; this raises an error when the plugin cannot be
      found.

      \lang{python}
      \doctest
      >>> try:
      ...     load_linsol("no_such_plugin")
      ... except RuntimeError as e:
      ...     print("no_such_plugin" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> try, load_linsol('no_such_plugin'); catch e, msg = e.message; end
      >> contains(msg, 'no_such_plugin')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> try
                 ca.load_linsol("no_such_plugin")
             catch e
                 println(occursin("no_such_plugin", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \identifier{2pn} */
  CASADI_EXPORT void load_linsol(const std::string& name);

  /** \brief Get the documentation string of a linear solver plugin

      The text includes the plugin-specific option table.

      \lang{python}
      \doctest
      >>> "QR" in doc_linsol("qr")
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> contains(doc_linsol('qr'), 'QR')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> occursin("QR", ca.doc_linsol("qr"))
      true
      \enddoctest
      \endlang

      \identifier{2po} */
  CASADI_EXPORT std::string doc_linsol(const std::string& name);

} // namespace casadi

#endif // CASADI_LINSOL_HPP
