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

#ifndef CASADI_FUNCTION_HPP
#define CASADI_FUNCTION_HPP

#include "sx_elem.hpp"
#include "mx.hpp"
#include "printable.hpp"
#include <exception>
#include <stack>

namespace casadi {

#ifndef SWIG
  /** Forward declaration of internal class */
  class FunctionInternal;
  class SerializingStream;
  class DeserializingStream;
#endif // SWIG

  /** \brief Function object

      A Function instance is a general multiple-input, multiple-output function
      where each input and output can be a sparse matrix.\n

      For an introduction to this class, see the CasADi user guide.\n

      Function is a reference counted and immutable class; copying a class instance
      is very cheap and its behavior (with some exceptions) is not affected by
      calling its member functions.\n

      \author Joel Andersson
      \date 2010-2017

      \identifier{1uw} */
  class CASADI_EXPORT Function :
    public SharedObject,
    public SWIG_IF_ELSE(PrintableCommon, Printable<Function>) {
  public:
    /** \brief Get type name

        \identifier{1ux} */
    static std::string type_name() {return "Function";}

    /** \brief Default constructor, null pointer

        \identifier{1uy} */
    Function();

    /** \brief Construct from a file

        \identifier{1uz} */
    Function(const std::string& fname);

    ///@{
    /** \brief Construct an SX function

        \identifier{1v0} */
    Function(const std::string& name,
             const std::vector<SX>& ex_in,
             const std::vector<SX>& ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name,
             const std::vector<SX>& ex_in,
             const std::vector<SX>& ex_out,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const std::map<std::string, SX>& dict,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    ///@}

    ///@{
    /** \brief Construct an MX function

        \identifier{1v1} */
    Function(const std::string& name,
             const std::vector<MX>& ex_in,
             const std::vector<MX>& ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name,
             const std::vector<MX>& ex_in,
             const std::vector<MX>& ex_out,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const std::map<std::string, MX>& dict,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    ///@}

    ///@{
    /** \brief To resolve ambiguity on some compilers

        \identifier{1v2} */
#ifndef SWIG
    Function(const std::string& name, SXIList ex_in,
      const SXVector& ex_out, const Dict& opts=Dict());
    Function(const std::string& name, const SXVector& ex_in,
      SXIList ex_out, const Dict& opts=Dict());
    Function(const std::string& name, SXIList ex_in,
      SXIList ex_out, const Dict& opts=Dict());
    Function(const std::string& name, SXIList ex_in, const SXVector& ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const SXVector& ex_in, SXIList ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, SXIList ex_in, SXIList ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, MXIList ex_in, const MXVector& ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const MXVector& ex_in, MXIList ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name, MXIList ex_in, MXIList ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name, MXIList ex_in, const MXVector& ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const MXVector& ex_in, MXIList ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, MXIList ex_in, MXIList ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
#endif // SWIG
    ///@}

    ///@{
    /** \brief Create a just-in-time compiled function from a C language string

     * The names and sparsity patterns of all the inputs and outputs must be provided.
     * If sparsities are not provided, all inputs and outputs are assumed to be scalar.
     * Only specify the function body, assuming that input and output nonzeros are
     * stored in arrays with the specified naming convension.
     * The data type used is 'casadi_real', which is typically equal to 'double` or
     * another data type with the same API as 'double'.
     *
     * Inputs may be null pointers. This means that the all entries are zero.
     * Outputs may be null points. This means that the corresponding result can be ignored.
     *
     * If an error occurs in the evaluation, issue "return 1;";
     *
     * The final generated function will have a structure similar to:
     *
     * casadi_int fname(const casadi_real** arg, casadi_real** res, casadi_int* iw,
                 casadi_real* w, void* mem) {
     *   const casadi_real *x1, *x2;
     *   casadi_real *r1, *r2;
     *   x1 = *arg++;
     *   x2 = *arg++;
     *   r1 = *res++;
     *   r2 = *res++;
     *   <FUNCTION_BODY>
     *   return 0;
     * }
     *
        \identifier{1v3} */
    static Function jit(const std::string& name, const std::string& body,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    static Function jit(const std::string& name, const std::string& body,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const std::vector<Sparsity>& sparsity_in,
             const std::vector<Sparsity>& sparsity_out,
             const Dict& opts=Dict());
    ///@}

    /** \brief  Destructor

        \identifier{1v4} */
    ~Function();

    /** \brief Expand a function to SX

        \identifier{1v5} */
    ///@{
    Function expand() const;
    Function expand(const std::string& name,
                    const Dict& opts=Dict()) const;
    ///@}

    /// \cond INTERNAL
#ifndef SWIG
    /** \brief  Create from node

        \identifier{1v6} */
    static Function create(FunctionInternal* node);

    /** \brief  Create from node and initialize

        \identifier{1v7} */
    static Function create(FunctionInternal* node, const Dict& opts);
#endif // SWIG
    /// \endcond

    /** \brief Get the number of function inputs

        \identifier{1v8} */
    casadi_int n_in() const;

    /** \brief Get the number of function outputs

        \identifier{1v9} */
    casadi_int n_out() const;

    ///@{
    /** \brief Get input dimension

        \identifier{1va} */
    casadi_int size1_in(casadi_int ind) const;
    casadi_int size1_in(const std::string& iname) const {return size1_in(index_in(iname));}
    casadi_int size2_in(casadi_int ind) const;
    casadi_int size2_in(const std::string& iname) const {return size2_in(index_in(iname));}
    std::pair<casadi_int, casadi_int> size_in(casadi_int ind) const;
    std::pair<casadi_int, casadi_int> size_in(const std::string& iname) const {
      return size_in(index_in(iname));
    }
    ///@}

    ///@{
    /** \brief Get output dimension

        \identifier{1vb} */
    casadi_int size1_out(casadi_int ind) const;
    casadi_int size1_out(const std::string& oname) const {return size1_out(index_out(oname));}
    casadi_int size2_out(casadi_int ind) const;
    casadi_int size2_out(const std::string& oname) const {return size2_out(index_out(oname));}
    std::pair<casadi_int, casadi_int> size_out(casadi_int ind) const;
    std::pair<casadi_int, casadi_int> size_out(const std::string& oname) const {
      return size_out(index_out(oname));
    }
    ///@}

    ///@{
    /** \brief  Get number of input nonzeros
     *
     * For a particular input or for all of the inputs

        \identifier{1vc} */
    casadi_int nnz_in() const;
    casadi_int nnz_in(casadi_int ind) const;
    casadi_int nnz_in(const std::string& iname) const {return nnz_in(index_in(iname));}
    ///@}

    ///@{
    /** \brief  Get number of output nonzeros
     *
     * For a particular output or for all of the outputs

        \identifier{1vd} */
    casadi_int nnz_out() const;
    casadi_int nnz_out(casadi_int ind) const;
    casadi_int nnz_out(const std::string& oname) const {return nnz_out(index_out(oname));}
    ///@}

    ///@{
    /** \brief  Get number of input elements
     *
     * For a particular input or for all of the inputs

        \identifier{1ve} */
    casadi_int numel_in() const;
    casadi_int numel_in(casadi_int ind) const;
    casadi_int numel_in(const std::string& iname) const {return numel_in(index_in(iname));}
    ///@}

    ///@{
    /** \brief  Get number of output elements
     *
     * For a particular output or for all of the outputs

        \identifier{1vf} */
    casadi_int numel_out() const;
    casadi_int numel_out(casadi_int ind) const;
    casadi_int numel_out(const std::string& oname) const {return numel_out(index_out(oname));}
    ///@}

    /** \brief Get input scheme

        \identifier{1vg} */
    const std::vector<std::string>& name_in() const;

    /** \brief Get output scheme

        \identifier{1vh} */
    const std::vector<std::string>& name_out() const;

    /** \brief Get input scheme name by index

        \identifier{1vi} */
    const std::string& name_in(casadi_int ind) const;

    /** \brief Get output scheme name by index

        \identifier{1vj} */
    const std::string& name_out(casadi_int ind) const;

    /** \brief Find the index for a string describing a particular entry of an input scheme
     *
     * example:  schemeEntry("x_opt")  -> returns  NLPSOL_X if FunctionInternal adheres to
     * SCHEME_NLPINput

        \identifier{1vk} */
    casadi_int index_in(const std::string &name) const;

    /** \brief Find the index for a string describing a particular entry of an output scheme
     *
     * example:  schemeEntry("x_opt")  -> returns  NLPSOL_X if FunctionInternal adheres to
     * SCHEME_NLPINput

        \identifier{1vl} */
    casadi_int index_out(const std::string &name) const;

    /** \brief Does the function have a particularly named input?

        \identifier{2c9} */
    bool has_in(const std::string &name) const;
    /** \brief Does the function have a particularly named output?

        \identifier{2ca} */
    bool has_out(const std::string &name) const;

    /** \brief Get default input value

        \identifier{1vm} */
    double default_in(casadi_int ind) const;

    /** \brief Get largest input value

        \identifier{1vn} */
    double max_in(casadi_int ind) const;

    /** \brief Get smallest input value

        \identifier{1vo} */
    double min_in(casadi_int ind) const;

    /** \brief Get nominal input value

        \identifier{1vp} */
    std::vector<double> nominal_in(casadi_int ind) const;

    /** \brief Get nominal output value

        \identifier{1vq} */
    std::vector<double> nominal_out(casadi_int ind) const;

    /** \brief Get sparsity of a given input

        \identifier{1vr} */
    /// @{
    const Sparsity& sparsity_in(casadi_int ind) const;
    const Sparsity& sparsity_in(const std::string& iname) const;
    /// @}

    /** \brief Get sparsity of a given output

        \identifier{1vs} */
    /// @{
    const Sparsity& sparsity_out(casadi_int ind) const;
    const Sparsity& sparsity_out(const std::string& iname) const;
    /// @}

    /** \brief Get differentiability of inputs/output

        \identifier{1vt} */
    /// @{
    bool is_diff_in(casadi_int ind) const;
    bool is_diff_out(casadi_int ind) const;
    std::vector<bool> is_diff_in() const;
    std::vector<bool> is_diff_out() const;
    /// @}

    // A linear combination of inputs
    typedef std::map<std::string, std::vector<std::string> > AuxOut;

    // Factory
    Function factory(const std::string& name,
                     const std::vector<std::string>& s_in,
                     const std::vector<std::string>& s_out,
                     const AuxOut& aux=AuxOut(),
                     const Dict& opts=Dict()) const;

    /** \brief Get oracle

        \identifier{1vu} */
    Function oracle() const;

    /** \brief Wrap in an Function instance consisting of only one MX call

        \identifier{1vv} */
    Function wrap() const;

    /** \brief Wrap in a Function with options

        \identifier{1vw} */
    Function wrap_as_needed(const Dict& opts) const;

    /** \brief Which variables enter with some order
    *
    * \param[in] order Only 1 (linear) and 2 (nonlinear) allowed
    * \param[in] tr   Flip the relationship. Return which expressions contain the variables

        \identifier{1vx} */
    std::vector<bool> which_depends(const std::string& s_in,
                                    const std::vector<std::string>& s_out,
                                    casadi_int order=1, bool tr=false) const;

    /** \brief Print dimensions of inputs and outputs

        \identifier{1vy} */
    void print_dimensions(std::ostream &stream=casadi::uout()) const;

    /** \brief Print options to a stream

        \identifier{1vz} */
    void print_options(std::ostream &stream=casadi::uout()) const;

    /** \brief Print all information there is to know about a certain option

        \identifier{1w0} */
    void print_option(const std::string &name, std::ostream &stream = casadi::uout()) const;

    /** \brief Does a particular option exist

        \identifier{1w1} */
    bool has_option(const std::string &option_name) const;

    /** \brief Change option after object creation for debugging

      * This is only possible for a selected number of options that do not change the numerical
      * results of the computation, e.g. to enable a more verbose output or saving to file.

        \identifier{1w2} */
    void change_option(const std::string& option_name, const GenericType& option_value);

    /** \brief Reset the counter used to name dump files

        \identifier{2dy} */
    void reset_dump_count();

    /** \brief Do the derivative functions need nondifferentiated outputs?

        \identifier{1w3} */
    bool uses_output() const;

#ifdef WITH_DEPRECATED_FEATURES
    /** \brief [DEPRECATED] Replaced by Function::factory.

        \identifier{1w4} */
    Function jacobian_old(casadi_int iind, casadi_int oind) const;

    /** \brief [DEPRECATED] Replaced by Function::factory.

        \identifier{1w5} */
    Function hessian_old(casadi_int iind, casadi_int oind) const;

    ///@{
    /// [DEPRECATED] Get, if necessary generate, the sparsity of a Jacobian block
    const Sparsity sparsity_jac(casadi_int iind, casadi_int oind,
                                bool compact=false, bool symmetric=false) const;
    const Sparsity sparsity_jac(const std::string &iind, casadi_int oind=0,
                                bool compact=false, bool symmetric=false) const {
      return sparsity_jac(index_in(iind), oind, compact, symmetric);
    }
    const Sparsity sparsity_jac(casadi_int iind, const std::string &oind,
                                bool compact=false, bool symmetric=false) const {
      return sparsity_jac(iind, index_out(oind), compact, symmetric);
    }
    const Sparsity sparsity_jac(const std::string &iind, const std::string &oind,
                                bool compact=false, bool symmetric=false) const {
      return sparsity_jac(index_in(iind), index_out(oind), compact, symmetric);
    }
    ///@}
#endif // WITH_DEPRECATED_FEATURES

  /** \brief Calculate all Jacobian blocks

    * Generates a function that takes all non-differentiated inputs and outputs
    * and calculates all Jacobian blocks.
    * Inputs that are not needed by the routine are all-zero sparse matrices
    * with the correct dimensions. Output blocks that are not calculated,
    * e.g. if the corresponding input or output is marked non-differentiated
    * are also all-zero sparse.
    * The Jacobian blocks are sorted starting by all the blocks for the first
    * output, then all the blocks for the second output and so on.
    * E.g. f : (x, y) -> (r, s) results in the function
    * jac_f : (x, y, out_r, out_s) -> (jac_r_x, jac_r_y, jac_s_x, jac_s_y)
    *
    * This function is cached.

      \identifier{1w6} */
    Function jacobian() const;

    ///@{
    /** \brief Evaluate the function symbolically or numerically

        \identifier{1w7} */
    void call(const std::vector<DM> &arg, std::vector<DM>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const std::vector<SX> &arg, std::vector<SX>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const std::vector<MX> &arg, std::vector<MX>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const DMDict& arg, DMDict& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const SXDict& arg, SXDict& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const MXDict& arg, MXDict& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    ///@}

#ifndef SWIG
    /// Check if same as another function
    bool operator==(const Function& f) const;

    ///@{
    /// Functor shorthand for evaluation
    std::vector<DM> operator()(const std::vector<DM>& arg) const;
    std::vector<SX> operator()(const std::vector<SX>& arg) const;
    std::vector<MX> operator()(const std::vector<MX>& arg) const;
    const DMDict operator()(const DMDict& arg) const;
    const SXDict operator()(const SXDict& arg) const;
    const MXDict operator()(const MXDict& arg) const;
    ///@}

    ///@{
    /** \brief Evaluate with temporary memory allocation

        \identifier{1w8} */
    void operator()(std::vector<const double*> arg, std::vector<double*> res) const;
    void operator()(std::vector<const bvec_t*> arg, std::vector<bvec_t*> res) const;
    void operator()(std::vector<const SXElem*> arg, std::vector<SXElem*> res) const;
    template<typename D> void call_gen(std::vector<const D*> arg, std::vector<D*> res) const;
    ///@}

    ///@{
    /** \brief Supported arguments for numerical evaluation and converters

        \identifier{1w9} */
    typedef const std::vector<std::vector<double>>& VecArg;
    std::vector<const double*> buf_in(VecArg arg) const;
    typedef std::vector<std::vector<double>>& VecRes;
    std::vector<double*> buf_out(VecRes res) const;
    typedef std::vector<std::vector<double>*> VPrRes;
    std::vector<double*> buf_out(VPrRes res) const;

    typedef const std::map<std::string, std::vector<double>>& MapArg;
    std::vector<const double*> buf_in(MapArg arg) const;
    typedef std::map<std::string, std::vector<double>>& MapRes;
    std::vector<double*> buf_out(MapRes res) const;
    typedef std::map<std::string, std::vector<double>*> MPrRes;
    std::vector<double*> buf_out(MPrRes res) const;
    ///@}

    ///@{
    /** \brief Numerical evaluation

        \identifier{1wa} */
    void operator()(VecArg arg, VecRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(VecArg arg, MapRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(VecArg arg, VPrRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(VecArg arg, MPrRes res) const { (*this)(buf_in(arg), buf_out(res)); }

    void operator()(MapArg arg, VecRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(MapArg arg, MapRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(MapArg arg, VPrRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(MapArg arg, MPrRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    ///@}

    ///@{
    /// Functor shorthand for evaluation, single argument (only C++)
    std::vector<DM> operator()(const DM& arg0) const {
      return operator()(std::vector<DM>{arg0});
    }
    std::vector<SX> operator()(const SX& arg0) const {
      return operator()(std::vector<SX>{arg0});
    }
    std::vector<MX> operator()(const MX& arg0) const {
      return operator()(std::vector<MX>{arg0});
    }
    ///@}

    /** \brief Evaluate memory-less, numerically

        \identifier{1wb} */
    int operator()(const double** arg, double** res,
        casadi_int* iw, double* w, int mem) const;

    /** \brief Evaluate numerically with checkout/release

        \identifier{1wc} */
    int operator()(const double** arg, double** res,
        casadi_int* iw, double* w) const;

    /** \brief Evaluate memory-less SXElem

        Same syntax as the double version, allowing use in templated code

        \identifier{1wd} */
    int operator()(const SXElem** arg, SXElem** res,
        casadi_int* iw, SXElem* w, int mem=0) const;

    /** \brief  Propagate sparsity forward

        \identifier{1we} */
    int operator()(const bvec_t** arg, bvec_t** res,
        casadi_int* iw, bvec_t* w, int mem=0) const;

    /** \brief  Propagate sparsity backward

        \identifier{1wf} */
    int rev(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, int mem=0) const;

    /** \brief Propagate sparsity backward with temporary memory allocation

        \identifier{1wg} */
    int rev(std::vector<bvec_t*> arg, std::vector<bvec_t*> res) const;

#endif // SWIG

    /** \brief  Evaluate symbolically in parallel and sum (matrix graph)

        \param parallelization Type of parallelization used: unroll|serial|openmp

        \identifier{1wh} */
    std::vector<MX> mapsum(const std::vector<MX > &x,
                           const std::string& parallelization="serial") const;

    ///@{
    /** \brief  Create a mapaccumulated version of this function

        Suppose the function has a signature of:
        \verbatim
           f: (x, u) -> (x_next , y )
        \endverbatim

        The the mapaccumulated version has the signature:
        \verbatim
           F: (x0, U) -> (X , Y )

            with
                U: horzcat([u0, u1, ..., u_(N-1)])
                X: horzcat([x1, x2, ..., x_N])
                Y: horzcat([y0, y1, ..., y_(N-1)])

            and
                x1, y0 <- f(x0, u0)
                x2, y1 <- f(x1, u1)
                ...
                x_N, y_(N-1) <- f(x_(N-1), u_(N-1))
        \endverbatim

        Mapaccum has the following benefits over writing an equivalent for-loop:
           - much faster at construction time
           - potentially much faster compilation times (for codegen)
           - offers a trade-off between memory and evaluation time

        The base (settable through the options dictionary, default 10),
        is used to create a tower of function calls,
        containing unrolled for-loops of length maximum base.

        This technique is much more scalable in terms of memory-usage,
        but slightly slower at evaluation, than a plain for-loop.
        The effect is similar to that of a for-loop with a check-pointing instruction
        after each chunk of iterations with size base.

        Set base to -1 to unroll all the way; no gains in memory efficiency here.

        \identifier{1wi} */
    Function mapaccum(const std::string& name, casadi_int N, const Dict& opts = Dict()) const;
    Function mapaccum(const std::string& name, casadi_int N, casadi_int n_accum,
                      const Dict& opts = Dict()) const;
    Function mapaccum(const std::string& name, casadi_int n,
                      const std::vector<casadi_int>& accum_in,
                      const std::vector<casadi_int>& accum_out,
                      const Dict& opts=Dict()) const;
    Function mapaccum(const std::string& name, casadi_int n,
                      const std::vector<std::string>& accum_in,
                      const std::vector<std::string>& accum_out,
                      const Dict& opts=Dict()) const;
    Function mapaccum(casadi_int N, const Dict& opts = Dict()) const;
    Function fold(casadi_int N, const Dict& opts = Dict()) const;
    ///@}

    /** \brief  Create a mapped version of this function

        Suppose the function has a signature of:
        \verbatim
           f: (a, p) -> ( s )
        \endverbatim

        The the mapped version has the signature:
        \verbatim
           F: (A, P) -> (S )

            with
                A: horzcat([a0, a1, ..., a_(N-1)])
                P: horzcat([p0, p1, ..., p_(N-1)])
                S: horzcat([s0, s1, ..., s_(N-1)])
            and
                s0 <- f(a0, p0)
                s1 <- f(a1, p1)
                ...
                s_(N-1) <- f(a_(N-1), p_(N-1))
        \endverbatim

        \param parallelization Type of parallelization used: unroll|serial|openmp

        \identifier{1wj} */
    Function map(casadi_int n, const std::string& parallelization="serial") const;
    Function map(casadi_int n, const std::string& parallelization,
      casadi_int max_num_threads) const;

    ///@{
    /** \brief Map with reduction

      A subset of the inputs are non-repeated and a subset of the outputs summed
      up.

        \identifier{1wk} */
    Function map(const std::string& name, const std::string& parallelization, casadi_int n,
      const std::vector<casadi_int>& reduce_in,
      const std::vector<casadi_int>& reduce_out,
      const Dict& opts=Dict()) const;
    Function map(const std::string& name, const std::string& parallelization, casadi_int n,
      const std::vector<std::string>& reduce_in,
      const std::vector<std::string>& reduce_out,
      const Dict& opts=Dict()) const;
    Function map(casadi_int n,
      const std::vector<bool>& reduce_in,
      const std::vector<bool>& reduce_out=std::vector<bool>(),
      const Dict& opts=Dict()) const;
    ///@}

    /** \brief returns a new function with a selection of inputs/outputs of the original

        \identifier{1wl} */
    Function slice(const std::string& name, const std::vector<casadi_int>& order_in,
                   const std::vector<casadi_int>& order_out, const Dict& opts=Dict()) const;

    /** \brief Constuct a switch function

        \identifier{1wm} */
    static Function conditional(const std::string& name, const std::vector<Function>& f,
                                const Function& f_def, const Dict& opts=Dict());

    /** \brief Conditional call to a function

        \identifier{1wn} */
    static Function conditional(const std::string& name,
      const Function& f, const Dict& opts=Dict());

    /** \brief BSpline evaluator function
     *
     *  Requires a known coefficient tensor

        \identifier{1wo} */
    static Function bspline(const std::string &name,
      const std::vector< std::vector<double> >& knots, const std::vector<double>& coeffs,
      const std::vector<casadi_int>& degree, casadi_int m=1, const Dict& opts=Dict());

    /** \brief Constructor (if-else)

        \identifier{1wp} */
    static Function if_else(const std::string& name, const Function& f_true,
                            const Function& f_false, const Dict& opts=Dict());

    /** \brief Get a function that calculates \a nfwd forward derivatives
     *
     *         Returns a function with <tt>n_in + n_out + n_in</tt> inputs
     *         and <tt>nfwd</tt> outputs.
     *         The first <tt>n_in</tt> inputs correspond to nondifferentiated inputs.
     *         The next <tt>n_out</tt> inputs correspond to nondifferentiated outputs.
     *         and the last <tt>n_in</tt> inputs correspond to forward seeds,
     *         stacked horizontally
     *         The  <tt>n_out</tt> outputs correspond to forward sensitivities,
     *         stacked horizontally.     *
     *         <tt>(n_in = n_in(), n_out = n_out())</tt>
     *
     *        The functions returned are cached, meaning that if called multiple timed
     *        with the same value, then multiple references to the same function will be returned.

        \identifier{1wq} */
    Function forward(casadi_int nfwd) const;

    /** \brief Get a function that calculates \a nadj adjoint derivatives
     *
     *         Returns a function with <tt>n_in + n_out + n_out</tt> inputs
     *         and <tt>n_in</tt> outputs.
     *         The first <tt>n_in</tt> inputs correspond to nondifferentiated inputs.
     *         The next <tt>n_out</tt> inputs correspond to nondifferentiated outputs.
     *         and the last <tt>n_out</tt> inputs correspond to adjoint seeds,
     *         stacked horizontally
     *         The  <tt>n_in</tt> outputs correspond to adjoint sensitivities,
     *         stacked horizontally.     *
     *         <tt>(n_in = n_in(), n_out = n_out())</tt>
     *
     *         <tt>(n_in = n_in(), n_out = n_out())</tt>
     *
     *        The functions returned are cached, meaning that if called multiple timed
     *        with the same value, then multiple references to the same function will be returned.

        \identifier{1wr} */
    Function reverse(casadi_int nadj) const;

    /** \brief Get, if necessary generate, the sparsity of all Jacobian blocks

        \identifier{1ws} */
    const std::vector<Sparsity>& jac_sparsity(bool compact = false) const;

    /** \brief Get, if necessary generate, the sparsity of a single Jacobian block

        \identifier{1wt} */
    Sparsity jac_sparsity(casadi_int oind, casadi_int iind, bool compact = false) const;

    /** \brief Export / Generate C code for the function

        \identifier{1wu} */
    std::string generate(const std::string& fname, const Dict& opts=Dict()) const;

    /** \brief Export / Generate C code for the function

        \identifier{1wv} */
    std::string generate(const Dict& opts=Dict()) const;

    /** \brief Export / Generate C code for the dependency function

        \identifier{1ww} */
    std::string generate_dependencies(const std::string& fname, const Dict& opts=Dict()) const;

    /** \brief Export an input file that can be passed to generate C code with a main
     *
     * \see generate_out
     * \see convert_in to convert between dict/map and vector

        \identifier{1wx} */
    /// @{
    void generate_in(const std::string& fname, const std::vector<DM>& arg);
    std::vector<DM> generate_in(const std::string& fname);
    /// @}

    /** \brief Export an output file that can be checked with generated C code output
     *
     * \see generate_in
     * \see convert_out to convert between dict/map and vector

        \identifier{1wy} */
    /// @{
    void generate_out(const std::string& fname, const std::vector<DM>& arg);
    std::vector<DM> generate_out(const std::string& fname);
    /// @}

    /** \brief Export function in specific language
     *
     * Only allowed for (a subset of) SX/MX Functions

        \identifier{1wz} */
    ///@{
    void export_code(const std::string& lang,
      const std::string &fname, const Dict& options=Dict()) const;

#ifndef SWIG
    /** \brief Serialize

        \identifier{1x0} */
    void serialize(std::ostream &stream, const Dict& opts=Dict()) const;

    /** \brief Serialize an object

        \identifier{1x1} */
    void serialize(SerializingStream &s) const;
#endif

    /** \brief Serialize

        \identifier{1x2} */
    std::string serialize(const Dict& opts=Dict()) const;

    /** \brief Save Function to a file

        \see load

        \identifier{240} */
    void save(const std::string &fname, const Dict& opts=Dict()) const;

    std::string export_code(const std::string& lang, const Dict& options=Dict()) const;
#ifndef SWIG
    void export_code(const std::string& lang,
      std::ostream &stream, const Dict& options=Dict()) const;
#endif // SWIG
    ///@}
#ifndef SWIG
    /// \cond INTERNAL
    /// Get a const pointer to the node
    FunctionInternal* get() const;

    /// Get a pointer and typecast
    template<typename T>
    T* get() const {
      T* ret = dynamic_cast<T*>(get());
      casadi_assert_dev(ret!=nullptr);
      return ret;
    }

    /** \brief  Const access functions of the node

        \identifier{1x3} */
    FunctionInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);
    /// \endcond
#endif // SWIG

    /// Get all statistics obtained at the end of the last evaluate call
    Dict stats(int mem=0) const;

    ///@{
    /** \brief Get symbolic primitives equivalent to the input expressions

     * There is no guarantee that subsequent calls return unique answers

        \identifier{1x4} */
    const SX sx_in(casadi_int iind) const;
    const SX sx_in(const std::string& iname) const {
      return sx_in(index_in(iname));
    }
    const std::vector<SX> sx_in() const;
    const MX mx_in(casadi_int ind) const;
    const MX mx_in(const std::string & iname) const {
      return mx_in(index_in(iname));
    }
    const std::vector<MX> mx_in() const;
    ///@}

    ///@{
    /** \brief Get symbolic primitives equivalent to the output expressions

    * There is no guarantee that subsequent calls return unique answers

        \identifier{1x5} */
    const SX sx_out(casadi_int oind) const;
    const SX sx_out(const std::string& oname) const {
      return sx_out(index_out(oname));
    }
    const std::vector<SX> sx_out() const;
    const MX mx_out(casadi_int ind) const;
    const MX mx_out(const std::string& oname) const {
      return mx_out(index_out(oname));
    }
    const std::vector<MX> mx_out() const;
    ///@}

    /** \brief Convert from/to flat vector of input/output nonzeros

        \identifier{1x6} */
    /// @{
    std::vector<double> nz_from_in(const std::vector<DM>& arg) const;
    std::vector<double> nz_from_out(const std::vector<DM>& arg) const;
    std::vector<DM> nz_to_in(const std::vector<double>& arg) const;
    std::vector<DM> nz_to_out(const std::vector<double>& arg) const;
    ///@}

    /** \brief Convert from/to input/output lists/map
    *
    * Will raise an error when an unknown key is used or a list has incorrect size.
    * Does not perform sparsity checking.

        \identifier{1x7} */
    /// @{
    DMDict convert_in(const std::vector<DM>& arg) const;
    std::vector<DM> convert_in(const DMDict& arg) const;
    DMDict convert_out(const std::vector<DM>& arg) const;
    std::vector<DM> convert_out(const DMDict& arg) const;
    SXDict convert_in(const std::vector<SX>& arg) const;
    std::vector<SX> convert_in(const SXDict& arg) const;
    SXDict convert_out(const std::vector<SX>& arg) const;
    std::vector<SX> convert_out(const SXDict& arg) const;
    MXDict convert_in(const std::vector<MX>& arg) const;
    std::vector<MX> convert_in(const MXDict& arg) const;
    MXDict convert_out(const std::vector<MX>& arg) const;
    std::vector<MX> convert_out(const MXDict& arg) const;
    /// @}

    /** \brief Does the function have free variables

        \identifier{1x8} */
    bool has_free() const;

    /** \brief Get free variables as a string

        \identifier{1x9} */
    std::vector<std::string> get_free() const;

    /** \brief Get all the free variables of the function

        \identifier{1xa} */
    std::vector<SX> free_sx() const;

    /** \brief Get all the free variables of the function

        \identifier{1xb} */
    std::vector<MX> free_mx() const;

    /** \brief Extract the functions needed for the Lifted Newton method

        \identifier{1xc} */
    void generate_lifted(Function& SWIG_OUTPUT(vdef_fcn),
                         Function& SWIG_OUTPUT(vinit_fcn)) const;

    /** \brief Number of nodes in the algorithm

        \identifier{1xd} */
    casadi_int n_nodes() const;

    /** \brief Number of instruction in the algorithm (SXFunction/MXFunction)

        \identifier{1xe} */
    casadi_int n_instructions() const;

    /** \brief Identifier index of the instruction (SXFunction/MXFunction)

        \identifier{1xf} */
    casadi_int instruction_id(casadi_int k) const;

    /** \brief Locations in the work vector for the inputs of the instruction

     * (SXFunction/MXFunction)

        \identifier{1xg} */
    std::vector<casadi_int> instruction_input(casadi_int k) const;

    /** \brief Get the floating point output argument of an instruction (SXFunction)

        \identifier{1xh} */
    double instruction_constant(casadi_int k) const;

    /** \brief Location in the work vector for the output of the instruction

     * (SXFunction/MXFunction)

        \identifier{1xi} */
    std::vector<casadi_int> instruction_output(casadi_int k) const;

    /** \brief Get the MX node corresponding to an instruction (MXFunction)

        \identifier{1xj} */
    MX instruction_MX(casadi_int k) const;

    /** \brief Get the SX node corresponding to all instructions (SXFunction)
     *
     * Note: input and output instructions have no SX representation.
     * This method returns nan for those instructions.

        \identifier{1xk} */
    SX instructions_sx() const;

    ///@{
    /** \brief  Is the class able to propagate seeds through the algorithm?

        \identifier{1xl} */
    bool has_spfwd() const;
    bool has_sprev() const;
    ///@}

    /** \brief Get required length of arg field

        \identifier{1xm} */
    size_t sz_arg() const;

    /** \brief Get required length of res field

        \identifier{1xn} */
    size_t sz_res() const;

    /** \brief Get required length of iw field

        \identifier{1xo} */
    size_t sz_iw() const;

    /** \brief Get required length of w field

        \identifier{1xp} */
    size_t sz_w() const;

#ifndef SWIG
    /** \brief Get number of temporary variables needed

        \identifier{1xq} */
    void sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const;

    /** \brief Set the (persistent) work vectors

        \identifier{1xr} */
    void set_work(const double**& arg, double**& res,
      casadi_int*& iw, double*& w, int mem=0) const;

    /** \brief Set the (temporary) work vectors

        \identifier{1xs} */
    void set_temp(const double** arg, double** res,
        casadi_int* iw, double* w, int mem=0) const;

    /** \brief Set the (persistent and temporary) work vectors

        \identifier{1xt} */
    void setup(const double** arg, double** res, casadi_int* iw, double* w, int mem=0) const;

    /** \brief Call using a map

        \identifier{1xu} */
    template<typename M>
    void call_gen(const std::map<std::string, M>& arg, std::map<std::string, M>& res,
               bool always_inline, bool never_inline) const;

    /** \brief List merge opportunitities

        \identifier{2b6} */
    void merge(const std::vector<MX>& arg,
        std::vector<MX>& subs_from, std::vector<MX>& subs_to) const;
#endif // SWIG
    /// \endcond

    /** \brief Name of the function

        \identifier{1xv} */
    const std::string& name() const;

    /** \brief Check if the function is of a particular type

        Optionally check if name matches one of the base classes (default true)

        \identifier{1xw} */
    bool is_a(const std::string& type, bool recursive=true) const;

    /** \brief Check if a string is a valid function name

     * Valid function names consist of
     *
     * - At least one character
     * - Upper and lower case letters: a-zA-Z
     * - Numbers 0-9, but never as first character
     * - Underscore, but never as first character and never next to another underscore
     *
     * May not be one of the following keywords: "null", "jac", "hess"

        \identifier{1xx} */
    static bool check_name(const std::string& name);

    /** \brief Turn a string into a valid function name as defined by "check_name"

     * Non-alphanumeric characters are converted into underscores and multiple
     * consecutive undercores are dropped

        \identifier{1xy} */
    static std::string fix_name(const std::string& name);

    /** \brief Build function from serialization

        \identifier{1xz} */
    static Function deserialize(std::istream& stream);

    /** \brief Build function from serialization

        \identifier{1y0} */
    static Function deserialize(const std::string& s);

    /** \brief Build function from serialization

        \identifier{1y1} */
    static Function load(const std::string& filename);

    /** \brief Build function from serialization

        \identifier{1y2} */
    static Function deserialize(DeserializingStream& s);

    /// Assert that an input dimension is equal so some given value
    void assert_size_in(casadi_int i, casadi_int nrow, casadi_int ncol) const;

    /// Assert that an output dimension is equal so some given value
    void assert_size_out(casadi_int i, casadi_int nrow, casadi_int ncol) const;

    /// Assert that an output sparsity is a multiple of some given sparsity
    void assert_sparsity_out(casadi_int i, const Sparsity& sp,
        casadi_int n = 1, bool allow_all_zero_sparse = true) const;

    /// Checkout a memory object
    casadi_int checkout() const;

    /// Release a memory object
    void release(int mem) const;

#ifndef SWIG
    /// Get memory object
    void* memory(int ind) const;

    static std::vector<SX> order(const std::vector<SX>& expr);
    static std::vector<MX> order(const std::vector<MX>& expr);
#endif // SWIG

    /** \brief Get all functions in the cache

        \identifier{26i} */
    Dict cache() const;

    /** \brief Get a list of all functions

        \identifier{1y3} */
    std::vector<std::string> get_function() const;

    /** \brief Get a dependency function

        \identifier{1y4} */
    Function get_function(const std::string &name) const;

    /** \brief Check if a particular dependency exists

        \identifier{1y5} */
    bool has_function(const std::string& fname) const;

    /** \brief Get all functions embedded in the expression graphs

      * \param[in] max_depth  Maximum depth - a negative number indicates no maximum

        \identifier{1y6} */
    std::vector<Function> find_functions(casadi_int max_depth = -1) const;

    /** \brief  Get a specific function embedded in the expression graphs

      * \param[in] name  Name of function needed
      * \param[in] max_depth  Maximum depth - a negative number indicates no maximum

        \identifier{1y7} */
    Function find_function(const std::string &name, casadi_int max_depth=-1) const;

    /** Obtain information about function */
    Dict info() const;

#ifndef SWIG
    protected:
    ///@{
    /** \brief Called by constructors

        \identifier{1y8} */
    void construct(const std::string& name,
                   const std::vector<SX>& ex_in, const std::vector<SX>& ex_out,
                   const std::vector<std::string>& name_in,
                   const std::vector<std::string>& name_out,
                   const Dict& opts);
    void construct(const std::string& name,
                   const std::vector<MX>& ex_in, const std::vector<MX>& ex_out,
                   const std::vector<std::string>& name_in,
                   const std::vector<std::string>& name_out,
                   const Dict& opts);
    template<typename M>
      void construct(const std::string& name, const std::map<std::string, M>& dict,
                     const std::vector<std::string>& name_in,
                     const std::vector<std::string>& name_out,
                     const Dict& opts);
    ///@}

    /// Helper function for parsing .casadi files
    static bool proceed_to(std::istream& file, const std::string& str);

    /// Helper function for mapaccum
    Function mapaccum(const std::string& name, const std::vector<Function>& chain,
                      casadi_int n_accum=1, const Dict& opts = Dict()) const;

#ifdef WITH_EXTRA_CHECKS
    public:
    // How many times have we passed through
    // operator()(const double** arg, double** res, casadi_int* iw, double* w, int mem)?
    static thread_local casadi_int call_depth_;
#endif


#endif // SWIG



  };


/** \brief Class to achieve minimal overhead function evaluations

    \identifier{1y9} */
class CASADI_EXPORT FunctionBuffer {
  Function f_;
  std::vector<double> w_;
  std::vector<casadi_int> iw_;
  std::vector<const double*> arg_;
  std::vector<double*> res_;
  FunctionInternal* f_node_;
  casadi_int mem_;
  void *mem_internal_;
  int ret_;
public:
  /** \brief Main constructor

      \identifier{1ya} */
  FunctionBuffer(const Function& f);
#ifndef SWIG
  ~FunctionBuffer();
  FunctionBuffer(const FunctionBuffer& f);
  FunctionBuffer& operator=(const FunctionBuffer& f);
#endif // SWIG

  /** \brief Set input buffer for input i

      mem.set_arg(0, memoryview(a))

      Note that CasADi uses 'fortran' order: column-by-column

      \identifier{1yb} */
  void set_arg(casadi_int i, const double* a, casadi_int size);

  /** \brief Set output buffer for ouput i

      mem.set_res(0, memoryview(a))

      Note that CasADi uses 'fortran' order: column-by-column

      \identifier{1yc} */
  void set_res(casadi_int i, double* a, casadi_int size);
  /// Get last return value
  int ret();
  void _eval();
  void* _self() { return this; }
  Dict stats() const;
};

void CASADI_EXPORT _function_buffer_eval(void* raw);


} // namespace casadi

#include "casadi_interrupt.hpp"
#include "runtime/shared.hpp"

#endif // CASADI_FUNCTION_HPP
