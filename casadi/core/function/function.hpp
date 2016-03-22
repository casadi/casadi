
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

#ifndef CASADI_FUNCTION_HPP
#define CASADI_FUNCTION_HPP

#include "../sx/sx_elem.hpp"
#include "../mx/mx.hpp"

#include <exception>

/*
  NOTE: The order of includes is as follows:
  1. Forward declaration of Function (in casadi_types.hpp)
  2. Declaration of Matrix class (in matrix.hpp), Function only as references or return values
  3. Definition of Function (this file), requires Matrix to be complete type
  4. Definition of Matrix class (in matrix_impl.hpp, included at the end of this file)
*/

namespace casadi {

#ifndef SWIG
  /** Forward declaration of internal class */
  class FunctionInternal;

#endif // SWIG

  /** \brief General function

      A general function \f$f\f$ in casadi can be multi-input, multi-output.\n
      Number of inputs:  \a nin    n_in()\n
      Number of outputs: \a nout   n_out()\n

      We can view this function as a being composed of a (\a nin, \a nout) grid of single-input,
      single-output primitive functions.\n
      Each such primitive function \f$f_ {i, j} \forall i \in [0, nin-1], j \in [0, nout-1]\f$ can
      map as \f$\mathbf {R}^{n, m}\to\mathbf{R}^{p, q}\f$,
      in which n, m, p, q can take different values for every (i, j) pair.\n

      When passing input, you specify which partition \f$i\f$ is active.
      You pass the numbers vectorized, as a vector of size \f$(n*m)\f$.\n
      When requesting output, you specify which partition \f$j\f$ is active.
      You get the numbers vectorized, as a vector of size \f$(p*q)\f$.\n

      To calculate Jacobians, you need to have \f$(m=1, q=1)\f$.

      Write the Jacobian as \f$J_ {i, j} = \nabla f_{i, j} =
      \frac {\partial f_{i, j}(\vec{x})}{\partial \vec{x}}\f$.

      We have the following relationships for function mapping from a row vector to a row vector:

      \f$ \vec {s}_f = \nabla f_{i, j} . \vec{v}\f$ \n
      \f$ \vec {s}_a = (\nabla f_{i, j})^T . \vec{w}\f$

      Some quantities in these formulas must be transposed: \n
      input  col: transpose \f$ \vec {v} \f$ and \f$\vec{s}_a\f$ \n
      output col: transpose \f$ \vec {w} \f$ and \f$\vec{s}_f\f$ \n

      NOTE: Functions are allowed to modify their input arguments when evaluating:
            implicitFunction, IDAS solver
      Further releases may disallow this.

      \internal
      \section Notes for developers

      Each function consists of 4 files:\n
      1. public class header file: imported in python\n
      2. public class implementation\n
      3. internal class header file: should only be used by derived classes\n
      4. internal class implementation\n

      python and c++ should be 1-to-1\n
      There should be no extra features in 1.\n
      All the functionality should exist in 1.\n
      If it means that c++ will be more "pythonic", so be it.
      \endinternal

      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT Function : public SharedObject {
  public:

    /** \brief Default constructor, null pointer */
    Function();

    /** \brief Construct from a file */
    Function(const std::string& fname);

    ///@{
    /** \brief Construct an SX function */
    Function(const std::string& name,
             const std::vector<SX>& arg, const std::vector<SX>& res,
             const Dict& opts=Dict());
    Function(const std::string& name,
             const std::vector<SX>& arg, const std::vector<SX>& res,
             const std::vector<std::string>& argn, const std::vector<std::string>& resn,
             const Dict& opts=Dict());
    Function(const std::string& name, const std::map<std::string, SX>& dict,
             const std::vector<std::string>& argn, const std::vector<std::string>& resn,
             const Dict& opts=Dict());
    ///@}

    ///@{
    /** \brief Construct an MX function */
    Function(const std::string& name,
             const std::vector<MX>& arg, const std::vector<MX>& res,
             const Dict& opts=Dict());
    Function(const std::string& name,
             const std::vector<MX>& arg, const std::vector<MX>& res,
             const std::vector<std::string>& argn, const std::vector<std::string>& resn,
             const Dict& opts=Dict());
    Function(const std::string& name, const std::map<std::string, MX>& dict,
             const std::vector<std::string>& argn, const std::vector<std::string>& resn,
             const Dict& opts=Dict());
    ///@}

    ///@{
    /** \brief To resolve ambiguity on some compilers */
#ifndef SWIG
    Function(const std::string& name, SXIList arg, const SXVector& res, const Dict& opts=Dict());
    Function(const std::string& name, const SXVector& arg, SXIList res, const Dict& opts=Dict());
    Function(const std::string& name, SXIList arg, SXIList res, const Dict& opts=Dict());
    Function(const std::string& name, SXIList arg, const SXVector& res,
             const StringVector& argn, const StringVector& resn, const Dict& opts=Dict());
    Function(const std::string& name, const SXVector& arg, SXIList res,
             const StringVector& argn, const StringVector& resn, const Dict& opts=Dict());
    Function(const std::string& name, SXIList arg, SXIList res,
             const StringVector& argn, const StringVector& resn, const Dict& opts=Dict());
    Function(const std::string& name, MXIList arg, const MXVector& res, const Dict& opts=Dict());
    Function(const std::string& name, const MXVector& arg, MXIList res, const Dict& opts=Dict());
    Function(const std::string& name, MXIList arg, MXIList res, const Dict& opts=Dict());
    Function(const std::string& name, MXIList arg, const MXVector& res,
             const StringVector& argn, const StringVector& resn, const Dict& opts=Dict());
    Function(const std::string& name, const MXVector& arg, MXIList res,
             const StringVector& argn, const StringVector& resn, const Dict& opts=Dict());
    Function(const std::string& name, MXIList arg, MXIList res,
             const StringVector& argn, const StringVector& resn, const Dict& opts=Dict());
#endif // SWIG
    ///@}

    /** \brief  Destructor */
    ~Function();

    /** \brief Expand a function to SX */
    ///@{
    Function expand() const;
    Function expand(const std::string& name,
                    const Dict& opts=Dict()) const;
    ///@}

    /// \cond INTERNAL
#ifndef SWIG
    /** \brief  Create from node */
    static Function create(FunctionInternal* node);
#endif // SWIG
    /// \endcond

    /** \brief Get the number of function inputs */
    int n_in() const;

    /** \brief Get the number of function outputs */
    int n_out() const;

    ///@{
    /** \brief Get input dimension */
    int size1_in(int ind) const;
    int size1_in(const std::string& iname) const {return size1_in(index_in(iname));}
    int size2_in(int ind) const;
    int size2_in(const std::string& iname) const {return size2_in(index_in(iname));}
    std::pair<int, int> size_in(int ind) const;
    std::pair<int, int> size_in(const std::string& iname) const {return size_in(index_in(iname));}
    ///@}

    ///@{
    /** \brief Get output dimension */
    int size1_out(int ind) const;
    int size1_out(const std::string& oname) const {return size1_out(index_out(oname));}
    int size2_out(int ind) const;
    int size2_out(const std::string& oname) const {return size2_out(index_out(oname));}
    std::pair<int, int> size_out(int ind) const;
    std::pair<int, int> size_out(const std::string& oname) const {
      return size_out(index_out(oname));
    }
    ///@}

    ///@{
    /** \brief  Get of number of input nonzeros
     * For a particular input or for all for all of the inputs
     */
    int nnz_in() const;
    int nnz_in(int ind) const;
    int nnz_in(const std::string& iname) const {return numel_in(index_in(iname));}
    ///@}

    ///@{
    /** \brief  Get of number of output nonzeros
     * For a particular output or for all for all of the outputs
     */
    int nnz_out() const;
    int nnz_out(int ind) const;
    int nnz_out(const std::string& oname) const {return numel_out(index_out(oname));}
    ///@}

    ///@{
    /** \brief  Get of number of input elements
     * For a particular input or for all for all of the inputs
     */
    int numel_in() const;
    int numel_in(int ind) const;
    int numel_in(const std::string& iname) const {return numel_in(index_in(iname));}
    ///@}

    ///@{
    /** \brief  Get of number of output elements
     * For a particular output or for all for all of the outputs
     */
    int numel_out() const;
    int numel_out(int ind) const;
    int numel_out(const std::string& oname) const {return numel_out(index_out(oname));}
    ///@}

    /** \brief Get input scheme */
    std::vector<std::string> name_in() const;

    /** \brief Get output scheme */
    std::vector<std::string> name_out() const;

    /** \brief Get input scheme name by index */
    std::string name_in(int ind) const;

    /** \brief Get output scheme name by index */
    std::string name_out(int ind) const;

    /** \brief Find the index for a string describing a particular entry of an input scheme
     *
     * example:  schemeEntry("x_opt")  -> returns  NLPSOL_X if FunctionInternal adheres to
     * SCHEME_NLPINput
     */
    int index_in(const std::string &name) const;

    /** \brief Find the index for a string describing a particular entry of an output scheme
     *
     * example:  schemeEntry("x_opt")  -> returns  NLPSOL_X if FunctionInternal adheres to
     * SCHEME_NLPINput
     */
    int index_out(const std::string &name) const;

    /** \brief Get default input value (NOTE: constant reference) */
    double default_in(int ind) const;

    /** \brief Get sparsity of a given input */
    /// @{
    const Sparsity& sparsity_in(int ind) const;
    const Sparsity& sparsity_in(const std::string& iname) const;
    /// @}

    /** \brief Get sparsity of a given output */
    /// @{
    const Sparsity& sparsity_out(int ind) const;
    const Sparsity& sparsity_out(const std::string& iname) const;
    /// @}

    /** \brief Print dimensions of inputs and outputs */
    void printDimensions(std::ostream &stream=casadi::userOut()) const;

    /** \brief Print options to a stream */
    void printOptions(std::ostream &stream=casadi::userOut()) const;

    /** \brief Print all information there is to know about a certain option */
    void printOption(const std::string &name, std::ostream &stream = casadi::userOut()) const;

    ///@{
    /** \brief Generate a Jacobian function of output \a oind with respect to input \a iind
     * \param iind The index of the input
     * \param oind The index of the output
     *
     * The default behavior of this class is defined by the derived class.
     * If compact is set to true, only the nonzeros of the input and output expressions are
     * considered.
     * If symmetric is set to true, the Jacobian being calculated is known to be symmetric
     * (usually a Hessian),
     * which can be exploited by the algorithm.
     *
     * The generated Jacobian has one more output than the calling function corresponding
     * to the Jacobian and the same number of inputs.
     *
     */
    Function jacobian(int iind=0, int oind=0, bool compact=false, bool symmetric=false);
    Function jacobian(const std::string& iind,  int oind=0, bool compact=false,
                      bool symmetric=false) {
        return jacobian(index_in(iind), oind, compact, symmetric);
    }
    Function jacobian(int iind, const std::string& oind, bool compact=false, bool symmetric=false) {
        return jacobian(iind, index_out(oind), compact, symmetric);
    }
    Function jacobian(const std::string& iind, const std::string& oind, bool compact=false,
                      bool symmetric=false) {
        return jacobian(index_in(iind), index_out(oind), compact, symmetric);
    }
    ///@}

    /** Set the Jacobian function of output \a oind with respect to input \a iind
     NOTE: Does _not_ take ownership, only weak references to the Jacobians are kept internally */
    void setJacobian(const Function& jac, int iind=0, int oind=0, bool compact=false);

    ///@{
    /** \brief Generate a gradient function of output \a oind with respect to input \a iind
     * \param iind The index of the input
     * \param oind The index of the output
     *
     * The default behavior of this class is defined by the derived class.
     * Note that the output must be scalar. In other cases, use the Jacobian instead.
     *
     */
    Function gradient(int iind=0, int oind=0);
    Function gradient(const std::string& iind, int oind=0) {
        return gradient(index_in(iind), oind);
    }
    Function gradient(int iind, const std::string& oind) {
        return gradient(iind, index_out(oind));
    }
    Function gradient(const std::string& iind, const std::string& oind) {
        return gradient(index_in(iind), index_out(oind));
    }
    ///@}

    ///@{
    /** \brief Generate a tangent function of output \a oind with respect to input \a iind
     * \param iind The index of the input
     * \param oind The index of the output
     *
     * The default behavior of this class is defined by the derived class.
     * Note that the input must be scalar. In other cases, use the Jacobian instead.
     *
     */
    Function tangent(int iind=0, int oind=0);
    Function tangent(const std::string& iind, int oind=0)
    { return tangent(index_in(iind), oind); }
    Function tangent(int iind, const std::string& oind)
    { return tangent(iind, index_out(oind)); }
    Function tangent(const std::string& iind, const std::string& oind)
    { return tangent(index_in(iind), index_out(oind)); }
    ///@}

    ///@{
    /** \brief Generate a Hessian function of output \a oind with respect to input \a iind
     * \param iind The index of the input
     * \param oind The index of the output
     *
     * The generated Hessian has two more outputs than the calling function corresponding
     * to the Hessian and the gradients.
     *
     */
    Function hessian(int iind=0, int oind=0);
    Function hessian(const std::string& iind, int oind=0)
    { return hessian(index_in(iind), oind); }
    Function hessian(int iind, const std::string& oind)
    { return hessian(iind, index_out(oind)); }
    Function hessian(const std::string& iind, const std::string& oind)
    { return hessian(index_in(iind), index_out(oind)); }
    ///@}

    /** \brief Generate a Jacobian function of all the inputs elements with respect to all
     * the output elements).
     */
    Function fullJacobian();

    /** Set the Jacobian of all the input nonzeros with respect to all output nonzeros
     NOTE: Does _not_ take ownership, only weak references to the Jacobian are kept internally */
    void setFullJacobian(const Function& jac);

    ///@{
    /** \brief Evaluate the function symbolically or numerically  */
    void call(const std::vector<DM> &arg, std::vector<DM>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false);
    void call(const std::vector<SX> &arg, std::vector<SX>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false);
    void call(const std::vector<MX> &arg, std::vector<MX>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false);
    void call(const DMDict& arg, DMDict& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false);
    void call(const SXDict& arg, SXDict& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false);
    void call(const MXDict& arg, MXDict& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false);
    ///@}

#ifndef SWIG
    ///@{
    /// Functor shorthand for evaluation
    std::vector<DM> operator()(const std::vector<DM>& arg);
    std::vector<SX> operator()(const std::vector<SX>& arg);
    std::vector<MX> operator()(const std::vector<MX>& arg);
    const DMDict operator()(const DMDict& arg);
    const SXDict operator()(const SXDict& arg);
    const MXDict operator()(const MXDict& arg);
    ///@}

    ///@{
    /** \brief Evaluate with temporary memory allocation */
    void operator()(std::vector<const double*> arg, std::vector<double*> res);
    void operator()(std::vector<const bvec_t*> arg, std::vector<bvec_t*> res);
    void operator()(std::vector<const SXElem*> arg, std::vector<SXElem*> res);
    template<typename D> void _call(std::vector<const D*> arg, std::vector<D*> res);
    ///@}

    ///@{
    /** \brief Supported arguments for numerical evaluation and converters */
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
    /** \brief Numerical evaluation */
    void operator()(VecArg arg, VecRes res) { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(VecArg arg, MapRes res) { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(VecArg arg, VPrRes res) { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(VecArg arg, MPrRes res) { (*this)(buf_in(arg), buf_out(res)); }

    void operator()(MapArg arg, VecRes res) { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(MapArg arg, MapRes res) { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(MapArg arg, VPrRes res) { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(MapArg arg, MPrRes res) { (*this)(buf_in(arg), buf_out(res)); }
    ///@}

    ///@{
    /// Functor shorthand for evaluation, single argument (only C++)
    std::vector<DM> operator()(const DM& arg0) {
      return operator()(std::vector<DM>{arg0});
    }
    std::vector<SX> operator()(const SX& arg0) {
      return operator()(std::vector<SX>{arg0});
    }
    std::vector<MX> operator()(const MX& arg0) {
      return operator()(std::vector<MX>{arg0});
    }
    ///@}

    /** \brief Evaluate memory-less, numerically */
    void operator()(const double** arg, double** res, int* iw, double* w, int mem=0) const;

    /** \brief Evaluate memory-less SXElem
        Same syntax as the double version, allowing use in templated code
     */
    void operator()(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem=0) const;

    /** \brief  Propagate sparsity forward */
    void operator()(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem=0) const;

    /** \brief  Propagate sparsity backward */
    void rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem=0);

    /** \brief Propagate sparsity backward with temporary memory allocation */
    void rev(std::vector<bvec_t*> arg, std::vector<bvec_t*> res);

#endif // SWIG

    /** \brief Create call to (cached) derivative function, forward mode  */
    void forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& fseed,
                 std::vector<std::vector<MX> >& SWIG_OUTPUT(fsens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    void reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& aseed,
                 std::vector<std::vector<MX> >& SWIG_OUTPUT(asens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, forward mode  */
    void forward(const std::vector<SX>& arg, const std::vector<SX>& res,
                 const std::vector<std::vector<SX> >& fseed,
                 std::vector<std::vector<SX> >& SWIG_OUTPUT(fsens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    void reverse(const std::vector<SX>& arg, const std::vector<SX>& res,
                 const std::vector<std::vector<SX> >& aseed,
                 std::vector<std::vector<SX> >& SWIG_OUTPUT(asens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, forward mode  */
    void forward(const std::vector<DM>& arg, const std::vector<DM>& res,
                 const std::vector<std::vector<DM> >& fseed,
                 std::vector<std::vector<DM> >& SWIG_OUTPUT(fsens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    void reverse(const std::vector<DM>& arg, const std::vector<DM>& res,
                 const std::vector<std::vector<DM> >& aseed,
                 std::vector<std::vector<DM> >& SWIG_OUTPUT(asens),
                 bool always_inline=false, bool never_inline=false);

    /// \cond INTERNAL
    ///@{
   /** \brief Evaluate the function symbolically or numerically with directional derivatives
     * The first two arguments are the nondifferentiated inputs and results of the evaluation,
     * the next two arguments are a set of forward directional seeds and the resulting forward
     * directional derivatives, the length of the vector being the number of forward directions.
     * The next two arguments are a set of adjoint directional seeds and the resulting adjoint
     * directional derivatives, the length of the vector being the number of adjoint directions.
     */
    void derivative(const DMVector& arg, DMVector& SWIG_OUTPUT(res),
                    const DMVectorVector& fseed, DMVectorVector& SWIG_OUTPUT(fsens),
                    const DMVectorVector& aseed, DMVectorVector& SWIG_OUTPUT(asens),
                    bool always_inline=false, bool never_inline=false);

    void derivative(const SXVector& arg, SXVector& SWIG_OUTPUT(res),
                    const SXVectorVector& fseed, SXVectorVector& SWIG_OUTPUT(fsens),
                    const SXVectorVector& aseed, SXVectorVector& SWIG_OUTPUT(asens),
                    bool always_inline=false, bool never_inline=false);

    void derivative(const MXVector& arg, MXVector& SWIG_OUTPUT(res),
                    const MXVectorVector& fseed, MXVectorVector& SWIG_OUTPUT(fsens),
                    const MXVectorVector& aseed, MXVectorVector& SWIG_OUTPUT(asens),
                    bool always_inline=false, bool never_inline=false);
    ///@}
    /// \endcond

    /** \brief  Evaluate symbolically in parallel (matrix graph)
        \param parallelization Type of parallelization used: unroll|serial|openmp
    */
    std::vector<MX> map(const std::vector<MX > &arg,
                        const std::string& parallelization="serial");

    /** \brief  Evaluate symbolically in parallel (matrix graph)
        \param parallelization Type of parallelization used: unroll|serial|openmp
    */
    std::map<std::string, MX> map(const std::map<std::string, MX> &arg,
                        const std::string& parallelization="serial");

    /** \brief  Evaluate symbolically in parallel and sum (matrix graph)
        \param parallelization Type of parallelization used: unroll|serial|openmp
    */
    std::vector<MX> mapsum(const std::vector<MX > &arg,
                           const std::string& parallelization="serial");

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


    */
    Function mapaccum(const std::string& name, int n, const Dict& opts = Dict());
    Function mapaccum(const std::string& name, int n,
                      const std::vector<int>& accum_in,
                      const std::vector<int>& accum_out,
                      const Dict& opts=Dict());
    ///@}


    ///@{
    /** \brief  Create a mapped version of this function

        Suppose the function has a signature of:
        \verbatim
           f: (a, p) -> ( s )
        \endverbatim

        The the mapaccumulated version has the signature:
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

    */

    Function map(const std::string& name, const std::string& parallelization, int n,
      const std::vector<int>& reduce_in, const std::vector<int>& reduce_out,
      const Dict& opts=Dict());
    Function map(const std::string& name, const std::string& parallelization, int n,
      const Dict& opts=Dict());
    ///@}

    /** \brief returns a new function with a selection of inputs/outputs of the original */
    Function slice(const std::vector<int>& order_in, const std::vector<int>& order_out,
                   const Dict& opts=Dict());

    /** \brief Constuct a switch function */
    static Function conditional(const std::string& name, const std::vector<Function>& f,
                                const Function& f_def, const Dict& opts=Dict());

    /** \brief Constructor (if-else) */
    static Function if_else(const std::string& name, const Function& f_true,
                            const Function& f_false, const Dict& opts=Dict());

    /** kernel_sum
        Consider a dense matrix V.

        KernelSum computes

        F(V,X)  = sum_i sum_j  f ( [i;j], V(i,j), X)

        with X: [x;y]

        where the summation is taken for all entries (i,j)
        that are a distance r away from X.

        This function assumes that V is fixed:
        sensitivities with respect to it are not computed.

        This allows for improved speed of evaluation.

        Having V fixed is a common use case:
        V may be a large bitmap (observation),
        onto which a kernel is fitted.

        \author Joris Gillis
        \date 2015
    */
    Function kernel_sum(const std::string& name,
                        const std::pair<int, int> & size,
                        double r, int n,
                        const Dict& opts=Dict()) const;

    /** \brief Get a function that calculates \a nfwd forward derivatives and nadj adjoint derivatives
     *         Legacy function: Use forward and reverse instead.
     *
     *         Returns a function with <tt>(1+nfwd)*n_in+nadj*n_out</tt> inputs
     *         and <tt>(1+nfwd)*n_out + nadj*n_in</tt> outputs.
     *         The first <tt>n_in</tt> inputs correspond to nondifferentiated inputs.
     *         The next <tt>nfwd*n_in</tt> inputs correspond to forward seeds,
     *         one direction at a time
     *         and the last <tt>nadj*n_out</tt> inputs correspond to adjoint seeds,
     *         one direction at a time.
     *         The first n_out outputs correspond to nondifferentiated outputs.
     *         The next <tt>nfwd*n_out</tt> outputs correspond to forward sensitivities,
     *         one direction at a time and the last <tt>nadj*n_in</tt> outputs corresponds to
     *         adjoint sensitivities, one direction at a time.
     *
     *         <tt>(n_in = n_in(), n_out = n_out())</tt>
     *
     */
    Function derivative(int nfwd, int nadj);

    /** \brief Get a function that calculates \a nfwd forward derivatives
     *
     *         Returns a function with <tt>n_in + n_out +nfwd*n_in</tt> inputs
     *         and <tt>nfwd*n_out</tt> outputs.
     *         The first <tt>n_in</tt> inputs correspond to nondifferentiated inputs.
     *         The next <tt>n_out</tt> inputs correspond to nondifferentiated outputs.
     *         and the last <tt>nfwd*n_in</tt> inputs correspond to forward seeds,
     *         one direction at a time
     *         The  <tt>nfwd*n_out</tt> outputs correspond to forward sensitivities,
     *         one direction at a time.     *
     *         <tt>(n_in = n_in(), n_out = n_out())</tt>
     *
     *        The functions returned are cached, meaning that if called multiple timed
     *        with the same value, then multiple references to the same function will be returned.
     */
    Function forward(int nfwd);

    /** \brief Get a function that calculates \a nadj adjoint derivatives
     *
     *         Returns a function with <tt>n_in + n_out +nadj*n_out</tt> inputs
     *         and <tt>nadj*n_in</tt> outputs.
     *         The first <tt>n_in</tt> inputs correspond to nondifferentiated inputs.
     *         The next <tt>n_out</tt> inputs correspond to nondifferentiated outputs.
     *         and the last <tt>nadj*n_out</tt> inputs correspond to adjoint seeds,
     *         one direction at a time
     *         The  <tt>nadj*n_in</tt> outputs correspond to adjoint sensitivities,
     *         one direction at a time.     *
     *         <tt>(n_in = n_in(), n_out = n_out())</tt>
     *
     *         <tt>(n_in = n_in(), n_out = n_out())</tt>
     *
     *        The functions returned are cached, meaning that if called multiple timed
     *        with the same value, then multiple references to the same function will be returned.
     */
    Function reverse(int nadj);

    /** \brief Set a function that calculates \a nfwd forward derivatives
        NOTE: Does _not_ take ownership, only weak references to the derivatives are kept internally */
    void set_forward(const Function& fcn, int nfwd);

    /** \brief Set a function that calculates \a nadj adjoint derivatives
        NOTE: Does _not_ take ownership, only weak references to the derivatives are kept internally */
    void set_reverse(const Function& fcn, int nadj);

    ///@{
    /// Get, if necessary generate, the sparsity of a Jacobian block
    const Sparsity sparsity_jac(int iind=0, int oind=0, bool compact=false, bool symmetric=false);
    const Sparsity sparsity_jac(const std::string &iind, int oind=0, bool compact=false,
                                bool symmetric=false) {
      return sparsity_jac(index_in(iind), oind, compact, symmetric);
    }
    const Sparsity sparsity_jac(int iind, const std::string &oind, bool compact=false,
                                bool symmetric=false) {
      return sparsity_jac(iind, index_out(oind), compact, symmetric);
    }
    const Sparsity sparsity_jac(const std::string &iind, const std::string &oind,
                          bool compact=false, bool symmetric=false) {
      return sparsity_jac(index_in(iind), index_out(oind), compact, symmetric);
    }
    ///@}

    ///@{
    /// Generate the sparsity of a Jacobian block
    void set_jac_sparsity(const Sparsity& sp, int iind, int oind, bool compact=false);
    void set_jac_sparsity(const Sparsity& sp, const std::string &iind, int oind,
                          bool compact=false) {
      set_jac_sparsity(sp, index_in(iind), oind, compact);
    }
    void set_jac_sparsity(const Sparsity& sp, int iind, const std::string &oind,
                          bool compact=false) {
      set_jac_sparsity(sp, iind, index_out(oind), compact);
    }
    void set_jac_sparsity(const Sparsity& sp, const std::string &iind, const std::string &oind,
                          bool compact=false) {
      set_jac_sparsity(sp, index_in(iind), index_out(oind), compact);
    }
    ///@}

    /** \brief Export / Generate C code for the function */
    void generate(const std::string& fname, const Dict& opts=Dict());

    /** \brief Export / Generate C code for the function */
    void generate(const Dict& opts=Dict());

    /** \brief Export / Generate C code for the dependency function */
    void generate_dependencies(const std::string& fname, const Dict& opts=Dict());

#ifndef SWIG
    /// \cond INTERNAL
    /// Get a const pointer to the node
    FunctionInternal* get() const;

    /** \brief  Const access functions of the node */
    FunctionInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectNode* ptr);
    /// \endcond
#endif // SWIG

    /// Get all statistics obtained at the end of the last evaluate call
    Dict stats(int mem=0) const;

    ///@{
    /** \brief Get symbolic primitives equivalent to the input expressions
     * There is no guarantee that subsequent calls return unique answers
     */
    const SX sx_in(int iind) const;
    const SX sx_in(const std::string& iname) const {
      return sx_in(index_in(iname));
    }
    const std::vector<SX> sx_in() const;
    const MX mx_in(int ind) const;
    const MX mx_in(const std::string & iname) const {
      return mx_in(index_in(iname));
    }
    const std::vector<MX> mx_in() const;
    ///@}

    ///@{
    /** \brief Get symbolic primitives equivalent to the output expressions
    * There is no guarantee that subsequent calls return unique answers
    */
    const SX sx_out(int oind) const;
    const SX sx_out(const std::string& oname) const {
      return sx_out(index_out(oname));
    }
    const std::vector<SX> sx_out() const;
    const MX mx_out(int ind) const;
    const MX mx_out(const std::string& oname) const {
      return mx_out(index_out(oname));
    }
    const std::vector<MX> mx_out() const;
    ///@}

    /** \brief Get all the free variables of the function */
    SX free_sx() const;

    /** \brief Get all the free variables of the function */
    std::vector<MX> free_mx() const;

    /** \brief Does the function have free variables */
    bool has_free() const;

    /** \brief Extract the functions needed for the Lifted Newton method */
    void generate_lifted(Function& SWIG_OUTPUT(vdef_fcn),
                                  Function& SWIG_OUTPUT(vinit_fcn));

    /** \brief Get the number of atomic operations */
    int getAlgorithmSize() const;

    /** \brief Get the length of the work vector */
    int getWorkSize() const;

    /** \brief Get an atomic operation operator index */
    int getAtomicOperation(int k) const;

    /** \brief Get the (integer) input arguments of an atomic operation */
    std::pair<int, int> getAtomicInput(int k) const;

    /** \brief Get the floating point output argument of an atomic operation */
    double getAtomicInputReal(int k) const;

    /** \brief Get the (integer) output argument of an atomic operation */
    int getAtomicOutput(int k) const;

    /** \brief Number of nodes in the algorithm */
    int n_nodes() const;

    /// \cond INTERNAL
    /** \brief Is the class able to propagate seeds through the algorithm?
     *
     * (for usage, see the example propagating_sparsity.cpp) */
    bool spCanEvaluate(bool fwd);

    /** \brief Get required length of arg field */
    size_t sz_arg() const;

    /** \brief Get required length of res field */
    size_t sz_res() const;

    /** \brief Get required length of iw field */
    size_t sz_iw() const;

    /** \brief Get required length of w field */
    size_t sz_w() const;

#ifndef SWIG
    /** \brief Get number of temporary variables needed */
    void sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const;

    /** \brief Set the (persistent) work vectors */
    void set_work(const double**& arg, double**& res, int*& iw, double*& w, int mem=0) const;

    /** \brief Set the (temporary) work vectors */
    void set_temp(const double** arg, double** res, int* iw, double* w, int mem=0) const;

    /** \brief Set the (persistent and temporary) work vectors */
    void setup(const double** arg, double** res, int* iw, double* w, int mem=0) const;
#endif // SWIG

    /// \endcond

    /** \brief Add modules to be monitored */
    void addMonitor(const std::string& mon);

    /** \brief Remove modules to be monitored */
    void removeMonitor(const std::string& mon);

    /// \cond INTERNAL
    /** \brief Check if the numerical values of the supplied bounds make sense */
    void checkInputs() const;
#ifndef SWIG
    /** \brief Call using a map */
    template<typename M>
    void _call(const std::map<std::string, M>& arg, std::map<std::string, M>& res,
               bool always_inline, bool never_inline);
#endif // SWIG
    /// \endcond

    /** \brief Name of the function */
    std::string name() const;

    /** \brief Get type name */
    std::string type_name() const;

    /** \brief Check if the function is of a particular type
        Optionally check if name matches one of the base classes (default true)
     */
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
     */
    static bool check_name(const std::string& name);

    /** \brief Turn a string into a valid function name as defined by "check_name"
     * Non-alphanumeric characters are converted into underscores and multiple
     * consecutive undercores are dropped
     */
    static std::string fix_name(const std::string& name);

    /// Checkout a memory object
    int checkout();

    /// Release a memory object
    void release(int mem);

    /// Create a solve node
    MX linsol_solve(const MX& A, const MX& B, bool tr=false);

#ifndef SWIG
    /// Get memory object
    void* memory(int ind) const;

    // Factorize linear system of equations
    void linsol_factorize(const double* A, int mem=0) const;

    // Solve factorized linear system of equations
    void linsol_solve(double* x, int nrhs=1, bool tr=false, int mem=0) const;

    ///@{
    /// Propagate sparsity through a linear solve
    void linsol_spsolve(bvec_t* X, const bvec_t* B, bool tr=false) const;
    void linsol_spsolve(DM& X, const DM& B, bool tr=false) const;
    ///@}

    /** \brief Solve the system of equations <tt>Lx = b</tt>
        Only when a Cholesky factorization is available
    */
    void linsol_solveL(double* x, int nrhs, bool tr, int mem=0) const;
#endif // SWIG

    /** \brief Obtain a symbolic Cholesky factorization
        Only for Cholesky solvers
    */
    Sparsity linsol_cholesky_sparsity(bool tr=false, int mem=0) const;

    /** \brief Obtain a numeric Cholesky factorization
        Only for Cholesky solvers
     */
    DM linsol_cholesky(bool tr=false, int mem=0) const;

    /// Access rhs function for a rootfinder
    Function rootfinder_fun();

    /// Access Jacobian of the ths function for a rootfinder
    Function rootfinder_jac();

    /// Access linear solver of a rootfinder
    Function rootfinder_linsol();

    /// Get the DAE for an integrator
    Function integrator_dae();

    /** Generate native code in the interfaced language for debugging */
    void qpsol_debug(const std::string &filename) const;

    /** Generate native code in the interfaced language for debugging */
    void qpsol_debug(std::ostream &file) const;

#ifndef SWIG
    protected:
    ///@{
    /** \brief Called by constructors */
    void construct(const std::string& name,
                   const std::vector<SX>& arg, const std::vector<SX>& res,
                   const Dict& opts);
    void construct(const std::string& name,
                   const std::vector<MX>& arg, const std::vector<MX>& res,
                   const Dict& opts);
    template<typename M>
      void construct(const std::string& name,
                     const std::vector<M>& arg, const std::vector<M>& res,
                     const std::vector<std::string>& argn, const std::vector<std::string>& resn,
                     const Dict& opts);
    template<typename M>
      void construct(const std::string& name, const std::map<std::string, M>& dict,
                     const std::vector<std::string>& argn,
                     const std::vector<std::string>& resn,
                     const Dict& opts);
    ///@}

    /// Helper function for parsing .casadi files
    static bool proceed_to(std::istream& file, const std::string& str);
#endif // SWIG
  };

} // namespace casadi

#include "../matrix_impl.hpp"

#endif // CASADI_FUNCTION_HPP
