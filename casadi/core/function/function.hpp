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

#include "io_interface.hpp"

namespace casadi {

  /** Forward declaration of internal class */
  class FunctionInternal;

  /** \brief General function

      A general function \f$f\f$ in casadi can be multi-input, multi-output.\n
      Number of inputs:  \a nin    nIn()\n
      Number of outputs: \a nout   nOut()\n

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
  class CASADI_EXPORT Function : public OptionsFunctionality, public IOInterface<Function>{
  public:

    /// \cond CLUTTER
    /** \brief  default constructor */
    Function();

    /** \brief  Destructor */
    ~Function();
    /// \endcond

    /// \cond INTERNAL
#ifndef SWIG
    /** \brief  Create from node */
    static Function create(FunctionInternal* node);
#endif // SWIG
    /// \endcond

    /** \brief Find the index for a string describing a particular entry of an input scheme
     *
     * example:  schemeEntry("x_opt")  -> returns  NLP_SOLVER_X if FunctionInternal adheres to
     * SCHEME_NLPINput
     */
    int inputIndex(const std::string &name) const;

    /** \brief Find the index for a string describing a particular entry of an output scheme
     *
     * example:  schemeEntry("x_opt")  -> returns  NLP_SOLVER_X if FunctionInternal adheres to
     * SCHEME_NLPINput
     */
    int outputIndex(const std::string &name) const;

    /** \brief Get input scheme name by index */
    std::string inputName(int ind) const;

    /** \brief Get output scheme name by index */
    std::string outputName(int ind) const;

    /** \brief Get input scheme description by index */
    std::string inputDescription(int ind) const;

    /** \brief Get output scheme description by index */
    std::string outputDescription(int ind) const;

    /** \brief Get default input value */
    double defaultInput(int ind) const;

    /** \brief Get sparsity of a given input */
    /// @{
    Sparsity inputSparsity(int ind=0) const;
    Sparsity inputSparsity(const std::string& iname) const;
    /// @}

    /** \brief Get sparsity of a given output */
    /// @{
    Sparsity outputSparsity(int ind=0) const;
    Sparsity outputSparsity(const std::string& iname) const;
    /// @}

#ifndef SWIG
    /// \cond UNSAFE
    /** \brief [UNSAFE] Obtain reference to inputs
     * \sa getInput, setInput
     */
    ///@{
    /// Access input argument
    const Matrix<double>& input(int i=0) const;
    const Matrix<double>& input(const std::string &iname) const;
    Matrix<double>& input(int i=0);
    Matrix<double>& input(const std::string &iname);
    ///@}

    /** \brief [UNSAFE] Obtain reference to outputs
     * \sa getOutput, getOutput
     */
    ///@{
    /// Access output argument
    const Matrix<double>& output(int i=0) const;
    const Matrix<double>& output(const std::string &oname) const;
    Matrix<double>& output(int i=0);
    Matrix<double>& output(const std::string &oname);
    ///@}
    /// \endcond
#endif

    /** \brief Get the number of function inputs */
    int nIn() const;

    /** \brief Get the number of function outputs */
    int nOut() const;

    /** \brief  Get total number of nonzeros in all of the matrix-valued inputs */
    int nnzIn() const;

    /** \brief  Get total number of nonzeros in all of the matrix-valued outputs */
    int nnzOut() const;

    /** \brief  Get total number of elements in all of the matrix-valued inputs */
    int numelIn() const;

    /** \brief  Get total number of elements in all of the matrix-valued outputs */
    int numelOut() const;

    /** \brief Get input scheme */
    std::vector<std::string> inputScheme() const;

    /** \brief Get output scheme */
    std::vector<std::string> outputScheme() const;

    /** \brief  Evaluate */
    void evaluate();

    /** \brief  Print dimensions of inputs and outputs */
    void printDimensions(std::ostream &stream=casadi::userOut()) const;

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
        return jacobian(inputIndex(iind), oind, compact, symmetric);
    }
    Function jacobian(int iind, const std::string& oind, bool compact=false, bool symmetric=false) {
        return jacobian(iind, outputIndex(oind), compact, symmetric);
    }
    Function jacobian(const std::string& iind, const std::string& oind, bool compact=false,
                      bool symmetric=false) {
        return jacobian(inputIndex(iind), outputIndex(oind), compact, symmetric);
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
        return gradient(inputIndex(iind), oind);
    }
    Function gradient(int iind, const std::string& oind) {
        return gradient(iind, outputIndex(oind));
    }
    Function gradient(const std::string& iind, const std::string& oind) {
        return gradient(inputIndex(iind), outputIndex(oind));
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
    { return tangent(inputIndex(iind), oind); }
    Function tangent(int iind, const std::string& oind)
    { return tangent(iind, outputIndex(oind)); }
    Function tangent(const std::string& iind, const std::string& oind)
    { return tangent(inputIndex(iind), outputIndex(oind)); }
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
    { return hessian(inputIndex(iind), oind); }
    Function hessian(int iind, const std::string& oind)
    { return hessian(iind, outputIndex(oind)); }
    Function hessian(const std::string& iind, const std::string& oind)
    { return hessian(inputIndex(iind), outputIndex(oind)); }
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
    void call(const std::vector<DMatrix> &arg, std::vector<DMatrix>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false);
    void call(const std::vector<SX> &arg, std::vector<SX>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false);
    void call(const std::vector<MX> &arg, std::vector<MX>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false);
    ///@}

    ///@{
    /// Functor shorthand for evaluation
    std::vector<DMatrix> operator()(const std::vector<DMatrix>& arg,
                                    bool always_inline=false, bool never_inline=false);
    std::vector<SX> operator()(const std::vector<SX>& arg,
                               bool always_inline=false, bool never_inline=false);
    std::vector<MX> operator()(const std::vector<MX>& arg,
                               bool always_inline=false, bool never_inline=false);
    const DMatrixDict operator()(const DMatrixDict& arg,
                                 bool always_inline=false, bool never_inline=false);
    const SXDict operator()(const SXDict& arg,
                            bool always_inline=false, bool never_inline=false);
    const MXDict operator()(const MXDict& arg,
                            bool always_inline=false, bool never_inline=false);
    ///@}

#ifndef SWIG
    ///@{
    /// Functor shorthand for evaluation, single argument (only C++)
    std::vector<DMatrix> operator()(const DMatrix& arg0) { return operator()(make_vector(arg0));}
    std::vector<SX> operator()(const SX& arg0) { return operator()(make_vector(arg0));}
    std::vector<MX> operator()(const MX& arg0) { return operator()(make_vector(arg0));}
    ///@}

    /** \brief Evaluate memory-less, numerically */
    void operator()(const double** arg, double** res, int* iw, double* w);

    /** \brief Evaluate memory-less SXElement
        Same syntax as the double version, allowing use in templated code
     */
    void operator()(const SXElement** arg, SXElement** res, int* iw, SXElement* w);

    /** \brief  Propagate sparsity forward */
    void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief Alias for spFwd
        Same syntax as the double and SXElement versions, allowing use in templated code
    */
    void operator()(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
      spFwd(arg, res, iw, w);
    }
#endif // SWIG

    /** \brief Create call to (cached) derivative function, forward mode  */
    void callForward(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& fseed,
                 std::vector<std::vector<MX> >& SWIG_OUTPUT(fsens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    void callReverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& aseed,
                 std::vector<std::vector<MX> >& SWIG_OUTPUT(asens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, forward mode  */
    void callForward(const std::vector<SX>& arg, const std::vector<SX>& res,
                 const std::vector<std::vector<SX> >& fseed,
                 std::vector<std::vector<SX> >& SWIG_OUTPUT(fsens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    void callReverse(const std::vector<SX>& arg, const std::vector<SX>& res,
                 const std::vector<std::vector<SX> >& aseed,
                 std::vector<std::vector<SX> >& SWIG_OUTPUT(asens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, forward mode  */
    void callForward(const std::vector<DMatrix>& arg, const std::vector<DMatrix>& res,
                 const std::vector<std::vector<DMatrix> >& fseed,
                 std::vector<std::vector<DMatrix> >& SWIG_OUTPUT(fsens),
                 bool always_inline=false, bool never_inline=false);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    void callReverse(const std::vector<DMatrix>& arg, const std::vector<DMatrix>& res,
                 const std::vector<std::vector<DMatrix> >& aseed,
                 std::vector<std::vector<DMatrix> >& SWIG_OUTPUT(asens),
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
    void callDerivative(const DMatrixVector& arg, DMatrixVector& SWIG_OUTPUT(res),
                        const DMatrixVectorVector& fseed, DMatrixVectorVector& SWIG_OUTPUT(fsens),
                        const DMatrixVectorVector& aseed, DMatrixVectorVector& SWIG_OUTPUT(asens),
                        bool always_inline=false, bool never_inline=false);

    void callDerivative(const SXVector& arg, SXVector& SWIG_OUTPUT(res),
                        const SXVectorVector& fseed, SXVectorVector& SWIG_OUTPUT(fsens),
                        const SXVectorVector& aseed, SXVectorVector& SWIG_OUTPUT(asens),
                        bool always_inline=false, bool never_inline=false);

    void callDerivative(const MXVector& arg, MXVector& SWIG_OUTPUT(res),
                        const MXVectorVector& fseed, MXVectorVector& SWIG_OUTPUT(fsens),
                        const MXVectorVector& aseed, MXVectorVector& SWIG_OUTPUT(asens),
                        bool always_inline=false, bool never_inline=false);
    ///@}
    /// \endcond

    /** \brief  Evaluate symbolically in parallel (matrix graph)
        \param parallelization Type of parallelization used: expand|serial|openmp
    */
    std::vector<std::vector<MX> > map(const std::vector<std::vector<MX> > &arg,
                                      const std::string& parallelization="serial");

    /** \brief  Evaluate symbolically in parallel (matrix graph)
        \param parallelization Type of parallelization used: expand|serial|openmp
    */
    std::vector<MX> map(const std::vector<MX > &arg,
                                      const std::string& parallelization="serial");

    /** \brief  Evaluate symbolically in parallel and sum (matrix graph)
        \param parallelization Type of parallelization used: expand|serial|openmp
    */
    std::vector<MX> mapsum(const std::vector<MX > &arg,
                                      const std::string& parallelization="serial");

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
    Function mapaccum(const std::string& name, int N, const Dict & options = Dict()) const;


    /** \brief  Create a mapped version of this function

        Suppose the function has a signature of:
        \verbatim
           f: (a, p) -> ( s )
        \endverbatim

        The the mapaccumulated version has the signature:
        \verbatim
           F: (A, P) -> (S )

            with
                a: horzcat([a0, a1, ..., a_(N-1)])
                p: horzcat([p0, p1, ..., p_(N-1)])
                s: horzcat([s0, s1, ..., s_(N-1)])
            and
                s0 <- f(a0, p0)
                s1 <- f(a1, p1)
                ...
                s_(N-1) <- f(a_(N-1), p_(N-1))
        \endverbatim


    */
    Function map(const std::string& name, int N,  const Dict & options = Dict()) const;

    /** \brief Get a function that calculates \a nfwd forward derivatives and nadj adjoint derivatives
     *         Legacy function: Use derForward and derReverse instead.
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
     *         <tt>(n_in = nIn(), n_out = nOut())</tt>
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
     *         <tt>(n_in = nIn(), n_out = nOut())</tt>
     *
     *        The functions returned are cached, meaning that if called multiple timed
     *        with the same value, then multiple references to the same function will be returned.
     */
    Function derForward(int nfwd);

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
     *         <tt>(n_in = nIn(), n_out = nOut())</tt>
     *
     *         <tt>(n_in = nIn(), n_out = nOut())</tt>
     *
     *        The functions returned are cached, meaning that if called multiple timed
     *        with the same value, then multiple references to the same function will be returned.
     */
    Function derReverse(int nadj);

    /** \brief Set a function that calculates \a nfwd forward derivatives
        NOTE: Does _not_ take ownership, only weak references to the derivatives are kept internally */
    void setDerForward(const Function& fcn, int nfwd);

    /** \brief Set a function that calculates \a nadj adjoint derivatives
        NOTE: Does _not_ take ownership, only weak references to the derivatives are kept internally */
    void setDerReverse(const Function& fcn, int nadj);

    ///@{
    /// Get, if necessary generate, the sparsity of a Jacobian block
    const Sparsity jacSparsity(int iind=0, int oind=0, bool compact=false, bool symmetric=false);
    const Sparsity jacSparsity(const std::string &iind, int oind=0, bool compact=false,
                          bool symmetric=false) {
        return jacSparsity(inputIndex(iind), oind, compact, symmetric); }
    const Sparsity jacSparsity(int iind, const std::string &oind, bool compact=false,
                          bool symmetric=false) {
        return jacSparsity(iind, outputIndex(oind), compact, symmetric); }
    const Sparsity jacSparsity(const std::string &iind, const std::string &oind,
                          bool compact=false, bool symmetric=false) {
        return jacSparsity(inputIndex(iind), outputIndex(oind), compact, symmetric); }
    ///@}

    ///@{
    /// Generate the sparsity of a Jacobian block
    void setJacSparsity(const Sparsity& sp, int iind, int oind, bool compact=false);
    void setJacSparsity(const Sparsity& sp, const std::string &iind, int oind, bool compact=false) {
        setJacSparsity(sp, inputIndex(iind), oind, compact); }
    void setJacSparsity(const Sparsity& sp, int iind, const std::string &oind, bool compact=false) {
        setJacSparsity(sp, iind, outputIndex(oind), compact); }
    void setJacSparsity(const Sparsity& sp, const std::string &iind, const std::string &oind,
                        bool compact=false) {
        setJacSparsity(sp, inputIndex(iind), outputIndex(oind), compact); }
    ///@}

    /** \brief Export / Generate C code for the function */
    void generate(const std::string& fname, const Dict& opts=Dict());

    /** \brief Export / Generate C code for the function */
    void generate(const Dict& opts=Dict());

    /// \cond INTERNAL
    /** \brief  Access functions of the node */
    FunctionInternal* operator->();

    /** \brief  Const access functions of the node */
    const FunctionInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
    /// \endcond

    /// Get all statistics obtained at the end of the last evaluate call
    const Dict& getStats() const;

    /// Get a single statistic obtained at the end of the last evaluate call
    GenericType getStat(const std::string& name) const;

    /** \brief  Get a vector of symbolic variables with the same dimensions as the inputs
     *
     * There is no guarantee that consecutive calls return identical objects
     */
    std::vector<MX> symbolicInput(bool unique=false) const;

    /** \brief Get a vector of symbolic variables with the same dimensions as the inputs, SX graph
     *
     * There is no guarantee that consecutive calls return identical objects
     */
    std::vector<SX> symbolicInputSX() const;

    /** \brief  Get a vector of symbolic variables with the same dimensions as the outputs
     *
     * There is no guarantee that consecutive calls return identical objects
     */
    std::vector<MX> symbolicOutput() const;

    /// \cond INTERNAL
    /** \brief Is the class able to propagate seeds through the algorithm?
     *
     * (for usage, see the example propagating_sparsity.cpp) */
    bool spCanEvaluate(bool fwd);

    /** \brief Reset the sparsity propagation
     *
     * (for usage, see the example propagating_sparsity.cpp) */
    void spInit(bool fwd);

    /** \brief Propagate the sparsity pattern through a set of directional
     *
     * derivatives forward or backward (for usage, see the example propagating_sparsity.cpp) */
    void spEvaluate(bool fwd);

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
    const std::map<std::string, M> callMap(const std::map<std::string, M>& arg,
                                       bool always_inline, bool never_inline);

    /** \brief Check if input arguments have correct length and dimensions
    *
    * \param hcat check if horizontal repetion of the function input is allowed
    */
    template<typename M>
    void checkArg(const std::vector<M>& arg, bool hcat=false) const;

    /** \brief Check if output arguments have correct length and dimensions */
    template<typename M>
    void checkRes(const std::vector<M>& res) const;

    /** \brief Check forward mode seeds dimensions */
    template<typename M>
    void checkFwdSeed(const std::vector<std::vector<M> >& fseed) const;

    /** \brief Check reverse mode seeds dimensions */
    template<typename M>
    void checkAdjSeed(const std::vector<std::vector<M> >& aseed) const;

    /** \brief Check if input arguments that needs to be replaced
    * 
    * \param hcat check if horizontal repetion of the function input is allowed
    */
    template<typename M>
    bool matchingArg(const std::vector<M>& arg, bool hcat=false) const;

    /** \brief Check if output arguments that needs to be replaced */
    template<typename M>
    bool matchingRes(const std::vector<M>& arg) const;

    /** \brief Check if there are 0-by-0 forward seeds that needs to be replaced */
    template<typename M>
    bool matchingFwdSeed(const std::vector<std::vector<M> >& fseed) const;

    /** \brief Check if there are 0-by-0 reverse seeds that needs to be replaced */
    template<typename M>
    bool matchingAdjSeed(const std::vector<std::vector<M> >& aseed) const;

    /** \brief Replace 0-by-0 inputs
    *
    * \param hcat check if horizontal repetion of the function input is allowed
    */
    template<typename M>
    std::vector<M> replaceArg(const std::vector<M>& arg, bool hcat=false) const;

    /** \brief Replace 0-by-0 outputs */
    template<typename M>
    std::vector<M> replaceRes(const std::vector<M>& res) const;

    /** \brief Replace 0-by-0 forward seeds */
    template<typename M>
    std::vector<std::vector<M> >
      replaceFwdSeed(const std::vector<std::vector<M> >& fseed) const;

    /** \brief Replace 0-by-0 reverse seeds */
    template<typename M>
    std::vector<std::vector<M> >
      replaceAdjSeed(const std::vector<std::vector<M> >& aseed) const;
#endif // SWIG
    /// \endcond

#ifndef SWIG
    // Create vector with 1 element
    inline friend std::vector<Function> make_vector(const Function& x0) {
      return std::vector<Function>(1, x0);
    }

    // Create vector with 2 elements
    inline friend std::vector<Function> make_vector(const Function& x0, const Function& x1) {
      Function x[] = {x0, x1};
      return std::vector<Function>(x, x+2);
    }

    // Create vector with 3 elements
    inline friend std::vector<Function> make_vector(const Function& x0, const Function& x1,
                                                   const Function& x2) {
      Function x[] = {x0, x1, x2};
      return std::vector<Function>(x, x+3);
    }

    // Create vector with 4 elements
    inline friend std::vector<Function> make_vector(const Function& x0, const Function& x1,
                                                   const Function& x2, const Function& x3) {
      Function x[] = {x0, x1, x2, x3};
      return std::vector<Function>(x, x+4);
    }

    // Create vector with 5 elements
    inline friend std::vector<Function> make_vector(const Function& x0, const Function& x1,
                                                   const Function& x2, const Function& x3,
                                                   const Function& x4) {
      Function x[] = {x0, x1, x2, x3, x4};
      return std::vector<Function>(x, x+5);
    }

    // Create vector with 6 elements
    inline friend std::vector<Function> make_vector(const Function& x0, const Function& x1,
                                                   const Function& x2, const Function& x3,
                                                   const Function& x4, const Function& x5) {
      Function x[] = {x0, x1, x2, x3, x4, x5};
      return std::vector<Function>(x, x+6);
    }
#endif // SWIG

    /** \brief get function name with all non alphanumeric characters converted to '_' */
    std::string getSanitizedName() const;

    /** \brief get function name with all non alphanumeric characters converted to '_' */
    static std::string sanitizeName(const std::string& name);
  };

} // namespace casadi

#endif // CASADI_FUNCTION_HPP
