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

 /*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
      Number of inputs:  \a nin    getNumInputs()\n
      Number of outputs: \a nout   getNumOutputs()\n

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

      Using \f$\vec {v} \in \mathbf{R}^n\f$ as a forward seed:  <tt>setFwdSeed(v, i)</tt>\n
      Retrieving \f$\vec {s}_f \in \mathbf{R}^p\f$ from:        <tt>getFwdSens(sf, j)</tt>\n

      Using \f$\vec {w} \in \mathbf{R}^p\f$ as a forward seed:  <tt>setAdjSeed(w, j)</tt>\n
      Retrieving \f$\vec {s}_a \in \mathbf{R}^n \f$ from:        <tt>getAdjSens(sa, i)</tt>\n

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
  class CASADI_CORE_EXPORT Function : public OptionsFunctionality, public IOInterface<Function>{

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

    /// \cond INTERNAL
    ///@{
    /** \brief Access input/output scheme */
    const casadi::IOScheme& inputScheme() const;
    const casadi::IOScheme& outputScheme() const;

    casadi::IOScheme& inputScheme();
    casadi::IOScheme& outputScheme();
    ///@}
    /// \endcond

    /// \cond INTERNAL
    ///@{
    /// Input/output structures of the function */
    const IOSchemeVector<DMatrix>& input_struct() const;
    const IOSchemeVector<DMatrix>& output_struct() const;
    IOSchemeVector<DMatrix>& input_struct();
    IOSchemeVector<DMatrix>& output_struct();
    ///@}
    /// \endcond

    /** \brief  Get total number of nonzeros in all of the matrix-valued inputs */
    int getNumInputNonzeros() const;

    /** \brief  Get total number of nonzeros in all of the matrix-valued outputs */
    int getNumOutputNonzeros() const;

    /** \brief  Get total number of elements in all of the matrix-valued inputs */
    int getNumInputElements() const;

    /** \brief  Get total number of elements in all of the matrix-valued outputs */
    int getNumOutputElements() const;

    /** \brief Set input scheme */
    void setInputScheme(const casadi::IOScheme &scheme);

    /** \brief Set output scheme */
    void setOutputScheme(const casadi::IOScheme &scheme);

    /** \brief Get input scheme */
    casadi::IOScheme getInputScheme() const;

    /** \brief Get output scheme */
    casadi::IOScheme getOutputScheme() const;

    /** \brief  Evaluate */
    void evaluate();

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
        return jacobian(inputSchemeEntry(iind), oind, compact, symmetric);
    }
    Function jacobian(int iind, const std::string& oind, bool compact=false, bool symmetric=false) {
        return jacobian(iind, outputSchemeEntry(oind), compact, symmetric);
    }
    Function jacobian(const std::string& iind, const std::string& oind, bool compact=false,
                      bool symmetric=false) {
        return jacobian(inputSchemeEntry(iind), outputSchemeEntry(oind), compact, symmetric);
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
        return gradient(inputSchemeEntry(iind), oind);
    }
    Function gradient(int iind, const std::string& oind) {
        return gradient(iind, outputSchemeEntry(oind));
    }
    Function gradient(const std::string& iind, const std::string& oind) {
        return gradient(inputSchemeEntry(iind), outputSchemeEntry(oind));
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
    { return tangent(inputSchemeEntry(iind), oind); }
    Function tangent(int iind, const std::string& oind)
    { return tangent(iind, outputSchemeEntry(oind)); }
    Function tangent(const std::string& iind, const std::string& oind)
    { return tangent(inputSchemeEntry(iind), outputSchemeEntry(oind)); }
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
    { return hessian(inputSchemeEntry(iind), oind); }
    Function hessian(int iind, const std::string& oind)
    { return hessian(iind, outputSchemeEntry(oind)); }
    Function hessian(const std::string& iind, const std::string& oind)
    { return hessian(inputSchemeEntry(iind), outputSchemeEntry(oind)); }
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
    std::vector<DMatrix> call(const std::vector<DMatrix> &arg, bool always_inline=false,
                              bool never_inline=false);
    std::vector<SX> call(const std::vector<SX> &arg, bool always_inline=false,
                         bool never_inline=false);
    std::vector<MX> call(const std::vector<MX> &arg, bool always_inline=false,
                         bool never_inline=false);
    ///@}

    ///@{
    /// Functor shorthand for evaluation
    std::vector<DMatrix> operator()(const std::vector<DMatrix>& arg) { return call(arg);}
    std::vector<SX> operator()(const std::vector<SX>& arg) { return call(arg);}
    std::vector<MX> operator()(const std::vector<MX>& arg) { return call(arg);}
    IOSchemeVector<DMatrix> operator()(const IOSchemeVector<DMatrix>& arg) {
      return outputScheme().fromVector(call(arg));
    }
    IOSchemeVector<SX> operator()(const IOSchemeVector<SX>& arg) {
      return outputScheme().fromVector(call(arg));
    }
    IOSchemeVector<MX> operator()(const IOSchemeVector<MX>& arg) {
      return outputScheme().fromVector(call(arg));
    }
    ///@}

#ifndef SWIG
    ///@{
    /// Functor shorthand for evaluation, single argument (only C++)
    std::vector<DMatrix> operator()(const DMatrix& arg0) { return call(toVector(arg0));}
    std::vector<SX> operator()(const SX& arg0) { return call(toVector(arg0));}
    std::vector<MX> operator()(const MX& arg0) { return call(toVector(arg0));}
    ///@}

    ///@{
    /// Functor shorthand for evaluation, two arguments (only C++)
    std::vector<DMatrix> operator()(const DMatrix& arg0, const DMatrix& arg1)
    { return call(toVector(arg0, arg1));}
    std::vector<SX> operator()(const SX& arg0, const SX& arg1) { return call(toVector(arg0, arg1));}
    std::vector<MX> operator()(const MX& arg0, const MX& arg1) { return call(toVector(arg0, arg1));}
    ///@}
#endif // SWIG

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
        \param paropt Set of options to be passed to the Parallelizer
    */
    std::vector<std::vector<MX> > callParallel(const std::vector<std::vector<MX> > &arg,
                                               const Dictionary& paropt=Dictionary());

    /** \brief Get a function that calculates nfwd forward derivatives and nadj adjoint derivatives
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
     *         <tt>(n_in = getNumInputs(), n_out = getNumOutputs())</tt>
     *
     *        The functions returned are cached, meaning that if called multiple timed
     *        with the same value, then multiple references to the same function will be returned.
     */
    Function derivative(int nfwd, int nadj);

    /** \brief Set a function that calculates \a nfwd forward derivatives
        and \a nadj adjoint derivatives

NOTE: Does _not_ take ownership, only weak references to the derivatives are kept internally */
    void setDerivative(const Function& fcn, int nfwd, int nadj);

    ///@{
    /// Get, if necessary generate, the sparsity of a Jacobian block
    Sparsity& jacSparsity(int iind=0, int oind=0, bool compact=false, bool symmetric=false);
    Sparsity& jacSparsity(const std::string &iind, int oind=0, bool compact=false,
                          bool symmetric=false) {
        return jacSparsity(inputSchemeEntry(iind), oind, compact, symmetric); }
    Sparsity& jacSparsity(int iind, const std::string &oind, bool compact=false,
                          bool symmetric=false) {
        return jacSparsity(iind, outputSchemeEntry(oind), compact, symmetric); }
    Sparsity& jacSparsity(const std::string &iind, const std::string &oind,
                          bool compact=false, bool symmetric=false) {
        return jacSparsity(inputSchemeEntry(iind), outputSchemeEntry(oind), compact, symmetric); }
    ///@}

    ///@{
    /// Generate the sparsity of a Jacobian block
    void setJacSparsity(const Sparsity& sp, int iind, int oind, bool compact=false);
    void setJacSparsity(const Sparsity& sp, const std::string &iind, int oind, bool compact=false) {
        setJacSparsity(sp, inputSchemeEntry(iind), oind, compact); }
    void setJacSparsity(const Sparsity& sp, int iind, const std::string &oind, bool compact=false) {
        setJacSparsity(sp, iind, outputSchemeEntry(oind), compact); }
    void setJacSparsity(const Sparsity& sp, const std::string &iind, const std::string &oind,
                        bool compact=false) {
        setJacSparsity(sp, inputSchemeEntry(iind), outputSchemeEntry(oind), compact); }
    ///@}

    /** \brief Export / Generate C code for the function */
    void generateCode(const std::string& filename, bool generate_main=false);

    /** \brief Generate C code for the function */
    std::string generateCode();

    /** \brief Generate C code for the function */
    void generateCode(std::ostream& filename, bool generate_main=false);

    /// \cond INTERNAL
    /** \brief  Access functions of the node */
    FunctionInternal* operator->();

    /** \brief  Const access functions of the node */
    const FunctionInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
    /// \endcond

    /// Get all statistics obtained at the end of the last evaluate call
    const Dictionary& getStats() const;

    /// Get a single statistic obtained at the end of the last evaluate call
    GenericType getStat(const std::string& name) const;

    /** \brief  Get a vector of symbolic variables with the same dimensions as the inputs
     *
     * There is no guarantee that consecutive calls return identical objects
     */
    std::vector<MX> symbolicInput() const;

    /** \brief Get a vector of symbolic variables with the same dimensions as the inputs, SX graph
     *
     * There is no guarantee that consecutive calls return identical objects
     */
    std::vector<SX> symbolicInputSX() const;

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

    /// \endcond

    /** \brief Add modules to be monitored */
    void addMonitor(const std::string& mon);

    /** \brief Remove modules to be monitored */
    void removeMonitor(const std::string& mon);

    /// \cond INTERNAL
    /** \brief Check if the numerical values of the supplied bounds make sense */
    void checkInputs() const;
    /// \endcond

    /** \brief get function name with all non alphanumeric characters converted to '_' */
    std::string getSanitizedName() const;
  };

} // namespace casadi

#ifdef SWIG
// Template instantiations
%template(FunctionVector)             std::vector<casadi::Function>;
#endif // SWIG


#endif // CASADI_FUNCTION_HPP
