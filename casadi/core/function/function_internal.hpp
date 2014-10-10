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


#ifndef CASADI_FUNCTION_INTERNAL_HPP
#define CASADI_FUNCTION_INTERNAL_HPP

#include "function.hpp"
#include "../functor.hpp"
#include "../weak_ref.hpp"
#include <set>
#include "code_generator.hpp"
#include "../matrix/sparse_storage.hpp"

// This macro is for documentation purposes
#define INPUTSCHEME(name)

// This macro is for documentation purposes
#define OUTPUTSCHEME(name)

/// \cond INTERNAL

namespace casadi {

  class MXFunction;

  /** \brief Internal class for Function
      \author Joel Andersson
      \date 2010
      A regular user should never work with any Node class. Use Function directly.
  */
  class CASADI_CORE_EXPORT FunctionInternal : public OptionsFunctionalityNode,
                                                  public IOInterface<FunctionInternal>{
    friend class Function;

  protected:
    /** \brief  Default constructor (accessible from the Function class and derived classes) */
    FunctionInternal();

  public:
    /** \brief  Destructor */
    virtual ~FunctionInternal() = 0;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Evaluate */
    virtual void evaluate() = 0;

    /** \brief Initialize
        Initialize and make the object ready for setting arguments and evaluation.
        This method is typically called after setting options but before evaluating.
        If passed to another class (in the constructor), this class should invoke
        this function when initialized. */
    virtual void init();

    /** \brief  Propagate the sparsity pattern through a set of directional
        derivatives forward or backward */
    virtual void spEvaluate(bool fwd);

    /** \brief  Propagate the sparsity pattern through a set of directional derivatives
        forward or backward, using the sparsity patterns */
    virtual void spEvaluateViaJacSparsity(bool fwd);

    /** \brief  Is the class able to propagate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd) { return false;}

    /** \brief  Reset the sparsity propagation */
    virtual void spInit(bool fwd) {}

    /** \brief  Evaluate symbolically, SXElement type, possibly nonmatching sparsity patterns */
    virtual void evalSX(const std::vector<SX>& arg, std::vector<SX>& res,
                        const std::vector<std::vector<SX> >& fseed,
                        std::vector<std::vector<SX> >& fsens,
                        const std::vector<std::vector<SX> >& aseed,
                        std::vector<std::vector<SX> >& asens);

    /** \brief  Evaluate symbolically, MX type */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res,
                        const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens,
                        const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens);

    /** \brief  Evaluate symbolically, SXElement type, matching sparsity patterns */
    virtual void evalSXsparse(const std::vector<SX>& arg, std::vector<SX>& res,
                              const std::vector<std::vector<SX> >& fseed,
                              std::vector<std::vector<SX> >& fsens,
                              const std::vector<std::vector<SX> >& aseed,
                              std::vector<std::vector<SX> >& asens);

    /** \brief  Create function call node */
    virtual void createCall(const std::vector<MX> &arg, std::vector<MX> &res,
                            const std::vector<std::vector<MX> > &fseed,
                            std::vector<std::vector<MX> > &fsens,
                            const std::vector<std::vector<MX> > &aseed,
                            std::vector<std::vector<MX> > &asens);

    /** \brief  Create derivative node */
    virtual void createCallDerivative(const std::vector<MX> &arg, std::vector<MX> &res,
                                      const std::vector<std::vector<MX> > &fseed,
                                      std::vector<std::vector<MX> > &fsens,
                                      const std::vector<std::vector<MX> > &aseed,
                                      std::vector<std::vector<MX> > &asens);

    /** \brief  Create a call to this */
    std::vector<MX> callSelf(const std::vector<MX> &arg);

    /** \brief Call a function, DMatrix type (overloaded) */
    void call(const DMatrixVector& arg, DMatrixVector& res,
              const DMatrixVectorVector& fseed, DMatrixVectorVector& fsens,
              const DMatrixVectorVector& aseed, DMatrixVectorVector& asens,
              bool always_inline, bool never_inline);

    /** \brief Call a function, MX type (overloaded) */
    void call(const MXVector& arg, MXVector& res,
              const MXVectorVector& fseed, MXVectorVector& fsens,
              const MXVectorVector& aseed, MXVectorVector& asens,
              bool always_inline, bool never_inline);

    /** \brief Call a function, SX type (overloaded) */
    void call(const SXVector& arg, SXVector& res,
              const SXVectorVector& fseed, SXVectorVector& fsens,
              const SXVectorVector& aseed, SXVectorVector& asens,
              bool always_inline, bool never_inline);

    ///@{
    /** \brief Return Hessian function */
    Function hessian(int iind, int oind);
    virtual Function getHessian(int iind, int oind);
    ///@}

    ///@{
    /** \brief Return gradient function */
    Function gradient(int iind, int oind);
    virtual Function getGradient(int iind, int oind);
    ///@}

    ///@{
    /** \brief Return tangent function */
    Function tangent(int iind, int oind);
    virtual Function getTangent(int iind, int oind);
    ///@}

    ///@{
    /** \brief Return Jacobian function */
    Function jacobian(int iind, int oind, bool compact, bool symmetric);
    void setJacobian(const Function& jac, int iind, int oind, bool compact);
    virtual Function getJacobian(int iind, int oind, bool compact, bool symmetric);
    virtual Function getNumericJacobian(int iind, int oind, bool compact, bool symmetric);
    ///@}

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements */
    Function fullJacobian();
    virtual Function getFullJacobian();
    ///@}

    ///@{
    /** \brief Return function that calculates forward derivatives
     *
     *    This method returns a cached instance if available,
     *    and calls <tt>Function getDerivative(int nfwd, int nadj)</tt>
     *    if no cached version is available.
     */
    Function derivative(int nfwd, int nadj);

    /** Set a function that calculates nfwd forward derivatives and nadj adjoint derivatives */
    void setDerivative(const Function& fcn, int nfwd, int nadj);

    /** \brief Constructs and returns a function that calculates forward derivatives */
    virtual Function getDerivative(int nfwd, int nadj);

    /** \brief Constructs and returns a function that calculates forward derivatives
     *
     *  by creating the Jacobian then multiplying */
    virtual Function getDerivativeViaJac(int nfwd, int nadj);

    ///@}


    /** \brief Create a helper MXFunction with some properties copied
    *
    * Copied properties:
    *
    *    input/outputscheme
    *    ad_mode
    *
    *  The function is not initialized
    */
    MXFunction wrapMXFunction();

    /** \brief  Print to a stream */
    virtual void generateCode(std::ostream &cfile, bool generate_main);

    /** \brief Generate code for function inputs and outputs */
    void generateIO(CodeGenerator& gen);

    /** \brief Generate code the function */
    virtual void generateFunction(std::ostream &stream, const std::string& fname,
                                  const std::string& input_type, const std::string& output_type,
                                  const std::string& type, CodeGenerator& gen) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(std::ostream &stream, const std::string& type,
                                      CodeGenerator& gen) const;

    /** \brief Generate code for the function body */
    virtual void generateBody(std::ostream &stream, const std::string& type,
                              CodeGenerator& gen) const;

    /** \brief  Print */
    virtual void print(std::ostream &stream) const;

    /** \brief  Print */
    virtual void repr(std::ostream &stream) const;

    /** \brief Check if the numerical values of the supplied bounds make sense */
    virtual void checkInputs() const {}

    /** \brief Get the unidirectional or bidirectional partition */
    void getPartition(int iind, int oind, Sparsity& D1, Sparsity& D2, bool compact, bool symmetric);

    /// Verbose mode?
    bool verbose() const;

    /// Is function fcn being monitored
    bool monitored(const std::string& mod) const;

    /** \brief  Get total number of nonzeros in all of the matrix-valued inputs */
    int getNumInputNonzeros() const;

    /** \brief  Get total number of nonzeros in all of the matrix-valued outputs */
    int getNumOutputNonzeros() const;

    /** \brief  Get total number of elements in all of the matrix-valued inputs */
    int getNumInputElements() const;

    /** \brief  Get total number of elements in all of the matrix-valued outputs */
    int getNumOutputElements() const;

    /// Get all statistics obtained at the end of the last evaluate call
    const Dictionary & getStats() const;

    /// Get single statistic obtained at the end of the last evaluate call
    GenericType getStat(const std::string & name) const;

    /// Generate the sparsity of a Jacobian block
    virtual Sparsity getJacSparsity(int iind, int oind, bool symmetric);

    /// A flavor of getJacSparsity without any magic
    Sparsity getJacSparsityPlain(int iind, int oind);

    /// A flavor of getJacSparsity that does hierarchical block structure recognition
    Sparsity getJacSparsityHierarchical(int iind, int oind);

    /** A flavor of getJacSparsity that does hierarchical block
    * structure recognition for symmetric Jacobians
    */
    Sparsity getJacSparsityHierarchicalSymm(int iind, int oind);

    /// Generate the sparsity of a Jacobian block
    void setJacSparsity(const Sparsity& sp, int iind, int oind, bool compact);

    /// Get, if necessary generate, the sparsity of a Jacobian block
    Sparsity& jacSparsity(int iind, int oind, bool compact, bool symmetric);

    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<MX> symbolicInput() const;

    /// Get a vector of symbolic variables corresponding to the outputs
    virtual std::vector<MX> symbolicOutput(const std::vector<MX>& arg);

    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<SX> symbolicInputSX() const;

    // Workaround helper functions: assign nonzeros but ignore all -1
    static void assignIgnore(MX& y, const MX& x, const std::vector<int>& nz);
    static void assignIgnore(SX& y, const SX& x, const std::vector<int>& nz);

    ///@{
    /** \brief Access input/output scheme */
    inline const IOScheme& inputScheme() const { return input_.scheme;}
    inline const IOScheme& outputScheme() const { return output_.scheme;}
    inline IOScheme& inputScheme() { return input_.scheme;}
    inline IOScheme& outputScheme() { return output_.scheme;}
    ///@}

    ///@{
    /// Input/output structures of the function */
    inline const IOSchemeVector<DMatrix>& input_struct() const { return input_;}
    inline const IOSchemeVector<DMatrix>& output_struct() const { return output_;}
    inline IOSchemeVector<DMatrix>& input_struct() { return input_;}
    inline IOSchemeVector<DMatrix>& output_struct() { return output_;}
    ///@}

    ///@{
    /// Input/output access without checking (faster, but unsafe)
    inline const Matrix<double>& inputNoCheck(int iind=0) const { return inputS<false>(iind);}
    inline const Matrix<double>& outputNoCheck(int oind=0) const { return outputS<false>(oind);}

    inline Matrix<double>& inputNoCheck(int iind=0) { return inputS<false>(iind);}
    inline Matrix<double>& outputNoCheck(int oind=0) { return outputS<false>(oind);}
    ///@}

    /** \brief  Log the status of the solver */
    void log(const std::string& msg) const;

    /** \brief  Log the status of the solver, function given */
    void log(const std::string& fcn, const std::string& msg) const;

    /// Codegen function
    Function dynamicCompilation(Function f, std::string fname, std::string fdescr,
                                std::string compiler);

    /// The following functions are called internally from EvaluateMX.
    /// For documentation, see the MXNode class
    ///@{
    virtual void evaluateD(MXNode* node, const DMatrixPtrV& arg, DMatrixPtrV& res,
                           std::vector<int>& itmp, std::vector<double>& rtmp);
    virtual void evaluateSX(MXNode* node, const SXPtrV& arg, SXPtrV& res,
                            std::vector<int>& itmp, std::vector<SXElement>& rtmp);
    virtual void evaluateMX(MXNode* node, const MXPtrV& arg, MXPtrV& res,
                            const MXPtrVV& fseed, MXPtrVV& fsens, const MXPtrVV& aseed,
                            MXPtrVV& asens, bool output_given);
    virtual void propagateSparsity(MXNode* node, DMatrixPtrV& input,
                                   DMatrixPtrV& output, std::vector<int>& itmp,
                                   std::vector<double>& rtmp, bool fwd);
    virtual void nTmp(MXNode* node, size_t& ni, size_t& nr);
    virtual void generateOperation(const MXNode* node, std::ostream &stream,
                                   const std::vector<std::string>& arg,
                                   const std::vector<std::string>& res, CodeGenerator& gen) const;
    virtual void printPart(const MXNode* node, std::ostream &stream, int part) const;
    ///@}

    /** \brief Prints out a human readable report about possible constraint violations
    *          - specific constraints
    *
    * Constraint visualizer strip:
    * \verbatim
    *  o-------=-------o   Indicates that the value is nicely inbetween the bounds
    *  o-=-------------o   Indicates that the value is closer to the lower bound
    *  X---------------o   Indicates that the lower bound is active
    *  8---------------o   Indicates that the lower bound is -infinity
    *  o------------=--o   Indicates that the value is closer to the upper bound
    *  o---------------X   Indicates that the upper bound is active
    *  o---------------8   Indicates that the upper bound is infinity
    *     VIOLATED         Indicates constraint violation
    * \endverbatim
    */
    static void reportConstraints(std::ostream &stream, const Matrix<double> &v,
                                  const Matrix<double> &lb, const Matrix<double> &ub,
                                  const std::string &name, double tol=1e-8);

    /** \brief  Inputs of the function */
    IOSchemeVector<DMatrix> input_;

    /** \brief  Output of the function */
    IOSchemeVector<DMatrix> output_;

    /** \brief  Verbose -- for debugging purposes */
    bool verbose_;

    /// Set of module names which are extra monitored
    std::set<std::string> monitors_;

    /** \brief  Dictionary of statistics (resulting from evaluate) */
    Dictionary stats_;

    /** \brief  Flag to indicate whether statistics must be gathered */
    bool gather_stats_;

    /// Cache for functions to evaluate directional derivatives
    std::vector<std::vector<WeakRef> > derivative_fcn_;

    /// Cache for full Jacobian
    WeakRef full_jacobian_;

    /// Cache for sparsities of the Jacobian blocks
    SparseStorage<Sparsity> jac_sparsity_, jac_sparsity_compact_;

    /// Cache for Jacobians
    SparseStorage<WeakRef> jac_, jac_compact_;

    /// User-set field
    void* user_data_;

    bool monitor_inputs_, monitor_outputs_;

    /// Errors are thrown when NaN is produced
    bool regularity_check_;

    /// Errors are thrown if numerical values of inputs look bad
    bool inputs_check_;

    /** \brief get function name with all non alphanumeric characters converted to '_' */
    std::string getSanitizedName() const;

  };


} // namespace casadi

/// \endcond

#endif // CASADI_FUNCTION_INTERNAL_HPP
