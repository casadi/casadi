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
#include "compiler.hpp"
#include "../matrix/sparse_storage.hpp"

// This macro is for documentation purposes
#define INPUTSCHEME(name)

// This macro is for documentation purposes
#define OUTPUTSCHEME(name)

/// \cond INTERNAL

namespace casadi {

  ///@{
  /** \brief  Function pointer types */
  typedef int (*sparsityPtr)(int i, int *n_row, int *n_col,
                             const int **colind, const int **row);
  typedef int (*workPtr)(int *n_iw, int *n_w);
  typedef int (*evalPtr)(const double** arg, double** res, int* iw, double* w);
  typedef void (*simplifiedPtr)(const double* arg, double* res);
  typedef int (*initPtr)(int *f_type, int *n_in, int *n_out, int *sz_arg, int* sz_res);
  ///@}

  class MXFunction;

  /** \brief Internal class for Function
      \author Joel Andersson
      \date 2010
      A regular user should never work with any Node class. Use Function directly.
  */
  class CASADI_EXPORT FunctionInternal : public OptionsFunctionalityNode,
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

    /** \brief  Evaluate using internal data structures */
    virtual void evaluate();

    /** \brief  Obtain solver name from Adaptor */
    virtual std::string getAdaptorSolverName() const { return ""; }

    /** \brief Initialize
        Initialize and make the object ready for setting arguments and evaluation.
        This method is typically called after setting options but before evaluating.
        If passed to another class (in the constructor), this class should invoke
        this function when initialized. */
    virtual void init();

    /** \brief Post-initialization
        This function, which visits the class hierarchy in reverse order is run after
        init() has been completed.
    */
    virtual void postinit();

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

    /** \brief  Evaluate numerically, possibly using just-in-time compilation */
    void eval(const double** arg, double** res, int* iw, double* w);

    /** \brief  Evaluate numerically, work vectors given */
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /** \brief Quickfix to avoid segfault, #1552 */
    virtual bool canEvalSX() const {return false;}

    /** \brief  Evaluate symbolically, SXElement type, possibly nonmatching sparsity patterns */
    virtual void evalSX(const SXElement** arg, SXElement** res,
                        int* iw, SXElement* w);

    /** \brief  Evaluate symbolically, SXElement type, possibly nonmatching sparsity patterns */
    virtual void evalSX(const std::vector<SX>& arg, std::vector<SX>& res);

    /** \brief  Evaluate symbolically, MX type */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Call a function, DMatrix type (overloaded) */
    void call(const DMatrixVector& arg, DMatrixVector& res,
              bool always_inline, bool never_inline);

    /** \brief Call a function, MX type (overloaded) */
    void call(const MXVector& arg, MXVector& res,
              bool always_inline, bool never_inline);

    /** \brief Call a function, SX type (overloaded) */
    void call(const SXVector& arg, SXVector& res,
              bool always_inline, bool never_inline);

    /** \brief Create call to (cached) derivative function, forward mode  */
    virtual void callForward(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens,
                         bool always_inline, bool never_inline);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    virtual void callReverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens,
                         bool always_inline, bool never_inline);

    /** \brief Create call to (cached) derivative function, forward mode  */
    virtual void callForward(const std::vector<SX>& arg, const std::vector<SX>& res,
                         const std::vector<std::vector<SX> >& fseed,
                         std::vector<std::vector<SX> >& fsens,
                         bool always_inline, bool never_inline);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    virtual void callReverse(const std::vector<SX>& arg, const std::vector<SX>& res,
                         const std::vector<std::vector<SX> >& aseed,
                         std::vector<std::vector<SX> >& asens,
                         bool always_inline, bool never_inline);

    /** \brief Create call to (cached) derivative function, forward mode  */
    virtual void callForward(const std::vector<DMatrix>& arg, const std::vector<DMatrix>& res,
                         const std::vector<std::vector<DMatrix> >& fseed,
                         std::vector<std::vector<DMatrix> >& fsens,
                         bool always_inline, bool never_inline);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    virtual void callReverse(const std::vector<DMatrix>& arg, const std::vector<DMatrix>& res,
                         const std::vector<std::vector<DMatrix> >& aseed,
                         std::vector<std::vector<DMatrix> >& asens,
                         bool always_inline, bool never_inline);
    ///@{
    /** \brief Return Hessian function */
    Function hessian(int iind, int oind);
    virtual Function getHessian(int iind, int oind);
    ///@}

    ///@{
    /** \brief Return gradient function */
    Function gradient(int iind, int oind);
    virtual Function getGradient(const std::string& name, int iind, int oind, const Dict& opts);
    ///@}

    ///@{
    /** \brief Return tangent function */
    Function tangent(int iind, int oind);
    virtual Function getTangent(const std::string& name, int iind, int oind, const Dict& opts);
    ///@}

    ///@{
    /** \brief Return Jacobian function */
    Function jacobian(int iind, int oind, bool compact, bool symmetric);
    void setJacobian(const Function& jac, int iind, int oind, bool compact);
    virtual Function getJacobian(const std::string& name, int iind, int oind,
                                 bool compact, bool symmetric, const Dict& opts);
    virtual Function getNumericJacobian(const std::string& name, int iind, int oind,
                                        bool compact, bool symmetric, const Dict& opts);
    ///@}

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements */
    Function fullJacobian();
    virtual bool hasFullJacobian() const;
    virtual Function getFullJacobian(const std::string& name, const Dict& opts);
    ///@}

    ///@{
    /** \brief Return function that calculates forward derivatives
     *    derForward(nfwd) returns a cached instance if available,
     *    and calls <tt>Function getDerForward(int nfwd)</tt>
     *    if no cached version is available.
     */
    Function derForward(int nfwd);
    virtual Function getDerForward(const std::string& name, int nfwd, Dict& opts);
    virtual int numDerForward() const { return hasSetOption("custom_forward") ? 64 : 0;}
    void setDerForward(const Function& fcn, int nfwd);
    ///@}

    ///@{
    /** \brief Return function that calculates adjoint derivatives
     *    derReverse(nadj) returns a cached instance if available,
     *    and calls <tt>Function getDerReverse(int nadj)</tt>
     *    if no cached version is available.
     */
    Function derReverse(int nadj);
    virtual Function getDerReverse(const std::string& name, int nadj, Dict& opts);
    virtual int numDerReverse() const { return hasSetOption("custom_reverse") ? 64 : 0;}
    void setDerReverse(const Function& fcn, int nadj);
    ///@}

    /** \brief Can derivatives be calculated in any way? */
    bool hasDerivative() const;

    /** \brief  Weighting factor for chosing forward/reverse mode */
    virtual double adWeight();

    /** \brief  Weighting factor for chosing forward/reverse mode,
        sparsity propagation */
    virtual double adWeightSp();

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

    /** \brief Generate function call */
    virtual void generate(CodeGenerator& g, const std::vector<int>& arg,
                          const std::vector<int>& res) const;

    /** \brief Generate code the function */
    virtual void generateFunction(CodeGenerator& g, const std::string& fname,
                                  bool decl_static) const;

    /** \brief Generate meta-information allowing a user to evaluate a generated function */
    void generateMeta(CodeGenerator& g, const std::string& fname) const;

    /** \brief Use simplified signature */
    virtual bool simplifiedCall() const { return false;}

    /** \brief Generate a call to a function (generic signature) */
    virtual std::string generateCall(const CodeGenerator& g,
                                     const std::string& arg, const std::string& res,
                                     const std::string& iw, const std::string& w) const;

    /** \brief Generate a call to a function (simplified signature) */
    virtual std::string generateCall(const CodeGenerator& g,
                                     const std::string& arg, const std::string& res) const;

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Generate code for the function body */
    virtual void generateBody(CodeGenerator& g) const;

    /** \brief  Print */
    virtual void print(std::ostream &stream) const;

    /** \brief  Print */
    virtual void repr(std::ostream &stream) const;

    /** \brief Check if the numerical values of the supplied bounds make sense */
    virtual void checkInputs() const {}

    /** \brief  Print dimensions of inputs and outputs */
    void printDimensions(std::ostream &stream) const;

    /** \brief Get the unidirectional or bidirectional partition */
    void getPartition(int iind, int oind, Sparsity& D1, Sparsity& D2, bool compact, bool symmetric);

    /// Verbose mode?
    bool verbose() const;

    /// Is function fcn being monitored
    bool monitored(const std::string& mod) const;

    /** \brief Get the number of function inputs */
    inline int nIn() const { return ibuf_.size();}

    /** \brief Get the number of function outputs */
    inline int nOut() const { return obuf_.size();}

    /** \brief  Get total number of nonzeros in all of the matrix-valued inputs */
    int nnzIn() const;

    /** \brief  Get total number of nonzeros in all of the matrix-valued outputs */
    int nnzOut() const;

    /** \brief  Get total number of elements in all of the matrix-valued inputs */
    int numelIn() const;

    /** \brief  Get total number of elements in all of the matrix-valued outputs */
    int numelOut() const;

    /// Get all statistics obtained at the end of the last evaluate call
    const Dict & getStats() const;

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

    /// Get a vector of symbolic variables with the same dimensions as the outputs
    virtual std::vector<MX> symbolicOutput() const;

    /// Get a vector of symbolic variables corresponding to the outputs
    virtual std::vector<MX> symbolicOutput(const std::vector<MX>& arg);

    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<SX> symbolicInputSX() const;

    /** \brief Get input scheme index by name */
    virtual int inputIndex(const std::string &name) const {
      const std::vector<std::string>& v=ischeme_;
      for (std::vector<std::string>::const_iterator i=v.begin(); i!=v.end(); ++i) {
        size_t col = i->find(':');
        if (i->compare(0, col, name)==0) return i-v.begin();
      }
      casadi_error("FunctionInternal::inputIndex: could not find entry \""
                   << name << "\". Available names are: " << v << ".");
      return -1;
    }

    /** \brief Get output scheme index by name */
    virtual int outputIndex(const std::string &name) const {
      const std::vector<std::string>& v=oscheme_;
      for (std::vector<std::string>::const_iterator i=v.begin(); i!=v.end(); ++i) {
        size_t col = i->find(':');
        if (i->compare(0, col, name)==0) return i-v.begin();
      }
      casadi_error("FunctionInternal::outputIndex: could not find entry \""
                   << name << "\". Available names are: " << v << ".");
      return -1;
    }

    /** \brief Get input scheme name by index */
    virtual std::string inputName(int ind) const {
      const std::string& s=ischeme_.at(ind);
      size_t col = s.find(':'); // Colon seprates name from description
      return s.substr(0, col);
    }

    /** \brief Get output scheme name by index */
    virtual std::string outputName(int ind) const {
      const std::string& s=oscheme_.at(ind);
      size_t col = s.find(':'); // Colon seprates name from description
      return s.substr(0, col);
    }

    /** \brief Get input scheme description by index */
    virtual std::string inputDescription(int ind) const {
      const std::string& s=ischeme_.at(ind);
      size_t col = s.find(':'); // Colon seprates name from description
      if (col==std::string::npos) {
        return std::string("No description available");
      } else {
        return s.substr(col+1);
      }
    }

    /** \brief Get output scheme description by index */
    virtual std::string outputDescription(int ind) const {
      const std::string& s=oscheme_.at(ind);
      size_t col = s.find(':'); // Colon seprates name from description
      if (col==std::string::npos) {
        return std::string("No description available");
      } else {
        return s.substr(col+1);
      }
    }

    /** \brief Get default input value */
    virtual double defaultInput(int ind) const { return 0;}

    /** \brief Get sparsity of a given input */
    /// @{
    inline Sparsity inputSparsity(int ind=0) const {
      return input(ind).sparsity();
    }
    inline Sparsity inputSparsity(const std::string& iname) const {
      return input(inputIndex(iname)).sparsity();
    }
    /// @}

    /** \brief Get sparsity of a given output */
    /// @{
    inline Sparsity outputSparsity(int ind=0) const {
      return output(ind).sparsity();
    }
    inline Sparsity outputSparsity(const std::string& iname) const {
      return output(outputIndex(iname)).sparsity();
    }
    /// @}


    /// Access input argument by index
    inline Matrix<double>& input(int i=0) {
      try {
        return ibuf_.at(i);
      } catch(std::out_of_range&) {
        std::stringstream ss;
        ss <<  "In function " << getOption("name")
           << ": input " << i << " not in interval [0, " << nIn() << ")";
        if (!isInit()) ss << std::endl << "Did you forget to initialize?";
        throw CasadiException(ss.str());
      }
    }

    /// Access input argument by name
    inline Matrix<double>& input(const std::string &iname) {
      return input(inputIndex(iname));
    }

    /// Const access input argument by index
    inline const Matrix<double>& input(int i=0) const {
      return const_cast<FunctionInternal*>(this)->input(i);
    }

    /// Const access input argument by name
    inline const Matrix<double>& input(const std::string &iname) const {
      return const_cast<FunctionInternal*>(this)->input(iname);
    }

    /// Access output argument by index
    inline Matrix<double>& output(int i=0) {
      try {
        return obuf_.at(i);
      } catch(std::out_of_range&) {
        std::stringstream ss;
        ss <<  "In function " << getOption("name")
           << ": output " << i << " not in interval [0, " << nOut() << ")";
        if (!isInit()) ss << std::endl << "Did you forget to initialize?";
        throw CasadiException(ss.str());
      }
    }

    /// Access output argument by name
    inline Matrix<double>& output(const std::string &oname) {
      return output(outputIndex(oname));
    }

    /// Const access output argument by index
    inline const Matrix<double>& output(int i=0) const {
      return const_cast<FunctionInternal*>(this)->output(i);
    }

    /// Const access output argument by name
    inline const Matrix<double>& output(const std::string &oname) const {
      return const_cast<FunctionInternal*>(this)->output(oname);
    }

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
    /** \brief  Propagate sparsity forward */
    virtual void spFwdSwitch(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdjSwitch(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief Get number of temporary variables needed */
    void sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const;

    /** \brief Get required length of arg field */
    size_t sz_arg() const { return sz_arg_;}

    /** \brief Get required length of res field */
    size_t sz_res() const { return sz_res_;}

    /** \brief Get required length of iw field */
    size_t sz_iw() const { return sz_iw_;}

    /** \brief Get required length of w field */
    size_t sz_w() const { return sz_w_;}

    /** \brief Ensure required length of arg field */
    void alloc_arg(size_t sz_arg);

    /** \brief Ensure required length of res field */
    void alloc_res(size_t sz_res);

    /** \brief Ensure required length of iw field */
    void alloc_iw(size_t sz_iw);

    /** \brief Ensure required length of w field */
    void alloc_w(size_t sz_w);

    /** \brief Ensure work vectors long enough to evaluate function */
    void alloc(const Function& f);

    /** \brief Update lengths of temporary work vectors */
    void alloc();
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

    ///@{
    /** \brief Calculate derivatives by multiplying the full Jacobian and multiplying */
    virtual bool fwdViaJac(int nfwd);
    virtual bool adjViaJac(int nadj);
    ///@}

    /// Input and output buffers
    std::vector<DMatrix> ibuf_, obuf_;

    /// Input and output scheme
    std::vector<std::string> ischeme_, oscheme_;

    /** \brief  Verbose -- for debugging purposes */
    bool verbose_;

    /** \brief  Use just-in-time compiler */
    bool jit_;
    evalPtr evalD_;

    /// Set of module names which are extra monitored
    std::set<std::string> monitors_;

    /** \brief  Dict of statistics (resulting from evaluate) */
    Dict stats_;

    /** \brief  Flag to indicate whether statistics must be gathered */
    bool gather_stats_;

    /// Cache for functions to evaluate directional derivatives (new)
    std::vector<WeakRef> derivative_fwd_, derivative_adj_;

    /// Cache for full Jacobian
    WeakRef full_jacobian_;

    /// Cache for sparsities of the Jacobian blocks
    SparseStorage<Sparsity> jac_sparsity_, jac_sparsity_compact_;

    /// Cache for Jacobians
    SparseStorage<WeakRef> jac_, jac_compact_;

    /// User-set field
    void* user_data_;

    /// Name
    std::string name_;

    /// Just-in-time compiler
    std::string compilerplugin_;
    Compiler compiler_;
    Dict jit_options_;

    bool monitor_inputs_, monitor_outputs_;

    /// Errors are thrown when NaN is produced
    bool regularity_check_;

    /// Errors are thrown if numerical values of inputs look bad
    bool inputs_check_;

    /** \brief get function name with all non alphanumeric characters converted to '_' */
    std::string getSanitizedName() const;

    /** \brief get function name with all non alphanumeric characters converted to '_' */
    static std::string sanitizeName(const std::string& name);

    /** \brief Can a derivative direction be skipped */
    template<typename MatType>
    static bool purgable(const std::vector<MatType>& seed);

    /** \brief Symbolic expressions for the forward seeds */
    template<typename MatType>
    std::vector<std::vector<MatType> > symbolicFwdSeed(int nfwd, const std::vector<MatType>& v);

    /** \brief Symbolic expressions for the adjoint seeds */
    template<typename MatType>
    std::vector<std::vector<MatType> > symbolicAdjSeed(int nadj, const std::vector<MatType>& v);

  protected:
    /** \brief  Temporary vector needed for the evaluation (integer) */
    std::vector<int> iw_tmp_;

    /** \brief  Temporary vector needed for the evaluation (real) */
    std::vector<double> w_tmp_;

  private:
    /** \brief Sizes of input and output buffers */
    size_t sz_arg_, sz_res_, sz_iw_, sz_w_;
  };

  // Template implementations
  template<typename MatType>
  bool FunctionInternal::purgable(const std::vector<MatType>& v) {
    for (typename std::vector<MatType>::const_iterator i=v.begin(); i!=v.end(); ++i) {
      if (!i->isZero()) return false;
    }
    return true;
  }

  template<typename MatType>
  std::vector<std::vector<MatType> >
  FunctionInternal::symbolicFwdSeed(int nfwd, const std::vector<MatType>& v) {
    std::vector<std::vector<MatType> > fseed(nfwd, v);
    for (int dir=0; dir<nfwd; ++dir) {
      // Replace symbolic inputs
      int iind=0;
      for (typename std::vector<MatType>::iterator i=fseed[dir].begin();
          i!=fseed[dir].end();
          ++i, ++iind) {
        // Name of the forward seed
        std::stringstream ss;
        ss << "f";
        if (nfwd>1) ss << dir;
        ss << "_";
        ss << iind;

        // Save to matrix
        *i = MatType::sym(ss.str(), i->sparsity());

      }
    }
    return fseed;
  }

  template<typename MatType>
  std::vector<std::vector<MatType> >
  FunctionInternal::symbolicAdjSeed(int nadj, const std::vector<MatType>& v) {
    std::vector<std::vector<MatType> > aseed(nadj, v);
    for (int dir=0; dir<nadj; ++dir) {
      // Replace symbolic inputs
      int oind=0;
      for (typename std::vector<MatType>::iterator i=aseed[dir].begin();
          i!=aseed[dir].end();
          ++i, ++oind) {
        // Name of the adjoint seed
        std::stringstream ss;
        ss << "a";
        if (nadj>1) ss << dir << "_";
        ss << oind;

        // Save to matrix
        *i = MatType::sym(ss.str(), i->sparsity());

      }
    }
    return aseed;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_FUNCTION_INTERNAL_HPP
