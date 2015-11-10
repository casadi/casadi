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
  template<typename T>
  std::vector<std::pair<std::string, T>> zip(const std::vector<std::string>& id,
                                             const std::vector<T>& mat) {
    casadi_assert(id.size()==mat.size());
    std::vector<std::pair<std::string, T>> r(id.size());
    for (unsigned int i=0; i<r.size(); ++i) r[i] = make_pair(id[i], mat[i]);
    return r;
  }

  ///@{
  /** \brief  Function pointer types */
  typedef int (*init_t)(void **mem, int *n_in, int *n_out, void *data);
  typedef int (*freemem_t)(void* mem);
  typedef int (*work_t)(void* mem, int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
  typedef int (*sparsity_t)(void* mem, int i, int *n_row, int *n_col,
                            const int **colind, const int **row);
  typedef int (*eval_t)(const double** arg, double** res, int* iw, double* w, void* mem);
  typedef void (*simple_t)(const double* arg, double* res);
  ///@}

  /** \brief Internal class for Function
      \author Joel Andersson
      \date 2010
      A regular user should never work with any Node class. Use Function directly.
  */
  class CASADI_EXPORT FunctionInternal : public OptionsFunctionalityNode,
                                         public IOInterface<FunctionInternal>{
    friend class Function;

  protected:
    /** \brief Constructor (accessible from the Function class and derived classes) */
    FunctionInternal(const std::string& name);

  public:
    /** \brief  Destructor */
    virtual ~FunctionInternal() = 0;

    /** \brief  Obtain solver name from Adaptor */
    virtual std::string getAdaptorSolverName() const { return ""; }

    /** \brief Initialize
        Initialize and make the object ready for setting arguments and evaluation.
        This method is typically called after setting options but before evaluating.
        If passed to another class (in the constructor), this class should invoke
        this function when initialized. */
    virtual void init();

    /** \brief Finalize the object creation
        This function, which visits the class hierarchy in reverse order is run after
        init() has been completed.
    */
    virtual void finalize();

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

    /** \brief  Evaluate numerically */
    virtual void eval(const double** arg, double** res, int* iw, double* w, void* mem);

    /** \brief  Evaluate symbolically, SXElem type, possibly nonmatching sparsity patterns */
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem);

    ///@{
    /** \brief  Evaluate symbolically, MX type */
    virtual void eval_mx(const MXVector& arg, MXVector& res, bool always_inline, bool never_inline);
    ///@}

    ///@{
    /** \brief Evaluate a function, overloaded */
    void _eval(const double** arg, double** res, int* iw, double* w, void* mem);
    void _eval(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem);
    void _eval(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem);
    ///@}

    ///@{
    /** \brief Call a function, overloaded */
    void _call(const MXVector& arg, MXVector& res, bool always_inline, bool never_inline) {
      eval_mx(arg, res, always_inline, never_inline);
    }
    template<typename D>
    void _call(const std::vector<Matrix<D> >& arg, std::vector<Matrix<D> >& res,
               bool always_inline, bool never_inline);
    ///@}

    /** \brief Call a function, templated */
    template<typename M>
      void call(const std::vector<M>& arg, std::vector<M>& res,
               bool always_inline, bool never_inline);

    /** Helper function */
    static bool checkMat(const Sparsity& arg, const Sparsity& inp, bool hcat=false);

    /** \brief Check if input arguments have correct length and dimensions
     *
     * \param hcat check if horizontal repetion of the function input is allowed
     */
    template<typename M>
      void checkArg(const std::vector<M>& arg, bool hcat=false) const;

    /** \brief Check if output arguments have correct length and dimensions */
    template<typename M>
      void checkRes(const std::vector<M>& res) const;

    /** \brief Check if input arguments that needs to be replaced
     * 
     * \param hcat check if horizontal repetion of the function input is allowed
     */
    template<typename M>
      bool matchingArg(const std::vector<M>& arg, bool hcat=false) const;

    /** \brief Check if output arguments that needs to be replaced */
    template<typename M>
      bool matchingRes(const std::vector<M>& arg) const;

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

    /** \brief Create call */
    virtual std::vector<MX> create_call(const std::vector<MX>& arg);

    ///@{
    /** \brief Forward mode AD, virtual functions overloaded in derived classes */
    virtual void forward_mx(const std::vector<MX>& arg, const std::vector<MX>& res,
                            const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens,
                            bool always_inline, bool never_inline);
    virtual void forward_sx(const std::vector<SX>& arg, const std::vector<SX>& res,
                            const std::vector<std::vector<SX> >& fseed,
                            std::vector<std::vector<SX> >& fsens,
                            bool always_inline, bool never_inline);
    virtual void forward_dm(const std::vector<DM>& arg, const std::vector<DM>& res,
                            const std::vector<std::vector<DM> >& fseed,
                            std::vector<std::vector<DM> >& fsens,
                            bool always_inline, bool never_inline);
    ///@}

    ///@{
    /** \brief Forward mode AD, overloaded */
    void _forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                  const std::vector<std::vector<MX> >& fseed,
                  std::vector<std::vector<MX> >& fsens,
                  bool always_inline, bool never_inline) {
      forward_mx(arg, res, fseed, fsens, always_inline, never_inline);
    }

    void _forward(const std::vector<SX>& arg, const std::vector<SX>& res,
                  const std::vector<std::vector<SX> >& fseed,
                  std::vector<std::vector<SX> >& fsens,
                  bool always_inline, bool never_inline) {
      forward_sx(arg, res, fseed, fsens, always_inline, never_inline);
    }

    void _forward(const std::vector<DM>& arg, const std::vector<DM>& res,
                  const std::vector<std::vector<DM> >& fseed,
                  std::vector<std::vector<DM> >& fsens,
                  bool always_inline, bool never_inline) {
      forward_dm(arg, res, fseed, fsens, always_inline, never_inline);
    }

    ///@}

    /** \brief Forward mode AD, templated */
    template<typename M>
      void forward(const std::vector<M>& arg, const std::vector<M>& res,
                   const std::vector<std::vector<M> >& fseed,
                   std::vector<std::vector<M> >& fsens,
                   bool always_inline, bool never_inline);

    ///@{
    /** \brief Reverse mode, virtual functions overloaded in derived classes */
    virtual void reverse_mx(const std::vector<MX>& arg, const std::vector<MX>& res,
                            const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens,
                            bool always_inline, bool never_inline);
    virtual void reverse_sx(const std::vector<SX>& arg, const std::vector<SX>& res,
                            const std::vector<std::vector<SX> >& aseed,
                            std::vector<std::vector<SX> >& asens,
                            bool always_inline, bool never_inline);
    virtual void reverse_dm(const std::vector<DM>& arg, const std::vector<DM>& res,
                            const std::vector<std::vector<DM> >& aseed,
                            std::vector<std::vector<DM> >& asens,
                            bool always_inline, bool never_inline);
    ///@}

    ///@{
    /** \brief Reverse mode AD, overloaded */
    void _reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                  const std::vector<std::vector<MX> >& aseed,
                  std::vector<std::vector<MX> >& asens,
                  bool always_inline, bool never_inline) {
      reverse_mx(arg, res, aseed, asens, always_inline, never_inline);
    }

    void _reverse(const std::vector<SX>& arg, const std::vector<SX>& res,
                  const std::vector<std::vector<SX> >& aseed,
                  std::vector<std::vector<SX> >& asens,
                  bool always_inline, bool never_inline) {
      reverse_sx(arg, res, aseed, asens, always_inline, never_inline);
    }

    void _reverse(const std::vector<DM>& arg, const std::vector<DM>& res,
                  const std::vector<std::vector<DM> >& aseed,
                  std::vector<std::vector<DM> >& asens,
                  bool always_inline, bool never_inline) {
      reverse_dm(arg, res, aseed, asens, always_inline, never_inline);
    }
    ///@}

    /** \brief Reverse mode AD, templated */
    template<typename M>
      void reverse(const std::vector<M>& arg, const std::vector<M>& res,
                   const std::vector<std::vector<M> >& aseed,
                   std::vector<std::vector<M> >& asens,
                   bool always_inline, bool never_inline);

    ///@{
    /** \brief Parallel evaluation */
    std::vector<std::vector<MX> > map_mx(const std::vector<std::vector<MX> > &arg,
                                         const std::string& parallelization);
    std::vector<MX> map_mx(const std::vector<MX > &arg, const std::string& parallelization);
    std::vector<MX> mapsum_mx(const std::vector<MX > &arg, const std::string& parallelization);
    ///@}

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
     *    forward(nfwd) returns a cached instance if available,
     *    and calls <tt>Function get_forward(int nfwd)</tt>
     *    if no cached version is available.
     */
    Function forward(int nfwd);
    virtual Function get_forward(const std::string& name, int nfwd, Dict& opts);
    virtual int get_n_forward() const { return 0;}
    void set_forward(const Function& fcn, int nfwd);
    ///@}

    ///@{
    /** \brief Return function that calculates adjoint derivatives
     *    reverse(nadj) returns a cached instance if available,
     *    and calls <tt>Function get_reverse(int nadj)</tt>
     *    if no cached version is available.
     */
    Function reverse(int nadj);
    virtual Function get_reverse(const std::string& name, int nadj, Dict& opts);
    virtual int get_n_reverse() const { return 0;}
    void set_reverse(const Function& fcn, int nadj);
    ///@}

    /** \brief Can derivatives be calculated in any way? */
    bool hasDerivative() const;

    /** \brief  Weighting factor for chosing forward/reverse mode */
    virtual double adWeight();

    /** \brief  Weighting factor for chosing forward/reverse mode,
        sparsity propagation */
    virtual double adWeightSp();

    /** \brief Gradient expression */
    virtual MX grad_mx(int iind=0, int oind=0);

    /** \brief Tangent expression */
    virtual MX tang_mx(int iind=0, int oind=0);

    /** \brief Jacobian expression */
    virtual MX jac_mx(int iind=0, int oind=0, bool compact=false, bool symmetric=false,
                      bool always_inline=true, bool never_inline=false);

    /** \brief Gradient expression */
    virtual SX grad_sx(int iind=0, int oind=0);

    /** \brief Tangent expression */
    virtual SX tang_sx(int iind=0, int oind=0);

    /** \brief Jacobian expression */
    virtual SX jac_sx(int iind=0, int oind=0, bool compact=false, bool symmetric=false,
                      bool always_inline=true, bool never_inline=false);

    /** \brief Hessian expression */
    virtual SX hess_sx(int iind=0, int oind=0);

    ///@{
    /** \brief Get function input(s) and output(s)  */
    virtual const SX sx_in(int ind) const;
    virtual const SX sx_out(int ind) const;
    virtual const std::vector<SX> sx_in() const;
    virtual const std::vector<SX> sx_out() const;
    virtual const MX mx_in(int ind) const;
    virtual const MX mx_out(int ind) const;
    virtual const std::vector<MX> mx_in() const;
    virtual const std::vector<MX> mx_out() const;
    ///@}

    /// Get free variables (MX)
    virtual std::vector<MX> free_mx() const;

    /// Get free variables (SX)
    virtual SX free_sx() const;

    /** \brief Extract the functions needed for the Lifted Newton method */
    virtual void generate_lifted(Function& vdef_fcn, Function& vinit_fcn);

    /** \brief Get the number of atomic operations */
    virtual int getAlgorithmSize() const;

    /** \brief Get the length of the work vector */
    virtual int getWorkSize() const;

    /** \brief Get an atomic operation operator index */
    virtual int getAtomicOperation(int k) const;

    /** \brief Get the (integer) input arguments of an atomic operation */
    virtual std::pair<int, int> getAtomicInput(int k) const;

    /** \brief Get the floating point output argument of an atomic operation */
    virtual double getAtomicInputReal(int k) const;

    /** \brief Get the (integer) output argument of an atomic operation */
    virtual int getAtomicOutput(int k) const;

    /** \brief Number of nodes in the algorithm */
    virtual int countNodes() const;

    /** \brief Create a helper MXFunction with some properties copied
    *
    * Copied properties:
    *
    *    input/outputscheme
    *    ad_mode
    *
    *  The function is not initialized
    */
    Function wrapMXFunction();

    /** \brief Generate function call */
    virtual void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const;

    /** \brief Generate code the function */
    virtual void generateFunction(CodeGenerator& g, const std::string& fname,
                                  bool decl_static) const;

    /** \brief Generate meta-information allowing a user to evaluate a generated function */
    void generateMeta(CodeGenerator& g, const std::string& fname) const;

    /** \brief Use simplified signature */
    virtual bool simplifiedCall() const { return false;}

    /** \brief Generate a call to a function (generic signature) */
    virtual std::string generic_call(const CodeGenerator& g, const std::string& arg,
                                     const std::string& res, const std::string& iw,
                                     const std::string& w, const std::string& mem) const;

    /** \brief Generate a call to a function (simplified signature) */
    virtual std::string simple_call(const CodeGenerator& g,
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

    ///@{
    /** \brief Number of function inputs and outputs */
    inline int n_in() const { return ibuf_.size();}
    virtual size_t get_n_in() const = 0;
    inline int n_out() const { return obuf_.size();}
    virtual size_t get_n_out() const = 0;
    ///@}

    ///@{
    /** \brief Number of input/output nonzeros */
    int nnz_in() const;
    int nnz_in(int ind) const { return sparsity_in(ind).nnz(); }
    int nnz_out() const;
    int nnz_out(int ind) const { return sparsity_out(ind).nnz(); }
    ///@}

    ///@{
    /** \brief Number of input/output elements */
    int numel_in() const;
    int numel_in(int ind) const { return input(ind).numel(); }
    int numel_out(int ind) const { return output(ind).numel(); }
    int numel_out() const;
    ///@}

    ///@{
    /** \brief Input/output dimensions */
    int size1_in(int ind) const { return input(ind).size1(); }
    int size2_in(int ind) const { return input(ind).size2(); }
    int size1_out(int ind) const { return output(ind).size1(); }
    int size2_out(int ind) const { return output(ind).size2(); }
    std::pair<int, int> size_in(int ind) const { return input(ind).size(); }
    std::pair<int, int> size_out(int ind) const { return output(ind).size(); }
    ///@}

    /** \brief Name of the function */
    const std::string& name() const { return name_;}

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
    void set_jac_sparsity(const Sparsity& sp, int iind, int oind, bool compact);

    /// Get, if necessary generate, the sparsity of a Jacobian block
    Sparsity& sparsity_jac(int iind, int oind, bool compact, bool symmetric);

    /// Get a vector of symbolic variables corresponding to the outputs
    virtual std::vector<MX> symbolicOutput(const std::vector<MX>& arg);

    /** \brief Get input scheme index by name */
    virtual int index_in(const std::string &name) const {
      const std::vector<std::string>& v=ischeme_;
      for (std::vector<std::string>::const_iterator i=v.begin(); i!=v.end(); ++i) {
        size_t col = i->find(':');
        if (i->compare(0, col, name)==0) return i-v.begin();
      }
      casadi_error("FunctionInternal::index_in: could not find entry \""
                   << name << "\". Available names are: " << v << ".");
      return -1;
    }

    /** \brief Get output scheme index by name */
    virtual int index_out(const std::string &name) const {
      const std::vector<std::string>& v=oscheme_;
      for (std::vector<std::string>::const_iterator i=v.begin(); i!=v.end(); ++i) {
        size_t col = i->find(':');
        if (i->compare(0, col, name)==0) return i-v.begin();
      }
      casadi_error("FunctionInternal::index_out: could not find entry \""
                   << name << "\". Available names are: " << v << ".");
      return -1;
    }

    /** \brief Get input scheme name by index */
    virtual std::string name_in(int ind) const {
      const std::string& s=ischeme_.at(ind);
      size_t col = s.find(':'); // Colon seprates name from description
      return s.substr(0, col);
    }

    /** \brief Get output scheme name by index */
    virtual std::string name_out(int ind) const {
      const std::string& s=oscheme_.at(ind);
      size_t col = s.find(':'); // Colon seprates name from description
      return s.substr(0, col);
    }

    /** \brief Get input scheme description by index */
    virtual std::string description_in(int ind) const {
      const std::string& s=ischeme_.at(ind);
      size_t col = s.find(':'); // Colon seprates name from description
      if (col==std::string::npos) {
        return std::string("No description available");
      } else {
        return s.substr(col+1);
      }
    }

    /** \brief Get output scheme description by index */
    virtual std::string description_out(int ind) const {
      const std::string& s=oscheme_.at(ind);
      size_t col = s.find(':'); // Colon seprates name from description
      if (col==std::string::npos) {
        return std::string("No description available");
      } else {
        return s.substr(col+1);
      }
    }

    /** \brief Get default input value */
    static const double& default_zero();
    static const double& default_inf();
    static const double& default_minf();
    virtual const double& default_in(int ind) const {
      return default_zero();
    }

    /** \brief Get sparsity of a given input */
    /// @{
    inline const Sparsity& sparsity_in(int ind) const {
      return input(ind).sparsity();
    }
    inline const Sparsity& sparsity_in(const std::string& iname) const {
      return sparsity_in(index_in(iname));
    }
    virtual Sparsity get_sparsity_in(int ind) const = 0;
    /// @}

    /** \brief Get sparsity of a given output */
    /// @{
    inline const Sparsity& sparsity_out(int ind) const {
      return output(ind).sparsity();
    }
    inline const Sparsity& sparsity_out(const std::string& iname) const {
      return sparsity_out(index_out(iname));
    }
    virtual Sparsity get_sparsity_out(int ind) const = 0;
    /// @}

    /// Access input argument by index
    inline Matrix<double>& input(int i=0) {
      try {
        return ibuf_.at(i);
      } catch(std::out_of_range&) {
        std::stringstream ss;
        ss <<  "In function " << name_
           << ": input " << i << " not in interval [0, " << n_in() << ")";
        throw CasadiException(ss.str());
      }
    }

    /// Access input argument by name
    inline Matrix<double>& input(const std::string &iname) {
      return input(index_in(iname));
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
        ss <<  "In function " << name_
           << ": output " << i << " not in interval [0, " << n_out() << ")";
        throw CasadiException(ss.str());
      }
    }

    /// Access output argument by name
    inline Matrix<double>& output(const std::string &oname) {
      return output(index_out(oname));
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
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem);

    /** \brief Get number of temporary variables needed */
    void sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const;

    /** \brief Get required length of arg field */
    size_t sz_arg() const { return sz_arg_per_ + sz_arg_tmp_;}

    /** \brief Get required length of res field */
    size_t sz_res() const { return sz_res_per_ + sz_res_tmp_;}

    /** \brief Get required length of iw field */
    size_t sz_iw() const { return sz_iw_per_ + sz_iw_tmp_;}

    /** \brief Get required length of w field */
    size_t sz_w() const { return sz_w_per_ + sz_w_tmp_;}

    /** \brief Ensure required length of arg field */
    void alloc_arg(size_t sz_arg, bool persistent=false);

    /** \brief Ensure required length of res field */
    void alloc_res(size_t sz_res, bool persistent=false);

    /** \brief Ensure required length of iw field */
    void alloc_iw(size_t sz_iw, bool persistent=false);

    /** \brief Ensure required length of w field */
    void alloc_w(size_t sz_w, bool persistent=false);

    /** \brief Ensure work vectors long enough to evaluate function */
    void alloc(const Function& f, bool persistent=false);

    /** \brief Update lengths of temporary work vectors */
    void alloc();
    ///@}

    ///@{
    /** \brief Calculate derivatives by multiplying the full Jacobian and multiplying */
    virtual bool fwdViaJac(int nfwd);
    virtual bool adjViaJac(int nadj);
    ///@}

    /// Input and output buffers
    std::vector<DM> ibuf_, obuf_;

    /// Input and output scheme
    std::vector<std::string> ischeme_, oscheme_;

    /** \brief  Verbose -- for debugging purposes */
    bool verbose_;

    /** \brief  Use just-in-time compiler */
    bool jit_;
    eval_t eval_;

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

    /** \brief Get type name */
    virtual std::string type_name() const;

    /** \brief Check if the function is of a particular type */
    virtual bool is_a(const std::string& type, bool recursive) const;

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

    /** \brief Allocate memory block */
    virtual void* alloc_mem() {return 0;}

    /** \brief Free allocated memory block */
    virtual void free_mem(void* mem) {}

    ///@{
    /// Linear solver specific (cf. Linsol class)
    virtual void linsol_prepare(const double** arg, double** res, int* iw, double* w, void* mem);
    virtual void linsol_solve(bool tr);
    virtual void linsol_solve(double* x, int nrhs, bool tr);
    virtual MX linsol_solve(const MX& A, const MX& B, bool tr);
    virtual void linsol_spsolve(bvec_t* X, const bvec_t* B, bool tr) const;
    virtual void linsol_spsolve(DM& X, const DM& B, bool tr) const;
    virtual void linsol_solveL(double* x, int nrhs, bool tr);
    virtual Sparsity linsol_cholesky_sparsity(bool tr) const;
    virtual DM linsol_cholesky(bool tr) const;
    virtual void linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem,
                               bool tr, int nrhs);
    virtual void linsol_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                                const std::vector<std::vector<MX> >& fseed,
                                std::vector<std::vector<MX> >& fsens, bool tr);
    virtual void linsol_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                                const std::vector<std::vector<MX> >& aseed,
                                std::vector<std::vector<MX> >& asens, bool tr);
    virtual void linsol_spFwd(const bvec_t** arg, bvec_t** res,
                              int* iw, bvec_t* w, void* mem, bool tr, int nrhs);
    virtual void linsol_spAdj(bvec_t** arg, bvec_t** res,
                              int* iw, bvec_t* w, void* mem, bool tr, int nrhs);
    ///@}

  protected:
    /** \brief  Temporary vector needed for the evaluation (integer) */
    std::vector<int> iw_tmp_;

    /** \brief  Temporary vector needed for the evaluation (real) */
    std::vector<double> w_tmp_;

  private:
    /** \brief Memory that is persistent during a call (but not between calls) */
    size_t sz_arg_per_, sz_res_per_, sz_iw_per_, sz_w_per_;

    /** \brief Temporary memory inside a function */
    size_t sz_arg_tmp_, sz_res_tmp_, sz_iw_tmp_, sz_w_tmp_;
  };

  // Template implementations
  template<typename MatType>
  bool FunctionInternal::purgable(const std::vector<MatType>& v) {
    for (typename std::vector<MatType>::const_iterator i=v.begin(); i!=v.end(); ++i) {
      if (!i->is_zero()) return false;
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

  template<typename M>
  void FunctionInternal::call(const std::vector<M>& arg, std::vector<M>& res,
                              bool always_inline, bool never_inline) {
    // Check if inputs need to be replaced
    if (!matchingArg(arg)) {
      return call(replaceArg(arg), res, always_inline, never_inline);
    }

    // Call the type-specific method
    _call(arg, res, always_inline, never_inline);
  }

  template<typename D>
  void FunctionInternal::_call(const std::vector<Matrix<D> >& arg, std::vector<Matrix<D> >& res,
                               bool always_inline, bool never_inline) {
    casadi_assert_message(!never_inline, "Call-nodes only possible in MX expressions");

    // Get the number of inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Check if matching input sparsity
    bool matching_sparsity = true;
    casadi_assert(arg.size()==n_in);
    for (int i=0; matching_sparsity && i<n_in; ++i)
      matching_sparsity = arg[i].sparsity()==sparsity_in(i);

    // Correct input sparsity if needed
    if (!matching_sparsity) {
      std::vector<Matrix<D> > arg2(arg);
      for (int i=0; i<n_in; ++i)
        if (arg2[i].sparsity()!=sparsity_in(i))
          arg2[i] = project(arg2[i], sparsity_in(i));
      return _call(arg2, res, always_inline, never_inline);
    }

    // Allocate results
    res.resize(n_out);
    for (int i=0; i<n_out; ++i)
      if (res[i].sparsity()!=sparsity_out(i))
        res[i] = Matrix<D>::zeros(sparsity_out(i));

    // Allocate temporary memory if needed
    iw_tmp_.resize(sz_iw());
    std::vector<D> w_tmp(sz_w());

    // Get pointers to input arguments
    std::vector<const D*> argp(sz_arg());
    for (int i=0; i<n_in; ++i) argp[i]=getPtr(arg[i]);

    // Get pointers to output arguments
    std::vector<D*> resp(sz_res());
    for (int i=0; i<n_out; ++i) resp[i]=getPtr(res[i]);

    // Call memory-less
    _eval(getPtr(argp), getPtr(resp), getPtr(iw_tmp_), getPtr(w_tmp), 0);
  }

  template<typename M>
  void FunctionInternal::
  forward(const std::vector<M>& arg, const std::vector<M>& res,
          const std::vector<std::vector<M> >& fseed,
          std::vector<std::vector<M> >& fsens,
          bool always_inline, bool never_inline) {
    // Check inputs
    checkArg(arg);
    checkRes(res);

    // Check forward mode seeds dimensions
    int n_in = this->n_in();
    for (int d=0; d<fseed.size(); ++d) {
      casadi_assert_message(fseed[d].size()==n_in,
                            "Incorrect number of forward seeds for direction " << d
                            << ": Expected " << n_in << ", got " << fseed[d].size());
      for (int i=0; i<n_in; ++i) {
        casadi_assert_message(checkMat(fseed[d][i].sparsity(), input(i).sparsity()),
                              "Forward seed " << i << " for direction " << d
                              << " has mismatching shape. Expected " << input(i).size()
                              << ", got " << fseed[d][i].size());
      }
    }

    // Check if forward mode seeds need to be replaced
    for (int d=0; d<fseed.size(); ++d) {
      if (!matchingArg(fseed[d])) {
        return forward(arg, res, replaceFwdSeed(fseed), fsens,
                       always_inline, never_inline);
      }
    }

    // Call the type-specific method
    _forward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  template<typename M>
  void FunctionInternal::
  reverse(const std::vector<M>& arg, const std::vector<M>& res,
          const std::vector<std::vector<M> >& aseed,
          std::vector<std::vector<M> >& asens,
          bool always_inline, bool never_inline) {
    // Check inputs
    checkArg(arg);
    checkRes(res);

    // Check reverse mode seeds dimensions
    int n_out = this->n_out();
    for (int d=0; d<aseed.size(); ++d) {
      casadi_assert_message(aseed[d].size()==n_out,
                            "Incorrect number of adjoint seeds for direction " << d
                            << ": Expected " << n_out << ", got " << aseed[d].size());
      for (int i=0; i<n_out; ++i) {
        casadi_assert_message(checkMat(aseed[d][i].sparsity(), output(i).sparsity()),
                              "Adjoint seed " << i << " for direction " << d
                              << " has mismatching shape. Expected " << output(i).size()
                              << ", got " << aseed[d][i].size());
      }
    }

    // Check if reverse mode seeds need to be replaced
    for (int d=0; d<aseed.size(); ++d) {
      if (!matchingRes(aseed[d])) {
        return reverse(arg, res, replaceAdjSeed(aseed), asens,
                       always_inline, never_inline);
      }
    }

    // Call the type-specific method
    _reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  template<typename M>
  void FunctionInternal::checkArg(const std::vector<M>& arg, bool hcat) const {
    int n_in = this->n_in();
    casadi_assert_message(arg.size()==n_in, "Incorrect number of inputs: Expected "
                          << n_in << ", got " << arg.size());
    for (int i=0; i<n_in; ++i) {
      casadi_assert_message(checkMat(arg[i].sparsity(), sparsity_in(i), hcat),
                            "Input " << i << " has mismatching shape. Expected "
                            << size_in(i) << ", got " << arg[i].size());
    }
  }

  template<typename M>
  void FunctionInternal::checkRes(const std::vector<M>& res) const {
    int n_out = this->n_out();
    casadi_assert_message(res.size()==n_out, "Incorrect number of outputs: Expected "
                          << n_out << ", got " << res.size());
    for (int i=0; i<n_out; ++i) {
      casadi_assert_message(checkMat(res[i].sparsity(), sparsity_out(i)),
                            "Output " << i << " has mismatching shape. Expected "
                            << size_out(i) << ", got " << res[i].size());
    }
  }

  template<typename M>
  bool FunctionInternal::matchingArg(const std::vector<M>& arg, bool hcat) const {
    checkArg(arg, hcat);
    int n_in = this->n_in();
    for (int i=0; i<n_in; ++i) {
      if (hcat) {
        if (arg.at(i).size1()!=size1_in(i)) return false;
        if (arg.at(i).size2() % size2_in(i)!=0 || arg.at(i).size2()==0) return false;
      } else {
        if (arg.at(i).size()!=size_in(i)) return false;
      }
    }
    return true;
  }

  template<typename M>
  bool FunctionInternal::matchingRes(const std::vector<M>& res) const {
    checkRes(res);
    int n_out = this->n_out();
    for (int i=0; i<n_out; ++i) {
      if (res.at(i).size()!=size_out(i)) return false;
    }
    return true;
  }

  template<typename M>
  M replaceMat(const M& arg, const Sparsity& inp, bool hcat=false) {
    if (arg.size()==inp.size()) {
      // Matching dimensions already
      return arg;
    } else if (hcat && arg.size1()==inp.size1() && arg.size2() % inp.size2()==0
               && arg.size2() >=0) {
      // Matching horzcat dimensions
      return arg;
    } else if (arg.is_empty()) {
      // Empty matrix means set zero
      return M(inp.size());
    } else if (arg.is_scalar()) {
      // Scalar assign means set all
      return M(inp, arg);
    } else {
      // Assign vector with transposing
      casadi_assert(arg.size1()==inp.size2() && arg.size2()==inp.size1()
                    && (arg.is_column() || inp.is_column()));
      return arg.T();
    }
  }

  template<typename M>
  std::vector<M> FunctionInternal::replaceArg(const std::vector<M>& arg, bool hcat) const {
    std::vector<M> r(arg.size());
    for (int i=0; i<r.size(); ++i) r[i] = replaceMat(arg[i], sparsity_in(i), hcat);
    return r;
  }

  template<typename M>
  std::vector<M> FunctionInternal::replaceRes(const std::vector<M>& res) const {
    std::vector<M> r(res.size());
    for (int i=0; i<r.size(); ++i) r[i] = replaceMat(res[i], sparsity_out(i));
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> >
  FunctionInternal::replaceFwdSeed(const std::vector<std::vector<M> >& fseed) const {
    std::vector<std::vector<M> > r(fseed.size());
    for (int d=0; d<r.size(); ++d) r[d] = replaceArg(fseed[d]);
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> >
  FunctionInternal::replaceAdjSeed(const std::vector<std::vector<M> >& aseed) const {
    std::vector<std::vector<M> > r(aseed.size());
    for (int d=0; d<r.size(); ++d) r[d] = replaceRes(aseed[d]);
    return r;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_FUNCTION_INTERNAL_HPP
