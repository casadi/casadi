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
#include <set>
#include <stack>
#include "code_generator.hpp"
#include "importer.hpp"
#include "sparse_storage.hpp"
#include "options.hpp"
#include "shared_object_internal.hpp"

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

  /// Combine two dictionaries, giving priority to first one
  Dict CASADI_EXPORT combine(const Dict& first, const Dict& second);

  /** \brief Internal class for Function
      \author Joel Andersson
      \date 2010-2015
  */
  class CASADI_EXPORT FunctionInternal : public SharedObjectInternal {
  public:
    /** \brief Constructor */
    FunctionInternal(const std::string& name);

    /** \brief  Destructor */
    ~FunctionInternal() override = 0;

    /** \brief  Obtain solver name from Adaptor */
    virtual std::string getAdaptorSolverName() const { return ""; }

    /** \brief Construct
        Prepares the function for evaluation
     */
    void construct(const Dict& opts);

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief Initialize
        Initialize and make the object ready for setting arguments and evaluation.
        This method is typically called after setting options but before evaluating.
        If passed to another class (in the constructor), this class should invoke
        this function when initialized. */
    virtual void init(const Dict& opts);

    /** \brief Finalize the object creation
        This function, which visits the class hierarchy in reverse order is run after
        init() has been completed.
    */
    virtual void finalize(const Dict& opts);

    /** \brief Get a public class instance */
    Function self() const { return shared_from_this<Function>();}

    // Factory
    virtual Function factory(const std::string& name,
                             const std::vector<std::string>& s_in,
                             const std::vector<std::string>& s_out,
                             const Function::AuxOut& aux,
                             const Dict& opts) const;

    // Get list of dependency functions
    virtual std::vector<std::string> get_function() const;

    // Get a dependency function
    virtual const Function& get_function(const std::string &name) const;

    // Check if a particular dependency exists
    virtual bool has_function(const std::string& fname) const {return false;}

    /** \brief Which variables enter with some order
    *
    * \param[in] order Only 1 (linear) and 2 (nonlinear) allowed
    * \param[in] tr   Flip the relationship. Return which expressions contain the variables
    */
    virtual std::vector<bool> which_depends(const std::string& s_in,
                                           const std::vector<std::string>& s_out,
                                           int order, bool tr=false) const;


    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(int i);
    virtual std::string get_name_out(int i);
    ///@}

    ///@{
    /** \brief  Is the class able to propagate seeds through the algorithm? */
    virtual bool has_spfwd() const { return false;}
    virtual bool has_sprev() const { return false;}
    ///@}

    ///@{
    /** \brief  Evaluate numerically */
    int _eval(const double** arg, double** res, int* iw, double* w, int mem) const;
    virtual int eval(const double** arg, double** res, int* iw, double* w, void* mem) const;
    ///@}

    /** \brief  Evaluate numerically, simplied syntax */
    virtual void simple(const double* arg, double* res) const;

    /** \brief  Evaluate with symbolic scalars */
    virtual int eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const;

    ///@{
    /** \brief  Evaluate with symbolic matrices */
    virtual void eval_mx(const MXVector& arg, MXVector& res,
                         bool always_inline, bool never_inline) const;
    ///@}

    ///@{
    /** \brief Evaluate a function, overloaded */
    int _eval(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const;
    int _eval(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const;
    ///@}

    ///@{
    /** \brief Call a function, overloaded */
    void _call(const MXVector& arg, MXVector& res,
               bool always_inline, bool never_inline) const {
      eval_mx(arg, res, always_inline, never_inline);
    }
    template<typename D>
    void _call(const std::vector<Matrix<D> >& arg, std::vector<Matrix<D> >& res,
               bool always_inline, bool never_inline) const;
    ///@}

    /** \brief Call a function, templated */
    template<typename M>
      void call(const std::vector<M>& arg, std::vector<M>& res,
               bool always_inline, bool never_inline) const;

    /** Helper function */
    static bool check_mat(const Sparsity& arg, const Sparsity& inp);

    /** \brief Check if input arguments have correct length and dimensions
     */
    template<typename M>
      void check_arg(const std::vector<M>& arg) const;

    /** \brief Check if output arguments have correct length and dimensions */
    template<typename M>
      void check_res(const std::vector<M>& res) const;

    /** \brief Check if input arguments that needs to be replaced
     */
    template<typename M>
      bool matching_arg(const std::vector<M>& arg) const;

    /** \brief Check if output arguments that needs to be replaced */
    template<typename M>
      bool matching_res(const std::vector<M>& arg) const;

    /** \brief Replace 0-by-0 inputs
     */
    template<typename M>
      std::vector<M> replace_arg(const std::vector<M>& arg) const;

    /** \brief Replace 0-by-0 outputs */
    template<typename M>
      std::vector<M> replace_res(const std::vector<M>& res) const;

    /** \brief Replace 0-by-0 forward seeds */
    template<typename M>
      std::vector<std::vector<M> >
      replace_fseed(const std::vector<std::vector<M> >& fseed) const;

    /** \brief Replace 0-by-0 reverse seeds */
    template<typename M>
      std::vector<std::vector<M> >
      replace_aseed(const std::vector<std::vector<M> >& aseed) const;

    ///@{
    /** \brief Forward mode AD, virtual functions overloaded in derived classes */
    virtual void call_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                            const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens,
                            bool always_inline, bool never_inline) const;
    virtual void call_forward(const std::vector<SX>& arg, const std::vector<SX>& res,
                            const std::vector<std::vector<SX> >& fseed,
                            std::vector<std::vector<SX> >& fsens,
                            bool always_inline, bool never_inline) const;
    ///@}

    ///@{
    /** \brief Reverse mode, virtual functions overloaded in derived classes */
    virtual void call_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                            const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens,
                            bool always_inline, bool never_inline) const;
    virtual void call_reverse(const std::vector<SX>& arg, const std::vector<SX>& res,
                            const std::vector<std::vector<SX> >& aseed,
                            std::vector<std::vector<SX> >& asens,
                            bool always_inline, bool never_inline) const;
    ///@}

    /** \brief Parallel evaluation */
    std::vector<MX> mapsum_mx(const std::vector<MX > &arg, const std::string& parallelization);

    /** \brief Do the derivative functions need nondifferentiated outputs? */
    virtual bool uses_output() const {return false;}

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements */
    Function jacobian() const;
    virtual bool has_jacobian() const { return false;}
    virtual Function get_jacobian(const std::string& name,
                                  const std::vector<std::string>& inames,
                                  const std::vector<std::string>& onames,
                                  const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return function that calculates forward derivatives
     *    forward(nfwd) returns a cached instance if available,
     *    and calls <tt>Function get_forward(int nfwd)</tt>
     *    if no cached version is available.
     */
    Function forward(int nfwd) const;
    virtual bool has_forward(int nfwd) const { return false;}
    virtual Function get_forward(int nfwd, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    /// \brief Get directional derivatives using finite differencing
    virtual Function get_fd(int nfwd, const std::string& name,
                            const std::vector<std::string>& inames,
                            const std::vector<std::string>& onames,
                            const Dict& opts) const;

    ///@{
    /** \brief Return function that calculates adjoint derivatives
     *    reverse(nadj) returns a cached instance if available,
     *    and calls <tt>Function get_reverse(int nadj)</tt>
     *    if no cached version is available.
     */
    Function reverse(int nadj) const;
    virtual bool has_reverse(int nadj) const { return false;}
    virtual Function get_reverse(int nadj, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    /** \brief returns a new function with a selection of inputs/outputs of the original */
    virtual Function slice(const std::string& name, const std::vector<int>& order_in,
                           const std::vector<int>& order_out, const Dict& opts) const;

    /** \brief Get oracle */
    virtual const Function& oracle() const;

    /** \brief Can derivatives be calculated in any way? */
    bool has_derivative() const;

    /** \brief  Weighting factor for chosing forward/reverse mode */
    virtual double ad_weight() const;

    /** \brief  Weighting factor for chosing forward/reverse mode,
        sparsity propagation */
    virtual double sp_weight() const;

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
    virtual std::vector<SX> free_sx() const;

    /** \brief Does the function have free variables */
    virtual bool has_free() const { return false;}

    /** \brief Extract the functions needed for the Lifted Newton method */
    virtual void generate_lifted(Function& vdef_fcn, Function& vinit_fcn) const;

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
    virtual int n_nodes() const;

    /** \brief Wrap in an Function instance consisting of only one MX call */
    Function wrap() const;

    /** \brief Generate code the function */
    virtual void codegen(CodeGenerator& g, const std::string& fname,
                                  bool decl_static) const;

    /** \brief Generate meta-information allowing a user to evaluate a generated function */
    void codegen_meta(CodeGenerator& g, const std::string& fname) const;

    /** \brief Use simplified signature */
    virtual bool simplified_call() const { return false;}

    /** \brief Get name of the evaluation function */
    std::string eval_name() const;

    /** \brief Get name in codegen */
    virtual std::string codegen_name(const CodeGenerator& g) const;

    /** \brief Add a dependent function */
    virtual void add_dependency(CodeGenerator& g) const;

    /** \brief Codegen incref for dependencies */
    virtual void codegen_incref(CodeGenerator& g) const {}

    /** \brief Codegen decref for dependencies */
    virtual void codegen_decref(CodeGenerator& g) const {}

    /** \brief Code generate the function  */
    std::string signature(const std::string& fname) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void codegen_declarations(CodeGenerator& g) const;

    /** \brief Generate code for the function body */
    virtual void codegen_body(CodeGenerator& g) const;

    /** \brief Export / Generate C code for the dependency function */
    virtual std::string generate_dependencies(const std::string& fname, const Dict& opts) const;

    /** \brief Is codegen supported? */
    virtual bool has_codegen() const { return false;}

    /** \brief Jit dependencies */
    virtual void jit_dependencies(const std::string& fname) {}

    /** \brief  Print */
    void print_long(std::ostream &stream) const override;

    /** \brief  Print */
    void print_short(std::ostream &stream) const override;

    /** \brief Get function signature: name:(inputs)->(outputs) */
    std::string definition() const;

    /** \brief Check if the numerical values of the supplied bounds make sense */
    virtual void checkInputs() const {}

    /** \brief Print dimensions of inputs and outputs */
    void print_dimensions(std::ostream &stream) const;

    /** \brief Print list of options */
    void print_options(std::ostream &stream) const;

    /** \brief Print all information there is to know about a certain option */
    void print_option(const std::string &name, std::ostream &stream) const;

    /** \brief Print free variables */
    virtual void print_free(std::ostream &stream) const;

    /** \brief Get the unidirectional or bidirectional partition */
    void get_partition(int iind, int oind, Sparsity& D1, Sparsity& D2,
                      bool compact, bool symmetric,
                      bool allow_forward, bool allow_reverse) const;

    /// Verbose mode?
    bool verbose() const;

    ///@{
    /** \brief Number of function inputs and outputs */
    inline int n_in() const { return isp_.size();}
    virtual size_t get_n_in();
    inline int n_out() const { return osp_.size();}
    virtual size_t get_n_out();
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
    int numel_in(int ind) const { return sparsity_in(ind).numel(); }
    int numel_out(int ind) const { return sparsity_out(ind).numel(); }
    int numel_out() const;
    ///@}

    ///@{
    /** \brief Input/output dimensions */
    int size1_in(int ind) const { return sparsity_in(ind).size1(); }
    int size2_in(int ind) const { return sparsity_in(ind).size2(); }
    int size1_out(int ind) const { return sparsity_out(ind).size1(); }
    int size2_out(int ind) const { return sparsity_out(ind).size2(); }
    std::pair<int, int> size_in(int ind) const { return sparsity_in(ind).size(); }
    std::pair<int, int> size_out(int ind) const { return sparsity_out(ind).size(); }
    ///@}

    ///@{
    /** \brief Are all inputs and outputs scalar */
    bool all_scalar() const;

    /** \brief Name of the function */
    const std::string& name() const { return name_;}

    /// Generate the sparsity of a Jacobian block
    virtual Sparsity getJacSparsity(int iind, int oind, bool symmetric) const;

    /// Get the sparsity pattern, forward mode
    template<bool fwd>
    Sparsity getJacSparsityGen(int iind, int oind, bool symmetric, int gr_i=1, int gr_o=1) const;

    /// A flavor of getJacSparsity that does hierarchical block structure recognition
    Sparsity getJacSparsityHierarchical(int iind, int oind) const;

    /** A flavor of getJacSparsity that does hierarchical block
    * structure recognition for symmetric Jacobians
    */
    Sparsity getJacSparsityHierarchicalSymm(int iind, int oind) const;

    /// Get, if necessary generate, the sparsity of a Jacobian block
    Sparsity& sparsity_jac(int iind, int oind, bool compact, bool symmetric) const;

    /// Get a vector of symbolic variables corresponding to the outputs
    virtual std::vector<MX> symbolic_output(const std::vector<MX>& arg) const;

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
      return ischeme_.at(ind);
    }

    /** \brief Get output scheme name by index */
    virtual std::string name_out(int ind) const {
      return oscheme_.at(ind);
    }

    /** \brief Get default input value */
    virtual double default_in(int ind) const {
      return 0;
    }

    /** \brief Get largest input value */
    virtual double max_in(int ind) const {
      return inf;
    }

    /** \brief Get smallest input value */
    virtual double min_in(int ind) const {
      return -inf;
    }

    /** \brief Get sparsity of a given input */
    /// @{
    inline const Sparsity& sparsity_in(int ind) const {
      return isp_.at(ind);
    }
    inline const Sparsity& sparsity_in(const std::string& iname) const {
      return sparsity_in(index_in(iname));
    }
    virtual Sparsity get_sparsity_in(int i);
    /// @}

    /** \brief Get sparsity of a given output */
    /// @{
    inline const Sparsity& sparsity_out(int ind) const {
      return osp_.at(ind);
    }
    inline const Sparsity& sparsity_out(const std::string& iname) const {
      return sparsity_out(index_out(iname));
    }
    virtual Sparsity get_sparsity_out(int i);
    /// @}

    /** \brief  Log the status of the solver */
    void log(const std::string& msg) const;

    /** \brief  Log the status of the solver, function given */
    void log(const std::string& fcn, const std::string& msg) const;

    /** \brief  Propagate sparsity forward */
    virtual int sp_forward(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const;

    /** \brief  Propagate sparsity backwards */
    virtual void sp_reverse(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const;

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

    /// Memory objects
    void* memory(int ind) const;

    /** \brief Create memory block */
    virtual void* alloc_memory() const {return 0;}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const {}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const;

    /** \brief Maximum number of memory objects */
    virtual int n_mem() const { return 0;}

    /** \brief Clear all memory (called from destructor) */
    void clear_memory();

    ///@{
    /// Get all statistics
    virtual Dict get_stats(void* mem) const { return Dict();}
    virtual Dict _get_stats(int mem) const { return get_stats(memory(mem));}
    ///@}

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const {}

    /** \brief Set the (temporary) work vectors */
    virtual void set_temp(void* mem, const double** arg, double** res,
                          int* iw, double* w) const {}

    /** \brief Set the (persistent and temporary) work vectors */
    void setup(void* mem, const double** arg, double** res, int* iw, double* w) const;

    ///@{
    /** \brief Calculate derivatives by multiplying the full Jacobian and multiplying */
    virtual bool fwdViaJac(int nfwd) const;
    virtual bool adjViaJac(int nadj) const;
    ///@}

    /// Checkout a memory object
    int checkout() const;

    /// Release a memory object
    void release(int mem) const;

    /// Input and output sparsity
    std::vector<Sparsity> isp_, osp_;

    /// Input and output scheme
    std::vector<std::string> ischeme_, oscheme_;

    /** \brief  Verbose -- for debugging purposes */
    bool verbose_;

    /** \brief  Use just-in-time compiler */
    bool jit_;

    /** \brief Numerical evaluation redirected to a C function */
    eval_t eval_;
    simple_t simple_;

    /** \brief Dict of statistics (resulting from evaluate) */
    Dict stats_;

    /** \brief Reference counting in codegen? */
    bool has_refcount_;

    /// Cache for functions to evaluate directional derivatives
    mutable std::vector<WeakRef> forward_, reverse_;

    /// Cache for full Jacobian
    mutable WeakRef jacobian_;

    /// Cache for sparsities of the Jacobian blocks
    mutable SparseStorage<Sparsity> jac_sparsity_, jac_sparsity_compact_;

    /// If the function is the derivative of another function
    Function derivative_of_;

    /// Wrapper function for indirect derivative calculation
    mutable WeakRef wrap_;

    /// User-set field
    void* user_data_;

    /// Name
    std::string name_;

    /// Just-in-time compiler
    std::string compilerplugin_;
    Importer compiler_;
    Dict jit_options_;

    /// Penalty factor for using a complete Jacobian to calculate directional derivatives
    double jac_penalty_;

    // Types of derivative calculation permitted
    bool enable_forward_, enable_reverse_, enable_jacobian_, enable_fd_;

    /// Weighting factor for derivative calculation and sparsity pattern calculation
    double ad_weight_, ad_weight_sp_;

    /// Maximum number of sensitivity directions
    int max_num_dir_;

    /// Errors are thrown when NaN is produced
    bool regularity_check_;

    /// Errors are thrown if numerical values of inputs look bad
    bool inputs_check_;

    // Print timing statistics
    bool print_time_;

    // Finite difference step
    Dict fd_options_;

    /** \brief Get type name */
    virtual std::string type_name() const = 0;

    /** \brief Check if the function is of a particular type */
    virtual bool is_a(const std::string& type, bool recursive) const;

    /** \brief Can a derivative direction be skipped */
    template<typename MatType>
    static bool purgable(const std::vector<MatType>& seed);

    /** \brief Symbolic expressions for the forward seeds */
    template<typename MatType>
    std::vector<std::vector<MatType> >
    fwd_seed(int nfwd) const;

    /** \brief Symbolic expressions for the adjoint seeds */
    template<typename MatType>
    std::vector<std::vector<MatType> >
    symbolicAdjSeed(int nadj, const std::vector<MatType>& v) const;

  protected:
    static void print_stats_line(int maxNameLen, std::string label, double n_call,
      double t_proc, double t_wall);
  private:
    /// Memory objects
    mutable std::vector<void*> mem_;

    /// Unused memory objects
    mutable std::stack<int> unused_;

    /** \brief Memory that is persistent during a call (but not between calls) */
    size_t sz_arg_per_, sz_res_per_, sz_iw_per_, sz_w_per_;

    /** \brief Temporary memory inside a function */
    size_t sz_arg_tmp_, sz_res_tmp_, sz_iw_tmp_, sz_w_tmp_;
  };

  // Template implementations
  template<typename MatType>
  bool FunctionInternal::purgable(const std::vector<MatType>& v) {
    for (auto i=v.begin(); i!=v.end(); ++i) {
      if (!i->is_zero()) return false;
    }
    return true;
  }

  template<typename MatType>
  std::vector<std::vector<MatType> >
  FunctionInternal::
  fwd_seed(int nfwd) const {
    std::vector<std::vector<MatType>> fseed(nfwd);
    int n_in = this->n_in();
    for (int dir=0; dir<nfwd; ++dir) {
      fseed[dir].resize(n_in);
      for (int iind=0; iind<n_in; ++iind) {
        std::string n = "f" + to_string(dir) + "_" +  name_in(iind);
        fseed[dir][iind] = MatType::sym(n, sparsity_in(iind));
      }
    }
    return fseed;
  }

  template<typename MatType>
  std::vector<std::vector<MatType> >
  FunctionInternal::
  symbolicAdjSeed(int nadj, const std::vector<MatType>& v) const {
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
                              bool always_inline, bool never_inline) const {
    // If all inputs are scalar ...
    if (all_scalar()) {
      // ... and some arguments are matrix-valued with matching dimensions ...
      bool matrix_call = false;
      std::pair<int, int> sz;
      for (auto&& a : arg) {
        if (!a.is_scalar() && !a.is_empty()) {
          if (!matrix_call) {
            // Matrix call
            matrix_call = true;
            sz = a.size();
          } else if (a.size()!=sz) {
            // Not same dimensions
            matrix_call = false;
            break;
          }
        }
      }

      // ... then, call multiple times
      if (matrix_call) {
        // Start with zeros
        res.resize(n_out());
        M z = M::zeros(sz);
        for (auto&& a : res) a = z;
        // Call multiple times
        std::vector<M> arg1 = arg, res1;
        for (int c=0; c<sz.second; ++c) {
          for (int r=0; r<sz.first; ++r) {
            // Get scalar arguments
            for (int i=0; i<arg.size(); ++i) {
              if (arg[i].size()==sz) arg1[i] = arg[i](r, c);
            }
            // Call recursively with scalar arguments
            call(arg1, res1, always_inline, never_inline);
            // Get results
            casadi_assert(res.size() == res1.size());
            for (int i=0; i<res.size(); ++i) res[i](r, c) = res1[i];
          }
        }
        // All elements assigned
        return;
      }
    }

    // Check if inputs need to be replaced
    if (!matching_arg(arg)) {
      return call(replace_arg(arg), res, always_inline, never_inline);
    }

    // Call the type-specific method
    _call(arg, res, always_inline, never_inline);
  }

  template<typename D>
  void FunctionInternal::_call(const std::vector<Matrix<D> >& arg, std::vector<Matrix<D> >& res,
                               bool always_inline, bool never_inline) const {
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
    for (int i=0; i<n_out; ++i) {
      if (res[i].sparsity()!=sparsity_out(i)) {
        res[i] = Matrix<D>::zeros(sparsity_out(i));
      }
    }

    // Allocate temporary memory if needed
    std::vector<int> iw_tmp(sz_iw());
    std::vector<D> w_tmp(sz_w());

    // Get pointers to input arguments
    std::vector<const D*> argp(sz_arg());
    for (int i=0; i<n_in; ++i) argp[i]=get_ptr(arg[i]);

    // Get pointers to output arguments
    std::vector<D*> resp(sz_res());
    for (int i=0; i<n_out; ++i) resp[i]=get_ptr(res[i]);

    // Call memory-less
    (void)_eval(get_ptr(argp), get_ptr(resp), get_ptr(iw_tmp), get_ptr(w_tmp), 0);
  }

  template<typename M>
  void FunctionInternal::check_arg(const std::vector<M>& arg) const {
    int n_in = this->n_in();
    casadi_assert_message(arg.size()==n_in, "Incorrect number of inputs: Expected "
                          << n_in << ", got " << arg.size());
    for (int i=0; i<n_in; ++i) {
      casadi_assert_message(check_mat(arg[i].sparsity(), sparsity_in(i)),
                            "Input " << i << " (" << name_in(i) << ") has mismatching shape. "
                            << "Expected " << size_in(i) << ", got " << arg[i].size());
    }
  }

  template<typename M>
  void FunctionInternal::check_res(const std::vector<M>& res) const {
    int n_out = this->n_out();
    casadi_assert_message(res.size()==n_out, "Incorrect number of outputs: Expected "
                          << n_out << ", got " << res.size());
    for (int i=0; i<n_out; ++i) {
      casadi_assert_message(check_mat(res[i].sparsity(), sparsity_out(i)),
                            "Output " << i << " (" << name_out(i) << ") has mismatching shape. "
                            "Expected " << size_out(i) << ", got " << res[i].size());
    }
  }

  template<typename M>
  bool FunctionInternal::matching_arg(const std::vector<M>& arg) const {
    check_arg(arg);
    int n_in = this->n_in();
    for (int i=0; i<n_in; ++i) {
      if (arg.at(i).size()!=size_in(i)) return false;
    }
    return true;
  }

  template<typename M>
  bool FunctionInternal::matching_res(const std::vector<M>& res) const {
    check_res(res);
    int n_out = this->n_out();
    for (int i=0; i<n_out; ++i) {
      if (res.at(i).size()!=size_out(i)) return false;
    }
    return true;
  }

  template<typename M>
  M replace_mat(const M& arg, const Sparsity& inp) {
    if (arg.size()==inp.size()) {
      // Matching dimensions already
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
  std::vector<M> FunctionInternal::replace_arg(const std::vector<M>& arg) const {
    std::vector<M> r(arg.size());
    for (int i=0; i<r.size(); ++i) r[i] = replace_mat(arg[i], sparsity_in(i));
    return r;
  }

  template<typename M>
  std::vector<M> FunctionInternal::replace_res(const std::vector<M>& res) const {
    std::vector<M> r(res.size());
    for (int i=0; i<r.size(); ++i) r[i] = replace_mat(res[i], sparsity_out(i));
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> >
  FunctionInternal::replace_fseed(const std::vector<std::vector<M> >& fseed) const {
    std::vector<std::vector<M> > r(fseed.size());
    for (int d=0; d<r.size(); ++d) r[d] = replace_arg(fseed[d]);
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> >
  FunctionInternal::replace_aseed(const std::vector<std::vector<M> >& aseed) const {
    std::vector<std::vector<M> > r(aseed.size());
    for (int d=0; d<r.size(); ++d) r[d] = replace_res(aseed[d]);
    return r;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_FUNCTION_INTERNAL_HPP
