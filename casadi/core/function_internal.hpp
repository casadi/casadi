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


#ifndef CASADI_FUNCTION_INTERNAL_HPP
#define CASADI_FUNCTION_INTERNAL_HPP

#include "function.hpp"
#include <set>
#include <stack>
#include "code_generator.hpp"
#include "importer.hpp"
#include "options.hpp"
#include "shared_object_internal.hpp"
#include "timing.hpp"
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

// This macro is for documentation purposes
#define INPUTSCHEME(name)

// This macro is for documentation purposes
#define OUTPUTSCHEME(name)

/// \cond INTERNAL

namespace casadi {
  template<typename T>
  std::vector<std::pair<std::string, T>> zip(const std::vector<std::string>& id,
                                             const std::vector<T>& mat) {
    casadi_assert_dev(id.size()==mat.size());
    std::vector<std::pair<std::string, T>> r(id.size());
    for (casadi_uint i=0; i<r.size(); ++i) r[i] = std::make_pair(id[i], mat[i]);
    return r;
  }

  /** \brief Function memory with temporary work vectors

      \identifier{ja} */
  struct CASADI_EXPORT ProtoFunctionMemory {
    // Function specific statistics
    std::map<std::string, FStats> fstats;

    // Short-hand for "total" fstats
    FStats* t_total;

    // Add a statistic
    void add_stat(const std::string& s) {
      bool added = fstats.insert(std::make_pair(s, FStats())).second;
      casadi_assert(added, "Duplicate stat: '" + s + "'");
    }
  };

  /** \brief Function memory with temporary work vectors

      \identifier{jb} */
  struct CASADI_EXPORT FunctionMemory : public ProtoFunctionMemory {
  };

  /** \brief Base class for FunctionInternal and LinsolInternal

    \author Joel Andersson
    \date 2017

      \identifier{jc} */
  class CASADI_EXPORT ProtoFunction : public SharedObjectInternal {
  public:
    /** \brief Constructor

        \identifier{jd} */
    ProtoFunction(const std::string& name);

    /** \brief  Destructor

        \identifier{je} */
    ~ProtoFunction() override = 0;

    /** \brief Construct

        Prepares the function for evaluation

        \identifier{jf} */
    void construct(const Dict& opts);

    ///@{
    /** \brief Options

        \identifier{jg} */
    static const Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /// Reconstruct options dict
    virtual Dict generate_options(const std::string& target) const;

    /** \brief Print list of options

        \identifier{jh} */
    void print_options(std::ostream &stream) const;

    /** \brief Print all information there is to know about a certain option

        \identifier{ji} */
    void print_option(const std::string &name, std::ostream &stream) const;

    /** \brief Does a particular option exist

        \identifier{jj} */
    bool has_option(const std::string &option_name) const;

    /** \brief Change option after object creation for debugging

        \identifier{jk} */
    virtual void change_option(const std::string& option_name, const GenericType& option_value);

    /** \brief Initialize

        Initialize and make the object ready for setting arguments and evaluation.
        This method is typically called after setting options but before evaluating.
        If passed to another class (in the constructor), this class should invoke
        this function when initialized.

        \identifier{jl} */
    virtual void init(const Dict& opts);

    /** \brief Finalize the object creation

        This function, which visits the class hierarchy in reverse order is run after
        init() has been completed.

        \identifier{jm} */
    virtual void finalize();

    /// Checkout a memory object
    int checkout() const;

    /// Release a memory object
    void release(int mem) const;

    /// Memory objects
    void* memory(int ind) const;

    /** \brief Create memory block

        \identifier{jn} */
    virtual void* alloc_mem() const { return new ProtoFunctionMemory(); }

    /** \brief Initalize memory block

        \identifier{jo} */
    virtual int init_mem(void* mem) const;

    /** \brief Free memory block

        \identifier{jp} */
    virtual void free_mem(void *mem) const { delete static_cast<ProtoFunctionMemory*>(mem); }

    /// Get all statistics
    virtual Dict get_stats(void* mem) const;

    /** \brief Clear all memory (called from destructor)

        \identifier{jq} */
    void clear_mem();

    /** \brief C-style formatted printing during evaluation

        \identifier{jr} */
    void print(const char* fmt, ...) const;

    /** \brief C-style formatted printing to string

        \identifier{js} */
    void sprint(char* buf, size_t buf_sz, const char* fmt, ...) const;

    /** \brief Format time in a fixed width 8 format

        \identifier{jt} */
    void format_time(char* buffer, double time) const;

    /** \brief Print timing statistics

        \identifier{ju} */
    void print_time(const std::map<std::string, FStats>& fstats) const;

    /** \brief Serialize an object

        \identifier{jv} */
    void serialize(SerializingStream &s) const;

    /** \brief Serialize an object without type information

        \identifier{jw} */
    virtual void serialize_body(SerializingStream &s) const;
    /** \brief Serialize type information

        \identifier{jx} */
    virtual void serialize_type(SerializingStream &s) const {}

    /** \brief String used to identify the immediate FunctionInternal subclass

        \identifier{jy} */
    virtual std::string serialize_base_function() const {
      return class_name();
    }

    /// Name
    std::string name_;

    /// Verbose printout
    bool verbose_;

    // Print timing statistics
    bool print_time_;

    // Print timing statistics
    bool record_time_;

    /// Errors are thrown when NaN is produced
    bool regularity_check_;

    /// Throw an exception on failure?
    bool error_on_fail_;

  protected:
    /** \brief Deserializing constructor

        \identifier{jz} */
    explicit ProtoFunction(DeserializingStream& s);

#ifdef CASADI_WITH_THREAD
    /// Mutex for thread safety
    mutable std::mutex mtx_;
#endif // CASADI_WITH_THREAD

  private:
    /// Memory objects
    mutable std::vector<void*> mem_;

    /// Unused memory objects
    mutable std::stack<int> unused_;
  };

  /** \brief Internal class for Function

      \author Joel Andersson
      \date 2010-2015

      \identifier{k0} */
  class CASADI_EXPORT FunctionInternal : public ProtoFunction {
    friend class Function;
  public:
    /** \brief Constructor

        \identifier{k1} */
    FunctionInternal(const std::string& name);

    /** \brief  Destructor

        \identifier{k2} */
    ~FunctionInternal() override = 0;

    /** \brief  Obtain solver name from Adaptor

        \identifier{k3} */
    virtual std::string getAdaptorSolverName() const { return ""; }

    ///@{
    /** \brief Options

        \identifier{k4} */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Reconstruct options dict
    Dict generate_options(const std::string& target) const override;

    /** \brief Change option after object creation for debugging

        \identifier{k5} */
    void change_option(const std::string& option_name, const GenericType& option_value) override;

    /** \brief Initialize

        \identifier{k6} */
    void init(const Dict& opts) override;

    /** \brief Finalize the object creation

        \identifier{k7} */
    void finalize() override;

    /** \brief Get a public class instance

        \identifier{k8} */
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

    // Add embedded function to map, helper function
    void add_embedded(std::map<FunctionInternal*, Function>& all_fun,
      const Function& dep, casadi_int max_depth) const;

    // Get all embedded functions, recursively
    virtual void find(std::map<FunctionInternal*, Function>& all_fun, casadi_int max_depth) const {}

    /** \brief Which variables enter with some order

    * \param[in] s_in Input name
    * \param[in] s_out Output name(s)
    * \param[in] order Only 1 (linear) and 2 (nonlinear) allowed
    * \param[in] tr   Flip the relationship. Return which expressions contain the variables

        \identifier{k9} */
    virtual std::vector<bool> which_depends(const std::string& s_in,
                                           const std::vector<std::string>& s_out,
                                           casadi_int order, bool tr=false) const;

    ///@{
    /** \brief  Is the class able to propagate seeds through the algorithm?

        \identifier{ka} */
    virtual bool has_spfwd() const { return false;}
    virtual bool has_sprev() const { return false;}
    ///@}

    ///@{
    /** \brief  Evaluate numerically

        \identifier{kb} */
    int eval_gen(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const;
    virtual int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const;
    ///@}

    /** \brief  Evaluate with symbolic scalars

        \identifier{kc} */
    virtual int eval_sx(const SXElem** arg, SXElem** res,
      casadi_int* iw, SXElem* w, void* mem) const;

    /** \brief  Evaluate with symbolic matrices

        \identifier{kd} */
    virtual void eval_mx(const MXVector& arg, MXVector& res,
                         bool always_inline, bool never_inline) const;

    ///@{
    /** \brief Evaluate with DM matrices

        \identifier{ke} */
    virtual std::vector<DM> eval_dm(const std::vector<DM>& arg) const;
    virtual bool has_eval_dm() const { return false;}
    ///@}

    ///@{
    /** \brief Evaluate a function, overloaded

        \identifier{kf} */
    int eval_gen(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w, void* mem) const {
      return eval_sx(arg, res, iw, w, mem);
    }
    int eval_gen(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
      return sp_forward(arg, res, iw, w, mem);
    }
    ///@}

    ///@{
    /** \brief Call a function, overloaded

        \identifier{kg} */
    void call_gen(const MXVector& arg, MXVector& res, casadi_int npar,
                  bool always_inline, bool never_inline) const;

    template<typename D>
    void call_gen(const std::vector<Matrix<D> >& arg, std::vector<Matrix<D> >& res,
                  casadi_int npar, bool always_inline, bool never_inline) const;
    ///@}

    /** \brief Call a function, templated

        \identifier{kh} */
    template<typename M>
      void call(const std::vector<M>& arg, std::vector<M>& res,
               bool always_inline, bool never_inline) const;

    ///@{
    /** Helper function
     *
     * \param npar[in] normal usage: 1, disallow pararallel calls: -1
     * \param npar[out] required number of parallel calls (or -1)
     */
    static bool check_mat(const Sparsity& arg, const Sparsity& inp, casadi_int& npar);
    ///@}

    ///@{
    /** \brief Check if input arguments have correct length and dimensions
     *
     * Raises errors.
     *
     * \param npar[in] normal usage: 1, disallow pararallel calls: -1
     * \param[out] npar: max number of horizontal repetitions across all arguments (or -1)

        \identifier{ki} */
    template<typename M>
    void check_arg(const std::vector<M>& arg, casadi_int& npar) const;
    ///@}

    ///@{
    /** \brief Check if output arguments have correct length and dimensions
     *
     * Raises errors.
     *
     * \param npar[in] normal usage: 1, disallow pararallel calls: -1
     * \param[out] npar: max number of horizontal repetitions across all arguments  (or -1)

        \identifier{kj} */
    template<typename M>
    void check_res(const std::vector<M>& res, casadi_int& npar) const;
    ///@}

    /** \brief Check if input arguments that needs to be replaced
     *
     * Raises errors
     *
     * \param npar[in] normal usage: 1, disallow pararallel calls: -1
     * \param[out] npar: max number of horizontal repetitions across all arguments  (or -1)

        \identifier{kk} */
    template<typename M> bool
    matching_arg(const std::vector<M>& arg, casadi_int& npar) const;

    /** \brief Check if output arguments that needs to be replaced
     *
     * Raises errors
     *
     * \param npar[in] normal usage: 1, disallow pararallel calls: -1
     * \param[out] npar: max number of horizontal repetitions across all arguments  (or -1)

        \identifier{kl} */
    template<typename M> bool
    matching_res(const std::vector<M>& arg, casadi_int& npar) const;

    /** \brief Replace 0-by-0 inputs

        \identifier{km} */
    template<typename M> std::vector<M>
    replace_arg(const std::vector<M>& arg, casadi_int npar) const;

    /** \brief Project sparsities
     *
        \identifier{kn} */
    template<typename M> std::vector<M>
    project_arg(const std::vector<M>& arg, casadi_int npar) const;

    /** \brief Project sparsities
     *
        \identifier{ko} */
    template<typename M> std::vector<M>
    project_res(const std::vector<M>& arg, casadi_int npar) const;

    /** \brief Replace 0-by-0 outputs

        \identifier{kp} */
    template<typename M> std::vector<M>
    replace_res(const std::vector<M>& res, casadi_int npar) const;

    /** \brief Replace 0-by-0 forward seeds

        \identifier{kq} */
    template<typename M> std::vector<std::vector<M>>
    replace_fseed(const std::vector<std::vector<M>>& fseed, casadi_int npar) const;

    /** \brief Replace 0-by-0 reverse seeds

        \identifier{kr} */
    template<typename M> std::vector<std::vector<M>>
    replace_aseed(const std::vector<std::vector<M>>& aseed, casadi_int npar) const;

    /** \brief Convert from/to input/output lists/map

        \identifier{ks} */
    /// @{
    template<typename M>
    std::map<std::string, M> convert_arg(const std::vector<M>& arg) const;
    template<typename M>
    std::vector<M> convert_arg(const std::map<std::string, M>& arg) const;
    template<typename M>
    std::map<std::string, M> convert_res(const std::vector<M>& res) const;
    template<typename M>
    std::vector<M> convert_res(const std::map<std::string, M>& res) const;
    /// @}

    /** \brief Convert from/to flat vector of input/output nonzeros

        \identifier{kt} */
    /// @{
    std::vector<double> nz_in(const std::vector<DM>& arg) const;
    std::vector<double> nz_out(const std::vector<DM>& res) const;
    std::vector<DM> nz_in(const std::vector<double>& arg) const;
    std::vector<DM> nz_out(const std::vector<double>& res) const;
    ///@}

    ///@{
    /** \brief Forward mode AD, virtual functions overloaded in derived classes

        \identifier{ku} */
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
    /** \brief Reverse mode, virtual functions overloaded in derived classes

        \identifier{kv} */
    virtual void call_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                            const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens,
                            bool always_inline, bool never_inline) const;
    virtual void call_reverse(const std::vector<SX>& arg, const std::vector<SX>& res,
                            const std::vector<std::vector<SX> >& aseed,
                            std::vector<std::vector<SX> >& asens,
                            bool always_inline, bool never_inline) const;
    ///@}

    /** \brief Parallel evaluation

        \identifier{kw} */
    std::vector<MX> mapsum_mx(const std::vector<MX > &arg, const std::string& parallelization);

    /** \brief Do the derivative functions need nondifferentiated outputs?

        \identifier{kx} */
    virtual bool uses_output() const {return false;}

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements

        \identifier{ky} */
    Function jacobian() const;
    virtual bool has_jacobian() const { return false;}
    virtual Function get_jacobian(const std::string& name,
                                  const std::vector<std::string>& inames,
                                  const std::vector<std::string>& onames,
                                  const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Get Jacobian sparsity

        \identifier{kz} */
    /// Get, if necessary generate, the sparsity of a Jacobian block
    Sparsity& jac_sparsity(casadi_int oind, casadi_int iind, bool compact, bool symmetric) const;
    virtual bool has_jac_sparsity(casadi_int oind, casadi_int iind) const { return false;}
    virtual Sparsity get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const;
    ///@}

    /// Helper function: Get name of forward derivative function
    static std::string forward_name(const std::string& fcn, casadi_int nfwd) {
      return "fwd" + str(nfwd) + "_" + fcn;
    }

    /// Determine prefix for differentiated functions
    std::string diff_prefix(const std::string& prefix) const;

    ///@{
    /** \brief Return function that calculates forward derivatives

     *    forward(nfwd) returns a cached instance if available,
     *    and calls <tt>Function get_forward(casadi_int nfwd)</tt>
     *    if no cached version is available.

        \identifier{l0} */
    Function forward(casadi_int nfwd) const;
    virtual bool has_forward(casadi_int nfwd) const { return false;}
    virtual Function get_forward(casadi_int nfwd, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    /// Helper function: Get name of adjoint derivative function
    static std::string reverse_name(const std::string& fcn, casadi_int nadj) {
      return "adj" + str(nadj) + "_" + fcn;
    }

    ///@{
    /** \brief Return function that calculates adjoint derivatives

     *    reverse(nadj) returns a cached instance if available,
     *    and calls <tt>Function get_reverse(casadi_int nadj)</tt>
     *    if no cached version is available.

        \identifier{l1} */
    Function reverse(casadi_int nadj) const;
    virtual bool has_reverse(casadi_int nadj) const { return false;}
    virtual Function get_reverse(casadi_int nadj, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    /** \brief Ensure that a matrix's sparsity is a horizontal multiple of another, or empty

        \identifier{26j} */
    template<typename MatType>
    static MatType ensure_stacked(const MatType& v, const Sparsity& sp, casadi_int n);

    /** \brief returns a new function with a selection of inputs/outputs of the original

        \identifier{l2} */
    virtual Function slice(const std::string& name, const std::vector<casadi_int>& order_in,
                           const std::vector<casadi_int>& order_out, const Dict& opts) const;

    /** \brief Get oracle

        \identifier{l3} */
    virtual const Function& oracle() const;

    /** \brief Can derivatives be calculated in any way?

        \identifier{l4} */
    bool has_derivative() const;

    /** \brief  Weighting factor for chosing forward/reverse mode

        \identifier{l5} */
    virtual double ad_weight() const;

    /** \brief  Weighting factor for chosing forward/reverse mode,

        sparsity propagation

        \identifier{l6} */
    virtual double sp_weight() const;

    ///@{
    /** \brief Get function input(s) and output(s)

        \identifier{l7} */
    virtual const SX sx_in(casadi_int ind) const;
    virtual const SX sx_out(casadi_int ind) const;
    virtual const std::vector<SX> sx_in() const;
    virtual const std::vector<SX> sx_out() const;
    virtual const MX mx_in(casadi_int ind) const;
    virtual const MX mx_out(casadi_int ind) const;
    virtual const std::vector<MX> mx_in() const;
    virtual const std::vector<MX> mx_out() const;
    const DM dm_in(casadi_int ind) const;
    const DM dm_out(casadi_int ind) const;
    const std::vector<DM> dm_in() const;
    const std::vector<DM> dm_out() const;
    ///@}

    /// Get free variables (MX)
    virtual std::vector<MX> free_mx() const;

    /// Get free variables (SX)
    virtual std::vector<SX> free_sx() const;

    /** \brief Does the function have free variables

        \identifier{l8} */
    virtual bool has_free() const { return false;}

    /** \brief Extract the functions needed for the Lifted Newton method

        \identifier{l9} */
    virtual void generate_lifted(Function& vdef_fcn, Function& vinit_fcn) const;

    /** \brief Get the number of atomic operations

        \identifier{la} */
    virtual casadi_int n_instructions() const;

    /** \brief Get an atomic operation operator index

        \identifier{lb} */
    virtual casadi_int instruction_id(casadi_int k) const;

    /** \brief Get the (integer) input arguments of an atomic operation

        \identifier{lc} */
    virtual std::vector<casadi_int> instruction_input(casadi_int k) const;

    /** \brief Get the floating point output argument of an atomic operation

        \identifier{ld} */
    virtual double instruction_constant(casadi_int k) const;

    /** \brief Get the (integer) output argument of an atomic operation

        \identifier{le} */
    virtual std::vector<casadi_int> instruction_output(casadi_int k) const;

    /** \brief Number of nodes in the algorithm

        \identifier{lf} */
    virtual casadi_int n_nodes() const;

    /** *\brief get MX expression associated with instruction

         \identifier{lg} */
    virtual MX instruction_MX(casadi_int k) const;

    /** *\brief get SX expression associated with instructions

         \identifier{lh} */
    virtual SX instructions_sx() const;

    /** \brief Wrap in an Function instance consisting of only one MX call

        \identifier{li} */
    Function wrap() const;

    /** \brief Wrap in an Function instance consisting of only one MX call

        \identifier{lj} */
    Function wrap_as_needed(const Dict& opts) const;

    /** \brief Get all functions in the cache

        \identifier{26g} */
    Dict cache() const;

    /** \brief Get function in cache

        \identifier{lk} */
    bool incache(const std::string& fname, Function& f, const std::string& suffix="") const;

    /** \brief Save function to cache

        \identifier{ll} */
    void tocache(const Function& f, const std::string& suffix="") const;

    /** \brief Generate code the function

        \identifier{lm} */
    void codegen(CodeGenerator& g, const std::string& fname) const;

    /** \brief Generate meta-information allowing a user to evaluate a generated function

        \identifier{ln} */
    void codegen_meta(CodeGenerator& g) const;

    /** \brief Codegen sparsities

        \identifier{lo} */
    void codegen_sparsities(CodeGenerator& g) const;

    /** \brief Get name in codegen

        \identifier{lp} */
    virtual std::string codegen_name(const CodeGenerator& g, bool ns=true) const;

    /** \brief Get thread-local memory object

        \identifier{lq} */
    std::string codegen_mem(CodeGenerator& g, const std::string& index="mem") const;

    /** \brief Codegen incref for dependencies

        \identifier{lr} */
    virtual void codegen_incref(CodeGenerator& g) const {}

    /** \brief Codegen decref for dependencies

        \identifier{ls} */
    virtual void codegen_decref(CodeGenerator& g) const {}

    /** \brief Codegen decref for alloc_mem

        \identifier{lt} */
    virtual void codegen_alloc_mem(CodeGenerator& g) const;

    /** \brief Codegen decref for init_mem

        \identifier{lu} */
    virtual void codegen_init_mem(CodeGenerator& g) const;

    /** \brief Codegen for free_mem

        \identifier{lv} */
    virtual void codegen_free_mem(CodeGenerator& g) const {}

    /** \brief Codegen for checkout

        \identifier{lw} */
    virtual void codegen_checkout(CodeGenerator& g) const;

    /** \brief Codegen for release

        \identifier{lx} */
    virtual void codegen_release(CodeGenerator& g) const;

    /** \brief Code generate the function

        \identifier{ly} */
    std::string signature(const std::string& fname) const;

    /** \brief Generate code for the declarations of the C function

        \identifier{lz} */
    virtual void codegen_declarations(CodeGenerator& g) const;

    /** \brief Generate code for the function body

        \identifier{m0} */
    virtual void codegen_body(CodeGenerator& g) const;

    /** \brief Thread-local memory object type

        \identifier{m1} */
    virtual std::string codegen_mem_type() const { return ""; }

    /** \brief Export / Generate C code for the dependency function

        \identifier{m2} */
    virtual std::string generate_dependencies(const std::string& fname, const Dict& opts) const;

    /** \brief Is codegen supported?

        \identifier{m3} */
    virtual bool has_codegen() const { return false;}

    /** \brief Jit dependencies

        \identifier{m4} */
    virtual void jit_dependencies(const std::string& fname) {}

    /** \brief Export function in a specific language

        \identifier{m5} */
    virtual void export_code(const std::string& lang,
      std::ostream &stream, const Dict& options) const;

    /** \brief Serialize type information

        \identifier{m6} */
    void serialize_type(SerializingStream &s) const override;

    /** \brief Serialize an object without type information

        \identifier{m7} */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Display object

        \identifier{m8} */
    void disp(std::ostream& stream, bool more) const override;

    /** \brief  Print more

        \identifier{m9} */
    virtual void disp_more(std::ostream& stream) const {}

    /** \brief Get function signature: name:(inputs)->(outputs)

        \identifier{ma} */
    std::string definition() const;

    /** \brief Print dimensions of inputs and outputs

        \identifier{mb} */
    void print_dimensions(std::ostream &stream) const;

    /** \brief Print free variables

        \identifier{mc} */
    virtual std::vector<std::string> get_free() const;

    /** \brief Get the unidirectional or bidirectional partition

        \identifier{md} */
    void get_partition(casadi_int iind, casadi_int oind, Sparsity& D1, Sparsity& D2,
                      bool compact, bool symmetric,
                      bool allow_forward, bool allow_reverse) const;

    ///@{
    /** \brief Number of input/output nonzeros

        \identifier{me} */
    casadi_int nnz_in() const;
    casadi_int nnz_in(casadi_int ind) const { return sparsity_in(ind).nnz(); }
    casadi_int nnz_out() const;
    casadi_int nnz_out(casadi_int ind) const { return sparsity_out(ind).nnz(); }
    ///@}

    ///@{
    /** \brief Number of input/output elements

        \identifier{mf} */
    casadi_int numel_in() const;
    casadi_int numel_in(casadi_int ind) const { return sparsity_in(ind).numel(); }
    casadi_int numel_out(casadi_int ind) const { return sparsity_out(ind).numel(); }
    casadi_int numel_out() const;
    ///@}

    ///@{
    /** \brief Input/output dimensions

        \identifier{mg} */
    casadi_int size1_in(casadi_int ind) const { return sparsity_in(ind).size1(); }
    casadi_int size2_in(casadi_int ind) const { return sparsity_in(ind).size2(); }
    casadi_int size1_out(casadi_int ind) const { return sparsity_out(ind).size1(); }
    casadi_int size2_out(casadi_int ind) const { return sparsity_out(ind).size2(); }
    std::pair<casadi_int, casadi_int> size_in(casadi_int ind) const {
      return sparsity_in(ind).size();
    }
    std::pair<casadi_int, casadi_int> size_out(casadi_int ind) const {
      return sparsity_out(ind).size();
    }
    ///@}

    ///@{
    /** \brief Input/output sparsity

        \identifier{mh} */
    const Sparsity& sparsity_in(casadi_int ind) const { return sparsity_in_.at(ind); }
    const Sparsity& sparsity_out(casadi_int ind) const { return sparsity_out_.at(ind); }
    ///@}

    ///@{
    /** \brief Are all inputs and outputs scalar

        \identifier{mi} */
    bool all_scalar() const;

    /// Is a Jacobian block known to be symmetric a priori?
    virtual bool jac_is_symm(casadi_int oind, casadi_int iind) const;

    /// Convert to compact Jacobian sparsity pattern
    Sparsity to_compact(casadi_int oind, casadi_int iind, const Sparsity& sp) const;

    /// Convert from compact Jacobian sparsity pattern
    Sparsity from_compact(casadi_int oind, casadi_int iind, const Sparsity& sp) const;

    /// Get the sparsity pattern via sparsity seed propagation
    template<bool fwd>
    Sparsity get_jac_sparsity_gen(casadi_int oind, casadi_int iind) const;

    /// A flavor of get_jac_sparsity_gen that does hierarchical block structure recognition
    Sparsity get_jac_sparsity_hierarchical(casadi_int oind, casadi_int iind) const;

    /** A flavor of get_jac_sparsity_gen that does hierarchical block
    * structure recognition for symmetric Jacobians
    */
    Sparsity get_jac_sparsity_hierarchical_symm(casadi_int oind, casadi_int iind) const;

    /// Get a vector of symbolic variables corresponding to the outputs
    virtual std::vector<MX> symbolic_output(const std::vector<MX>& arg) const;

    ///@{
    /** \brief Number of function inputs and outputs

        \identifier{mj} */
    virtual size_t get_n_in();
    virtual size_t get_n_out();
    ///@}


    ///@{
    /** \brief Names of function input and outputs

        \identifier{mk} */
    virtual std::string get_name_in(casadi_int i);
    virtual std::string get_name_out(casadi_int i);
    ///@}

    /** \brief Get default input value

        \identifier{ml} */
    virtual double get_default_in(casadi_int ind) const {
      return 0;
    }

    /** \brief Get largest input value

        \identifier{mm} */
    virtual double get_max_in(casadi_int ind) const {
      return inf;
    }

    /** \brief Get smallest input value

        \identifier{mn} */
    virtual double get_min_in(casadi_int ind) const {
      return -inf;
    }

    virtual std::vector<double> get_nominal_in(casadi_int ind) const {
      return std::vector<double>(nnz_in(ind), 1.);
    }

    virtual std::vector<double> get_nominal_out(casadi_int ind) const {
      return std::vector<double>(nnz_out(ind), 1.);
    }

    /** \brief Get relative tolerance

        \identifier{mo} */
    virtual double get_reltol() const {
      return eps;
    }

    /** \brief Get absolute tolerance

        \identifier{mp} */
    virtual double get_abstol() const {
      return eps;
    }

    /** \brief Get sparsity of a given input

        \identifier{mq} */
    virtual Sparsity get_sparsity_in(casadi_int i);

    /** \brief Get sparsity of a given output

        \identifier{mr} */
    virtual Sparsity get_sparsity_out(casadi_int i);

    /** \brief Which inputs are differentiable

        \identifier{ms} */
    virtual bool get_diff_in(casadi_int i) { return true; }

    /** \brief Which outputs are differentiable

        \identifier{mt} */
    virtual bool get_diff_out(casadi_int i) { return true; }

    /** \brief Get input scheme index by name

        \identifier{mu} */
    casadi_int index_in(const std::string &name) const {
      for (casadi_int i=0; i<name_in_.size(); ++i) {
        if (name_in_[i]==name) return i;
      }
      casadi_error("FunctionInternal::index_in: could not find entry \""
                   + name + "\". Available names are: " + str(name_in_) + ".");
      return -1;
    }

    /** \brief Get output scheme index by name

        \identifier{mv} */
    casadi_int index_out(const std::string &name) const {
      for (casadi_int i=0; i<name_out_.size(); ++i) {
        if (name_out_[i]==name) return i;
      }
      casadi_error("FunctionInternal::index_out: could not find entry \""
                   + name + "\". Available names are: " + str(name_out_) + ".");
      return -1;
    }

    /** \brief  Propagate sparsity forward

        \identifier{mw} */
    virtual int sp_forward(const bvec_t** arg, bvec_t** res,
                            casadi_int* iw, bvec_t* w, void* mem) const;

    /** \brief  Propagate sparsity forward, specific block

        \identifier{mx} */
    virtual int sp_forward_block(const bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem, casadi_int oind, casadi_int iind) const;

    /** \brief  Propagate sparsity backwards

        \identifier{my} */
    virtual int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const;

    /** \brief Get number of temporary variables needed

        \identifier{mz} */
    void sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const;

    /** \brief Get required length of arg field

        \identifier{n0} */
    size_t sz_arg() const { return sz_arg_per_ + sz_arg_tmp_;}

    /** \brief Get required length of res field

        \identifier{n1} */
    size_t sz_res() const { return sz_res_per_ + sz_res_tmp_;}

    /** \brief Get required length of iw field

        \identifier{n2} */
    size_t sz_iw() const { return sz_iw_per_ + sz_iw_tmp_;}

    /** \brief Get required length of w field

        \identifier{n3} */
    size_t sz_w() const { return sz_w_per_ + sz_w_tmp_;}

    /** \brief Ensure required length of arg field

        \identifier{n4} */
    void alloc_arg(size_t sz_arg, bool persistent=false);

    /** \brief Ensure required length of res field

        \identifier{n5} */
    void alloc_res(size_t sz_res, bool persistent=false);

    /** \brief Ensure required length of iw field

        \identifier{n6} */
    void alloc_iw(size_t sz_iw, bool persistent=false);

    /** \brief Ensure required length of w field

        \identifier{n7} */
    void alloc_w(size_t sz_w, bool persistent=false);

    /** \brief Ensure work vectors long enough to evaluate function

        \identifier{n8} */
    void alloc(const Function& f, bool persistent=false, int num_threads=1);

    /** \brief Set the (persistent) work vectors

        \identifier{n9} */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const {}

    /** \brief Set the (temporary) work vectors

        \identifier{na} */
    virtual void set_temp(void* mem, const double** arg, double** res,
                          casadi_int* iw, double* w) const {}

    /** \brief Set the (persistent and temporary) work vectors

        \identifier{nb} */
    void setup(void* mem, const double** arg, double** res, casadi_int* iw, double* w) const;

    ///@{
    /** \brief Calculate derivatives by multiplying the full Jacobian and multiplying

        \identifier{nc} */
    virtual bool fwdViaJac(casadi_int nfwd) const;
    virtual bool adjViaJac(casadi_int nadj) const;
    ///@}

    /** Obtain information about function */
    virtual Dict info() const;

    /** \brief Generate/retrieve cached serial map

        \identifier{nd} */
    Function map(casadi_int n, const std::string& parallelization) const;

    /** \brief Export an input file that can be passed to generate C code with a main

        \identifier{ne} */
    void generate_in(const std::string& fname, const double** arg) const;
    void generate_out(const std::string& fname, double** res) const;

    bool always_inline_, never_inline_;

    /// Number of inputs and outputs
    size_t n_in_, n_out_;

    /// Are inputs and outputs differentiable?
    std::vector<bool> is_diff_in_, is_diff_out_;

    /// Input and output sparsity
    std::vector<Sparsity> sparsity_in_, sparsity_out_;

    /// Input and output scheme
    std::vector<std::string> name_in_, name_out_;

    /** \brief  Use just-in-time compiler

        \identifier{nf} */
    bool jit_;

    /** \brief  Cleanup jit source file

        \identifier{ng} */
    bool jit_cleanup_;

    /** \brief  Serialize behaviour

        \identifier{nh} */
    std::string jit_serialize_;

    /** \brief  Name if jit source file

        \identifier{ni} */
    std::string jit_name_;

    std::string jit_base_name_;

    /** \brief Use a temporary name

        \identifier{nj} */
    bool jit_temp_suffix_;

    /** \brief Numerical evaluation redirected to a C function

        \identifier{nk} */
    eval_t eval_;

    /** \brief Checkout redirected to a C function

        \identifier{nl} */
    casadi_checkout_t checkout_;

   /** \brief Release redirected to a C function

       \identifier{nm} */
    casadi_release_t release_;

    /** \brief Dict of statistics (resulting from evaluate)

        \identifier{nn} */
    Dict stats_;

    /** \brief Reference counting in codegen?

        \identifier{no} */
    bool has_refcount_;

    /** \brief Values to prepopulate the function cache with

        \identifier{26h} */
    Dict cache_init_;

    /// Function cache
    mutable std::map<std::string, WeakRef> cache_;

    /// Cache for sparsities of the Jacobian blocks
    mutable std::vector<Sparsity> jac_sparsity_[2];

    /// If the function is the derivative of another function
    Function derivative_of_;

    /// User-set field
    void* user_data_;

    /// Just-in-time compiler
    std::string compiler_plugin_;
    Importer compiler_;
    Dict jit_options_;

    /// Penalty factor for using a complete Jacobian to calculate directional derivatives
    double jac_penalty_;

    // Types of derivative calculation permitted
    bool enable_forward_, enable_reverse_, enable_jacobian_, enable_fd_;
    bool enable_forward_op_, enable_reverse_op_, enable_jacobian_op_, enable_fd_op_;

    /// Weighting factor for derivative calculation and sparsity pattern calculation
    double ad_weight_, ad_weight_sp_;

    /// Maximum number of sensitivity directions
    casadi_int max_num_dir_;

    /// Errors are thrown if numerical values of inputs look bad
    bool inputs_check_;

    // Finite difference step
    Dict fd_options_;

    // Finite difference step size
    double fd_step_;

    // Finite difference method
    std::string fd_method_;

    // Print input/output
    bool print_in_;
    bool print_out_;

    // Warn when number of inputs or outputs exceed this value
    casadi_int max_io_;

    // Dump input/output
    bool dump_in_, dump_out_, dump_;

    // Directory to dump to
    std::string dump_dir_;

    // Format to dump with
    std::string dump_format_;

    // Forward/reverse/Jacobian options
    Dict forward_options_, reverse_options_, jacobian_options_, der_options_;

    // Store a reference to a custom Jacobian
    Function custom_jacobian_;

    // Counter for unique names for dumping inputs and output
    mutable casadi_int dump_count_ = 0;
#ifdef CASADI_WITH_THREAD
    mutable std::mutex dump_count_mtx_;
#endif // CASADI_WITH_THREAD

    /** \brief Check if the function is of a particular type

        \identifier{np} */
    virtual bool is_a(const std::string& type, bool recursive) const;

    /** \brief Can a derivative direction be skipped

        \identifier{nq} */
    template<typename MatType>
    static bool purgable(const std::vector<MatType>& seed);

    /** \brief Symbolic expressions for the forward seeds

        \identifier{nr} */
    template<typename MatType>
    std::vector<std::vector<MatType> >
    fwd_seed(casadi_int nfwd) const;

    /** \brief Symbolic expressions for the adjoint seeds

        \identifier{ns} */
    template<typename MatType>
    std::vector<std::vector<MatType> >
    symbolicAdjSeed(casadi_int nadj, const std::vector<MatType>& v) const;

    static std::string string_from_UnifiedReturnStatus(UnifiedReturnStatus status);

    /** \brief Deserializing constructor

        \identifier{nt} */
    explicit FunctionInternal(DeserializingStream& e);

    /** \brief Deserialize with type disambiguation

        \identifier{nu} */
    static Function deserialize(DeserializingStream& s);
    static std::map<std::string, ProtoFunction* (*)(DeserializingStream&)> deserialize_map;

    /** \brief Print inputs

        \identifier{nv} */
    void print_in(std::ostream &stream, const double** arg, bool truncate) const;

    /** \brief Print outputs

        \identifier{nw} */
    void print_out(std::ostream &stream, double** res, bool truncate) const;

  protected:
    /** \brief Populate jac_sparsity_ and jac_sparsity_compact_ during initialization

        \identifier{nx} */
    void set_jac_sparsity(casadi_int oind, casadi_int iind, const Sparsity& sp);

  private:
    // @{
    /// Dumping functionality
    casadi_int get_dump_id() const;
    void dump_in(casadi_int id, const double** arg) const;
    void dump_out(casadi_int id, double** res) const;
    void dump() const;
    // @}

    /** \brief Memory that is persistent during a call (but not between calls)

        \identifier{ny} */
    size_t sz_arg_per_, sz_res_per_, sz_iw_per_, sz_w_per_;

    /** \brief Temporary memory inside a function

        \identifier{nz} */
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
  fwd_seed(casadi_int nfwd) const {
    std::vector<std::vector<MatType>> fseed(nfwd);
    for (casadi_int dir=0; dir<nfwd; ++dir) {
      fseed[dir].resize(n_in_);
      for (casadi_int iind=0; iind<n_in_; ++iind) {
        std::string n = "f" + str(dir) + "_" +  name_in_[iind];
        Sparsity sp = is_diff_in_[iind] ? sparsity_in(iind) : Sparsity(size_in(iind));
        fseed[dir][iind] = MatType::sym(n, sp);
      }
    }
    return fseed;
  }

  template<typename MatType>
  std::vector<std::vector<MatType> >
  FunctionInternal::
  symbolicAdjSeed(casadi_int nadj, const std::vector<MatType>& v) const {
    std::vector<std::vector<MatType> > aseed(nadj, v);
    for (casadi_int dir=0; dir<nadj; ++dir) {
      // Replace symbolic inputs
      casadi_int oind=0;
      for (typename std::vector<MatType>::iterator i=aseed[dir].begin();
          i!=aseed[dir].end();
          ++i, ++oind) {
        // Name of the adjoint seed
        std::stringstream ss;
        ss << "a";
        if (nadj>1) ss << dir << "_";
        ss << oind;

        // Save to matrix
        *i = MatType::sym(ss.str(), is_diff_out_[oind] ? i->sparsity() : Sparsity(i->size()));

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
      std::pair<casadi_int, casadi_int> sz;
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
        res.resize(n_out_);
        M z = M::zeros(sz);
        for (auto&& a : res) a = z;
        // Call multiple times
        std::vector<M> arg1 = arg, res1;
        for (casadi_int c=0; c<sz.second; ++c) {
          for (casadi_int r=0; r<sz.first; ++r) {
            // Get scalar arguments
            for (casadi_int i=0; i<arg.size(); ++i) {
              if (arg[i].size()==sz) arg1[i] = arg[i](r, c);
            }
            // Call recursively with scalar arguments
            call(arg1, res1, always_inline, never_inline);
            // Get results
            casadi_assert_dev(res.size() == res1.size());
            for (casadi_int i=0; i<res.size(); ++i) res[i](r, c) = res1[i];
          }
        }
        // All elements assigned
        return;
      }
    }

    // Check if inputs need to be replaced
    casadi_int npar = 1;
    if (!matching_arg(arg, npar)) {
      return call(replace_arg(arg, npar), res, always_inline, never_inline);
    }

    // Call the type-specific method
    call_gen(arg, res, npar, always_inline, never_inline);
  }

  template<typename M>
  std::vector<M> FunctionInternal::
  project_arg(const std::vector<M>& arg, casadi_int npar) const {
    casadi_assert_dev(arg.size()==n_in_);

    // Which arguments require mapped evaluation
    std::vector<bool> mapped(n_in_);
    for (casadi_int i=0; i<n_in_; ++i) {
      mapped[i] = arg[i].size2()!=size2_in(i);
    }

    // Check if matching input sparsity
    std::vector<bool> matching(n_in_);
    bool any_mismatch = false;
    for (casadi_int i=0; i<n_in_; ++i) {
      if (mapped[i]) {
        matching[i] = arg[i].sparsity().is_stacked(sparsity_in(i), npar);
      } else {
        matching[i] = arg[i].sparsity()==sparsity_in(i);
      }
      any_mismatch = any_mismatch || !matching[i];
    }

    // Correct input sparsity
    if (any_mismatch) {
      std::vector<M> arg2(arg);
      for (casadi_int i=0; i<n_in_; ++i) {
        if (!matching[i]) {
          if (mapped[i]) {
            arg2[i] = project(arg2[i], repmat(sparsity_in(i), 1, npar));
          } else {
            arg2[i] = project(arg2[i], sparsity_in(i));
          }
        }
      }
      return arg2;
    }
    return arg;
  }

  template<typename M>
    std::vector<M> FunctionInternal::
    project_res(const std::vector<M>& arg, casadi_int npar) const {
      return arg;
  }

  template<typename D>
  void FunctionInternal::
  call_gen(const std::vector<Matrix<D> >& arg, std::vector<Matrix<D> >& res,
           casadi_int npar, bool always_inline, bool never_inline) const {
    casadi_assert(!never_inline, "Call-nodes only possible in MX expressions");
    std::vector< Matrix<D> > arg2 = project_arg(arg, npar);

    // Which arguments require mapped evaluation
    std::vector<bool> mapped(n_in_);
    for (casadi_int i=0; i<n_in_; ++i) {
      mapped[i] = arg[i].size2()!=size2_in(i);
    }

    // Allocate results
    res.resize(n_out_);
    for (casadi_int i=0; i<n_out_; ++i) {
      if (!res[i].sparsity().is_stacked(sparsity_out(i), npar)) {
        res[i] = Matrix<D>::zeros(repmat(sparsity_out(i), 1, npar));
      }
    }

    // Allocate temporary memory if needed
    std::vector<casadi_int> iw_tmp(sz_iw());
    std::vector<D> w_tmp(sz_w());

    // Get pointers to input arguments
    std::vector<const D*> argp(sz_arg());
    for (casadi_int i=0; i<n_in_; ++i) argp[i]=get_ptr(arg2[i]);

    // Get pointers to output arguments
    std::vector<D*> resp(sz_res());
    for (casadi_int i=0; i<n_out_; ++i) resp[i]=get_ptr(res[i]);

    // For all parallel calls
    for (casadi_int p=0; p<npar; ++p) {
      // Call memory-less
      if (eval_gen(get_ptr(argp), get_ptr(resp),
                   get_ptr(iw_tmp), get_ptr(w_tmp), memory(0))) {
        if (error_on_fail_) casadi_error("Evaluation failed");
      }
      // Update offsets
      if (p==npar-1) break;
      for (casadi_int i=0; i<n_in_; ++i) if (mapped[i]) argp[i] += nnz_in(i);
      for (casadi_int i=0; i<n_out_; ++i) resp[i] += nnz_out(i);
    }
  }

  template<typename M>
  void FunctionInternal::check_arg(const std::vector<M>& arg, casadi_int& npar) const {
    casadi_assert(arg.size()==n_in_, "Incorrect number of inputs: Expected "
                          + str(n_in_) + ", got " + str(arg.size()));
    for (casadi_int i=0; i<n_in_; ++i) {
      if (!check_mat(arg[i].sparsity(), sparsity_in(i), npar)) {
        // Dimensions
        std::string d_arg = str(arg[i].size1()) + "-by-" + str(arg[i].size2());
        std::string d_in = str(size1_in(i)) + "-by-" + str(size2_in(i));
        std::string e = "Input " + str(i) + " (" + name_in_[i] + ") has mismatching shape. "
                     "Got " + d_arg + ". Allowed dimensions, in general, are:\n"
                     " - The input dimension N-by-M (here " + d_in + ")\n"
                     " - A scalar, i.e. 1-by-1\n"
                     " - M-by-N if N=1 or M=1 (i.e. a transposed vector)\n"
                     " - N-by-M1 if K*M1=M for some K (argument repeated horizontally)\n";
        if (npar!=-1) {
          e += " - N-by-P*M, indicating evaluation with multiple arguments (P must be a "
                     "multiple of " + str(npar) + " for consistency with previous inputs)";
        }
        casadi_error(e);
      }
    }
  }

  template<typename M>
  void FunctionInternal::check_res(const std::vector<M>& res, casadi_int& npar) const {
    casadi_assert(res.size()==n_out_, "Incorrect number of outputs: Expected "
                          + str(n_out_) + ", got " + str(res.size()));
    for (casadi_int i=0; i<n_out_; ++i) {
      casadi_assert(check_mat(res[i].sparsity(), sparsity_out(i), npar),
                    "Output " + str(i) + " (" + name_out_[i] + ") has mismatching shape. "
                    "Expected " + str(size_out(i)) + ", got " + str(res[i].size()));
    }
  }

  template<typename M>
  bool FunctionInternal::matching_arg(const std::vector<M>& arg, casadi_int& npar) const {
    check_arg(arg, npar);
    for (casadi_int i=0; i<n_in_; ++i) {
      if (arg.at(i).size1()!=size1_in(i)) return false;
      if (arg.at(i).size2()!=size2_in(i) && arg.at(i).size2()!=npar*size2_in(i)) return false;
    }
    return true;
  }

  template<typename M>
  bool FunctionInternal::matching_res(const std::vector<M>& res, casadi_int& npar) const {
    check_res(res, npar);
    for (casadi_int i=0; i<n_out_; ++i) {
      if (res.at(i).size1()!=size1_out(i)) return false;
      if (res.at(i).size2()!=size2_out(i) && res.at(i).size2()!=npar*size2_out(i)) return false;
    }
    return true;
  }

  template<typename M>
  M replace_mat(const M& arg, const Sparsity& inp, casadi_int npar) {
    if (arg.size()==inp.size()) {
      // Matching dimensions already
      return arg;
    } else if (arg.is_empty()) {
      // Empty matrix means set zero
      return M(inp.size());
    } else if (arg.is_scalar()) {
      // Scalar assign means set all
      return M(inp, arg);
    } else if (arg.is_vector() && inp.size()==std::make_pair(arg.size2(), arg.size1())) {
      // Transpose vector
      return arg.T();
    } else if (arg.size1()==inp.size1() && arg.size2()>0 && inp.size2()>0
               && inp.size2()%arg.size2()==0) {
      // Horizontal repmat
      return repmat(arg, 1, inp.size2()/arg.size2());
    } else {
      casadi_assert_dev(npar!=-1);
      // Multiple evaluation
      return repmat(arg, 1, (npar*inp.size2())/arg.size2());
    }
  }

  template<typename M>
  std::vector<M> FunctionInternal::
  replace_arg(const std::vector<M>& arg, casadi_int npar) const {
    std::vector<M> r(arg.size());
    for (casadi_int i=0; i<r.size(); ++i) r[i] = replace_mat(arg[i], sparsity_in(i), npar);
    return r;
  }

  template<typename M>
  std::vector<M> FunctionInternal::
  replace_res(const std::vector<M>& res, casadi_int npar) const {
    std::vector<M> r(res.size());
    for (casadi_int i=0; i<r.size(); ++i) r[i] = replace_mat(res[i], sparsity_out(i), npar);
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> > FunctionInternal::
  replace_fseed(const std::vector<std::vector<M> >& fseed, casadi_int npar) const {
    std::vector<std::vector<M> > r(fseed.size());
    for (casadi_int d=0; d<r.size(); ++d) r[d] = replace_arg(fseed[d], npar);
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> > FunctionInternal::
  replace_aseed(const std::vector<std::vector<M> >& aseed, casadi_int npar) const {
    std::vector<std::vector<M> > r(aseed.size());
    for (casadi_int d=0; d<r.size(); ++d) r[d] = replace_res(aseed[d], npar);
    return r;
  }

  template<typename M>
  std::map<std::string, M> FunctionInternal::
  convert_arg(const std::vector<M>& arg) const {
    casadi_assert(arg.size()==n_in_, "Incorrect number of inputs: Expected "
                      + str(n_in_) + ", got " + str(arg.size()));
    std::map<std::string, M> ret;
    for (casadi_int i=0;i<n_in_;++i) {
      ret[name_in_[i]] = arg[i];
    }
    return ret;
  }

  template<typename M>
  std::vector<M> FunctionInternal::
  convert_arg(const std::map<std::string, M>& arg) const {
    // Get default inputs
    std::vector<M> arg_v(n_in_);
    for (casadi_int i=0; i<arg_v.size(); ++i) {
      arg_v[i] = get_default_in(i);
    }

    // Assign provided inputs
    for (auto&& e : arg) {
      arg_v.at(index_in(e.first)) = e.second;
    }

    return arg_v;
  }

  template<typename M>
  std::map<std::string, M> FunctionInternal::
  convert_res(const std::vector<M>& res) const {
    casadi_assert(res.size()==n_out_, "Incorrect number of outputs: Expected "
                      + str(n_out_) + ", got " + str(res.size()));
    std::map<std::string, M> ret;
    for (casadi_int i=0;i<n_out_;++i) {
      ret[name_out_[i]] = res[i];
    }
    return ret;
  }

  template<typename M>
  std::vector<M> FunctionInternal::
  convert_res(const std::map<std::string, M>& res) const {
    // Get default inputs
    std::vector<M> res_v(n_out_);
    for (casadi_int i=0; i<res_v.size(); ++i) {
      res_v[i] = std::numeric_limits<double>::quiet_NaN();
    }

    // Assign provided inputs
    for (auto&& e : res) {
      M a = e.second;
      res_v.at(index_out(e.first)) = a;
    }
    return res_v;
  }

  template<typename MatType>
  MatType FunctionInternal::ensure_stacked(const MatType& v, const Sparsity& sp, casadi_int n) {
    // Check dimensions
    if (v.size1() == sp.size1() && v.size2() == n * sp.size2()) {
      // Ensure that sparsity is a horizontal multiple of original input, or has no entries
      if (v.nnz() != 0 && !v.sparsity().is_stacked(sp, n)) {
        return project(v, repmat(sp, 1, n));
      }
    } else {
      // Correct empty sparsity
      casadi_assert_dev(v.is_empty());
      return MatType(sp.size1(), sp.size2() * n);
    }
    // No correction needed
    return v;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_FUNCTION_INTERNAL_HPP
