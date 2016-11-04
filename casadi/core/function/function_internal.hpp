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
#include <stack>
#include "code_generator.hpp"
#include "importer.hpp"
#include "../sparse_storage.hpp"
#include "../options.hpp"

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
  class CASADI_EXPORT FunctionInternal : public SharedObjectNode {
  public:
    /** \brief Constructor */
    FunctionInternal(const std::string& name);

    /** \brief  Destructor */
    virtual ~FunctionInternal() = 0;

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

#ifdef WITH_DEPRECATED_FEATURES
    /** \brief [DEPRECATED] Which variables enter nonlinearly
    *
    * Use which_depends instead.
    */
    virtual std::vector<bool> nl_var(const std::string& s_in,
                             const std::vector<std::string>& s_out) const;
#endif

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
    void _eval(const double** arg, double** res, int* iw, double* w, int mem);
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;
    ///@}

    /** \brief  Evaluate numerically, simplied syntax */
    virtual void simple(const double* arg, double* res);

    /** \brief  Evaluate with symbolic scalars */
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    ///@{
    /** \brief  Evaluate with symbolic matrices */
    virtual void eval_mx(const MXVector& arg, MXVector& res, bool always_inline, bool never_inline);
    ///@}

    ///@{
    /** \brief Evaluate a function, overloaded */
    void _eval(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);
    void _eval(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);
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
    virtual Function get_forward_old(const std::string& name, int nfwd, Dict& opts);
    virtual Function get_forward(const std::string& name, int nfwd,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_forward() const { return 0;}
    ///@}

    ///@{
    /** \brief Return function that calculates adjoint derivatives
     *    reverse(nadj) returns a cached instance if available,
     *    and calls <tt>Function get_reverse(int nadj)</tt>
     *    if no cached version is available.
     */
    Function reverse(int nadj);
    virtual Function get_reverse_old(const std::string& name, int nadj, Dict& opts);
    virtual Function get_reverse(const std::string& name, int nadj,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names, const Dict& opts);
    virtual int get_n_reverse() const { return 0;}
    ///@}

    /** \brief returns a new function with a selection of inputs/outputs of the original */
    virtual Function slice(const std::string& name, const std::vector<int>& order_in,
                           const std::vector<int>& order_out, const Dict& opts) const;

    /** \brief Get oracle */
    virtual const Function& oracle() const;

    /** \brief Propagate options */
    virtual Dict derived_options() const;

    /** \brief Can derivatives be calculated in any way? */
    bool hasDerivative() const;

    /** \brief  Weighting factor for chosing forward/reverse mode */
    virtual double ad_weight() const;

    /** \brief  Weighting factor for chosing forward/reverse mode,
        sparsity propagation */
    virtual double sp_weight() const;

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
    virtual std::vector<SX> free_sx() const;

    /** \brief Does the function have free variables */
    virtual bool has_free() const { return false;}

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
    virtual int n_nodes() const;

    /** \brief Wrap in an Function instance consisting of only one MX call */
    Function wrap() const;

    /** \brief Generate code the function */
    virtual void generateFunction(CodeGenerator& g, const std::string& fname,
                                  bool decl_static) const;

    /** \brief Generate meta-information allowing a user to evaluate a generated function */
    void generateMeta(CodeGenerator& g, const std::string& fname) const;

    /** \brief Use simplified signature */
    virtual bool simplifiedCall() const { return false;}

    /** \brief Generate shorthand macro */
    void addShorthand(CodeGenerator& g, const std::string& name) const;

    /** \brief Get name of the evaluation function */
    std::string eval_name() const;

    /** \brief Get name in codegen */
    virtual std::string codegen_name(const CodeGenerator& g) const;

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const;

    /** \brief Codegen incref for dependencies */
    virtual void codegen_incref(CodeGenerator& g) const {}

    /** \brief Codegen decref for dependencies */
    virtual void codegen_decref(CodeGenerator& g) const {}

    /** \brief Code generate the function  */
    std::string signature(const std::string& fname) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Generate code for the function body */
    virtual void generateBody(CodeGenerator& g) const;

    /** \brief Export / Generate C code for the dependency function */
    virtual std::string generate_dependencies(const std::string& fname, const Dict& opts);

    /** \brief Is codegen supported? */
    virtual bool has_codegen() const { return false;}

    /** \brief Jit dependencies */
    virtual void jit_dependencies(const std::string& fname) {}

    /** \brief  Print */
    virtual void print(std::ostream &stream) const;

    /** \brief  Print */
    virtual void repr(std::ostream &stream) const;

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
    void getPartition(int iind, int oind, Sparsity& D1, Sparsity& D2, bool compact, bool symmetric);

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

    /** \brief Name of the function */
    const std::string& name() const { return name_;}

    /// Generate the sparsity of a Jacobian block
    virtual Sparsity getJacSparsity(int iind, int oind, bool symmetric);

    /// Get the sparsity pattern, forward mode
    template<bool fwd>
    Sparsity getJacSparsityGen(int iind, int oind, bool symmetric, int gr_i=1, int gr_o=1);

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
    virtual void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

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

    ///@{
    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const {}
    virtual void _set_work(const double**& arg, double**& res,
                           int*& iw, double*& w, int mem) const {
      set_work(memory(mem), arg, res, iw, w);
    }
    ///@}

    ///@{
    /** \brief Set the (temporary) work vectors */
    virtual void set_temp(void* mem, const double** arg, double** res,
                          int* iw, double* w) const {}
    virtual void _set_temp(const double** arg, double** res,
                          int* iw, double* w, int mem) const {
      set_temp(memory(mem), arg, res, iw, w);
    }
    ///@}

    ///@{
    /** \brief Set the (persistent and temporary) work vectors */
    void setup(void* mem, const double** arg, double** res, int* iw, double* w) const;
    void _setup(const double** arg, double** res, int* iw, double* w, int mem) const {
      setup(memory(mem), arg, res, iw, w);
    }
    ///@}

    ///@{
    /** \brief Calculate derivatives by multiplying the full Jacobian and multiplying */
    virtual bool fwdViaJac(int nfwd);
    virtual bool adjViaJac(int nadj);
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
    std::vector<WeakRef> forward_, reverse_;

    /// Cache for full Jacobian
    WeakRef full_jacobian_;

    /// Cache for sparsities of the Jacobian blocks
    SparseStorage<Sparsity> jac_sparsity_, jac_sparsity_compact_;

    /// Cache for Jacobians
    SparseStorage<WeakRef> jac_, jac_compact_;

    /// If the function is the derivative of another function
    Function derivative_of_;

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

    /** \brief Get type name */
    virtual std::string type_name() const;

    /** \brief Check if the function is of a particular type */
    virtual bool is_a(const std::string& type, bool recursive) const;

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
  FunctionInternal::symbolicFwdSeed(int nfwd, const std::vector<MatType>& v) {
    std::vector<std::vector<MatType> > fseed(nfwd, v);
    for (int dir=0; dir<nfwd; ++dir) {
      // Replace symbolic inputs
      int iind=0;
      for (auto i=fseed[dir].begin(); i!=fseed[dir].end(); ++i, ++iind) {
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
    _eval(get_ptr(argp), get_ptr(resp), get_ptr(iw_tmp), get_ptr(w_tmp), 0);
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
        casadi_assert_message(checkMat(fseed[d][i].sparsity(), sparsity_in(i)),
                              "Forward seed " << i << "(" << name_in(i) << ") for direction " << d
                              << " has mismatching shape. Expected " << size_in(i)
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
        casadi_assert_message(checkMat(aseed[d][i].sparsity(), sparsity_out(i)),
                              "Adjoint seed " << i << " (" << name_out(i) << ") for direction " << d
                              << " has mismatching shape. Expected " << size_out(i)
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
                            "Input " << i << " (" << name_in(i) << ") has mismatching shape. "
                            << "Expected " << size_in(i) << ", got " << arg[i].size());
    }
  }

  template<typename M>
  void FunctionInternal::checkRes(const std::vector<M>& res) const {
    int n_out = this->n_out();
    casadi_assert_message(res.size()==n_out, "Incorrect number of outputs: Expected "
                          << n_out << ", got " << res.size());
    for (int i=0; i<n_out; ++i) {
      casadi_assert_message(checkMat(res[i].sparsity(), sparsity_out(i)),
                            "Output " << i << " (" << name_out(i) << ") has mismatching shape. "
                            "Expected " << size_out(i) << ", got " << res[i].size());
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
