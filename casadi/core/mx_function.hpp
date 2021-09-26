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


#ifndef CASADI_MX_FUNCTION_HPP
#define CASADI_MX_FUNCTION_HPP

#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "x_function.hpp"
#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {

#ifndef SWIG
  /** \brief  An element of the algorithm, namely an MX node

      \identifier{1x} */
  struct MXAlgEl {
    /// Operator index
    casadi_int op;

    /// Data associated with the operation
    MX data;

    /// Work vector indices of the arguments
    std::vector<casadi_int> arg;

    /// Work vector indices of the results
    std::vector<casadi_int> res;

    /// Data type
    std::string data_type;
  };
#endif // SWIG

  /** \brief  Internal node class for MXFunction

      \author Joel Andersson
      \date 2010-2015

      \identifier{1y} */
  class CASADI_EXPORT MXFunction :
        public XFunction<MXFunction, MX, MXNode>{
  public:
    /** \brief  An element of the algorithm, namely an MX node

        \identifier{1z} */
    typedef MXAlgEl AlgEl;

    /** \brief  All the runtime elements in the order of evaluation

        \identifier{20} */
    std::vector<AlgEl> algorithm_;

    /** \brief Offsets for elements in the w_ vector

        \identifier{21} */
    std::vector<casadi_int> workloc_;
    std::vector<casadi_int> workloc_sz_self_;
    std::vector<std::string> worktype_;

    /** \brief Offset for the 'extra' working memory (non-io) */
    size_t w_extra_offset_;

    /// Free variables
    std::vector<MX> free_vars_;

    /// Default input values
    std::vector<double> default_in_;

    /// Offset for copy elision
    std::vector<casadi_int> ce_off_;

    /// Number elided locations
    casadi_int n_ce_;

    /// Copy elision used?
    bool ce_active_;

    /// 0-buffer for copy elision
    std::vector<double> zero_array_;
    /// 0-buffer for copy elision
    std::vector<bvec_t> zero_array_bvec_t_;

    /// Live variables?
    bool live_variables_;

    /// Print instructions during evaluation
    bool print_instructions_;

    /** \brief Constructor

        \identifier{22} */
    MXFunction(const std::string& name,
      const std::vector<MX>& input, const std::vector<MX>& output,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out);

    /** \brief  Destructor

        \identifier{23} */
    ~MXFunction() override;

    /** \brief  Evaluate numerically, work vectors given

        \identifier{24} */
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /** \brief  Print description

        \identifier{25} */
    void disp_more(std::ostream& stream) const override;

    /** \brief Get type name

        \identifier{26} */
    std::string class_name() const override {return "MXFunction";}

    /** \brief Check if the function is of a particular type

        \identifier{27} */
    bool is_a(const std::string& type, bool recursive) const override;

    ///@{
    /** \brief Options

        \identifier{28} */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /// Reconstruct options dict
    Dict generate_options(bool is_temp, bool keep_dim) const override;

    /** \brief  Initialize

        \identifier{29} */
    void init(const Dict& opts) override;

    /** \brief Generate code for the declarations of the C function

        \identifier{2a} */
    void codegen_declarations(CodeGenerator& g, const Instance& inst) const override;

    /** \brief Codegen incref for dependencies

        \identifier{2b} */
    void codegen_incref(CodeGenerator& g, const Instance& inst) const override;

    /** \brief Codegen decref for dependencies

        \identifier{2c} */
    void codegen_decref(CodeGenerator& g, const Instance& inst) const override;

    /** \brief Generate code for the body of the C function

        \identifier{2d} */
    void codegen_body(CodeGenerator& g, const Instance& inst) const override;

    /** \brief Serialize an object without type information

        \identifier{2e} */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation

        \identifier{2f} */
    static ProtoFunction* deserialize(DeserializingStream& s);

    /** \brief Extract the residual function G and the modified function Z out of an expression

     * (see Albersmeyer2010 paper)

        \identifier{2g} */
    void generate_lifted(Function& vdef_fcn, Function& vinit_fcn) const override;

    /** Inline calls? */
    bool should_inline(bool always_inline, bool never_inline) const override;

    /** \brief Evaluate symbolically, SX type

        \identifier{2h} */
    int eval_sx(const SXElem** arg, SXElem** res,
                casadi_int* iw, SXElem* w, void* mem) const override;

    /** \brief Evaluate symbolically, MX type

        \identifier{2i} */
    void eval_mx(const MXVector& arg, MXVector& res,
                 bool always_inline, bool never_inline) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{2j} */
    void ad_forward(const std::vector<std::vector<MX> >& fwdSeed,
                        std::vector<std::vector<MX> >& fwdSens) const;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{2k} */
    void ad_reverse(const std::vector<std::vector<MX> >& adjSeed,
                        std::vector<std::vector<MX> >& adjSens) const;

    /// Get a vector of symbolic variables corresponding to the outputs
    std::vector<MX> symbolic_output(const std::vector<MX>& arg) const override;

    /** \brief  Propagate sparsity forward

        \identifier{2l} */
    int sp_forward(const bvec_t** arg, bvec_t** res,
                  casadi_int* iw, bvec_t* w, void* mem) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{2m} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const override;

    // print an element of an algorithm
    std::string print(const AlgEl& el) const;

    // Print the input arguments of an instruction
    void print_arg(std::ostream &stream, casadi_int k, const AlgEl& el, const double** arg) const;

    // Print the output arguments of an instruction
    void print_res(std::ostream &stream, casadi_int k, const AlgEl& el, double** res) const;

    ///@{
    /** \brief Get function input(s) and output(s)

        \identifier{2n} */
    const MX mx_in(casadi_int ind) const override;
    const std::vector<MX> mx_in() const override;
    ///@}

    /// Get free variables (MX)
    std::vector<MX> free_mx() const override {return free_vars_;}

    /** \brief Does the function have free variables

        \identifier{2o} */
    bool has_free() const override { return !free_vars_.empty();}

    /** \brief Print free variables

        \identifier{2p} */
    std::vector<std::string> get_free() const override {
      std::vector<std::string> ret;
      for (auto&& e : free_vars_) ret.push_back(e.name());
      return ret;
    }

    // Get list of dependency functions
    std::vector<std::string> get_function() const override;

    // Get a dependency function
    const Function& get_function(const std::string &name) const override;

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Layout get_layout_in(casadi_int i) override { return in_.at(i).layout();}
    Layout get_layout_out(casadi_int i) override { return out_.at(i).layout();}
    /// @}

    /** \brief Number of nodes in the algorithm

        \identifier{2q} */
    casadi_int n_nodes() const override { return algorithm_.size();}

    casadi_int n_instructions() const override { return algorithm_.size();}

    /** *\brief get MX expression associated with instruction

         \identifier{2r} */
    MX instruction_MX(casadi_int k) const override;

    /** \brief Get an atomic operation operator index

        \identifier{2s} */
    casadi_int instruction_id(casadi_int k) const override { return algorithm_.at(k).op;}

    /** \brief Get default input value

        \identifier{2t} */
    double get_default_in(casadi_int ind) const override { return default_in_.at(ind);}

    /** \brief Get the (integer) input arguments of an atomic operation

        \identifier{2u} */
    std::vector<casadi_int> instruction_input(casadi_int k) const override;

    /** \brief Get the (integer) output argument of an atomic operation

        \identifier{2v} */
    std::vector<casadi_int> instruction_output(casadi_int k) const override;

    /** \brief Export function in a specific language

        \identifier{2w} */
    void export_code_body(const std::string& lang,
      std::ostream &stream, const Dict& options) const override;

    /// Substitute inplace, internal implementation
    void substitute_inplace(std::vector<MX>& vdef, std::vector<MX>& ex) const;

    // Get all embedded functions, recursively
    void find(std::map<FunctionInternal*, Function>& all_fun, casadi_int max_depth) const override;

    /** \brief Change option after object creation for debugging

        \identifier{2x} */
    void change_option(const std::string& option_name, const GenericType& option_value) override;

  protected:
    /** \brief Deserializing constructor

        \identifier{2y} */
    explicit MXFunction(DeserializingStream& s);
  };

} // namespace casadi
/// \endcond

#endif // CASADI_MX_FUNCTION_HPP
