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

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "x_function.hpp"
#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {

#ifndef SWIG
  /** \brief  An element of the algorithm, namely an MX node */
  struct MXAlgEl {
    /// Operator index
    int op;

    /// Data associated with the operation
    MX data;

    /// Work vector indices of the arguments
    std::vector<int> arg;

    /// Work vector indices of the results
    std::vector<int> res;
  };
#endif // SWIG

  /** \brief  Internal node class for MXFunction
      \author Joel Andersson
      \date 2010-2015
  */
  class CASADI_EXPORT MXFunction :
        public XFunction<MXFunction, MX, MXNode>{
  public:
    /** \brief  An element of the algorithm, namely an MX node */
    typedef MXAlgEl AlgEl;

    /** \brief  All the runtime elements in the order of evaluation */
    std::vector<AlgEl> algorithm_;

    /** \brief Offsets for elements in the w_ vector */
    std::vector<int> workloc_;

    /// Free variables
    std::vector<MX> free_vars_;

    /// Default input values
    std::vector<double> default_in_;

    /** \brief Constructor */
    MXFunction(const std::string& name, const std::vector<MX>& input,
                       const std::vector<MX>& output);

    /** \brief  Destructor */
    ~MXFunction() override;

    /** \brief  Evaluate numerically, work vectors given */
    void eval(void* mem, const double** arg, double** res, int* iw, double* w) const override;

    /** \brief  Print description */
    void print(std::ostream &stream) const override;

    /** \brief Get type name */
    std::string type_name() const override {return "mxfunction";}

    /** \brief Check if the function is of a particular type */
    bool is_a(const std::string& type, bool recursive) const override;

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Generate code for the declarations of the C function */
    void generateDeclarations(CodeGenerator& g) const override;

    /** \brief Codegen incref for dependencies */
    void codegen_incref(CodeGenerator& g) const override;

    /** \brief Codegen decref for dependencies */
    void codegen_decref(CodeGenerator& g) const override;

    /** \brief Generate code for the body of the C function */
    void generateBody(CodeGenerator& g) const override;

    /** \brief Extract the residual function G and the modified function Z out of an expression
     * (see Albersmeyer2010 paper) */
    void generate_lifted(Function& vdef_fcn, Function& vinit_fcn) const override;

    /** Inline calls? */
    bool should_inline(bool always_inline, bool never_inline) const override {
      return always_inline;
    }

    /** \brief Evaluate symbolically, SX type*/
    void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const override;

    /** \brief Evaluate symbolically, MX type */
    void eval_mx(const MXVector& arg, MXVector& res,
                         bool always_inline, bool never_inline) const override;

    /** \brief Calculate forward mode directional derivatives */
    void eval_forward(const std::vector<std::vector<MX> >& fwdSeed,
                        std::vector<std::vector<MX> >& fwdSens) const;

    /** \brief Calculate reverse mode directional derivatives */
    void eval_reverse(const std::vector<std::vector<MX> >& adjSeed,
                        std::vector<std::vector<MX> >& adjSens) const;

    /** \brief Expand the matrix valued graph into a scalar valued graph */
    Function expand(const std::vector<SX>& inputv);

    /// Get a vector of symbolic variables corresponding to the outputs
    std::vector<MX> symbolicOutput(const std::vector<MX>& arg) override;

    /** \brief  Propagate sparsity forward */
    void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

    /** \brief  Propagate sparsity backwards */
    void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

    // print an element of an algorithm
    std::string print(const AlgEl& el) const;

    /** \brief Gradient expression */
    MX grad_mx(int iind=0, int oind=0) override;

    /** \brief Tangent expression */
    MX tang_mx(int iind=0, int oind=0) override;

    /** \brief Jacobian expression */
    MX jac_mx(int iind=0, int oind=0, const Dict& opts = Dict()) override;

    ///@{
    /** \brief Get function input(s) and output(s)  */
    const MX mx_in(int ind) const override;
    const std::vector<MX> mx_in() const override;
    ///@}

    /// Get free variables (MX)
    std::vector<MX> free_mx() const override {return free_vars_;}

    /** \brief Does the function have free variables */
    bool has_free() const override { return !free_vars_.empty();}

    /** \brief Print free variables */
    void print_free(std::ostream &stream) const override {
      stream << free_vars_;
    }

    /** \brief Number of nodes in the algorithm */
    int n_nodes() const override { return algorithm_.size();}

    /** \brief Get default input value */
    double default_in(int ind) const override { return default_in_.at(ind);}

    /// Substitute inplace, internal implementation
    void substitute_inplace(std::vector<MX>& vdef, std::vector<MX>& ex) const;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_MX_FUNCTION_HPP
