/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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


#ifndef CASADI_CALL_HPP
#define CASADI_CALL_HPP

#include "multiple_output.hpp"
#include "function.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Embeds a function call in an expression graph
      \author Joel Andersson
      \date 2010-2015
  */
  class CASADI_EXPORT Call : public MultipleOutput {
  public:
    /** \brief  Create function call node

        \identifier{6j} */
    static std::vector<MX> create(const Function& fcn, const std::vector<MX>& arg);

    /** \brief  Create function call node

        \identifier{287} */
    static MX create_call(const Function& fcn, const std::vector<MX>& arg);

    /** \brief  Get an output

        \identifier{2bt} */
    MX get_output(casadi_int oind) const override;

    /** \brief  Destructor

        \identifier{6k} */
    ~Call() override {}

    /** \brief Project a function input to a particular sparsity

        \identifier{6l} */
    static MX projectArg(const MX& x, const Sparsity& sp, casadi_int i);

    /** \brief  Print expression

        \identifier{6m} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Add a dependent function

        \identifier{6n} */
    void add_dependency(CodeGenerator& g) const override;

    /** \brief Is reference counting needed in codegen?

        \identifier{6o} */
    bool has_refcount() const override;

    /** \brief Codegen incref

        \identifier{6p} */
    void codegen_incref(CodeGenerator& g, std::set<void*>& added) const override;

    /** \brief Codegen decref

        \identifier{6q} */
    void codegen_decref(CodeGenerator& g, std::set<void*>& added) const override;

    /** \brief Generate code for the operation

        \identifier{6r} */
    void generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res,
                          const std::vector<bool>& arg_is_ref,
                          std::vector<bool>& res_is_ref) const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{6s} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{6t} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{6u} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward

        \identifier{6v} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{6w} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Get called function

        \identifier{6x} */
    const Function& which_function() const override { return fcn_;}

    /** \brief  Get function output

        \identifier{6y} */
    casadi_int which_output() const override { return -1;}

    /** \brief  Number of outputs

        \identifier{6z} */
    casadi_int nout() const override;

    /** \brief  Get the sparsity of output oind

        \identifier{70} */
    const Sparsity& sparsity(casadi_int oind) const override;

    /** \brief Get the operation

        \identifier{71} */
    casadi_int op() const override { return OP_CALL;}

    /** \brief Get required length of arg field

        \identifier{72} */
    size_t sz_arg() const override;

    /** \brief Get required length of res field

        \identifier{73} */
    size_t sz_res() const override;

    /** \brief Get required length of iw field

        \identifier{74} */
    size_t sz_iw() const override;

    /** \brief Get required length of w field

        \identifier{75} */
    size_t sz_w() const override;

    /** \brief Serialize an object without type information

        \identifier{76} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information

        \identifier{77} */
    static MXNode* deserialize(DeserializingStream& s) { return new Call(s); }

  protected:
    /** \brief  Constructor (should not be used directly)

        \identifier{78} */
    explicit Call(const Function& fcn, const std::vector<MX>& arg);

    /** \brief Deserializing constructor

        \identifier{79} */
    explicit Call(DeserializingStream& s);

    /** \brief Find a common conditional argument for all seeds

        \identifier{7a} */
    static MX common_cond(const std::vector<std::vector<MX>>& seed);

    // Function to be evaluated
    Function fcn_;

    /// Output node cache
    mutable WeakCache<casadi_int, MX> cache_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_CALL_HPP
