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


#ifndef CASADI_CALL_HPP
#define CASADI_CALL_HPP

#include "multiple_output.hpp"
#include "function.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Base class for nodes involving function calls
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT GenericCall : public MultipleOutput {
  public:

    /** \brief Constructor */
    GenericCall() {}

    /** \brief Destructor */
    ~GenericCall() override {}

    /** \brief  Number of functions */
    int numFunctions() const override = 0;

    /** \brief  Get function reference */
    const Function& getFunction(int i) const override = 0;

    /** \brief Project a function input to a particular sparsity */
    static MX projectArg(const MX& x, const Sparsity& sp, int i);
  };

  /** Embeds a function call in an expression graph
      \author Joel Andersson
      \date 2010-2015
  */
  class CASADI_EXPORT Call : public GenericCall {
  public:
    /** \brief  Create function call node */
    static std::vector<MX> create(const Function& fcn, const std::vector<MX>& arg);

    /** \brief  Destructor */
    ~Call() override {}

    /** \brief  Print expression */
    std::string print(const std::vector<std::string>& arg) const override;

    /** \brief Add a dependent function */
    void addDependency(CodeGenerator& g) const override;

    /** \brief Is reference counting needed in codegen? */
    bool has_refcount() const override;

    /** \brief Codegen incref */
    void codegen_incref(CodeGenerator& g, std::set<void*>& added) const override;

    /** \brief Codegen decref */
    void codegen_decref(CodeGenerator& g, std::set<void*>& added) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const override;

    /// Evaluate the function numerically
    void eval(const double** arg, double** res, int* iw, double* w, int mem) const override;

    /// Evaluate the function symbolically (SX)
    void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void eval_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void eval_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward */
    void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

    /** \brief  Propagate sparsity backwards */
    void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

    /** \brief  Number of functions */
    int numFunctions() const override {return 1;}

    /** \brief  Get function reference */
    const Function& getFunction(int i) const override { return fcn_;}

    /** \brief  Get function input */
    int getFunction_input() const override { return -1;}

    /** \brief  Get function output */
    int getFunctionOutput() const override { return -1;}

    /** \brief  Number of outputs */
    int nout() const override;

    /** \brief  Get the sparsity of output oind */
    const Sparsity& sparsity(int oind) const override;

    /** \brief Get the operation */
    int op() const override { return OP_CALL;}

    /** \brief Get required length of arg field */
    size_t sz_arg() const override;

    /** \brief Get required length of res field */
    size_t sz_res() const override;

    /** \brief Get required length of iw field */
    size_t sz_iw() const override;

    /** \brief Get required length of w field */
    size_t sz_w() const override;

  protected:
    /** \brief  Constructor (should not be used directly) */
    explicit Call(const Function& fcn, const std::vector<MX>& arg);

    // Function to be evaluated
    Function fcn_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_CALL_HPP
