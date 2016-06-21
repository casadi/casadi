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
#include "../function/function.hpp"

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
    virtual ~GenericCall() {}

    /** \brief  Number of functions */
    virtual int numFunctions() const = 0;

    /** \brief  Get function reference */
    virtual const Function& getFunction(int i) const = 0;

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
    virtual ~Call() {}

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const;

    /** \brief Is reference counting needed in codegen? */
    virtual bool has_refcount() const;

    /** \brief Codegen incref */
    virtual void codegen_incref(CodeGenerator& g, std::set<void*>& added) const;

    /** \brief Codegen decref */
    virtual void codegen_decref(CodeGenerator& g, std::set<void*>& added) const;

    /** \brief Generate code for the operation */
    virtual void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const;

    /// Evaluate the function numerically
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem) const;

    /// Evaluate the function symbolically (SX)
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Propagate sparsity forward */
    virtual void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Number of functions */
    virtual int numFunctions() const {return 1;}

    /** \brief  Get function reference */
    virtual const Function& getFunction(int i) const { return fcn_;}

    /** \brief  Get function input */
    virtual int getFunction_input() const { return -1;}

    /** \brief  Get function output */
    virtual int getFunctionOutput() const { return -1;}

    /** \brief  Number of outputs */
    virtual int nout() const;

    /** \brief  Get the sparsity of output oind */
    virtual const Sparsity& sparsity(int oind) const;

    /** \brief Get the operation */
    virtual int op() const { return OP_CALL;}

    /** \brief Get required length of arg field */
    virtual size_t sz_arg() const;

    /** \brief Get required length of res field */
    virtual size_t sz_res() const;

    /** \brief Get required length of iw field */
    virtual size_t sz_iw() const;

    /** \brief Get required length of w field */
    virtual size_t sz_w() const;

  protected:
    /** \brief  Constructor (should not be used directly) */
    explicit Call(const Function& fcn, const std::vector<MX>& arg);

    // Function to be evaluated
    Function fcn_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_CALL_HPP
