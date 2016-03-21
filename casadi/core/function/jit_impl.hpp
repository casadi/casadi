
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


#ifndef CASADI_JIT_IMPL_HPP
#define CASADI_JIT_IMPL_HPP

#include "jit.hpp"
#include "function_internal.hpp"

namespace casadi {

  class CASADI_EXPORT Jit : public FunctionInternal {
  public:
    /** \brief Constructor */
    Jit(const std::string& name, int n_in, int n_out,
        const std::string& body, const Dict& opts);

    /** \brief Get type name */
    virtual std::string type_name() const { return "jit";}

    /** \brief Destructor */
    virtual ~Jit();

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief Initialize */
    virtual void init(const Dict& opts);

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return n_in_;}
    virtual size_t get_n_out() { return n_out_;}
    ///@}

    /// @{
    /** \brief All inputs and outputs are scalars */
    virtual Sparsity get_sparsity_in(int i) { return Sparsity::scalar();}
    virtual Sparsity get_sparsity_out(int i) { return Sparsity::scalar();}
    /// @}

    /** \brief Use simplified signature */
    virtual bool simplifiedCall() const { return true;}

    /** \brief Generate code for the function body */
    virtual void generateBody(CodeGenerator& g) const;

    ///@{
    /** \brief Jacobian of all outputs with respect to all inputs */
    bool hasFullJacobian() const;
    virtual Function getFullJacobian(const std::string& name, const Dict& opts);
    ///@}

  private:
    // Number of inputs and outputs
    int n_in_, n_out_;

    // Function body
    std::string body_;

    // Jacobian function body
    std::string jac_body_;

    // Hessian function body
    std::string hess_body_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_JIT_IMPL_HPP
