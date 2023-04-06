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

#ifndef CASADI_CALLBACK_INTERNAL_HPP
#define CASADI_CALLBACK_INTERNAL_HPP

#include "callback.hpp"
#include "function_internal.hpp"

namespace casadi {

  class CASADI_EXPORT CallbackInternal : public FunctionInternal {
    friend class CallbackFunction;
  public:

    /** \brief Constructor

        \identifier{17w} */
    explicit CallbackInternal(const std::string& name, Callback* self);

    /** \brief Destructor

        \identifier{17x} */
    ~CallbackInternal() override;

    /** \brief Get type name

        \identifier{17y} */
    std::string class_name() const override {return "CallbackInternal";}

    ///@{
    /** \brief Number of function inputs and outputs

        \identifier{17z} */
    size_t get_n_in() override;
    size_t get_n_out() override;
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs

        \identifier{180} */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs

        \identifier{181} */
    std::string get_name_in(casadi_int i) override;
    std::string get_name_out(casadi_int i) override;
    /// @}

    /** \brief  Initialize

        \identifier{182} */
    void init(const Dict& opts) override;

    /** \brief Finalize the object creation

        \identifier{183} */
    void finalize() override;

    ///@{
    /** \brief Evaluate with DM matrices

        \identifier{184} */
    std::vector<DM> eval_dm(const std::vector<DM>& arg) const override;
    bool has_eval_dm() const override { return !has_eval_buffer_;}
    ///@}

    /** \brief  Evaluate numerically

        \identifier{185} */
    virtual int eval(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const override;
    bool has_eval_buffer() const;

    /** \brief Do the derivative functions need nondifferentiated outputs?

        \identifier{186} */
    bool uses_output() const override;

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements

        \identifier{187} */
    bool has_jacobian() const override;
    Function get_jacobian(const std::string& name,
                          const std::vector<std::string>& inames,
                          const std::vector<std::string>& onames,
                          const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Return sparsity of Jacobian of an output respect to an input

        \identifier{188} */
    bool has_jac_sparsity(casadi_int oind, casadi_int iind) const override;
    Sparsity get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const override;
    ///@}

    ///@{
    /** \brief Return function that calculates forward derivatives

        \identifier{189} */
    bool has_forward(casadi_int nfwd) const override;
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Return function that calculates adjoint derivatives

        \identifier{18a} */
    bool has_reverse(casadi_int nadj) const override;
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    /** \brief Pointer to the public class

        \identifier{18b} */
    Callback* self_;

    // For buffered evaluation
    std::vector<casadi_int> sizes_arg_, sizes_res_;
    bool has_eval_buffer_;
  };

} // namespace casadi

#endif // CASADI_CALLBACK_INTERNAL_HPP
