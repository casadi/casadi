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


#ifndef CASADI_DERIVATIVE_HPP
#define CASADI_DERIVATIVE_HPP

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Evaluate in parallel
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT Derivative : public FunctionInternal {
  public:
    // Create function (use instead of constructor)
    static Function create(const std::string& name, const Dict& opts);

    /** \brief Destructor */
    ~Derivative() override;

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(int i) override;
    Sparsity get_sparsity_out(int i) override;
    /// @}

    /** \brief Get default input value */
    double default_in(int ind) const override;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override;
    size_t get_n_out() override;
    ///@}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(int i) override;
    std::string get_name_out(int i) override;
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    // Evaluate numerically
    void eval(void* mem, const double** arg, double** res, int* iw, double* w) const override;

    /** \brief  Number of perturbed function calls */
    virtual int n_calls() const = 0;

    /** \brief Function to be called */
    virtual const Function& f() const { return derivative_of_;}

    /** \brief  Calculate perturbed function inputs */
    virtual void perturb(const double** f_arg, double* f_arg_pert) const = 0;

    /** \brief Calculate the finite difference approximation */
    virtual void finalize(const double** f_res, const double* f_res_pert, double* jac) const = 0;

  protected:
    // Constructor (protected, use create function)
    Derivative(const std::string& name, double h);

    // Perturbation
    double h_, h2_;
  };

  // Forward differences, first order
  class CASADI_EXPORT Forward : public Derivative {
  public:
    // Constructor
    Forward(const std::string& name, double h) : Derivative(name, h) { }

    /** \brief Destructor */
    ~Forward() override {}

    /** \brief Get type name */
    std::string type_name() const override {return "forward";}

    /** \brief  Number of function calls */
    int n_calls() const override { return 2;}

    /** \brief Is the scheme using the (nondifferentiated) output? */
    bool uses_output() const override {return false;}

    /** \brief  Calculate perturbed function inputs */
    void perturb(const double** f_arg, double* f_arg_pert) const override;

    /** \brief Calculate the finite difference approximation */
    void finalize(const double** f_res, const double* f_res_pert, double* jac) const override;
  };

  // Central differences, first order
  class CASADI_EXPORT Central : public Derivative {
  public:
    // Constructor
    Central(const std::string& name, double h) : Derivative(name, h) { }

    /** \brief Destructor */
    ~Central() override {}

    /** \brief Get type name */
    std::string type_name() const override {return "central";}

    /** \brief  Number of function calls */
    int n_calls() const override { return 2;}

    /** \brief Is the scheme using the (nondifferentiated) output? */
    bool uses_output() const override {return false;}

    /** \brief  Calculate perturbed function inputs */
    void perturb(int i, const double** f_arg, double* f_arg_pert) const override;

    /** \brief Calculate the finite difference approximation */
    void finalize(const double** f_res, const double* f_res_pert, double* jac) const override;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_DERIVATIVE_HPP
