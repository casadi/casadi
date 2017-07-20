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
    static Function create(const std::string& name, const Function& f, int n,
                           const Dict& opts);

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
    size_t get_n_in() override { return f_.n_in() + f_.n_out() + f_.n_in();}
    size_t get_n_out() override { return f_.n_out();}
    ///@}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(int i) override;
    std::string get_name_out(int i) override;
    /// @}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    // Evaluate numerically
    void eval(void* mem, const double** arg, double** res, int* iw, double* w) const override;

    /** \brief  Number of perturbed function calls */
    virtual int n_calls() const = 0;

    /** \brief  Calculate perturbed function inputs */
    virtual void perturb(const double** f_arg, double* f_arg_pert, const double** seed) const = 0;

    /** \brief Calculate the finite difference approximation */
    virtual void finalize(const double** f_res, const double* f_res_pert, double** sens) const = 0;

  protected:
    // Constructor (protected, use create function)
    Derivative(const std::string& name, const Function& f, int n, double h);

    // The function which is to be differentiated
    Function f_;

    // Number of directional derivatives
    int n_;

    // Perturbation
    double h_;
  };

  // Forward differences
  class CASADI_EXPORT Forward1 : public Derivative {
  public:
    // Constructor
    Forward1(const std::string& name, const Function& f, int n, double h)
             : Derivative(name, f, n, h) { }

    /** \brief Get type name */
    std::string type_name() const override {return "forward1";}

    /** \brief  Number of function calls */
    int n_calls() const override { return 1;}

    /** \brief Is the scheme using the (nondifferentiated) output? */
    bool uses_output() const override {return true;}

    /** \brief  Calculate perturbed function inputs */
    void perturb(const double** f_arg, double* f_arg_pert, const double** seed) const override;

    /** \brief Calculate the finite difference approximation */
    void finalize(const double** f_res, const double* f_res_pert, double** sens) const override;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_DERIVATIVE_HPP
