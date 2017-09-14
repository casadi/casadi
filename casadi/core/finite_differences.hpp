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


#ifndef CASADI_FINITE_DIFFERENCES_HPP
#define CASADI_FINITE_DIFFERENCES_HPP

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Calculate derivative using finite differences
    * \author Joel Andersson
    * \date 2017
  */
  class CASADI_EXPORT FiniteDiff : public FunctionInternal {
  public:
    // Constructor (protected, use create function)
    FiniteDiff(const std::string& name, int n, double h);

    /** \brief Destructor */
    ~FiniteDiff() override;

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
    int eval(const double** arg, double** res, int* iw, double* w, void* mem) const override;

    /** \brief Is the scheme using the (nondifferentiated) output? */
    bool uses_output() const override {return true;}

    /** \brief Is codegen supported? */
    bool has_codegen() const override { return true;}

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Generate code for the body of the C function */
    void codegen_body(CodeGenerator& g) const override;

  protected:
    // Number of function evaluations needed
    virtual int n_pert() const = 0;

    // Get perturbation expression
    virtual std::string pert(const std::string& k) const = 0;

    // Get perturbation expression
    virtual double pert(int k) const = 0;

    // Calculate finite difference approximation
    virtual void calc_fd(double** yk, double* y0, double* J) const = 0;

    // Codegen finite difference approximation
    virtual void calc_fd(CodeGenerator& g, const std::string& yk,
                         const std::string& y0, const std::string& J) const = 0;

    // Number of directional derivatives
    int n_;

    // Perturbation
    double h_, h2_;

    // Dimensions
    int n_z_, n_y_;

    // Maximum allowed precision
    double h_max_;

    // Input precision
    double eps_in_;

    // Output precision
    double eps_out_;

    // Ratio of roundoff error to truncation error
    double u_aim_;
  };

  /** Calculate derivative using forward differences
    * \author Joel Andersson
    * \date 2017
  */
  class CASADI_EXPORT ForwardDiff : public FiniteDiff {
  public:
    // Constructor
    ForwardDiff(const std::string& name, int n, double h) : FiniteDiff(name, n, h) {}

    /** \brief Destructor */
    ~ForwardDiff() override {}

    /** \brief Get type name */
    std::string class_name() const override {return "ForwardDiff";}

    // Number of function evaluations needed
    int n_pert() const override {return 1;};

    // Get perturbation expression
    std::string pert(const std::string& k) const override {
      return str(h_);
    }

    // Get perturbation expression
    double pert(int k) const override {
      return h_;
    }

    // Calculate finite difference approximation
    void calc_fd(double** yk, double* y0, double* J) const override;

    // Codegen finite difference approximation
    void calc_fd(CodeGenerator& g, const std::string& yk,
                 const std::string& y0, const std::string& J) const override;
  };

  /** Calculate derivative using central differences
    * \author Joel Andersson
    * \date 2017
  */
  class CASADI_EXPORT CentralDiff : public FiniteDiff {
  public:
    // Constructor
    CentralDiff(const std::string& name, int n, double h) : FiniteDiff(name, n, h) {}

    /** \brief Destructor */
    ~CentralDiff() override {}

    /** \brief Get type name */
    std::string class_name() const override {return "CentralDiff";}

    // Number of function evaluations needed
    int n_pert() const override {return 2;};

    // Get perturbation expression
    std::string pert(const std::string& k) const override {
      return "(2*" + k + "-1)*" + str(h_);
    }

    // Get perturbation expression
    double pert(int k) const override {
      return (2*k-1)*h_;
    }

    // Calculate finite difference approximation
    void calc_fd(double** yk, double* y0, double* J) const override;

    // Codegen finite difference approximation
    void calc_fd(CodeGenerator& g, const std::string& yk,
                 const std::string& y0, const std::string& J) const override;

    ///@{
    /** \brief Second order derivatives */
    bool has_forward(int nfwd) const override { return true;}
    Function get_forward(int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_FINITE_DIFFERENCES_HPP
