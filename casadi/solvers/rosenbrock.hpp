/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl,
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


#ifndef CASADI_ROSENBROCK_HPP
#define CASADI_ROSENBROCK_HPP

#include "casadi/core/runtime_integrator.hpp"
#include <casadi/solvers/casadi_integrator_rosenbrock_export.h>

/** \defgroup plugin_Integrator_rosenbrock Title
    \par

    Rosenbrock-Wanner integrator for stiff ODEs and semi-explicit index-1 DAEs, with
    adaptive or fixed steps.

    An 8-stage, stiffly accurate, L-stable method of order 5 with an embedded method of
    order 4. Each step needs one Jacobian and one sparse QR factorization; the method is
    linearly implicit, without Newton iterations. Steps cross the output times within a
    control group, which are interpolated by a dense output of order 4. Step sizes are
    chosen from the nominal solution only, so forward sensitivities are the exact
    derivatives of the discrete solution (internal numerical differentiation) and adjoint
    sensitivities are the exact discrete adjoint of the taped steps. With step_size, each
    run of steps over a control group is divided into equal steps instead, and the tape has
    a known bound.

    Zero-crossing events (DAE output zero, option transition) are located on the dense
    output of the step in which a function becomes negative; integration restarts from the
    state after the transition. Sensitivities, forward, adjoint and forward-over-adjoint,
    include those of the event times.

    \identifier{2ka} */
/** \pluginsection{Integrator,rosenbrock} */

/// \cond INTERNAL
namespace casadi {

  // Memory
  struct CASADI_INTEGRATOR_ROSENBROCK_EXPORT RosenbrockMemory : public RuntimeIntegratorMemory {
    // Method data, and of the nominal integration whose steps are replayed
    casadi_rosenbrock_data<double> r, r_nom;
  };

  /** \brief \pluginbrief{Integrator,rosenbrock}

      @copydoc plugin_Integrator_rosenbrock
  */
  class CASADI_INTEGRATOR_ROSENBROCK_EXPORT Rosenbrock : public RuntimeIntegrator {
   public:

    /// Constructor
    Rosenbrock(const std::string& name, const Function& dae, double t0,
      const std::vector<double>& tout);

    /** \brief  Create a new integrator */
    static Integrator* creator(const std::string& name, const Function& dae,
        double t0, const std::vector<double>& tout) {
      return new Rosenbrock(name, dae, t0, tout);
    }

    /// Destructor
    ~Rosenbrock() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "rosenbrock";}

    // Get name of the class
    std::string class_name() const override { return "Rosenbrock";}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Initialize
    void init(const Dict& opts) override;

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
      casadi_int*& iw, double*& w) const override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new RosenbrockMemory();}

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<RosenbrockMemory*>(mem);}

    /// A documentation string
    static const std::string meta_doc;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new Rosenbrock(s); }

   protected:

    /** \brief Deserializing constructor */
    explicit Rosenbrock(DeserializingStream& s);

    // The iteration matrix needs the Jacobian
    bool needs_jacobian() const override { return true;}

    // Call the kernel
    int kernel(RuntimeIntegratorMemory* m) const override;

    // Generate code for the problem and data of the method
    void codegen_method(CodeGenerator& g) const override;

   private:
    // Sparsity of the iteration matrix and symbolic QR factorizations, from the Jacobian blocks
    void set_iteration_matrix();

    // Method coefficients in transformed form
    void set_coefficients();

    // Problem structure of the method, and of the nominal integration if replayed
    void set_rosenbrock_prob();

    // Generate code for the problem structure of the method in variable pe, shared part p
    void codegen_rosenbrock_prob(CodeGenerator& g, const std::string& pe,
      const std::string& p) const;

    // Method coefficients, transformed form, beyond the shared ones
    std::vector<double> c_, gamma_sum_;

    // Problem structure of the method, and of the nominal integration
    casadi_rosenbrock_prob<double> pr_, pr_nom_;

    // Options
    double init_tol_;
    casadi_int init_max_iter_;

    // Iteration matrix: sparsity, nonzero positions of the Jacobian blocks and x-diagonal
    Sparsity sp_iter_;
    std::vector<casadi_int> iter_map_;

    // Sparse QR of the iteration matrix and of d(alg)/dz
    Sparsity sp_iter_v_, sp_iter_r_, sp_alg_v_, sp_alg_r_;
    std::vector<casadi_int> iter_prinv_, iter_pc_, alg_prinv_, alg_pc_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_ROSENBROCK_HPP
