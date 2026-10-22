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


#ifndef CASADI_ERK_HPP
#define CASADI_ERK_HPP

#include "casadi/core/runtime_integrator.hpp"
#include <casadi/solvers/casadi_integrator_erk_export.h>

/** \defgroup plugin_Integrator_erk Title
    \par

    Explicit Runge-Kutta integrator for ODEs: a pair of order 5(4) with adaptive or fixed
    steps, or the classical method of order 4 with fixed steps.

    Steps are independent of the output times. Adaptive steps are chosen by error control on
    the nominal states; with fixed steps, each run of steps over a control group (or up to
    the end) is divided into equal steps of at most step_size. Outputs inside a step come
    from a dense output (order 4 and 3 respectively) at no extra cost. Forward sensitivities
    are the exact derivatives of the discrete solution and adjoint sensitivities are the
    exact discrete adjoint of the taped steps. With fixed steps the tape has a known bound,
    so adjoints can be code-generated; with adaptive steps, set max_tape_steps.

    Zero-crossing events (DAE output zero, option transition) are located on the dense
    output of the step in which a function becomes negative; integration restarts from the
    state after the transition. Sensitivities, forward, adjoint and forward-over-adjoint,
    include those of the event times.

    \identifier{2kb} */
/** \pluginsection{Integrator,erk} */

/// \cond INTERNAL
namespace casadi {

  // Memory
  struct CASADI_INTEGRATOR_ERK_EXPORT ErkMemory : public RuntimeIntegratorMemory {
    // Method data, and of the nominal integration whose steps are replayed
    casadi_erk_data<double> e, e_nom;
  };

  /** \brief \pluginbrief{Integrator,erk}

      @copydoc plugin_Integrator_erk
  */
  class CASADI_INTEGRATOR_ERK_EXPORT Erk : public RuntimeIntegrator {
   public:

    /// Constructor
    Erk(const std::string& name, const Function& dae, double t0,
      const std::vector<double>& tout);

    /** \brief  Create a new integrator */
    static Integrator* creator(const std::string& name, const Function& dae,
        double t0, const std::vector<double>& tout) {
      return new Erk(name, dae, t0, tout);
    }

    /// Destructor
    ~Erk() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "erk";}

    // Get name of the class
    std::string class_name() const override { return "Erk";}

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
    void* alloc_mem() const override { return new ErkMemory();}

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<ErkMemory*>(mem);}

    /// A documentation string
    static const std::string meta_doc;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new Erk(s); }

   protected:

    /** \brief Deserializing constructor */
    explicit Erk(DeserializingStream& s);

    // Explicit stages: no linearization
    bool needs_jacobian() const override { return false;}

    // Call the kernel
    int kernel(RuntimeIntegratorMemory* m) const override;

    // Generate code for the problem and data of the method
    void codegen_method(CodeGenerator& g) const override;

   private:
    // Method coefficients
    void set_coefficients();

    // Problem structure of the method
    void set_erk_prob();
    casadi_erk_prob<double> pe_;

    // Options
    std::string method_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_ERK_HPP
