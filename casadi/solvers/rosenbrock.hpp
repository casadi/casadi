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

#include "casadi/core/integrator_impl.hpp"
#include <casadi/solvers/casadi_integrator_rosenbrock_export.h>

/** \defgroup plugin_Integrator_rosenbrock Title
    \par

    Adaptive Rosenbrock-Wanner integrator for stiff ODEs and semi-explicit index-1 DAEs.

    An 8-stage, stiffly accurate, L-stable method of order 5 with an embedded
    method of order 4. Each step needs one Jacobian and one sparse QR factorization;
    the method is linearly implicit, without Newton iterations. Steps cross output
    times, which are interpolated by a dense output of order 4, unless the
    controls change there. Step sizes are chosen from the nominal solution
    only, so forward sensitivities are the exact derivatives of the discrete solution
    (internal numerical differentiation) and adjoint sensitivities are the exact
    discrete adjoint of the taped steps.

*/
/** \pluginsection{Integrator,rosenbrock} */

/// \cond INTERNAL
namespace casadi {

  // Memory
  struct CASADI_INTEGRATOR_ROSENBROCK_EXPORT RosenbrockMemory : public IntegratorMemory {
    // Reverse communication data, and of the nominal integration whose steps are replayed
    casadi_rosenbrock_data<double> d, d_nom;
    // Tape of unbounded size, if the work vector has none
    std::vector<double> tape;
  };

  /** \brief \pluginbrief{Integrator,rosenbrock}

      @copydoc plugin_Integrator_rosenbrock
  */
  class CASADI_INTEGRATOR_ROSENBROCK_EXPORT Rosenbrock : public Integrator {
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

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<RosenbrockMemory*>(mem);}

    /// Reset the forward problem at the start or after a change in the controls
    void reset(IntegratorMemory* mem, bool first_call) const override;

    /// Advance solution in time
    int advance_noevent(IntegratorMemory* mem) const override;

    /// Reset the backward problem
    void resetB(IntegratorMemory* mem) const override;

    /// Introduce an impulse into the backwards integration at the current time
    void impulseB(IntegratorMemory* mem,
      const double* adj_x, const double* adj_z, const double* adj_q) const override;

    /// Next output interval at which the nominal controls change
    casadi_int next_stop(casadi_int k, const double* u) const override;

    /// Retreat solution in time
    void retreat(IntegratorMemory* mem, const double* u,
      double* adj_x, double* adj_p, double* adj_u) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    ///@{
    /** \brief Code generation */
    bool has_codegen() const override;
    void codegen_declarations(CodeGenerator& g) const override;
    void codegen_body(CodeGenerator& g) const override;
    ///@}

    /// Print solver statistics
    void print_stats(IntegratorMemory* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new Rosenbrock(s); }

   protected:

    /** \brief Deserializing constructor */
    explicit Rosenbrock(DeserializingStream& s);

   private:
    // Set up the method coefficients and the reverse communication problem
    void set_rosenbrock_prob();

    // Integrator at the start of the chain of derivatives
    const Rosenbrock* root() const;

    // Problem of the nominal integration whose steps are replayed
    void set_nominal_prob();

    // Carry out the tasks requested by the kernel
    void run(RosenbrockMemory* m) const;

    // Source buffers of function arguments and results, see slot_ptr and slot_expr
    enum Src {SRC_T, SRC_Y, SRC_P, SRC_U, SRC_RES, SRC_W, SRC_JAC, SRC_JAC_FWD, SRC_JAC_ADJ,
      SRC_JAC_ADJ_FWD, SRC_VJP_X, SRC_VJP_Z, SRC_VJP_P, SRC_VJP_U};

    // Argument or result i of a function call: buffer src at offset off
    struct Slot {
      casadi_int i;
      Src src;
      casadi_int off;
    };

    // Function call of a task: arguments and results not listed are null
    struct TaskCall {
      std::string fname;
      casadi_int n_in, n_out;
      std::vector<Slot> arg, res;
    };

    // For each task of the kernel, the function calls (function names with a suffix)
    std::vector<std::vector<TaskCall>> task_calls(const std::string& suffix) const;

    // Carry out a task of data d with a task table
    int eval_task(RosenbrockMemory* m, casadi_rosenbrock_data<double>* d,
      const std::vector<std::vector<TaskCall>>& calls) const;

    // Buffer of a slot in the virtual machine, and in generated code
    double* slot_ptr(RosenbrockMemory* m, casadi_rosenbrock_data<double>* d,
      const Slot& e) const;
    static std::string slot_expr(const std::string& d, const Slot& e);

    // Name of a task in generated code
    static std::string task_name(casadi_int task);

    // Reverse communication: evaluate the function requested by the kernel
    int task(RosenbrockMemory* m) const;

    // Generate code for the tasks of data variable d, as eval_task
    void codegen_task(CodeGenerator& g, const std::string& d,
      const std::vector<std::vector<TaskCall>>& calls) const;

    // Generate code for the task loop of the kernel
    void codegen_loop(CodeGenerator& g) const;

    // Generate code for the problem structure of the kernel, variable p
    void codegen_prob(CodeGenerator& g, const std::string& p) const;

    // Capacity of the tape in the work vector (0: unbounded, in the virtual machine only)
    casadi_int tape_capacity() const;

    // Method coefficients, transformed form
    std::vector<double> a_, c_, alpha_, gamma_sum_, m_, m_err_, m_dense_;

    // Problem structure for the reverse communication kernel
    casadi_rosenbrock_prob<double> p_;

    // Steps replayed from a nominal integration with the functions of the root integrator
    bool replay_;
    casadi_rosenbrock_prob<double> p_nom_;

    // Function calls of the tasks, and of the tasks of the nominal integration
    std::vector<std::vector<TaskCall>> calls_, calls_nom_;

    // Options
    double abstol_, reltol_, step0_, max_step_size_, min_step_size_, init_tol_;
    casadi_int max_num_steps_, fixed_steps_, init_max_iter_, max_tape_steps_;
    bool quad_err_con_;

    // Dimensions of the nominal problem, subject to error control
    casadi_int nx_err_, nz_err_, nq_err_, nu_err_;

    // Sparsity of the Jacobian blocks, see casadi_rosenbrock_jac_t
    std::vector<Sparsity> sp_jac_;

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
