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


#ifndef CASADI_RUNTIME_INTEGRATOR_HPP
#define CASADI_RUNTIME_INTEGRATOR_HPP

#include "integrator_impl.hpp"

/// \cond INTERNAL

namespace casadi {

struct CASADI_EXPORT RuntimeIntegratorMemory : public IntegratorMemory {
  // Reverse communication data, and of the nominal integration whose steps are replayed
  casadi_ivp_data<double> d, d_nom;
  // Tape of unbounded size, if the work vector has none
  std::vector<double> tape;
};

/** \brief Base class of one-step integrators driven by a reverse communication kernel

    The kernel (casadi_ivp.hpp, with the step of a method) requests function evaluations;
    stages, dense output, stops, tape and the discrete adjoint are shared, both in the
    virtual machine and in generated code.

    \identifier{2kc} */
class CASADI_EXPORT RuntimeIntegrator : public Integrator {
 public:

  /// Constructor
  RuntimeIntegrator(const std::string& name, const Function& dae, double t0,
    const std::vector<double>& tout);

  /// Destructor
  ~RuntimeIntegrator() override;

  ///@{
  /** \brief Options

      \identifier{2kd} */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /// Initialize: functions of the kernel tasks
  void init(const Dict& opts) override;

  /** \brief Set the (persistent) work vectors

      \identifier{2ke} */
  void set_work(void* mem, const double**& arg, double**& res,
    casadi_int*& iw, double*& w) const override;

  /** \brief Initalize memory block

      \identifier{2kf} */
  int init_mem(void* mem) const override;

  /// Reset the forward problem at the start or after a change in the controls
  void reset(IntegratorMemory* mem, bool first_call) const override;

  /// Advance solution in time
  int advance_noevent(IntegratorMemory* mem) const override;

  /// Adjoint sensitivities, also with events
  bool has_reverse(casadi_int nadj) const override { return true;}

  /// Reset the backward problem
  void resetB(IntegratorMemory* mem) const override;

  /// Introduce an impulse into the backwards integration at the current time
  void impulseB(IntegratorMemory* mem,
    const double* adj_x, const double* adj_z, const double* adj_q) const override;

  ///@{
  /// Controls and their adjoints: a column per control group
  Sparsity get_sparsity_in(casadi_int i) override;
  Sparsity get_sparsity_out(casadi_int i) override;
  ///@}

  /// Evaluate numerically: a run of steps per control group, events in the kernel
  int eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const override;

  ///@{
  /// Propagate sparsity, with the controls of a group at each of its output intervals
  int sp_forward(const bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w, void* mem) const override;
  int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w,
    void* mem) const override;
  ///@}

  /// Retreat solution in time
  void retreat(IntegratorMemory* mem, const double* u,
    double* adj_x, double* adj_p, double* adj_u) const override;

  /// Get all statistics
  Dict get_stats(void* mem) const override;

  /// Print solver statistics
  void print_stats(IntegratorMemory* mem) const override;

  ///@{
  /** \brief Code generation

      \identifier{2kg} */
  bool has_codegen() const override;
  void codegen_declarations(CodeGenerator& g) const override;
  void codegen_body(CodeGenerator& g) const override;
  ///@}

  /** \brief Serialize an object without type information

      \identifier{2kh} */
  void serialize_body(SerializingStream &s) const override;

 protected:

  /** \brief Deserializing constructor

      \identifier{2ki} */
  explicit RuntimeIntegrator(DeserializingStream& s);

  /// Does the method linearize (tasks EVAL_JAC, EVAL_JAC_ADJ)?
  virtual bool needs_jacobian() const = 0;

  /// Call the kernel of the method: returns 1 while a task is requested
  virtual int kernel(RuntimeIntegratorMemory* m) const = 0;

  /// Generate code for the problem and data of the method, variables pe and e
  virtual void codegen_method(CodeGenerator& g) const = 0;

  /// Integrator at the start of the chain of derivatives
  const RuntimeIntegrator* root() const;

  /// Shared problem structure, from the method coefficients
  void set_ivp_prob();

  /// State after an event as a function of the state before it, see casadi_ivp_event.hpp
  Function jump_function() const;

  /// Problem of the nominal integration whose steps are replayed, if any
  void set_nominal_prob();

  /// Number of control groups, columns of the controls
  casadi_int ngroup() const { return group_start_.size() - 1;}

  /// Work vectors of the shared kernel data
  void alloc_ivp_work();

  /// Carry out the tasks requested by the kernel
  void run(RuntimeIntegratorMemory* m) const;

  /// Capacity of the tape in the work vector (0: unbounded, in the virtual machine only)
  casadi_int tape_capacity() const;

  // Method coefficients: stage arguments, nodes, weights, dense output weights, weights of the
  // error estimate (empty: fixed steps only) and its order
  std::vector<double> a_, alpha_, m_, m_dense_, m_err_;
  casadi_int n_dense_, order_err_;

  // First output interval of each control group, followed by the number of output times
  std::vector<casadi_int> group_start_;

  // Leading entries of x, z, q subject to error control
  casadi_int nx_err_, nz_err_, nq_err_;

  // Options
  double step_size_, abstol_, reltol_, step0_, max_step_size_, min_step_size_;
  casadi_int max_num_steps_, max_tape_steps_;
  bool quad_err_con_;

  // Sparsity of the Jacobian blocks, see casadi_ivp_jac_t (empty: no linearization)
  std::vector<Sparsity> sp_jac_;

  // Shared problem structure
  casadi_ivp_prob<double> p_;

  // Steps replayed from a nominal integration with the functions of the root integrator
  bool replay_;
  casadi_ivp_prob<double> p_nom_;

 private:
  // Source buffers of function arguments and results, see slot_ptr and slot_expr
  enum Src {SRC_T, SRC_Y, SRC_P, SRC_U, SRC_RES, SRC_W, SRC_JAC, SRC_JAC_FWD, SRC_JAC_ADJ,
    SRC_JAC_ADJ_FWD, SRC_VJP_X, SRC_VJP_Z, SRC_VJP_P, SRC_VJP_U, SRC_EV_ONE, SRC_EV_W, SRC_YDOT,
    SRC_EDOT};

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
  int eval_task(RuntimeIntegratorMemory* m, casadi_ivp_data<double>* d,
    const std::vector<std::vector<TaskCall>>& calls) const;

  // Buffer of a slot in the virtual machine, and in generated code
  double* slot_ptr(RuntimeIntegratorMemory* m, casadi_ivp_data<double>* d, const Slot& e) const;
  static std::string slot_expr(const std::string& d, const Slot& e);

  // Reverse communication: evaluate the function requested by the kernel
  int task(RuntimeIntegratorMemory* m) const;

  // Generate code for the tasks of data variable d, as eval_task
  void codegen_task(CodeGenerator& g, const std::string& d,
    const std::vector<std::vector<TaskCall>>& calls, bool in_event = false) const;

  // Generate code for the task loop of the kernel
  void codegen_loop(CodeGenerator& g) const;

  // Generate code for the shared problem structure in variable p
  void codegen_ivp_prob(CodeGenerator& g, const std::string& p, casadi_int tape_size) const;

  // Function calls of the tasks, and of the tasks of the nominal integration
  std::vector<std::vector<TaskCall>> calls_, calls_nom_;
};

} // namespace casadi
/// \endcond

#endif // CASADI_RUNTIME_INTEGRATOR_HPP
