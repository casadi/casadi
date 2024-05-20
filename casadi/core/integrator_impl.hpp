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


#ifndef CASADI_INTEGRATOR_IMPL_HPP
#define CASADI_INTEGRATOR_IMPL_HPP

#include "integrator.hpp"
#include "oracle_function.hpp"
#include "plugin_interface.hpp"
#include "casadi_enum.hpp"

/// \cond INTERNAL

namespace casadi {

/** \brief Integrator memory

    \identifier{1lp} */
struct CASADI_EXPORT IntegratorMemory : public OracleMemory {
  // Work vectors, forward problem
  double *q, *x, *z, *p, *u, *e, *edot, *old_e, *xdot, *zdot;
  // Work vectors, backward problem
  double *adj_x, *adj_z, *adj_p, *adj_q;
  // Temporary work vectors of length max(nx + nz, nrx, nrz)
  double *tmp1, *tmp2;
  // Current control interval
  casadi_int k;
  // Current time
  double t;
  // Next time to be visited by the integrator
  double t_next;
  // Time not to be exceeded by during integrator integration
  double t_stop;
  // Time at the beginning of the current control interval
  double t_start;
  // Next output time
  double t_next_out;
  // Next stop time due to step change in input
  double t_step;
  // Which events have been triggered
  casadi_int *event_triggered;
  // Do we need to reset the solver?
  bool reset_solver;
  // Number of root-finding iterations for a single event
  casadi_int event_iter;
  // Number of root-finding iterations for all events during an interval
  casadi_int interval_event_iter;
};

/// Memory struct, forward sparsity pattern propagation
struct CASADI_EXPORT SpForwardMem {
  const bvec_t** arg;
  bvec_t** res;
  casadi_int* iw;
  bvec_t* w;
};

/// Memory struct, backward sparsity pattern propagation
struct CASADI_EXPORT SpReverseMem {
  bvec_t** arg;
  bvec_t** res;
  casadi_int* iw;
  bvec_t* w;
};

/** \brief Internal storage for integrator related data


    \author Joel Andersson
    \date 2010

    \identifier{1lq} */
class CASADI_EXPORT
Integrator : public OracleFunction, public PluginInterface<Integrator> {
 public:
  /** \brief  Constructor

      \identifier{1lr} */
  Integrator(const std::string& name, const Function& oracle,
    double t0, const std::vector<double>& tout);

  /** \brief  Destructor

      \identifier{1ls} */
  ~Integrator() override=0;

  ///@{
  /** \brief Number of function inputs and outputs

      \identifier{1lt} */
  size_t get_n_in() override { return INTEGRATOR_NUM_IN;}
  size_t get_n_out() override { return INTEGRATOR_NUM_OUT;}
  ///@}

  /// @{
  /** \brief Sparsities of function inputs and outputs

      \identifier{1lu} */
  Sparsity get_sparsity_in(casadi_int i) override;
  Sparsity get_sparsity_out(casadi_int i) override;
  /// @}

  ///@{
  /** \brief Names of function input and outputs

      \identifier{1lv} */
  std::string get_name_in(casadi_int i) override { return integrator_in(i);}
  std::string get_name_out(casadi_int i) override { return integrator_out(i);}
  /// @}

  /** \brief Initalize memory block

      \identifier{1lw} */
  int init_mem(void* mem) const override;

  ///@{
  /** \brief Options

      \identifier{1lx} */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /** \brief  Initialize

      \identifier{1ly} */
  void init(const Dict& opts) override;

  /** \brief Set the (persistent) work vectors */
  void set_work(void* mem, const double**& arg, double**& res,
    casadi_int*& iw, double*& w) const override;

  /** Helper for a more powerful 'integrator' factory */
  virtual Function create_advanced(const Dict& opts);

  virtual MX algebraic_state_init(const MX& x0, const MX& z0) const { return z0; }
  virtual MX algebraic_state_output(const MX& Z) const { return Z; }

  // Set the quadrature states
  void set_q(IntegratorMemory* m, const double* q) const;

  // Set the differential states
  void set_x(IntegratorMemory* m, const double* x) const;

  // Set the algebraic variables
  void set_z(IntegratorMemory* m, const double* z) const;

  // Set the parameters
  void set_p(IntegratorMemory* m, const double* p) const;

  // Set the controls
  void set_u(IntegratorMemory* m, const double* u) const;

  // Get the quadrature states
  void get_q(IntegratorMemory* m, double* q) const;

  // Get the differential states
  void get_x(IntegratorMemory* m, double* x) const;

  // Get the algebraic variables
  void get_z(IntegratorMemory* m, double* z) const;

  /** \brief  Reset the forward solver at the start or after an event */
  virtual void reset(IntegratorMemory* mem, bool first_call) const {}

  /** \brief  Find next stop time

      \identifier{25b} */
  casadi_int next_stop(casadi_int k, const double* u) const;

  /** \brief Linearize the zero crossing function */
  int calc_edot(IntegratorMemory* m) const;

  /** \brief Predict next event time */
  int predict_events(IntegratorMemory* m) const;

  /** \brief Handle events, if any */
  int handle_events(IntegratorMemory* m) const;

  /** \brief Trigger an event */
  int trigger_event(IntegratorMemory* m, casadi_int* ind) const;

  /** \brief  Advance solution in time

      \identifier{25c} */
  virtual void advance(IntegratorMemory* mem) const = 0;

  /** \brief Reset the backward problem

      \identifier{25d} */
  virtual void resetB(IntegratorMemory* mem) const = 0;

  /** \brief  Find next stop time

      \identifier{25e} */
  casadi_int next_stopB(casadi_int k, const double* u) const;

  /** \brief Introduce an impulse into the backwards integration at the current time

      \identifier{25f} */
  virtual void impulseB(IntegratorMemory* mem,
    const double* adj_x, const double* adj_z, const double* adj_q) const = 0;

  /** \brief  Retreat solution in time

      \identifier{25g} */
  virtual void retreat(IntegratorMemory* mem, const double* u,
    double* adj_x, double* adj_p, double* adj_u) const = 0;

  /** \brief  evaluate

      \identifier{1m3} */
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

  /** \brief  Print solver statistics

      \identifier{1m4} */
  virtual void print_stats(IntegratorMemory* mem) const {}

  /// Forward sparsity pattern propagation through DAE, forward problem
  int fdae_sp_forward(SpForwardMem* m, const bvec_t* x,
    const bvec_t* p, const bvec_t* u, bvec_t* ode, bvec_t* alg) const;

  /// Forward sparsity pattern propagation through quadratures, forward problem
  int fquad_sp_forward(SpForwardMem* m, const bvec_t* x, const bvec_t* z,
    const bvec_t* p, const bvec_t* u, bvec_t* quad) const;

  /// Forward sparsity pattern propagation through DAE, backward problem
  int bdae_sp_forward(SpForwardMem* m, const bvec_t* x, const bvec_t* z,
    const bvec_t* p, const bvec_t* u, const bvec_t* adj_ode, const bvec_t* adj_quad,
    bvec_t* adj_x, bvec_t* adj_z) const;

  /// Forward sparsity pattern propagation through quadratures, backward problem
  int bquad_sp_forward(SpForwardMem* m, const bvec_t* x, const bvec_t* z,
    const bvec_t* p, const bvec_t* u, const bvec_t* adj_ode, const bvec_t* adj_alg,
    const bvec_t* adj_quad, bvec_t* adj_p, bvec_t* adj_u) const;

  /** \brief  Propagate sparsity forward

      \identifier{1m5} */
  int sp_forward(const bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w, void* mem) const override;

  /// Reverse sparsity pattern propagation through DAE, forward problem
  int fdae_sp_reverse(SpReverseMem* m, bvec_t* x,
    bvec_t* p, bvec_t* u, bvec_t* ode, bvec_t* alg) const;

  /// Reverse sparsity pattern propagation through quadratures, forward problem
  int fquad_sp_reverse(SpReverseMem* m, bvec_t* x, bvec_t* z,
    bvec_t* p, bvec_t* u, bvec_t* quad) const;

  /// Reverse sparsity pattern propagation through DAE, backward problem
  int bdae_sp_reverse(SpReverseMem* m, bvec_t* x, bvec_t* z,
    bvec_t* p, bvec_t* u, bvec_t* adj_ode, bvec_t* adj_quad,
    bvec_t* adj_x, bvec_t* adj_z) const;

  /// Reverse sparsity pattern propagation through quadratures, backward problem
  int bquad_sp_reverse(SpReverseMem* m, bvec_t* x, bvec_t* z,
    bvec_t* p, bvec_t* u, bvec_t* adj_ode, bvec_t* adj_alg, bvec_t* adj_quad,
    bvec_t* adj_p, bvec_t* adj_u) const;

  /** \brief  Propagate sparsity backwards

      \identifier{1m6} */
  int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const override;

  ///@{
  /// Is the class able to propagate seeds through the algorithm?
  bool has_spfwd() const override { return true;}
  bool has_sprev() const override { return true;}
  ///@}

  ///@{
  /** \brief Generate a function that calculates \a nfwd forward derivatives

      \identifier{1m7} */
  Function get_forward(casadi_int nfwd, const std::string& name,
                        const std::vector<std::string>& inames,
                        const std::vector<std::string>& onames,
                        const Dict& opts) const override;
  bool has_forward(casadi_int nfwd) const override { return true;}
  ///@}

  ///@{
  /** \brief Generate a function that calculates \a nadj adjoint derivatives

      \identifier{1m8} */
  Function get_reverse(casadi_int nadj, const std::string& name,
                        const std::vector<std::string>& inames,
                        const std::vector<std::string>& onames,
                        const Dict& opts) const override;
  bool has_reverse(casadi_int nadj) const override { return true;}
  ///@}

  /** \brief Set solver specific options to generated augmented integrators

      \identifier{1ma} */
  virtual Dict getDerivativeOptions(bool fwd) const;

  ///@{
  /** \brief Generate the augmented DAE system

      \identifier{261} */
  template<typename MatType> Function get_forward_dae(const std::string& name) const;
  Function augmented_dae() const;
  ///@}

  /// Helper function: Vector has only zeros?
  static bool all_zero(const double* v, casadi_int n);

  /// Helper function, get augmented system Jacobian
  Sparsity sp_jac_aug(const Sparsity& J, const Sparsity& J1) const;

  /// Create sparsity pattern of the extended Jacobian (forward problem)
  Sparsity sp_jac_dae();

  /// Create sparsity pattern of the extended Jacobian (backward problem)
  Sparsity sp_jac_rdae();

  /// Sparsity pattern of the extended Jacobians
  Sparsity sp_jac_dae_, sp_jac_rdae_;

  /// Number of output times
  inline casadi_int nt() const { return tout_.size();}

  ///@{
  /** \brief IO conventions for continuous time dynamics

      \identifier{260} */
  enum DaeOut { DAE_ODE, DAE_ALG, DAE_NUM_OUT};
  static std::vector<std::string> dae_out() { return {"ode", "alg"}; }
  enum QuadOut { QUAD_QUAD, QUAD_NUM_OUT};
  static std::vector<std::string> quad_out() { return {"quad"}; }
  enum BDynIn { BDYN_T, BDYN_X, BDYN_Z, BDYN_P, BDYN_U,
    BDYN_OUT_ODE, BDYN_OUT_ALG, BDYN_OUT_QUAD, BDYN_OUT_ZERO,
    BDYN_ADJ_ODE, BDYN_ADJ_ALG, BDYN_ADJ_QUAD, BDYN_ADJ_ZERO, BDYN_NUM_IN};
  static std::string bdyn_in(casadi_int i);
  static std::vector<std::string> bdyn_in();
  enum BDynOut { BDYN_ADJ_T, BDYN_ADJ_X, BDYN_ADJ_Z, BDYN_ADJ_P, BDYN_ADJ_U, BDYN_NUM_OUT};
  static std::string bdyn_out(casadi_int i);
  static std::vector<std::string> bdyn_out();
  enum DAEBOut { BDAE_ADJ_X, BDAE_ADJ_Z, BDAE_NUM_OUT};
  static std::vector<std::string> bdae_out() { return {"adj_x", "adj_z"}; }
  enum QuadBOut { BQUAD_ADJ_P, BQUAD_ADJ_U, BQUAD_NUM_OUT};
  static std::vector<std::string> bquad_out() { return {"adj_p", "adj_u"}; }
  ///@}

  /// Initial time
  double t0_;

  /// Output time grid
  std::vector<double> tout_;

  /// Number of sensitivities
  casadi_int nfwd_, nadj_;

  /// Backwards DAE function
  Function rdae_;

  /// Number of states for the forward integration
  casadi_int nx_, nz_, nq_, nx1_, nz1_, nq1_;

  /// Number of states for the backward integration
  casadi_int nrx_, nrz_, nrq_, nuq_, nrx1_, nrz1_, nrq1_, nuq1_;

  /// Number of forward and backward parameters
  casadi_int np_, nrp_, np1_, nrp1_;

  /// Number of controls
  casadi_int nu_, nu1_;

  /// Number of of zero-crossing functions
  casadi_int ne_;

  /// Length of the tmp1, tmp2 vectors
  casadi_int ntmp_;

  // Nominal values for states
  std::vector<double> nom_x_, nom_z_;

  /// Augmented user option
  Dict augmented_options_;

  /// Copy of the options
  Dict opts_;

  /// Options
  bool print_stats_;

  /// Function to be called at state events
  Function event_transition_;

  /// Maximum number of event iterations for a single event
  casadi_int max_event_iter_;

  /// Maximum number of total event iterations during an interval
  casadi_int max_interval_event_iter_;

  /// Termination tolerance for the event iteration
  double event_tol_;

  /// Acceptable tolerance for the event iteration
  double event_acceptable_tol_;

  // Creator function for internal class
  typedef Integrator* (*Creator)(const std::string& name, const Function& oracle,
    double t0, const std::vector<double>& tout);

  // No static functions exposed
  struct Exposed{ };

  /// Collection of solvers
  static std::map<std::string, Plugin> solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

  /// Infix
  static const std::string infix_;

  /// Convert dictionary to Problem
  template<typename XType>
  static Function map2oracle(const std::string& name, const std::map<std::string, XType>& d);

  /** \brief Serialize an object without type information

      \identifier{1md} */
  void serialize_body(SerializingStream &s) const override;
  /** \brief Serialize type information

      \identifier{1me} */
  void serialize_type(SerializingStream &s) const override;

  /** \brief Deserialize into MX

      \identifier{1mf} */
  static ProtoFunction* deserialize(DeserializingStream& s);

  /** \brief String used to identify the immediate FunctionInternal subclass

      \identifier{1mg} */
  std::string serialize_base_function() const override { return "Integrator"; }

  /// Is an input repeated for each grid point?
  static bool grid_in(casadi_int i);

  /// Is an output repeated for each grid point?
  static bool grid_out(casadi_int i);

  /// Which output is used to calculate a given input in adjoint sensitivity analysis
  static casadi_int adjmap_out(casadi_int i);

 protected:
  /** \brief Deserializing constructor

      \identifier{1mh} */
  explicit Integrator(DeserializingStream& s);
};

/// Input arguments of a forward stepping function
enum StepIn {
  /// Current time
  STEP_T,
  /// Step size
  STEP_H,
  /// State vector
  STEP_X0,
  /// Dependent variables
  STEP_V0,
  /// Parameter
  STEP_P,
  /// Controls
  STEP_U,
  /// Number of arguments
  STEP_NUM_IN
};

/// Output arguments of a forward stepping function
enum StepOut {
  /// State vector at next time
  STEP_XF,
  /// Dependent variables at next time
  STEP_VF,
  /// Quadrature state contribution
  STEP_QF,
  /// Number of arguments
  STEP_NUM_OUT
};

/// Input arguments of a backward stepping function
enum BStepIn {
  BSTEP_T,
  BSTEP_H,
  BSTEP_X0,
  BSTEP_V0,
  BSTEP_P,
  BSTEP_U,
  BSTEP_OUT_XF,
  BSTEP_OUT_VF,
  BSTEP_OUT_QF,
  BSTEP_ADJ_XF,
  BSTEP_ADJ_VF,
  BSTEP_ADJ_QF,
  BSTEP_NUM_IN
};

/// Output arguments of a backward stepping function
enum BStepOut {
  BSTEP_ADJ_T,
  BSTEP_ADJ_H,
  BSTEP_ADJ_X0,
  BSTEP_ADJ_V0,
  BSTEP_ADJ_P,
  BSTEP_ADJ_U,
  BSTEP_NUM_OUT
};

struct CASADI_EXPORT FixedStepMemory : public IntegratorMemory {
  /// Work vectors, forward problem
  double *v, *v_prev, *q_prev;

  /// Work vectors, backward problem
  double *rv, *adj_u, *adj_p_prev, *adj_u_prev;

  /// State and dependent variables at all times
  double *x_tape, *v_tape;
};

class CASADI_EXPORT FixedStepIntegrator : public Integrator {
 public:

  /// Constructor
  explicit FixedStepIntegrator(const std::string& name, const Function& dae,
    double t0, const std::vector<double>& tout);

  /// Destructor
  ~FixedStepIntegrator() override;

  ///@{
  /** \brief Options

      \identifier{1mi} */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /// Initialize stage
  void init(const Dict& opts) override;

  /** \brief Set the (persistent) work vectors

      \identifier{25h} */
  void set_work(void* mem, const double**& arg, double**& res,
    casadi_int*& iw, double*& w) const override;

  /** Helper for a more powerful 'integrator' factory */
  Function create_advanced(const Dict& opts) override;

  /** \brief Create memory block

      \identifier{1mj} */
  void* alloc_mem() const override { return new FixedStepMemory();}

  /** \brief Initalize memory block

      \identifier{1mk} */
  int init_mem(void* mem) const override;

  /** \brief Free memory block

      \identifier{1ml} */
  void free_mem(void *mem) const override { delete static_cast<FixedStepMemory*>(mem);}

  /// Setup step functions
  virtual void setup_step() = 0;

  /** \brief  Reset the forward solver at the start or after an event */
  void reset(IntegratorMemory* mem, bool first_call) const override;

  /** \brief  Advance solution in time

      \identifier{25j} */
  void advance(IntegratorMemory* mem) const override;

  /// Reset the backward problem and take time to tf
  void resetB(IntegratorMemory* mem) const override;

  /// Introduce an impulse into the backwards integration at the current time
  void impulseB(IntegratorMemory* mem,
    const double* adj_x, const double* adj_z, const double* adj_q) const override;

  /** \brief Retreat solution in time

      \identifier{25k} */
  void retreat(IntegratorMemory* mem, const double* u,
    double* adj_x, double* adj_p, double* adj_u) const override;

  /// Take integrator step forward
  void stepF(FixedStepMemory* m, double t, double h,
    const double* x0, const double* v0, double* xf, double* vf, double* qf) const;

  /// Take integrator step backward
  void stepB(FixedStepMemory* m, double t, double h,
    const double* x0, const double* xf, const double* vf,
    const double* adj_xf, const double* rv0,
    double* adj_x0, double* adj_p, double* adj_u) const;

  // Target number of finite elements
  casadi_int nk_target_;

  // Number of steps per control interval
  std::vector<casadi_int> disc_;

  /// Number of dependent variables in the discrete time integration
  casadi_int nv_, nv1_, nrv_, nrv1_;

  /** \brief Serialize an object without type information

      \identifier{1mp} */
  void serialize_body(SerializingStream &s) const override;

protected:
  /** \brief Deserializing constructor

      \identifier{1mq} */
  explicit FixedStepIntegrator(DeserializingStream& s);
};

class CASADI_EXPORT ImplicitFixedStepIntegrator : public FixedStepIntegrator {
 public:

  /// Constructor
  explicit ImplicitFixedStepIntegrator(const std::string& name, const Function& dae,
    double t0, const std::vector<double>& tout);

  /// Destructor
  ~ImplicitFixedStepIntegrator() override;

  ///@{
  /** \brief Options

      \identifier{1mr} */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /// Initialize stage
  void init(const Dict& opts) override;

  /** \brief Serialize an object without type information

      \identifier{1ms} */
  void serialize_body(SerializingStream &s) const override;

protected:
  /** \brief Deserializing constructor

      \identifier{1mt} */
  explicit ImplicitFixedStepIntegrator(DeserializingStream& s);
};

} // namespace casadi
/// \endcond

#endif // CASADI_INTEGRATOR_IMPL_HPP
