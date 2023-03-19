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
  // Current control interval
  casadi_int k;
  // Current time
  double t;
  // Next time to be visited by the integrator
  double t_next;
  // Next stop time due to step change in input, continuous
  double t_stop;
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

  /** Helper for a more powerful 'integrator' factory */
  virtual Function create_advanced(const Dict& opts);

  virtual MX algebraic_state_init(const MX& x0, const MX& z0) const { return z0; }
  virtual MX algebraic_state_output(const MX& Z) const { return Z; }

  /** \brief Reset the forward problem

      \identifier{25a} */
  virtual void reset(IntegratorMemory* mem,
    const double* x, const double* z, const double* p) const = 0;

  /** \brief  Find next stop time

      \identifier{25b} */
  casadi_int next_stop(casadi_int k, const double* u) const;

  /** \brief  Advance solution in time

      \identifier{25c} */
  virtual void advance(IntegratorMemory* mem,
    const double* u, double* x, double* z, double* q) const = 0;

  /** \brief Reset the backward problem

      \identifier{25d} */
  virtual void resetB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const = 0;

  /** \brief  Find next stop time

      \identifier{25e} */
  casadi_int next_stopB(casadi_int k, const double* u) const;

  /** \brief Introduce an impulse into the backwards integration at the current time

      \identifier{25f} */
  virtual void impulseB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const = 0;

  /** \brief  Retreat solution in time

      \identifier{25g} */
  virtual void retreat(IntegratorMemory* mem, const double* u,
    double* rx, double* rz, double* rq, double* uq) const = 0;

  /** \brief  evaluate

      \identifier{1m3} */
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

  /** \brief  Print solver statistics

      \identifier{1m4} */
  virtual void print_stats(IntegratorMemory* mem) const {}

  /** \brief  Propagate sparsity forward

      \identifier{1m5} */
  int sp_forward(const bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w, void* mem) const override;

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

  /** \brief Generate a augmented DAE system with \a nfwd forward sensitivities

      \identifier{1mb} */
  template<typename MatType> std::map<std::string, MatType> aug_fwd(casadi_int nfwd) const;

  /** \brief Generate a augmented DAE system with \a nadj adjoint sensitivities

      \identifier{1mc} */
  template<typename MatType> std::map<std::string, MatType> aug_adj(casadi_int nadj) const;

  /// Create sparsity pattern of the extended Jacobian (forward problem)
  Sparsity sp_jac_dae();

  /// Create sparsity pattern of the extended Jacobian (backward problem)
  Sparsity sp_jac_rdae();

  // Sparsity pattern of the extended Jacobians
  Sparsity sp_jac_dae_, sp_jac_rdae_;

  /// New oracle, to replace existing oracle
  Function nonaug_oracle_;

  ///@{
  // Shorthands
  const Sparsity&  t() const { return oracle_.sparsity_in(DYN_T);}
  const Sparsity&  x() const { return oracle_.sparsity_in(DYN_X);}
  const Sparsity&  z() const { return oracle_.sparsity_in(DYN_Z);}
  const Sparsity&  p() const { return oracle_.sparsity_in(DYN_P);}
  const Sparsity&  u() const { return oracle_.sparsity_in(DYN_U);}
  const Sparsity&  q() const { return oracle_.sparsity_out(DYN_QUAD);}
  const Sparsity& rx() const { return oracle_.sparsity_in(DYN_RX);}
  const Sparsity& rz() const { return oracle_.sparsity_in(DYN_RZ);}
  const Sparsity& rp() const { return oracle_.sparsity_in(DYN_RP);}
  const Sparsity& rq() const { return oracle_.sparsity_out(DYN_RQUAD);}
  const Sparsity& uq() const { return oracle_.sparsity_out(DYN_UQUAD);}
  inline casadi_int nt() const { return tout_.size();}
  ///@}

  ///@{
  // Shorthands (new oracle definition)
  const Sparsity&  x1() const { return nonaug_oracle_.sparsity_in(DYN_X);}
  const Sparsity&  z1() const { return nonaug_oracle_.sparsity_in(DYN_Z);}
  const Sparsity&  p1() const { return nonaug_oracle_.sparsity_in(DYN_P);}
  const Sparsity&  u1() const { return nonaug_oracle_.sparsity_in(DYN_U);}
  const Sparsity&  q1() const { return nonaug_oracle_.sparsity_out(DYN_QUAD);}
  const Sparsity& rx1() const { return nonaug_oracle_.sparsity_in(DYN_RX);}
  const Sparsity& rz1() const { return nonaug_oracle_.sparsity_in(DYN_RZ);}
  const Sparsity& rp1() const { return nonaug_oracle_.sparsity_in(DYN_RP);}
  const Sparsity& rq1() const { return nonaug_oracle_.sparsity_out(DYN_RQUAD);}
  const Sparsity& uq1() const { return nonaug_oracle_.sparsity_out(DYN_UQUAD);}
  ///@}

  // Initial time
  double t0_;

  // Output time grid
  std::vector<double> tout_;

  /// Number of sensitivities
  casadi_int nfwd_;

  /// Number of states for the forward integration
  casadi_int nx_, nz_, nq_, nx1_, nz1_, nq1_;

  /// Number of states for the backward integration
  casadi_int nrx_, nrz_, nrq_, nuq_, nrx1_, nrz1_, nrq1_, nuq1_;

  /// Number of forward and backward parameters
  casadi_int np_, nrp_, np1_, nrp1_;

  /// Number of controls
  casadi_int nu_, nu1_;

  /// Number of sensitivities
  casadi_int ns_;

  // Augmented user option
  Dict augmented_options_;

  // Copy of the options
  Dict opts_;

  /// Options
  bool print_stats_;

  // Creator function for internal class
  typedef Integrator* (*Creator)(const std::string& name, const Function& oracle,
    double t0, const std::vector<double>& tout);

  // No static functions exposed
  struct Exposed{ };

  /// Collection of solvers
  static std::map<std::string, Plugin> solvers_;

  /// Infix
  static const std::string infix_;

  /// Convert dictionary to Problem
  template<typename XType>
    static Function map2oracle(const std::string& name,
      const std::map<std::string, XType>& d, const Dict& opts=Dict());


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

  /// Which input is used to calculate a given output in adjoint sensitivity analysis
  static casadi_int adjmap_in(casadi_int i);

  /// Which output is used to calculate a given input in adjoint sensitivity analysis
  static casadi_int adjmap_out(casadi_int i);

 protected:
  /** \brief Deserializing constructor

      \identifier{1mh} */
  explicit Integrator(DeserializingStream& s);
};

/// Input arguments of a forward stepping function
enum FStepIn {
  /// Current time
  FSTEP_T0,
  /// Step size
  FSTEP_H,
  /// State vector
  FSTEP_X0,
  /// Dependent variables
  FSTEP_V0,
  /// Parameter
  FSTEP_P,
  /// Controls
  FSTEP_U,
  /// Number of arguments
  FSTEP_NUM_IN
};

/// Output arguments of a forward stepping function
enum FStepOut {
  /// State vector at next time
  FSTEP_XF,
  /// Dependent variables at next time
  FSTEP_VF,
  /// Quadrature state contribution
  FSTEP_QF,
  /// Number of arguments
  FSTEP_NUM_OUT
};

/// Input arguments of a backward stepping function
enum BStepIn {
  /// Current time
  BSTEP_T0,
  /// Step size
  BSTEP_H,
  /// State vector for backward problem
  BSTEP_RX0,
  /// Dependent variables for backward problem
  BSTEP_RV0,
  /// Parameter vector for backward problem
  BSTEP_RP,
  /// State vector for forward problem
  BSTEP_X,
  /// Dependent variables for forward problem
  BSTEP_V,
  /// Parameter vector for forward problem
  BSTEP_P,
  /// Controls
  BSTEP_U,
  /// Number of arguments
  BSTEP_NUM_IN
};

/// Output arguments of a backward stepping function
enum BStepOut {
  /// State vector for backward problem at the next time
  BSTEP_RXF,
  /// Dependent variables for backward problem at the next time
  BSTEP_RVF,
  /// Quadrature state contribution for backward problem, summing
  BSTEP_RQF,
  /// Quadrature state contribution for backward problem, non-summing
  BSTEP_UQF,
  /// Number of arguments
  BSTEP_NUM_OUT
};

struct CASADI_EXPORT FixedStepMemory : public IntegratorMemory {
  // Work vectors, allocated in base class
  double *x, *z, *rx, *rz, *rq, *x_prev, *rx_prev;

  /// Work vectors, forward problem
  double *v, *p, *u, *q, *v_prev, *q_prev;

  /// Work vectors, backward problem
  double *rv, *rp, *uq, *rv_prev, *rq_prev, *uq_prev;

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

  /** \brief Reset the forward problem

      \identifier{25i} */
  void reset(IntegratorMemory* mem,
    const double* x, const double* z, const double* p) const override;

  /** \brief  Advance solution in time

      \identifier{25j} */
  void advance(IntegratorMemory* mem,
    const double* u, double* x, double* z, double* q) const override;

  /// Reset the backward problem and take time to tf
  void resetB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const override;

  /// Introduce an impulse into the backwards integration at the current time
  void impulseB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const override;

  /** \brief Retreat solution in time

      \identifier{25k} */
  void retreat(IntegratorMemory* mem, const double* u,
    double* rx, double* rz, double* rq, double* uq) const override;

  // Target number of finite elements
  casadi_int nk_target_;

  // Number of steps per control interval
  std::vector<casadi_int> disc_;

  /// Number of dependent variables in the discrete time integration
  casadi_int nv_, nrv_;

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
