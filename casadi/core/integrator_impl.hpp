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

/// \cond INTERNAL

namespace casadi {

/** \brief Integrator memory

    \identifier{1lp} */
struct CASADI_EXPORT IntegratorMemory : public OracleMemory {
};

/** \brief Internal storage for integrator related data

    @copydoc DAE_doc
    \author Joel Andersson
    \date 2010

    \identifier{1lq} */
class CASADI_EXPORT
Integrator : public OracleFunction, public PluginInterface<Integrator> {
 public:
  /** \brief  Constructor

      \identifier{1lr} */
  Integrator(const std::string& name, const Function& oracle);

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

      \identifier{1lz} */
  virtual void reset(IntegratorMemory* mem, double t,
                      const double* x, const double* z, const double* p) const = 0;

  /** \brief  Advance solution in time

      \identifier{1m0} */
  virtual void advance(IntegratorMemory* mem, double t,
                        double* x, double* z, double* q) const = 0;

  /** \brief Reset the backward problem

      \identifier{1m1} */
  virtual void resetB(IntegratorMemory* mem, double t,
                      const double* rx, const double* rz, const double* rp) const = 0;

  /** \brief  Retreat solution in time

      \identifier{1m2} */
  virtual void retreat(IntegratorMemory* mem, double t,
                        double* rx, double* rz, double* rq) const = 0;

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

  /** \brief  Set stop time for the integration

      \identifier{1m9} */
  virtual void setStopTime(IntegratorMemory* mem, double tf) const;

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

  ///@{
  // Shorthands
  const Sparsity&  t() const { return oracle_.sparsity_in(DE_T);}
  const Sparsity&  x() const { return oracle_.sparsity_in(DE_X);}
  const Sparsity&  z() const { return oracle_.sparsity_in(DE_Z);}
  const Sparsity&  p() const { return oracle_.sparsity_in(DE_P);}
  const Sparsity&  q() const { return oracle_.sparsity_out(DE_QUAD);}
  const Sparsity& rx() const { return oracle_.sparsity_in(DE_RX);}
  const Sparsity& rz() const { return oracle_.sparsity_in(DE_RZ);}
  const Sparsity& rp() const { return oracle_.sparsity_in(DE_RP);}
  const Sparsity& rq() const { return oracle_.sparsity_out(DE_RQUAD);}
  ///@}

  /// Number of states for the forward integration
  casadi_int nx_, nz_, nq_, nx1_, nz1_, nq1_;

  /// Number of states for the backward integration
  casadi_int nrx_, nrz_, nrq_, nrx1_, nrz1_, nrq1_;

  /// Number of forward and backward parameters
  casadi_int np_, nrp_, np1_, nrp1_;

  /// Number of sensitivities
  casadi_int ns_;

  // Time grid
  std::vector<double> grid_;
  casadi_int ngrid_;

  // Augmented user option
  Dict augmented_options_;

  // Copy of the options
  Dict opts_;

  /// One step
  Function onestep_;

  /// Options
  bool print_stats_;

  /// Output the state at the initial time
  bool output_t0_;
  casadi_int ntout_;

  // Creator function for internal class
  typedef Integrator* (*Creator)(const std::string& name, const Function& oracle);

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

 protected:
  /** \brief Deserializing constructor

      \identifier{1mh} */
  explicit Integrator(DeserializingStream& s);
};

struct CASADI_EXPORT FixedStepMemory : public IntegratorMemory {
  // Current time
  double t;

  // Discrete time
  casadi_int k;

  // Current state
  std::vector<double> x, z, p, q, rx, rz, rp, rq;

  // Previous state
  std::vector<double> x_prev, Z_prev, q_prev, rx_prev, RZ_prev, rq_prev;

  /// Algebraic variables for the discrete time integration
  std::vector<double> Z, RZ;

  // Tape
  std::vector<std::vector<double> > x_tape, Z_tape;
};

class CASADI_EXPORT FixedStepIntegrator : public Integrator {
 public:

  /// Constructor
  explicit FixedStepIntegrator(const std::string& name, const Function& dae);

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

  /// Setup F and G
  virtual void setupFG() = 0;

  /** \brief Reset the forward problem

      \identifier{1mm} */
  void reset(IntegratorMemory* mem, double t,
    const double* x, const double* z, const double* p) const override;

  /** \brief  Advance solution in time

      \identifier{1mn} */
  void advance(IntegratorMemory* mem, double t,
    double* x, double* z, double* q) const override;

  /// Reset the backward problem and take time to tf
  void resetB(IntegratorMemory* mem, double t,
    const double* rx, const double* rz, const double* rp) const override;

  /** \brief  Retreat solution in time

      \identifier{1mo} */
  void retreat(IntegratorMemory* mem, double t,
    double* rx, double* rz, double* rq) const override;

  /// Get explicit dynamics
  virtual const Function& getExplicit() const { return F_;}

  /// Get explicit dynamics (backward problem)
  virtual const Function& getExplicitB() const { return G_;}

  // Discrete time dynamics
  Function F_, G_;

  // Number of finite elements
  casadi_int nk_;

  // Time step size
  double h_;

  /// Number of algebraic variables for the discrete time integration
  casadi_int nZ_, nRZ_;

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
  explicit ImplicitFixedStepIntegrator(const std::string& name, const Function& dae);

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

  /// Get explicit dynamics
  const Function& getExplicit() const override { return rootfinder_;}

  /// Get explicit dynamics (backward problem)
  const Function& getExplicitB() const override { return backward_rootfinder_;}

  // Implicit function solver
  Function rootfinder_, backward_rootfinder_;

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
