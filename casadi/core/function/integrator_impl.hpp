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
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief Integrator memory */
  struct CASADI_EXPORT IntegratorMemory : public WorkMemory {
  };

  /** \brief Internal storage for integrator related data

      @copydoc DAE_doc
      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT
  Integrator : public FunctionInternal, public PluginInterface<Integrator> {
  public:
    /** \brief  Constructor */
    Integrator(const std::string& name, Oracle* dae);

    /** \brief  Destructor */
    virtual ~Integrator()=0;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return INTEGRATOR_NUM_IN;}
    virtual size_t get_n_out() { return INTEGRATOR_NUM_OUT;}
    ///@}

   /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i);
    virtual Sparsity get_sparsity_out(int i);
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(int i) { return integrator_in(i);}
    virtual std::string get_name_out(int i) { return integrator_out(i);}
    /// @}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /** \brief Set the work vectors */
    virtual void set_temp(void* mem, const double** arg, double** res,
                          int* iw, double* w) const;

    /** \brief Reset the forward problem */
    virtual void reset(IntegratorMemory* mem, double t,
                       const double* x, const double* z, const double* p) const = 0;

    /** \brief  Advance solution in time */
    virtual void advance(IntegratorMemory* mem, double t,
                         double* x, double* z, double* q) const = 0;

    /** \brief Reset the backward problem */
    virtual void resetB(IntegratorMemory* mem, double t,
                        const double* rx, const double* rz, const double* rp) const = 0;

    /** \brief  Retreat solution in time */
    virtual void retreat(IntegratorMemory* mem, double t,
                         double* rx, double* rz, double* rq) const = 0;

    /** \brief  evaluate */
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief  Print solver statistics */
    virtual void printStats(IntegratorMemory* mem, std::ostream &stream) const {}

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /// Is the class able to propagate seeds through the algorithm?
    virtual bool spCanEvaluate(bool fwd) { return true;}

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function get_forward(const std::string& name, int nfwd, Dict& opts);
    virtual int get_n_forward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function get_reverse(const std::string& name, int nadj, Dict& opts);
    virtual int get_n_reverse() const { return 64;}
    ///@}

    /** \brief  Set stop time for the integration */
    virtual void setStopTime(IntegratorMemory* mem, double tf) const;

    // Helper structure
    struct AugOffset {
      std::vector<int> x, z, q, p, rx, rz, rq, rp;
    };

    /** \brief Set solver specific options to generated augmented integrators */
    virtual Dict getDerivativeOptions(bool fwd);

    /** \brief Generate a augmented DAE system with \a nfwd forward sensitivities  */
    template<typename MatType>
      std::map<std::string, MatType> aug_fwd(int nfwd, AugOffset& offset);

    /** \brief Generate a augmented DAE system with \a nfwd forward sensitivities
    * and \a nadj adjoint sensitivities */
    template<typename MatType>
      std::map<std::string, MatType> aug_adj(int nadj, AugOffset& offset);

    /// Get offsets in augmented problem
    AugOffset getAugOffset(int nfwd, int nadj);

    /// Create sparsity pattern of the extended Jacobian (forward problem)
    Sparsity spJacF();

    /// Create sparsity pattern of the extended Jacobian (backward problem)
    Sparsity spJacG();

    // Sparities
    Sparsity t_, x_, z_, p_, q_, rx_, rz_, rp_, rq_;

    ///@{
    // Shorthands
    const Sparsity&  t() { return t_;}
    const Sparsity&  x() { return x_;}
    const Sparsity&  z() { return z_;}
    const Sparsity&  p() { return p_;}
    const Sparsity&  q() { return q_;}
    const Sparsity& rx() { return rx_;}
    const Sparsity& rz() { return rz_;}
    const Sparsity& rp() { return rp_;}
    const Sparsity& rq() { return rq_;}
    ///@}

    /// Number of states for the forward integration
    int nx_, nz_, nq_;

    /// Number of states for the backward integration
    int nrx_, nrz_, nrq_;

    /// Number of forward and backward parameters
    int np_, nrp_;

    // Time grid
    std::vector<double> grid_;
    int ngrid_;

    // Augmented user option
    Dict augmented_options_;

    // Copy of the options
    Dict opts_;

    // Dae
    Oracle* dae_;

    /// One step
    Function onestep_;

    /// ODE/DAE forward integration function
    Function f_;

    /// ODE/DAE backward integration function, if any
    Function g_;

    /// Integrator for sparsity pattern propagation
    Function linsol_f_, linsol_g_;

    /// Options
    bool print_stats_;

    /// Output the state at the initial time
    bool output_t0_;
    int ntout_;

    // Creator function for internal class
    typedef Integrator* (*Creator)(const std::string& name, Oracle* dae);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Convert dictionary to Problem
    template<typename XType>
      static Oracle* map2problem(const std::map<std::string, XType>& d);

    /// Get the (legacy) dae forward function
    template<typename XType>
      static Oracle* fun2problem(Function f, Function g=Function());
  };

  template<typename XType>
  Oracle* Integrator::map2problem(const std::map<std::string, XType>& d) {
    std::vector<XType> de_in(DE_NUM_IN), de_out(DE_NUM_OUT);
    for (auto&& i : d) {
      if (i.first=="t") {
        de_in[DE_T]=i.second;
      } else if (i.first=="x") {
        de_in[DE_X]=i.second;
      } else if (i.first=="z") {
        de_in[DE_Z]=i.second;
      } else if (i.first=="p") {
        de_in[DE_P]=i.second;
      } else if (i.first=="rx") {
        de_in[DE_RX]=i.second;
      } else if (i.first=="rz") {
        de_in[DE_RZ]=i.second;
      } else if (i.first=="rp") {
        de_in[DE_RP]=i.second;
      } else if (i.first=="ode") {
        de_out[DE_ODE]=i.second;
      } else if (i.first=="alg") {
        de_out[DE_ALG]=i.second;
      } else if (i.first=="quad") {
        de_out[DE_QUAD]=i.second;
      } else if (i.first=="rode") {
        de_out[DE_RODE]=i.second;
      } else if (i.first=="ralg") {
        de_out[DE_RALG]=i.second;
      } else if (i.first=="rquad") {
        de_out[DE_RQUAD]=i.second;
      } else {
        casadi_error("No such field: " + i.first);
      }
    }
    return Oracle::construct(de_in, de_out, DE_INPUTS, DE_OUTPUTS);
  }

  template<typename XType>
  Oracle* Integrator::fun2problem(Function f, Function g) {
    std::vector<XType> dae_in(DE_NUM_IN), dae_out(DE_NUM_OUT);
    std::vector<XType> v = XType::get_input(f), vf=v, vg=v;
    dae_in[DE_T] = v[DAE_T];
    dae_in[DE_X] = v[DAE_X];
    dae_in[DE_Z] = v[DAE_Z];
    dae_in[DE_P] = v[DAE_P];
    v = f(v);
    dae_out[DE_ODE] = v[DAE_ODE];
    dae_out[DE_ALG] = v[DAE_ALG];
    dae_out[DE_QUAD] = v[DAE_QUAD];
    if (!g.is_null()) {
      v = XType::get_input(g);
      dae_in[DE_RX] = v[RDAE_RX];
      dae_in[DE_RZ] = v[RDAE_RZ];
      dae_in[DE_RP] = v[RDAE_RP];
      vg[DAE_T] = v[RDAE_T];
      vg[DAE_X] = v[RDAE_X];
      vg[DAE_Z] = v[RDAE_Z];
      vg[DAE_P] = v[RDAE_P];
      v = substitute(g(v), vg, vf);
      dae_out[DE_RODE] = v[RDAE_ODE];
      dae_out[DE_RALG] = v[RDAE_ALG];
      dae_out[DE_RQUAD] = v[RDAE_QUAD];
    }
    return Oracle::construct(dae_in, dae_out, DE_INPUTS, DE_OUTPUTS);
  }

  struct CASADI_EXPORT FixedStepMemory : public IntegratorMemory {
    // Current time
    double t;

    // Discrete time
    int k;

    // Current state
    std::vector<double> x, z, p, q, rx, rz, rp, rq;

    // Previous state
    std::vector<double> x_prev, Z_prev, q_prev, rx_prev, RZ_prev, rq_prev;

    /// Algebraic variables for the discrete time integration
    DM Z, RZ;

    // Tape
    std::vector<std::vector<double> > x_tape, Z_tape;
  };

  class CASADI_EXPORT FixedStepIntegrator : public Integrator {
  public:

    /// Constructor
    explicit FixedStepIntegrator(const std::string& name, Oracle* dae);

    /// Destructor
    virtual ~FixedStepIntegrator();

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /// Initialize stage
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new FixedStepMemory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<FixedStepMemory*>(mem);}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    /// Setup F and G
    virtual void setupFG() = 0;

    /** \brief Reset the forward problem */
    virtual void reset(IntegratorMemory* mem, double t,
                       const double* x, const double* z, const double* p) const;

    /** \brief  Advance solution in time */
    virtual void advance(IntegratorMemory* mem, double t,
                         double* x, double* z, double* q) const;

    /// Reset the backward problem and take time to tf
    virtual void resetB(IntegratorMemory* mem, double t,
                        const double* rx, const double* rz, const double* rp) const;

    /** \brief  Retreat solution in time */
    virtual void retreat(IntegratorMemory* mem, double t,
                         double* rx, double* rz, double* rq) const;

    /// Get explicit dynamics
    virtual const Function& getExplicit() const { return F_;}

    /// Get explicit dynamics (backward problem)
    virtual const Function& getExplicitB() const { return G_;}

    // Discrete time dynamics
    Function F_, G_;

    // Number of finite elements
    int nk_;

    // Time step size
    double h_;

    /// Number of algebraic variables for the discrete time integration
    int nZ_, nRZ_;
  };

  class CASADI_EXPORT ImplicitFixedStepIntegrator : public FixedStepIntegrator {
  public:

    /// Constructor
    explicit ImplicitFixedStepIntegrator(const std::string& name, Oracle* dae);

    /// Destructor
    virtual ~ImplicitFixedStepIntegrator();

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /// Initialize stage
    virtual void init(const Dict& opts);

    /// Get explicit dynamics
    virtual const Function& getExplicit() const { return rootfinder_;}

    /// Get explicit dynamics (backward problem)
    virtual const Function& getExplicitB() const { return backward_rootfinder_;}

    // Implicit function solver
    Function rootfinder_, backward_rootfinder_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_INTEGRATOR_IMPL_HPP
