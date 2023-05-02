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


#ifndef CASADI_SUNDIALS_INTERFACE_HPP
#define CASADI_SUNDIALS_INTERFACE_HPP

#include <casadi/interfaces/sundials/casadi_sundials_common_export.h>
#include "casadi/core/integrator_impl.hpp"

#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_types.h>

#include <ctime>

/// \cond INTERNAL
namespace casadi {

  // IdasMemory
  struct CASADI_SUNDIALS_COMMON_EXPORT SundialsMemory : public IntegratorMemory {
    // N-vectors for the forward integration
    N_Vector xz, xzdot, q;

    // N-vectors for the backward integration
    N_Vector rxz, rxzdot, ruq;

    // Initialize or reinitialize?
    bool first_callB;

    // Parameters
    double *p, *rp;

    // Controls
    double *u;

    /// Jacobian memory blocks
    double *jac_ode_x, *jac_alg_x, *jac_ode_z, *jac_alg_z;

    // Jacobian
    double *jacF;

    /// Stats, forward integration
    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;
    long nniters, nncfails;

    /// Stats, backward integration
    long nstepsB, nfevalsB, nlinsetupsB, netfailsB;
    int qlastB, qcurB;
    double hinusedB, hlastB, hcurB, tcurB;
    long nnitersB, nncfailsB;

    /// Offsets for stats in backward integration
    long nstepsB_off, nfevalsB_off, nlinsetupsB_off, netfailsB_off;
    long nnitersB_off, nncfailsB_off;

    // Temporaries for [x;z] or [rx;rz]
    double *v1, *v2;

    /// number of checkpoints stored so far
    int ncheck;

    // Absolute tolerance
    N_Vector abstolv;

    /// Linear solver memory objects
    int mem_linsolF;

    /// Constructor
    SundialsMemory();

    /// Destructor
    ~SundialsMemory();
  };

  class CASADI_SUNDIALS_COMMON_EXPORT SundialsInterface : public Integrator {
  public:
    /** \brief  Constructor */
    SundialsInterface(const std::string& name, const Function& dae,
      double t0, const std::vector<double>& tout);

    /** \brief  Destructor */
    ~SundialsInterface() override=0;

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
      casadi_int*& iw, double*& w) const override;

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Get relative tolerance */
    double get_reltol() const override { return reltol_;}

    /** \brief Get absolute tolerance */
    double get_abstol() const override { return abstol_;}

    // DAE right-hand-side, forward problem
    int calc_daeF(SundialsMemory* m, double t, const double* x, const double* z,
      double* ode, double* alg) const;

    // DAE right-hand-side, forward problem
    int calc_daeB(SundialsMemory* m, double t, const double* x, const double* z,
      const double* rx, const double* rz, const double* rp, double* adj_x, double* adj_z) const;

    // Quadrature right-hand-side, forward problem
    int calc_quadF(SundialsMemory* m, double t, const double* x, const double* z,
      double* quad) const;

    // Quadrature right-hand-side, backward problem
    int calc_quadB(SundialsMemory* m, double t, const double* x, const double* z,
      const double* rx, const double* rz, double* adj_p, double* adj_u) const;

    // Jacobian of DAE-times-vector function, forward problem
    int calc_jtimesF(SundialsMemory* m, double t, const double* x, const double* z,
      const double* fwd_x, const double* fwd_z, double* fwd_ode, double* fwd_alg) const;

    // Jacobian of DAE right-hand-side function, forward problem
    int calc_jacF(SundialsMemory* m, double t, const double* x, const double* z,
      double* jac_ode_x, double* jac_alg_x, double* jac_ode_z, double* jac_alg_z) const;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief  Print solver statistics */
    void print_stats(IntegratorMemory* mem) const override;

    /** \brief  Reset the forward problem and bring the time back to t0 */
    void reset(IntegratorMemory* mem, const double* x,
      const double* z, const double* p) const override;

    /** \brief  Reset the backward problem and take time to tf */
    void resetB(IntegratorMemory* mem) const override;

    /** \brief Introduce an impulse into the backwards integration at the current time */
    void impulseB(IntegratorMemory* mem,
      const double* rx, const double* rz, const double* rp) const override;

    /** \brief Reset stats */
    void reset_stats(SundialsMemory* m) const;

    /** \brief Save stats offsets before reset */
    void save_offsets(SundialsMemory* m) const;

    /** \brief Add stats offsets to stats */
    void add_offsets(SundialsMemory* m) const;

    /** \brief Cast to memory object */
    static SundialsMemory* to_mem(void *mem) {
      SundialsMemory* m = static_cast<SundialsMemory*>(mem);
      casadi_assert_dev(m);
      return m;
    }

    ///@{
    /** \brief IO conventions for continuous time dynamics */
    enum JtimesFIn { JTIMESF_T, JTIMESF_X, JTIMESF_Z, JTIMESF_P, JTIMESF_U, JTIMESF_FWD_X,
      JTIMESF_FWD_Z, JTIMESF_NUM_IN};
    enum JtimesFOut { JTIMESF_FWD_ODE, JTIMESF_FWD_ALG, JTIMESF_NUM_OUT};
    enum JacFOut {JACF_ODE_X, JACF_ALG_X, JACF_ODE_Z, JACF_ALG_Z, JACF_NUM_OUT};
    ///@}

    ///@{
    /// Options
    double abstol_, reltol_;
    casadi_int max_num_steps_;
    bool stop_at_end_;
    bool quad_err_con_;
    casadi_int steps_per_checkpoint_;
    bool disable_internal_warnings_;
    casadi_int max_multistep_order_;
    std::string linear_solver_;
    Dict linear_solver_options_;
    casadi_int max_krylov_;
    bool use_precon_;
    bool second_order_correction_;
    double step0_;
    double max_step_size_;
    double nonlin_conv_coeff_;
    casadi_int max_order_;
    bool scale_abstol_;
    ///@}

    /// Linear solver
    Linsol linsolF_;

    /// Supported iterative solvers in Sundials
    enum NewtonScheme {SD_DIRECT, SD_GMRES, SD_BCGSTAB, SD_TFQMR} newton_scheme_;

    // Supported interpolations in Sundials
    enum InterpType {SD_POLYNOMIAL, SD_HERMITE} interp_;

    /// Linear solver data (dense) -- what is this?
    struct LinSolDataDense {};

    // Print a variable
    static void printvar(const std::string& id, double v) {
      uout() << id << " = " << v << std::endl;
    }
    // Print an N_Vector
    static void printvar(const std::string& id, N_Vector v) {
      std::vector<double> tmp(NV_DATA_S(v), NV_DATA_S(v)+NV_LENGTH_S(v));
      uout() << id << " = " << tmp << std::endl;
    }

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

  protected:
    /** \brief Deserializing constructor */
    explicit SundialsInterface(DeserializingStream& s);
  };

  // Check if N_Vector is regular
  inline bool is_regular(N_Vector v) {
    std::vector<double> tmp(NV_DATA_S(v), NV_DATA_S(v)+NV_LENGTH_S(v));
    return is_regular(tmp);
  }

} // namespace casadi

/// \endcond
#endif // CASADI_SUNDIALS_INTERFACE_HPP
