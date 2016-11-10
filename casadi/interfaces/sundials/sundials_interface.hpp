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


#ifndef CASADI_SUNDIALS_INTERFACE_HPP
#define CASADI_SUNDIALS_INTERFACE_HPP

#include <casadi/interfaces/sundials/casadi_sundials_common_export.h>
#include "casadi/core/function/integrator_impl.hpp"

#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_types.h>

#include <ctime>

/// \cond INTERNAL
namespace casadi {

  // IdasMemory
  struct CASADI_SUNDIALS_COMMON_EXPORT SundialsMemory : public IntegratorMemory {
    // Current time
    double t;

    // N-vectors for the forward integration
    N_Vector xz, xzdot, q;

    // N-vectors for the backward integration
    N_Vector rxz, rxzdot, rq;

    // Initialize or reinitialize?
    bool first_callB;

    // Parameters
    double *p, *rp;

    // Jacobian
    double *jac, *jacB;

    /// Stats
    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;

    long nstepsB, nfevalsB, nlinsetupsB, netfailsB;
    int qlastB, qcurB;
    double hinusedB, hlastB, hcurB, tcurB;

    // Temporary for max(x,rx)
    double *xtmp;

    // Temporary for max(z,rz)
    double *ztmp;

    /// number of checkpoints stored so far
    int ncheck;

    /// Constructor
    SundialsMemory();

    /// Destructor
    ~SundialsMemory();
  };

  class CASADI_SUNDIALS_COMMON_EXPORT SundialsInterface : public Integrator {
  public:
    /** \brief  Constructor */
    SundialsInterface(const std::string& name, const Function& dae);

    /** \brief  Destructor */
    virtual ~SundialsInterface()=0;

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    // Get system Jacobian
    virtual Function getJ(bool backward) const = 0;

    /// Get all statistics
    virtual Dict get_stats(void* mem) const;

    /** \brief  Print solver statistics */
    virtual void print_stats(IntegratorMemory* mem, std::ostream &stream) const;

    /** \brief  Reset the forward problem and bring the time back to t0 */
    virtual void reset(IntegratorMemory* mem, double t, const double* x,
                       const double* z, const double* p) const;

    /** \brief  Reset the backward problem and take time to tf */
    virtual void resetB(IntegratorMemory* mem, double t, const double* rx,
                        const double* rz, const double* rp) const;

    /** \brief Cast to memory object */
    static SundialsMemory* to_mem(void *mem) {
      SundialsMemory* m = static_cast<SundialsMemory*>(mem);
      casadi_assert(m);
      return m;
    }

    ///@{
    /// Options
    double abstol_, reltol_;
    int max_num_steps_;
    bool stop_at_end_;
    bool quad_err_con_;
    int steps_per_checkpoint_;
    bool disable_internal_warnings_;
    int max_multistep_order_;
    std::string linear_solver_;
    Dict linear_solver_options_;
    bool use_iterative_solver_;
    ///@}

    /// Supported iterative solvers in Sundials
    enum IterativeSolverType {SD_GMRES, SD_BCGSTAB, SD_TFQMR} itsol_;

    // Supported interpolations in Sundials
    enum InterpType {SD_POLYNOMIAL, SD_HERMITE} interp_;

    /// Linear solver data (dense)
    struct LinSolDataDense {};

    /// Max Krylov size
    int max_krylov_;

    /// Use preconditioning
    bool use_precon_;

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;

    // Print a variable
    static void printvar(const std::string& id, double v) {
      userOut() << id << " = " << v << std::endl;
    }
    // Print an N_Vector
    static void printvar(const std::string& id, N_Vector v) {
      std::vector<double> tmp(NV_DATA_S(v), NV_DATA_S(v)+NV_LENGTH_S(v));
      userOut() << id << " = " << tmp << std::endl;
    }
  };

  // Check if N_Vector is regular
  inline bool is_regular(N_Vector v) {
    std::vector<double> tmp(NV_DATA_S(v), NV_DATA_S(v)+NV_LENGTH_S(v));
    return is_regular(tmp);
  }

} // namespace casadi

/// \endcond
#endif // CASADI_SUNDIALS_INTERFACE_HPP
