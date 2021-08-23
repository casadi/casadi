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


#ifndef CASADI_CVODES_SIMULATOR_HPP
#define CASADI_CVODES_SIMULATOR_HPP

#include <casadi/interfaces/sundials/casadi_simulator_cvodes_export.h>
#include "sundials_simulator.hpp"
#include <cvodes/cvodes.h>            /* prototypes for CVode fcts. and consts. */
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_impl.h> /* Needed for the provided linear solver */
#include <ctime>

/** \defgroup plugin_Simulator_cvodes

      Interface to CVodes from the Sundials suite.

      A call to evaluate will integrate to the end.

      You can retrieve the entire state trajectory as follows, after the evaluate call:
      Call reset. Then call integrate(t_i) and getOuput for a series of times t_i.
*/

/** \pluginsection{Simulator,cvodes} */

/// \cond INTERNAL
namespace casadi {
  // Forward declaration
  class CvodesSimulator;

  // CvodesSimMemory
  struct CASADI_SIMULATOR_CVODES_EXPORT CvodesSimMemory : public SundialsSimMemory {
    /// Function object
    const CvodesSimulator& self;

    // CVodes memory block
    void* mem;

    // Ids of backward problem
    int whichB;

    // Remember the gamma and gammaB from last factorization
    double gamma, gammaB;

    /// Constructor
    CvodesSimMemory(const CvodesSimulator& s);

    /// Destructor
    ~CvodesSimMemory();
  };

  /** \brief \pluginbrief{Simulator,cvodes}

      @copydoc DAE_doc
      @copydoc plugin_Simulator_cvodes

  */
  class CASADI_SIMULATOR_CVODES_EXPORT CvodesSimulator : public SundialsSimulator {
  public:
    /** \brief  Constructor */
    explicit CvodesSimulator(const std::string& name, const Function& dae,
      const std::vector<double>& grid);

    /** \brief  Create a new simulator */
    static Simulator* creator(const std::string& name, const Function& dae,
        const std::vector<double>& grid) {
      return new CvodesSimulator(name, dae, grid);
    }

    /** \brief  Destructor */
    ~CvodesSimulator() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "cvodes";}

    // Get name of the class
    std::string class_name() const override { return "CvodesSimulator";}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    double min_step_size_;
    ///@}


    /** \brief  Initialize stage */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new CvodesSimMemory(*this);}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<CvodesSimMemory*>(mem);}

    /** \brief  Reset the forward problem and bring the time back to t0 */
    void reset(SimulatorMemory* mem, double t, const double* x,
                       const double* z, const double* p) const override;

    /** \brief  Advance solution in time */
    void advance(SimulatorMemory* mem, double t, double* x,
                         double* z, double* q) const override;

    /** \brief  Reset the backward problem and take time to tf */
    void resetB(SimulatorMemory* mem, double t,
                        const double* rx, const double* rz, const double* rp) const override;

    /** \brief  Retreat solution in time */
    void retreat(SimulatorMemory* mem, double t, double* rx,
                         double* rz, double* rq) const override;

    /** \brief  Set the stop time of the forward integration */
    void setStopTime(SimulatorMemory* mem, double tf) const override;

    /** \brief Cast to memory object */
    static CvodesSimMemory* to_mem(void *mem) {
      CvodesSimMemory* m = static_cast<CvodesSimMemory*>(mem);
      casadi_assert_dev(m);
      return m;
    }

    ///@{
    // Get system Jacobian
    Function getJ(bool backward) const override;
    template<typename MatType> Function getJ(bool backward) const;
    ///@}

    /// A documentation string
    static const std::string meta_doc;

  protected:

    // Sundials callback functions
    static int rhs(double t, N_Vector x, N_Vector xdot, void *user_data);
    static void ehfun(int error_code, const char *module, const char *function, char *msg,
                      void *user_data);
    static int rhsQ(double t, N_Vector x, N_Vector qdot, void *user_data);
    static int rhsB(double t, N_Vector x, N_Vector xB, N_Vector xdotB, void *user_data);
    static int rhsQB(double t, N_Vector x, N_Vector xB, N_Vector qdotB, void *user_data);
    static int jtimes(N_Vector v, N_Vector Jv, double t, N_Vector x, N_Vector xdot,
                      void *user_data, N_Vector tmp);
    static int jtimesB(N_Vector vB, N_Vector JvB, double t, N_Vector x, N_Vector xB,
                       N_Vector xdotB, void *user_data , N_Vector tmpB);
    static int psolve(double t, N_Vector x, N_Vector xdot, N_Vector r, N_Vector z,
                      double gamma, double delta, int lr, void *user_data, N_Vector tmp);
    static int psolveB(double t, N_Vector x, N_Vector xB, N_Vector xdotB, N_Vector rvecB,
                       N_Vector zvecB, double gammaB, double deltaB,
                       int lr, void *user_data, N_Vector tmpB);
    static int psetup(double t, N_Vector x, N_Vector xdot, booleantype jok,
                      booleantype *jcurPtr, double gamma, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int psetupB(double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                       booleantype jokB, booleantype *jcurPtrB, double gammaB,
                       void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    static int lsetup(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                      booleantype *jcurPtr,
                      N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    static int lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector x,
                      N_Vector xdot);
    static int lsetupB(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                       booleantype *jcurPtr,
                       N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    static int lsolveB(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                       N_Vector x, N_Vector xdot);

    // Throw error
    static void cvodes_error(const char* module, int flag);

    casadi_int lmm_; // linear multistep method
    casadi_int iter_; // nonlinear solver iteration
  };

} // namespace casadi

/// \endcond
#endif // CASADI_CVODES_SIMULATOR_HPP
