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


#ifndef CASADI_IDAS_INTERFACE_HPP
#define CASADI_IDAS_INTERFACE_HPP

#include <casadi/interfaces/sundials/casadi_integrator_idas_export.h>
#include "sundials_interface.hpp"
#include <idas/idas.h>            /* prototypes for CVODE fcts. and consts. */
#include <idas/idas_dense.h>
#include <idas/idas_band.h>
#include <idas/idas_spgmr.h>
#include <idas/idas_spbcgs.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_impl.h> /* Needed for the provided linear solver */
#include <ctime>

/** \defgroup plugin_Integrator_idas Title
    \par

      Interface to IDAS from the Sundials suite.

    \identifier{225} */

/** \pluginsection{Integrator,idas} */

/// \cond INTERNAL

namespace casadi {

// Forward declaration
class IdasInterface;

// IdasMemory
struct CASADI_INTEGRATOR_IDAS_EXPORT IdasMemory : public SundialsMemory {
  /// Function object
  const IdasInterface& self;

  /// Idas memory block
  void* mem;

  /// Ids of backward problem
  int whichB;

  /// cj used in the last factorization
  double cj_last;

  /// Constructor
  IdasMemory(const IdasInterface& s);

  /// Destructor
  ~IdasMemory();
};

/** \brief \pluginbrief{Integrator,idas}

    @copydoc IdasIntegrator_doc
    @copydoc plugin_Integrator_idas

    \author Joel Andersson
    \date 2010
*/
class CASADI_INTEGRATOR_IDAS_EXPORT IdasInterface : public SundialsInterface {
 public:

  /** \brief  Constructor */
  IdasInterface(const std::string& name, const Function& dae,
    double t0, const std::vector<double>& tout);

  /** \brief  Create a new integrator */
  static Integrator* creator(const std::string& name, const Function& dae,
      double t0, const std::vector<double>& tout) {
    return new IdasInterface(name, dae, t0, tout);
  }

  /** \brief  Destructor */
  ~IdasInterface() override;

  // Get name of the plugin
  const char* plugin_name() const override { return "idas";}

  // Get name of the class
  std::string class_name() const override { return "IdasInterface";}

  ///@{
  /** \brief Options */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /** \brief  Initialize */
  void init(const Dict& opts) override;

  /** \brief Create memory block */
  void* alloc_mem() const override { return new IdasMemory(*this);}

  /** \brief Initalize memory block */
  int init_mem(void* mem) const override;

  /** \brief Free memory block */
  void free_mem(void *mem) const override { delete static_cast<IdasMemory*>(mem);}

  /** \brief  Reset the forward problem and bring the time back to t0 */
  void reset(IntegratorMemory* mem,
    const double* x, const double* z, const double* p) const override;

  /** \brief  Advance solution in time */
  void advance(IntegratorMemory* mem,
    const double* u, double* x, double* z, double* q) const override;

  /** \brief  Reset the backward problem and take time to tf */
  void resetB(IntegratorMemory* mem) const override;

  /** \brief Propagate impulse from rz to rx */
  void z_impulseB(IdasMemory* m, const double* rz) const;

  /** \brief Solve transposed linear system */
  int solve_transposed(IdasMemory* m, double t, const double* xz, const double* rxz,
    const double* rhs, double* sol) const;

  /** \brief Introduce an impulse into the backwards integration at the current time */
  void impulseB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const override;

  /** \brief  Retreat solution in time */
  void retreat(IntegratorMemory* mem, const double* u,
    double* rx, double* rq, double* uq) const override;

  /** \brief Cast to memory object */
  static IdasMemory* to_mem(void *mem) {
    IdasMemory* m = static_cast<IdasMemory*>(mem);
    casadi_assert_dev(m);
    return m;
  }

  /// A documentation string
  static const std::string meta_doc;

 protected:

  // Sundials callback functions
  static int resF(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, void *user_data);
  static int resB(double t, N_Vector xz, N_Vector xzdot, N_Vector rxz, N_Vector rxzdot,
    N_Vector rr, void *user_data);
  static void ehfun(int error_code, const char *module, const char *function, char *msg,
    void *eh_data);
  static int jtimesF(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector v,
    N_Vector Jv, double cj, void *user_data, N_Vector tmp1, N_Vector tmp2);
  static int jtimesB(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
    N_Vector resvalB, N_Vector vB, N_Vector JvB, double cjB,
    void *user_data, N_Vector tmp1B, N_Vector tmp2B);
  static int rhsQF(double t, N_Vector xz, N_Vector xzdot, N_Vector qdot, void *user_data);
  static int rhsQB(double t, N_Vector xz, N_Vector xzdot, N_Vector rxz,
    N_Vector rxzdot, N_Vector ruqdot, void *user_data);
  static int psolveF(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector rvec,
    N_Vector zvec, double cj, double delta, void *user_data, N_Vector tmp);
  static int psetupF(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, double cj,
    void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int psolveB(double t, N_Vector xz, N_Vector xzdot, N_Vector rxz, N_Vector rxzdot,
    N_Vector resvalB, N_Vector rvecB, N_Vector zvecB, double cjB,
    double deltaB, void *user_dataB, N_Vector tmpB);
  static int psetupB(double t, N_Vector xz, N_Vector xzdot, N_Vector rxz, N_Vector rxzdot,
    N_Vector resvalB, double cjB, void *user_dataB,
    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
  static int lsetupF(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp,
    N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  static int lsolveF(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur,
    N_Vector xzdotcur, N_Vector rescur);
  static int lsetupB(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp,
    N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  static int lsolveB(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur,
    N_Vector xzdotcur, N_Vector rescur);
 public:

  // Throw error
  static void idas_error(const char* module, int flag);

  // Options
  bool cj_scaling_;
  bool calc_ic_;
  bool calc_icB_;
  bool suppress_algebraic_;
  std::vector<double> abstolv_;
  double first_time_;

  // Constraints on decision variables
  std::vector<casadi_int> y_c_;

  //  Initial values for \p xdot
  std::vector<double> init_xdot_;

  /** \brief Serialize an object without type information */
  void serialize_body(SerializingStream &s) const override;

  /** \brief Deserialize into MX */
  static ProtoFunction* deserialize(DeserializingStream& s) { return new IdasInterface(s); }

 protected:
  /** \brief Deserializing constructor */
  explicit IdasInterface(DeserializingStream& s);
};

} // namespace casadi

/// \endcond
#endif // CASADI_IDAS_INTERFACE_HPP
