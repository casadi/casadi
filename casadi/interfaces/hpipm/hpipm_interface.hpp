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


#ifndef CASADI_HPIPM_INTERFACE_HPP
#define CASADI_HPIPM_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include "casadi/core/linsol.hpp"
#include <casadi/interfaces/hpipm/casadi_conic_hpipm_export.h>

/** \defgroup plugin_Conic_hpipm
Interface to HMPC Solver


In order to use this interface, you must:

 - Decision variables must only by state and control,
   and the variable ordering must be [x0 u0 x1 u1 ...]
 - The constraints must be in order: [ gap0 lincon0 gap1 lincon1  ]

    gap: Ak+1 = Ak xk + Bk uk
    lincon: yk= Ck xk + Dk uk

    \verbatim
       A0 B0 -I
       C0 D0
              A1 B1 -I
              C1 D1
    \endverbatim

   where I must be a diagonal sparse matrix
 - Either supply all of N, nx, ng, nu options or rely on automatic detection


*/

#include <blasfeo_d_aux_ext_dep.h>

#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_utils.h>

/** \pluginsection{Conic,hpipm} */

/// \cond INTERNAL
namespace casadi {

  // Forward declaration
  class HpipmInterface;

  struct CASADI_CONIC_HPIPM_EXPORT HpipmMemory : public ConicMemory {

    std::vector<double> A, B, b, b2, Q, S, R, q, r, lb, ub, C, D, lg, ug;

    std::vector<double> I;
    std::vector<double> x, u, pi, lam, lamg;
    std::vector<double*> xs, us, pis, lams, lamgs;

    std::vector<int> iidxb, iidxs;
    std::vector<int*> hidxb, hidxs;

    std::vector<double*> hA, hB, hb, hQ, hR, hS, hq, hr;
    std::vector<double*> hlb, hub, hC, hD, hlg, hug, hZl, hZu, hzl, hzu, hlls, hlus;
    std::vector<double*> hu_guess, hx_guess;
    std::vector<double*> hI;

    // Number of states in each interval
    std::vector<int> nx;
    // Number of controls in each interval
    std::vector<int> nu;
    // Number of non-dynamic constraints in each interval
    std::vector<int> ng;
    // Number of bounds on variables
    std::vector<int> nbx;
    // Number of bounds on controls
    std::vector<int> nbu;
    // Number of soft constraints
    std::vector<int> ns;

    std::vector<double> stats;

    std::vector<char> workspace;

    std::vector<double> pv;

    int iter_count;
    int return_status;
    double res_stat;
	  double res_eq;
	  double res_ineq;
	  double res_comp;

    /// Constructor
    HpipmMemory();

    /// Destructor
    ~HpipmMemory();
  };

  /** \brief \pluginbrief{Conic,hpipm}
   *
   * @copydoc QPSolver_doc
   * @copydoc plugin_Conic_hpipm
   *
   * \author Joris Gillis
   * \date 2016
   *
   * */
  class CASADI_CONIC_HPIPM_EXPORT HpipmInterface : public Conic {
  public:
    /** \brief  Constructor */
    explicit HpipmInterface();

    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new HpipmInterface(name, st);
    }

    /** \brief  Create a new Solver */
    explicit HpipmInterface(const std::string& name,
                              const std::map<std::string, Sparsity>& st);

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief  Destructor */
    ~HpipmInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "hpipm";}

    // Get name of the class
    std::string class_name() const override { return "HpipmInterface";}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new HpipmMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<HpipmMemory*>(mem);}

    /** \brief  Evaluate numerically */
    int solve(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const override;

    /** \brief Helper function */
    static void mproject(double factor, const double* x, const casadi_int* sp_x,
                         double* y, const casadi_int* sp_y, double* w);

    /** Dense transfer: y(y_sp).nonzeros() <- x(x_sp).nonzeros()
     (length >= max(number of rows, nnz)) */
    static void dense_transfer(double factor, const double* x, const casadi_int* sp_x, double* y,
                               const casadi_int* sp_y, double* w);

    /// A documentation string
    static const std::string meta_doc;

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new HpipmInterface(s); }

  protected:
    explicit HpipmInterface(DeserializingStream& s);
    
    struct Block {
      casadi_int offset_r;
      casadi_int offset_c;
      casadi_int rows;
      casadi_int cols;
    };

    static Sparsity blocksparsity(casadi_int rows, casadi_int cols,
                                   const std::vector<Block>& blocks, bool eye=false);
    static void blockptr(std::vector<double *>& vs, std::vector<double>& v,
      const std::vector<Block>& blocks, bool eye=false);
    Sparsity Asp_, Bsp_, Csp_, Dsp_, Isp_, Rsp_, Ssp_, Qsp_, bsp_, lugsp_, usp_, xsp_;
    Sparsity theirs_xsp_, theirs_usp_, theirs_Xsp_, theirs_Usp_;
    Sparsity lamg_gapsp_, lamg_csp_, lam_ulsp_, lam_uusp_, lam_xlsp_, lam_xusp_, lam_clsp_;
    Sparsity lam_cusp_, pisp_;

    std::vector< Block > R_blocks, S_blocks, Q_blocks;
    std::vector< Block > b_blocks, lug_blocks;
    std::vector< Block > u_blocks, x_blocks;
    std::vector< Block > lam_ul_blocks, lam_xl_blocks, lam_uu_blocks, lam_xu_blocks, lam_cl_blocks;
    std::vector< Block > lam_cu_blocks, A_blocks, B_blocks, C_blocks, D_blocks, I_blocks;

    std::vector<casadi_int> nxs_;
    std::vector<casadi_int> nus_;
    std::vector<casadi_int> ngs_;
    casadi_int N_;
    casadi_int print_level_;

    bool warm_start_;
    double inf_;

    d_ocp_qp_ipm_arg hpipm_options_;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_HPIPM_INTERFACE_HPP
