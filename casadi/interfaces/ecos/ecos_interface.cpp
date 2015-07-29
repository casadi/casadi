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


#include "ecos_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_SOCPSOLVER_ECOS_EXPORT
  casadi_register_socpsolver_ecos(SocpSolverInternal::Plugin* plugin) {
    plugin->creator = EcosInterface::creator;
    plugin->name = "ecos";
    plugin->doc = EcosInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_SOCPSOLVER_ECOS_EXPORT casadi_load_socpsolver_ecos() {
    SocpSolverInternal::registerPlugin(casadi_register_socpsolver_ecos);
  }

  EcosInterface* EcosInterface::clone() const {
    // Return a deep copy
    EcosInterface* node = new EcosInterface(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  EcosInterface::EcosInterface(const std::vector<Sparsity> &st) : SocpSolverInternal(st) {
    // Define options
    addOption("gamma"         , OT_REAL     , GAMMA,
                "Scaling the final step length");
    addOption("delta"         , OT_REAL     , DELTA,
                "Regularization parameter");
    addOption("eps"           , OT_REAL     , EPS,
                "Regularization threshold");
    addOption("feastol"       , OT_REAL     , FEASTOL,
                "Primal/dual infeasibility tolerance");
    addOption("abstol"        , OT_REAL     , ABSTOL,
                 "Absolute tolerance on duality gap");
    addOption("reltol"        , OT_REAL     , RELTOL,
                 "Relative tolerance on duality gap");
    addOption("feastol_inacc" , OT_REAL     , FTOL_INACC,
                "Primal/dual infeasibility relaxed tolerance");
    addOption("abstol_inacc"  , OT_REAL     , ATOL_INACC,
                "Absolute relaxed tolerance on duality gap");
    addOption("reltol_inacc"  , OT_REAL     , RTOL_INACC,
                "Relative relaxed tolerance on duality gap");
    addOption("nitref"        , OT_INTEGER  , NITREF,
                "Number of iterative refinement steps");
    addOption("maxit"         , OT_INTEGER  , MAXIT,
                "Maximum number of iterations");
    addOption("verbose"       , OT_INTEGER  , VERBOSE,
                "Verbosity bool for PRINTLEVEL < 3");
  }

  EcosInterface::~EcosInterface() {
    if (is_init_) {
      // Destructor
    }
  }

  void EcosInterface::init() {
    // Initialize the base classes
    SocpSolverInternal::init();
    log("EcosInterface::init", "Enter");

    // Initialize b
    ecos_b_vec_.resize(n_);

    // Initialize c
    ecos_c_vec_.resize(m_+N_+2*nc_+2*n_);

    // Initialize q
    for (int i=0;i<m_;++i) ecos_q_vec_.push_back(ni_[i]+1);

    // Initalize h
    for (int i=0;i<m_+N_+2*nc_+2*n_;++i) ecos_h_vec_.push_back(0.0);

    // Initialize G
    ecos_Gpr_vec_ = std::vector<pfloat>(m_+N_+2*nc_+2*n_, -1.0);
    ecos_Gjc_vec_.resize(m_+N_+2*nc_+2*n_+1);
    ecos_Gir_vec_.resize(m_+N_+2*nc_+2*n_);
    for (int i=0;i<m_+N_+2*nc_+2*n_+1;++i) ecos_Gjc_vec_[i] = i;

    // Initialize A
    int sizeof_A = input(SOCP_SOLVER_E).size() + input(SOCP_SOLVER_G).size() +
                2*input(SOCP_SOLVER_A).size()+2*n_;
    ecos_Ajc_vec_.reserve(m_+N_+2*nc_+2*n_+1);
    ecos_Air_vec_.reserve(sizeof_A);

  }

  void EcosInterface::evaluate() {
    if (inputs_check_) checkInputs();
    if (print_problem_) printProblem();

    // ECOS expects the SOCP in conic form, obtain this form by formulating dual SOCP
    convertToDualSocp();

    /** Assemble input for ECOS solver
             n (idxint)    : number of variables
             m (idxint)    : number of inequality constraints
             p (idxint)    : number of equality constraints
             l (idxint)    : dimension of positive orthant
        ncones (idxint)    : number of second order cones
             q (idxint[])  : dimension of each cone
           Gpr (pfloat[])  : matrix G in CCS format, data
           Gjc (idxint[])  : matrix G in CCS format, column indices
           Gir (idxint[])  : matrix G in CCS format, row index arrays
           Apr (pfloat[])  : matrix A in CCS format, data
           Ajc (idxint[])  : matrix A in CCS format, column indices
           Air (idxint[])  : matrix A in CCS format, row index arrays
             c (pfloat[])  : vector c (dense)
             h (pfloat[])  : vector h (dense)
             b (pfloat[])  : vector b (dense)
       
        For problem:
          min c'*x
         s.t. A*x = b
              G*x <=_K h
       
        where the last inequality is generalized, i.e. h - G*x belongs to the cone K. ECOS supports the
        positive orthant R_+ and second-order cones Q_n defined as
             Q_n = { (t,x) | t >= || x ||_2 }
       
        In the definition above, t is a scalar and x is in R_{n-1}. The cone K is therefore a direct product of
        the positive orthant and second-order cones:
             K = R_+ x Q_n1 x ... x Q_nN
    */

    // Update A-matrix
    ecos_Ajc_vec_.resize(0);
    ecos_Air_vec_.resize(0);
    for (int i=0;i<dual_A_colind_.size();++i) ecos_Ajc_vec_.push_back(dual_A_colind_[i]);
    for (int i=0;i<dual_A_row_.size();++i) ecos_Air_vec_.push_back(dual_A_row_[i]);
    std::transform(dual_c_.begin(), dual_c_.end(), ecos_c_vec_.begin(), std::negate<pfloat>());
    std::transform(dual_b_.begin(), dual_b_.end(), ecos_b_vec_.begin(), std::negate<pfloat>());
    // Update G-matrix (change row-indices)
    int sum_ni = 0;
    int dim_pos_orthant = dual_n_-m_-N_;
    for (int i=0;i<m_;++i) {
      ecos_Gir_vec_[i] = dim_pos_orthant + sum_ni + i;
      for (int ii=0;ii<ni_[i];++ii) {
        ecos_Gir_vec_[m_+sum_ni+ii] = dim_pos_orthant + sum_ni + i + ii + 1;
      }
      sum_ni += ni_[i];
    }
    for (int i=0;i<dim_pos_orthant;++i) ecos_Gir_vec_[m_+N_+i] = i;

    idxint n        = dual_n_;
    idxint m        = dual_n_;
    idxint p        = dual_nc_;
    idxint l        = dim_pos_orthant;
    idxint ncones   = m_;
    idxint* q       = &ecos_q_vec_[0];
    idxint nex      = 0;
    pfloat* Gpr     = &ecos_Gpr_vec_[0];
    idxint* Gjc     = &ecos_Gjc_vec_[0];
    idxint* Gir     = &ecos_Gir_vec_[0];
    pfloat* Apr     = &dual_A_data_[0];
    idxint* Ajc     = &ecos_Ajc_vec_[0];
    idxint* Air     = &ecos_Air_vec_[0];
    pfloat* c       = &ecos_c_vec_[0];
    pfloat* h       = &ecos_h_vec_[0];
    pfloat* b       = &ecos_b_vec_[0];


    // Setup ECOS for new problem based on problem definition. We should be able to place
    // this in init(), it requires ECOS to do a re-initialize.
#ifdef EXPCONE
    pwork* ecos_work =  ECOS_setup(n, m, p, l, ncones, q, nex, Gpr, Gjc, Gir, 
                Apr, Ajc, Air, c, h, b);
#else // EXPCONE
    pwork* ecos_work =  ECOS_setup(n, m, p, l, ncones, q, Gpr, Gjc, Gir, 
                Apr, Ajc, Air, c, h, b);
#endif // EXPCONE

    // Pass all options to ECOS
    int temp_nitref = getOption("nitref");
    int temp_maxit = getOption("maxit");
    int temp_verbose = getOption("verbose");
    ecos_work->stgs->gamma = getOption("gamma");
    ecos_work->stgs->delta = getOption("delta");
    ecos_work->stgs->eps = getOption("eps");
    ecos_work->stgs->feastol = getOption("feastol");
    ecos_work->stgs->abstol = getOption("abstol");
    ecos_work->stgs->reltol = getOption("reltol");
    ecos_work->stgs->feastol_inacc = getOption("feastol_inacc");
    ecos_work->stgs->abstol_inacc = getOption("abstol_inacc");
    ecos_work->stgs->reltol_inacc = getOption("reltol_inacc");
    ecos_work->stgs->nitref = temp_nitref;
    ecos_work->stgs->maxit = temp_maxit;
    ecos_work->stgs->verbose = temp_verbose;


    // Solve problem with ECOS
    idxint ret = ECOS_solve(ecos_work);

    // Extract solution from ECOS
    double* dual_sol = ecos_work->x;
    double* primal_sol = ecos_work->y;
    double* output_primal = &output(SOCP_SOLVER_X).data()[0];
    double* output_lam_cone = &output(SOCP_SOLVER_LAM_CONE).data()[0];
    std::transform(primal_sol, primal_sol + n_, output_primal, std::negate<pfloat>());
    // Primal and dual cost
    output(SOCP_SOLVER_COST).set(-ecos_work->info->dcost);
    output(SOCP_SOLVER_DUAL_COST).set(-ecos_work->info->pcost);
    // Fill lagrange multipliers on A*x and x with zeros
    std::fill(output(SOCP_SOLVER_LAM_A).begin(), output(SOCP_SOLVER_LAM_A).end(), 0);
    std::fill(output(SOCP_SOLVER_LAM_X).begin(), output(SOCP_SOLVER_LAM_X).end(), 0);
    // Cycle through solution vector to attain Lagrange multipliers
    int idx = m_+N_;
    int k;
    for (int i=0;i<primal_idx_lba_.size();++i) {
      k = primal_idx_lba_[i];
      output(SOCP_SOLVER_LAM_A).data()[k] = -dual_sol[idx];
      idx += 1;
    }
    for (int i=0;i<primal_idx_uba_.size();++i) {
      k = primal_idx_uba_[i];
      if (std::abs(output(SOCP_SOLVER_LAM_A).data()[k]) < dual_sol[idx])
                output(SOCP_SOLVER_LAM_A).data()[k] += dual_sol[idx];
      idx += 1;
    }
    for (int i=0;i<primal_idx_lbx_.size();++i) {
      k = primal_idx_lbx_[i];
      output(SOCP_SOLVER_LAM_X).data()[k] = -dual_sol[idx];
      idx += 1;
    }
    for (int i=0;i<primal_idx_ubx_.size();++i) {
      k = primal_idx_ubx_[i];
      if (std::abs(output(SOCP_SOLVER_LAM_X).data()[k]) < dual_sol[idx])
                output(SOCP_SOLVER_LAM_X).data()[k] += dual_sol[idx];
      idx += 1;
    }

    // Write solution statistics
    stats_["exit_code"] = exitCode(ret);
    stats_["exit_code_explanation"] = exitCodeLong(ret);
    stats_["pcost"] = ecos_work->info->pcost;
    stats_["dcost"] = ecos_work->info->dcost;
    stats_["pres"] = ecos_work->info->pres;
    stats_["dres"] = ecos_work->info->dres;
    stats_["pinf"] = ecos_work->info->pinf;
    stats_["dinf"] = ecos_work->info->dinf;
    stats_["pinfres"] = ecos_work->info->pinfres;
    stats_["dinfres"] = ecos_work->info->dinfres;
    stats_["gap"] = ecos_work->info->gap;
    stats_["relgap"] = ecos_work->info->relgap;
    stats_["sigma"] = ecos_work->info->sigma;
    stats_["mu"] = ecos_work->info->mu;
    stats_["step"] = ecos_work->info->step;
    stats_["step_aff"] = ecos_work->info->step_aff;
    stats_["kapovert"] = ecos_work->info->kapovert;
    stats_["iter"] = static_cast<int>(ecos_work->info->iter);
    stats_["nitref1"] = static_cast<int>(ecos_work->info->nitref1);
    stats_["nitref2"] = static_cast<int>(ecos_work->info->nitref2);
    stats_["nitref3"] = static_cast<int>(ecos_work->info->nitref3);

    // Cleanup memory
    ECOS_cleanup(ecos_work, 0);

  }

  const std::string EcosInterface::exitCode(const idxint& retval) const {
    std::string exit_code;
    switch (retval) {
      case ECOS_OPTIMAL:
        exit_code = "ECOS_OPTIMAL";
        break;
      case ECOS_PINF:
        exit_code = "ECOS_PINF";
        break;
      case ECOS_DINF:
        exit_code = "ECOS_DINF";
        break;
      case ECOS_INACC_OFFSET:
        exit_code = "ECOS_INACC_OFFSET";
        break;
      case ECOS_MAXIT:
        exit_code = "ECOS_MAXIT";
        break;
      case ECOS_NUMERICS:
        exit_code = "ECOS_NUMERICS";
        break;
      case ECOS_OUTCONE:
        exit_code = "ECOS_OUTCONE";
        break;
      case ECOS_SIGINT:
        exit_code = "ECOS_SIGINT";
        break;
      case ECOS_FATAL:
        exit_code = "ECOS_FATAL";
        break;
      default:
        exit_code = "Unknown ECOS exit code, contact CasADi developers.";
    }
    return exit_code;
  }

  const std::string EcosInterface::exitCodeLong(const idxint& retval) const {
    std::string exit_code;
    switch (retval) {
      case ECOS_OPTIMAL:
        exit_code = "Problem solved to optimality";
        break;
      case ECOS_PINF:
        exit_code = "Found certificate of primal infeasibility";
        break;
      case ECOS_DINF:
        exit_code = "Found certificate of dual infeasibility";
        break;
      case ECOS_INACC_OFFSET:
        exit_code = "Offset exitflag at inaccurate results";
        break;
      case ECOS_MAXIT:
        exit_code = "Maximum number of iterations reached";
        break;
      case ECOS_NUMERICS:
        exit_code = "Search direction unreliable";
        break;
      case ECOS_OUTCONE:
        exit_code = "s or z got outside the cone, numerics?";
        break;
      case ECOS_SIGINT:
        exit_code = "Solver interrupted by a signal/ctrl-c ";
        break;
      case ECOS_FATAL:
        exit_code = "Unknown problem in solver";
        break;
      default:
        exit_code = "Unknown ECOS exit code, contact CasADi developers.";
    }
    return exit_code;
  }


} // namespace casadi
