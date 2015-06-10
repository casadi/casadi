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


#include "mosek_socp_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

using namespace std;
namespace casadi {

  // Stream printer for MOSEK
  static void MSKAPI printstr(void *handle, MSKCONST char str[]) {
    std::cout << str;
  }

  extern "C"
  int CASADI_SOCPSOLVER_MOSEK_EXPORT
  casadi_register_socpsolver_mosek(SocpSolverInternal::Plugin* plugin) {
    plugin->creator = MosekSocpInterface::creator;
    plugin->name = "mosek";
    plugin->doc = MosekSocpInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_SOCPSOLVER_MOSEK_EXPORT casadi_load_socpsolver_mosek() {
    SocpSolverInternal::registerPlugin(casadi_register_socpsolver_mosek);
  }

  MosekSocpInterface* MosekSocpInterface::clone() const {
    // Return a deep copy
    MosekSocpInterface* node = new MosekSocpInterface(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  MosekSocpInterface::MosekSocpInterface(const std::vector<Sparsity> &st) : SocpSolverInternal(st) {
    // Introduce temporary task
    MSKenv_t temp_env;
    MSKtask_t temp_task;
    MSK_makeenv(&temp_env, NULL);
    MSK_maketask(temp_env, 0, 0, &temp_task);

    // Some variables needed to harvest parameters
    int num_param;
    int str_len;
    char str_buffer[100] = "";
    double dou_buffer;
    int int_buffer;

    // Harvest string parameters MOSEK
    MSK_getparammax(temp_task, MSK_PAR_STR_TYPE, &num_param);
    for (int i=0;i<num_param;++i) {
      MSK_getparamname(temp_task, MSK_PAR_STR_TYPE, i, str_buffer);
      std::string par_name(str_buffer);
      MSK_getnastrparam(temp_task, str_buffer, 100, &str_len, str_buffer);
      std::string par_value(str_buffer);
      addOption(par_name, OT_STRING, par_value, "Consult MOSEK manual.");
    }

    // Harvest double parameters for MOSEK
    MSK_getparammax(temp_task, MSK_PAR_DOU_TYPE, &num_param);
    for (int i=0;i<num_param;++i) {
      MSK_getparamname(temp_task, MSK_PAR_DOU_TYPE, i, str_buffer);
      std::string par_name(str_buffer);
      MSK_getnadouparam(temp_task, str_buffer, &dou_buffer);
      addOption(par_name, OT_REAL, dou_buffer, "Consult MOSEK manual.");
    }

    // Harvest integer parameters for MOSEK, passed as string to accept enumeration options
    MSK_getparammax(temp_task, MSK_PAR_INT_TYPE, &num_param);
    for (int i=0;i<num_param;++i) {
      MSK_getparamname(temp_task, MSK_PAR_INT_TYPE, i, str_buffer);
      std::string par_name(str_buffer);
      MSK_getnaintparam(temp_task, str_buffer, &int_buffer);
      std::stringstream ss;
      ss << int_buffer;
      std::string par_value_str = ss.str();
      addOption(par_name, OT_STRING, par_value_str, "Consult MOSEK manual.");
    }

    MSK_deletetask(&temp_task);
    MSK_deleteenv(&temp_env);
  }

  MosekSocpInterface::~MosekSocpInterface() {
    if (is_init_) {
      MSK_deletetask(&mosek_task_);
      MSK_deleteenv(&mosek_env_);
    }
  }

  const char* MosekSocpInterface::terminationReason(int flag) {

  }

  const char* MosekSocpInterface::solutionType(int flag) {

  }

  void MosekSocpInterface::init() {
    // Initialize the base classes
    SocpSolverInternal::init();
    log("MosekSocpInterface::init", "Enter");

    // Create the MOSEK environment
    MSK_makeenv(&mosek_env_, NULL);
    // Create empty MOSEK task
    MSK_maketask(mosek_env_, n_, m_+N_, &mosek_task_);
    // Link MOSEK task to stream printer
    MSK_linkfunctotaskstream(mosek_task_, MSK_STREAM_LOG, NULL, printstr);
    // Pre-allocate memory for MOSEK
    MSK_putmaxnumvar(mosek_task_, m_+N_+2*nc_+2*n_);
    // Append variables and constraints
    MSK_appendcons(mosek_task_, n_);
    MSK_appendvars(mosek_task_, m_+N_);
    // Set that we desire to minimize SOCP
    MSK_putobjsense(mosek_task_, MSK_OBJECTIVE_SENSE_MAXIMIZE);
    // Set known variable bounds
    for (int i=0;i<m_+N_;++i) {
      MSK_putvarbound(mosek_task_, i, MSK_BK_FR, -MSK_INFINITY, +MSK_INFINITY);
    }

    // Pass all the options to mosek
    const std::vector<string> options = getOptionNames();
    for (std::vector<string>::const_iterator it=options.begin(); it!=options.end(); ++it)
      if (hasSetOption(*it) && ((*it).find("MSK_")!=std::string::npos)) {
        GenericType op = getOption(*it);
        std::stringstream ss;
        ss << op;
        std::string optionAsString = ss.str();
        MSK_putparam(mosek_task_, (*it).c_str(), optionAsString.c_str());
      }

    // Set cone constraints (fixed)
    if (m_ > 0) {
      int largest_cone = *std::max_element(ni_.begin(), ni_.end());
      std::vector<int> submem;
      submem.resize(largest_cone);
      int sum_ni = 0;
      for (int i=0;i<m_;++i) {
        submem[0] = i;
        for (int ii=0;ii<ni_[i];++ii) submem[ii+1] = m_ + sum_ni + ii;
        MSK_appendcone(mosek_task_, MSK_CT_QUAD, 0.0, ni_[i]+1, &submem[0]);
        sum_ni += ni_[i];
      }
    }

    // Pre-allocate memory for vector indices of relevant bounds
    primal_idx_lba_.reserve(nc_);
    primal_idx_uba_.reserve(nc_);
    primal_idx_lbx_.reserve(n_);
    primal_idx_ubx_.reserve(n_);

    // Start to construct base of dual_c_
    std::vector<DMatrix> dual_c_v;
    dual_c_v.push_back(input(SOCP_SOLVER_F));
    dual_c_v.push_back(input(SOCP_SOLVER_H));
    dual_c_ = vertcat(dual_c_v).data();
    // Reserve maximum size of dual_c_
    int sizeof_c = m_+N_+2*nc_+2*n_;
    dual_c_.reserve(sizeof_c);

    // Start to construct base of dual_A_
    std::vector<DMatrix> dual_A_v;
    DMatrix dual_A_DMatrix;
    dual_A_v.push_back(input(SOCP_SOLVER_E));
    dual_A_v.push_back(input(SOCP_SOLVER_G));
    dual_A_DMatrix = horzcat(dual_A_v);
    dual_A_data_ = dual_A_DMatrix.data();
    dual_A_row_ = dual_A_DMatrix.sparsity().getRow();
    dual_A_colind_ = dual_A_DMatrix.sparsity().getColind();
    // Reserve memory for maximum size of dual_A_
    int sizeof_A = input(SOCP_SOLVER_E).size() + input(SOCP_SOLVER_G).size() +
                     2*input(SOCP_SOLVER_A).size()+2*n_;
    dual_A_data_.reserve(sizeof_A);
    dual_A_row_.reserve(sizeof_A);
    dual_A_colind_.reserve(m_+N_+2*nc_+2*n_+1);

    // Start to construct base of dual_b_
    dual_b_ = input(SOCP_SOLVER_C).data();
  }

  void MosekSocpInterface::evaluate() {
    if (inputs_check_) checkInputs();
    if (print_problem_) printProblem();

    // MOSEK expects the SOCP in conic form, obtain this form by formulating dual SOCP
    convertToDualSocp();

    // Remove or append variables
    int numvar_old;
    MSK_getnumvar(mosek_task_, &numvar_old);
    int numvar_new = m_ + N_ + primal_idx_lba_.size() + primal_idx_uba_.size() +
                       primal_idx_lbx_.size() + primal_idx_ubx_.size();
    int num_vars_to_remove = numvar_old - numvar_new;
    if (num_vars_to_remove < 0) {
      MSK_appendvars(mosek_task_, -num_vars_to_remove);
    } else if (num_vars_to_remove > 0) {
      std::vector<int> vars_to_remove;
      vars_to_remove.resize(num_vars_to_remove);
      for (int i=0;i<num_vars_to_remove;++i) {
        vars_to_remove[i] = numvar_new + num_vars_to_remove - i - 1;
      }
      MSK_removevars(mosek_task_, num_vars_to_remove, &vars_to_remove[0]);
    }

    // Add objective function
    std::vector<int> subj;
    subj.resize(numvar_new);
    double* c_val = &dual_c_[0];
    for (int i=0;i<numvar_new;++i) subj[i] = i;
    MSK_putclist(mosek_task_, numvar_new, &subj[0], c_val);
    for (int i=m_+N_;i<numvar_new;++i) {
      MSK_putvarbound(mosek_task_, i, MSK_BK_LO, 0.0, +MSK_INFINITY);
    }

    // Add equality constraints
    int* ptrb = &dual_A_colind_[0];
    int* ptre = &dual_A_colind_[1];
    int* asub = &dual_A_row_[0];
    double* aval = &dual_A_data_[0];
    MSK_putacollist(mosek_task_, numvar_new, &subj[0], ptrb, ptre, asub, aval);
    for (int i=0;i<n_;++i) {
      MSK_putconbound(mosek_task_, i, MSK_BK_FX, -dual_b_[i], -dual_b_[i]);
    }

    // Solve SOCP
    MSKrescodee r;
    MSKrescodee trmcode;

    r = MSK_optimizetrm(mosek_task_, &trmcode);

    casadi_assert_message(r==MSK_RES_OK, "MosekSocpInterface failed");

    MSK_solutionsummary(mosek_task_, MSK_STREAM_MSG);

    // Extract solution from MOSEK
    double primal_objective;
    double dual_objective;
    double* primal_solution = &output(SOCP_SOLVER_X).data()[0];
    std::vector<double> dual_solution;
    dual_solution.resize(numvar_new);
    MSK_getdualobj(mosek_task_, MSK_SOL_ITR, &primal_objective);
    MSK_getprimalobj(mosek_task_, MSK_SOL_ITR, &dual_objective);
    MSK_gety(mosek_task_, MSK_SOL_ITR, primal_solution);
    MSK_getxx(mosek_task_, MSK_SOL_ITR, &dual_solution[0]);
    output(SOCP_SOLVER_COST).set(primal_objective);
    output(SOCP_SOLVER_DUAL_COST).set(dual_objective);
    // Change sign of primal solution
    std::for_each(output(SOCP_SOLVER_X).data().begin(), output(SOCP_SOLVER_X).data().end(),
                    [](double& in){in=-in;});
    // Interpret dual solution
    std::copy(dual_solution.begin(), dual_solution.begin()+m_,
                output(SOCP_SOLVER_LAM_CONE).data().begin());
    // Fill lagrange multipliers on A*x and x with zeros
    std::fill(output(SOCP_SOLVER_LAM_A).begin(), output(SOCP_SOLVER_LAM_A).end(), 0);
    std::fill(output(SOCP_SOLVER_LAM_X).begin(), output(SOCP_SOLVER_LAM_X).end(), 0);
    // Cycle through solution vector to attain Lagrange multipliers
    int idx = m_+N_;
    for (int i : primal_idx_lba_) {
      output(SOCP_SOLVER_LAM_A).data()[i] = -dual_solution[idx];
      idx += 1;
    }
    for (int i : primal_idx_uba_) {
      if (std::abs(output(SOCP_SOLVER_LAM_A).data()[i]) < dual_solution[idx])
            output(SOCP_SOLVER_LAM_A).data()[i] = dual_solution[idx];
      idx += 1;
    }
    for (int i : primal_idx_lbx_) {
      output(SOCP_SOLVER_LAM_X).data()[i] = -dual_solution[idx];
      idx += 1;
    }
    for (int i : primal_idx_ubx_) {
      if (std::abs(output(SOCP_SOLVER_LAM_X).data()[i]) < dual_solution[idx])
            output(SOCP_SOLVER_LAM_X).data()[i] = dual_solution[idx];
      idx += 1;
    }

    // Solution status
    MSKsolstae solsta;
    MSK_getsolsta(mosek_task_, MSK_SOL_ITR, &solsta);
    stats_["solution_status"] = solutionStatus(solsta);

    // Problem status
    MSKprostae prosta;
    MSK_getprosta(mosek_task_, MSK_SOL_ITR, &prosta);
    stats_["problem_status"] = problemStatus(prosta);

    // Termination reason
    stats_["termination_reason"] = trmcode;

  }

  std::string MosekSocpInterface::solutionStatus(MSKsolstae& solsta) {
    std::string solution_status;
    switch (solsta) {
      case MSK_SOL_STA_OPTIMAL:
        solution_status = "MSK_SOL_STA_OPTIMAL";
        break;
      case MSK_SOL_STA_NEAR_OPTIMAL:
        solution_status = "MSK_SOL_STA_NEAR_OPTIMAL";
        break;
      case MSK_SOL_STA_DUAL_INFEAS_CER:
        solution_status = "MSK_SOL_STA_DUAL_INFEAS_CER";
        break;
      case MSK_SOL_STA_PRIM_INFEAS_CER:
        solution_status = "MSK_SOL_STA_PRIM_INFEAS_CER";
        break;
      case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
        solution_status = "MSK_SOL_STA_NEAR_DUAL_INFEAS_CER";
        break;
      case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
        solution_status = "MSK_SOL_STA_NEAR_PRIM_INFEAS_CER";
        break;
      case MSK_SOL_STA_UNKNOWN:
        solution_status = "MSK_SOL_STA_UNKNOWN";
        break;
      default:
        solution_status = "Unknown MOSEK solution status, contact CasADi developers.";
    }
    return solution_status;
  }

  std::string MosekSocpInterface::problemStatus(MSKprostae& prosta) {
    std::string problem_status;
    switch (prosta) {
      case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:
        problem_status = "MSK_PRO_STA_PRIM_AND_DUAL_FEAS";
        break;
      case MSK_PRO_STA_PRIM_FEAS:
        problem_status = "MSK_PRO_STA_PRIM_FEAS";
        break;
      case MSK_PRO_STA_DUAL_FEAS:
        problem_status = "MSK_PRO_STA_DUAL_FEAS";
        break;
      case MSK_PRO_STA_PRIM_INFEAS:
        problem_status = "MSK_PRO_STA_PRIM_INFEAS";
        break;
      case MSK_PRO_STA_DUAL_INFEAS:
        problem_status = "MSK_PRO_STA_DUAL_INFEAS";
        break;
      case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:
        problem_status = "MSK_PRO_STA_PRIM_AND_DUAL_INFEAS";
        break;
      case MSK_PRO_STA_ILL_POSED:
        problem_status = "MSK_PRO_STA_ILL_POSED";
        break;
      case MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS:
        problem_status = "MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS";
        break;
      case MSK_PRO_STA_NEAR_PRIM_FEAS:
        problem_status = "MSK_PRO_STA_NEAR_PRIM_FEAS";
        break;
      case MSK_PRO_STA_NEAR_DUAL_FEAS:
        problem_status = "MSK_PRO_STA_NEAR_DUAL_FEAS";
        break;
      case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
        problem_status = "MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED";
        break;
      case MSK_PRO_STA_UNKNOWN:
        problem_status = "MSK_PRO_STA_UNKNOWN";
        break;
      default:
        problem_status = "Unknown MOSEK problem status, contact CasADi developers.";
    }
    return problem_status;
  }

} // namespace casadi
