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


#include "mosek_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

using namespace std;
namespace casadi {

  // Stream printer for MOSEK
  static void MSKAPI printstr(void *handle, MSKCONST char str[]){
    std::cout << str;
  }

  extern "C"
  int CASADI_SOCPSOLVER_MOSEK_EXPORT
  casadi_register_socpsolver_mosek(SocpSolverInternal::Plugin* plugin) {
    plugin->creator = MosekInterface::creator;
    plugin->name = "mosek";
    plugin->doc = MosekInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_SOCPSOLVER_MOSEK_EXPORT casadi_load_socpsolver_mosek() {
    SocpSolverInternal::registerPlugin(casadi_register_socpsolver_mosek);
  }

  MosekInterface* MosekInterface::clone() const {
    // Return a deep copy
    MosekInterface* node = new MosekInterface(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  MosekInterface::MosekInterface(const std::vector<Sparsity> &st) : SocpSolverInternal(st) {
    // Introduce temporary task
    MSKenv_t temp_env;
    MSKtask_t temp_task;
    MSK_makeenv(&temp_env,NULL);
    MSK_maketask(temp_env,0,0,&temp_task);
    
    // Some variables needed to harvest parameters
    int num_param;
    int str_len;
    char str_buffer [100] = "";
    double dou_buffer;
    int int_buffer;

    // Harvest string parameters MOSEK
    MSK_getparammax(temp_task,MSK_PAR_STR_TYPE,&num_param);
    for (int i=0;i<num_param;++i) {
      MSK_getparamname(temp_task,MSK_PAR_STR_TYPE,i,str_buffer);
      std::string par_name(str_buffer);
      MSK_getnastrparam(temp_task,str_buffer,100,&str_len,str_buffer);
      std::string par_value(str_buffer);
      addOption(par_name, OT_STRING, par_value, "Consult MOSEK manual.");
    }

    // Harvest double parameters for MOSEK
    MSK_getparammax(temp_task,MSK_PAR_DOU_TYPE,&num_param);
    for (int i=0;i<num_param;++i) {
      MSK_getparamname(temp_task,MSK_PAR_DOU_TYPE,i,str_buffer);
      std::string par_name(str_buffer);
      MSK_getnadouparam(temp_task,str_buffer,&dou_buffer);
      addOption(par_name, OT_REAL, dou_buffer, "Consult MOSEK manual.");
    }

    // Harvest integer parameters for MOSEK, passed as string to accept enumeration options
    MSK_getparammax(temp_task,MSK_PAR_INT_TYPE,&num_param);
    for (int i=0;i<num_param;++i) {
      MSK_getparamname(temp_task,MSK_PAR_INT_TYPE,i,str_buffer);
      std::string par_name(str_buffer);
      MSK_getnaintparam(temp_task,str_buffer,&int_buffer);
      std::stringstream ss;
      ss << int_buffer;
      std::string par_value_str = ss.str();
      addOption(par_name, OT_STRING, par_value_str, "Consult MOSEK manual.");
    }

    MSK_deletetask(&temp_task);
    MSK_deleteenv(&temp_env);
  }

  MosekInterface::~MosekInterface() {
    if (is_init_) {
      MSK_deletetask(&mosek_task_);
      MSK_deleteenv(&mosek_env_);
    }
  }

  const char* MosekInterface::terminationReason(int flag) {

  }

  const char* MosekInterface::solutionType(int flag) {

  }

  void MosekInterface::init() {
    // Initialize the base classes
    SocpSolverInternal::init();
    log("MosekInterface::init", "Enter");
	   
	  // Create the MOSEK environment
	  MSK_makeenv(&mosek_env_,NULL);
    // Create empty MOSEK task
    MSK_maketask(mosek_env_,n_,m_+N_,&mosek_task_);
    // Link MOSEK task to stream printer
    MSK_linkfunctotaskstream(mosek_task_,MSK_STREAM_LOG,NULL,printstr);
    // Pre-allocate memory for MOSEK
    MSK_putmaxnumvar(mosek_task_,m_+N_+2*nc_+2*n_);
    // Append variables and constraints
    MSK_appendcons(mosek_task_,n_);
    MSK_appendvars(mosek_task_,m_+N_);
    // Set that we desire to minimize SOCP
    MSK_putobjsense(mosek_task_, MSK_OBJECTIVE_SENSE_MAXIMIZE);
    // Set known variable bounds
    for (int i=0;i<m_+N_;++i){
      MSK_putvarbound(mosek_task_,i,MSK_BK_FR,-MSK_INFINITY,+MSK_INFINITY);
    }

    // Pass all the options to mosek
    const std::vector<string> options = getOptionNames();
    for (std::vector<string>::const_iterator it=options.begin(); it!=options.end(); ++it)
      if (hasSetOption(*it) && ((*it).find("MSK_")!=std::string::npos)) {
        GenericType op = getOption(*it);
        std::stringstream ss;
        ss << op;
        std::string optionAsString = ss.str();
        MSK_putparam(mosek_task_,(*it).c_str(), optionAsString.c_str());
      }

    // Set cone constraints (fixed)
    int largest_cone = *std::max_element(ni_.begin(),ni_.end());
    int submem [largest_cone];
    int sum_ni = 0;
    for (int i=0;i<m_;++i) {
      submem [0] = i;
      for(int ii=0;ii<ni_[i];++ii) submem[ii+1] = m_ + sum_ni + ii;
      MSK_appendcone(mosek_task_,MSK_CT_QUAD,0.0,ni_[i]+1,submem);
      sum_ni += ni_[i];
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
    int sizeof_A = input(SOCP_SOLVER_E).size()+input(SOCP_SOLVER_G).size()+2*input(SOCP_SOLVER_A).size()+2*n_;
    dual_A_data_.reserve(sizeof_A);
    dual_A_row_.reserve(sizeof_A);
    dual_A_colind_.reserve(m_+N_+2*nc_+2*n_+1);

    // Start to construct base of dual_b_
    dual_b_ = input(SOCP_SOLVER_C).data();

    // Start to construct base of dual_yi_
    for (int i=m_;i<m_+N_;++i) dual_yi_.push_back(i);

    // Start to construct base of dual_ti_
    for (int i=0;i<m_;++i) dual_ti_.push_back(i);
  }

  void MosekInterface::evaluate() {
    if (inputs_check_) checkInputs();
    if (print_problem_) printProblem();

    // MOSEK expects the SOCP in conic form, obtain this form by formulating dual SOCP
    convertToDualSocp();

    // Remove or append variables
    int numvar_old;
    MSK_getnumvar(mosek_task_,&numvar_old);
    int numvar_new = m_ + N_ + primal_idx_lba_.size() + primal_idx_uba_.size() + primal_idx_lbx_.size() + primal_idx_ubx_.size();
    int num_vars_to_remove = numvar_old - numvar_new;
    if (num_vars_to_remove < 0){
      MSK_appendvars(mosek_task_,-num_vars_to_remove);
    } else if (num_vars_to_remove > 0){
      int vars_to_remove [num_vars_to_remove];
      for (int i=0;i<num_vars_to_remove;++i) vars_to_remove[i] = numvar_new + num_vars_to_remove - i - 1;
      MSK_removevars(mosek_task_,num_vars_to_remove,vars_to_remove);      
    }

    // Add objective function
    int subj [numvar_new];
    double* c_val = &dual_c_[0];
    for (int i=0;i<numvar_new;++i) subj [i] = i;
    MSK_putclist(mosek_task_,numvar_new,subj,c_val);
    for (int i=m_+N_;i<numvar_new;++i) {
      MSK_putvarbound(mosek_task_,i,MSK_BK_LO,0.0,+MSK_INFINITY);
    }

    // Add equality constraints
    int* ptrb = &dual_A_colind_[0];
    int* ptre = &dual_A_colind_[1];
    int* asub = &dual_A_row_[0];
    double* aval = &dual_A_data_[0];
    MSK_putacollist(mosek_task_,numvar_new,subj,ptrb,ptre,asub,aval);
    for (int i=0;i<n_;++i) {
      MSK_putconbound(mosek_task_,i,MSK_BK_FX,-dual_b_[i],-dual_b_[i]);
    }

    // Solve SOCP
    MSKrescodee r;
    MSKrescodee trmcode;

    r = MSK_optimizetrm(mosek_task_,&trmcode);

    casadi_assert_message(r==MSK_RES_OK, "MosekInterface failed");

    MSK_solutionsummary(mosek_task_,MSK_STREAM_MSG);

    // Extract solution from MOSEK
    double primal_objective;
    double dual_objective;
    double* primal_solution = &output(SOCP_SOLVER_X).data()[0];
    double dual_solution [numvar_new];
    MSK_getdualobj(mosek_task_,MSK_SOL_ITR,&primal_objective);
    MSK_getprimalobj(mosek_task_,MSK_SOL_ITR,&dual_objective);
    MSK_gety(mosek_task_,MSK_SOL_ITR,primal_solution);
    MSK_getxx(mosek_task_,MSK_SOL_ITR,dual_solution);
    output(SOCP_SOLVER_COST).set(primal_objective);
    output(SOCP_SOLVER_DUAL_COST).set(dual_objective);
    // Change sign of primal solution
    std::for_each(output(SOCP_SOLVER_X).data().begin(),output(SOCP_SOLVER_X).data().end(),[](double& in){in=-in;});
    // Interpret dual solution
    std::copy(dual_solution,dual_solution+m_,output(SOCP_SOLVER_LAM_CONE).data().begin());
    // Fill lagrange multipliers on A*x and x with zeros
    std::fill(output(SOCP_SOLVER_LAM_A).begin(),output(SOCP_SOLVER_LAM_A).end(),0);
    std::fill(output(SOCP_SOLVER_LAM_X).begin(),output(SOCP_SOLVER_LAM_X).end(),0);
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
    MSK_getsolsta(mosek_task_,MSK_SOL_ITR,&solsta);
    stats_["solution_status"] = solutionStatus(solsta);

    // Problem status
    MSKprostae prosta;
    MSK_getprosta(mosek_task_,MSK_SOL_ITR,&prosta);
    stats_["problem_status"] = problemStatus(prosta);    

    // Termination reason
    stats_["termination_reason"] = trmcode;

  }

  void MosekInterface::convertToDualSocp() {
    /**

    Convert Second Order Cone Programming (SOCP) problem in standard form to one in a conic form.

    Dual:

    max          c' x
    x
    subject to
    || yi ||_2 <= ti  i = 1..m
    A x + b == 0
    lbx <= x
     
    Dimensions                       | Meaning in terms of primal variables
    ------------------------------------------------------------------------------------------------------------
    with x ( nx x 1)                 | [lag_ti' lag_yi' lag_lba' lag_uba' lag_lbx' lag_ubx']'
    c dense ( nx x 1 )               | [-fi' hi' LBA' -UBA' LBX' -UBX']'
    yi dense ( ni )                  | Lagrange multipliers for equality constraints: yi = Gi' x + hi
    ti dense ( m )                   | Lagrange multipliers for equality constraint: ti = ei' x + fi
    A  sparse ( n x nx )             | [-ei Gi -Alba' Auba' -Ilbx Iubx]
    b  dense ( n x 1 )               | [c]
    lbx dense ( nx x 1 )             | [-inf' 0']'
    nx = m + N + nlba + nuba + nlbx + nubx

    */

    /*****************************************************************
     * Define aliases                                                *
     *****************************************************************/
    DMatrix& primal_C           = input(SOCP_SOLVER_C);
    DMatrix& primal_G           = input(SOCP_SOLVER_G);
    DMatrix& primal_H           = input(SOCP_SOLVER_H);
    DMatrix& primal_E           = input(SOCP_SOLVER_E);
    DMatrix& primal_F           = input(SOCP_SOLVER_F);
    DMatrix& primal_A           = input(SOCP_SOLVER_A);
    DMatrix& primal_LBA         = input(SOCP_SOLVER_LBA);
    DMatrix& primal_UBA         = input(SOCP_SOLVER_UBA);
    DMatrix& primal_LBX         = input(SOCP_SOLVER_LBX);
    DMatrix& primal_UBX         = input(SOCP_SOLVER_UBX);
    const Sparsity& primal_A_sparse   = primal_A.sparsity();
    const Sparsity& primal_G_sparse   = primal_G.sparsity();

    /*****************************************************************
     * Categorize number of optimization variables for dual problem  *
     *****************************************************************/
    
    // Empty list of indices 
    primal_idx_lba_.resize(0);
    primal_idx_uba_.resize(0);
    primal_idx_lbx_.resize(0);
    primal_idx_ubx_.resize(0);

    // Loop over linear equality constraints
    for (int i=0;i<nc_;++i) {
      if (primal_LBA.at(i) != -std::numeric_limits<double>::infinity()){
        primal_idx_lba_.push_back(i);
      }
      if (primal_UBA.at(i) != std::numeric_limits<double>::infinity()){
        primal_idx_uba_.push_back(i);
      }
    }

    // Loop over simple bounds
    for (int i=0;i<n_;++i) {
      if (primal_LBX.at(i) != -std::numeric_limits<double>::infinity()){
        primal_idx_lbx_.push_back(i);
      }
      if (primal_UBX.at(i) != std::numeric_limits<double>::infinity()){
        primal_idx_ubx_.push_back(i);
      }
    }

    /*****************************************************************
     * Set up dual problem                                           *
     *****************************************************************/

    // c-vector: [fi' -hi' -LBA' UBA' -LBX' UBX']'
    dual_c_.resize(m_+N_);
    int start_idx = 0;
    int end_idx = m_;
    std::copy(primal_F.data().begin(),primal_F.data().end(),dual_c_.begin()+start_idx);
    std::for_each(dual_c_.begin()+start_idx,dual_c_.begin()+end_idx,[](double& in){in=-in;});
    start_idx = end_idx;
    end_idx += N_;
    std::copy(primal_H.data().begin(),primal_H.data().end(),dual_c_.begin()+start_idx);
    for (int i : primal_idx_lba_) dual_c_.push_back(primal_LBA.data()[i]);
    for (int i : primal_idx_uba_) dual_c_.push_back(-primal_UBA.data()[i]);
    for (int i : primal_idx_lbx_) dual_c_.push_back(primal_LBX.data()[i]);
    for (int i : primal_idx_ubx_) dual_c_.push_back(-primal_UBX.data()[i]);

    // A-matrix: [-ei Gi -Alba' Auba' -Ilbx Iubx]
    int begin_colind;
    int end_colind;
    dual_A_data_.resize(primal_E.size()+primal_G.size());
    dual_A_row_.resize(primal_E.size()+primal_G.size());
    dual_A_colind_.resize(m_+N_+1);

    // TODO: replace T()-call by casadi_trans()
    DMatrix primal_A_T = primal_A.T();
    const Sparsity& primal_A_T_sparse = primal_A_T.sparsity();
    start_idx = 0;
    end_idx = primal_E.size();
    std::copy(primal_E.data().begin(),primal_E.data().end(),dual_A_data_.begin()+start_idx);
    std::for_each(dual_A_data_.begin()+start_idx,dual_A_data_.begin()+end_idx,[](double& in){in=-in;});
    start_idx = end_idx;
    end_idx += primal_G.size();
    std::copy(primal_G.data().begin(),primal_G.data().end(),dual_A_data_.begin()+start_idx);
    for (int i : primal_idx_lba_) {
      begin_colind = primal_A_T_sparse.colind(i);
      end_colind = primal_A_T_sparse.colind(i+1);
      for (int ii=begin_colind; ii<end_colind;++ii){
        dual_A_data_.push_back(-primal_A_T.data()[ii]);
        dual_A_row_.push_back(primal_A_T_sparse.getRow()[ii]);
      }
      dual_A_colind_.push_back(dual_A_colind_.back()+end_colind-begin_colind);
    }
    for (int i : primal_idx_uba_) {
      begin_colind = primal_A_T_sparse.colind(i);
      end_colind = primal_A_T_sparse.colind(i+1);
      for (int ii=begin_colind; ii<end_colind;++ii){
        dual_A_data_.push_back(primal_A_T.data()[ii]);
        dual_A_row_.push_back(primal_A_T_sparse.getRow()[ii]);
      }
      dual_A_colind_.push_back(dual_A_colind_.back()+end_colind-begin_colind);
    }
    for (int i : primal_idx_lbx_) {
      dual_A_data_.push_back(-1);
      dual_A_row_.push_back(i);
      dual_A_colind_.push_back(dual_A_colind_.back()+1);
    }
    for (int i : primal_idx_ubx_) {
      dual_A_data_.push_back(1);
      dual_A_row_.push_back(i);
      dual_A_colind_.push_back(dual_A_colind_.back()+1);
    }

    // b-vector
    std::copy(primal_C.data().begin(),primal_C.data().end(),dual_b_.begin());

  }

  std::string MosekInterface::solutionStatus(MSKsolstae& solsta) {
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

  std::string MosekInterface::problemStatus(MSKprostae& prosta) {
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
