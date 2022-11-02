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


#include "knitro_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.thread.h>
#else // CASADI_WITH_THREAD_MINGW
#include <thread>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_KNITRO_EXPORT
  casadi_register_nlpsol_knitro(Nlpsol::Plugin* plugin) {
    plugin->creator = KnitroInterface::creator;
    plugin->name = "knitro";
    plugin->doc = KnitroInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &KnitroInterface::options_;
    plugin->deserialize = &KnitroInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_KNITRO_EXPORT casadi_load_nlpsol_knitro() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_knitro);
  }

  KnitroInterface::KnitroInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }


  KnitroInterface::~KnitroInterface() {
    clear_mem();
  }

  const Options KnitroInterface::options_
  = {{&Nlpsol::options_},
     {{"knitro",
       {OT_DICT,
        "Options to be passed to KNITRO"}},
      {"detect_linear_constraints",
       {OT_BOOL,
        "Detect type of constraints"}},
      {"contype",
       {OT_INTVECTOR,
        "Type of constraint"}},
      {"complem_variables",
       {OT_INTVECTORVECTOR,
        "List of complementary constraints on simple bounds. "
        "Pair (i, j) encodes complementarity between the bounds on variable i and variable j."}}
     }
  };

  void KnitroInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);
    bool detect_linear_constraints = true;
    std::vector< std::vector<casadi_int> > complem_variables;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="knitro") {
        opts_ = op.second;
      } else if (op.first=="contype") {
        contype_ = op.second;
      } else if (op.first=="detect_linear_constraints") {
        detect_linear_constraints = op.second;
      } else if (op.first=="complem_variables") {
        complem_variables = op.second;
      }
    }

    // Type of constraints, general by default
    if (contype_.empty()) {
      contype_.resize(ng_, KN_CONTYPE_GENERAL);
      if (detect_linear_constraints) {
        std::vector<bool> nl_g = oracle_.which_depends("x", {"g"}, 2, true);
        for (casadi_int i=0;i<ng_;++i)
          contype_[i] = nl_g[i] ? KN_CONTYPE_GENERAL : KN_CONTYPE_LINEAR;
      }
    }

    casadi_assert_dev(contype_.size()==ng_);


    comp_type_.resize(complem_variables.size(), KN_CCTYPE_VARVAR);
    comp_i1_.reserve(complem_variables.size());
    comp_i2_.reserve(complem_variables.size());
    for (auto && e : complem_variables) {
      casadi_assert(e.size()==2, "Complementary constraints must come in pairs.");
      casadi_assert(e[0]>=0, "Invalid variable index.");
      casadi_assert(e[1]>=0, "Invalid variable index.");
      casadi_assert(e[0]<nx_, "Invalid variable index.");
      casadi_assert(e[1]<nx_, "Invalid variable index.");
      comp_i1_.push_back(e[0]);
      comp_i2_.push_back(e[1]);
    }

    // Setup NLP functions
    create_function("nlp_fg", {"x", "p"}, {"f", "g"});
    Function gf_jg_fcn = create_function("nlp_gf_jg", {"x", "p"}, {"grad:f:x", "jac:g:x"});
    jacg_sp_ = gf_jg_fcn.sparsity_out(1);
    Function hess_l_fcn = create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                                  {"hess:gamma:x:x"},
                                  {{"gamma", {"f", "g"}}});
    hesslag_sp_ = hess_l_fcn.sparsity_out(0);

    unsigned int hc = 0;
    #ifdef CASADI_WITH_THREAD
    //may return 0 when not able to detect. If it's the case, then return 8.
    hc = std::thread::hardware_concurrency();
    #endif
    int processor_count = hc ? hc : 8;
    //Obtain maximum number of threads needed
    int ms_numthreads = 1;
    int findiff_numthreads = 1;
    int numthreads = 1;
    int mip_numthreads = 1;
    for (auto&& op : opts_) {
      if (op.first=="ms_numthreads") {
        ms_numthreads = op.second;
      }
      if (op.first=="findiff_numthreads") {
        findiff_numthreads = op.second;
      }
      if (op.first=="numthreads") {
        numthreads = op.second;
      }
      if (op.first=="mip_numthreads") {
        mip_numthreads = op.second;
      }
    }
    max_num_threads_ = std::max({processor_count, ms_numthreads, findiff_numthreads, numthreads, mip_numthreads});

    // Allocate persistent memory
    alloc_w(nx_, true); // wlbx_
    alloc_w(nx_, true); // wubx_
    alloc_w(ng_, true); // wlbg_
    alloc_w(ng_, true); // wubg_
  }

  int KnitroInterface::init_mem(void* mem) const {
    return Nlpsol::init_mem(mem);
    //auto m = static_cast<KnitroMemory*>(mem);

    // Commented out since I have not found out how to change the bounds
    // Allocate KNITRO memory block
    //  m.kc = KN_new();
  }

  void KnitroInterface::set_work(void* mem, const double**& arg, double**& res,
                                 casadi_int*& iw, double*& w) const {
    auto m = static_cast<KnitroMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Copy inputs to temporary arrays
    m->wlbx = w; w += nx_;
    m->wubx = w; w += nx_;
    m->wlbg = w; w += ng_;
    m->wubg = w; w += ng_;
  }

  int casadi_KN_puts(const char * const str, void * const userParams) {
    std::string s(str);
    uout() << s << std::flush;
    return s.size();
  }

  int KnitroInterface::solve(void* mem) const {
    auto m = static_cast<KnitroMemory*>(mem);
    auto d_nlp = &m->d_nlp;
    casadi_int status;

    // Allocate KNITRO memory block (move back to init!)
    casadi_assert_dev(m->kc==nullptr);
    status = KN_new(&m->kc);
    casadi_assert_dev(m->kc!=nullptr);

    status = KN_set_puts_callback(m->kc, casadi_KN_puts, nullptr);
    casadi_assert(status == 0, "KN_set_puts_callback failed");

    // Jacobian sparsity
    vector<int> Jcol, Jrow;
    if (!jacg_sp_.is_null()) {
      assign_vector(jacg_sp_.get_col(), Jcol);
      assign_vector(jacg_sp_.get_row(), Jrow);
    }

    // Hessian sparsity
    casadi_int nnzH = hesslag_sp_.is_null() ? 0 : hesslag_sp_.nnz();
    vector<int> Hcol, Hrow;
    if (nnzH>0) {
      assign_vector(hesslag_sp_.get_col(), Hcol);
      assign_vector(hesslag_sp_.get_row(), Hrow);
      status = KN_set_int_param(m->kc, KN_PARAM_HESSOPT, KN_HESSOPT_EXACT);
      casadi_assert(status==0, "KN_set_int_param failed");
    } else {
      status = KN_set_int_param(m->kc, KN_PARAM_HESSOPT, KN_HESSOPT_LBFGS);
      casadi_assert(status==0, "KN_set_int_param failed");
    }

    // Pass user set options
    for (auto&& op : opts_) {
      int param_id;
      casadi_assert(KN_get_param_id(m->kc, op.first.c_str(), &param_id)==0,
        "Unknown parameter '" + op.first + "'.");

      int param_type;
      casadi_assert(!KN_get_param_type(m->kc, param_id, &param_type),
        "Error when setting option '" + op.first + "'.");

      switch (param_type) {
        case KN_PARAMTYPE_INTEGER:
          casadi_assert(!KN_set_int_param(m->kc, param_id, op.second),
            "Error when setting option '" + op.first + "'.");
          continue;
        case KN_PARAMTYPE_FLOAT:
          casadi_assert(!KN_set_double_param(m->kc, param_id, op.second),
            "Error when setting option '" + op.first + "'.");
          continue;
        case KN_PARAMTYPE_STRING:
          {
            string str = op.second.to_string();
            casadi_assert(!KN_set_char_param(m->kc, param_id, str.c_str()),
              "Error when setting option '" + op.first + "'.");
          }
          continue;
        default:
          casadi_error("Error when setting option '" + op.first + "'.");
      }
    }

    // "Correct" upper and lower bounds
    casadi_copy(d_nlp->lbz, nx_, m->wlbx);
    casadi_copy(d_nlp->ubz, nx_, m->wubx);
    casadi_copy(d_nlp->lbz+nx_, ng_, m->wlbg);
    casadi_copy(d_nlp->ubz+nx_, ng_, m->wubg);
    for (casadi_int i=0; i<nx_; ++i) if (isinf(m->wlbx[i])) m->wlbx[i] = -KN_INFINITY;
    for (casadi_int i=0; i<nx_; ++i) if (isinf(m->wubx[i])) m->wubx[i] =  KN_INFINITY;
    for (casadi_int i=0; i<ng_; ++i) if (isinf(m->wlbg[i])) m->wlbg[i] = -KN_INFINITY;
    for (casadi_int i=0; i<ng_; ++i) if (isinf(m->wubg[i])) m->wubg[i] =  KN_INFINITY;

    vector<int> xindex(nx_);
    vector<int> gindex(ng_);
    iota (begin(xindex), end(xindex), 0);
    iota (begin(gindex), end(gindex), 0);

    status = KN_add_vars(m->kc, nx_, get_ptr(xindex));
    casadi_assert(status==0, "KN_add_vars failed");
    status = KN_set_obj_goal(m->kc, KN_OBJGOAL_MINIMIZE);
    casadi_assert(status==0, "KN_set_obj_goal failed");
    status = KN_set_var_lobnds_all(m->kc, m->wlbx);
    casadi_assert(status==0, "KN_set_var_lobnds failed");
    status = KN_set_var_upbnds_all(m->kc, m->wubx);
    casadi_assert(status==0, "KN_set_var_upbnds failed");
    status = KN_add_cons(m->kc, ng_, get_ptr(gindex));
    casadi_assert(status==0, "KN_add_cons failed");
    status = KN_set_con_lobnds_all(m->kc, m->wlbg);
    casadi_assert(status==0, "KN_set_con_lobnds failed");
    status = KN_set_con_upbnds_all(m->kc, m->wubg);
    casadi_assert(status==0, "KN_set_con_upbnds failed");
    status = KN_set_var_primal_init_values_all(m->kc, d_nlp->z);
    casadi_assert(status==0, "KN_set_var_primal_init_values failed");
    if (mi_) {
        // Types of variables
        vector<int> vtype;
        vtype.reserve(nx_);
        for (auto&& e : discrete_) {
            vtype.push_back(e ? KN_VARTYPE_INTEGER : KN_VARTYPE_CONTINUOUS);
        }
        status = KN_set_var_types_all(m->kc, get_ptr(vtype));
        casadi_assert(status==0, "KN_set_var_types failed");
    }

    // Complementarity constraints
    status = KN_set_compcons(m->kc, comp_i1_.size(),
      get_ptr(comp_type_), get_ptr(comp_i1_), get_ptr(comp_i2_));
    casadi_assert(status==0, "KN_set_compcons failed");

    // Register callback functions
    status = KN_add_eval_callback(m->kc, true, ng_, get_ptr(gindex), &callback, &m->cb);
    casadi_assert(status==0, "KN_add_eval_callback failed");

    status = KN_set_cb_grad(m->kc, m->cb, nx_, get_ptr(xindex), Jcol.size(), get_ptr(Jrow), get_ptr(Jcol), &callback);
    casadi_assert(status==0, "KN_set_cb_grad failed");

    if (nnzH>0) {
      status = KN_set_cb_hess(m->kc, m->cb, nnzH, get_ptr(Hrow), get_ptr(Hcol), &callback);
      casadi_assert(status==0, "KN_set_cb_hess failed");
    }

    status = KN_set_cb_user_params(m->kc, m->cb, static_cast<void*>(m));
    casadi_assert(status==0, "KN_set_cb_user_params failed");

    // NumThreads to 1 to prevent segmentation fault
    status = KN_set_int_param(m->kc, KN_PARAM_NUMTHREADS, 1);
    casadi_assert(status==0, "KN_set_cb_user_params failed");

    // Lagrange multipliers
    vector<double> lambda(nx_+ng_);

    // objective solution
    double objSol;

    // Solve NLP
    status = KN_solve(m->kc);
    int statusKnitro = int(status);

    m->return_status = return_codes(status);
    m->success = status==KN_RC_OPTIMAL_OR_SATISFACTORY ||
                 status==KN_RC_NEAR_OPT;
    if (status==KN_RC_ITER_LIMIT_FEAS  ||
        status==KN_RC_TIME_LIMIT_FEAS  ||
        status==KN_RC_FEVAL_LIMIT_FEAS ||
        status==KN_RC_ITER_LIMIT_INFEAS  ||
        status==KN_RC_TIME_LIMIT_INFEAS  ||
        status==KN_RC_FEVAL_LIMIT_INFEAS)
      m->unified_return_status = SOLVER_RET_LIMITED;

    // Output optimal cost
    casadi_int error;
    error = KN_get_solution(m->kc, &statusKnitro, &objSol, d_nlp->z, get_ptr(lambda));
    casadi_assert(error == 0, "KN_get_solution failed");
    // Output dual solution
    casadi_copy(get_ptr(lambda), ng_, d_nlp->lam + nx_);
    casadi_copy(get_ptr(lambda)+ng_, nx_, d_nlp->lam);

    d_nlp->f = objSol;
    
    // Calculate constraints
    if (ng_>0) {
      m->arg[0] = d_nlp->z;
      m->arg[1] = d_nlp->p;
      m->res[0] = nullptr;
      m->res[1] = d_nlp->z + nx_;
      calc_function(m, "nlp_fg");
    }

    // Free memory (move to destructor!)
    status = KN_free(&m->kc);
    casadi_assert(status == 0, "KN_free failed");
    m->kc = nullptr;
    return 0;
  }

  int KnitroInterface::callback(KN_context_ptr kc,
                   CB_context_ptr             cb,
                   KN_eval_request_ptr const  evalRequest,
                   KN_eval_result_ptr  const  evalResult,
                   void             *  const  userParams) {
    try {
      int thread_id = evalRequest->threadID;
      // Get a pointer to the calling object
      auto m = static_cast<KnitroMemory*>(userParams);
      // Get thread local memory
      auto ml = m->thread_local_mem.at(thread_id);
      auto d_nlp = &m->d_nlp;
      const double *x;
      // Direct to the correct function
      switch (evalRequest->type) {
      case KN_RC_EVALFC:
      double *obj;
      double *c;
      x = evalRequest->x;
      obj = evalResult->obj;
      c = evalResult->c;
      ml->arg[0] = x;
      ml->arg[1] = d_nlp->p;
      ml->res[0] = obj;
      ml->res[1] = c;
      if (m->self.calc_function(m, "nlp_fg", nullptr, thread_id)) return KN_RC_EVAL_ERR;
      break;
      case KN_RC_EVALGA:
      double *objGrad;
      double *jac;
      x = evalRequest->x;
      objGrad = evalResult->objGrad;
      jac = evalResult->jac;
      ml->arg[0] = x;
      ml->arg[1] = d_nlp->p;
      ml->res[0] = objGrad;
      ml->res[1] = jac;
      if (m->self.calc_function(m, "nlp_gf_jg", nullptr, thread_id)) return KN_RC_EVAL_ERR;
      break;
      case KN_RC_EVALH_NO_F:
      case KN_RC_EVALH:
      const double *lambda;
      double sigma;
      double *hess;
      x = evalRequest->x;
      hess = evalResult->hess;
      lambda = evalRequest->lambda;
      sigma = *(evalRequest->sigma);
      ml->arg[0] = x;
      ml->arg[1] = d_nlp->p;
      ml->arg[2] = &sigma;
      ml->arg[3] = lambda;
      ml->res[0] = hess;
      if (m->self.calc_function(m, "nlp_hess_l", nullptr, thread_id)) {casadi_error("calc_hess_l failed");}
      break;
      default:
        casadi_error("KnitroInterface::callback: unknown method");
      }

      return 0;
    } catch(KeyboardInterruptException& ex) {
      return KN_RC_USER_TERMINATION;
    } catch(exception& ex) {
      uerr() << "KnitroInterface::callback caught exception: "
                               << ex.what() << endl;
      return -1;
    }

  }

  const char* KnitroInterface::return_codes(int flag) {
    switch (flag) {
    case KN_RC_OPTIMAL_OR_SATISFACTORY: return "KN_RC_OPTIMAL_OR_SATISFACTORY";
    case KN_RC_NEAR_OPT: return "KN_RC_NEAR_OPT";
    case KN_RC_FEAS_XTOL: return "KN_RC_FEAS_XTOL";
    case KN_RC_FEAS_NO_IMPROVE: return "KN_RC_FEAS_NO_IMPROVE";
    case KN_RC_FEAS_FTOL: return "KN_RC_FEAS_FTOL";
    case KN_RC_INFEASIBLE: return "KN_RC_INFEASIBLE";
    case KN_RC_INFEAS_XTOL: return "KN_RC_INFEAS_XTOL";
    case KN_RC_INFEAS_NO_IMPROVE: return "KN_RC_INFEAS_NO_IMPROVE";
    case KN_RC_INFEAS_MULTISTART: return "KN_RC_INFEAS_MULTISTART";
    case KN_RC_INFEAS_CON_BOUNDS: return "KN_RC_INFEAS_CON_BOUNDS";
    case KN_RC_INFEAS_VAR_BOUNDS: return "KN_RC_INFEAS_VAR_BOUNDS";
    case KN_RC_UNBOUNDED: return "KN_RC_UNBOUNDED";
    case KN_RC_ITER_LIMIT_FEAS: return "KN_RC_ITER_LIMIT_FEAS";
    case KN_RC_TIME_LIMIT_FEAS: return "KN_RC_TIME_LIMIT_FEAS";
    case KN_RC_FEVAL_LIMIT_FEAS: return "KN_RC_FEVAL_LIMIT_FEAS";
    case KN_RC_MIP_EXH_FEAS: return "KN_RC_MIP_EXH_FEAS";
    case KN_RC_MIP_TERM_FEAS: return "KN_RC_MIP_TERM_FEAS";
    case KN_RC_MIP_SOLVE_LIMIT_FEAS: return "KN_RC_MIP_SOLVE_LIMIT_FEAS";
    case KN_RC_MIP_NODE_LIMIT_FEAS: return "KN_RC_MIP_NODE_LIMIT_FEAS";
    case KN_RC_ITER_LIMIT_INFEAS: return "KN_RC_ITER_LIMIT_INFEAS";
    case KN_RC_TIME_LIMIT_INFEAS: return "KN_RC_TIME_LIMIT_INFEAS";
    case KN_RC_FEVAL_LIMIT_INFEAS: return "KN_RC_FEVAL_LIMIT_INFEAS";
    case KN_RC_MIP_EXH_INFEAS: return "KN_RC_MIP_EXH_INFEAS";
    case KN_RC_MIP_SOLVE_LIMIT_INFEAS: return "KN_RC_MIP_SOLVE_LIMIT_INFEAS";
    case KN_RC_MIP_NODE_LIMIT_INFEAS: return "KN_RC_MIP_NODE_LIMIT_INFEAS";
    case KN_RC_CALLBACK_ERR: return "KN_RC_CALLBACK_ERR";
    case KN_RC_LP_SOLVER_ERR: return "KN_RC_LP_SOLVER_ERR";
    case KN_RC_EVAL_ERR: return "KN_RC_EVAL_ERR";
    case KN_RC_OUT_OF_MEMORY: return "KN_RC_OUT_OF_MEMORY";
    case KN_RC_USER_TERMINATION: return "KN_RC_USER_TERMINATION";
    case KN_RC_OPEN_FILE_ERR: return "KN_RC_OPEN_FILE_ERR";
    case KN_RC_BAD_N_OR_F: return "KN_RC_BAD_N_OR_F";
    case KN_RC_BAD_CONSTRAINT: return "KN_RC_BAD_CONSTRAINT";
    case KN_RC_BAD_JACOBIAN: return "KN_RC_BAD_JACOBIAN";
    case KN_RC_BAD_HESSIAN: return "KN_RC_BAD_HESSIAN";
    case KN_RC_BAD_CON_INDEX: return "KN_RC_BAD_CON_INDEX";
    case KN_RC_BAD_JAC_INDEX: return "KN_RC_BAD_JAC_INDEX";
    case KN_RC_BAD_HESS_INDEX: return "KN_RC_BAD_HESS_INDEX";
    case KN_RC_BAD_CON_BOUNDS: return "KN_RC_BAD_CON_BOUNDS";
    case KN_RC_BAD_VAR_BOUNDS: return "KN_RC_BAD_VAR_BOUNDS";
    case KN_RC_ILLEGAL_CALL: return "KN_RC_ILLEGAL_CALL";
    case KN_RC_BAD_KCPTR: return "KN_RC_BAD_KCPTR";
    case KN_RC_NULL_POINTER: return "KN_RC_NULL_POINTER";
    case KN_RC_BAD_INIT_VALUE: return "KN_RC_BAD_INIT_VALUE";
    case KN_RC_BAD_PARAMINPUT: return "KN_RC_BAD_PARAMINPUT";
    case KN_RC_LINEAR_SOLVER_ERR: return "KN_RC_LINEAR_SOLVER_ERR";
    case KN_RC_DERIV_CHECK_FAILED: return "KN_RC_DERIV_CHECK_FAILED";
    case KN_RC_DERIV_CHECK_TERMINATE: return "KN_RC_DERIV_CHECK_TERMINATE";
    case KN_RC_INTERNAL_ERROR: return "KN_RC_INTERNAL_ERROR";
    }
    return nullptr;
  }

  Dict KnitroInterface::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<KnitroMemory*>(mem);
    stats["return_status"] = m->return_status;

    return stats;
  }

  KnitroInterface::KnitroInterface(DeserializingStream& s) : Nlpsol(s) {
    s.version("KnitroInterface", 1);
    s.unpack("KnitroInterface::contype", contype_);
    s.unpack("KnitroInterface::comp_type", comp_type_);
    s.unpack("KnitroInterface::comp_i1", comp_i1_);
    s.unpack("KnitroInterface::comp_i2", comp_i2_);
    s.unpack("KnitroInterface::opts", opts_);
    s.unpack("KnitroInterface::jacg_sp", jacg_sp_);
    s.unpack("KnitroInterface::hesslag_sp", hesslag_sp_);
  }

  void KnitroInterface::serialize_body(SerializingStream &s) const {
    Nlpsol::serialize_body(s);
    s.version("KnitroInterface", 1);
    s.pack("KnitroInterface::contype", contype_);
    s.pack("KnitroInterface::comp_type", comp_type_);
    s.pack("KnitroInterface::comp_i1", comp_i1_);
    s.pack("KnitroInterface::comp_i2", comp_i2_);
    s.pack("KnitroInterface::opts", opts_);
    s.pack("KnitroInterface::jacg_sp", jacg_sp_);
    s.pack("KnitroInterface::hesslag_sp", hesslag_sp_);
  }

  KnitroMemory::KnitroMemory(const KnitroInterface& self) : self(self) {
    this->kc = nullptr;
  }

  KnitroMemory::~KnitroMemory() {
    // Currently no persistent memory since KNITRO requires knowledge of nature of bounds
    // if (this->kc) {
    //   KTR_free(&this->kc);
    //}
  }

} // namespace casadi
