//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

// Ipopt runtime without the Ipopt C interface: compiled into the virtual machine as well as
// emitted into generated code (with option "stats"), so both name return statuses and write
// stats records alike

// C-REPLACE "case SOLVER_RET_SUCCESS" "case 0"
// C-REPLACE "case SOLVER_RET_LIMITED" "case 2"
// C-REPLACE "case SOLVER_RET_NAN" "case 3"
// C-REPLACE "= SOLVER_RET_UNKNOWN;" "= 1;"

// SYMBOL "ipopt_return_status_string"
inline
const char* casadi_ipopt_return_status_string(casadi_int status) {
  // Name of an ApplicationReturnStatus; values as in IpReturnCodes_inc.h, as literals so that
  // no Ipopt version guards are needed
  switch (status) {
    case 0: return "Solve_Succeeded";
    case 1: return "Solved_To_Acceptable_Level";
    case 2: return "Infeasible_Problem_Detected";
    case 3: return "Search_Direction_Becomes_Too_Small";
    case 4: return "Diverging_Iterates";
    case 5: return "User_Requested_Stop";
    case 6: return "Feasible_Point_Found";
    case -1: return "Maximum_Iterations_Exceeded";
    case -2: return "Restoration_Failed";
    case -3: return "Error_In_Step_Computation";
    case -4: return "Maximum_CpuTime_Exceeded";
    case -5: return "Maximum_WallTime_Exceeded";
    case -10: return "Not_Enough_Degrees_Of_Freedom";
    case -11: return "Invalid_Problem_Definition";
    case -12: return "Invalid_Option";
    case -13: return "Invalid_Number_Detected";
    case -100: return "Unrecoverable_Exception";
    case -101: return "NonIpopt_Exception_Thrown";
    case -102: return "Insufficient_Memory";
    case -199: return "Internal_Error";
  }
  return "Unknown";
}

// FILTER-MACROS OFF
// The fields of an Ipopt iteration, in the order casadi_ipopt_stats_declare_fields gives them
enum casadi_ipopt_field {
  CASADI_IPOPT_ALG_MOD,
  CASADI_IPOPT_OBJ,
  CASADI_IPOPT_INF_PR,
  CASADI_IPOPT_INF_DU,
  CASADI_IPOPT_MU,
  CASADI_IPOPT_D_NORM,
  CASADI_IPOPT_REGULARIZATION_SIZE,
  CASADI_IPOPT_ALPHA_DU,
  CASADI_IPOPT_ALPHA_PR,
  CASADI_IPOPT_LS_TRIALS,
  CASADI_IPOPT_N_FIELDS
};
// FILTER-MACROS ON

// SYMBOL "ipopt_stats_declare_fields"
inline
void casadi_ipopt_stats_declare_fields(struct casadi_stats_sink* s, casadi_int call) {
  // Once per solve; same order as enum casadi_ipopt_field
  static const char* names[CASADI_IPOPT_N_FIELDS] = {"alg_mod", "obj", "inf_pr", "inf_du", "mu",
    "d_norm", "regularization_size", "alpha_du", "alpha_pr", "ls_trials"};
  casadi_stats_declare_fields(s, call, names, CASADI_IPOPT_N_FIELDS);
}

// SYMBOL "ipopt_stats_next_iteration"
template<typename T1>
casadi_int casadi_ipopt_stats_next_iteration(struct casadi_stats_sink* s, casadi_int call,
    casadi_int scope, casadi_int iter, casadi_int alg_mod, T1 obj, T1 inf_pr, T1 inf_du, T1 mu,
    T1 d_norm, T1 regularization_size, T1 alpha_du, T1 alpha_pr, casadi_int ls_trials) {
  // Close scope (pre, or the previous iteration) and open iteration iter: the new scope
  casadi_int it;
  if (!s) return scope;
  casadi_stats_end_scope(s, scope);
  it = casadi_stats_begin_iteration(s, call, iter);
  casadi_stats_set_field_text(s, it, CASADI_IPOPT_ALG_MOD, alg_mod ? "restoration" : "regular");
  casadi_stats_set_field_real(s, it, CASADI_IPOPT_OBJ, obj);
  casadi_stats_set_field_real(s, it, CASADI_IPOPT_INF_PR, inf_pr);
  casadi_stats_set_field_real(s, it, CASADI_IPOPT_INF_DU, inf_du);
  casadi_stats_set_field_real(s, it, CASADI_IPOPT_MU, mu);
  casadi_stats_set_field_real(s, it, CASADI_IPOPT_D_NORM, d_norm);
  casadi_stats_set_field_real(s, it, CASADI_IPOPT_REGULARIZATION_SIZE, regularization_size);
  casadi_stats_set_field_real(s, it, CASADI_IPOPT_ALPHA_DU, alpha_du);
  casadi_stats_set_field_real(s, it, CASADI_IPOPT_ALPHA_PR, alpha_pr);
  casadi_stats_set_field_int(s, it, CASADI_IPOPT_LS_TRIALS, ls_trials);
  return it;
}

// SYMBOL "ipopt_stats_set_outcome"
inline
void casadi_ipopt_stats_set_outcome(struct casadi_stats_sink* s, casadi_int call,
    casadi_int status, casadi_int unified_return_status, casadi_int success,
    casadi_int iter_count) {
  // status: Ipopt's ApplicationReturnStatus; each status as label and as enum value
  const char* unified;
  if (!s) return;
  switch (unified_return_status) {
    case SOLVER_RET_SUCCESS: unified = "SOLVER_RET_SUCCESS"; break;
    case SOLVER_RET_LIMITED: unified = "SOLVER_RET_LIMITED"; break;
    case SOLVER_RET_NAN: unified = "SOLVER_RET_NAN"; break;
    default:
      unified = "SOLVER_RET_UNKNOWN";
      unified_return_status = SOLVER_RET_UNKNOWN;
  }
  casadi_stats_put_text(s, call, CASADI_STATS_SET, "return_status",
    casadi_ipopt_return_status_string(status));
  casadi_stats_put_int(s, call, CASADI_STATS_SET, "return_status_enum", status);
  casadi_stats_put_text(s, call, CASADI_STATS_SET, "unified_return_status", unified);
  casadi_stats_put_int(s, call, CASADI_STATS_SET, "unified_return_status_enum",
    unified_return_status);
  casadi_stats_put_bool(s, call, CASADI_STATS_SET, "success", success);
  casadi_stats_put_int(s, call, CASADI_STATS_SET, "iter_count", iter_count);
}
