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

// Ipopt runtime independent of the Ipopt C interface, shared by the VM and generated code


// ipmc runtime independent of the ipmc C API, shared by the VM and generated code

// C-REPLACE "case SOLVER_RET_SUCCESS" "case 0"
// C-REPLACE "case SOLVER_RET_LIMITED" "case 2"
// C-REPLACE "case SOLVER_RET_NAN" "case 3"
// C-REPLACE "case SOLVER_RET_INFEASIBLE" "case 4"
// C-REPLACE "= SOLVER_RET_UNKNOWN;" "= 1;"

// SYMBOL "ipmc_return_status_string"
inline
const char* casadi_ipmc_return_status_string(casadi_int status) {
  // Literals from enum IpmcStatus in ipmc.h
  switch (status) {
    case 0: return "IPMC_SOLVED";
    case 1: return "IPMC_FAILED";
    case 2: return "IPMC_INFEASIBLE";
    case 3: return "IPMC_INTERRUPTED";
    case 100: return "IPMC_RETURNED_FROM_RESTO";
  }
  return "Unknown";
}

// SYMBOL "ipmc_stats_set_outcome"
inline
void casadi_ipmc_stats_set_outcome(struct casadi_stats_sink* s, casadi_int call,
    casadi_int status, casadi_int unified_return_status, casadi_int success,
    casadi_int iter_count) {
  const char* unified;
  if (!s) return;
  switch (unified_return_status) {
    case SOLVER_RET_SUCCESS: unified = "SOLVER_RET_SUCCESS"; break;
    case SOLVER_RET_LIMITED: unified = "SOLVER_RET_LIMITED"; break;
    case SOLVER_RET_NAN: unified = "SOLVER_RET_NAN"; break;
    case SOLVER_RET_INFEASIBLE: unified = "SOLVER_RET_INFEASIBLE"; break;
    default:
      unified = "SOLVER_RET_UNKNOWN";
      unified_return_status = SOLVER_RET_UNKNOWN;
  }
  casadi_stats_set_text(s, call, "return_status", casadi_ipmc_return_status_string(status));
  casadi_stats_set_int(s, call, "return_status_enum", status);
  casadi_stats_set_text(s, call, "unified_return_status", unified);
  casadi_stats_set_int(s, call, "unified_return_status_enum", unified_return_status);
  casadi_stats_set_bool(s, call, "success", success);
  casadi_stats_set_int(s, call, "iter_count", iter_count);
}
