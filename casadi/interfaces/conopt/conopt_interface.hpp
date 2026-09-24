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

#ifndef CASADI_CONOPT_INTERFACE_HPP
#define CASADI_CONOPT_INTERFACE_HPP

#include "casadi/core/nlpsol_impl.hpp"
#include <casadi/interfaces/conopt/casadi_nlpsol_conopt_export.h>
#include <conopt.h>
#include <atomic>
#include <map>
#include <vector>
#include <string>
#include <utility>

namespace casadi {
  class ConoptInterface;

  enum class ConoptModelStatus {
    Unset              = 0,
    Optimal            = 1,
    LocallyOptimal     = 2,
    Unbounded          = 3,
    Infeasible         = 4,
    LocallyInfeasible  = 5,
    IntermediateInfeas = 6,
    IntermediateNonOpt = 7,
    UnknownError       = 12,
    ErrorNoSolution    = 13
  };

  // CONOPT's TYPE row-type convention (used for both TYPEX[] passed to CONOPT
  // and the interface's own conopt_type bookkeeping).
  enum class ConoptRowType {
    Equality     = 0,  // an equality constraint
    GreaterEqual = 1,  // a greater than or equal constraint
    LessEqual    = 2,  // a less than or equal constraint
    Free         = 3   // a free row
  };

  // CONOPT's IniStat=2 basis-status convention, shared by VSTA (variables) and
  // ESTA (constraint slacks) in the ReadMatrix callback.
  enum class ConoptBasisStatus {
    AtLower   = 0,  // initialized at lower bound
    AtUpper   = 1,  // initialized at upper bound
    Basic     = 2,  // initialized basic
    SuperBasic = 3  // initialized superbasic
  };

  enum class ConoptSolverStatus {
    Unset              = 0,
    NormalCompletion   = 1,
    IterationLimit     = 2,
    TimeLimit          = 3,
    TerminatedBySolver = 4,
    EvalErrorLimit     = 5,
    UserInterrupt      = 8,
    SetupFailure       = 9,
    MajorSolverError   = 10,
    MajorSolverErrorFeas = 11,
    SystemError        = 13,
    QuickModeTermination = 15
  };

  // One instruction of an SX function's evaluation tape (copied from
  // SXFunction::algorithm_; d is only meaningful for OP_CONST).
  struct ConoptInstr {
    int op, i0, i1, i2;
    double d;
  };

  // Evaluation tape of an SX function, built in init() so that FDEvalIni can
  // execute only the instructions the requested rows depend on.
  struct ConoptTape {
    std::vector<ConoptInstr> instr;
    // Instruction that produced each operand of an instruction (-1: none)
    std::vector<int> dep1, dep2;
    // OP_OUTPUT instruction of each output nonzero, per output
    std::vector<std::vector<int>> out_instr;
    casadi_int sz_w = 0;
  };

  // Outcome of collecting the instructions for one set of row slots: either
  // too much of the tape (evaluate in full) or the sorted instruction list.
  struct ConoptRowSet {
    bool heavy;
    std::vector<int> instr;
  };

  struct CASADI_NLPSOL_CONOPT_EXPORT ConoptMemory : public NlpsolMemory {
    const ConoptInterface& self;
    coiHandle_t cntvect;

    ConoptModelStatus modsta;
    ConoptSolverStatus solsta;
    int iter;
    std::string return_status;

    // Caching state for the evaluation block
    std::vector<double> cached_x;
    double cached_f;
    std::vector<double> cached_grad_f;
    std::vector<double> cached_g;
    std::vector<double> cached_jac_g;
    std::atomic<bool> cache_valid{false};
    std::atomic<bool> cache_valid_jac{false};  // true only when cached_jac_g was computed
    // Whether the cache holds values/Jacobian for cached_x from a previous
    // single-row FDEvalIni. Unlike cache_valid*, these survive FDEvalEnd so a
    // repeated single-row request at the same point can skip re-evaluation.
    bool stored_fun = false;
    bool stored_jac = false;
    bool nan_encountered;

    // Row-subset evaluation state (only used when the interface has tapes).
    // Row slots are CasADi constraint indices, with slot ng_ for the objective.
    // A slot's values (derivatives) are cached for cached_x when its stamp_fun
    // (stamp_jac) entry equals cur_stamp; bumping cur_stamp invalidates all.
    std::vector<unsigned> stamp_fun;
    std::vector<unsigned> stamp_jac;
    std::vector<unsigned> queued;  // dedupes slots within one FDEvalIni
    unsigned cur_stamp = 0;
    unsigned cur_queue = 0;
    std::vector<int> pending_rows;
    std::vector<double> tape_w;
    std::vector<char> tape_mark;   // all zero between evaluations
    std::vector<int> tape_stack;
    std::vector<int> tape_list;
    // Row sets seen so far, keyed by sorted slots, per tape (0: fg, 1: fg_jac).
    // Structural, so kept across solves; cleared when row_set_ints exceeds a cap.
    std::map<std::vector<int>, ConoptRowSet> row_sets[2];
    size_t row_set_ints = 0;
    // Statistics
    casadi_int n_eval_subset = 0;   // FDEvalIni calls evaluated on a row subset
    casadi_int n_eval_full = 0;     // FDEvalIni calls evaluated in full
    casadi_int n_eval_reused = 0;   // FDEvalIni calls answered from the cache
    casadi_int n_eval_on_demand = 0;  // FDEval rows missing from the FDEvalIni list
    double subset_instr_frac = 0;   // sum over subset evaluations of executed/total

    // Invalidates the row-subset cache
    void bump_stamp();

    // Options handling
    std::vector<std::pair<std::string, GenericType>> custom_options;

    // Constant Jacobian values for linear (NLFLAG=0) entries (populated in solve())
    std::vector<double> const_jac_vals;
    // Constant objective gradient values for linear (NLFLAG=0) entries
    std::vector<double> gradf_const_vals;
    // Scratch buffer for the linear part of G at x0 per row.
    std::vector<double> linear_at_x0;

    // Range-constraint expansion state (recomputed each solve)
    int ng_expanded;
    int numnz_expanded;
    std::vector<int> row_nnz;            // nnz per CasADi row, persistent across solves
    std::vector<int> conopt_to_casadi;   // CONOPT constraint row (0-indexed) → CasADi index
    std::vector<int> casadi_to_conopt_lb_row;  // CasADi index → CONOPT row for lb (or only) side
    // CasADi index -> CONOPT row for ub side (range only); -1 otherwise
    std::vector<int> casadi_to_conopt_ub_row;
    std::vector<ConoptRowType> conopt_type;  // CONOPT TYPE for each expanded row
    std::vector<double> conopt_rhs;      // CONOPT RHS for each expanded row

    // Multiplier buffer for Hessian callback (avoids per-call allocation)
    std::vector<double> hess_lam_g_;

    // Stored constant objective value when the interface is in feasibility mode
    double obj_const_;
    // Per-row constant absorbed into conopt_rhs, restored in the reported g.
    std::vector<double> row_const_;
    // Affine objective constant, restored in the reported f.
    double obj_const_lin_;

    ConoptMemory(const ConoptInterface& interface);
    ~ConoptMemory();
  };

  class CASADI_NLPSOL_CONOPT_EXPORT ConoptInterface : public Nlpsol {
  public:
    explicit ConoptInterface(const std::string& name, const Function& nlp);
    ~ConoptInterface() override;

    const char* plugin_name() const override { return "conopt"; }
    std::string class_name() const override { return "ConoptInterface"; }

    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new ConoptInterface(name, nlp);
    }

    static const Options options_;
    const Options& get_options() const override { return options_; }

    /// A documentation string
    static const std::string meta_doc;

    void init(const Dict& opts) override;
    void* alloc_mem() const override { return new ConoptMemory(*this); }
    int init_mem(void* mem) const override;
    void free_mem(void* mem) const override;
    void set_work(void* mem, const double**& arg, double**& res,
                  casadi_int*& iw, double*& w) const override;
    int solve(void* mem) const override;
    Dict get_stats(void* mem) const override;

    // Grows conopt_to_casadi/conopt_type/conopt_rhs (which always stay the same
    // size as each other) when they are full, used by solve()'s row-expansion loop.
    void ensure_row_capacity(ConoptMemory* m, casadi_int remaining_rows) const;

    // Clears jacg_nlflag_/gradf_nlflag_ for columns absent from hesslag_sp_
    // (CONOPT rejects nonlinear columns missing from the Hessian structure).
    void refine_nlflags_with_hessian();

    // Builds tape_fg_/tape_fg_jac_ when subset_eval_ is set and the evaluation
    // functions are SX without calls; sets has_tape_.
    void build_tapes();
    static bool build_tape(const Function& f, ConoptTape& t);

    // Evaluates the row slots in m->pending_rows at m->cached_x using the tapes.
    // Returns 0 on success, 1 on an evaluation error (NaN/Inf), and -1 when the
    // rows depend on so much of the tape that a full evaluation is cheaper.
    int eval_rows(ConoptMemory* m, bool need_jac) const;

    // Full evaluation at m->cached_x into the cache (nlp_fg or nlp_fg_jac)
    int eval_full(ConoptMemory* m, bool need_jac) const;

    // FDEvalIni when row-subset evaluation is active
    int fdevalini_subset(ConoptMemory* m, const double X[], const int ROWLIST[],
                         int MODE, int LISTSIZE) const;

    // Evaluate only the rows CONOPT requests (see the 'subset_eval' option)
    bool subset_eval_;
    bool has_tape_ = false;
    ConoptTape tape_fg_;      // nlp_fg: f, g
    ConoptTape tape_fg_jac_;  // nlp_fg_jac: f, grad:f:x, g, jac:g:x

    // Sparsities for the problem components
    Sparsity gradf_sp_;
    Sparsity jacg_sp_;
    Sparsity hesslag_sp_;
    bool exact_hessian_;
    bool warm_start_;
    bool debug_;
    Dict opts_; // CONOPT specific options
    std::string optfile_; // Path to CONOPT option file (for string-valued CR-cells)

    // License passed via the 'license' option; takes precedence over the
    // CONOPT_LICENSE_* environment variables. Not serialized.
    bool has_license_ = false;
    int license_int_[3] = {0, 0, 0};
    std::string license_text_;

    // Per-column flag: true if the objective gradient has a nonzero in that column
    std::vector<bool> gradf_col_flag_;

    // Row-indexed structure for fast Jacobian scatter in cb_fd_eval (nonlinear entries only)
    std::vector<int> jacg_rowstart_;  // size ng_+1; nonlinear nnz prefix sums per row
    std::vector<int> jacg_nzidx_;    // CCS indices of nonlinear entries, in row-major order
    std::vector<int> jacg_col_;      // column (variable) index of each entry in jacg_nzidx_

    // Per-nonzero linearity flag (0 = constant/linear, 1 = nonlinear), CCS order
    std::vector<int> jacg_nlflag_;
    bool has_linear_jac_;

    // Objective gradient linearity: flag per nonzero of gradf_sp_
    std::vector<int> gradf_nlflag_;
    bool has_linear_gradf_;
    // Maps variable index → nonzero index in gradf_sp_ (-1 if absent)
    std::vector<casadi_int> gradf_col_to_nz_;

    // Serialization and Deserialization
    void serialize_body(SerializingStream &s) const override;
    static ProtoFunction* deserialize(DeserializingStream& s) { return new ConoptInterface(s); }

    // CONOPT Mandatory Callbacks
    static int COI_CALLCONV cb_read_matrix(double LOWER[], double CURR[], double UPPER[],
                                            int VSTA[], int TYPEX[], double RHS[],
                                            int ESTA[], int COLSTA[], int ROWNO[],
                                            double VALUE[], int NLFLAG[], int NUMVAR,
                                            int NUMCON, int NUMNZ, void* USRMEM);
    static int COI_CALLCONV cb_fdevalini(const double X[], const int ROWLIST[], int MODE,
                                          int LISTSIZE, int NUMTHREAD, int IGNERR,
                                          int* ERRCNT, int NUMVAR, void* USRMEM);
    static int COI_CALLCONV cb_fd_eval(const double X[], double* G, double JAC[],
                                        int ROWNO, const int JACNUM[], int MODE, int IGNERR,
                                        int* ERRCNT, int NUMVAR, int NUMJAC, int THREAD,
                                        void* USRMEM);
    static int COI_CALLCONV cb_fdevalend(int IGNERR, int* ERRCNT, void* USRMEM);
    static int COI_CALLCONV cb_status(int MODSTA, int SOLSTA, int ITER, double OBJVAL,
                                       void* USRMEM);
    static int COI_CALLCONV cb_solution(const double XVAL[], const double XMAR[],
                                         const int XBAS[], const int XSTA[],
                                         const double YVAL[], const double YMAR[],
                                         const int YBAS[], const int YSTA[],
                                         int NUMVAR, int NUMCON, void* USRMEM);

    // Logging and Messages
    static int COI_CALLCONV cb_message(int SMSG, int DMSG, int NMSG, char* MSGV[],
                                        void* USRMEM);
    static int COI_CALLCONV cb_errmsg(int ROWNO, int COLNO, int POSNO, const char* MSG,
                                       void* USRMEM);

    // Options and Progress
    static int COI_CALLCONV cb_option(int NCALL, double* RVAL, int* IVAL, int* LVAL,
                                       char* NAME, void* USRMEM);
    static int COI_CALLCONV cb_progress(int LEN_INT, const int INTX[], int LEN_RL,
                                         const double RL[], const double X[], void* USRMEM);

    // CONOPT 2nd Order Callbacks
    static int COI_CALLCONV cb_2dlagrstr(int HSRW[], int HSCL[], int* NODRV, int NUMVAR,
                                          int NUMCON, int NHESS, void* USRMEM);
    static int COI_CALLCONV cb_2dlagrval(const double X[], const double U[],
                                          const int HSRW[], const int HSCL[], double HSVL[],
                                          int* NODRV, int NUMVAR, int NUMCON, int NHESS,
                                          void* USRMEM);

  protected:
    explicit ConoptInterface(DeserializingStream& s);
  };
}  // namespace casadi
#endif // CASADI_CONOPT_INTERFACE_HPP
