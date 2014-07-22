/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include <stdio.h>
#include <string.h>
#include <ctime>
#include <utility>
#include <string>
#include <algorithm>
#include <vector>
#include <iomanip>

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/matrix/sparsity_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

#include "snopt_internal.hpp"
#include "wsnopt.hpp"

namespace casadi {

  extern "C"
  int CASADI_NLPSOLVER_SNOPT_EXPORT
  casadi_register_nlpsolver_snopt(NlpSolverInternal::Plugin* plugin) {
    plugin->creator = SnoptInternal::creator;
    plugin->name = "snopt";
    plugin->doc = SnoptInternal::meta_doc.c_str();
    plugin->version = 20;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOLVER_SNOPT_EXPORT casadi_load_nlpsolver_snopt() {
    NlpSolverInternal::registerPlugin(casadi_register_nlpsolver_snopt);
  }

  SnoptInternal::SnoptInternal(const Function& nlp) : NlpSolverInternal(nlp) {
    addOption("detect_linear", OT_BOOLEAN, true,
              "Make an effort to treat linear constraints and linear variables specially.");

    // Monitors
    addOption("monitor", OT_STRINGVECTOR, GenericType(),  "", "eval_nlp|setup_nlp", true);

    addOption("print_time", OT_BOOLEAN, true, "print information about execution time");
    addOption("_start",  OT_STRING, "Cold",  "", "Cold|Warm");
    addOption("_iprint",  OT_INTEGER, 0);
    addOption("_isumm",   OT_INTEGER, 6);

    // Snopt options
    typedef std::pair<opt_type, std::string> optPair;
    // optionsmap_["_major_print_level"]     = optPair(OT_INTEGER, "Major print level");
    // optionsmap_["_minor_print_level"]     = optPair(OT_INTEGER, "Minor print level");
    optionsmap_["_verify_level"]          = optPair(OT_INTEGER, "Verify level");
    optionsmap_["_major_iteration_limit"] = optPair(OT_INTEGER, "Major iteration limit");
    optionsmap_["_minor_iteration_limit"] = optPair(OT_INTEGER, "Minor iteration limit");
    optionsmap_["_feasibility_tolerance"] = optPair(OT_REAL,    "Feasibility tolerance");
    optionsmap_["_optimality_tolerance"]  = optPair(OT_REAL,    "Optimality tolerance");
    optionsmap_["_scale_option"]          = optPair(OT_INTEGER, "Scale option");

    // Add the Snopt Options
    for (OptionsMap::const_iterator it = optionsmap_.begin(); it != optionsmap_.end(); it++) {
      addOption(it->first, it->second.first, GenericType(), it->second.second);
    }
  }

  SnoptInternal::~SnoptInternal() {
  }

  SnoptInternal* SnoptInternal::clone() const {
    // Use default copy routine
    SnoptInternal* node = new SnoptInternal(*this);

    return node;
  }

  void SnoptInternal::init() {
    // Read in casadi options
    detect_linear_ = getOption("detect_linear");

    // Call the init method of the base class
    NlpSolverInternal::init();

    // Get/generate required functions
    jacF();
    jacG();

    // A large part of what follows is about classiyning variables
    //  and building a mapping

    // Classify the decision variables into (2: nonlinear, 1: linear, 0:zero)
    // according to the objective.
    // x_type_f_:  original variable index -> category w.r.t f
    // x_type_g_:  original variable index -> category w.r.t g
    // x_type_f_:  original constraint index -> category
    x_type_f_.resize(nx_);
    x_type_g_.resize(nx_);
    g_type_.resize(ng_);

    if (detect_linear_) {
      jacF_.spInit(true);
      // Detect dependencies w.r.t. gradF
      // Dependency seeds
      bvec_t* input_v_x =  get_bvec_t(jacF_->inputNoCheck(GRADF_X).data());
      bvec_t* input_v_p =  get_bvec_t(jacF_->inputNoCheck(GRADF_P).data());
      // Make a column with all variables active
      std::fill(input_v_x, input_v_x+nx_, bvec_t(1));
      std::fill(input_v_p, input_v_p+np_, bvec_t(0));
      bvec_t* output_v = get_bvec_t(jacF_->outputNoCheck().data());
      // Perform a single dependency sweep
      jacF_.spEvaluate(true);

      // Harvest the results
      for (int j = 0; j < nx_; ++j) {
        if (jacF_.output().colind(j) == jacF_.output().colind(j+1)) {
          x_type_f_[j] = 0;
        } else {
          x_type_f_[j] = output_v[jacF_.output().colind(j)]?  2 : 1;
        }
      }

      if (!jacG_.isNull()) {  // Detect dependencies w.r.t. jacG
        // Dependency seeds
        bvec_t* input_v_x =  get_bvec_t(jacG_->inputNoCheck(JACG_X).data());
        bvec_t* input_v_p =  get_bvec_t(jacG_->inputNoCheck(JACG_P).data());
        // Make a column with all variables active
        std::fill(input_v_x, input_v_x+nx_, bvec_t(1));
        std::fill(input_v_p, input_v_p+np_, bvec_t(0));
        bvec_t* output_v = get_bvec_t(jacG_->outputNoCheck().data());
        // Perform a single dependency sweep
        jacG_.spEvaluate(true);

        DMatrix out_trans = jacG_.output().T();
        bvec_t* output_v_trans = get_bvec_t(out_trans.data());

        for (int j = 0; j < nx_; ++j) {  // Harvest the results
          if (jacG_.output().colind(j) == jacG_.output().colind(j+1)) {
            x_type_g_[j] = 0;
          } else {
            bool linear = true;
            for (int k = jacG_.output().colind(j); k < jacG_.output().colind(j+1); ++k) {
              linear = linear && !output_v[k];
            }
            x_type_g_[j] = linear? 1 : 2;
          }
        }
        for (int j = 0; j < ng_; ++j) {  // Harvest the results
          if (out_trans.colind(j) == out_trans.colind(j+1)) {
            g_type_[j] = 0;
          } else {
            bool linear = true;
            for (int k = out_trans.colind(j); k < out_trans.colind(j+1); ++k) {
              linear = linear && !output_v_trans[k];
            }
            g_type_[j] = linear? 1 : 2;
          }
        }
      } else {  // Assume all non-linear
        std::fill(x_type_g_.begin(), x_type_g_.end(), 1);
        std::fill(g_type_.begin(), g_type_.end(), 1);
      }

    } else {  // Assume all non-linear variables
      std::fill(x_type_f_.begin(), x_type_f_.end(), 2);
      std::fill(x_type_g_.begin(), x_type_g_.end(), 2);
      std::fill(g_type_.begin(), g_type_.end(), 2);
    }

    if (monitored("setup_nlp")) {
      std::cout << "Variable classification (obj): " << x_type_f_ << std::endl;
      std::cout << "Variable classification (con): " << x_type_g_ << std::endl;
      std::cout << "Constraint classification: " << g_type_ << std::endl;
    }

    // An encoding of the desired sorting pattern
    // Digits xy  with x correspodning to x_type_f_ and y corresponding to x_type_g_
    std::vector<int> order_template;
    order_template.reserve(9);
    order_template.push_back(22);
    order_template.push_back(12);
    order_template.push_back(2);
    order_template.push_back(21);
    order_template.push_back(20);
    order_template.push_back(11);
    order_template.push_back(10);
    order_template.push_back(1);
    order_template.push_back(0);

    // prepare the mapping for variables
    x_order_.resize(0);
    x_order_.reserve(nx_);

    // Populate the order vector
    std::vector<int> x_order_count;  // Cumulative index into x_order_
    for (int p = 0; p < order_template.size(); ++p) {
      for (int k = 0; k < nx_; ++k) {
        if (x_type_f_[k]*10+x_type_g_[k] == order_template[p]) {
          x_order_.push_back(k);
        }
      }
      // Save a cumulative index
      x_order_count.push_back(x_order_.size());
    }

    // prepare the mapping for constraints
    g_order_.resize(0);
    g_order_.reserve(ng_);
    std::vector<int> g_order_count;  // Cumulative index into g_order_
    for (int p = 2; p >= 0; --p) {
      for (int k = 0; k < ng_; ++k) {
        if (g_type_[k] == p) {
          g_order_.push_back(k);
        }
      }
      g_order_count.push_back(g_order_.size());
    }
    nnJac_ = x_order_count[2];
    nnObj_ = x_order_count[4];
    nnCon_ = g_order_count[0];

    if (monitored("setup_nlp")) {
      for (int p = 0; p < order_template.size(); ++p) {
        int start_k = (p > 0 ?x_order_count[p-1]:0);
        std::cout << "Variables (" << order_template[p]/10 << ", "
                  << order_template[p]%10 << ") - " << x_order_count[p]-start_k << ":"
                  << std::vector<int>(x_order_.begin()+start_k,
                                      x_order_.begin()+std::min(x_order_count[p], 200+start_k))
                  << std::endl;
      }

      std::cout << "Variable order:" << x_order_ << std::endl;
      std::cout << "Constraint order:" << g_order_ << std::endl;
      std::cout << "nnJac:" << nnJac_ << std::endl;
      std::cout << "nnObj:" << nnObj_ << std::endl;
      std::cout << "nnCon:" << nnCon_ << std::endl;
    }

    // Here follows the core of the mapping
    //  Two integer matrices are constructed:
    //  one with gradF sparsity, and one with jacG sparsity
    //  the integer values denote the nonzero locations into the original gradF/jacG
    //  but with a special encoding: entries of gradF are encoded "-1-i" and
    //  entries of jacG are encoded "1+i"
    //  "0" is to be interpreted not as an index but as a literal zero

    IMatrix mapping_jacG  = IMatrix::sparse(0, nx_);
    IMatrix mapping_gradF = IMatrix(jacF_.output().sparsity(),
                                    range(-1, -1-jacF_.output().size(), -1));

    if (!jacG_.isNull()) {
      mapping_jacG = IMatrix(jacG_.output().sparsity(), range(1, jacG_.output().size()+1));
    }

    // First, remap jacG
    A_structure_ = mapping_jacG(g_order_, x_order_);

    m_ = ng_;

    // Construct the linear objective row
    IMatrix d = mapping_gradF(Slice(0), x_order_);

    std::vector<int> ii = mapping_gradF.sparsity().getCol();
    for (int j = 0; j < nnObj_; ++j) {
      if (d.colind(j) != d.colind(j+1)) {
        int k = d.colind(j);
        int i = d.data()[k];  // Nonzero original index into gradF
        if (x_type_f_[ii[-1-i]] == 2) {
          d[k] = 0;
        }
      }
    }

    // Make it as sparse as you can
    d = sparse(d);

    jacF_row_ = d.size() != 0;
    if (jacF_row_) {  // We need an objective gradient row
      A_structure_ = vertcat(A_structure_, d);
      m_ +=1;
    }
    iObj_ = jacF_row_ ? m_ : 0;

    // Is the A matrix completely empty?
    dummyrow_ = A_structure_.size() == 0;  // Then we need a dummy row
    if (dummyrow_) {
      IMatrix dummyrow = IMatrix::sparse(1, nx_);
      dummyrow(0, 0) = 0;
      A_structure_ = vertcat(A_structure_, dummyrow);
      m_+=1;
    }

    // We don't need a dummy row if a linear objective row is present
    casadi_assert(!(dummyrow_ && jacF_row_));

    if (monitored("setup_nlp")) {
      std::cout << "Objective gradient row presence: " << jacF_row_ << std::endl;
      std::cout << "Dummy row presence: " << dummyrow_ << std::endl;
      std::cout << "iObj: " << iObj_ << std::endl;
    }

    // The following section deals with querying SNOPT upfront
    // how much memory it needs

    int iPrint = 0;
    int iSumm = 0;

    // Manual says 500 is a bare minimum, which we will correct later
    int clen = 500; int rlen = 500; int ilen = 500;
    snopt_cw_.resize(clen*8, 'a');
    snopt_iw_.resize(ilen, 0);
    snopt_rw_.resize(rlen, 0);

    snInit(iPrint, iSumm);

    int mincw = 0;
    int miniw = 0;
    int minrw = 0;
    int INFO = 0;
    // Number of nonzeros in entire A
    int neA = A_structure_.size();
    // Number of nonzeros in nonlinear part of A
    int negCon = std::max(A_structure_(Slice(0, nnCon_), Slice(0, nnJac_)).size(), 1);
    snMemb(&INFO, &m_, &nx_, &neA, &negCon, &nnCon_, &nnJac_, &nnObj_, &mincw, &miniw, &minrw);

    casadi_assert(INFO == 104);

    snopt_cw_.resize(mincw*8, 'a');
    snopt_iw_.resize(miniw, 0);
    snopt_rw_.resize(minrw, 0);
    clen = mincw;
    ilen = miniw;
    rlen = minrw;

    iPrint = getOption("_iprint");
    iSumm = getOption("_isumm");;
    snInit(iPrint, iSumm);

    for (OptionsMap::const_iterator it = optionsmap_.begin(); it != optionsmap_.end(); it++) {
      int Error = 0;
      const std::string & snopt_name = it->second.second;
      //int bufferlen = snopt_name.size();
      if (hasSetOption(it->first)) {
        switch (it->second.first) {
          case OT_INTEGER: {
            int value = getOption(it->first);
            snSeti(snopt_name, value);
          }; break;
          case OT_REAL: {
            double value = getOption(it->first);
            snSetr(snopt_name, value);
          }; break;
          case OT_STRING: {
            std::string value = getOption(it->first);
            snSet(snopt_name, value);
            }; break;
          default:
            casadi_error("Unkown type " << it->second.first);
        }
        casadi_assert_message(Error == 0, "snopt error setting option \"" + snopt_name + "\"")
      }
      // else {
      //   switch (it->second.first) {
      //     case OT_INTEGER: {
      //       int value = 0;
      //       snopt_geti(snopt_name.c_str(), &bufferlen, &value, &Error,
      //                  getPtr(snopt_cw_), &clen,
      //                  getPtr(snopt_iw_), &ilen,
      //                  getPtr(snopt_rw_), &rlen);
      //       setOption(it->first, value);
      //     }; break;
      //     case OT_REAL: {
      //       double value = 0;
      //       snopt_getr(snopt_name.c_str(), &bufferlen, &value, &Error,
      //                  getPtr(snopt_cw_), &clen,
      //                  getPtr(snopt_iw_), &ilen,
      //                  getPtr(snopt_rw_), &rlen);
      //       setOption(it->first, value);
      //     }; break;
      //     case OT_STRING: {
      //       char value[8];
      //       snopt_getc(snopt_name.c_str(), &bufferlen, value, &Error,
      //                  getPtr(snopt_cw_), &clen,
      //                  getPtr(snopt_iw_), &ilen,
      //                  getPtr(snopt_rw_), &rlen);
      //       setOption(it->first, std::string(value));
      //       }; break;
      //     default:
      //       casadi_error("Unkown type " << it->second.first);
      //   }
      //   casadi_assert_message(Error == 0, "snopt error getting option \"" + snopt_name + "\"")
      // }
    }

    // Allocate data structures needed in evaluate
    bl_.resize(nx_+m_);
    bu_.resize(nx_+m_);
    hs_.resize(nx_+m_);
    x_.resize(nx_+m_);
    pi_.resize(m_);
    rc_.resize(nx_+m_);
    A_data_.resize(A_structure_.size());

    // Reset the counters
    t_eval_grad_f_ = t_eval_jac_g_ = t_callback_fun_ = t_mainloop_ = 0;
    n_eval_grad_f_ = n_eval_jac_g_ = n_callback_fun_ = n_iter_ = 0;
  }

  void SnoptInternal::reset() {
  }

  void SnoptInternal::setQPOptions() {
  }

  void SnoptInternal::passOptions() {
  }

  std::string SnoptInternal::formatStatus(int status) const {
    if (status_.find(status) == status_.end()) {
      std::stringstream ss;
      ss << "Unknown status: " << status;
      return ss.str();
    } else {
      return (*status_.find(status)).second;
    }
  }

  void SnoptInternal::evaluate() {
    log("SnoptInternal::evaluate");

    // Initial checks
    if (inputs_check_) checkInputs();
    checkInitialBounds();

    if (gather_stats_) {
      Dictionary iterations;
      iterations["inf_pr"] = std::vector<double>();
      iterations["inf_du"] = std::vector<double>();
      iterations["merit"] = std::vector<double>();
      iterations["step_size"] = std::vector<double>();
      iterations["pen_norm"] = std::vector<double>();
      iterations["cond_H"] = std::vector<double>();
      iterations["qp_num_iter"] = std::vector<int>();
      stats_["iterations"] = iterations;
    }

    std::string start = getOption("_start");

    // Evaluate gradF and jacG at initial value
    if (!jacG_.isNull()) {
      jacG_.setInput(input(NLP_SOLVER_X0), JACG_X);
      jacG_.setInput(input(NLP_SOLVER_P), JACG_P);
      jacG_.evaluate();
    }
    jacF_.setInput(input(NLP_SOLVER_X0), GRADF_X);
    jacF_.setInput(input(NLP_SOLVER_P), GRADF_P);

    jacF_.evaluate();

    // perform the mapping:
    // populate A_data_ (the nonzeros of A)
    // with numbers pulled from jacG and gradF
    for (int k = 0; k < A_structure_.size(); ++k) {
      int i = A_structure_.data()[k];
      if (i == 0) {
        A_data_[k] = 0;
      } else if (i > 0) {
        A_data_[k] = jacG_.output().data()[i-1];
      } else {
        A_data_[k] = jacF_.output().data()[-i-1];
      }
    }

    // Obtain constraint offsets for linear constraints
    if (detect_linear_) {
      nlp_.setInput(0.0, NL_X);
      // Setting the zero might actually be problematic
      nlp_.setInput(input(NLP_SOLVER_P), NL_P);
      nlp_.evaluate();
    }

    int n = nx_;
    int nea = A_structure_.size();
    double ObjAdd = 0;

    // Problem name
    std::string prob = "  CasADi";

    // Obtain sparsity pattern of A (fortan is Index-1 based)
    std::vector<int> row(A_structure_.size());
    for (int k = 0; k < A_structure_.size(); ++k) {
      row[k] = 1+A_structure_.row()[k];
    }

    std::vector<int> col(nx_+1);
    for (int k = 0; k < nx_+1; ++k) {
      col[k] = 1+A_structure_.colind()[k];
    }

    // Obtain initial guess and bounds through the mapping
    for (int k = 0; k < nx_; ++k) {
      int kk = x_order_[k];
      bl_[k] = input(NLP_SOLVER_LBX).data()[kk];
      bu_[k] = input(NLP_SOLVER_UBX).data()[kk];
      x_[k] = input(NLP_SOLVER_X0).data()[kk];
    }
    for (int k = 0; k < ng_; ++k) {
      int kk = g_order_[k];
      if (g_type_[kk] < 2) {
        // casadi_error("woops");
        bl_[nx_+k] = input(NLP_SOLVER_LBG).data()[kk]-nlp_.output(NL_G).data()[kk];
        bu_[nx_+k] = input(NLP_SOLVER_UBG).data()[kk]-nlp_.output(NL_G).data()[kk];
      } else {
        bl_[nx_+k] = input(NLP_SOLVER_LBG).data()[kk];
        bu_[nx_+k] = input(NLP_SOLVER_UBG).data()[kk];
      }
      x_[nx_+k] = input(NLP_SOLVER_LAM_G0).data()[kk];
    }

    // Objective row / dummy row should be unbounded
    if (dummyrow_ || jacF_row_) {
      bl_.back() = -1e20;  // -std::numeric_limits<double>::infinity();
      bu_.back() = 1e20;  // std::numeric_limits<double>::infinity();
    }

    int nS = 0;
    int info = 0;

    // Outputs
    int miniw = 0, minrw = 0, mincw = 0;
    int nInf = 0;
    double Obj = 0;
    double sInf = 0;

    casadi_assert(m_ > 0);
    casadi_assert(n > 0);
    casadi_assert(nea > 0);
    casadi_assert(row.size() == nea);
    casadi_assert(hs_.size() == n+m_);
    casadi_assert(col.size() == n+1);
    casadi_assert(A_structure_.size() == nea);
    casadi_assert(bl_.size() == n+m_);
    casadi_assert(bu_.size() == n+m_);
    casadi_assert(pi_.size() == m_);
    casadi_assert(x_.size() == n+m_);
    casadi_assert(col.at(0) == 1);
    casadi_assert(col.at(n) == nea+1);

    // Pointer magic, courtesy of Greg
    int iulen = 8;
    std::vector<int> iu(iulen);
    SnoptInternal* source = this;
    memcpy(&(iu[0]), &source, sizeof(SnoptInternal*));

    casadi_assert_message(!jacF_.isNull(), "blaasssshc");

    if (monitored("setup_nlp")) {
      std::cout << "indA:" << row << std::endl;
      std::cout << "locA:" << col << std::endl;
      std::cout << "colA:" << A_data_ << std::endl;
      A_structure_.sparsity().spy();
      std::cout << "A:" << DMatrix(A_structure_.sparsity(), A_data_) << std::endl;
      std::cout << "n:" << n << std::endl;
      std::cout << "m:" << m_ << std::endl;
      std::cout << "nea:" << nea << std::endl;
      std::cout << "bl_:" << bl_ << std::endl;
      std::cout << "bu_:" << bu_ << std::endl;
    }

    // Run SNOPT
    double time0 = clock();
    snoptC(start, &m_, &n, &nea,
           &nnCon_, &nnObj_, &nnJac_, &iObj_, &ObjAdd,
           prob, userfunPtr, snStopPtr,
           A_data_, row, col, bl_, bu_,
           // Initial values
           getPtr(hs_), getPtr(x_), getPtr(pi_), getPtr(rc_),
           // Outputs
           &info, &mincw, &miniw, &minrw, &nS, &nInf, &sInf, &Obj,
           iu);
    t_mainloop_ = static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;

    stats_["return_status"] = info;

    // Store results into output
    for (int k = 0; k < nx_; ++k) {
      int kk =  x_order_[k];
      output(NLP_SOLVER_X).data()[kk] = x_[k];
      output(NLP_SOLVER_LAM_X).data()[kk] = -rc_[k];
    }

    setOutput(Obj+ (jacF_row_? x_[nx_+ng_] : 0), NLP_SOLVER_F);

    for (int k = 0; k < ng_; ++k) {
      int kk = g_order_[k];
      output(NLP_SOLVER_LAM_G).data()[kk] = -rc_[nx_+k];
      output(NLP_SOLVER_G).data()[kk] = x_[nx_+k];  // TODO(Joris): this is not quite right
      // mul_no_alloc
    }


    // print timing information
    // save state
    std::ios state(NULL);
    state.copyfmt(std::cout);
    const int w_time = 7;
    const int p_time = 3;
    const int w_ms = 7;
    const int p_ms = 2;
    const int w_n = 5;
    if (hasOption("print_time") && static_cast<bool>(getOption("print_time"))) {
      std::cout << std::endl;

      std::cout << "time spent in eval_grad_f       "
                << std::fixed << std::setw(w_time) << std::setprecision(p_time)
                << t_eval_grad_f_ << " s.";
      if (n_eval_grad_f_>0)
        std::cout << " (" << std::setw(w_n) << n_eval_grad_f_ << " calls, "
                  << std::setw(w_ms) << std::setprecision(p_ms)
                  << (t_eval_grad_f_/n_eval_grad_f_)*1000 << " ms average)";
      std::cout << std::endl;

      std::cout << "time spent in eval_jac_g        "
                << std::fixed << std::setw(w_time) << std::setprecision(p_time)
                << t_eval_jac_g_ << " s.";
      if (n_eval_jac_g_>0)
        std::cout << " (" << std::setw(w_n) << n_eval_jac_g_ << " calls, "
                  << std::setw(w_ms) << std::setprecision(p_ms)
                  << (t_eval_jac_g_/n_eval_jac_g_)*1000 << " ms average)";
      std::cout << std::endl;

      std::cout << "time spent in callback function "
                << std::fixed << std::setw(w_time) << std::setprecision(p_time)
                << t_callback_fun_ << " s.";
      if (n_callback_fun_>0)
        std::cout << " (" << std::setw(w_n) << n_callback_fun_ << " calls, "
                  << std::setw(w_ms) << std::setprecision(p_ms)
                  << (t_callback_fun_/n_callback_fun_)*1000 << " ms average)";
      std::cout << std::endl;

      std::cout << "time spent in main loop         "
                << std::setw(w_time) << std::setprecision(p_time)
                << t_mainloop_ << " s." << std::endl;
    }
    // restore state
    std::cout.copyfmt(state);

    // set timing information
    stats_["t_eval_grad_f"] = t_eval_grad_f_;
    stats_["t_eval_jac_g"] = t_eval_jac_g_;
    stats_["t_mainloop"] = t_mainloop_;
    stats_["t_callback_fun"] = t_callback_fun_;
    stats_["n_eval_grad_f"] = n_eval_grad_f_;
    stats_["n_eval_jac_g"] = n_eval_jac_g_;
    stats_["n_callback_fun"] = n_callback_fun_;
    stats_["iter_count"] = n_iter_-1;

    // Reset the counters
    t_eval_grad_f_ = t_eval_jac_g_ = t_callback_fun_ = t_mainloop_ = 0;
    n_eval_grad_f_ = n_eval_jac_g_ = n_callback_fun_ = n_iter_ = 0;
  }

  void SnoptInternal::setOptionsFromFile(const std::string & file) {
  }

  void SnoptInternal::userfun(
      int* mode, int nnObj, int nnCon, int nnJac, int nnL, int neJac,
      double* x, double* fObj, double*gObj, double* fCon, double* gCon,
      int nState, char* cu, int lencu, int* iu, int leniu, double* ru, int lenru) {
    try {
      double time0 = clock();

      casadi_assert_message(nnCon_ == nnCon, "Con " << nnCon_ << " <-> " << nnCon);
      casadi_assert_message(nnObj_ == nnObj, "Obj " << nnObj_ << " <-> " << nnObj);
      casadi_assert_message(nnJac_ == nnJac, "Jac " << nnJac_ << " <-> " << nnJac);

      // Evaluate gradF with the linear variables put to zero
      jacF_.setInput(0.0, NL_X);
      jacF_.setInput(input(NLP_SOLVER_P), NL_P);
      for (int k = 0; k < nnObj; ++k) {
        if (x_type_f_[x_order_[k]] == 2) {
          jacF_.input(NL_X)[x_order_[k]] = x[k];
        }
      }

      if (monitored("eval_nlp")) {
        std::cout << "mode: " << *mode << std::endl;
        std::cout << "A before we touch it:"
                  << DMatrix(A_structure_.sparsity(), A_data_) << std::endl;
        std::cout << "x (obj - sorted indices   - all elements present):"
                  << std::vector<double>(x, x+nnObj) << std::endl;
        std::cout << "x (obj - original indices - linear elements zero):"
                  << jacF_.input(NL_X) << std::endl;
      }

      jacF_.evaluate();

      // provide objective (without linear contributions) to SNOPT
      *fObj = jacF_.output(1).at(0);

      // provide nonlinear part of objective gradient to SNOPT
      for (int k = 0; k < nnObj; ++k) {
        int i = x_order_[k];
        if (x_type_f_[i] == 2) {
          int el = jacF_.output().colind(i);
          if (jacF_.output().colind(i+1) > el) {
            gObj[k] = jacF_.output().data()[el];
          } else {
            gObj[k] = 0;
          }
        } else {
          gObj[k] = 0;
        }
      }

      jacF_.output().sparsity().sanityCheck(true);
      jacF_.output().sparsity().sanityCheck(false);

      // timing and counters
      t_eval_grad_f_ += static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;
      n_eval_grad_f_ += 1;


      if (monitored("eval_nlp")) {
        std::cout << "fObj:" << *fObj << std::endl;
        std::cout << "gradF:" << jacF_.output() << std::endl;
        std::cout << "gObj:" << std::vector<double>(gObj, gObj+nnObj) << std::endl;
      }

      time0 = clock();
      if (!jacG_.isNull()) {
        // Evaluate jacG with the linear variabes put to zero
        jacG_.setInput(0.0, JACG_X);
        for (int k = 0; k < nnJac; ++k) {
          jacG_.input(JACG_X)[x_order_[k]] = x[k];
        }
        if (monitored("eval_nlp")) {
          std::cout << "x (con - sorted indices   - all elements present):"
                    << std::vector<double>(x, x+nnJac) << std::endl;
          std::cout << "x (con - original indices - linear elements zero):"
                    << jacG_.input(JACG_X) << std::endl;
        }
        jacG_.setInput(input(NLP_SOLVER_P), JACG_P);

        // Evaluate the function
        jacG_.evaluate();

        // provide nonlinear part of constraint jacobian to SNOPT
        int kk = 0;
        for (int j = 0; j < nnJac; ++j) {
          for (int k = A_structure_.colind(j); k < A_structure_.sparsity().colind(j+1); ++k) {
            if (A_structure_.row(k) >= nnCon) break;
            int i = A_structure_.data()[k];
            if (i > 0) {
              gCon[kk++] = jacG_.output().data()[i-1];
            }
          }
        }

        casadi_assert(kk == 0 || kk == neJac);

        if (monitored("eval_nlp")) {
          std::cout << jacG_.output(GRADF_G) << std::endl;
        }

        // provide nonlinear part of objective to SNOPT
        DMatrix g = jacG_.output();
        for (int k = 0; k < nnCon; ++k) {
          fCon[k] = jacG_.output(GRADF_G).data()[g_order_[k]];
        }

        // timing and counters
        t_eval_jac_g_ += static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;
        n_eval_jac_g_ += 1;

        if (monitored("eval_nlp")) {
          std::cout << "fCon:" << std::vector<double>(fCon, fCon+nnCon) << std::endl;
          std::cout << "gCon:" << std::vector<double>(gCon, gCon+neJac) << std::endl;
        }
      }

    } catch(std::exception& ex) {
      std::cerr << "eval_nlp failed: " << ex.what() << std::endl;
      *mode = -1;  // Reduce step size - we've got problems
      return;
    }
  }

  void SnoptInternal::callback(
      int* iAbort, int* info, int HQNType, int* KTcond, int MjrPrt, int minimz,
      int m, int maxS, int n, int nb, int nnCon0, int nnCon, int nnObj0, int nnObj, int nS,
      int itn, int nMajor, int nMinor, int nSwap,
      double condHz, int iObj, double sclObj, double ObjAdd,
      double fMrt,  double PenNrm,  double step,
      double prInf,  double duInf,  double vimax,  double virel, int* hs,
      int ne, int nlocJ, int* locJ, int* indJ, double* Jcol, int negCon,
      double* Ascale, double* bl, double* bu, double* fCon, double* gCon, double* gObj,
      double* yCon, double* pi, double* rc, double* rg, double* x,
      double*  cu, int lencu, int* iu, int leniu, double* ru, int lenru,
      char*   cw, int lencw,  int* iw, int leniw, double* rw, int lenrw) {
    try {
      n_iter_+=1;
      if (gather_stats_) {
        Dictionary & iterations = stats_["iterations"];
        static_cast<std::vector<double> &>(iterations["inf_pr"]).push_back(prInf);
        static_cast<std::vector<double> &>(iterations["inf_du"]).push_back(duInf);
        static_cast<std::vector<double> &>(iterations["merit"]).push_back(fMrt);
        static_cast<std::vector<double> &>(iterations["step_size"]).push_back(step);
        static_cast<std::vector<double> &>(iterations["pen_norm"]).push_back(PenNrm);
        static_cast<std::vector<double> &>(iterations["cond_H"]).push_back(condHz);
        static_cast<std::vector<int> &>(
        iterations["qp_num_iter"]).push_back(nMinor);
      }
      if (!callback_.isNull()) {
        double time0 = clock();
        for (int k = 0; k < nx_; ++k) {
          int kk = x_order_[k];
          output(NLP_SOLVER_X).data()[kk] = x_[k];
          // output(NLP_SOLVER_LAM_X).data()[kk] = -rc_[k];
        }

        // setOutput(Obj+ (jacF_row_? x_[nx_+ng_] : 0), NLP_SOLVER_F);
        for (int k = 0; k < ng_; ++k) {
          int kk =  g_order_[k];
          // output(NLP_SOLVER_LAM_G).data()[kk] = -rc_[nx_+k];
          output(NLP_SOLVER_G).data()[kk] = x_[nx_+k];
        }

        *iAbort = callback_(ref_, user_data_);
        t_callback_fun_ += static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;
        n_callback_fun_ += 1;
      }
    } catch(std::exception& ex) {
      if (getOption("iteration_callback_ignore_errors")) {
        std::cerr << "callback: " << ex.what() << std::endl;
      } else {
        throw ex;
      }
    }
  }

  void SnoptInternal::userfunPtr(
      int * mode, int* nnObj, int * nnCon, int *nJac, int *nnL, int * neJac,
      double *x, double *fObj, double *gObj, double * fCon, double* gCon,
      int* nState,
      char* cu, int* lencu, int* iu, int* leniu, double* ru, int *lenru) {
    SnoptInternal* interface;  // = reinterpret_cast<SnoptInternal*>(iu);
    memcpy(&interface, &(iu[0]), sizeof(SnoptInternal*));

    interface->userfun(mode, *nnObj, *nnCon, *nJac, *nnL, *neJac,
                       x, fObj, gObj, fCon, gCon, *nState,
                       cu, *lencu, iu, *leniu, ru, *lenru);
  }

  void SnoptInternal::snStopPtr(
    int* iAbort, int* info, int* HQNType, int* KTcond, int* MjrPrt, int* minimz,
    int* m, int* maxS, int* n, int* nb,
    int* nnCon0, int* nnCon, int* nnObj0, int* nnObj, int* nS,
    int* itn, int* nMajor, int* nMinor, int* nSwap,
    double * condHz, int* iObj, double * sclObj,  double *ObjAdd,
    double * fMrt,  double * PenNrm,  double * step,
    double *prInf,  double *duInf,  double *vimax,  double *virel, int* hs,
    int* ne, int* nlocJ, int* locJ, int* indJ, double* Jcol, int* negCon,
    double* Ascale, double* bl, double* bu, double* fCon, double* gCon, double* gObj,
    double* yCon, double* pi, double* rc, double* rg, double* x,
    double*  cu, int * lencu, int* iu, int* leniu, double* ru, int *lenru,
    char*   cw, int* lencw,  int* iw, int *leniw, double* rw, int* lenrw) {
    SnoptInternal* interface;  // = reinterpret_cast<SnoptInternal*>(iu);
    memcpy(&interface, &(iu[0]), sizeof(SnoptInternal*));

    interface->callback(iAbort, info, *HQNType, KTcond, *MjrPrt, *minimz,
    *m, *maxS, *n, *nb, *nnCon0, *nnCon, *nnObj0, *nnObj, *nS,
    *itn, *nMajor, *nMinor, *nSwap,
    *condHz, *iObj, *sclObj, *ObjAdd,  *fMrt,  *PenNrm,  *step,
    *prInf,  *duInf, *vimax, *virel, hs,
    *ne, *nlocJ, locJ, indJ, Jcol, *negCon,
    Ascale, bl, bu, fCon, gCon, gObj,
    yCon, pi, rc, rg, x,
     cu, *lencu, iu, *leniu, ru, *lenru,
    cw, *lencw,  iw, *leniw, rw, *lenrw);
  }

  void SnoptInternal::snMemb(int *INFO, const int *m_, const int *nx_,
                             const int *neA, const int * negCon, const int * nnCon_,
                             const int *nnJac_, const int *nnObj_,
                             int *mincw, int *miniw, int *minrw) {
      int clen = snopt_cw_.size()/8;
      int ilen = snopt_iw_.size();
      int rlen = snopt_rw_.size();
      snmemb_(INFO, m_, nx_, neA, negCon, nnCon_, nnJac_, nnObj_,
              mincw, miniw, minrw,
              getPtr(snopt_cw_), &clen,
              getPtr(snopt_iw_), &ilen,
              getPtr(snopt_rw_), &rlen,
              clen*8);
      //snopt_memb(INFO, m_, nx_, neA, negCon, nnCon_, nnJac_, nnObj_,
      //           mincw, miniw, minrw,
      //           getPtr(snopt_cw_), &clen,
      //           getPtr(snopt_iw_), &ilen,
      //           getPtr(snopt_rw_), &rlen);
  }
  void SnoptInternal::snInit(int iPrint, int iSumm) {
      int clen = snopt_cw_.size()/8;
      int ilen = snopt_iw_.size();
      int rlen = snopt_rw_.size();
      sninit_(&iPrint, &iSumm,
              getPtr(snopt_cw_), &clen,
              getPtr(snopt_iw_), &ilen,
              getPtr(snopt_rw_), &rlen,
              clen*8);
      //snopt_init(&iPrint, &iSumm,
      //           getPtr(snopt_cw_), &clen,
      //           getPtr(snopt_iw_), &ilen,
      //           getPtr(snopt_rw_), &rlen);
  }

  void SnoptInternal::snSeti(const std::string &snopt_name, int value) {
      int bufferlen = snopt_name.size();
      int clen = snopt_cw_.size()/8;
      int ilen = snopt_iw_.size();
      int rlen = snopt_rw_.size();
      int iSumm = getOption("_isumm");
      int iPrint = getOption("_iprint");
      int Error = 0;
      snseti_(snopt_name.c_str(), &value, &iPrint, &iSumm, &Error,
              getPtr(snopt_cw_), &clen,
              getPtr(snopt_iw_), &ilen,
              getPtr(snopt_rw_), &rlen,
              bufferlen, clen*8);
      //snopt_seti(snopt_name.c_str(), &bufferlen, &value, &iPrint, &iSumm, &Error,
      //           getPtr(snopt_cw_), &clen,
      //           getPtr(snopt_iw_), &ilen,
      //           getPtr(snopt_rw_), &rlen);
      casadi_assert_message(Error == 0, "snopt error setting option \"" + snopt_name + "\"")
  }

  void SnoptInternal::snSetr(const std::string &snopt_name, double value) {
      int bufferlen = snopt_name.size();
      int clen = snopt_cw_.size()/8;
      int ilen = snopt_iw_.size();
      int rlen = snopt_rw_.size();
      int iSumm = getOption("_isumm");
      int iPrint = getOption("_iprint");
      int Error = 0;
      snsetr_(snopt_name.c_str(), &value, &iPrint, &iSumm, &Error,
              getPtr(snopt_cw_), &clen,
              getPtr(snopt_iw_), &ilen,
              getPtr(snopt_rw_), &rlen,
              bufferlen, clen*8);
      //snopt_setr(snopt_name.c_str(), &bufferlen, &value, &iPrint, &iSumm, &Error,
      //           getPtr(snopt_cw_), &clen,
      //           getPtr(snopt_iw_), &ilen,
      //           getPtr(snopt_rw_), &rlen);
      casadi_assert_message(Error == 0, "snopt error setting option \"" + snopt_name + "\"")
  }

  void SnoptInternal::snSet(const std::string &snopt_name, const std::string &value0) {
      std::string value = value0;
      assert(value.size() <= 8);
      value.append(8-value.size(), ' ');
      std::string buffer = snopt_name;
      buffer.append(" = ");
      buffer.append(value);
      int bufferlen = buffer.size();

      int clen = snopt_cw_.size()/8;
      int ilen = snopt_iw_.size();
      int rlen = snopt_rw_.size();
      int iSumm = getOption("_isumm");
      int iPrint = getOption("_iprint");
      int Error = 0;

      snset_(buffer.c_str(), &iPrint, &iSumm, &Error,
             getPtr(snopt_cw_), &clen,
             getPtr(snopt_iw_), &ilen,
             getPtr(snopt_rw_), &rlen,
             bufferlen, clen*8);
      //snopt_set(buffer.c_str(), &bufferlen, &iPrint, &iSumm, &Error,
      //          getPtr(snopt_cw_), &clen,
      //          getPtr(snopt_iw_), &ilen,
      //          getPtr(snopt_rw_), &rlen);

      casadi_assert_message(Error == 0, "snopt error setting option \"" + snopt_name + "\"")
  }

  void SnoptInternal::snoptC(
         const std::string & start, const int * m_, const int * n, const int * neA,
         const int *nnCon, const int *nnObj, const int *nnJac,
         const int *iObj, const double *ObjAdd,
         const std::string & prob, UserFun userfunPtr, snStop snStopPtr,
         const std::vector<double>& A_data_, const std::vector<int>& row,
         const std::vector<int>& col,
         std::vector<double>& bl_, std::vector<double>& bu_,
         // Initial values
         int* hs, double* x, double* pi, double * rc,
         // Outputs
         int *info, int* mincw, int* miniw, int* minrw, int * nS,
         int* nInf, double* sInf, double* Obj,
         std::vector<int>&iu) {
      const int lenstart = start.size();
      const int clen = snopt_cw_.size()/8;
      const int rlen = snopt_rw_.size();
      const int ilen = snopt_iw_.size();

      const int iulen = iu.size();
      const int nName = 1;

      const int prob_len = prob.size();
      casadi_assert_message(prob_len==8, "snopt bindings: problem name must be 8 chars, your name: "
                            "\"" << prob << "\"");
//      snopt_c(
//        start.c_str(), &lenstart, m_, n, neA, &nName,
//        nnCon, nnObj, nnJac, iObj, ObjAdd,
//        prob.c_str(), userfunPtr,
//        getPtr(A_data_), getPtr(row), getPtr(col), getPtr(bl_), getPtr(bu_),
//        0,
//        // Initial values
//        getPtr(hs_), getPtr(x_), getPtr(pi_), getPtr(rc_),
//        // Outputs
//        info, mincw, miniw, minrw, nS, nInf, sInf, Obj,
//        // Working spaces for usrfun
//        getPtr(snopt_cw_), &clen,
//        getPtr(iu), &iulen,
//        getPtr(snopt_rw_), &rlen,
//        // Working spaces for SNOPT
//        getPtr(snopt_cw_), &clen,
//        getPtr(snopt_iw_), &ilen,
//        getPtr(snopt_rw_), &rlen);
      snkerc_(
          start.c_str(), m_, n, neA, &nName,
          nnCon, nnObj, nnJac, iObj, ObjAdd,
          prob.c_str(), userfunPtr, snlog_, snlog2_, sqlog_, snStopPtr,
          getPtr(A_data_), getPtr(row), getPtr(col), getPtr(bl_), getPtr(bu_),
          0,
          // Initial values
          getPtr(hs_), getPtr(x_), getPtr(pi_), getPtr(rc_),
          // Outputs
          info, mincw, miniw, minrw, nS, nInf, sInf, Obj,
          // Working spaces for usrfun
          getPtr(snopt_cw_), &clen,
          getPtr(iu), &iulen,
          getPtr(snopt_rw_), &rlen,
          // Working spaces for SNOPT
          getPtr(snopt_cw_), &clen,
          getPtr(snopt_iw_), &ilen,
          getPtr(snopt_rw_), &rlen,
          lenstart, prob_len, 0,
          clen*8, clen*8);
  }

}  // namespace casadi
