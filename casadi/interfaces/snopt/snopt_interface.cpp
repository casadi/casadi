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

#include "snopt_interface.hpp"
#include "casadi/core/std_vector_tools.hpp"

#include <stdio.h>
#include <string.h>
#include <ctime>
#include <utility>
#include <algorithm>
#include <iomanip>

namespace casadi {

  extern "C"
  int CASADI_NLPSOL_SNOPT_EXPORT
  casadi_register_nlpsol_snopt(Nlpsol::Plugin* plugin) {
    plugin->creator = SnoptInterface::creator;
    plugin->name = "snopt";
    plugin->doc = SnoptInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_SNOPT_EXPORT casadi_load_nlpsol_snopt() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_snopt);
  }

  SnoptInterface::SnoptInterface(const std::string& name, const XProblem& nlp)
    : Nlpsol(name, nlp) {

    addOption("detect_linear", OT_BOOLEAN, true,
              "Make an effort to treat linear constraints and linear variables specially.");

    // casadi options
    addOption("print_time", OT_BOOLEAN, true, "print information about execution time");

    // options which are handled seperately
    addOption("start",  OT_STRING, "Cold",  "", "Cold|Basis|Warm");
    addOption("print file",  OT_STRING);
    addOption("specs file",  OT_STRING);
    addOption("summary", OT_BOOLEAN, true);

    // Snopt options
    typedef std::pair<std::string, std::string> spair;

    // Printing
    intOpts_["Major print level"] = "1 * 1-line major iteration log";
    intOpts_["Minor print level"] = "1 * 1-line minor iteration log";
    // special om["Print file"] = OT_S; //  * specified by subroutine sn_init
    // special om["Summary file"] = OT_S; //  * specified by subroutine sn_init
    intOpts_["Print frequency"] = "100 * minor iterations log on Print file";
    intOpts_["Summary frequency"] = "100 * minor iterations log on Summary file";
    strOpts_["Solution"] = spair("Yes|No|If Optimal|If Infeasible|If Unbounded",
                                 "Yes * on the Print file");

    // * Suppress options listing * options are normally listed
    strOpts_["System information"] =
      spair("No|Yes", "No * Yes prints more system information");

    // * Problem specification
    // special Minimize * (opposite of Maximize)

    // * Feasible point * (alternative to Max or Min)
    // Objective row 1 * has precedence over ObjRow (snOptA)
    // Infinite bound 1.0e+20 *

    // * Convergence Tolerances
    realOpts_["Major feasibility tolerance"] = "1.0e-6 * target nonlinear constraint violation";
    realOpts_["Major optimality tolerance"] = "1.0e-6 * target complementarity gap";
    realOpts_["Minor feasibility tolerance"] = "1.0e-6 * for satisfying the QP bounds";

    // * Derivative checking
    intOpts_["Verify level"] = "0 * cheap check on gradients";
    // Start objective check at col 1 * NOT ALLOWED IN snOptA
    // Stop objective check at col n'1 * NOT ALLOWED IN snOptA
    // Start constraint check at col 1 * NOT ALLOWED IN snOptA
    // Stop constraint check at col n''1 * NOT ALLOWED IN snOptA

    // * Scaling
    intOpts_["Scale option"] = "1 * linear constraints and variables";
    realOpts_["Scale tolerance"] = "0.9";
    // SPECIAL // * Scale Print * default: scales are not printed

    // * Other Tolerances
    realOpts_["Crash tolerance"] = "0.1";
    realOpts_["Linesearch tolerance"] = "0.9 * smaller for more accurate search";
    realOpts_["Pivot tolerance"] = "3.7e-11 * e^2/3";

    // * QP subproblems
    strOpts_["QPSolver"] = spair("Cholesky|CG|QN", "Cholesky * default");
    intOpts_["Crash option"] = "3 * first basis is essentially triangular";
    realOpts_["Elastic weight"] = "1.0e+4 * used only during elastic mode";
    intOpts_["Iterations limit"] = "10000 * or 20m if that is more";
    intOpts_["Partial price"] = "1 * 10 for large LPs";

    //* SQP method
    // special * Cold start * has precedence over argument start
    // special * Warm start * (alternative to a cold start)
    intOpts_["Major iterations limit"] = "1000 * or m if that is more";
    intOpts_["Minor iterations limit"] = "500 * or 3m if that is more";
    realOpts_["Major step limit"] = "2.0";
    intOpts_["Superbasics limit"] = "n1 + 1 * n1 = number of nonlinear variables";
    intOpts_["Derivative level"] = "3";
    // Derivative option 1 * ONLY FOR snOptA
    //?????    om["Derivative linesearch"] *

    // * Nonderivative linesearch *
    realOpts_["Function precision"] = "3.0e-13 * e^0.8 (almost full accuracy)";
    realOpts_["Difference interval"] = "5.5e-7 * (Function precision)^1/2";
    realOpts_["Central difference interval"] = "6.7e-5 * (Function precision)^1/3";
    intOpts_["New superbasics limit"] = "99 * controls early termination of QPs";
    // Objective row ObjRow * row number of objective in F(x)
    realOpts_["Penalty parameter"] = "0.0 * initial penalty parameter";
    intOpts_["Proximal point method"] = "1 * satisfies linear constraints near x0";
    intOpts_["Reduced Hessian dimension"] = "2000 * or Superbasics limit if that is less";
    realOpts_["Violation limit"] = "10.0 * unscaled constraint violation limit";
    realOpts_["Unbounded step size"] = "1.0e+18";
    realOpts_["Unbounded objective"] = "1.0e+15";

    // * Hessian approximation
    strOpts_["Hessian"] = spair("full memory|limited memory",
                                "   full memory * default if n1 ≤ 75\n"
                                "limited memory * default if n1 > 75");
    intOpts_["Hessian frequency"] = "999999 * for full Hessian (never reset)";
    intOpts_["Hessian updates"] = "10 * for limited memory Hessian";
    intOpts_["Hessian flush"] = "999999 * no flushing";

    // * Frequencies
    intOpts_["Check frequency"] = "60 * test row residuals kAx - sk";
    intOpts_["Expand frequency"] = "10000 * for anti-cycling procedure";
    intOpts_["Factorization frequency"] = "50 * 100 for LPs";
    intOpts_["Save frequency"] = "100 * save basis map";

    // * LUSOL options
    realOpts_["LU factor tolerance"] = "3.99 * for NP (100.0 for LP)";
    realOpts_["LU update tolerance"] = "3.99 * for NP ( 10.0 for LP)";
    realOpts_["LU singularity tolerance"] = "3.2e-11";
    strOpts_["LU"] = spair("partial pivoting|rook pivoting|complete pivoting",
                           "LU partial pivoting * default threshold pivoting strategy\n"
                           "LU rook pivoting * threshold rook pivoting\n"
                           "LU complete pivoting * threshold complete pivoting");

    // * Basis files
    intOpts_["Old basis file"] = "0 * input basis map";
    intOpts_["New basis file"] = "0 * output basis map";
    intOpts_["Backup basis file"] = "0 * output extra basis map";
    intOpts_["Insert file"] = "0 * input in industry format";
    intOpts_["Punch file"] = "0 * output Insert data";
    intOpts_["Load file"] = "0 * input names and values";
    intOpts_["Dump file"] = "0 * output Load data";
    intOpts_["Solution file"] = "0 * different from printed solution";

    // * Partitions of cw, iw, rw
    // Total character workspace lencw *
    // Total integer workspace leniw *
    // Total real workspace lenrw *
    // User character workspace 500 *
    // User integer workspace 500 *
    // User real workspace 500 *

    // * Miscellaneous
    intOpts_["Debug level"] = "0 * for developers";
    strOpts_["Sticky parameters"] = spair("No|Yes", "No * Yes makes parameter values persist");
    intOpts_["Timing level"] = "3 * print cpu times";

    // Add the Snopt Options
    for (std::map<std::string, std::string>::const_iterator it = intOpts_.begin();
         it != intOpts_.end(); ++it)
      addOption(it->first, OT_INTEGER, GenericType(), it->second);
    for (std::map<std::string, std::string>::const_iterator it = realOpts_.begin();
         it != realOpts_.end(); ++it)
      addOption(it->first, OT_REAL, GenericType(), it->second);
    for (std::map<std::string, std::pair<std::string, std::string> >::const_iterator
           it = strOpts_.begin(); it != strOpts_.end(); ++it)
      addOption(it->first, OT_STRING, GenericType(), it->second.second, it->second.first);
  }

  SnoptInterface::~SnoptInterface() {
  }

  void SnoptInterface::init() {
    // Read in casadi options
    detect_linear_ = option("detect_linear");

    // Call the init method of the base class
    Nlpsol::init();

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
      // Detect dependencies w.r.t. gradF
      std::vector<const bvec_t*> arg(jacF_.sz_arg(), 0);
      std::vector<bvec_t> input_v_x(nx_, 1);
      arg[GRADF_X] = getPtr(input_v_x);
      std::vector<bvec_t*> res(jacF_.sz_res(), 0);
      std::vector<bvec_t> output_v(jacF_.nnz_out(0), 0);
      res[0] = getPtr(output_v);

      // Perform a single dependency sweep
      jacF_(arg, res);

      // Harvest the results
      for (int j = 0; j < nx_; ++j) {
        if (jacF_.sparsity_out(0).colind(j) == jacF_.sparsity_out(0).colind(j+1)) {
          x_type_f_[j] = 0;
        } else {
          x_type_f_[j] = output_v[jacF_.sparsity_out(0).colind(j)]?  2 : 1;
        }
      }

      if (!jacG_.isNull()) {
        // Detect dependencies w.r.t. jacG
        std::vector<const bvec_t*> arg(jacG_.sz_arg(), 0);
        std::vector<bvec_t> input_v_x(nx_, 1);
        arg[JACG_X] = getPtr(input_v_x);
        std::vector<bvec_t*> res(jacG_.sz_res(), 0);
        std::vector<bvec_t> output_v(jacG_.nnz_out(0), 0);
        res[0] = getPtr(output_v);

        // Perform a single dependency sweep
        jacG_(arg, res);

        Sparsity out_trans = jacG_.sparsity_out(0).T();
        std::vector<bvec_t> output_v_trans(output_v.size());
        std::vector<int> iw(jacG_.size2_out(0));
        casadi_trans(getPtr(output_v), jacG_.sparsity_out(0),
                     getPtr(output_v_trans), out_trans, getPtr(iw));

        for (int j = 0; j < nx_; ++j) {  // Harvest the results
          if (jacG_.sparsity_out(0).colind(j) == jacG_.sparsity_out(0).colind(j+1)) {
            x_type_g_[j] = 0;
          } else {
            bool linear = true;
            for (int k = jacG_.sparsity_out(0).colind(j);
                 k<jacG_.sparsity_out(0).colind(j+1); ++k) {
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

    // An encoding of the desired sorting pattern
    // Digits xy  with x correspodning to x_type_f_ and y corresponding to x_type_g_
    std::vector<int> order_template = {22, 12, 2, 21, 20, 11, 10, 1, 0};

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

    // Here follows the core of the mapping
    //  Two integer matrices are constructed:
    //  one with gradF sparsity, and one with jacG sparsity
    //  the integer values denote the nonzero locations into the original gradF/jacG
    //  but with a special encoding: entries of gradF are encoded "-1-i" and
    //  entries of jacG are encoded "1+i"
    //  "0" is to be interpreted not as an index but as a literal zero

    IM mapping_jacG  = IM(0, nx_);
    IM mapping_gradF = IM(jacF_.sparsity_out(0),
                          range(-1, -1-jacF_.nnz_out(0), -1));

    if (!jacG_.isNull()) {
      mapping_jacG = IM(jacG_.sparsity_out(0), range(1, jacG_.nnz_out(0)+1));
    }

    // First, remap jacG
    A_structure_ = mapping_jacG(g_order_, x_order_);

    m_ = ng_;

    // Construct the linear objective row
    IM d = mapping_gradF(Slice(0), x_order_);

    std::vector<int> ii = mapping_gradF.sparsity().get_col();
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
    d = sparsify(d);

    jacF_row_ = d.nnz() != 0;
    if (jacF_row_) {  // We need an objective gradient row
      A_structure_ = vertcat(A_structure_, d);
      m_ +=1;
    }
    iObj_ = jacF_row_ ? (m_ - 1) : -1;

    // Is the A matrix completely empty?
    dummyrow_ = A_structure_.nnz() == 0;  // Then we need a dummy row
    if (dummyrow_) {
      IM dummyrow = IM(1, nx_);
      dummyrow(0, 0) = 0;
      A_structure_ = vertcat(A_structure_, dummyrow);
      m_+=1;
    }

    // We don't need a dummy row if a linear objective row is present
    casadi_assert(!(dummyrow_ && jacF_row_));

    // Allocate data structures needed in evaluate
    bl_.resize(nx_+m_);
    bu_.resize(nx_+m_);
    hs_.resize(nx_+m_);
    xk_.resize(nx_+m_);
    pi_.resize(m_);
    rc_.resize(nx_+m_);
    A_data_.resize(A_structure_.nnz());

    // Reset the counters
    t_eval_grad_f_ = t_eval_jac_g_ = t_callback_fun_ = t_mainloop_ = 0;
    n_eval_grad_f_ = n_eval_jac_g_ = n_callback_fun_ = n_iter_ = 0;

    // Allocate temporary memory
    alloc_w(nx_, true); // xk2_
    alloc_w(ng_, true); // lam_gk_
    alloc_w(nx_, true); // lam_xk_
    alloc_w(ng_, true); // gk_
    alloc_w(jacF_.nnz_out(0), true); // jac_fk_
    alloc_w(jacG_.nnz_out(0), true); // jac_gk_
  }

  void SnoptInterface::passOptions(snoptProblemC &probC) {
    for (std::map<std::string, std::string>::const_iterator it = intOpts_.begin();
         it != intOpts_.end(); ++it)
      if (hasSetOption(it->first)) {
        int value = option(it->first);
        casadi_assert(it->first.size() <= 55);
        int Error = probC.setIntParameter(it->first.c_str(), value);
        casadi_assert_message(0 == Error, "snopt error setting option \"" + it->first + "\"");
      }
    for (std::map<std::string, std::string>::const_iterator it = realOpts_.begin();
         it != realOpts_.end(); ++it)
      if (hasSetOption(it->first)) {
        double value = option(it->first);
        casadi_assert(it->first.size() <= 55);
        int Error = probC.setRealParameter(it->first.c_str(), value);
        casadi_assert_message(0 == Error, "snopt error setting option \"" + it->first + "\"");
      }
    for (std::map<std::string, std::pair<std::string, std::string> >::const_iterator
           it = strOpts_.begin(); it != strOpts_.end(); ++it)
      if (hasSetOption(it->first)) {
        std::string value = option(it->first);
        std::string buffer = it->first;
        buffer.append(" ");
        buffer.append(value);
        casadi_assert(buffer.size() <= 72);
        int Error = probC.setParameter(buffer.c_str());
        casadi_assert_message(0 == Error, "snopt error setting option \"" + it->first + "\"");
      }
  } // passOptions()

  std::string SnoptInterface::formatStatus(int status) const {
    if (status_.find(status) == status_.end()) {
      std::stringstream ss;
      ss << "Unknown status: " << status;
      return ss.str();
    } else {
      return (*status_.find(status)).second;
    }
  }

  void SnoptInterface::reset(void* mem, const double**& arg, double**& res, int*& iw, double*& w) {
    // Reset the base classes
    Nlpsol::reset(mem, arg, res, iw, w);

    // Work vectors
    xk2_ = w; w += nx_;
    lam_gk_ = w; w += ng_;
    lam_xk_ = w; w += nx_;
    gk_ = w; w += ng_;
    jac_fk_ = w; w += jacF_.nnz_out(0);
    jac_gk_ = w; w += jacG_.nnz_out(0);
  }

  void SnoptInterface::solve(void* mem) {
    for (int i=0; i<NLPSOL_NUM_IN; ++i) {
      const double *v;
      switch (i) {
      case NLPSOL_X0: v = x0_; break;
      case NLPSOL_P: v = p_; break;
      case NLPSOL_LBX: v = lbx_; break;
      case NLPSOL_UBX: v = ubx_; break;
      case NLPSOL_LBG: v = lbg_; break;
      case NLPSOL_UBG: v = ubg_; break;
      case NLPSOL_LAM_X0: v = lam_x0_; break;
      case NLPSOL_LAM_G0: v = lam_g0_; break;
      default: casadi_assert(0);
      }
      casadi_copy(v, nnz_in(i), input(i).ptr());
    }

    // Check the provided inputs
    checkInputs(mem);

    // Evaluate gradF and jacG at initial value
    if (!jacG_.isNull()) {
      std::fill_n(arg_, jacG_.n_in(), nullptr);
      arg_[JACG_X] = input(NLPSOL_X0).ptr();
      arg_[JACG_P] = input(NLPSOL_P).ptr();
      std::fill_n(res_, jacG_.n_out(), nullptr);
      res_[0] = jac_gk_;
      jacG_(arg_, res_, iw_, w_, 0);
    }
    std::fill_n(arg_, jacF_.n_in(), nullptr);
    arg_[GRADF_X] = input(NLPSOL_X0).ptr();
    arg_[GRADF_P] = input(NLPSOL_P).ptr();
    std::fill_n(res_, jacF_.n_out(), nullptr);
    res_[0] = jac_fk_;
    jacF_(arg_, res_, iw_, w_, 0);

    // perform the mapping:
    // populate A_data_ (the nonzeros of A)
    // with numbers pulled from jacG and gradF
    for (int k = 0; k < A_structure_.nnz(); ++k) {
      int i = A_structure_.data()[k];
      if (i == 0) {
        A_data_[k] = 0;
      } else if (i > 0) {
        A_data_[k] = jac_gk_[i-1];
      } else {
        A_data_[k] = jac_fk_[-i-1];
      }
    }

    // Obtain constraint offsets for linear constraints
    if (detect_linear_) {
      std::fill_n(arg_, nlp_.n_in(), nullptr);
      arg_[NL_P] = input(NLPSOL_P).ptr();
      std::fill_n(res_, nlp_.n_out(), nullptr);
      res_[NL_F] = &fk_;
      res_[NL_G] = gk_;
      nlp_(arg_, res_, iw_, w_, 0);
    }

    // Obtain sparsity pattern of A (Fortran is Index-1 based, but the C++ wrappers are Index-0)
    std::vector<int> row(A_structure_.nnz());
    for (int k = 0; k < A_structure_.nnz(); ++k) {
      row[k] = A_structure_.row()[k];
    }

    std::vector<int> col(nx_+1);
    for (int k = 0; k < nx_+1; ++k) {
      col[k] = A_structure_.colind()[k];
    }

    // Obtain initial guess and bounds through the mapping
    for (int k = 0; k < nx_; ++k) {
      int kk = x_order_[k];
      bl_[k] = input(NLPSOL_LBX).data()[kk];
      bu_[k] = input(NLPSOL_UBX).data()[kk];
      xk_[k] = input(NLPSOL_X0).data()[kk];
    }
    for (int k = 0; k < ng_; ++k) {
      int kk = g_order_[k];
      if (g_type_[kk] < 2) {
        // casadi_error("woops");
        bl_[nx_+k] = input(NLPSOL_LBG).data()[kk]-gk_[kk];
        bu_[nx_+k] = input(NLPSOL_UBG).data()[kk]-gk_[kk];
      } else {
        bl_[nx_+k] = input(NLPSOL_LBG).data()[kk];
        bu_[nx_+k] = input(NLPSOL_UBG).data()[kk];
      }
      xk_[nx_+k] = input(NLPSOL_LAM_G0).data()[kk];
    }

    // Objective row / dummy row should be unbounded
    if (dummyrow_ || jacF_row_) {
      bl_.back() = -1e20;  // -std::numeric_limits<double>::infinity();
      bu_.back() = 1e20;  // std::numeric_limits<double>::infinity();
    }

    int n = nx_;
    int nea = A_structure_.nnz();
    double ObjAdd = 0;

    casadi_assert(m_ > 0);
    casadi_assert(n > 0);
    casadi_assert(nea > 0);
    casadi_assert(row.size() == nea);
    casadi_assert(hs_.size() == n+m_);
    casadi_assert(col.size() == n+1);
    casadi_assert(A_structure_.nnz() == nea);
    casadi_assert(bl_.size() == n+m_);
    casadi_assert(bu_.size() == n+m_);
    casadi_assert(pi_.size() == m_);
    casadi_assert(xk_.size() == n+m_);
    casadi_assert(col.at(0) == 0);
    casadi_assert(col.at(n) == nea);

    // Pointer magic, courtesy of Greg
    casadi_assert_message(!jacF_.isNull(), "blaasssshc");

    // Outputs
    double Obj = 0; // TODO(Greg): get this from snopt

    if (option("summary"))
      summaryOn();
    else
      summaryOff();
    snoptProblemC snoptProbC = snoptProblemC();
    if (hasSetOption("specs file")) {
      std::string specs_file = option("specs file");
      snoptProbC.setSpecsFile(specs_file.c_str());
    }
    if (hasSetOption("print file")) {
      std::string print_file = option("print file");
      snoptProbC.setPrintFile(print_file.c_str());
    }

    snoptProbC.setProblemSize(m_, nx_, nnCon_, nnJac_, nnObj_);
    snoptProbC.setObjective(iObj_, ObjAdd);
    snoptProbC.setJ(nea, getPtr(A_data_), getPtr(row), getPtr(col));
    snoptProbC.setX(getPtr(bl_), getPtr(bu_), getPtr(xk_), getPtr(pi_), getPtr(rc_), getPtr(hs_));
    snoptProbC.setUserFun(userfunPtr);
    snoptProbC.setSTOP(snStopPtr);
    passOptions(snoptProbC);

    // user data
    int iulen = 8;
    std::vector<int> iu(iulen);
    SnoptInterface* source = this;
    memcpy(&(iu[0]), &source, sizeof(SnoptInterface*));
    snoptProbC.setUserI(getPtr(iu), iulen);

    // Run SNOPT
    double time0 = clock();
    int info = snoptProbC.solve(optionEnumValue("start"));
    casadi_assert_message(99 != info, "snopt problem set up improperly");
    t_mainloop_ = static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;

    stats_["return_status"] = info;

    // Store results into output
    for (int k = 0; k < nx_; ++k) {
      int kk =  x_order_[k];
      output(NLPSOL_X).data()[kk] = xk_[k];
      output(NLPSOL_LAM_X).data()[kk] = -rc_[k];
    }

    setOutput(Obj+ (jacF_row_? xk_[nx_+ng_] : 0), NLPSOL_F);

    for (int k = 0; k < ng_; ++k) {
      int kk = g_order_[k];
      output(NLPSOL_LAM_G).data()[kk] = -rc_[nx_+k];
      output(NLPSOL_G).data()[kk] = xk_[nx_+k];  // TODO(Joris): this is not quite right
      // mul_no_alloc
    }

    // todo(Greg): get these from snopt
    // we overwrite F and G for now because the current snopt interface
    // doesn't provide F, and the above comment suggests that G is wrong
    std::fill_n(arg_, nlp_.n_in(), nullptr);
    arg_[NL_X] = output(NLPSOL_X).ptr();
    arg_[NL_P] = input(NLPSOL_P).ptr();
    std::fill_n(res_, nlp_.n_out(), nullptr);
    res_[NL_F] = &fk_;
    res_[NL_G] = gk_;
    nlp_(arg_, res_, iw_, w_, 0);

    setOutputNZ(&fk_, NLPSOL_F);
    setOutputNZ(gk_, NLPSOL_G);

    // print timing information
    // save state
    std::ios state(NULL);
    state.copyfmt(userOut());
    const int w_time = 7;
    const int p_time = 3;
    const int w_ms = 7;
    const int p_ms = 2;
    const int w_n = 5;

    // restore state
    userOut().copyfmt(state);

    // Reset the counters
    t_eval_grad_f_ = t_eval_jac_g_ = t_callback_fun_ = t_mainloop_ = 0;
    n_eval_grad_f_ = n_eval_jac_g_ = n_callback_fun_ = n_iter_ = 0;

    for (int i=0; i<NLPSOL_NUM_OUT; ++i) {
      double **v;
      switch (i) {
      case NLPSOL_X: v = &x_; break;
      case NLPSOL_F: v = &f_; break;
      case NLPSOL_G: v = &g_; break;
      case NLPSOL_LAM_X: v = &lam_x_; break;
      case NLPSOL_LAM_G: v = &lam_g_; break;
      case NLPSOL_LAM_P: v = &lam_p_; break;
      default: casadi_assert(0);
      }
      casadi_copy(output(i).ptr(), nnz_out(i), *v);
    }
  }

  void SnoptInterface::setOptionsFromFile(const std::string & file) {
  }

  void SnoptInterface::
  userfun(int* mode, int nnObj, int nnCon, int nnJac, int nnL, int neJac,
          double* x, double* fObj, double*gObj, double* fCon, double* gCon,
          int nState, char* cu, int lencu, int* iu, int leniu, double* ru,
          int lenru) {
    try {
      double time0 = clock();

      casadi_assert_message(nnCon_ == nnCon, "Con " << nnCon_ << " <-> " << nnCon);
      casadi_assert_message(nnObj_ == nnObj, "Obj " << nnObj_ << " <-> " << nnObj);
      casadi_assert_message(nnJac_ == nnJac, "Jac " << nnJac_ << " <-> " << nnJac);

      // Get reduced decision variables
      casadi_fill(xk2_, nx_, 0.);
      for (int k = 0; k < nnObj; ++k) {
        if (x_type_f_[x_order_[k]] == 2) {
          xk2_[x_order_[k]] = x[k];
        }
      }

      // Evaluate gradF with the linear variables put to zero
      std::fill_n(arg_, jacF_.n_in(), nullptr);
      arg_[NL_X] = xk2_;
      arg_[NL_P] = input(NLPSOL_P).ptr();
      std::fill_n(res_, jacF_.n_out(), nullptr);
      res_[0] = jac_fk_;
      res_[1] = fObj;
      jacF_(arg_, res_, iw_, w_, 0);

      // provide nonlinear part of objective gradient to SNOPT
      for (int k = 0; k < nnObj; ++k) {
        int i = x_order_[k];
        if (x_type_f_[i] == 2) {
          int el = jacF_.sparsity_out(0).colind(i);
          if (jacF_.sparsity_out(0).colind(i+1) > el) {
            gObj[k] = jac_fk_[el];
          } else {
            gObj[k] = 0;
          }
        } else {
          gObj[k] = 0;
        }
      }

      jacF_.sparsity_out(0).sanity_check(true);
      jacF_.sparsity_out(0).sanity_check(false);

      // timing and counters
      t_eval_grad_f_ += static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;
      n_eval_grad_f_ += 1;


      time0 = clock();
      if (!jacG_.isNull()) {
        // Get reduced decision variables
        casadi_fill(xk2_, nx_, 0.);
        for (int k = 0; k < nnJac; ++k) {
          xk2_[x_order_[k]] = x[k];
        }

        // Evaluate jacG with the linear variabes put to zero
        std::fill_n(arg_, jacG_.n_in(), nullptr);
        arg_[JACG_X] = xk2_;
        arg_[JACG_P] = input(NLPSOL_P).ptr();
        std::fill_n(res_, jacG_.n_out(), nullptr);
        res_[0] = jac_gk_;
        res_[GRADF_G] = gk_;
        jacG_(arg_, res_, iw_, w_, 0);

        // provide nonlinear part of constraint jacobian to SNOPT
        int kk = 0;
        for (int j = 0; j < nnJac; ++j) {
          for (int k = A_structure_.colind(j); k < A_structure_.sparsity().colind(j+1); ++k) {
            if (A_structure_.row(k) >= nnCon) break;
            int i = A_structure_.data()[k];
            if (i > 0) {
              gCon[kk++] = jac_gk_[i-1];
            }
          }
        }

        casadi_assert(kk == 0 || kk == neJac);

        // provide nonlinear part of objective to SNOPT
        for (int k = 0; k < nnCon; ++k) {
          fCon[k] = gk_[g_order_[k]];
        }

        // timing and counters
        t_eval_jac_g_ += static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;
        n_eval_jac_g_ += 1;
      }

    } catch(std::exception& ex) {
      userOut<true, PL_WARN>() << "eval_nlp failed: " << ex.what() << std::endl;
      *mode = -1;  // Reduce step size - we've got problems
      return;
    }
  }

  void SnoptInterface::
  userfunPtr(int * mode, int* nnObj, int * nnCon, int *nJac,
             int *nnL, int * neJac, double *x, double *fObj,
             double *gObj, double * fCon, double* gCon,
             int* nState, char* cu, int* lencu, int* iu,
             int* leniu, double* ru, int *lenru) {
    SnoptInterface* interface;
    memcpy(&interface, &(iu[0]), sizeof(SnoptInterface*));
    interface->userfun(mode, *nnObj, *nnCon, *nJac, *nnL, *neJac,
                       x, fObj, gObj, fCon, gCon, *nState,
                       cu, *lencu, iu, *leniu, ru, *lenru);
  }

  void SnoptInterface::
  snStopPtr(int *iAbort, int KTcond[], int *MjrPrt, int *minimz,
            int *m, int *maxS, int *n, int *nb,
            int *nnCon0, int *nnCon, int *nnObj0, int *nnObj, int *nS,
            int *itn, int *nMajor, int *nMinor, int *nSwap,
            double *condHz, int *iObj, double *sclObj, double *ObjAdd,
            double *fMrt, double *PenNrm, double *step,
            double *prInf, double *duInf, double *vimax, double *virel, int hs[],
            int *ne, int *nlocJ, int locJ[], int indJ[], double Jcol[], int *negCon,
            double Ascale[], double bl[], double bu[], double fCon[], double gCon[], double gObj[],
            double yCon[], double pi[], double rc[], double rg[], double x[],
            char cu[], int *lencu, int iu[], int *leniu, double ru[], int *lenru,
            char cw[], int *lencw, int iw[], int *leniw, double rw[], int *lenrw ) {
    SnoptInterface* self;
    memcpy(&self, &(iu[0]), sizeof(SnoptInterface*));
    try {
      self->n_iter_+=1;
      if (!self->fcallback_.isNull()) {
        double time0 = clock();
        for (int k = 0; k < self->nx_; ++k) {
          int kk = self->x_order_[k];
          self->output(NLPSOL_X).data()[kk] = self->xk_[k];
        }
        for (int k = 0; k < self->ng_; ++k) {
          int kk =  self->g_order_[k];
          self->output(NLPSOL_G).data()[kk] = self->xk_[self->nx_+k];
        }

        for (int i=0; i<NLPSOL_NUM_OUT; ++i) {
          self->fcallback_.setInput(self->output(i), i);
        }
        self->fcallback_.evaluate();
        double ret_double;
        self->fcallback_.getOutput(ret_double);
        *iAbort = static_cast<int>(ret_double);
        self->t_callback_fun_ += static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;
        self->n_callback_fun_ += 1;
      }
    } catch(std::exception& ex) {
      if (self->option("iteration_callback_ignore_errors")) {
        userOut<true, PL_WARN>() << "callback: " << ex.what() << std::endl;
      } else {
        throw ex;
      }
    }
  }
}  // namespace casadi
