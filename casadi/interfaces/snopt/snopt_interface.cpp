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

    // Printing
    intOpts_["Major print level"] = "1 * 1-line major iteration log";
    intOpts_["Minor print level"] = "1 * 1-line minor iteration log";
    // special om["Print file"] = OT_S; //  * specified by subroutine sn_init
    // special om["Summary file"] = OT_S; //  * specified by subroutine sn_init
    intOpts_["Print frequency"] = "100 * minor iterations log on Print file";
    intOpts_["Summary frequency"] = "100 * minor iterations log on Summary file";
    strOpts_["Solution"] = "Yes * on the Print file";

    // * Suppress options listing * options are normally listed
    strOpts_["System information"] = "No * Yes prints more system information";

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
    strOpts_["QPSolver"] = "Cholesky * default";
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
    strOpts_["Hessian"] = "full memory * default if n1 ≤ 75\n"
      "limited memory * default if n1 > 75";
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
    strOpts_["LU"] = "LU partial pivoting * default threshold pivoting strategy\n"
      "LU rook pivoting * threshold rook pivoting\n"
      "LU complete pivoting * threshold complete pivoting";

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
    strOpts_["Sticky parameters"] = "No * Yes makes parameter values persist";
    intOpts_["Timing level"] = "3 * print cpu times";

    // Add the Snopt Options
    for (auto it = intOpts_.begin(); it != intOpts_.end(); ++it)
      addOption(it->first, OT_INTEGER, GenericType(), it->second);
    for (auto it = realOpts_.begin(); it != realOpts_.end(); ++it)
      addOption(it->first, OT_REAL, GenericType(), it->second);
    for (auto it = strOpts_.begin(); it != strOpts_.end(); ++it)
      addOption(it->first, OT_STRING, GenericType(), it->second);
  }

  SnoptInterface::~SnoptInterface() {
  }

  void SnoptInterface::init() {
    // Read in casadi options
    detect_linear_ = option("detect_linear");

    // Call the init method of the base class
    Nlpsol::init();

    // Get/generate required functions
    setup_jac_f();
    setup_jac_g();

    // prepare the mapping for constraints
    nnJac_ = nx_;
    nnObj_ = nx_;
    nnCon_ = ng_;

    // Here follows the core of the mapping
    //  Two integer matrices are constructed:
    //  one with gradF sparsity, and one with jacG sparsity
    //  the integer values denote the nonzero locations into the original gradF/jacG
    //  but with a special encoding: entries of gradF are encoded "-1-i" and
    //  entries of jacG are encoded "1+i"
    //  "0" is to be interpreted not as an index but as a literal zero

    IM mapping_jacG  = IM(0, nx_);
    IM mapping_gradF = IM(jac_f_fcn_.sparsity_out(1),
                          range(-1, -1-jac_f_fcn_.nnz_out(1), -1));

    if (!jac_g_fcn_.isNull()) {
      mapping_jacG = IM(jacg_sp_, range(1, jacg_sp_.nnz()+1));
    }

    // First, remap jacG
    A_structure_ = mapping_jacG;

    m_ = ng_;

    // Construct the linear objective row
    IM d = mapping_gradF(Slice(0), Slice());

    std::vector<int> ii = mapping_gradF.sparsity().get_col();
    for (int j = 0; j < nnObj_; ++j) {
      if (d.colind(j) != d.colind(j+1)) {
        int k = d.colind(j);
        d[k] = 0;
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

    // Allocate temporary memory
    alloc_w(nx_, true); // xk2_
    alloc_w(ng_, true); // lam_gk_
    alloc_w(nx_, true); // lam_xk_
    alloc_w(ng_, true); // gk_
    alloc_w(jac_f_fcn_.nnz_out(1), true); // jac_fk_
    if (!jacg_sp_.isNull()) {
      alloc_w(jacg_sp_.nnz(), true); // jac_gk_
    }
  }

  void SnoptInterface::init_memory(Memory& mem) const {
    Nlpsol::init_memory(mem);
    SnoptMemory& m = dynamic_cast<SnoptMemory&>(mem);

    // Allocate data structures needed in evaluate
    m.A_data.resize(A_structure_.nnz());

    // Reset the counters
    m.t_eval_grad_f = m.t_eval_jac_g = m.t_callback_fun = m.t_mainloop = 0;
    m.n_eval_grad_f = m.n_eval_jac_g = m.n_callback_fun = m.n_iter = 0;
  }

  std::string SnoptInterface::formatStatus(int status) const {
    if (status_.find(status) == status_.end()) {
      std::stringstream ss;
      ss << "Unknown status: " << status;
      return ss.str();
    } else {
      return (*status_.find(status)).second;
    }
  }

  void SnoptInterface::set_work(Memory& mem, const double**& arg, double**& res,
                                int*& iw, double*& w) const {
    SnoptMemory& m = dynamic_cast<SnoptMemory&>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Work vectors
    m.xk2 = w; w += nx_;
    m.lam_gk = w; w += ng_;
    m.lam_xk = w; w += nx_;
    m.gk = w; w += ng_;
    m.jac_fk = w; w += jac_f_fcn_.nnz_out(1);
    if (!jacg_sp_.isNull()) {
      m.jac_gk = w; w += jacg_sp_.nnz();
    }
  }

  void SnoptInterface::solve(Memory& mem) const {
    SnoptMemory& m = dynamic_cast<SnoptMemory&>(mem);

    // Check the provided inputs
    checkInputs(mem);

    // Memory object
    snProblem prob;

    // Evaluate gradF and jacG at initial value
    if (!jac_g_fcn_.isNull()) {
      std::fill_n(m.arg, jac_g_fcn_.n_in(), nullptr);
      m.arg[0] = m.x0;
      m.arg[1] = m.p;
      std::fill_n(m.res, jac_g_fcn_.n_out(), nullptr);
      m.res[0] = 0;
      m.res[1] = m.jac_gk;
      jac_g_fcn_(m.arg, m.res, m.iw, m.w, 0);
    }
    calc_jac_f(m, m.x0, m.p, 0, m.jac_fk);

    // perform the mapping:
    // populate A_data_ (the nonzeros of A)
    // with numbers pulled from jacG and gradF
    for (int k = 0; k < A_structure_.nnz(); ++k) {
      int i = A_structure_.data()[k];
      if (i == 0) {
        m.A_data[k] = 0;
      } else if (i > 0) {
        m.A_data[k] = m.jac_gk[i-1];
      } else {
        m.A_data[k] = m.jac_fk[-i-1];
      }
    }

    int n = nx_;
    int nea = A_structure_.nnz();
    double ObjAdd = 0;

    casadi_assert(m_ > 0);
    casadi_assert(n > 0);
    casadi_assert(nea > 0);
    casadi_assert(A_structure_.nnz() == nea);

    // Pointer magic, courtesy of Greg
    casadi_assert_message(!jac_f_fcn_.isNull(), "blaasssshc");

    // Outputs
    double Obj = 0; // TODO(Greg): get this from snopt

    // if (option("summary")) {
    //   summaryOn();
    // } else {
    //   summaryOff();
    // }

    // snInit must be called first.
    //   9, 6 are print and summary unit numbers (for Fortran).
    //   6 == standard out
    int iprint = 9;
    int isumm = 6;
    std::string outname = name_ + ".out";
    snInit(&prob, const_cast<char*>(name_.c_str()),
           const_cast<char*>(outname.c_str()), iprint, isumm);

    // snoptProblemC snoptProbC = snoptProblemC();
    // if (hasSetOption("specs file")) {
    //   std::string specs_file = option("specs file");
    //   snoptProbC.setSpecsFile(specs_file.c_str());
    // }
    // if (hasSetOption("print file")) {
    //   std::string print_file = option("print file");
    //   snoptProbC.setPrintFile(print_file.c_str());
    // }

    // Set the problem size and other data.
    // This will allocate arrays inside snProblem struct.
    setProblemSize(&prob, m_, nx_, nea, nnCon_, nnJac_, nnObj_);
    setObjective(&prob, iObj_, ObjAdd);
    setUserfun(&prob, userfunPtr);

    // user data
    prob.leniu = 1;
    prob.iu = &m.memind;

    // Pass bounds
    casadi_copy(m.lbx, nx_, prob.bl);
    casadi_copy(m.ubx, nx_, prob.bu);
    casadi_copy(m.lbg, ng_, prob.bl + nx_);
    casadi_copy(m.ubg, ng_, prob.bu + nx_);

    // Initialize states and slack
    casadi_fill(prob.hs, ng_ + nx_, 0);
    casadi_copy(m.x0, nx_, prob.x);
    casadi_fill(prob.x + nx_, ng_, 0.);

    // Initialize multipliers
    casadi_copy(m.lam_g0, ng_, prob.pi);

    // Set up Jacobian matrix
    casadi_copy(A_structure_.colind(), A_structure_.size2()+1, prob.locJ);
    casadi_copy(A_structure_.row(), A_structure_.nnz(), prob.indJ);
    casadi_copy(get_ptr(m.A_data), A_structure_.nnz(), prob.valJ);

    for (auto&& pp : intOpts_) {
      if (hasSetOption(pp.first)) {
        int value = option(pp.first);
        casadi_assert(pp.first.size() <= 55);
        int Error = setIntParameter(&prob, const_cast<char*>(pp.first.c_str()), value);
        casadi_assert_message(0 == Error, "snopt error setting option \"" + pp.first + "\"");
      }
    }
    for (auto&& pp : realOpts_) {
      if (hasSetOption(pp.first)) {
        double value = option(pp.first);
        casadi_assert(pp.first.size() <= 55);
        int Error = setRealParameter(&prob, const_cast<char*>(pp.first.c_str()), value);
        casadi_assert_message(0 == Error, "snopt error setting option \"" + pp.first + "\"");
      }
    }
    for (auto&& pp : strOpts_) {
      if (hasSetOption(pp.first)) {
        std::string value = option(pp.first);
        std::string buffer = pp.first;
        buffer.append(" ");
        buffer.append(value);
        casadi_assert(buffer.size() <= 72);
        int Error = setParameter(&prob, const_cast<char*>(buffer.c_str()));
        casadi_assert_message(0 == Error, "snopt error setting option \"" + pp.first + "\"");
      }
    }

    // Run SNOPT
    int Cold = 0;
    int info = solveC(&prob, Cold, &m.fk);
    casadi_assert_message(99 != info, "snopt problem set up improperly");

    // Negate rc to match CasADi's definition
    casadi_scal(nx_ + ng_, -1., prob.rc);

    // Get primal solution
    casadi_copy(prob.x, nx_, m.x);

    // Get dual solution
    casadi_copy(prob.rc, nx_, m.lam_x);
    casadi_copy(prob.rc+nx_, ng_, m.lam_g);

    // Copy optimal cost to output
    if (m.f) *m.f = m.fk;

    // Copy optimal constraint values to output
    casadi_copy(m.gk, ng_, m.g);

    // Free memory
    deleteSNOPT(&prob);
  }

  void SnoptInterface::
  userfun(SnoptMemory &m, int* mode, int nnObj, int nnCon, int nnJac, int nnL, int neJac,
          double* x, double* fObj, double*gObj, double* fCon, double* gCon,
          int nState, char* cu, int lencu, int* iu, int leniu, double* ru,
          int lenru) const {
    try {
      double time0 = clock();

      casadi_assert_message(nnCon_ == nnCon, "Con " << nnCon_ << " <-> " << nnCon);
      casadi_assert_message(nnObj_ == nnObj, "Obj " << nnObj_ << " <-> " << nnObj);
      casadi_assert_message(nnJac_ == nnJac, "Jac " << nnJac_ << " <-> " << nnJac);

      // Get reduced decision variables
      casadi_fill(m.xk2, nx_, 0.);
      for (int k = 0; k < nnObj; ++k) m.xk2[k] = x[k];

      // Evaluate gradF with the linear variables put to zero
      calc_jac_f(m, m.xk2, m.p, fObj, m.jac_fk);

      // provide nonlinear part of objective gradient to SNOPT
      for (int k = 0; k < nnObj; ++k) {
        int el = jac_f_fcn_.sparsity_out(1).colind(k);
        if (jac_f_fcn_.sparsity_out(1).colind(k+1) > el) {
          gObj[k] = m.jac_fk[el];
        } else {
          gObj[k] = 0;
        }
      }

      jac_f_fcn_.sparsity_out(1).sanity_check(true);
      jac_f_fcn_.sparsity_out(1).sanity_check(false);

      // timing and counters
      m.t_eval_grad_f += static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;
      m.n_eval_grad_f += 1;

      time0 = clock();
      if (!jac_g_fcn_.isNull()) {
        // Get reduced decision variables
        casadi_fill(m.xk2, nx_, 0.);
        for (int k = 0; k < nnJac; ++k) {
          m.xk2[k] = x[k];
        }

        // Evaluate jacG with the linear variabes put to zero
        std::fill_n(m.arg, jac_g_fcn_.n_in(), nullptr);
        m.arg[0] = m.xk2;
        m.arg[1] = m.p;
        std::fill_n(m.res, jac_g_fcn_.n_out(), nullptr);
        m.res[0] = m.gk;
        m.res[1] = m.jac_gk;
        jac_g_fcn_(m.arg, m.res, m.iw, m.w, 0);

        // provide nonlinear part of constraint jacobian to SNOPT
        int kk = 0;
        for (int j = 0; j < nnJac; ++j) {
          for (int k = A_structure_.colind(j); k < A_structure_.sparsity().colind(j+1); ++k) {
            if (A_structure_.row(k) >= nnCon) break;
            int i = A_structure_.data()[k];
            if (i > 0) {
              gCon[kk++] = m.jac_gk[i-1];
            }
          }
        }

        casadi_assert(kk == 0 || kk == neJac);

        // provide nonlinear part of objective to SNOPT
        for (int k = 0; k < nnCon; ++k) {
          fCon[k] = m.gk[k];
        }

        // timing and counters
        m.t_eval_jac_g += static_cast<double>(clock()-time0)/CLOCKS_PER_SEC;
        m.n_eval_jac_g += 1;
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
    SnoptMemory &m = *SnoptMemory::mempool.at(iu[0]);
    m.self.userfun(m, mode, *nnObj, *nnCon, *nJac, *nnL, *neJac,
                   x, fObj, gObj, fCon, gCon, *nState,
                   cu, *lencu, iu, *leniu, ru, *lenru);
  }

  SnoptMemory::SnoptMemory(const SnoptInterface& self) : self(self) {
    // Put in memory pool
    auto mem_it = std::find(mempool.begin(), mempool.end(), nullptr);
    if (mem_it==mempool.end()) {
      // Append to end
      memind = mempool.size();
      mempool.push_back(this);
    } else {
      // Reuse freed element
      memind = mem_it - mempool.begin();
      *mem_it = this;
    }
  }

  SnoptMemory::~SnoptMemory() {
    // Remove from memory pool
    auto mem_it = std::find(mempool.begin(), mempool.end(), this);
    if (mem_it==mempool.end()) {
      casadi_warning("SNOPT memory pool failure");
    } else {
      *mem_it = nullptr;
    }
  }

  std::vector<SnoptMemory*> SnoptMemory::mempool;

}  // namespace casadi
