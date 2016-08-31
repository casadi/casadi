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
    plugin->version = 30;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_SNOPT_EXPORT casadi_load_nlpsol_snopt() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_snopt);
  }

  SnoptInterface::SnoptInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }

  SnoptInterface::~SnoptInterface() {
    clear_memory();
  }

  Options SnoptInterface::options_
  = {{&Nlpsol::options_},
     {{"snopt",
       {OT_DICT,
        "Options to be passed to SNOPT"}},
      {"start",
       {OT_STRING,
        "Warm-start options for Worhp: cold|warm|hot"}}
     }
  };

  void SnoptInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Default: cold start
    Cold_ = 0;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="snopt") {
        opts_ = op.second;
      } else if (op.first=="start") {
        std::string start = op.second.to_string();
        if (start=="cold") {
          Cold_ = 0;
        } else if (start=="warm") {
          Cold_ = 1;
        } else if (start=="hot") {
          Cold_ = 2;
        } else {
          casadi_error("Unknown start option: " + start);
        }
      }
    }

    // Get/generate required functions
    jac_f_fcn_ = create_function("nlp_jac_f", {"x", "p"}, {"f", "jac:f:x"});
    jac_g_fcn_ = create_function("nlp_jac_g", {"x", "p"}, {"g", "jac:g:x"});
    jacg_sp_ = jac_g_fcn_.sparsity_out(1);

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

    if (!jac_g_fcn_.is_null()) {
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
        d.nz(k) = 0;
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
    if (!jacg_sp_.is_null()) {
      alloc_w(jacg_sp_.nnz(), true); // jac_gk_
    }
  }

  void SnoptInterface::init_memory(void* mem) const {
    Nlpsol::init_memory(mem);
    auto m = static_cast<SnoptMemory*>(mem);

    // Allocate data structures needed in evaluate
    m->A_data.resize(A_structure_.nnz());

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

  void SnoptInterface::set_work(void* mem, const double**& arg, double**& res,
                                int*& iw, double*& w) const {
    auto m = static_cast<SnoptMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Work vectors
    m->xk2 = w; w += nx_;
    m->lam_gk = w; w += ng_;
    m->lam_xk = w; w += nx_;
    m->gk = w; w += ng_;
    m->jac_fk = w; w += jac_f_fcn_.nnz_out(1);
    if (!jacg_sp_.is_null()) {
      m->jac_gk = w; w += jacg_sp_.nnz();
    }
  }

  void SnoptInterface::solve(void* mem) const {
    auto m = static_cast<SnoptMemory*>(mem);

    // Check the provided inputs
    checkInputs(mem);

    m->fstats.at("mainloop").tic();

    // Memory object
    snProblem prob;

    // Evaluate gradF and jacG at initial value
    const double** arg = m->arg;
    *arg++ = m->x0;
    *arg++ = m->p;
    double** res = m->res;
    *res++ = 0;
    *res++ = m->jac_gk;
    calc_function(m, "nlp_jac_g");
    res = m->res;
    *res++ = 0;
    *res++ = m->jac_fk;
    calc_function(m, "nlp_jac_f");

    // perform the mapping:
    // populate A_data_ (the nonzeros of A)
    // with numbers pulled from jacG and gradF
    for (int k = 0; k < A_structure_.nnz(); ++k) {
      int i = A_structure_.nonzeros()[k];
      if (i == 0) {
        m->A_data[k] = 0;
      } else if (i > 0) {
        m->A_data[k] = m->jac_gk[i-1];
      } else {
        m->A_data[k] = m->jac_fk[-i-1];
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
    casadi_assert_message(!jac_f_fcn_.is_null(), "blaasssshc");

    // Outputs
    //double Obj = 0; // TODO(Greg): get this from snopt

    // snInit must be called first.
    //   9, 6 are print and summary unit numbers (for Fortran).
    //   6 == standard out
    int iprint = 9;
    int isumm = 6;
    std::string outname = name_ + ".out";
    snInit(&prob, const_cast<char*>(name_.c_str()),
           const_cast<char*>(outname.c_str()), iprint, isumm);

    // Set the problem size and other data.
    // This will allocate arrays inside snProblem struct.
    setProblemSize(&prob, m_, nx_, nea, nnCon_, nnJac_, nnObj_);
    setObjective(&prob, iObj_, ObjAdd);
    setUserfun(&prob, userfunPtr);

    // user data
    prob.leniu = 1;
    prob.iu = &m->memind;

    // Pass bounds
    casadi_copy(m->lbx, nx_, prob.bl);
    casadi_copy(m->ubx, nx_, prob.bu);
    casadi_copy(m->lbg, ng_, prob.bl + nx_);
    casadi_copy(m->ubg, ng_, prob.bu + nx_);

    // Initialize states and slack
    casadi_fill(prob.hs, ng_ + nx_, 0);
    casadi_copy(m->x0, nx_, prob.x);
    casadi_fill(prob.x + nx_, ng_, 0.);

    // Initialize multipliers
    casadi_copy(m->lam_g0, ng_, prob.pi);

    // Set up Jacobian matrix
    casadi_copy(A_structure_.colind(), A_structure_.size2()+1, prob.locJ);
    casadi_copy(A_structure_.row(), A_structure_.nnz(), prob.indJ);
    casadi_copy(get_ptr(m->A_data), A_structure_.nnz(), prob.valJ);

    for (auto&& op : opts_) {
      // Replace underscores with spaces
      std::string opname = op.first;
      std::replace(opname.begin(), opname.end(), '_', ' ');

      // Try integer
      if (op.second.can_cast_to(OT_INT)) {
        casadi_assert(opname.size() <= 55);
        int flag = setIntParameter(&prob, const_cast<char*>(opname.c_str()),
                                   op.second.to_int());
        if (flag==0) continue;
      }

      // Try double
      if (op.second.can_cast_to(OT_DOUBLE)) {
        casadi_assert(opname.size() <= 55);
        int flag = setRealParameter(&prob, const_cast<char*>(opname.c_str()),
                                    op.second.to_double());
        if (flag==0) continue;
      }

      // try string
      if (op.second.can_cast_to(OT_STRING)) {
        std::string buffer = opname + " " + op.second.to_string();
        casadi_assert(buffer.size() <= 72);
        int flag = setParameter(&prob, const_cast<char*>(buffer.c_str()));
        if (flag==0) continue;
      }

      // Error if reached this point
      casadi_error("SNOPT error setting option \"" + opname + "\"");
    }

    m->fstats.at("mainloop").toc();

    // Run SNOPT
    int info = solveC(&prob, Cold_, &m->fk);
    casadi_assert_message(99 != info, "snopt problem set up improperly");

    // Negate rc to match CasADi's definition
    casadi_scal(nx_ + ng_, -1., prob.rc);

    // Get primal solution
    casadi_copy(prob.x, nx_, m->x);

    // Get dual solution
    casadi_copy(prob.rc, nx_, m->lam_x);
    casadi_copy(prob.rc+nx_, ng_, m->lam_g);

    // Copy optimal cost to output
    if (m->f) *m->f = m->fk;

    // Copy optimal constraint values to output
    casadi_copy(m->gk, ng_, m->g);

    // Free memory
    deleteSNOPT(&prob);
  }

  void SnoptInterface::
  userfun(SnoptMemory* m, int* mode, int nnObj, int nnCon, int nnJac, int nnL, int neJac,
          double* x, double* fObj, double*gObj, double* fCon, double* gCon,
          int nState, char* cu, int lencu, int* iu, int leniu, double* ru,
          int lenru) const {
    try {

      casadi_assert_message(nnCon_ == nnCon, "Con " << nnCon_ << " <-> " << nnCon);
      casadi_assert_message(nnObj_ == nnObj, "Obj " << nnObj_ << " <-> " << nnObj);
      casadi_assert_message(nnJac_ == nnJac, "Jac " << nnJac_ << " <-> " << nnJac);

      // Get reduced decision variables
      casadi_fill(m->xk2, nx_, 0.);
      for (int k = 0; k < nnObj; ++k) m->xk2[k] = x[k];

      // Evaluate gradF with the linear variables put to zero
      const double** arg = m->arg;
      *arg++ = m->xk2;
      *arg++ = m->p;
      double** res = m->res;
      *res++ = fObj;
      *res++ = m->jac_fk;
      calc_function(m, "nlp_jac_f");

      // provide nonlinear part of objective gradient to SNOPT
      for (int k = 0; k < nnObj; ++k) {
        int el = jac_f_fcn_.sparsity_out(1).colind(k);
        if (jac_f_fcn_.sparsity_out(1).colind(k+1) > el) {
          gObj[k] = m->jac_fk[el];
        } else {
          gObj[k] = 0;
        }
      }

      jac_f_fcn_.sparsity_out(1).sanity_check(true);
      jac_f_fcn_.sparsity_out(1).sanity_check(false);

      if (!jac_g_fcn_.is_null()) {
        // Get reduced decision variables
        casadi_fill(m->xk2, nx_, 0.);
        for (int k = 0; k < nnJac; ++k) {
          m->xk2[k] = x[k];
        }

        // Evaluate jacG with the linear variabes put to zero
        const double** arg = m->arg;
        *arg++ = m->xk2;
        *arg++ = m->p;
        double** res = m->res;
        *res++ = m->gk;
        *res++ = m->jac_gk;
        calc_function(m, "nlp_jac_g");

        // provide nonlinear part of constraint jacobian to SNOPT
        int kk = 0;
        for (int j = 0; j < nnJac; ++j) {
          for (int k = A_structure_.colind(j); k < A_structure_.sparsity().colind(j+1); ++k) {
            if (A_structure_.row(k) >= nnCon) break;
            int i = A_structure_.nonzeros()[k];
            if (i > 0) {
              gCon[kk++] = m->jac_gk[i-1];
            }
          }
        }

        casadi_assert(kk == 0 || kk == neJac);

        // provide nonlinear part of objective to SNOPT
        for (int k = 0; k < nnCon; ++k) {
          fCon[k] = m->gk[k];
        }
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
    auto m = SnoptMemory::mempool.at(iu[0]);
    m->self.userfun(m, mode, *nnObj, *nnCon, *nJac, *nnL, *neJac,
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
