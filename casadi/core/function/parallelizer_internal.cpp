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


#include "parallelizer_internal.hpp"
#include "mx_function.hpp"
#include <algorithm>
#ifdef WITH_OPENMP
#include <omp.h>
#endif //WITH_OPENMP

using namespace std;

namespace casadi {

  ParallelizerInternal::ParallelizerInternal(const Function& f, int n) : f_(f), n_(n) {
    addOption("parallelization", OT_STRING, "serial", "", "serial|openmp|mpi");
  }

  ParallelizerInternal::~ParallelizerInternal() {
  }

  void ParallelizerInternal::init() {
    // Get mode
    if (getOption("parallelization")=="serial") {
      mode_ = SERIAL;
    } else if (getOption("parallelization")=="openmp") {
      mode_ = OPENMP;
    } else if (getOption("parallelization")=="mpi") {
      mode_ = MPI;
    } else {
      casadi_error("Parallelization mode " << getOption("parallelization") << " unknown.");
    }

    // Switch to serial mode if OPENMP is not supported
#ifndef WITH_OPENMP
    if (mode_ == OPENMP) {
      casadi_warning("OpenMP parallelization is not available, switching to serial mode. "
                     "Recompile CasADi setting the option WITH_OPENMP to ON.");
      mode_ = SERIAL;
    }
#endif // WITH_OPENMP

    // Initialize the function
    f_.init();

    // Inputs
    int f_num_in = f_.getNumInputs();
    setNumInputs(n_ * f_num_in);
    int k=0;
    for (int i=0; i<n_; ++i)
      for (int j=0; j<f_num_in; ++j)
        input(k++) = f_.input(j);

    // Outputs
    int f_num_out = f_.getNumOutputs();
    setNumOutputs(n_ * f_num_out);
    k=0;
    for (int i=0; i<n_; ++i)
      for (int j=0; j<f_num_out; ++j)
        output(k++) = f_.output(j);

    // Clear the indices
    inind_.clear();   inind_.push_back(0);
    outind_.clear();  outind_.push_back(0);

    // Add the inputs and outputs
    for (int i=0; i<n_; ++i) {
      inind_.push_back(inind_.back()+f_.getNumInputs());
      outind_.push_back(outind_.back()+f_.getNumOutputs());
    }


    // Call the init function of the base class
    FunctionInternal::init();

    // Allocate work vectors
    size_t ni, nr;
    f_.nTmp(ni, nr);
    if (mode_ == OPENMP) {
      itmp_.resize(n_ * ni);
      rtmp_.resize(n_ * nr);
    } else {
      itmp_.resize(ni);
      rtmp_.resize(nr);
    }
  }

  void ParallelizerInternal::evaluate() {
    // Let the first call (which may contain memory allocations) be serial when using OpenMP
    if (mode_== SERIAL) {
      for (int task=0; task<n_; ++task) {
        evaluateTask(task);
      }
    } else if (mode_== OPENMP) {
#ifdef WITH_OPENMP
#pragma omp parallel for
      for (int task=0; task<n_; ++task) {
        evaluateTask(task);
      }
#endif //WITH_OPENMP
#ifndef WITH_OPENMP
      casadi_error("ParallelizerInternal::evaluate: OPENMP support was not available "
                   "during CasADi compilation");
#endif //WITH_OPENMP
    } else if (mode_ == MPI) {
      casadi_error("ParallelizerInternal::evaluate: MPI not implemented");
    }
  }

  void ParallelizerInternal::evaluateTask(int task) {
    // Copy inputs to functions
    int f_num_in = f_.getNumInputs();
    for (int j=0; j<f_num_in; ++j) {
      f_.input(j).set(input(task*f_num_in + j));
    }

    // Evaluate
    f_.evaluate();

    // Get the results
    int f_num_out = f_.getNumOutputs();
    for (int j=0; j<f_num_out; ++j) {
      f_.output(j).get(output(task*f_num_out + j));
    }
  }

  Sparsity ParallelizerInternal::getJacSparsity(int iind, int oind, bool symmetric) {
    int f_num_in = f_.getNumInputs();
    int f_num_out = f_.getNumOutputs();

    // Get the local sparsity patterm
    if (iind / f_num_in == oind / f_num_out) {
      // Same task
      return f_.jacSparsity(iind % f_num_in, oind % f_num_out);
    } else {
      // Different tasks: All-zero jacobian
      return Sparsity();
    }
  }

  Function ParallelizerInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
    int f_num_in = f_.getNumInputs();
    int f_num_out = f_.getNumOutputs();

    // Get the local sparsity patterm
    if (iind / f_num_in == oind / f_num_out) {
      // Same task
      return f_.jacobian(iind % f_num_in, oind % f_num_out, compact, symmetric);
    } else {
      // Different tasks: All-zero jacobian
      return Function();
    }
  }

  void ParallelizerInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    FunctionInternal::deepCopyMembers(already_copied);
    f_ = deepcopy(f_, already_copied);
  }

  void ParallelizerInternal::spInit(bool use_fwd) {
    f_.spInit(use_fwd);
  }

  void ParallelizerInternal::spEvaluate(bool use_fwd) {
    // This function can be parallelized. Move logic in "evaluate" to a template function.
    for (int task=0; task<n_; ++task) {
      spEvaluateTask(use_fwd, task);
    }
  }

  void ParallelizerInternal::spEvaluateTask(bool use_fwd, int task) {
    if (use_fwd) {
      // Set input influence
      for (int j=inind_[task]; j<inind_[task+1]; ++j) {
        int nv = input(j).nnz();
        const bvec_t* p_v = get_bvec_t(input(j).data());
        bvec_t* f_v = get_bvec_t(f_.input(j-inind_[task]).data());
        copy(p_v, p_v+nv, f_v);
      }

      // Propagate
      f_.spEvaluate(use_fwd);

      // Get output dependence
      for (int j=outind_[task]; j<outind_[task+1]; ++j) {
        int nv = output(j).nnz();
        bvec_t* p_v = get_bvec_t(output(j).data());
        const bvec_t* f_v = get_bvec_t(f_.output(j-outind_[task]).data());
        copy(f_v, f_v+nv, p_v);
      }

    } else {

      // Set output influence
      for (int j=outind_[task]; j<outind_[task+1]; ++j) {
        int nv = output(j).nnz();
        const bvec_t* p_v = get_bvec_t(output(j).data());
        bvec_t* f_v = get_bvec_t(f_.output(j-outind_[task]).data());
        copy(p_v, p_v+nv, f_v);
      }

      // Propagate
      f_.spEvaluate(use_fwd);

      // Get input dependence
      for (int j=inind_[task]; j<inind_[task+1]; ++j) {
        int nv = input(j).nnz();
        bvec_t* p_v = get_bvec_t(input(j).data());
        const bvec_t* f_v = get_bvec_t(f_.input(j-inind_[task]).data());
        copy(f_v, f_v+nv, p_v);
      }
    }
  }

  Function ParallelizerInternal::getDerForward(int nfwd) {
    // Generate derivative expressions
    Function der_f = f_.derForward(nfwd);

    // Create a new parallelizer for the derivatives
    Parallelizer par(der_f, n_);

    // Set options and initialize
    par.setOption(dictionary());
    par.init();

    // Create a function call to the parallelizer
    vector<MX> par_arg = par.symbolicInput();
    vector<MX> par_res = par(par_arg);

    // Get the offsets: copy to allow being used as a counter
    std::vector<int> par_inind = par->inind_;
    std::vector<int> par_outind = par->outind_;

    // Arguments and results of the return function
    vector<MX> ret_arg, ret_res;
    ret_arg.reserve(par_arg.size());
    ret_res.reserve(par_res.size());

    // Nondifferentiated inputs and outputs
    for (int i=0; i<n_; ++i) {
      for (int j=inind_[i]; j<inind_[i+1]; ++j) ret_arg.push_back(par_arg[par_inind[i]++]);
      for (int j=outind_[i]; j<outind_[i+1]; ++j) ret_arg.push_back(par_arg[par_inind[i]++]);
    }

    // Forward seeds/sensitivities
    for (int dir=0; dir<nfwd; ++dir) {
      for (int i=0; i<n_; ++i) {
        for (int j=inind_[i]; j<inind_[i+1]; ++j) ret_arg.push_back(par_arg[par_inind[i]++]);
        for (int j=outind_[i]; j<outind_[i+1]; ++j) ret_res.push_back(par_res[par_outind[i]++]);
      }
    }

    // Assemble the return function
    return MXFunction(ret_arg, ret_res);
  }

  Function ParallelizerInternal::getDerReverse(int nadj) {
    // Generate derivative expressions
    Function der_f = f_.derReverse(nadj);

    // Create a new parallelizer for the derivatives
    Parallelizer par(der_f, n_);

    // Set options and initialize
    par.setOption(dictionary());
    par.init();

    // Create a function call to the parallelizer
    vector<MX> par_arg = par.symbolicInput();
    vector<MX> par_res = par(par_arg);

    // Get the offsets: copy to allow being used as a counter
    std::vector<int> par_inind = par->inind_;
    std::vector<int> par_outind = par->outind_;

    // Arguments and results of the return function
    vector<MX> ret_arg, ret_res;
    ret_arg.reserve(par_arg.size());
    ret_res.reserve(par_res.size());

    // Nondifferentiated inputs and outputs
    for (int i=0; i<n_; ++i) {
      for (int j=inind_[i]; j<inind_[i+1]; ++j) ret_arg.push_back(par_arg[par_inind[i]++]);
      for (int j=outind_[i]; j<outind_[i+1]; ++j) ret_arg.push_back(par_arg[par_inind[i]++]);
    }

    // Adjoint seeds/sensitivities
    for (int dir=0; dir<nadj; ++dir) {
      for (int i=0; i<n_; ++i) {
        for (int j=outind_[i]; j<outind_[i+1]; ++j) ret_arg.push_back(par_arg[par_inind[i]++]);
        for (int j=inind_[i]; j<inind_[i+1]; ++j) ret_res.push_back(par_res[par_outind[i]++]);
      }
    }

    // Assemble the return function
    return MXFunction(ret_arg, ret_res);
  }

} // namespace casadi

