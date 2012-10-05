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

#include "evaluation_mx.hpp"
#include "../fx/fx_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace CasADi {

EvaluationMX::EvaluationMX(const FX& fcn, const std::vector<MX> &arg,
    const std::vector<std::vector<MX> > &fseed,
    const std::vector<std::vector<MX> > &aseed, bool output_given) :
    fcn_(fcn), output_given_f_(output_given) {

  // Number inputs and outputs
  int num_in = fcn.getNumInputs();
  int num_out = fcn.getNumOutputs();

  // Number of directional derivatives
  nfwd_f_ = fseed.size();
  nadj_f_ = aseed.size();

  // All dependencies of the function
  vector<MX> d;
  d.reserve(num_in + nfwd_f_ * num_in + nadj_f_ * num_out);

  // Add arguments
  d.insert(d.end(), arg.begin(), arg.end());
  d.resize(d.size() + num_in - arg.size()); // Add empty matrices

  // Add forward seeds
  for (int dir = 0; dir < nfwd_f_; ++dir) {
    d.insert(d.end(), fseed[dir].begin(), fseed[dir].end());
    d.resize(d.size() + num_in - fseed[dir].size()); // Add empty matrices
  }

  // Add adjoint seeds
  for (int dir = 0; dir < nadj_f_; ++dir) {
    d.insert(d.end(), aseed[dir].begin(), aseed[dir].end());
    d.resize(d.size() + num_out - aseed[dir].size()); // Add empty matrices
  }

  setDependencies(d);
  setSparsity(CRSSparsity(1, 1, true));
}

EvaluationMX* EvaluationMX::clone() const {
  return new EvaluationMX(*this);
}

void EvaluationMX::printPart(std::ostream &stream, int part) const {
  if (part == 0) {
    stream << fcn_ << ".call([";
  } else if (part == ndep()) {
    stream << "])";
  } else {
    stream << ",";
  }
}

void EvaluationMX::evaluateD(const DMatrixPtrV& arg, DMatrixPtrV& res,
    const DMatrixPtrVV& fseed, DMatrixPtrVV& fsens,
    const DMatrixPtrVV& aseed, DMatrixPtrVV& asens) {
  
  // Check if the function is differentiated
  bool is_diff = nfwd_f_ > 0 || nadj_f_ > 0;

  // Number of inputs and outputs
  int num_in = fcn_.getNumInputs();
  int num_out = fcn_.getNumOutputs();

  // Number of derivative directions to calculate
  int nfwd = is_diff ? nfwd_f_ : fsens.size();
  int nadj = is_diff ? nadj_f_ : aseed.size();

  // Number of derivative directions supported by the function
  int max_nfwd = fcn_->nfdir_;
  int max_nadj = fcn_->nadir_;

  // Current forward and adjoint direction
  int offset_fwd = 0, offset_adj = 0;

  // Has the function been evaluated once
  bool fcn_evaluated = false;

  // Pass the inputs to the function
  for (int i = 0; i < num_in; ++i) {
    if (arg[i] != 0 && arg[i]->size() != 0) {
      fcn_.setInput(*arg[i], i);
    }
  }

  // Evaluate until everything has been determinated
  while (!fcn_evaluated || offset_fwd < nfwd || offset_adj < nadj) {

    // Number of forward and adjoint directions in the current "batch"
    int nfwd_f_batch = std::min(nfwd - offset_fwd, max_nfwd);
    int nadj_f_batch = std::min(nadj - offset_adj, max_nadj);

    // Pass the forward seeds to the function
    for (int d = 0; d < nfwd_f_batch; ++d) {
      int dir = offset_fwd + d;
      for (int i = 0; i < num_in; ++i) {
        DMatrix *a =
            is_diff ?
                arg[num_in + dir * num_in + i] : fseed[dir][i];
        if (a != 0 && a->size() != 0) {
          fcn_.setFwdSeed(a->data(), i, d);
        }
      }
    }

    // Pass the adjoint seed to the function
    for (int d = 0; d < nadj_f_batch; ++d) {
      int dir = offset_adj + d;
      for (int i = 0; i < num_out; ++i) {
        DMatrix *a =
            is_diff ?
                arg[num_in + nfwd * num_in + dir * num_out + i] :
                aseed[dir][i];
        if (a != 0 && a->size() != 0) {
          fcn_.setAdjSeed(a->data(), i, d);
        }
      }
    }

    // Evaluate
    fcn_.evaluate(nfwd_f_batch, nadj_f_batch);

    // Get the outputs if first evaluation
    if (!fcn_evaluated) {
      for (int i = 0; i < num_out; ++i) {
        if (res[i] != 0 && res[i]->size() != 0) {
          fcn_.getOutput(res[i]->data(), i);
        }
      }
    }

    // Marked as evaluated
    fcn_evaluated = true;

    // Get the forward sensitivities
    for (int d = 0; d < nfwd_f_batch; ++d) {
      int dir = offset_fwd + d;
      for (int i = 0; i < num_out; ++i) {
        DMatrix *a =
            is_diff ?
                res[num_out + dir * num_out + i] :
                fsens[dir][i];
        if (a != 0 && a->size() != 0) {
          fcn_.getFwdSens(a->data(), i, d);
        }
      }
    }

    // Get the adjoint sensitivities
    for (int d = 0; d < nadj_f_batch; ++d) {
      int dir = offset_adj + d;
      for (int i = 0; i < num_in; ++i) {
        DMatrix *a =
            is_diff ?
                res[num_out + nfwd * num_out + dir * num_in + i] :
                asens[dir][i];
        if (a != 0 && a->size() != 0) {
          transform(a->begin(), a->end(), fcn_.adjSens(i, d).begin(),
              a->begin(), plus<double>());
        }
      }
    }

    // Update direction offsets
    offset_fwd += nfwd_f_batch;
    offset_adj += nadj_f_batch;
  }
}

int EvaluationMX::getNumOutputs() const {
  int num_in = fcn_.getNumInputs();
  int num_out = fcn_.getNumOutputs();
  return num_out + nfwd_f_ * num_out + nadj_f_ * num_in;
}

const CRSSparsity& EvaluationMX::sparsity(int oind) {
  int num_in = fcn_.getNumInputs();
  int num_out = fcn_.getNumOutputs();
  int num_out_and_fsens = num_out + num_out * nfwd_f_; // number of outputs _and_ forward sensitivities combined
  if (oind < num_out_and_fsens) { // Output or forward sensitivity
    return fcn_.output(oind % num_out).sparsity();
  } else { // Adjoint sensitivity
    return fcn_.input((oind - num_out_and_fsens) % num_in).sparsity();
  }
}

FX& EvaluationMX::getFunction() {
  return fcn_;
}

void EvaluationMX::evaluateSX(const SXMatrixPtrV& arg, SXMatrixPtrV& res,
    const SXMatrixPtrVV& fseed, SXMatrixPtrVV& fsens,
    const SXMatrixPtrVV& aseed, SXMatrixPtrVV& asens) {
  // Create input arguments
  vector<SXMatrix> argv = getVector(arg);

  // Evaluate symbolically
  vector<SXMatrix> resv = fcn_.eval(argv);

  // Collect the result
  for (int i = 0; i < res.size(); ++i) {
    if (res[i] != 0)
      *res[i] = resv[i];
  }
}

void EvaluationMX::evaluateMX(const MXPtrV& arg, MXPtrV& res, const MXPtrVV& fseed, MXPtrVV& fsens, const MXPtrVV& aseed, MXPtrVV& asens, bool output_given) {
  
  // Number of sensitivity directions
  int nfwd = fsens.size();
  casadi_assert(nfwd==0 || fcn_.spCanEvaluate(true));
  int nadj = aseed.size();
  casadi_assert(nadj==0 || fcn_.spCanEvaluate(false));

  // Get/generate the derivative function
  FX d = fcn_.derivative(nfwd, nadj);

  // Temporary
  vector<MX> tmp;

  // Assemble inputs
  vector<MX> d_arg;
  d_arg.reserve(d.getNumInputs());

  // Nondifferentiated inputs
  tmp = getVector(arg);
  d_arg.insert(d_arg.end(), tmp.begin(), tmp.end());
  for (MXPtrVV::const_iterator i = fseed.begin(); i != fseed.end(); ++i) {
    tmp = getVector(*i);
    d_arg.insert(d_arg.end(), tmp.begin(), tmp.end());
  }
  for (MXPtrVV::const_iterator i = aseed.begin(); i != aseed.end(); ++i) {
    tmp = getVector(*i);
    d_arg.insert(d_arg.end(), tmp.begin(), tmp.end());
  }

  // Evaluate symbolically
  vector<MX> d_res = d.call(d_arg);
  vector<MX>::const_iterator d_res_it = d_res.begin();

  // Collect the nondifferentiated results
  for (MXPtrV::iterator i = res.begin(); i != res.end(); ++i, ++d_res_it) {
    if (!output_given && *i) **i = *d_res_it;
  }

  // Collect the forward sensitivities
  for (MXPtrVV::iterator j = fsens.begin(); j != fsens.end(); ++j) {
    for (MXPtrV::iterator i = j->begin(); i != j->end(); ++i, ++d_res_it) {
      if (*i) **i = *d_res_it;
    }
  }

  // Collect the adjoint sensitivities
  for (MXPtrVV::iterator j = asens.begin(); j != asens.end(); ++j) {
    for (MXPtrV::iterator i = j->begin(); i != j->end(); ++i, ++d_res_it) {
      if(*i && !d_res_it->isNull()){
        **i += *d_res_it;
      }
    }
  }

  // Make sure that we've got to the end of the outputs
  casadi_assert(d_res_it==d_res.end());
}

void EvaluationMX::deepCopyMembers(
  std::map<SharedObjectNode*, SharedObject>& already_copied) {
  MXNode::deepCopyMembers(already_copied);
  fcn_ = deepcopy(fcn_, already_copied);
}

void EvaluationMX::propagateSparsity(DMatrixPtrV& arg, DMatrixPtrV& res,bool use_fwd) {
  if (fcn_.spCanEvaluate(use_fwd)) {
    // Propagating sparsity pattern supported
    
    // Pass/clear forward seeds/adjoint sensitivities
    for (int iind = 0; iind < fcn_.getNumInputs(); ++iind) {
      // Input vector
      vector<double> &v = fcn_.input(iind).data();
      if (v.empty()) continue; // FIXME: remove?

      if (arg[iind] == 0) {
        // Set to zero if not used
        fill_n(get_bvec_t(v), v.size(), bvec_t(0));
      } else {
        // Copy output
        fcn_.input(iind).sparsity().set(
          get_bvec_t(fcn_.input(iind).data()),
          get_bvec_t(arg[iind]->data()),
          arg[iind]->sparsity());
      }
    }

    // Pass/clear adjoint seeds/forward sensitivities
    for (int oind = 0; oind < fcn_.getNumOutputs(); ++oind) {
      // Output vector
      vector<double> &v = fcn_.output(oind).data();
      if (v.empty()) continue; // FIXME: remove?

      if (res[oind] == 0) {
        // Set to zero if not used
        fill_n(get_bvec_t(v), v.size(), bvec_t(0));
      } else {
        // Copy output
        fcn_.output(oind).sparsity().set(
            get_bvec_t(fcn_.output(oind).data()),
            get_bvec_t(res[oind]->data()),
            res[oind]->sparsity());
      }
    }

    // Propagate seedsfcn_.
    fcn_.spInit(use_fwd); // NOTE: should only be done once
    fcn_.spEvaluate(use_fwd);

    // Get the sensitivities
    if (use_fwd) {
      for (int oind = 0; oind < res.size(); ++oind) {
        if (res[oind] != 0) {
          res[oind]->sparsity().set(
              get_bvec_t(res[oind]->data()),
              get_bvec_t(fcn_.output(oind).data()),
              fcn_.output(oind).sparsity());
        }
      }
    } else {
      for (int iind = 0; iind < arg.size(); ++iind) {
        if (arg[iind] != 0) {
          arg[iind]->sparsity().set(
              get_bvec_t(arg[iind]->data()),
              get_bvec_t(fcn_.input(iind).data()),
              fcn_.input(iind).sparsity());
        }
      }
    }

    // Clear seeds and sensitivities
    for (int iind = 0; iind < arg.size(); ++iind) {
      vector<double> &v = fcn_.input(iind).data();
      fill(v.begin(), v.end(), 0);
    }
    for (int oind = 0; oind < res.size(); ++oind) {
      vector<double> &v = fcn_.output(oind).data();
      fill(v.begin(), v.end(), 0);
    }

  } else {
    // Propagating sparsity pattern not supported

    if (use_fwd) {
      // Clear the outputs
      for (int oind = 0; oind < res.size(); ++oind) {
        // Skip of not used
        if (res[oind] == 0)
          continue;

        // Get data array for output and clear it
        bvec_t *outputd = get_bvec_t(res[oind]->data());
        fill_n(outputd, res[oind]->size(), 0);
      }
    }

    // Loop over inputs
    for (int iind = 0; iind < arg.size(); ++iind) {
      // Skip of not used
      if (arg[iind] == 0)
        continue;

      // Skip if no seeds
      if (use_fwd && arg[iind]->empty())
        continue;

      // Get data array for input
      bvec_t *inputd = get_bvec_t(arg[iind]->data());

      // Loop over outputs
      for (int oind = 0; oind < res.size(); ++oind) {

        // Skip of not used
        if (res[oind] == 0)
          continue;

        // Skip if no seeds
        if (!use_fwd && res[oind]->empty())
          continue;

        // Get the sparsity of the Jacobian block
        CRSSparsity& sp = fcn_.jacSparsity(iind, oind, true);
        if (sp.isNull() || sp.size() == 0)
          continue; // Skip if zero
        const int d1 = sp.size1();
        //const int d2 = sp.size2();
        const vector<int>& rowind = sp.rowind();
        const vector<int>& col = sp.col();

        // Get data array for output
        bvec_t *outputd = get_bvec_t(res[oind]->data());

        // Carry out the sparse matrix-vector multiplication
        for (int i = 0; i < d1; ++i) {
          for (int el = rowind[i]; el < rowind[i + 1]; ++el) {
            // Get column
            int j = col[el];

            // Propagate dependencies
            if (use_fwd) {
              outputd[i] |= inputd[j];
            } else {
              inputd[j] |= outputd[i];
            }
          }
        }
      }
    }
  }
}

void EvaluationMX::create(const FX& fcn, const std::vector<MX> &arg,
    std::vector<MX> &res, const std::vector<std::vector<MX> > &fseed,
    std::vector<std::vector<MX> > &fsens,
    const std::vector<std::vector<MX> > &aseed,
    std::vector<std::vector<MX> > &asens, bool output_given) {

  // Number inputs and outputs
  int num_in = fcn.getNumInputs();
  int num_out = fcn.getNumOutputs();

  // Number of directional derivatives
  int nfwd = fseed.size();
  int nadj = aseed.size();

  // Create the evaluation node
  MX ev;
  ev.assignNode(new EvaluationMX(fcn, arg, fseed, aseed, output_given));

  // Output index
  int ind = 0;

  // Create the output nodes corresponding to the nondifferented function
  res.resize(num_out);
  for (int i = 0; i < num_out; ++i, ++ind) {
    if (fcn.output(i).size1() > 0 && fcn.output(i).size2() > 0) {
      res[i].assignNode(new OutputNode(ev, ind));
    } else {
      res[i] = MX();
    }
  }

  // Forward sensitivities
  fsens.resize(nfwd);
  for (int dir = 0; dir < nfwd; ++dir) {
    fsens[dir].resize(num_out);
    for (int i = 0; i < num_out; ++i, ++ind) {
      if (fcn.output(i).size1() > 0 && fcn.output(i).size2() > 0) {
        fsens[dir][i].assignNode(new OutputNode(ev, ind));
      } else {
        fsens[dir][i] = MX();
      }
    }
  }

  // Adjoint sensitivities
  asens.resize(nadj);
  for (int dir = 0; dir < nadj; ++dir) {
    asens[dir].resize(num_in);
    for (int i = 0; i < num_in; ++i, ++ind) {
      if (fcn.output(i).size1() > 0 && fcn.output(i).size2() > 0) {
        asens[dir][i].assignNode(new OutputNode(ev, ind));
      } else {
        asens[dir][i] = MX();
      }
    }
  }
}

} // namespace CasADi
