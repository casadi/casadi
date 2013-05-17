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
#include "../fx/derivative.hpp"

using namespace std;

namespace CasADi {

  EvaluationMX::EvaluationMX(const FX& fcn, std::vector<MX> arg) : fcn_(fcn) {

    // Number inputs and outputs
    int num_in = fcn.getNumInputs();
    casadi_assert(arg.size()<=num_in);

    // Add arguments if needed
    arg.resize(num_in);

    // Replace nulls with zeros of the right dimension
    for(int i=0; i<arg.size(); ++i){
      if(arg[i].isNull()) arg[i] = MX::zeros(fcn_.input(i).sparsity());
    }

    setDependencies(arg);
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
  
    // Number of inputs and outputs
    int num_in = fcn_.getNumInputs();
    int num_out = fcn_.getNumOutputs();

    // Number of derivative directions to calculate
    int nfdir = fsens.size();
    int nadir = aseed.size();

    // Number of derivative directions supported by the function
    int max_nfdir = fcn_.numAllocFwd();
    int max_nadir = fcn_.numAllocAdj();

    // Current forward and adjoint direction
    int offset_nfdir = 0, offset_nadir = 0;

    // Has the function been evaluated once
    bool fcn_evaluated = false;

    // Pass the inputs to the function
    for (int i = 0; i < num_in; ++i) {
      DMatrix *a = arg[i];
      if(a != 0){
        fcn_.setInput(*a, i);
      } else {
        fcn_.setInput(0., i);
      }
    }
  
    // Evaluate until everything has been determinated
    while (!fcn_evaluated || offset_nfdir < nfdir || offset_nadir < nadir) {

      // Number of forward and adjoint directions in the current "batch"
      int nfdir_f_batch = std::min(nfdir - offset_nfdir, max_nfdir);
      int nadir_f_batch = std::min(nadir - offset_nadir, max_nadir);

      // Pass the forward seeds to the function
      for(int d = 0; d < nfdir_f_batch; ++d){
        for(int i = 0; i < num_in; ++i){
          DMatrix *a = fseed[offset_nfdir + d][i];
          if(a != 0){
            fcn_.setFwdSeed(*a, i, d);
          } else {
            fcn_.setFwdSeed(0., i, d);
          }
        }
      }

      // Pass the adjoint seed to the function
      for(int d = 0; d < nadir_f_batch; ++d){
        for(int i = 0; i < num_out; ++i) {
          DMatrix *a = aseed[offset_nadir + d][i];
          if(a != 0){
            fcn_.setAdjSeed(*a, i, d);
          } else {
            fcn_.setAdjSeed(0., i, d);
          }
        }
      }

      // Evaluate
      fcn_.evaluate(nfdir_f_batch, nadir_f_batch);
    
      // Get the outputs if first evaluation
      if(!fcn_evaluated){
        for(int i = 0; i < num_out; ++i) {
          if(res[i] != 0) fcn_.getOutput(*res[i], i);
        }
      }

      // Marked as evaluated
      fcn_evaluated = true;

      // Get the forward sensitivities
      for(int d = 0; d < nfdir_f_batch; ++d){
        for(int i = 0; i < num_out; ++i) {
          DMatrix *a = fsens[offset_nfdir + d][i];
          if(a != 0) fcn_.getFwdSens(*a, i, d);
        }
      }

      // Get the adjoint sensitivities
      for (int d = 0; d < nadir_f_batch; ++d) {
        for (int i = 0; i < num_in; ++i) {
          DMatrix *a = asens[offset_nadir + d][i];
          if(a != 0){
            a->sparsity().add(a->ptr(),fcn_.adjSens(i,d).ptr(),fcn_.adjSens(i,d).sparsity());
          }
        }
      }

      // Update direction offsets
      offset_nfdir += nfdir_f_batch;
      offset_nadir += nadir_f_batch;
    }

    // Clear adjoint seeds
    clearVector(aseed);
  }

  int EvaluationMX::getNumOutputs() const {
    return fcn_.getNumOutputs();
  }

  const CRSSparsity& EvaluationMX::sparsity(int oind) const{
    return fcn_.output(oind).sparsity();
  }

  FX& EvaluationMX::getFunction() {
    return fcn_;
  }

  void EvaluationMX::evaluateSX(const SXMatrixPtrV& arg, SXMatrixPtrV& res,
                                const SXMatrixPtrVV& fseed, SXMatrixPtrVV& fsens,
                                const SXMatrixPtrVV& aseed, SXMatrixPtrVV& asens) {
  
    // Create input arguments
    vector<SXMatrix> argv(arg.size());
    for(int i=0; i<arg.size(); ++i){
      argv[i] = SXMatrix(fcn_.input(i).sparsity(),0.);
      if(arg[i] != 0)
        argv[i].set(*arg[i]);
    }

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
    int nfdir = fsens.size();
    casadi_assert(nfdir==0 || fcn_.spCanEvaluate(true));
    int nadir = aseed.size();
    casadi_assert(nadir==0 || fcn_.spCanEvaluate(false));

    // Get/generate the derivative function
    FX d = fcn_.derivative(nfdir, nadir);

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
      for (MXPtrV::const_iterator j = i->begin(); j != i->end(); ++j){
        if(*j!=0) **j = MX();
      }
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

  void EvaluationMX::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    fcn_ = deepcopy(fcn_, already_copied);
  }

  void EvaluationMX::propagateSparsity(DMatrixPtrV& arg, DMatrixPtrV& res, bool use_fwd) {
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
          if(!use_fwd) fill_n(get_bvec_t(res[oind]->data()),res[oind]->size(),bvec_t(0));
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
            arg[iind]->sparsity().bor(
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
      if(!use_fwd){
        for(int oind=0; oind<res.size(); ++oind){
          if(res[oind]!=0){
            vector<double> &w = res[oind]->data();
            fill_n(get_bvec_t(w),w.size(),bvec_t(0));
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
    int nfdir = fseed.size();
    int nadir = aseed.size();

    // Create the evaluation node
    MX ev;
    if(nfdir>0 || nadir>0){
      // Create derivative function
      Derivative dfcn(fcn,nfdir,nadir);
      stringstream ss;
      ss << "der_" << fcn.getOption("name") << "_" << nfdir << "_" << nadir;
      dfcn.setOption("verbose",fcn.getOption("verbose"));
      dfcn.setOption("name",ss.str());
      dfcn.init();
    
      // All inputs
      vector<MX> darg;
      darg.reserve(num_in*(1+nfdir) + num_out*nadir);
      darg.insert(darg.end(),arg.begin(),arg.end());
    
      // Forward seeds
      for(int dir=0; dir<nfdir; ++dir){
        darg.insert(darg.end(),fseed[dir].begin(),fseed[dir].end());
      }
    
      // Adjoint seeds
      for(int dir=0; dir<nadir; ++dir){
        darg.insert(darg.end(),aseed[dir].begin(),aseed[dir].end());
      }
    
      ev.assignNode(new EvaluationMX(dfcn, darg));
    } else {
      ev.assignNode(new EvaluationMX(fcn, arg));
    }

    // Output index
    int ind = 0;

    // Create the output nodes corresponding to the nondifferented function
    res.resize(num_out);
    for (int i = 0; i < num_out; ++i, ++ind) {
      if(!output_given){
        if(!fcn.output(i).empty()){
          res[i].assignNode(new OutputNode(ev, ind));
        } else {
          res[i] = MX();
        }
      }
    }

    // Forward sensitivities
    fsens.resize(nfdir);
    for(int dir = 0; dir < nfdir; ++dir){
      fsens[dir].resize(num_out);
      for (int i = 0; i < num_out; ++i, ++ind) {
        if (!fcn.output(i).empty()){
          fsens[dir][i].assignNode(new OutputNode(ev, ind));
        } else {
          fsens[dir][i] = MX();
        }
      }
    }

    // Adjoint sensitivities
    asens.resize(nadir);
    for (int dir = 0; dir < nadir; ++dir) {
      asens[dir].resize(num_in);
      for (int i = 0; i < num_in; ++i, ++ind) {
        if (!fcn.input(i).empty()) {
          asens[dir][i].assignNode(new OutputNode(ev, ind));
        } else {
          asens[dir][i] = MX();
        }
      }
    }
  }

  void EvaluationMX::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
  
    // Running index of the temporary used
    int nr=0;

    // Copy arguments with nonmatching sparsities to the temp vector
    vector<string> arg_mod = arg;
    for(int i=0; i<fcn_.getNumInputs(); ++i){
      if(dep(i).sparsity()!=fcn_.input(i).sparsity()){
        arg_mod[i] = "rrr+" + CodeGenerator::numToString(nr);
        nr += fcn_.input(i).size();
        
        // Codegen "copy sparse"
        gen.addAuxiliary(CodeGenerator::AUX_COPY_SPARSE);
        
        int sp_arg = gen.getSparsity(dep(i).sparsity());
        int sp_input = gen.addSparsity(fcn_.input(i).sparsity());
        stream << "  casadi_copy_sparse(" << arg[i] << ",s" << sp_arg << "," << arg_mod[i] << ",s" << sp_input << ");" << std::endl;
      }
    }

    // Get the index of the function
    int f = gen.getDependency(fcn_);
    stream << "  f" << f << "(";
  
    // Pass inputs to the function input buffers
    for(int i=0; i<arg.size(); ++i){
      stream << arg_mod.at(i);
      if(i+1<arg.size()+res.size()) stream << ",";
    }

    // Separate arguments and results with an extra space
    stream << " ";

    // Pass results to the function input buffers
    for(int i=0; i<res.size(); ++i){
      stream << res.at(i);
      if(i+1<res.size()) stream << ",";
    }
  
    // Finalize the function call
    stream << ");" << endl;  
  }
  
  void EvaluationMX::nTmp(size_t& ni, size_t& nr){
    // Start with no extra memory
    ni=0;
    nr=0;

    // Add memory for all inputs with nonmatching sparsity
    for(int i=0; i<fcn_.getNumInputs(); ++i){
      if(dep(i).isNull() || dep(i).sparsity()!=fcn_.input(i).sparsity()){
        nr += fcn_.input(i).size();
      }
    }
  }

} // namespace CasADi
