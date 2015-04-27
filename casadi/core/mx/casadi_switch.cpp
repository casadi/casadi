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


#include "casadi_switch.hpp"
#include "../function/function_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace casadi {

  Switch::Switch(const MX& ind, const std::vector<MX>& arg,
                 const std::vector<Function>& f, const Function& f_def)
    : f_(f), f_def_(f_def) {

    // Consitency check
    casadi_assert(!f_.empty());

    // Make sure ind scalar and dense (constructor should not have been called if sparse)
    casadi_assert(ind.isScalar(true));

    // Get input and output sparsities
    int num_in = -1, num_out=-1;
    for (int k=0; k<=f_.size(); ++k) {
      Function& fk = k<f_.size() ? f_[k] : f_def_;
      if (fk.isNull()) continue;
      fk.assertInit();
      if (num_in<0) {
        // Number of inputs and outputs
        num_in=fk.getNumInputs();
        num_out=fk.getNumOutputs();
        // Output sparsity
        sp_out_.resize(num_out);
        for (int i=0; i<num_out; ++i) sp_out_[i] = fk.output(i).sparsity();
        // Input sparsity
        sp_in_.resize(num_in);
        for (int i=0; i<num_in; ++i) sp_in_[i] = fk.input(i).sparsity();
      } else {
        // Assert matching number of inputs and outputs
        casadi_assert(num_in==fk.getNumInputs());
        casadi_assert(num_out==fk.getNumOutputs());
        // Intersect with output sparsity
        for (int i=0; i<num_out; ++i) {
          sp_out_[i] = sp_out_[i].patternIntersection(fk.output(i).sparsity());
        }
        // Intersect with input sparsity
        for (int i=0; i<num_in; ++i) {
          sp_out_[i] = sp_out_[i].patternIntersection(fk.input(i).sparsity());
        }
      }
    }

    // Illegal to pass only "null" functions
    casadi_assert_message(num_in>=0, "All functions are null");

    // Check number of arguments
    casadi_assert_message(arg.size()==num_in, "Argument list length (" << arg.size()
                          << ") does not match number of inputs (" << num_in << ")");

    // Create arguments of the right dimensions and sparsity
    vector<MX> arg1(num_in+1);
    arg1[0] = ind;
    for (int i=0; i<num_in; ++i) {
      arg1[i+1] = projectArg(arg[i], sp_in_[i], i);
    }
    setDependencies(arg1);
    setSparsity(Sparsity::scalar());
  }

  Switch* Switch::clone() const {
    return new Switch(*this);
  }

  void Switch::evalD(cp_double* arg, p_double* res,
                           int* itmp, double* rtmp) {
    // Get conditional
    int k = static_cast<int>(*arg[0]);

    // Get the function to be evaluated
    Function& fk = k>=0 && k<f_.size() ? f_[k] : f_def_;

    // Input buffers
    for (int i=0; i<sp_in_.size(); ++i) {
      const Sparsity& s = fk.input(i).sparsity();
      casadi_assert_message(s==sp_in_[i], "Not implemented");
    }

    // Output buffers
    for (int i=0; i<sp_out_.size(); ++i) {
      const Sparsity& s = fk.output(i).sparsity();
      casadi_assert_message(s==sp_out_[i], "Not implemented");
    }

    // Evaluate case
    fk->evalD(arg+1, res, itmp, rtmp);
  }

  int Switch::nout() const {
    return sp_out_.size();
  }

  const Sparsity& Switch::sparsity(int oind) const {
    return sp_out_.at(oind);
  }

  void Switch::evalSX(cp_SXElement* arg, p_SXElement* res, int* itmp, SXElement* rtmp) {
    //    fcn_->evalSX(arg, res, itmp, rtmp);
  }

  void Switch::evalMX(const vector<MX>& arg, vector<MX>& res) {
    //    res = fcn_->createCall(arg);
  }

  void Switch::evalFwd(const vector<vector<MX> >& fseed,
                     vector<vector<MX> >& fsens) {
    /*
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Call the cached functions
    fcn_.callForward(arg, res, fseed, fsens);
    */
  }

  void Switch::evalAdj(const vector<vector<MX> >& aseed,
                     vector<vector<MX> >& asens) {
    /*
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Call the cached functions
    vector<vector<MX> > v;
    fcn_.callReverse(arg, res, aseed, v);
    for (int i=0; i<v.size(); ++i) {
      for (int j=0; j<v[i].size(); ++j) {
        if (!v[i][j].isEmpty()) { // TODO(@jaeandersson): Hack
          asens[i][j] += v[i][j];
        }
      }
    }
    */
  }

  void Switch::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    //    fcn_ = deepcopy(fcn_, already_copied);
  }

  void Switch::spFwd(cp_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    //    fcn_.spFwd(arg, res, itmp, rtmp);
  }

  void Switch::spAdj(p_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    //    fcn_.spAdj(arg, res, itmp, rtmp);
  }

  void Switch::generate(const vector<int>& arg, const vector<int>& res,
                      CodeGenerator& g) const {
    // Argument without the case index
    vector<int> arg1(arg.begin()+1, arg.end());

    // Codegen switch
    g.body << "  switch(int(" << g.workel(arg[0], dep(0).nnz()) << ")) {" << endl;
    for (int k=0; k<f_.size(); ++k) {
      g.body << "  case " << k << ": {" << endl;
      f_[k]->generate(arg1, res, g);
      g.body << "  }" << endl;
    }

    // Codegen default case
    g.body << "  default: {" << endl;
    f_def_->generate(arg1, res, g);
    g.body << "  }" << endl;

    // End switch
    g.body << "  }" << endl;
  }

  void Switch::nTmp(size_t& ni, size_t& nr) {
    ni=0;
    nr=0;
    for (int k=0; k<=f_.size(); ++k) {
      Function& fk = k<f_.size() ? f_[k] : f_def_;
      if (fk.isNull()) continue;

      // Get local work vector sizes
      size_t ni_k, nr_k;
      fk.nTmp(ni_k, nr_k);

      // Add size for input buffers
      for (int i=0; i<sp_in_.size(); ++i) {
        const Sparsity& s = fk.input(i).sparsity();
        if (s!=sp_in_[i]) nr_k += s.nnz();
      }

      // Add size for output buffers
      for (int i=0; i<sp_out_.size(); ++i) {
        const Sparsity& s = fk.output(i).sparsity();
        if (s!=sp_out_[i]) nr_k += s.nnz();
      }

      // Find the largest
      ni = max(ni, ni_k);
      nr = max(nr, nr_k);
    }
  }

  Function& Switch::getFunction(int i) {
    if (i<f_.size()) {
      return f_[i];
    } else {
      casadi_assert(i==f_.size());
      return f_def_;
    }
  }

  std::string Switch::printArg(const std::vector<std::string>& arg) {
    stringstream ss;
    ss << "[";
    for (int i=1; i<arg.size(); ++i) {
      if (i>1) ss << ", ";
      ss << arg[i];
    }
    ss << "]";
    return ss.str();
  }

  std::string Switch::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    if (f_.size()==1) {
      // Print as if-then-else
      ss << "if_then_else(" << arg[0] << ", " << printArg(arg) << ", "
         << f_def_.getOption("name") << ", " << f_[0].getOption("name") << ")";
    } else {
      // Print generic
      ss << "conditional(" << arg[0] << ", " << printArg(arg) << ", [";
      for (int k=0; k<f_.size(); ++k) {
        if (k!=0) ss << ", ";
        ss << f_[k].getOption("name");
      }
      ss << "], " << f_def_.getOption("name") << ")";
    }
    return ss.str();
  }

  std::vector<MX> Switch::create(const MX& ind, const std::vector<MX>& arg,
                                 const std::vector<Function>& f, const Function& f_def) {
    return MX::createMultipleOutput(new Switch(ind, arg, f, f_def));
  }

} // namespace casadi
