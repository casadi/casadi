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


#include "switch_internal.hpp"
#include "mx_function.hpp"

using namespace std;

namespace casadi {

  SwitchInternal::SwitchInternal(const std::string& name,
                                 const std::vector<Function>& f, const Function& f_def)
    : FunctionInternal(name), f_(f), f_def_(f_def) {

    // Consitency check
    casadi_assert(!f_.empty());
  }

  SwitchInternal::~SwitchInternal() {
  }

  void SwitchInternal::init() {
    // Initialize the functions, get input and output sparsities
    // Input and output sparsities
    std::vector<Sparsity> sp_in, sp_out;
    int num_in = -1, num_out=-1;
    for (int k=0; k<=f_.size(); ++k) {
      Function& fk = k<f_.size() ? f_[k] : f_def_;
      if (fk.isNull()) continue;
      fk.init(false);
      if (num_in<0) {
        // Number of inputs and outputs
        num_in=fk.n_in();
        num_out=fk.n_out();
        // Output sparsity
        sp_out.resize(num_out);
        for (int i=0; i<num_out; ++i) sp_out[i] = fk.output(i).sparsity();
        // Input sparsity
        sp_in.resize(num_in);
        for (int i=0; i<num_in; ++i) sp_in[i] = fk.input(i).sparsity();
      } else {
        // Assert matching number of inputs and outputs
        casadi_assert(num_in==fk.n_in());
        casadi_assert(num_out==fk.n_out());
        // Intersect with output sparsity
        for (int i=0; i<num_out; ++i) {
          sp_out[i] = sp_out[i].patternIntersection(fk.output(i).sparsity());
        }
        // Intersect with input sparsity
        for (int i=0; i<num_in; ++i) {
          sp_in[i] = sp_in[i].patternIntersection(fk.input(i).sparsity());
        }
      }
    }

    // Illegal to pass only "null" functions
    casadi_assert_message(num_in>=0, "All functions are null");

    // Allocate input and output buffers
    ibuf_.resize(1+num_in);
    input(0) = 0; // conditional
    for (int i=0; i<num_in; ++i) input(i+1) = DMatrix::zeros(sp_in[i]);
    obuf_.resize(num_out);
    for (int i=0; i<num_out; ++i) output(i) = DMatrix::zeros(sp_out[i]);

    // Call the initialization method of the base class
    FunctionInternal::init();

    // Get required work
    for (int k=0; k<=f_.size(); ++k) {
      const Function& fk = k<f_.size() ? f_[k] : f_def_;
      if (fk.isNull()) continue;

      // Get local work vector sizes
      alloc(fk);
      size_t sz_w = fk.sz_w();

      // Add size for input buffers
      for (int i=1; i<n_in(); ++i) {
        const Sparsity& s = fk.input(i-1).sparsity();
        if (s!=input(i).sparsity()) sz_w += s.nnz();
      }

      // Add size for output buffers
      for (int i=0; i<n_out(); ++i) {
        const Sparsity& s = fk.output(i).sparsity();
        if (s!=output(i).sparsity()) sz_w += s.nnz();
      }

      // Make sure enough work for this
      alloc_w(sz_w);
    }
  }

  void SwitchInternal::evalD(const double** arg, double** res, int* iw, double* w) {
    // Get conditional
    int k = static_cast<int>(*arg[0]);

    // Get the function to be evaluated
    Function& fk = k<f_.size() ? f_[k] : f_def_;

    // Input buffers
    for (int i=1; i<n_in(); ++i) {
      const Sparsity& s = fk.input(i-1).sparsity();
      casadi_assert_message(s==input(i).sparsity(), "Not implemented");
    }

    // Output buffers
    for (int i=0; i<n_out(); ++i) {
      const Sparsity& s = fk.output(i).sparsity();
      casadi_assert_message(s==output(i).sparsity(), "Not implemented");
    }

    // Evaluate the corresponding function
    fk->eval(arg+1, res, iw, w);
  }

  Function SwitchInternal
  ::getDerForward(const std::string& name, int nfwd, Dict& opts) {
    // Derivative of each case
    vector<Function> der(f_.size());
    for (int k=0; k<f_.size(); ++k) {
      if (!f_[k].isNull()) der[k] = f_[k].derForward(nfwd);
    }

    // Default case
    Function der_def;
    if (!f_def_.isNull()) der_def = f_def_.derForward(nfwd);

    // New Switch for derivatives
    stringstream ss;
    ss << "fwd" << nfwd << "_" << name_;
    Switch sw(ss.str(), der, der_def);

    // Construct wrapper inputs and arguments for calling sw
    vector<MX> arg = mx_in();
    vector<MX> res = symbolicOutput();
    vector<vector<MX> > seed = symbolicFwdSeed(nfwd, arg);
    // Wrapper input being constructed
    vector<MX> w_in = arg;
    w_in.insert(w_in.end(), res.begin(), res.end());
    // Arguments for calling sw being constructed
    vector<MX> v;
    v.insert(v.end(), arg.begin(), arg.end());
    v.insert(v.end(), res.begin(), res.end());
    for (int d=0; d<nfwd; ++d) {
      // ignore seed for ind
      seed[d][0] = MX::sym(seed[d][0].getName(), Sparsity::scalar(false));
      // Add to wrapper input
      w_in.insert(w_in.end(), seed[d].begin(), seed[d].end());
      // Add to sw argument vector
      v.insert(v.end(), seed[d].begin()+1, seed[d].end());
    }

    // Create wrapper
    casadi_assert(v.size()==sw.n_in());
    return MXFunction(name, w_in, sw(v), opts);
  }

  Function SwitchInternal
  ::getDerReverse(const std::string& name, int nadj, Dict& opts) {
    // Derivative of each case
    vector<Function> der(f_.size());
    for (int k=0; k<f_.size(); ++k) {
      if (!f_[k].isNull()) der[k] = f_[k].derReverse(nadj);
    }

    // Default case
    Function der_def;
    if (!f_def_.isNull()) der_def = f_def_.derReverse(nadj);

    // New Switch for derivatives
    stringstream ss;
    ss << "adj" << nadj << "_" << name_;
    Switch sw(ss.str(), der, der_def);

    // Construct wrapper inputs and arguments for calling sw
    vector<MX> arg = mx_in();
    vector<MX> res = symbolicOutput();
    vector<vector<MX> > seed = symbolicAdjSeed(nadj, res);
    vector<MX> w_in = arg;
    w_in.insert(w_in.end(), res.begin(), res.end());
    // Arguments for calling sw being constructed
    vector<MX> v;
    v.push_back(arg.at(0)); // index
    for (int d=0; d<nadj; ++d) {
      // Add to wrapper input
      w_in.insert(w_in.end(), seed[d].begin(), seed[d].end());
      // Add to sw argument vector
      v.insert(v.end(), arg.begin()+1, arg.end());
      v.insert(v.end(), res.begin(), res.end());
      v.insert(v.end(), seed[d].begin(), seed[d].end());
    }

    // Construct wrapper outputs
    casadi_assert(v.size()==sw.n_in());
    v = sw(v);
    vector<MX> w_out;
    MX ind_sens = MX(1, 1); // no dependency on index
    vector<MX>::const_iterator v_it = v.begin(), v_it_next;
    for (int d=0; d<nadj; ++d) {
      w_out.push_back(ind_sens);
      v_it_next = v_it + (n_in()-1);
      w_out.insert(w_out.end(), v_it, v_it_next);
      v_it = v_it_next;
    }

    // Create wrapper
    return MXFunction(name, w_in, w_out, opts);
  }

  inline string name(const Function& f) {
    if (f.isNull()) {
      return "NULL";
    } else {
      return f.name();
    }
  }

  void SwitchInternal::print(ostream &stream) const {
    if (f_.size()==1) {
      // Print as if-then-else
      stream << "Switch(" << name(f_def_) << ", " << name(f_[0]) << ")";
    } else {
      // Print generic
      stream << "Switch([";
      for (int k=0; k<f_.size(); ++k) {
        if (k!=0) stream << ", ";
        stream << name(f_[k]);
      }
      stream << "], " << name(f_def_) << ")";
    }
  }

  void SwitchInternal::generateDeclarations(CodeGenerator& g) const {
    for (int k=0; k<=f_.size(); ++k) {
      const Function& fk = k<f_.size() ? f_[k] : f_def_;
      fk->addDependency(g);
    }
  }

  void SwitchInternal::generateBody(CodeGenerator& g) const {
    // Codegen as if-else
    bool if_else = f_.size()==1;

    // Codegen condition
    g.body << "  " << (if_else ? "if" : "switch")  << " (to_int(arg[0][0])) {" << endl;

    // Loop over cases/functions
    for (int k=0; k<=f_.size(); ++k) {

      // For if,  reverse order
      int k1 = if_else ? 1-k : k;

      if (!if_else) {
        // Codegen cases
        if (k1<f_.size()) {
          g.body << "  case " << k1 << ":" << endl;
        } else {
          g.body << "  default:" << endl;
        }
      } else if (k1==0) {
        // Else
        g.body << "  } else {" << endl;
      }

      // Get the function:
      const Function& fk = k1<f_.size() ? f_[k1] : f_def_;
      if (fk.isNull()) {
        g.body << "    return 1;" << endl;
      } else {
        // Call function
        if (g.simplifiedCall(fk)) {
          casadi_error("Not implemented.");
        } else {
          g.body << "    if (" << g.call(fk, "arg+1", "res", "iw", "w") << ") return 1;" << endl;
          if (!if_else)
            g.body << "    break;" << endl;
        }
      }
    }

    // End switch/else
    g.body << "  }" << endl;
  }

} // namespace casadi
