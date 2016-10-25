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


#include "switch.hpp"

using namespace std;

namespace casadi {

  Switch::Switch(const std::string& name,
                 const std::vector<Function>& f, const Function& f_def)
    : FunctionInternal(name), f_(f), f_def_(f_def) {

    // Consitency check
    casadi_assert(!f_.empty());
  }

  Switch::~Switch() {
  }

  size_t Switch::get_n_in() {
    for (auto&& i : f_) if (!i.is_null()) return 1+i.n_in();
    casadi_assert(!f_def_.is_null());
    return 1+f_def_.n_in();
  }

  size_t Switch::get_n_out() {
    for (auto&& i : f_) if (!i.is_null()) return i.n_out();
    casadi_assert(!f_def_.is_null());
    return f_def_.n_out();
  }

  Sparsity Switch::get_sparsity_in(int i) {
    if (i==0) {
      return Sparsity::scalar();
    } else {
      Sparsity ret;
      for (auto&& fk : f_) {
        if (!fk.is_null()) {
          const Sparsity& s = fk.sparsity_in(i-1);
          ret = ret.is_null() ? s : ret.unite(s);
        }
      }
      casadi_assert(!f_def_.is_null());
      const Sparsity& s = f_def_.sparsity_in(i-1);
      ret = ret.is_null() ? s : ret.unite(s);
      return ret;
    }
  }

  Sparsity Switch::get_sparsity_out(int i) {
    Sparsity ret;
    for (auto&& fk : f_) {
      if (!fk.is_null()) {
        const Sparsity& s = fk.sparsity_out(i);
        ret = ret.is_null() ? s : ret.unite(s);
      }
    }
    casadi_assert(!f_def_.is_null());
    const Sparsity& s = f_def_.sparsity_out(i);
    ret = ret.is_null() ? s : ret.unite(s);
    return ret;
  }

  void Switch::init(const Dict& opts) {
    // Initialize the functions, get input and output sparsities
    // Input and output sparsities
    std::vector<Sparsity> sp_in, sp_out;
    int num_in = -1, num_out=-1;
    for (int k=0; k<=f_.size(); ++k) {
      Function& fk = k<f_.size() ? f_[k] : f_def_;
      if (fk.is_null()) continue;
      if (num_in<0) {
        // Number of inputs and outputs
        num_in=fk.n_in();
        num_out=fk.n_out();
        // Output sparsity
        sp_out.resize(num_out);
        for (int i=0; i<num_out; ++i) sp_out[i] = fk.sparsity_out(i);
        // Input sparsity
        sp_in.resize(num_in);
        for (int i=0; i<num_in; ++i) sp_in[i] = fk.sparsity_in(i);
      } else {
        // Assert matching number of inputs and outputs
        casadi_assert(num_in==fk.n_in());
        casadi_assert(num_out==fk.n_out());
        // Intersect with output sparsity
        for (int i=0; i<num_out; ++i) {
          sp_out[i] = sp_out[i].intersect(fk.sparsity_out(i));
        }
        // Intersect with input sparsity
        for (int i=0; i<num_in; ++i) {
          sp_in[i] = sp_in[i].intersect(fk.sparsity_in(i));
        }
      }
    }

    // Illegal to pass only "null" functions
    casadi_assert_message(num_in>=0, "All functions are null");

    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Get required work
    for (int k=0; k<=f_.size(); ++k) {
      const Function& fk = k<f_.size() ? f_[k] : f_def_;
      if (fk.is_null()) continue;

      // Get local work vector sizes
      alloc(fk);
      size_t sz_w = fk.sz_w();

      // Add size for input buffers
      for (int i=1; i<n_in(); ++i) {
        const Sparsity& s = fk.sparsity_in(i-1);
        if (s!=sparsity_in(i)) sz_w += s.nnz();
      }

      // Add size for output buffers
      for (int i=0; i<n_out(); ++i) {
        const Sparsity& s = fk.sparsity_out(i);
        if (s!=sparsity_out(i)) sz_w += s.nnz();
      }

      // Make sure enough work for this
      alloc_w(sz_w);
    }
  }

  void Switch::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    // Get conditional
    int k = static_cast<int>(*arg[0]);

    // Get the function to be evaluated
    const Function& fk = k<f_.size() ? f_[k] : f_def_;

    // Input buffers
    for (int i=1; i<n_in(); ++i) {
      const Sparsity& s = fk.sparsity_in(i-1);
      casadi_assert_message(s==sparsity_in(i), "Not implemented");
    }

    // Output buffers
    for (int i=0; i<n_out(); ++i) {
      const Sparsity& s = fk.sparsity_out(i);
      casadi_assert_message(s==sparsity_out(i), "Not implemented");
    }

    // Evaluate the corresponding function
    fk(arg+1, res, iw, w, 0);
  }

  Function Switch
  ::get_forward(const std::string& name, int nfwd, Dict& opts) {
    // Derivative of each case
    vector<Function> der(f_.size());
    for (int k=0; k<f_.size(); ++k) {
      if (!f_[k].is_null()) der[k] = f_[k].forward_new(nfwd);
    }

    // Default case
    Function der_def;
    if (!f_def_.is_null()) der_def = f_def_.forward_new(nfwd);

    // New Switch for derivatives
    stringstream ss;
    ss << "fwd" << nfwd << "_" << name_;
    Function sw = Function::conditional(ss.str(), der, der_def);

    // Construct wrapper inputs and arguments for calling sw
    vector<MX> arg = sw.mx_in();

    // ignore seed for ind
    vector<MX>::iterator ind_seed = arg.begin() + n_in() + n_out();
    *ind_seed = MX(ind_seed->size());

    // Get result expression
    vector<MX> res = sw(arg);

    // Remove index seed from return function input
    arg.erase(ind_seed);

    // Create wrapper
    return Function(name, arg, res, opts);
  }

  Function Switch
  ::get_reverse(const std::string& name, int nadj, Dict& opts) {
    // Derivative of each case
    vector<Function> der(f_.size());
    for (int k=0; k<f_.size(); ++k) {
      if (!f_[k].is_null()) der[k] = f_[k].reverse_new(nadj);
    }

    // Default case
    Function der_def;
    if (!f_def_.is_null()) der_def = f_def_.reverse_new(nadj);

    // New Switch for derivatives
    stringstream ss;
    ss << "adj" << nadj << "_" << name_;
    Function sw = Function::conditional(ss.str(), der, der_def);

    // Construct wrapper inputs and arguments for calling sw
    vector<MX> arg = sw.mx_in();

    // ignore seed for ind
    vector<MX>::iterator ind_seed = arg.begin() + n_in() + n_out();
    *ind_seed = MX(ind_seed->size());

    // Get result expression
    vector<MX> res = sw(arg);

    // Remove index seed from return function input
    arg.erase(ind_seed);

    // Create wrapper
    return Function(name, arg, res, opts);
  }

  void Switch::print(ostream &stream) const {
    if (f_.size()==1) {
      // Print as if-then-else
      stream << "Switch(" << f_def_.name() << ", " << f_[0].name() << ")";
    } else {
      // Print generic
      stream << "Switch([";
      for (int k=0; k<f_.size(); ++k) {
        if (k!=0) stream << ", ";
        stream << f_[k].name();
      }
      stream << "], " << f_def_.name() << ")";
    }
  }

  void Switch::generateDeclarations(CodeGenerator& g) const {
    for (int k=0; k<=f_.size(); ++k) {
      const Function& fk = k<f_.size() ? f_[k] : f_def_;
      fk->addDependency(g);
    }
  }

  void Switch::generateBody(CodeGenerator& g) const {
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
      if (fk.is_null()) {
        g.body << "    return 1;" << endl;
      } else {
        // Call function
        if (g.simplifiedCall(fk)) {
          casadi_error("Not implemented.");
        } else {
          g.body << "    if (" << g(fk, "arg+1", "res", "iw", "w") << ") return 1;" << endl;
          if (!if_else)
            g.body << "    break;" << endl;
        }
      }
    }

    // End switch/else
    g.body << "  }" << endl;
  }

} // namespace casadi
