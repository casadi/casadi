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
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Buffer for mismatching sparsities
    size_t sz_buf=0;

    // Keep track of sparsity projections
    project_in_ = project_out_ = false;

    // Get required work
    for (int k=0; k<=f_.size(); ++k) {
      const Function& fk = k<f_.size() ? f_[k] : f_def_;
      if (fk.is_null()) continue;

      // Memory for evaluation
      alloc(fk);

      // Required work vectors
      size_t sz_buf_k=0;

      // Add size for input buffers
      for (int i=1; i<n_in(); ++i) {
        const Sparsity& s = fk.sparsity_in(i-1);
        if (s!=sparsity_in(i)) {
          project_in_ = true;
          alloc_w(s.size1()); // for casadi_project
          sz_buf_k += s.nnz();
        }
      }

      // Add size for output buffers
      for (int i=0; i<n_out(); ++i) {
        const Sparsity& s = fk.sparsity_out(i);
        if (s!=sparsity_out(i)) {
          project_out_ = true;
          alloc_w(s.size1()); // for casadi_project
          sz_buf_k += s.nnz();
        }
      }

      // Only need the largest of these work vectors
      sz_buf = max(sz_buf, sz_buf_k);
    }

    // Memory for the work vectors
    alloc_w(sz_buf, true);
  }

  void Switch::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    // Shorthands
    int n_in=this->n_in()-1, n_out=this->n_out();

    // Get the function to be evaluated
    int k = arg[0] ? static_cast<int>(*arg[0]) : 0;
    const Function& fk = k>=0 && k<f_.size() ? f_[k] : f_def_;

    // Project arguments with different sparsity
    const double** arg1;
    if (project_in_) {
      // Project one or more argument
      arg1 = arg + 1 + n_in;
      for (int i=0; i<n_in; ++i) {
        const Sparsity& f_sp = fk.sparsity_in(i);
        const Sparsity& sp = sparsity_in(i+1);
        arg1[i] = arg[i+1];
        if (arg1[i] && f_sp!=sp) {
          casadi_project(arg1[i], sp, w, f_sp, w + f_sp.nnz());
          arg1[i] = w; w += f_sp.nnz();
        }
      }
    } else {
      // No inputs projected
      arg1 = arg + 1;
    }

    // Temporary memory for results with different sparsity
    double** res1;
    if (project_out_) {
      // Project one or more results
      res1 = res + n_out;
      for (int i=0; i<n_out; ++i) {
        const Sparsity& f_sp = fk.sparsity_out(i);
        const Sparsity& sp = sparsity_out(i);
        res1[i] = res[i];
        if (res1[i] && f_sp!=sp) {
          res1[i] = w;
          w += f_sp.nnz();
        }
      }
    } else {
      // No outputs projected
      res1 = res;
    }

    // Evaluate the corresponding function
    fk(arg1, res1, iw, w, 0);

    // Project results with different sparsity
    if (project_out_) {
      for (int i=0; i<n_out; ++i) {
        const Sparsity& f_sp = fk.sparsity_out(i);
        const Sparsity& sp = sparsity_out(i);
        if (res[i] && f_sp!=sp) {
          casadi_project(res1[i], f_sp, res[i], sp, w);
        }
      }
    }
  }

  Function Switch
  ::get_forward(int nfwd, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Derivative of each case
    vector<Function> der(f_.size());
    for (int k=0; k<f_.size(); ++k) {
      if (!f_[k].is_null()) der[k] = f_[k].forward(nfwd);
    }

    // Default case
    Function der_def;
    if (!f_def_.is_null()) der_def = f_def_.forward(nfwd);

    // New Switch for derivatives
    Function sw = Function::conditional("switch_" + name, der, der_def);

    // Get expressions for the derivative switch
    vector<MX> arg = sw.mx_in();
    vector<MX> res = sw(arg);

    // Ignore seed for ind
    arg.insert(arg.begin() + n_in() + n_out(), MX(1, nfwd));

    // Create wrapper
    return Function(name, arg, res, inames, onames, opts);
  }

  Function Switch
  ::get_reverse(const std::string& name, int nadj,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Derivative of each case
    vector<Function> der(f_.size());
    for (int k=0; k<f_.size(); ++k) {
      if (!f_[k].is_null()) der[k] = f_[k].reverse(nadj);
    }

    // Default case
    Function der_def;
    if (!f_def_.is_null()) der_def = f_def_.reverse(nadj);

    // New Switch for derivatives
    Function sw = Function::conditional("switch_" + name, der, der_def);

    // Get expressions for the derivative switch
    vector<MX> arg = sw.mx_in();
    vector<MX> res = sw(arg);

    // No derivatives with respect to index
    res.insert(res.begin(), MX(1, nadj));

    // Create wrapper
    return Function(name, arg, res, inames, onames, opts);
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

  void Switch::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const {

    // Shorthands
    int n_in=this->n_in()-1, n_out=this->n_out();

    // Input and output buffers
    const SXElem** arg1 = arg + 1 + n_in;
    SXElem** res1 = res + n_out;

    // Extra memory needed for chaining if_else calls
    std::vector<SXElem> w_extra(nnz_out());
    std::vector<SXElem*> res_tempv(n_out);
    SXElem** res_temp = get_ptr(res_tempv);

    for (int k=0; k<f_.size()+1; ++k) {

      // Local work vector
      SXElem* wl = w;

      // Local work vector
      SXElem* wll = get_ptr(w_extra);

      if (k==0) {
        // For the default case, redirect the temporary results to res
        copy_n(res, n_out, res_temp);
      } else {
        // For the other cases, store the temporary results
        for (int i=0; i<n_out; ++i) {
          res_temp[i] = wll;
          wll += nnz_out(i);
        }
      }

      copy_n(arg+1, n_in, arg1);
      copy_n(res_temp, n_out, res1);

      const Function& fk = k==0 ? f_def_ : f_[k-1];

      // Project arguments with different sparsity
      for (int i=0; i<n_in; ++i) {
        if (arg1[i]) {
          const Sparsity& f_sp = fk.sparsity_in(i);
          const Sparsity& sp = sparsity_in(i+1);
          if (f_sp!=sp) {
            SXElem *t = wl; wl += f_sp.nnz(); // t is non-const
            casadi_project(arg1[i], sp, t, f_sp, wl);
            arg1[i] = t;
          }
        }
      }

      // Temporary memory for results with different sparsity
      for (int i=0; i<n_out; ++i) {
        if (res1[i]) {
          const Sparsity& f_sp = fk.sparsity_out(i);
          const Sparsity& sp = sparsity_out(i);
          if (f_sp!=sp) { res1[i] = wl; wl += f_sp.nnz();}
        }
      }

      // Evaluate the corresponding function
      fk(arg1, res1, iw, wl, 0);

      // Project results with different sparsity
      for (int i=0; i<n_out; ++i) {
        if (res1[i]) {
          const Sparsity& f_sp = fk.sparsity_out(i);
          const Sparsity& sp = sparsity_out(i);
          if (f_sp!=sp) casadi_project(res1[i], f_sp, res_temp[i], sp, wl);
        }
      }

      if (k>0) { // output the temporary results via an if_else
        SXElem cond = k-1==arg[0][0];
        for (int i=0; i<n_out; ++i) {
          if (res[i]) {
            for (int j=0; j<nnz_out(i); ++j) {
              res[i][j] = if_else(cond, res_temp[i][j], res[i][j], true);
            }
          }
        }
      }

    }

  }

  void Switch::generateBody(CodeGenerator& g) const {
    // Shorthands
    int n_in=this->n_in()-1, n_out=this->n_out();

    // Project arguments with different sparsity
    if (project_in_) {
      // Project one or more argument
      g.local("i", "int");
      g << "const real_t** arg1 = arg + " << (1 + n_in) << ";\n"
        << "for(i=0; i<" << n_in << "; ++i) arg1[i]=arg[i+1];\n";
    }

    // Temporary memory for results with different sparsity
    if (project_out_) {
      // Project one or more results
      g.local("i", "int");
      g << "real_t** res1 = res + " << n_out << ";\n"
        << "for (i=0; i<" << n_out << "; ++i) res1[i]=res[i];\n";
    }

    // Codegen condition
    bool if_else = f_.size()==1;
    g << (if_else ? "if" : "switch")  << " (arg[0] ? to_int(*arg[0]) : 0) {\n";

    // Loop over cases/functions
    for (int k=0; k<=f_.size(); ++k) {

      // For if,  reverse order
      int k1 = if_else ? 1-k : k;

      if (!if_else) {
        // Codegen cases
        if (k1<f_.size()) {
          g << "case " << k1 << ":\n";
        } else {
          g << "default:\n";
        }
      } else if (k1==0) {
        // Else
        g << "} else {\n";
      }

      // Get the function:
      const Function& fk = k1<f_.size() ? f_[k1] : f_def_;
      if (fk.is_null()) {
        g << "return 1;\n";
      } else if (g.simplifiedCall(fk)) {
        casadi_error("Not implemented.");
      } else {
        // Project arguments with different sparsity
        for (int i=0; i<n_in; ++i) {
          const Sparsity& f_sp = fk.sparsity_in(i);
          const Sparsity& sp = sparsity_in(i+1);
          if (f_sp!=sp) {
            if (f_sp.nnz()==0) {
              g << "arg1[" << i << "]=0;\n";
            } else {
              g.local("t", "real_t", "*");
              g << "t=w, w+=" << f_sp.nnz() << ";\n"
                << g.project("arg1[" + to_string(i) + "]", sp, "t", f_sp, "w") << "\n"
                << "arg1[" << i << "]=t;\n";
            }
          }
        }

        // Temporary memory for results with different sparsity
        for (int i=0; i<n_out; ++i) {
          const Sparsity& f_sp = fk.sparsity_out(i);
          const Sparsity& sp = sparsity_out(i);
          if (f_sp!=sp) {
            if (f_sp.nnz()==0) {
              g << "res1[" << i << "]=0;\n";
            } else {
              g << "res1[" << i << "]=w, w+=" << f_sp.nnz() << ";\n";
            }
          }
        }

        // Function call
        g << "if (" << g(fk, project_in_ ? "arg1" : "arg+1",
                         project_out_ ? "res1" : "res",
                         "iw", "w") << ") return 1;\n";

        // Project results with different sparsity
        for (int i=0; i<n_out; ++i) {
          const Sparsity& f_sp = fk.sparsity_out(i);
          const Sparsity& sp = sparsity_out(i);
          if (f_sp!=sp) {
            g << g.project("res1[" + to_string(i) + "]", f_sp,
                           "res[" + to_string(i) + "]", sp, "w") << "\n";
          }
        }

        // Break (if switch)
        if (!if_else) g << "break;\n";
      }
    }

    // End switch/else
    g << "}\n";
  }

} // namespace casadi
