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


#include "map.hpp"

using namespace std;

namespace casadi {

  Function Map::create(const std::string& parallelization, const Function& f, int n) {
    // Create instance of the right class
    string name = f.name() + "_" + to_string(n);
    if (parallelization == "serial") {
      return Function::create(new Map(name, f, n), Dict());
    } else if (parallelization== "openmp") {
      return Function::create(new MapOmp(name, f, n), Dict());
    } else {
      casadi_error("Unknown parallelization: " + parallelization);
    }
  }

  Map::Map(const std::string& name, const Function& f, int n)
    : FunctionInternal(name), f_(f), n_(n) {
  }

  Map::~Map() {
  }

  void Map::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Allocate sufficient memory for serial evaluation
    alloc_arg(f_.sz_arg());
    alloc_res(f_.sz_res());
    alloc_w(f_.sz_w());
    alloc_iw(f_.sz_iw());
  }

  template<typename T>
  int Map::evalGen(const T** arg, T** res, int* iw, T* w) const {
    int n_in = this->n_in(), n_out = this->n_out();
    const T** arg1 = arg+n_in;
    copy_n(arg, n_in, arg1);
    T** res1 = res+n_out;
    copy_n(res, n_out, res1);
    for (int i=0; i<n_; ++i) {
      if (f_(arg1, res1, iw, w, 0)) return 1;
      for (int j=0; j<n_in; ++j) {
        if (arg1[j]) arg1[j] += f_.nnz_in(j);
      }
      for (int j=0; j<n_out; ++j) {
        if (res1[j]) res1[j] += f_.nnz_out(j);
      }
    }
    return 0;
  }

  int Map::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const {
    return evalGen(arg, res, iw, w);
  }

  int Map::sp_forward(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const {
    return evalGen(arg, res, iw, w);
  }

  void Map::sp_reverse(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const {
    int n_in = this->n_in(), n_out = this->n_out();
    bvec_t** arg1 = arg+n_in;
    copy_n(arg, n_in, arg1);
    bvec_t** res1 = res+n_out;
    copy_n(res, n_out, res1);
    for (int i=0; i<n_; ++i) {
      f_->sp_reverse(arg1, res1, iw, w, 0);
      for (int j=0; j<n_in; ++j) {
        if (arg1[j]) arg1[j] += f_.nnz_in(j);
      }
      for (int j=0; j<n_out; ++j) {
        if (res1[j]) res1[j] += f_.nnz_out(j);
      }
    }
  }

  void Map::codegen_declarations(CodeGenerator& g) const {
    f_->add_dependency(g);
  }

  void Map::codegen_body(CodeGenerator& g) const {
    int n_in = this->n_in(), n_out = this->n_out();
    g << "int i;\n";
    // Input buffer
    g << "const casadi_real** arg1 = arg+" << n_in << ";\n"
      << "for (i=0; i<" << n_in << "; ++i) arg1[i]=arg[i];\n";
    // Output buffer
    g << "casadi_real** res1 = res+" << n_out << ";\n"
      << "for (i=0; i<" << n_out << "; ++i) res1[i]=res[i];\n"
      << "for (i=0; i<" << n_ << "; ++i) {\n";
    // Evaluate
    g << "if (" << g(f_, "arg1", "res1", "iw", "w") << ") return 1;\n";
    // Update input buffers
    for (int j=0; j<n_in; ++j) {
      g << "if (arg1[" << j << "]) arg1[" << j << "]+=" << f_.nnz_in(j) << ";\n";
    }
    // Update output buffers
    for (int j=0; j<n_out; ++j) {
      g << "if (res1[" << j << "]) res1[" << j << "]+=" << f_.nnz_out(j) << ";\n";
    }
    g << "}\n";
  }

  Function Map
  ::get_forward(int nfwd, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Shorthands
    int n_in = this->n_in(), n_out = this->n_out();

    // Generate map of derivative
    Function df = f_.forward(nfwd);
    Function dm = df.map(n_, parallelization());

    // Input expressions
    vector<MX> arg = dm.mx_in();

    // Need to reorder sensitivity inputs
    vector<MX> res = arg;
    vector<MX>::iterator it=res.begin()+n_in+n_out;
    vector<int> ind;
    for (int i=0; i<n_in; ++i, ++it) {
      int sz = f_.size2_in(i);
      ind.clear();
      for (int k=0; k<n_; ++k) {
        for (int d=0; d<nfwd; ++d) {
          for (int j=0; j<sz; ++j) {
            ind.push_back((d*n_ + k)*sz + j);
          }
        }
      }
      *it = (*it)(Slice(), ind);
    }

    // Get output expressions
    res = dm(res);

    // Reorder sensitivity outputs
    it = res.begin();
    for (int i=0; i<n_out; ++i, ++it) {
      int sz = f_.size2_out(i);
      ind.clear();
      for (int d=0; d<nfwd; ++d) {
        for (int k=0; k<n_; ++k) {
          for (int j=0; j<sz; ++j) {
            ind.push_back((k*nfwd + d)*sz + j);
          }
        }
      }
      *it = (*it)(Slice(), ind);
    }

    // Construct return function
    return Function(name, arg, res, inames, onames, opts);
  }

  Function Map
  ::get_reverse(int nadj, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Shorthands
    int n_in = this->n_in(), n_out = this->n_out();

    // Generate map of derivative
    Function df = f_.reverse(nadj);
    Function dm = df.map(n_, parallelization());

    // Input expressions
    vector<MX> arg = dm.mx_in();

    // Need to reorder sensitivity inputs
    vector<MX> res = arg;
    vector<MX>::iterator it=res.begin()+n_in+n_out;
    vector<int> ind;
    for (int i=0; i<n_out; ++i, ++it) {
      int sz = f_.size2_out(i);
      ind.clear();
      for (int k=0; k<n_; ++k) {
        for (int d=0; d<nadj; ++d) {
          for (int j=0; j<sz; ++j) {
            ind.push_back((d*n_ + k)*sz + j);
          }
        }
      }
      *it = (*it)(Slice(), ind);
    }

    // Get output expressions
    res = dm(res);

    // Reorder sensitivity outputs
    it = res.begin();
    for (int i=0; i<n_in; ++i, ++it) {
      int sz = f_.size2_in(i);
      ind.clear();
      for (int d=0; d<nadj; ++d) {
        for (int k=0; k<n_; ++k) {
          for (int j=0; j<sz; ++j) {
            ind.push_back((k*nadj + d)*sz + j);
          }
        }
      }
      *it = (*it)(Slice(), ind);
    }

    // Construct return function
    return Function(name, arg, res, inames, onames, opts);
  }

  int Map::eval(const double** arg, double** res, int* iw, double* w, void* mem) const {
    return evalGen(arg, res, iw, w);
  }

  MapOmp::~MapOmp() {
  }

  int MapOmp::eval(const double** arg, double** res, int* iw, double* w, void* mem) const {
#ifndef WITH_OPENMP
    return Map::eval(arg, res, iw, w, mem);
#else // WITH_OPENMP
    int n_in = this->n_in(), n_out = this->n_out();
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);

    // Checkout memory objects
    int* ind = iw; iw += n_;
    for (int i=0; i<n_; ++i) ind[i] = f_.checkout();

    // Evaluate in parallel
#pragma omp parallel for
    for (int i=0; i<n_; ++i) {
      // Input buffers
      const double** arg1 = arg + n_in + i*sz_arg;
      for (int j=0; j<n_in; ++j) {
        arg1[j] = arg[j] ? arg[j] + i*f_.nnz_in(j) : 0;
      }

      // Output buffers
      double** res1 = res + n_out + i*sz_res;
      for (int j=0; j<n_out; ++j) {
        res1[j] = res[j] ? res[j] + i*f_.nnz_out(j) : 0;
      }

      // Evaluation
      f_(arg1, res1, iw + i*sz_iw, w + i*sz_w, ind[i]);
    }
    // Release memory objects
    for (int i=0; i<n_; ++i) f_.release(ind[i]);
#endif  // WITH_OPENMP
  }

  void MapOmp::codegen_body(CodeGenerator& g) const {
    int n_in = this->n_in(), n_out = this->n_out();
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);
    g << "int i;\n"
      << "const double** arg1;\n"
      << "double** res1;\n"
      << "#pragma omp parallel for private(i,arg1,res1)\n"
      << "for (i=0; i<" << n_ << "; ++i) {\n"
      << "arg1 = arg + " << n_in << "+i*" << sz_arg << ";\n";
    for (int j=0; j<n_in; ++j) {
      g << "arg1[" << j << "] = arg[" << j << "] ? "
        << "arg[" << j << "]+i*" << f_.nnz_in(j) << ": 0;\n";
    }
    g << "res1 = res + " <<  n_out << "+i*" <<  sz_res << ";\n";
    for (int j=0; j<n_out; ++j) {
      g << "res1[" << j << "] = res[" << j << "] ?"
        << "res[" << j << "]+i*" << f_.nnz_out(j) << ": 0;\n";
    }
    g << g(f_, "arg1", "res1",
           "iw+i*" + to_string(sz_iw), "w+i*" + to_string(sz_w)) << ";\n"
      << "}\n";
  }

  void MapOmp::init(const Dict& opts) {
    // Call the initialization method of the base class
    Map::init(opts);

    // Allocate memory for holding memory object references
    alloc_iw(n_, true);

    // Allocate sufficient memory for parallel evaluation
    alloc_arg(f_.sz_arg() * n_);
    alloc_res(f_.sz_res() * n_);
    alloc_w(f_.sz_w() * n_);
    alloc_iw(f_.sz_iw() * n_);
  }

} // namespace casadi
