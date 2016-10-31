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

  Function Map::create(const std::string& name,
                          const std::string& parallelization, Function& f, int n,
                          const Dict& opts) {
    // Create instance of the right class
    Function ret;
    if (parallelization == "serial") {
      ret.assignNode(new Map(name, f, n));
    } else if (parallelization== "openmp") {
      ret.assignNode(new MapOmp(name, f, n));
    } else {
      casadi_error("Unknown parallelization: " + parallelization);
    }
    // Finalize creation
    ret->construct(opts);
    return ret;
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
  void Map::evalGen(const T** arg, T** res, int* iw, T* w) const {
    int n_in = this->n_in(), n_out = this->n_out();
    const T** arg1 = arg+n_in;
    T** res1 = res+n_out;
    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in; ++j) {
        arg1[j] = arg[j] ? arg[j]+i*f_.nnz_in(j): 0;
      }
      for (int j=0; j<n_out; ++j) {
        res1[j]= res[j] ? res[j]+i*f_.nnz_out(j): 0;
      }
      f_(arg1, res1, iw, w, 0);
    }
  }

  void Map::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen(arg, res, iw, w);
  }

  void Map::sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    evalGen(arg, res, iw, w);
  }

  void Map::sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    int n_in = this->n_in(), n_out = this->n_out();
    bvec_t** arg1 = arg+n_in;
    bvec_t** res1 = res+n_out;
    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in; ++j) {
        arg1[j] = arg[j] ? arg[j]+i*f_.nnz_in(j): 0;
      }
      for (int j=0; j<n_out; ++j) {
        res1[j]= res[j] ? res[j]+i*f_.nnz_out(j): 0;
      }
      f_->sp_rev(arg1, res1, iw, w, 0);
    }
  }

  void Map::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void Map::generateBody(CodeGenerator& g) const {
    int n_in = this->n_in(), n_out = this->n_out();

    g.body << "  const real_t** arg1 = arg+" << n_in << ";"<< endl;
    g.body << "  real_t** res1 = res+" << n_out << ";" << endl;

    g.body << "  int i;" << endl;
    g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;
    for (int j=0; j<n_in; ++j) {
      g.body << "    arg1[" << j << "] = arg[" << j << "]? " <<
        "arg[" << j << "]+i*" << f_.nnz_in(j) << " : 0;" << endl;
    }
    for (int j=0; j<n_out; ++j) {
      g.body << "    res1[" << j << "] = res[" << j << "]? " <<
        "res[" << j << "]+i*" << f_.nnz_out(j) << " : 0;" << endl;
    }
    g.body << "    if (" << g(f_, "arg1", "res1", "iw", "w") << ") return 1;" << endl;
    g.body << "  }" << std::endl;
  }

  Function Map
  ::get_forward_old(const std::string& name, int nfwd, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.forward(nfwd);

    // Construct and return
    return df.map(name, parallelization(), n_, opts);
  }

  Function Map
  ::get_forward(const std::string& name, int nfwd, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.forward_new(nfwd);

    // Construct and return
    Function dm = df.map(name + "_map", parallelization(), n_, derived_options());

    // Input expressions
    vector<MX> arg = dm.mx_in();

    // Need to reorder sensitivity inputs
    vector<MX> parg = arg;
    for (int i=0; i<n_in(); ++i) {
      vector<MX> v = horzsplit(arg[n_in()+n_out()+i], f_.size2_in(i));
      casadi_assert(v.size()==n_*nfwd);
      vector<MX> w(v.size());
      for (int k=0; k<n_; ++k) {
        for (int d=0; d<nfwd; ++d) {
          w[k*nfwd + d] = v[d*n_ + k];
        }
      }
      parg[n_in()+n_out()+i] = horzcat(w);
    }

    // Get output expressions
    vector<MX> res = dm(parg);

    // Reorder sensitivity outputs
    for (int i=0; i<n_out(); ++i) {
      vector<MX> v = horzsplit(res[i], f_.size2_out(i));
      casadi_assert(v.size()==n_*nfwd);
      vector<MX> w(v.size());
      for (int k=0; k<n_; ++k) {
        for (int d=0; d<nfwd; ++d) {
          w[d*n_ + k] = v[k*nfwd + d];
        }
      }
      res[i] = horzcat(w);
    }

    // Construct return function
    return Function(name, arg, res, opts);
  }

  Function Map
  ::get_reverse(const std::string& name, int nadj, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.reverse_new(nadj);

    // Construct and return
    Function dm = df.map(name + "_map", parallelization(), n_, derived_options());

    // Input expressions
    vector<MX> arg = dm.mx_in();

    // Need to reorder sensitivity inputs
    vector<MX> parg = arg;
    for (int i=0; i<n_out(); ++i) {
      vector<MX> v = horzsplit(arg[n_in()+n_out()+i], f_.size2_out(i));
      casadi_assert(v.size()==n_*nadj);
      vector<MX> w(v.size());
      for (int k=0; k<n_; ++k) {
        for (int d=0; d<nadj; ++d) {
          w[k*nadj + d] = v[d*n_ + k];
        }
      }
      parg[n_in()+n_out()+i] = horzcat(w);
    }

    // Get output expressions
    vector<MX> res = dm(parg);

    // Reorder sensitivity outputs
    for (int i=0; i<n_in(); ++i) {
      vector<MX> v = horzsplit(res[i], f_.size2_in(i));
      casadi_assert(v.size()==n_*nadj);
      vector<MX> w(v.size());
      for (int k=0; k<n_; ++k) {
        for (int d=0; d<nadj; ++d) {
          w[d*n_ + k] = v[k*nadj + d];
        }
      }
      res[i] = horzcat(w);
    }

    // Construct return function
    return Function(name, arg, res, opts);
  }

  Function Map
  ::get_reverse_old(const std::string& name, int nadj, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.reverse(nadj);

    // Construct and return
    return df.map(name, parallelization(), n_, opts);
  }

  void Map::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    evalGen(arg, res, iw, w);
  }

  MapOmp::~MapOmp() {
  }

  void MapOmp::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
#ifdef WITH_OPENMP
    return Map::eval(mem, arg, res, iw, w);
#else // WITH_OPENMP
    int n_in = this->n_in(), n_out = this->n_out();
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);
#pragma omp parallel for
    for (int i=0; i<n_; ++i) {
      const double** arg_i = arg + n_in + sz_arg*i;
      for (int j=0; j<n_in; ++j) {
        arg_i[j] = arg[j]+i*f_.nnz_in(j);
      }
      double** res_i = res + n_out + sz_res*i;
      for (int j=0; j<n_out; ++j) {
        res_i[j] = res[j]? res[j]+i*f_.nnz_out(j) : 0;
      }
      int* iw_i = iw + i*sz_iw;
      double* w_i = w + i*sz_w;
      f_->eval(0, arg_i, res_i, iw_i, w_i);
    }
#endif  // WITH_OPENMP
  }

  void MapOmp::generateBody(CodeGenerator& g) const {
    int n_in = this->n_in(), n_out = this->n_out();
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);

    g.body << "  int i;" << endl;
    g.body << "#pragma omp parallel for" << endl;
    g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;
    g.body << "    const double** arg_i = arg + " << n_in << "+" << sz_arg << "*i;" << endl;
    for (int j=0; j<n_in; ++j) {
      g.body << "    arg_i[" << j << "] = arg[" << j << "]+i*" << f_.nnz_in(j) << ";" << endl;
    }
    g.body << "    double** res_i = res + " <<  n_out << "+" <<  sz_res << "*i;" << endl;
    for (int j=0; j<n_out; ++j) {
      g.body << "    res_i[" << j << "] = res[" << j << "] ?" <<
                "res[" << j << "]+i*" << f_.nnz_out(j) << ": 0;" << endl;
    }
    g.body << "    int* iw_i = iw + i*" << sz_iw << ";" << endl;
    g.body << "    double* w_i = w + i*" << sz_w << ";" << endl;
    g.body << "    " << g(f_, "arg_i", "res_i", "iw_i", "w_i") << ";" << endl;
    g.body << "  }" << std::endl;
  }

  void MapOmp::init(const Dict& opts) {
    // Call the initialization method of the base class
    Map::init(opts);

    // Allocate sufficient memory for parallel evaluation
    alloc_arg(f_.sz_arg() * n_);
    alloc_res(f_.sz_res() * n_);
    alloc_w(f_.sz_w() * n_);
    alloc_iw(f_.sz_iw() * n_);
  }

} // namespace casadi
