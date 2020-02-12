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
#include "serializing_stream.hpp"

#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.thread.h>
#else // CASADI_WITH_THREAD_MINGW
#include <thread>
#endif // CASADI_WITH_THREAD_MINGW
#endif // CASADI_WITH_THREAD

using namespace std;

namespace casadi {

  Function Map::create(const std::string& parallelization, const Function& f, casadi_int n) {
    // Create instance of the right class
    string suffix = str(n) + "_" + f.name();
    if (parallelization == "serial") {
      return Function::create(new Map("map" + suffix, f, n), Dict());
    } else if (parallelization== "openmp") {
      return Function::create(new OmpMap("ompmap" + suffix, f, n), Dict());
    } else if (parallelization== "thread") {
      return Function::create(new ThreadMap("threadmap" + suffix, f, n), Dict());
    } else {
      casadi_error("Unknown parallelization: " + parallelization);
    }
  }

  Map::Map(const std::string& name, const Function& f, casadi_int n)
    : FunctionInternal(name), f_(f), n_(n) {
  }

  void Map::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);
    s.pack("Map::f", f_);
    s.pack("Map::n", n_);
  }

  void Map::serialize_type(SerializingStream &s) const {
    FunctionInternal::serialize_type(s);
    s.pack("Map::class_name", class_name());
  }

  Map::Map(DeserializingStream& s) : FunctionInternal(s) {
    s.unpack("Map::f", f_);
    s.unpack("Map::n", n_);
  }

  ProtoFunction* Map::deserialize(DeserializingStream& s) {
    std::string class_name;
    s.unpack("Map::class_name", class_name);
    if (class_name=="Map") {
      return new Map(s);
    } else if (class_name=="OmpMap") {
      return new OmpMap(s);
    } else if (class_name=="ThreadMap") {
      return new ThreadMap(s);
    } else {
      casadi_error("class name '" + class_name + "' unknown.");
    }
  }

  Map::~Map() {
    clear_mem();
  }

  void Map::init(const Dict& opts) {
    is_diff_in_ = f_.is_diff_in();
    is_diff_out_ = f_.is_diff_out();
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    align_w_ = f_.align_w();

    // Allocate sufficient memory for serial evaluation
    alloc_arg(f_.sz_arg());
    alloc_res(f_.sz_res());
    alloc_w(f_.sz_w());
    alloc_iw(f_.sz_iw());
  }

  template<typename T>
  int Map::eval_gen(const T** arg, T** res, casadi_int* iw, T* w, int mem) const {
    const T** arg1 = arg+n_in_;
    copy_n(arg, n_in_, arg1);
    T** res1 = res+n_out_;
    copy_n(res, n_out_, res1);
    for (casadi_int i=0; i<n_; ++i) {
      if (f_(arg1, res1, iw, w, mem)) return 1;
      for (casadi_int j=0; j<n_in_; ++j) {
        if (arg1[j]) arg1[j] += f_.nnz_in(j);
      }
      for (casadi_int j=0; j<n_out_; ++j) {
        if (res1[j]) res1[j] += f_.nnz_out(j);
      }
    }
    return 0;
  }

  int Map::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w, void* mem) const {
    return eval_gen(arg, res, iw, w);
  }

  int Map::sp_forward(const bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    return eval_gen(arg, res, iw, w);
  }

  int Map::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    bvec_t** arg1 = arg+n_in_;
    copy_n(arg, n_in_, arg1);
    bvec_t** res1 = res+n_out_;
    copy_n(res, n_out_, res1);
    for (casadi_int i=0; i<n_; ++i) {
      if (f_.rev(arg1, res1, iw, w)) return 1;
      for (casadi_int j=0; j<n_in_; ++j) {
        if (arg1[j]) arg1[j] += f_.nnz_in(j);
      }
      for (casadi_int j=0; j<n_out_; ++j) {
        if (res1[j]) res1[j] += f_.nnz_out(j);
      }
    }
    return 0;
  }

  void Map::codegen_declarations(CodeGenerator& g, const Instance& inst) const {
    uout() << "inst" << inst.arg_null << std::endl;
    g.add_dependency(f_, inst);
  }

  void Map::codegen_body(CodeGenerator& g,
      const Instance& inst) const {
    g.local("i", "casadi_int");
    g.local("arg1[" + str(f_.sz_arg()) + "]", "const casadi_real*");
    g.local("res1[" + str(f_.sz_res()) + "]", "casadi_real*");

    // Input buffer
    g << "for (i=0; i<" << n_in_ << "; ++i) arg1[i]=arg[i];\n";
    // Output buffer
    if (str(f_).find("SXFunction")!= std::string::npos) {
      g << "for (i=0; i<" << n_out_ << "; ++i) res1[i]=res[i];\n"
        << "#pragma omp simd\n"
        << "for (i=0; i<" << n_ << "; ++i) {\n";
    } else {
      g << "for (i=0; i<" << n_out_ << "; ++i) res1[i]=res[i];\n"
        << "for (i=0; i<" << n_ << "; ++i) {\n";
    }
    // Evaluate
    if (str(f_).find("SXFunction")!= std::string::npos) {
      g << g(f_, "arg1", "res1", "iw", "w", inst) << ";\n";
    } else {
      g << "if (" << g(f_, "arg1", "res1", "iw", "w", inst) << ") return 1;\n";
    }

    // Update input buffers
    for (casadi_int j=0; j<n_in_; ++j) {
      if (f_.nnz_in(j)) {
        if (inst.arg_null.empty()) {
          g << "if (arg1[" << j << "]) arg1[" << j << "]+=" << f_.nnz_in(j) << ";\n";
        } else {
          if (!inst.arg_null[j]) g << "arg1[" << j << "]+=" << f_.nnz_in(j) << ";\n";
        }
      }
    }
    // Update output buffers
    for (casadi_int j=0; j<n_out_; ++j) {
      if (f_.nnz_out(j)) {
        if (inst.res_null.empty()) {
          g << "if (res1[" << j << "]) res1[" << j << "]+=" << f_.nnz_out(j) << ";\n";
        } else {
          if (!inst.res_null[j]) g << "res1[" << j << "]+=" << f_.nnz_out(j) << ";\n";
        }
      }
    }
    g << "}\n";
  }

  Function Map
  ::get_forward(casadi_int nfwd, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Generate map of derivative
    Function df = f_.forward(nfwd);
    Function dm = df.map(n_, parallelization());

    // Input expressions
    vector<MX> arg = dm.mx_in();

    // Need to reorder sensitivity inputs
    vector<MX> res = arg;
    vector<MX>::iterator it=res.begin()+n_in_+n_out_;
    vector<casadi_int> ind;
    for (casadi_int i=0; i<n_in_; ++i, ++it) {
      casadi_int sz = f_.size2_in(i);
      ind.clear();
      for (casadi_int k=0; k<n_; ++k) {
        for (casadi_int d=0; d<nfwd; ++d) {
          for (casadi_int j=0; j<sz; ++j) {
            ind.push_back((d*n_ + k)*sz + j);
          }
        }
      }
      *it = (*it)(Slice(), ind); // NOLINT
    }

    // Get output expressions
    res = dm(res);

    // Reorder sensitivity outputs
    it = res.begin();
    for (casadi_int i=0; i<n_out_; ++i, ++it) {
      casadi_int sz = f_.size2_out(i);
      ind.clear();
      for (casadi_int d=0; d<nfwd; ++d) {
        for (casadi_int k=0; k<n_; ++k) {
          for (casadi_int j=0; j<sz; ++j) {
            ind.push_back((k*nfwd + d)*sz + j);
          }
        }
      }
      *it = (*it)(Slice(), ind); // NOLINT
    }

    // Construct return function
    return Function(name, arg, res, inames, onames, opts);
  }

  Function Map
  ::get_reverse(casadi_int nadj, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Generate map of derivative
    Function df = f_.reverse(nadj);
    Function dm = df.map(n_, parallelization());

    // Input expressions
    vector<MX> arg = dm.mx_in();

    // Need to reorder sensitivity inputs
    vector<MX> res = arg;
    vector<MX>::iterator it=res.begin()+n_in_+n_out_;
    vector<casadi_int> ind;
    for (casadi_int i=0; i<n_out_; ++i, ++it) {
      casadi_int sz = f_.size2_out(i);
      ind.clear();
      for (casadi_int k=0; k<n_; ++k) {
        for (casadi_int d=0; d<nadj; ++d) {
          for (casadi_int j=0; j<sz; ++j) {
            ind.push_back((d*n_ + k)*sz + j);
          }
        }
      }
      *it = (*it)(Slice(), ind); // NOLINT
    }

    // Get output expressions
    res = dm(res);

    // Reorder sensitivity outputs
    it = res.begin();
    for (casadi_int i=0; i<n_in_; ++i, ++it) {
      casadi_int sz = f_.size2_in(i);
      ind.clear();
      for (casadi_int d=0; d<nadj; ++d) {
        for (casadi_int k=0; k<n_; ++k) {
          for (casadi_int j=0; j<sz; ++j) {
            ind.push_back((k*nadj + d)*sz + j);
          }
        }
      }
      *it = (*it)(Slice(), ind); // NOLINT
    }

    // Construct return function
    return Function(name, arg, res, inames, onames, opts);
  }

  int Map::eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    // This checkout/release dance is an optimization.
    // Could also use the thread-safe variant f_(arg1, res1, iw, w)
    // in Map::eval_gen
    scoped_checkout<Function> m(f_);
    return eval_gen(arg, res, iw, w, m);
  }

  OmpMap::~OmpMap() {
    clear_mem();
  }

  int OmpMap::eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
#ifndef WITH_OPENMP
    return Map::eval(arg, res, iw, w, mem);
#else // WITH_OPENMP
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);

    // Error flag
    casadi_int flag = 0;

    // Checkout memory objects
    std::vector< scoped_checkout<Function> > ind; ind.reserve(n_);
    for (casadi_int i=0; i<n_; ++i) ind.emplace_back(f_);

    // Evaluate in parallel
#pragma omp parallel for reduction(||:flag)
    for (casadi_int i=0; i<n_; ++i) {
      // Input buffers
      const double** arg1 = arg + n_in_ + i*sz_arg;
      for (casadi_int j=0; j<n_in_; ++j) {
        arg1[j] = arg[j] ? arg[j] + i*f_.nnz_in(j) : 0;
      }

      // Output buffers
      double** res1 = res + n_out_ + i*sz_res;
      for (casadi_int j=0; j<n_out_; ++j) {
        res1[j] = res[j] ? res[j] + i*f_.nnz_out(j) : 0;
      }

      // Evaluation
      try {
        flag = f_(arg1, res1, iw + i*sz_iw, w + i*sz_w, ind[i]) || flag;
      } catch (std::exception& e) {
        flag = 1;
        casadi_warning("Exception raised: " + std::string(e.what()));
      } catch (...) {
        flag = 1;
        casadi_warning("Uncaught exception.");
      }

    }

    // Return error flag
    return flag;
#endif  // WITH_OPENMP
  }

  void OmpMap::codegen_body(CodeGenerator& g,
      const Instance& inst) const {
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);
    g << "casadi_int i;\n"
      << "const double** arg1;\n"
      << "double** res1;\n"
      << "casadi_int flag = 0;\n"
      << "#pragma omp parallel for private(i,arg1,res1) reduction(||:flag)\n"
      << "for (i=0; i<" << n_ << "; ++i) {\n"
      << "arg1 = arg + " << n_in_ << "+i*" << sz_arg << ";\n";
    for (casadi_int j=0; j<n_in_; ++j) {
      g << "arg1[" << j << "] = arg[" << j << "] ? "
        << g.arg(j) << "+i*" << f_.nnz_in(j) << ": 0;\n";
    }
    g << "res1 = res + " <<  n_out_ << "+i*" <<  sz_res << ";\n";
    for (casadi_int j=0; j<n_out_; ++j) {
      g << "res1[" << j << "] = res[" << j << "] ?"
        << g.res(j) << "+i*" << f_.nnz_out(j) << ": 0;\n";
    }
    g << "flag = "
      << g(f_, "arg1", "res1", "iw+i*" + str(sz_iw), "w+i*" + str(sz_w), inst) << " || flag;\n"
      << "}\n"
      << "if (flag) return 1;\n";
  }

  void OmpMap::init(const Dict& opts) {
#ifndef CASADI_WITH_THREAD
    casadi_warning("CasADi was not compiled with WITH_THREAD=ON. "
                   "Falling back to serial evaluation.");
#endif // CASADI_WITH_THREAD
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


  ThreadMap::~ThreadMap() {
    clear_mem();
  }

  void ThreadsWork(const Function& f, casadi_int i,
      const double** arg, double** res,
      casadi_int* iw, double* w,
      casadi_int ind, int& ret) {

    // Function dimensions
    casadi_int n_in = f.n_in();
    casadi_int n_out = f.n_out();

    // Function work sizes
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f.sz_work(sz_arg, sz_res, sz_iw, sz_w);

    // Input buffers
    const double** arg1 = arg + n_in + i*sz_arg;
    for (casadi_int j=0; j<n_in; ++j) {
      arg1[j] = arg[j] ? arg[j] + i*f.nnz_in(j) : nullptr;
    }

    // Output buffers
    double** res1 = res + n_out + i*sz_res;
    for (casadi_int j=0; j<n_out; ++j) {
      res1[j] = res[j] ? res[j] + i*f.nnz_out(j) : nullptr;
    }

    try {
      ret = f(arg1, res1, iw + i*sz_iw, w + i*sz_w, ind);
    } catch (std::exception& e) {
      ret = 1;
      casadi_warning("Exception raised: " + std::string(e.what()));
    } catch (...) {
      ret = 1;
      casadi_warning("Uncaught exception.");
    }
  }

  int ThreadMap::eval(const double** arg, double** res, casadi_int* iw, double* w,
      void* mem) const {

#ifndef CASADI_WITH_THREAD
    return Map::eval(arg, res, iw, w, mem);
#else // CASADI_WITH_THREAD
    // Checkout memory objects
    std::vector< scoped_checkout<Function> > ind; ind.reserve(n_);
    for (casadi_int i=0; i<n_; ++i) ind.emplace_back(f_);

    // Allocate space for return values
    std::vector<int> ret_values(n_);

    // Spawn threads
    std::vector<std::thread> threads;
    for (casadi_int i=0; i<n_; ++i) {
      // Why the lambda function?
      // Because it was the first iteration to pass tests on MingGW
      // using mingw-std-threads.
      threads.emplace_back(
        [i](const Function& f, const double** arg, double** res,
            casadi_int* iw, double* w, casadi_int ind, int& ret) {
              ThreadsWork(f, i, arg, res, iw, w, ind, ret);
            },
        std::ref(f_), arg, res, iw, w, casadi_int(ind[i]), std::ref(ret_values[i]));
    }

    // Join threads
    for (auto && th : threads) th.join();

    // Anticipate success
    int ret = 0;

    // Compute aggregate return value
    for (int e : ret_values) ret = ret || e;

    return ret;
#endif // CASADI_WITH_THREAD
  }

  void ThreadMap::codegen_body(CodeGenerator& g,
      const Instance& inst) const {
    Map::codegen_body(g, inst);
  }

  void ThreadMap::init(const Dict& opts) {
#ifndef CASADI_WITH_THREAD
    casadi_warning("CasADi was not compiled with WITH_THREAD=ON. "
                   "Falling back to serial evaluation.");
#endif // CASADI_WITH_THREAD
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
