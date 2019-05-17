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


#include "mapsum.hpp"
#include "serializing_stream.hpp"

using namespace std;

namespace casadi {

  Function MapSum::create(const std::string& name, const std::string& parallelization,
                          const Function& f, casadi_int n,
                          const std::vector<bool>& reduce_in,
                          const std::vector<bool>& reduce_out,
                          const Dict& opts) {
    if (reduce_out.empty()) return create(name, parallelization, f, n,
                                     reduce_in, std::vector<bool>(f.n_out(), false));
    casadi_assert(reduce_in.size()==f.n_in(), "Dimension mismatch");
    casadi_assert(reduce_out.size()==f.n_out(), "Dimension mismatch");

    if (parallelization == "serial") {
      string suffix = str(reduce_in)+str(reduce_out);
      Function ret;
      if (!f->incache(name, ret, suffix)) {
        // Create new serial map
        ret = Function::create(new MapSum(name, f, n, reduce_in, reduce_out), opts);
        casadi_assert_dev(ret.name()==name);
        // Save in cache
        f->tocache(ret, suffix);
      }
      return ret.wrap_as_needed(opts);
    } else {
      casadi_error("Unknown parallelization: " + parallelization);
    }
  }

  MapSum::MapSum(const std::string& name, const Function& f, casadi_int n,
                 const std::vector<bool>& reduce_in,
                 const std::vector<bool>& reduce_out)
    : FunctionInternal(name), f_(f), n_(n), reduce_in_(reduce_in), reduce_out_(reduce_out) {
    casadi_assert_dev(reduce_in.size()==f.n_in());
    casadi_assert_dev(reduce_out.size()==f.n_out());
  }

  void MapSum::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);
    s.pack("MapSum::f", f_);
    s.pack("MapSum::n", n_);
    s.pack("MapSum::reduce_in", reduce_in_);
    s.pack("MapSum::reduce_out", reduce_out_);
  }

  void MapSum::serialize_type(SerializingStream &s) const {
    FunctionInternal::serialize_type(s);
    s.pack("MapSum::class_name", class_name());
  }

  MapSum::MapSum(DeserializingStream& s) : FunctionInternal(s) {
    s.unpack("MapSum::f", f_);
    s.unpack("MapSum::n", n_);
    s.unpack("MapSum::reduce_in", reduce_in_);
    s.unpack("MapSum::reduce_out", reduce_out_);
  }

  ProtoFunction* MapSum::deserialize(DeserializingStream& s) {
    std::string class_name;
    s.unpack("MapSum::class_name", class_name);
    if (class_name=="MapSum") {
      return new MapSum(s);
    } else {
      casadi_error("class name '" + class_name + "' unknown.");
    }
  }

  MapSum::~MapSum() {
    clear_mem();
  }

  void MapSum::init(const Dict& opts) {
    is_diff_in_ = f_.is_diff_in();
    is_diff_out_ = f_.is_diff_out();

    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Allocate sufficient memory for serial evaluation
    alloc_arg(f_.sz_arg());
    alloc_res(f_.sz_res());
    alloc_w(f_.sz_w(), true);
    alloc_iw(f_.sz_iw());

    // Allocate scratch space for dummping result of reduced outputs
    for (casadi_int j=0;j<n_out_;++j) {
      if (reduce_out_[j]) alloc_w(f_.nnz_out(j), true);
    }
  }

  template<typename T1>
  void casadi_add(casadi_int n, const T1* x, T1* y) {
    casadi_int i;
    if (!x || !y) return;
    for (i=0; i<n; ++i) *y++ += *x++;
  }

  template<>
  void casadi_add(casadi_int n, const bvec_t* x, bvec_t* y) {
    casadi_int i;
    if (!x || !y) return;
    for (i=0; i<n; ++i) *y++ |= *x++;
  }

  template<typename T>
  int MapSum::eval_gen(const T** arg, T** res, casadi_int* iw, T* w, int mem) const {
    const T** arg1 = arg+n_in_;
    copy_n(arg, n_in_, arg1);
    T** res1 = res+n_out_;

    T* w_scratch = w + f_.sz_w();
    for (casadi_int j=0;j<n_out_;++j) {
      if (res[j] && reduce_out_[j]) {
        casadi_clear(res[j], f_.nnz_out(j)); // clear sums
        res1[j] = w_scratch; // Make the function dump result in scratch space
        w_scratch += f_.nnz_out(j);
      } else {
        res1[j] = res[j];
      }
    }
    for (casadi_int i=0; i<n_; ++i) {
      if (f_(arg1, res1, iw, w, mem)) return 1;
      for (casadi_int j=0; j<n_in_; ++j) {
        if (arg1[j] && !reduce_in_[j]) arg1[j] += f_.nnz_in(j);
      }
      for (casadi_int j=0; j<n_out_; ++j) {
        if (res1[j]) {
          if (reduce_out_[j]) {
            casadi_add(f_.nnz_out(j), res1[j], res[j]); // Perform sum
          } else {
            res1[j] += f_.nnz_out(j);
          }
        }
      }
    }
    return 0;
  }

  int MapSum::eval_sx(const SXElem** arg, SXElem** res,
      casadi_int* iw, SXElem* w, void* mem) const {
    return eval_gen(arg, res, iw, w);
  }

  int MapSum::sp_forward(const bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    return eval_gen(arg, res, iw, w);
  }

  int MapSum::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    // Note: f_.rev(arg,res,iw,w)
    // has a side effect of clearing res
    // Reduced outputs should not be cleared;
    // they must influence each iteration

    // Store reduced res in scratch space
    bvec_t* w_scratch = w + f_.sz_w();
    for (casadi_int j=0;j<n_out_;++j) {
      if (res[j] && reduce_out_[j]) {
        casadi_copy(res[j], f_.nnz_out(j), w_scratch);
        w_scratch += f_.nnz_out(j);
      }
    }
    bvec_t** arg1 = arg+n_in_;
    copy_n(arg, n_in_, arg1);
    bvec_t** res1 = res+n_out_;
    copy_n(res, n_out_, res1);
    for (casadi_int i=0; i<n_; ++i) {
      // Restore res1[j] from scratch space
      w_scratch = w + f_.sz_w();
      for (casadi_int j=0;j<n_out_;++j) {
        if (res[j] && reduce_out_[j]) {
          casadi_copy(w_scratch, f_.nnz_out(j), res1[j]);
          w_scratch += f_.nnz_out(j);
        }
      }
      if (f_.rev(arg1, res1, iw, w)) return 1;
      for (casadi_int j=0; j<n_in_; ++j) {
        if (arg1[j] && !reduce_in_[j]) arg1[j] += f_.nnz_in(j);
      }
      for (casadi_int j=0; j<n_out_; ++j) {
        if (res1[j] && !reduce_out_[j]) res1[j] += f_.nnz_out(j);
      }
    }
    return 0;
  }

  void MapSum::codegen_declarations(CodeGenerator& g) const {
    g.add_dependency(f_);
  }

  void MapSum::codegen_body(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_CLEAR);
    g.local("i", "casadi_int");
    g.local("arg1", "const casadi_real*", "*");
    g.local("res1", "casadi_real*", "*");
    g.local("w_scratch", "casadi_real*", "*");
    // Input buffer
    g << "arg1 = arg+" << n_in_ << ";\n"
      << "for (i=0; i<" << n_in_ << "; ++i) arg1[i]=arg[i];\n";
    // Output buffer
    g << "res1 = res+" << n_out_ << ";\n";
    g << "w_scratch = w+" << f_.sz_w() << ";\n";
    for (casadi_int j=0;j<n_out_;++j) {
      if (reduce_out_[j]) {
        g << "if (res[" << j << "]) {\n";
        g << "casadi_clear(res[" << j << "], " << f_.nnz_out(j) << ");\n";
        g << "res1[" << j << "] = w_scratch;\n";
        g << "w_scratch+=" << f_.nnz_out(j) << ";\n";
        g << "} else {\n";
        g << "res1[" << j << "] = res[" << j << "];\n";
        g << "}\n";
      } else {
        g << "res1[" << j << "] = res[" << j << "];\n";
      }
    }

    g << "for (i=0; i<" << n_ << "; ++i) {\n";
    // Evaluate
    g << "if (" << g(f_, "arg1", "res1", "iw", "w") << ") return 1;\n";
    // Update input buffers
    for (casadi_int j=0; j<n_in_; ++j) {
      if (!reduce_in_[j] && f_.nnz_in(j)) {
        g << "if (arg1[" << j << "]) arg1[" << j << "]+=" << f_.nnz_in(j) << ";\n";
      }
    }
    // Update output buffers
    for (casadi_int j=0; j<n_out_; ++j) {
      if (reduce_out_[j]) {
        g << "if (res1[" << j << "]) ";
        g << g.axpy(f_.nnz_out(j), "1.0", "res1[" + str(j) + "]", "res[" + str(j) + "]") << "\n";
      } else {
        if (f_.nnz_out(j)) {
          g << "if (res1[" << j << "]) ";
          g << "res1[" << j << "]+=" << f_.nnz_out(j) << ";\n";
        }
      }
    }
    g << "}\n";
  }

  Function MapSum
  ::get_forward(casadi_int nfwd, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Generate map of derivative
    Function df = f_.forward(nfwd);

    for (casadi_int i=0;i<n_out_;++i) {
      if (reduce_out_[i]) casadi_assert(df.nnz_in(n_in_+i)==0, "Case not implemented");
    }

    std::vector<bool> reduce_in = join(reduce_in_, reduce_out_, reduce_in_);
    Function dm = MapSum::create("mapsum" + str(n_) + "_" + df.name(), parallelization(),
      df, n_, reduce_in, reduce_out_);

    // Input expressions
    vector<MX> arg = dm.mx_in();

    // Need to reorder sensitivity inputs
    vector<MX> res = arg;
    vector<MX>::iterator it=res.begin()+n_in_+n_out_;
    vector<casadi_int> ind;
    for (casadi_int i=0; i<n_in_; ++i, ++it) {
      if (reduce_in_[i]) continue;
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
      if (reduce_out_[i]) continue;
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
    Dict custom_opts = opts;
    custom_opts["always_inline"] = true;
    return Function(name, arg, res, inames, onames, custom_opts);
  }

  Function MapSum
  ::get_reverse(casadi_int nadj, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Generate map of derivative
    Function df = f_.reverse(nadj);

    for (casadi_int i=0;i<n_out_;++i) {
      if (reduce_out_[i]) casadi_assert(df.nnz_in(n_in_+i)==0, "Case not implemented");
    }

    std::vector<bool> reduce_in = join(reduce_in_, reduce_out_, reduce_out_);
    Function dm = MapSum::create("mapsum" + str(n_) + "_" + df.name(), parallelization(),
      df, n_, reduce_in, reduce_in_);

    // Input expressions
    vector<MX> arg = dm.mx_in();

    // Need to reorder sensitivity inputs
    vector<MX> res = arg;
    vector<MX>::iterator it=res.begin()+n_in_+n_out_;
    vector<casadi_int> ind;
    for (casadi_int i=0; i<n_out_; ++i, ++it) {
      if (reduce_out_[i]) continue;
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
      if (reduce_in_[i]) continue;
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
    Dict custom_opts = opts;
    custom_opts["always_inline"] = true;
    return Function(name, arg, res, inames, onames, custom_opts);
  }

  int MapSum::eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    // This checkout/release dance is an optimization.
    // Could also use the thread-safe variant f_(arg1, res1, iw, w)
    // in Map::eval_gen
    scoped_checkout<Function> m(f_);
    return eval_gen(arg, res, iw, w, m);
  }

} // namespace casadi
