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
#include "global_options.hpp"

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

    Function ret;
    if (parallelization == "serial") {
      string suffix = str(reduce_in)+str(reduce_out);
      if (!f->incache(name, ret, suffix)) {
        // Create new serial map
        ret = Function::create(new MapSum(name, f, n, reduce_in, reduce_out), opts);
        casadi_assert_dev(ret.name()==name);
        // Save in cache
        f->tocache(ret, suffix);
      }
    } else {
      casadi_error("Unknown parallelization: " + parallelization);
    }

    if (!vectorize_f(f, n)) return ret.wrap_as_needed(opts);

    const MapSum& m = *static_cast<const MapSum*>(ret.get());

    // Input expressions
    vector<MX> arg = ret.mx_in(); // change to layout of reinterpret_layout

    vector<MX> res = m.permute_out(ret(m.permute_in(arg)));

    // Construct return function
    Dict custom_opts = opts;
    custom_opts["always_inline"] = true;
    uout() << "opts" << opts << std::endl;
    return Function(ret.name(), arg, res, custom_opts);

  }

  MapSum::MapSum(const std::string& name, const Function& f, casadi_int n,
                 const std::vector<bool>& reduce_in,
                 const std::vector<bool>& reduce_out)
    : FunctionInternal(name), f_(f), n_(n), reduce_in_(reduce_in), reduce_out_(reduce_out) {
    casadi_assert_dev(reduce_in.size()==f.n_in());
    casadi_assert_dev(reduce_out.size()==f.n_out());
    uout() << "init It's map sum!" << std::endl;

    f_orig_ = f;

    if (vectorize_f(f, n_)) {
      std::vector<casadi_int> stride_in(f_.n_in());
      std::vector<casadi_int> stride_out(f_.n_out());
      for (casadi_int j=0; j<f_.n_in(); ++j) {
        stride_in[j] = vectorize_f(f, n_) && !reduce_in[j] ? n_padded(n) : 1;
      }
      for (casadi_int j=0; j<f_.n_out(); ++j) {
        stride_out[j] = vectorize_f(f, n_) ? (reduce_out_[j] ? GlobalOptions::vector_width_real : n_padded(n)) : 1;
      }

      Dict opts;
      opts["stride_in"] = stride_in;
      opts["stride_out"] = stride_out;
      f_ = f_->with_options(opts);
    }
  }

  Layout MapSum::get_layout_in(casadi_int i) {
    if (!vectorize_f()) return Layout();
    return permute_in(i).target();
  }

  Layout MapSum::get_layout_out(casadi_int i) {
    if (!vectorize_f()) return Layout();
    return permute_out(i).source();
  }

  void MapSum::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);
    s.pack("MapSum::f", f_);
    s.pack("MapSum::f_orig", f_orig_);
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
    s.unpack("MapSum::f_orig", f_orig_);
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

 std::vector<std::string> MapSum::get_function() const {
    return {"f"};
  }

  const Function& MapSum::get_function(const std::string &name) const {
    casadi_assert(has_function(name),
      "No function \"" + name + "\" in " + name_ + ". " +
      "Available functions: " + join(get_function()) + ".");
    return f_;
  }

  bool MapSum::has_function(const std::string& fname) const {
    return fname=="f";
  }

  void MapSum::init(const Dict& opts) {
    is_diff_in_ = f_.is_diff_in();
    is_diff_out_ = f_.is_diff_out();

    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Allocate sufficient memory for serial evaluation
    alloc_arg(f_.sz_arg());
    alloc_res(f_.sz_res());
    alloc_w(f_.sz_w()*GlobalOptions::vector_width_real, true);
    alloc_iw(f_.sz_iw());

    // Allocate scratch space for dummping result of reduced outputs
    for (casadi_int j=0;j<n_out_;++j) {
      if (reduce_out_[j]) alloc_w(f_.nnz_out(j), true);
    }

    dump_in_ = true;
    dump_out_ = true;

  }

  casadi_int MapSum::n_padded() const {
    if (vectorize_f()) {
      return n_padded(n_);
    }
    return n_;
  }

  casadi_int MapSum::n_padded(casadi_int n) {
    return n;
    casadi_int rem = n % GlobalOptions::vector_width_real;
    if (rem==0) return n;
    return n + GlobalOptions::vector_width_real - rem;
  }

  vector<MX> MapSum::permute_in(const vector<MX> & arg, bool invert) const {
    vector<MX> ret = arg;
    for (casadi_int i=0;i<f_.n_in();++i) {
      ret[i] = permute_layout(ret[i], permute_in(i, invert));
    }
    return ret;
  }

  vector<MX> MapSum::permute_out(const vector<MX> & res, bool invert) const {
    vector<MX> ret = res;
    for (casadi_int i=0;i<f_.n_out();++i) {
      ret[i] = permute_layout(ret[i], permute_out(i, invert));
    }
    return ret;
  }

  Relayout MapSum::permute_in(casadi_int i, bool invert) const {
    if (reduce_in_[i]) return Relayout();
    Layout source({f_.nnz_in(i), n_});
    Layout target({n_, f_.nnz_in(i)}, {n_padded(), f_.nnz_in(i)});
    Relayout ret = Relayout(source, {1, 0}, target);
    if (invert) return ret.invert();
    return ret;
  }

  Relayout MapSum::permute_out(casadi_int i, bool invert) const {
    if (reduce_out_[i]) return Relayout();
    Layout source({n_, f_.nnz_out(i)}, {n_padded(), f_.nnz_out(i)});
    Layout target({f_.nnz_out(i), n_});
    Relayout ret = Relayout(source, {1, 0}, target);
    if (invert) return ret.invert();
    return ret;
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
  int MapSum::local_eval_gen(const T** arg, T** res, casadi_int* iw, T* w, int mem) const {


    uout() << "local_eval_gen:" << name_ << ":" << dump_in_ << std::endl;
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
        if (arg1[j] && !reduce_in_[j]) arg1[j] += vectorize_f() ? 1 : f_.nnz_in(j);
      }
      for (casadi_int j=0; j<n_out_; ++j) {
        if (res1[j]) {
          if (reduce_out_[j]) {
            casadi_add(f_.nnz_out(j), res1[j], res[j]); // Perform sum
          } else {
            res1[j] += vectorize_f() ? 1 : f_.nnz_out(j);
          }
        }
      }
    }
    return 0;
  }

  int MapSum::eval_sx(const SXElem** arg, SXElem** res,
      casadi_int* iw, SXElem* w, void* mem) const {
    return local_eval_gen(arg, res, iw, w);
  }

  int MapSum::sp_forward(const bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    return local_eval_gen(arg, res, iw, w);
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

  void MapSum::codegen_declarations(CodeGenerator& g, const Instance& inst) const {
    Instance local = inst;
    local.stride_in.resize(f_.n_in());
    for (casadi_int j=0; j<n_in_; ++j) {
      local.stride_in[j] = (vectorize_f() && !reduce_in_[j]) ? n_padded() : 1;
    }
    local.stride_out.resize(f_.n_out());
    for (casadi_int j=0; j<n_out_; ++j) {
      local.stride_out[j] = vectorize_f() ? (reduce_out_[j] ? GlobalOptions::vector_width_real : n_padded()) : 1;
    }
    g.add_dependency(f_, local);
  }

  void MapSum::codegen_body(CodeGenerator& g, const Instance& inst) const {
    Instance local = inst;
    local.stride_in.resize(f_.n_in());
    for (casadi_int j=0; j<n_in_; ++j) {
      local.stride_in[j] = (vectorize_f() && !reduce_in_[j]) ? n_padded() : 1;
    }
    local.stride_out.resize(f_.n_out());
    for (casadi_int j=0; j<n_out_; ++j) {
      local.stride_out[j] = vectorize_f() ? (reduce_out_[j] ? GlobalOptions::vector_width_real : n_padded()) : 1;
    }

    bool any_reduce_out = any(reduce_out_);

    // Outer loop is needed to make reduction on outputs
    bool outer_loop = vectorize_f() && any_reduce_out;

    g.add_dependency(f_, local);
    g.add_auxiliary(CodeGenerator::AUX_CLEAR);
    g.local("i", "casadi_int");
    g.local("arg1[" + str(f_.sz_arg()) + "]", "const casadi_real*");

    g.local("res1[" + str(f_.sz_res()) + "]", "casadi_real*");
    if (any_reduce_out) {
      g.local("w_scratch", "casadi_real*", "*");
    }



    // loop bug: if inner body contains branches, and outer_loop is active, gcc-8 gets confused.
    // about data dependency (gcc-9 not)
    bool loop_bug = outer_loop;
    // Input buffer
    g  << "for (i=0; i<" << n_in_ << "; ++i) arg1[i]=arg[i];\n";
    if (loop_bug) {
      g  << "for (i=0; i<" << n_in_ << "; ++i) arg2[i]=arg[i];\n";
    }
    // Output buffer
    if (any_reduce_out) {
      g << "w_scratch = w+" << f_.sz_w() << ";\n";
    }
    for (casadi_int j=0;j<n_out_;++j) {
      if (reduce_out_[j]) {
        g << "if (res[" << j << "]) {\n";
        g << "casadi_clear(res[" << j << "], " << f_.nnz_out(j) << ");\n";
        //if (vectorize_f()) {
        //  g << "res1[" << j << "] = res[" << j << "];\n";
        //} else {
          g << "res"<< (loop_bug? "2" : "1") <<"[" << j << "] = w_scratch;\n";
          g << "w_scratch+=" << f_.nnz_out(j)*GlobalOptions::vector_width_real << ";\n";
        //}
        g << "} else {\n";
        g << "res"<< (loop_bug? "2" : "1") <<"[" << j << "] = res[" << j << "];\n";
        g << "}\n";
        if (loop_bug) g << "res1[" << j << "] = res2[" << j << "];\n";
      } else {
        g << "res"<< (loop_bug? "2" : "1") <<"[" << j << "] = res[" << j << "];\n";
      }
      
    }

    if (loop_bug) {
      g.local("arg2[" + str(n_in_) + "]", "const casadi_real*", "*");
      g.local("res2[" + str(n_out_) + "]", "casadi_real*", "*");
    }

    if (vectorize_f()) {
      g.local("j", "casadi_int");
      if (outer_loop) {
        g << "for (i=0; i<" << n_padded() / GlobalOptions::vector_width_real << "; ++i) {\n";
        if (loop_bug) {
          for (casadi_int j=0; j<n_in_; ++j) {
            if (!reduce_in_[j] && f_.nnz_in(j)) {
              g << "arg1[" << j << "]=arg2[" << j << "];\n";
            }
          }
          for (casadi_int j=0; j<n_out_; ++j) {
            if (f_.nnz_out(j)) {
              g << "res1[" << j << "]=res2[" << j << "];\n";
            }
          }
        }
        g << "#pragma omp simd\n";
        g << "for (j=0; j<" << GlobalOptions::vector_width_real << "; ++j) {\n";
      } else {
        g << "#pragma omp simd\n";
        g << "for (j=0; j<" << n_padded() << "; ++j) {\n";
      }
    } else {
      g << "for (i=0; i<" << n_ << "; ++i) {\n";
    }

    // Evaluate
    if (str(f_).find("SXFunction")!= std::string::npos) {
      g << g(f_, "arg1", "res1", "iw", "w", local) << ";\n";
    } else {
      g << "if (" << g(f_, "arg1", "res1", "iw", "w", local) << ") return 1;\n";
    }

    uout() << "debug" << name_ << std::endl;
    uout() << "n_in_" << n_in_ << std::endl;
    uout() << "reduce_in_" << reduce_in_ << std::endl;
    uout() << "inst.arg_null" << inst.arg_null << std::endl;
    uout() << "inst.arg_null" << inst.arg_null << std::endl;

    // Update input buffers
    for (casadi_int j=0; j<n_in_; ++j) {
      uout() << "j nnz" << j << ":" << f_.nnz_in(j) << std::endl;
      if (!reduce_in_[j] && f_.nnz_in(j)) {
        if (inst.arg_null.empty()) {
          g << "if (arg1[" << j << "]) arg1[" << j << "]+=" << (vectorize_f() ? 1 : f_.nnz_in(j)) << ";\n";
        } else {
          if (!inst.arg_null[j]) g << "arg1[" << j << "]+=" << (vectorize_f() ? 1 : f_.nnz_in(j)) << ";\n";
        }
      }
    }
    // Update output buffers
    for (casadi_int j=0; j<n_out_; ++j) {
      if (reduce_out_[j]) {
        if (vectorize_f()) {
          g << "res1[" << j << "]+=1;\n";
        } else {
          g << "if (res1[" << j << "]) ";
          g << g.axpy(f_.nnz_out(j), "1.0", "res1[" + str(j) + "]", "res[" + str(j) + "]") << "\n";
        }
      } else {
        if (f_.nnz_out(j)) {
          if (inst.res_null.empty()) {
            g << "if (res1[" << j << "]) ";
            g << "res1[" << j << "]+=" << (vectorize_f() ? 1 : f_.nnz_out(j)) << ";\n";
          } else {
            if (!inst.res_null[j]) g << "res1[" << j << "]+=" << (vectorize_f() ? 1 : f_.nnz_out(j)) << ";\n";
          }
        }
      }
    }

    bool has_rem = n_ != n_padded();
    casadi_int rem = GlobalOptions::vector_width_real-(n_padded()-n_);
    casadi_int n_outer_it = n_padded() / GlobalOptions::vector_width_real; 
    if (outer_loop) {
      g << "}\n";
      if (loop_bug) {
        for (casadi_int j=0; j<n_in_; ++j) {
          if (!reduce_in_[j] && f_.nnz_in(j)) {
            if (inst.arg_null.empty()) {
              g << "if (arg1[" << j << "]) arg2[" << j << "]+=" << GlobalOptions::vector_width_real << ";\n";
            } else {
              if (!inst.arg_null[j]) g << "arg2[" << j << "]+=" << GlobalOptions::vector_width_real << ";\n";
            }
          }
        }
      }
      for (casadi_int j=0; j<n_out_; ++j) {
        if (reduce_out_[j]) {
          if (!loop_bug)
            g << "res1[" << j << "]-= " << GlobalOptions::vector_width_real << " ;\n";
          g.local("sum"+str(j), "casadi_real", "*");
          g << "sum" << j << " = res[" << j << "];\n";
          if (has_rem && n_outer_it>1) {
            // avoid adding nans
            g << "if (i<" << ((n_padded() / GlobalOptions::vector_width_real)-1) << ") {\n";
          }
          if (!has_rem || n_outer_it>1) {
            g.local("k", "casadi_int");
            g.local("j", "casadi_int");
            g << "#pragma omp simd reduction(+:sum" << j << "[:" << f_.nnz_out(j) << "])\n";
            g << "for (k=0; k<" << f_.nnz_out(j) << "; ++k) {\n";
            g << "for (j=0; j<" << GlobalOptions::vector_width_real << "; ++j) ";
            std::string res = loop_bug ? "res2" : "res1";
            g << "sum" << j << "[k]+=" << res << " [" << j << "][j+" << GlobalOptions::vector_width_real << "*k];\n";
            g << "}\n";
          }
          if (has_rem && n_outer_it>1) {
            g << "}\n";
          }
        } else {
          if (loop_bug) {
            if (f_.nnz_out(j)) {
              g << "res2[" << j << "]+=" << GlobalOptions::vector_width_real << ";\n";
            }
          }
        }
      }

    }
    g << "}\n";
    if (has_rem && outer_loop) {
      for (casadi_int j=0; j<n_out_; ++j) {
        if (reduce_out_[j]) {
          g.local("k", "casadi_int");
          g.local("k", "casadi_int");
          g << "for (k=0; k<" << f_.nnz_out(j) << "; ++k) {\n";
          // you might be adding stray numbers
          g << "for (j=0; j<" << rem << "; ++j) ";
          std::string res = loop_bug ? "res2" : "res1";
          g << "sum" << j << "[k]+=" << res << " [" << j << "][j+" << GlobalOptions::vector_width_real << "*k];\n";
          g << "}\n";
        }
      }
    }

  }

  Function MapSum
  ::get_forward(casadi_int nfwd, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Generate map of derivative
    Function df = f_orig_.forward(nfwd);

    for (casadi_int i=0;i<n_out_;++i) {
      if (reduce_out_[i]) casadi_assert(df.nnz_in(n_in_+i)==0, "Case not implemented");
    }

    std::vector<bool> reduce_in = join(reduce_in_, reduce_out_, reduce_in_);
    Function dm = MapSum::create("mapsum" + str(n_) + "_" + df.name(), parallelization(),
      df, n_, reduce_in, reduce_out_);

    // Strip permuting layer when vectorized
    if (dm.is_a("MXFunction")) {
      dm = dm.get_function(dm.get_function()[0]);
    }

    // Input expressions
    vector<MX> arg = dm.mx_in(); // change to layout of reinterpret_layout

    vector<MX> res = arg;
    for (casadi_int i=0;i<f_.n_in();++i) {
      MX& x = res[f_.n_in()+f_.n_out()+i];
      if (!vectorize_f() || reduce_in_[i]) {
        Layout source({df.nnz_in(f_.n_in()+f_.n_out()+i)/nfwd,reduce_in_[i] ? 1 : n_,nfwd});
        Layout target({df.nnz_in(f_.n_in()+f_.n_out()+i)/nfwd, nfwd, reduce_in_[i] ? 1 : n_});
        x = permute_layout(x, Relayout(source, {0, 2, 1}, target, "get_forward_in_"));
      }
    }

    // Get output expressions
    res = dm(res);

    for (casadi_int i=0;i<f_.n_out();++i) {
      if (!vectorize_f() || reduce_out_[i]) {
        Layout source({df.nnz_out(i)/nfwd,nfwd,reduce_out_[i] ? 1 : n_});
        Layout target({df.nnz_out(i)/nfwd,reduce_out_[i] ? 1 : n_,nfwd});
        res[i] = permute_layout(res[i],Relayout(source, {0, 2, 1}, target,"get_forward_out"));
      }
    }

    // Construct return function
    Dict custom_opts = opts;
    custom_opts["always_inline"] = true;
    return Function(name, arg, res, inames, onames, custom_opts);
  }

  bool MapSum::vectorize_f() const {
    return vectorize_f(f_, n_);
  }

  bool MapSum::vectorize_f(const Function& f, casadi_int n) {
    return GlobalOptions::vector_width_real>1 && f.is_a("SXFunction") && n>=GlobalOptions::vector_width_real;
  }

  Function MapSum
  ::get_reverse(casadi_int nadj, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Generate map of derivative
    Function df = f_orig_.reverse(nadj);

    for (casadi_int i=0;i<n_out_;++i) {
      if (reduce_out_[i]) casadi_assert(df.nnz_in(n_in_+i)==0, "Case not implemented");
    }

    std::vector<bool> reduce_in = join(reduce_in_, reduce_out_, reduce_out_);
    Dict options;
    options["dump_in"] = true;
    options["dump_out"] = true;
    options["print_in"] = true;
    Function dm = MapSum::create("mapsum" + str(n_) + "_" + df.name(), parallelization(),
      df, n_, reduce_in, reduce_in_, options);

    // Strip permuting layer when vectorized
    if (dm.is_a("MXFunction")) {
      dm = dm.get_function(dm.get_function()[0]);
    }

    // Input expressions
    vector<MX> arg = dm.mx_in();

    vector<MX> res = arg;
    for (casadi_int i=0;i<f_.n_out();++i) {
      MX& x = res[f_.n_in()+f_.n_out()+i];
      if (!vectorize_f() || reduce_out_[i]) {
        Layout source({df.nnz_in(f_.n_in()+f_.n_out()+i)/nadj,reduce_out_[i] ? 1 : n_,nadj});
        Layout target({df.nnz_in(f_.n_in()+f_.n_out()+i)/nadj,nadj,reduce_out_[i] ? 1 : n_});
        x = permute_layout(x, Relayout(source, {0, 2, 1}, target, "get_forward_in_"));
      }
    }

    // Get output expressions
    res = dm(res);

    for (casadi_int i=0;i<f_.n_in();++i) {
      MX& x = res[i];
      if (!vectorize_f() || reduce_in_[i]) {
        Layout source({df.nnz_out(i)/nadj,nadj,reduce_in_[i] ? 1 : n_});
        Layout target({df.nnz_out(i)/nadj,reduce_in_[i] ? 1 : n_,nadj});
        x = permute_layout(x,Relayout(source, {0, 2, 1}, target,"get_forward_out"));
      }
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
    return local_eval_gen(arg, res, iw, w, m);
  }

} // namespace casadi
