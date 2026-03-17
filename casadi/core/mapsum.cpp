/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

namespace casadi {

  // Build the index vector for a user-to-SIMD perfect shuffle of an (nnz x N) tensor.
  //   User layout:  offset(a, k) = a + k*nnz   (nnz-element a fast, iteration k slow)
  //   SIMD layout:  offset(a, k) = k + a*N      (iteration k fast, nnz-element a slow)
  // Returns perm such that simd[perm[j]] = user[j], i.e.
  //   x.nz(perm) reads from SIMD positions in user order.
  static std::vector<casadi_int> simd_user_perm(casadi_int nnz, casadi_int N) {
    std::vector<casadi_int> perm(nnz * N);
    for (casadi_int k = 0; k < N; ++k)
      for (casadi_int a = 0; a < nnz; ++a)
        perm[a + k * nnz] = k + a * N;
    return perm;
  }

  // Reorder the flat nnz array of x from SIMD to user layout, preserving sparsity.
  // Uses GetNonzeros + sparsity_cast; AD is correct scatter/gather (no 3D decomposition).
  static MX simd_to_user_nz(const MX& x, casadi_int nnz, casadi_int N) {
    return sparsity_cast(x.nz(simd_user_perm(nnz, N)), x.sparsity());
  }

  // Reorder the flat nnz array of x from user to SIMD layout (inverse of above).
  static MX user_to_simd_nz(const MX& x, casadi_int nnz, casadi_int N) {
    auto fwd = simd_user_perm(nnz, N);
    std::vector<casadi_int> inv(fwd.size());
    for (casadi_int j = 0; j < static_cast<casadi_int>(fwd.size()); ++j) inv[fwd[j]] = j;
    return sparsity_cast(x.nz(inv), x.sparsity());
  }

  Function MapSum::create_bare(const std::string& name, const std::string& parallelization,
                          const Function& f, casadi_int n,
                          const std::vector<bool>& reduce_in,
                          const std::vector<bool>& reduce_out,
                          const Dict& opts) {
    if (reduce_out.empty()) return create_bare(name, parallelization, f, n,
                                     reduce_in, std::vector<bool>(f.n_out(), false), opts);
    casadi_assert(reduce_in.size()==f.n_in(), "Dimension mismatch");
    casadi_assert(reduce_out.size()==f.n_out(), "Dimension mismatch");

    Function ret;
    if (parallelization == "serial") {
      std::string suffix = str(reduce_in)+str(reduce_out)+str(opts);
      if (!f->incache(name, ret, suffix)) {
        // Create new serial map
        ret = Function::create(new MapSum(name, f, n, reduce_in, reduce_out), opts);
        casadi_assert_dev(ret.name()==name);
        // Save in cache
        f->tocache_if_missing(ret, suffix);
      }
    } else {
      casadi_error("Unknown parallelization: " + parallelization);
    }

    return ret;
  }

  Function MapSum::create(const std::string& name, const std::string& parallelization,
                          const Function& f, casadi_int n,
                          const std::vector<bool>& reduce_in,
                          const std::vector<bool>& reduce_out,
                          const Dict& opts) {
    Function ret = create_bare(name, parallelization, f, n, reduce_in, reduce_out, opts);

    if (!vectorize_f(f, n)) return ret.wrap_as_needed(ret.name(), opts);

    const MapSum& m = *static_cast<const MapSum*>(ret.get());

    // Input expressions
    std::vector<MX> arg = ret.mx_in();

    // Wrap: user->SIMD on inputs, evaluate, SIMD->user on outputs.
    // nz-reorder (GetNonzeros) instead of permute_layout: AD through GetNonzeros
    // is a well-tested scatter/gather, avoiding the buggy 3D decomposition in
    // PermuteLayout::eval_mx that scrambles block-structured Jacobians.
    std::vector<MX> call_args = arg;
    for (casadi_int i = 0; i < ret.n_in(); ++i)
      if (!m.reduce_in_[i])
        call_args[i] = user_to_simd_nz(arg[i], f.nnz_in(i), n);
    std::vector<MX> res = ret(call_args);
    for (casadi_int i = 0; i < ret.n_out(); ++i)
      if (!m.reduce_out_[i])
        res[i] = simd_to_user_nz(res[i], f.nnz_out(i), n);

    // Construct return function
    Dict custom_opts = opts;
    custom_opts["always_inline"] = true;
    //REMOVE uout() << "opts" << opts << std::endl;
    return Function(ret.name(), arg, res, custom_opts);

  }

  MapSum::MapSum(const std::string& name, const Function& f, casadi_int n,
                 const std::vector<bool>& reduce_in,
                 const std::vector<bool>& reduce_out)
    : FunctionInternal(name), f_(f), n_(n), reduce_in_(reduce_in), reduce_out_(reduce_out) {
    casadi_assert_dev(reduce_in.size()==f.n_in());
    casadi_assert_dev(reduce_out.size()==f.n_out());
    //REMOVE uout() << "init It's map sum!" << std::endl;

    f_orig_ = f;

    if (vectorize_f(f, n_)) {
      std::vector<casadi_int> stride_in(f_.n_in());
      std::vector<casadi_int> stride_out(f_.n_out());
      for (casadi_int j=0; j<f_.n_in(); ++j) {
        stride_in[j] = vectorize_f(f, n_) && !reduce_in[j] ? n_ : 1;
        if (reduce_in_[j]) stride_in[j]*= -1;
      }
      for (casadi_int j=0; j<f_.n_out(); ++j) {
        stride_out[j] = vectorize_f(f, n_) ? (reduce_out_[j] ? GlobalOptions::vector_width_real : n_) : 1;
      }

      Dict opts;
      opts["stride_in"] = stride_in;
      opts["stride_out"] = stride_out;
      f_ = f_->with_options(opts);
    }
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
    return {"f", "f_orig"};
  }

  const Function& MapSum::get_function(const std::string &name) const {
    casadi_assert(has_function(name),
      "No function \"" + name + "\" in " + name_ + ". " +
      "Available functions: " + join(get_function()) + ".");
    if (name=="f") {
      return f_;
    } else {
      return f_orig_;
    }
  }

  void MapSum::find(std::map<FunctionInternal*, std::pair<Function, size_t>> & all_fun,
      casadi_int max_depth) const {
    // Call to base class
    FunctionInternal::find(all_fun, max_depth);
    add_embedded(all_fun, f_, max_depth);
  }

  bool MapSum::has_function(const std::string& fname) const {
    return fname=="f" || fname=="f_orig" ;
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

    // Allocate scratch space for reduced outputs
    // When vectorize_f(), f_ writes at strided offsets (0, vw, 2*vw, ...)
    // so scratch needs nnz_out * vw elements, not just nnz_out
    for (casadi_int j=0;j<n_out_;++j) {
      if (reduce_out_[j]) {
        casadi_int vw = vectorize_f() ? GlobalOptions::vector_width_real : 1;
        alloc_w(f_.nnz_out(j) * vw, true);
      }
    }

    has_refcount_ = f_->has_refcount_;
  }


  // Apply layout permutation to each input: user {nnz,n} -> SIMD {n,nnz}
  std::vector<MX> MapSum::permute_in(const std::vector<MX> & arg, bool invert) const {
    std::vector<MX> ret = arg;
    for (casadi_int i=0;i<f_.n_in();++i) {
      ret[i] = permute_layout(ret[i], permute_in(i, invert)).set_meta("MapSum::permute_in at " + CASADI_WHERE);
    }
    return ret;
  }

  // Apply layout permutation to each output: SIMD {n,nnz} -> user {nnz,n}
  std::vector<MX> MapSum::permute_out(const std::vector<MX> & res, bool invert) const {
    std::vector<MX> ret = res;
    for (casadi_int i=0;i<f_.n_out();++i) {
      ret[i] = permute_layout(ret[i], permute_out(i, invert)).set_meta("MapSum::permute_out at " + CASADI_WHERE);
    }
    return ret;
  }

  // Relayout descriptor for input i: transpose {nnz,n} <-> {n,nnz}
  // Reduced inputs are shared across iterations, so no permutation needed.
  Relayout MapSum::permute_in(casadi_int i, bool invert) const {
    if (reduce_in_[i]) return Relayout(); // reduced: same data for all n iterations
    Layout source({f_.nnz_in(i), n_});
    Layout target({n_, f_.nnz_in(i)});
    Relayout ret = Relayout(source, {1, 0}, target);
    if (invert) return ret.invert();
    return ret;
  }

  // Relayout descriptor for output i: transpose {n,nnz} <-> {nnz,n}
  // Reduced outputs are summed across iterations, so no permutation needed.
  Relayout MapSum::permute_out(casadi_int i, bool invert) const {
    if (reduce_out_[i]) return Relayout(); // reduced: accumulated scalar per output element
    Layout source({n_, f_.nnz_out(i)});
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

  template<typename T1>
  void casadi_inplace_add(T1& y, const T1& x) {
    y += x;
  }

  inline void casadi_inplace_add(bvec_t& y, const bvec_t& x) {
    y |= x;
  }

  template<typename T>
  int MapSum::local_eval_gen(const T** arg, T** res, casadi_int* iw, T* w, int mem) const {

    const T** arg1 = arg+n_in_;
    std::copy_n(arg, n_in_, arg1);
    T** res1 = res+n_out_;

    T* w_scratch = w + f_.sz_w();
    casadi_int vw = vectorize_f() ? GlobalOptions::vector_width_real : 1;
    for (casadi_int j=0;j<n_out_;++j) {
      if (res[j] && reduce_out_[j]) {
        casadi_clear(res[j], f_.nnz_out(j)); // clear sums
        res1[j] = w_scratch; // Make the function dump result in scratch space
        w_scratch += f_.nnz_out(j) * vw;
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
            // Strided add: f_ writes at offsets 0, vw, 2*vw, ...
            for (casadi_int k = 0; k < f_.nnz_out(j); k++) {
              casadi_inplace_add(res[j][k], res1[j][k * vw]);
            }
          } else {
            res1[j] += vectorize_f() ? 1 : f_.nnz_out(j);
          }
        }
      }
    }
    return 0;
  }

  int MapSum::eval_sx(const SXElem** arg, SXElem** res,
      casadi_int* iw, SXElem* w, void* mem,
      bool always_inline, bool never_inline) const {
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

    casadi_int vw = vectorize_f() ? GlobalOptions::vector_width_real : 1;

    bvec_t** arg1 = arg+n_in_;
    std::copy_n(arg, n_in_, arg1);
    bvec_t** res1 = res+n_out_;
    std::copy_n(res, n_out_, res1);

    // For reduced outputs, redirect res1[j] to strided scratch space
    // (f_ may use strided i2 indices; res[j] is non-strided from the caller)
    bvec_t* w_scratch = w + f_.sz_w();
    for (casadi_int j=0;j<n_out_;++j) {
      if (res[j] && reduce_out_[j]) {
        res1[j] = w_scratch;
        w_scratch += f_.nnz_out(j) * vw;
      }
    }

    for (casadi_int i=0; i<n_; ++i) {
      // Scatter non-strided seeds from res[j] into strided scratch res1[j]
      for (casadi_int j=0;j<n_out_;++j) {
        if (res[j] && reduce_out_[j]) {
          casadi_clear(res1[j], f_.nnz_out(j) * vw);
          for (casadi_int k = 0; k < f_.nnz_out(j); k++) {
            res1[j][k * vw] = res[j][k];
          }
        }
      }
      if (f_.rev(arg1, res1, iw, w)) return 1;
      for (casadi_int j=0; j<n_in_; ++j) {
        if (arg1[j] && !reduce_in_[j]) arg1[j] += vectorize_f() ? 1 : f_.nnz_in(j);
      }
      for (casadi_int j=0; j<n_out_; ++j) {
        if (res1[j] && !reduce_out_[j]) res1[j] += vectorize_f() ? 1 : f_.nnz_out(j);
      }
    }
    // Clear consumed seeds (f_.rev cleared the scratch, not the original res[j])
    for (casadi_int j=0;j<n_out_;++j) {
      if (res[j] && reduce_out_[j]) {
        casadi_clear(res[j], f_.nnz_out(j));
      }
    }
    return 0;
  }

  void MapSum::codegen_declarations(CodeGenerator& g, const Instance& inst) const {
    Instance local = inst;
    local.stride_in.resize(f_.n_in());
    for (casadi_int j=0; j<n_in_; ++j) {
      local.stride_in[j] = (vectorize_f() && !reduce_in_[j]) ? n_ : 1;
      if (reduce_in_[j]) local.stride_in[j]*= -1;
    }
    local.stride_out.resize(f_.n_out());
    for (casadi_int j=0; j<n_out_; ++j) {
      local.stride_out[j] = vectorize_f() ? (reduce_out_[j] ? GlobalOptions::vector_width_real : n_) : 1;
    }
    g.add_dependency(f_, local);
  }

  void MapSum::codegen_body(CodeGenerator& g, const Instance& inst) const {
    Instance local = inst;
    local.stride_in.resize(f_.n_in());
    for (casadi_int j=0; j<n_in_; ++j) {
      local.stride_in[j] = (vectorize_f() && !reduce_in_[j]) ? n_ : 1;
      if (reduce_in_[j]) local.stride_in[j]*= -1;
    }
    local.stride_out.resize(f_.n_out());
    for (casadi_int j=0; j<n_out_; ++j) {
      local.stride_out[j] = vectorize_f() ? (reduce_out_[j] ? GlobalOptions::vector_width_real : n_) : 1;
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
          g << "w_scratch+=" << f_.nnz_out(j)* (vectorize_f() ? GlobalOptions::vector_width_real : 1) << ";\n";
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
        // Declare sum pointers for reduced output accumulation
        for (casadi_int j=0; j<n_out_; ++j) {
          if (reduce_out_[j]) {
            g.local("sum"+str(j), "casadi_real", "*");
            g << "sum" << j << " = res[" << j << "];\n";
          }
        }
        g << "for (i=0; i<" << n_ / GlobalOptions::vector_width_real << "; ++i) {\n";
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
        g << "#pragma omp simd safelen(" << n_ << ")\n";
        g << "for (j=0; j<" << n_ << "; ++j) {\n";
      }
    } else {
      g << "for (i=0; i<" << n_ << "; ++i) {\n";
    }

    bool local_increment = vectorize_f() && str(f_).find("SXFunction")!= std::string::npos;

    // Evaluate
    if (local_increment) {
      g << g(f_, "arg1", "res1", "iw", "w", "1", local, "j") << ";\n";
    } else {
      g << "if (" << g(f_, "arg1", "res1", "iw", "w", "1", local) << ") return 1;\n";
    }

   //REMOVE  uout() << "debug" << name_ << std::endl;
    //REMOVE uout() << "n_in_" << n_in_ << std::endl;
    //REMOVE uout() << "reduce_in_" << reduce_in_ << std::endl;
    //REMOVE uout() << "inst.arg_null" << inst.arg_null << std::endl;
   //REMOVE  uout() << "inst.arg_null" << inst.arg_null << std::endl;

    // Update input buffers
    for (casadi_int j=0; j<n_in_; ++j) {
      //REMOVE uout() << "j nnz" << j << ":" << f_.nnz_in(j) << std::endl;
      if (!reduce_in_[j] && f_.nnz_in(j) && !local_increment) {
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
          if (!local_increment) g << "res1[" << j << "]+=1;\n";
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
            if (!local_increment) if (!inst.res_null[j]) g << "res1[" << j << "]+=" << (vectorize_f() ? 1 : f_.nnz_out(j)) << ";\n";
          }
        }
      }
    }

    if (outer_loop) {
      g << "}\n"; // close inner SIMD loop
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
          if (f_.nnz_out(j)==0) continue;
          if (!loop_bug)
            g << "res1[" << j << "]-= " << GlobalOptions::vector_width_real << " ;\n";
          g << "if (sum" << j << ") {\n";
          // Full reduction — all vw lanes are valid (no remainder in this loop)
          g.local("k", "casadi_int");
          g.local("j", "casadi_int");
          g << "#pragma omp simd reduction(+:sum" << j << "[:" << f_.nnz_out(j) << "])\n";
          g << "for (k=0; k<" << f_.nnz_out(j) << "; ++k) {\n";
          g << "for (j=0; j<" << GlobalOptions::vector_width_real << "; ++j) ";
          std::string res = loop_bug ? "res2" : "res1";
          g << "sum" << j << "[k]+=" << res << "[" << j << "][j+" << GlobalOptions::vector_width_real << "*k];\n";
          g << "}\n";
          g << "}\n"; // close if(sum)
        } else {
          if (loop_bug) {
            if (f_.nnz_out(j)) {
              g << "res2[" << j << "]+=" << GlobalOptions::vector_width_real << ";\n";
            }
          }
        }
      }
    }
    g << "}\n"; // close outer/main loop

    // Scalar tail for remainder iterations
    casadi_int rem = n_ % GlobalOptions::vector_width_real;
    if (rem > 0 && outer_loop) {
      // Sync pointers after outer loop (loop_bug: arg1/res1 are stale copies of arg2/res2)
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
      g << "for (j=0; j<" << rem << "; ++j) {\n";
      // Kernel call (same as SIMD path, but scalar)
      if (local_increment) {
        g << g(f_, "arg1", "res1", "iw", "w", "1", local, "j") << ";\n";
      } else {
        g << "if (" << g(f_, "arg1", "res1", "iw", "w", "1", local) << ") return 1;\n";
      }
      // Update input buffers
      for (casadi_int j=0; j<n_in_; ++j) {
        if (!reduce_in_[j] && f_.nnz_in(j) && !local_increment) {
          if (inst.arg_null.empty()) {
            g << "if (arg1[" << j << "]) arg1[" << j << "]+=1;\n";
          } else {
            if (!inst.arg_null[j]) g << "arg1[" << j << "]+=1;\n";
          }
        }
      }
      // Update output buffers
      for (casadi_int j=0; j<n_out_; ++j) {
        if (reduce_out_[j]) {
          if (!local_increment) g << "res1[" << j << "]+=1;\n";
        } else {
          if (f_.nnz_out(j)) {
            if (inst.res_null.empty()) {
              g << "if (res1[" << j << "]) res1[" << j << "]+=1;\n";
            } else {
              if (!local_increment) if (!inst.res_null[j]) g << "res1[" << j << "]+=1;\n";
            }
          }
        }
      }
      g << "}\n"; // close scalar tail loop
      // Reduce remainder lanes from scratch to sum
      for (casadi_int j=0; j<n_out_; ++j) {
        if (reduce_out_[j]) {
          g << "if (sum" << j << ") {\n";
          g.local("k", "casadi_int");
          g << "for (k=0; k<" << f_.nnz_out(j) << "; ++k) {\n";
          g << "for (j=0; j<" << rem << "; ++j) ";
          std::string res = loop_bug ? "res2" : "res1";
          g << "sum" << j << "[k]+=" << res << "[" << j << "][j+" << GlobalOptions::vector_width_real << "*k];\n";
          g << "}\n";
          g << "}\n";
        }
      }
    }

  }


  void MapSum::codegen_incref(CodeGenerator& g, const Instance& inst) const {
    auto i = g.incref_added_.insert(f_.get());
    if (i.second) { // prevent duplicate calls
      g << f_->codegen_name(g, inst) << "_incref();\n";
    }
  }

  void MapSum::codegen_decref(CodeGenerator& g, const Instance& inst) const {
    auto i = g.decref_added_.insert(f_.get());
    if (i.second) { // prevent duplicate calls
      g << f_->codegen_name(g, inst) << "_decref();\n";
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
    Function dm = MapSum::create_bare("mapsum" + str(n_) + "_" + df.name(), parallelization(),
      df, n_, reduce_in, reduce_out_);

    // Input expressions
    std::vector<MX> arg = dm.mx_in(); // change to layout of reinterpret_layout

    std::vector<MX> res = arg;
    for (casadi_int i=0;i<f_.n_in();++i) {
      MX& x = res[f_.n_in()+f_.n_out()+i];
      if (!vectorize_f() || reduce_in_[i]) {
        Layout source({df.nnz_in(f_.n_in()+f_.n_out()+i)/nfwd,reduce_in_[i] ? 1 : n_,nfwd});
        Layout target({df.nnz_in(f_.n_in()+f_.n_out()+i)/nfwd, nfwd, reduce_in_[i] ? 1 : n_});
        x = permute_layout(x, Relayout(source, {0, 2, 1}, target, "get_forward_in_")).set_meta("MapSum::get_forward fwd_in at " + CASADI_WHERE);
      }
    }

    // Get output expressions
    res = dm(res);

    for (casadi_int i=0;i<f_.n_out();++i) {
      if (!vectorize_f() || reduce_out_[i]) {
        Layout source({df.nnz_out(i)/nfwd,nfwd,reduce_out_[i] ? 1 : n_});
        Layout target({df.nnz_out(i)/nfwd,reduce_out_[i] ? 1 : n_,nfwd});
        res[i] = permute_layout(res[i],Relayout(source, {0, 2, 1}, target,"get_forward_out")).set_meta("MapSum::get_forward fwd_out at " + CASADI_WHERE);
      }
    }

    // Construct return function
    Dict custom_opts = opts;
    custom_opts["always_inline"] = true;
    custom_opts["allow_duplicate_io_names"] = true;
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
    Function dm = MapSum::create_bare("mapsum" + str(n_) + "_" + df.name(), parallelization(),
      df, n_, reduce_in, reduce_in_, options);

    // Input expressions
    std::vector<MX> arg = dm.mx_in();

    std::vector<MX> res = arg;
    for (casadi_int i=0;i<f_.n_out();++i) {
      MX& x = res[f_.n_in()+f_.n_out()+i];
      if (!vectorize_f() || reduce_out_[i]) {
        Layout source({df.nnz_in(f_.n_in()+f_.n_out()+i)/nadj,reduce_out_[i] ? 1 : n_,nadj});
        Layout target({df.nnz_in(f_.n_in()+f_.n_out()+i)/nadj,nadj,reduce_out_[i] ? 1 : n_});
        x = permute_layout(x, Relayout(source, {0, 2, 1}, target, "get_forward_in_")).set_meta("MapSum::get_reverse adj_in at " + CASADI_WHERE);
      }
    }

    // Get output expressions
    res = dm(res);

    for (casadi_int i=0;i<f_.n_in();++i) {
      MX& x = res[i];
      if (!vectorize_f() || reduce_in_[i]) {
        Layout source({df.nnz_out(i)/nadj,nadj,reduce_in_[i] ? 1 : n_});
        Layout target({df.nnz_out(i)/nadj,reduce_in_[i] ? 1 : n_,nadj});
        x = permute_layout(x,Relayout(source, {0, 2, 1}, target,"get_forward_out")).set_meta("MapSum::get_reverse adj_out at " + CASADI_WHERE);
      }
    }

    // Construct return function
    options = opts;
    options["always_inline"] = true;
    options["allow_duplicate_io_names"] = true;
    return Function(name, arg, res, inames, onames, options);
  }

  Function MapSum::get_jacobian(const std::string& name,
                      const std::vector<std::string>& inames,
                      const std::vector<std::string>& onames,
                      const Dict& opts) const {

    Dict options = opts;
    // Matching get_forward/get_reverse: use f_orig_ (no strides) + create_bare
    // (no wrapping). The outer create()'s nz-reorder differentiates cleanly
    // via GetNonzeros scatter/gather.
    Function Jf = f_orig_.jacobian();

    std::vector<bool> reduce_in = reduce_in_;
    for (size_t oind = 0; oind < n_out_; ++oind) { // Nominal outputs
      reduce_in.push_back(reduce_out_[oind]);
    }
    std::vector<bool> reduce_out;
    reduce_out.reserve(n_in_ * n_out_);
    for (size_t oind = 0; oind < n_out_; ++oind) {
      for (size_t iind = 0; iind < n_in_; ++iind) {
        reduce_out.push_back(reduce_out_[oind] && reduce_in_[iind]);
      }
    }

    Function Jmap = MapSum::create_bare("mapsum" + str(n_) + "_" + Jf.name(), parallelization(),
      Jf, n_, reduce_in, reduce_out);

    // Input expressions
    std::vector<MX> arg = Jmap.mx_in();

    std::vector<MX> res = Jmap(arg);

    // Helper: build the inverse of simd_user_perm (user->SIMD index vector)
    // for use with r(perm, ...) row/column selection.
    auto inv_perm = [](casadi_int nnz, casadi_int N) {
      auto fwd = simd_user_perm(nnz, N);
      std::vector<casadi_int> inv(fwd.size());
      for (size_t j = 0; j < fwd.size(); ++j) inv[fwd[j]] = j;
      return inv;
    };

    size_t i = 0;
    // Post-process each Jacobian block J_{oind,iind}.
    //
    // The bare Jmap (from create_bare) stores non-reduced blocks in SIMD nnz
    // layout.  Three things may need fixing:
    //
    //   (a) simd_to_user_nz: fix nnz so structural operations (sparsity_cast,
    //       horzsplit) assign values to correct block positions.
    //
    //   (b) sparsity_cast / horzsplit: build block-diagonal or stacked matrix
    //       (only for non-reduced FUNCTION outputs; same as non-vectorized path).
    //
    //   (c) Row/column reindex to SIMD order: the outer create()'s nz-reorder
    //       chain rule applies simd_to_user to tangent vectors.  We reindex
    //       the dimensions that have a corresponding nz-reorder in create()
    //       so the chain rule's conversion yields the correct final user order.
    //       - Rows:    only for non-reduced outputs  (!reduce_out_[oind])
    //       - Columns: only for non-reduced inputs   (!reduce_in_[iind])
    //
    for (size_t oind = 0; oind < n_out_; ++oind) {
      for (size_t iind = 0; iind < n_in_; ++iind) {
        MX& r = res[i];

        // (a) SIMD -> user nnz layout (for any non-reduced Jmap block)
        if (vectorize_f() && !reduce_out[i])
          r = simd_to_user_nz(r, Jf.nnz_out(i), n_);

        // (b) Build block-diagonal structure (non-reduced function outputs only)
        if (!reduce_out_[oind]) {
          if (reduce_in_[iind]) {
            r = vertcat(horzsplit(r, Jf.size2_out(i)));
          } else {
            r = sparsity_cast(r, Sparsity::kron(Sparsity::diag(n_), Jf.sparsity_out(i)));
          }
        }

        // (c) Reindex dimensions that have outer nz-reorder to SIMD order
        if (vectorize_f() && !reduce_out[i]) {
          bool need_rows = !reduce_out_[oind];
          bool need_cols = !reduce_in_[iind];
          if (need_rows && need_cols) {
            r = r(inv_perm(f_orig_.nnz_out(oind), n_),
                  inv_perm(f_orig_.nnz_in(iind), n_));
          } else if (need_rows) {
            r = r(inv_perm(f_orig_.nnz_out(oind), n_), Slice());
          } else if (need_cols) {
            r = r(Slice(), inv_perm(f_orig_.nnz_in(iind), n_));
          }
        }

        i++;
      }
    }

    // Construct return function
    Dict custom_opts = options;
    custom_opts["always_inline"] = true;
    custom_opts["allow_duplicate_io_names"] = true;
    return Function(name, arg, res, inames, onames, custom_opts);
  }

  int MapSum::eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    // This checkout/release dance is an optimization.
    // Could also use the thread-safe variant f_(arg1, res1, iw, w)
    // in Map::eval_gen
    setup(mem, arg, res, iw, w);
    scoped_checkout<Function> m(f_);
    return local_eval_gen(arg, res, iw, w, m);
  }

  Function MapSum::pull_out(const std::vector<casadi_int>& in, Function& outer) const {
    uout() << "pull_out" << name_ << ":" << in << n_in_ << reduce_in_ << vector_slice(reduce_in_, in) << std::endl;
    //casadi_assert_dev(all(vector_slice(reduce_in_, in)));
    Function f_core = f_.pull_out(in, outer);


    std::vector<bool> reduce_in = vector_slice(reduce_in_, in, true);
    reduce_in.push_back(false);

    return f_core.map(n_, reduce_in, reduce_out_);
  }

} // namespace casadi
