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


#include "map.hpp"
#include "serializing_stream.hpp"
#include "global_options.hpp"

#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.thread.h>
#else // CASADI_WITH_THREAD_MINGW
#include <thread>
#endif // CASADI_WITH_THREAD_MINGW
#endif // CASADI_WITH_THREAD

namespace casadi {

  Function Map::create(const std::string& parallelization, const Function& f, casadi_int n) {
    Function ret;
    // Create instance of the right class
    std::string suffix = str(n) + "_" + f.name();
    if (parallelization == "serial") {
      ret = Function::create(new Map("map" + suffix, f, n), Dict());
    } else if (parallelization== "openmp") {
      ret = Function::create(new OmpMap("ompmap" + suffix, f, n), Dict());
    } else if (parallelization== "thread") {
      ret = Function::create(new ThreadMap("threadmap" + suffix, f, n), Dict());
    } else {
      casadi_error("Unknown parallelization: " + parallelization);
    }

    if (!vectorize_f(f, n)) return ret;

    const Map& m = *static_cast<const Map*>(ret.get());


    // Input expressions
    std::vector<MX> arg = ret.mx_in(); // change to layout of reinterpret_layout

    std::vector<MX> res = m.permute_out(ret(m.permute_in(arg)));

    // Construct return function
    Dict custom_opts;
    custom_opts["always_inline"] = true;
    Function retf(ret.name(), arg, res, custom_opts);
    static std::vector<Function> leaking;
    leaking.push_back(retf);
    return retf;
  }

  std::vector<MX> Map::permute_in(const std::vector<MX> & arg, bool invert) const {
    std::vector<MX> ret(arg.size());
    for (casadi_int i=0;i<f_.n_in();++i) {
      ret[i] = permute_layout(arg[i], permute_in(i, invert));
    }
    //REMOVE uout() << "testje" << ret << std::endl;
    return ret;
  }

  std::vector<MX> Map::permute_out(const std::vector<MX> & res, bool invert) const {
    std::vector<MX> ret(res.size());
    for (casadi_int i=0;i<f_.n_out();++i) {
      ret[i] = permute_layout(res[i], permute_out(i, invert));
    }
    return ret;
  }

  Relayout Map::permute_in(casadi_int i, bool invert) const {
    Layout source({f_.nnz_in(i), n_});
    Layout target({n_, f_.nnz_in(i)}, {n_padded(), f_.nnz_in(i)});
    Relayout ret = Relayout(source, {1, 0}, target);
    if (invert) return ret.invert();
    return ret;
  }

  Relayout Map::permute_out(casadi_int i, bool invert) const {
    Layout source({n_, f_.nnz_out(i)}, {n_padded(), f_.nnz_out(i)});
    Layout target({f_.nnz_out(i), n_});
    Relayout ret = Relayout(source, {1, 0}, target);
    if (invert) return ret.invert();
    return ret;
  }

  Map::Map(const std::string& name, const Function& f, casadi_int n)
    : FunctionInternal(name), f_(f), n_(n) {
    //REMOVE uout() << "It's map!" << std::endl;

    f_orig_ = f;
    if (vectorize_f(f, n_)) {

      std::vector<casadi_int> stride_in(f_.n_in());
      std::vector<casadi_int> stride_out(f_.n_out());
      for (casadi_int j=0; j<f_.n_in(); ++j) {
        stride_in[j] = vectorize_f(f, n_) ? n_padded(n) : 1;
      }
      for (casadi_int j=0; j<f_.n_out(); ++j) {
        stride_out[j] = vectorize_f(f, n_) ? n_padded(n) : 1;
      }

      Dict opts;
      opts["stride_in"] = stride_in;
      opts["stride_out"] = stride_out;
      f_ = f_->with_options(opts);

      // Options
      /*Dict my_opts = f->generate_options();
      my_opts["stride_in"] = stride_in;
      my_opts["stride_out"] = stride_out;
      // Wrap the function
      std::vector<SX> arg = f.sx_in();
      std::vector<SX> res = f(arg);
      f_ = Function(f.name(), arg, res, f.name_in(), f.name_out(), my_opts);*/
    }
  }

  bool Map::is_a(const std::string& type, bool recursive) const {
    return type=="Map"
      || (recursive && FunctionInternal::is_a(type, recursive));
  }

  bool OmpMap::is_a(const std::string& type, bool recursive) const {
    return type=="OmpMap"
      || (recursive && Map::is_a(type, recursive));
  }

  bool ThreadMap::is_a(const std::string& type, bool recursive) const {
    return type=="ThreadMap"
      || (recursive && Map::is_a(type, recursive));
  }

 std::vector<std::string> Map::get_function() const {
    return {"f"};
  }

  const Function& Map::get_function(const std::string &name) const {
    casadi_assert(has_function(name),
      "No function \"" + name + "\" in " + name_ + ". " +
      "Available functions: " + join(get_function()) + ".");
    return f_;
  }

  bool Map::has_function(const std::string& fname) const {
    return fname=="f";
  }

  bool Map::vectorize_f() const {
    return vectorize_f(f_, n_);
  }

  bool Map::vectorize_f(const Function& f, casadi_int n) {
    return GlobalOptions::vector_width_real>1 && f.is_a("SXFunction") && n>=GlobalOptions::vector_width_real;
  }

  Layout Map::get_layout_in(casadi_int i) {
    if (!vectorize_f()) return Layout();
    return Layout({f_.nnz_in(i), n_}, {f_.nnz_in(i), n_padded()});
  }

  Layout Map::get_layout_out(casadi_int i) {
    if (!vectorize_f()) return Layout();
    // Just here to allocate enough sz_self for outputs
    return Layout({f_.nnz_out(i), n_}, {f_.nnz_out(i), n_padded()});
  }

  void Map::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);
    s.pack("Map::f", f_);
    s.pack("Map::f_orig", f_orig_);
    s.pack("Map::n", n_);
  }

  void Map::serialize_type(SerializingStream &s) const {
    FunctionInternal::serialize_type(s);
    s.pack("Map::class_name", class_name());
  }

  Map::Map(DeserializingStream& s) : FunctionInternal(s) {
    s.unpack("Map::f", f_);
    s.unpack("Map::f_orig", f_orig_);
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
    std::copy_n(arg, n_in_, arg1);
    T** res1 = res+n_out_;
    std::copy_n(res, n_out_, res1);
    for (casadi_int i=0; i<n_; ++i) {
      if (f_(arg1, res1, iw, w, mem)) return 1;
      for (casadi_int j=0; j<n_in_; ++j) {
        casadi_int stride = vectorize_f() ? 1 : f_.nnz_in(j);
        if (arg1[j]) arg1[j] += stride;
      }
      for (casadi_int j=0; j<n_out_; ++j) {
        casadi_int stride = vectorize_f() ? 1 : f_.nnz_out(j);
        if (res1[j]) res1[j] += stride;
      }
    }
    return 0;
  }

  int Map::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w, void* mem,
      bool always_inline, bool never_inline) const {
    return eval_gen(arg, res, iw, w);
  }

  int Map::sp_forward(const bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    return eval_gen(arg, res, iw, w);
  }

  int Map::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    bvec_t** arg1 = arg+n_in_;
    std::copy_n(arg, n_in_, arg1);
    bvec_t** res1 = res+n_out_;
    std::copy_n(res, n_out_, res1);
    for (casadi_int i=0; i<n_; ++i) {
      if (f_.rev(arg1, res1, iw, w)) return 1;
      for (casadi_int j=0; j<n_in_; ++j) {
        casadi_int stride = vectorize_f() ? 1 : f_.nnz_in(j);
        if (arg1[j]) arg1[j] += stride;
      }
      for (casadi_int j=0; j<n_out_; ++j) {
        casadi_int stride = vectorize_f() ? 1 : f_.nnz_out(j);
        if (res1[j]) res1[j] += stride;
      }
    }
    return 0;
  }

  void Map::codegen_declarations(CodeGenerator& g, const Instance& inst) const {
    Instance local = inst;
    local.stride_in.resize(f_.n_in());
    for (casadi_int j=0; j<n_in_; ++j) {
      local.stride_in[j] = vectorize_f() ? n_padded() : 1;
    }
    local.stride_out.resize(f_.n_out());
    for (casadi_int j=0; j<n_out_; ++j) {
      local.stride_out[j] = vectorize_f() ? n_padded() : 1;
    }
    g.add_dependency(f_, local);
    //REMOVE uout() << "codegen_dec" << local.arg_null << f_ << std::endl;
  }

  void Map::codegen_body(CodeGenerator& g,
      const Instance& inst) const {
    Instance local = inst;
    local.stride_in.resize(f_.n_in());
    for (casadi_int j=0; j<n_in_; ++j) {
      local.stride_in[j] = vectorize_f() ? n_padded() : 1;
    }
    local.stride_out.resize(f_.n_out());
    for (casadi_int j=0; j<n_out_; ++j) {
      local.stride_out[j] = vectorize_f() ? n_padded() : 1;
    }
    //REMOVE uout() << "codegen_body" << local.arg_null << f_ << std::endl;
    g.local("i", "casadi_int");
    g.local("arg1[" + str(f_.sz_arg()) + "]", "const casadi_real*");
    g.local("res1[" + str(f_.sz_res()) + "]", "casadi_real*");

    // Input buffer
    g << "for (i=0; i<" << n_in_ << "; ++i) arg1[i]=arg[i];\n";
    // Output buffer
    g << "for (i=0; i<" << n_out_ << "; ++i) res1[i]=res[i];\n";
    if (vectorize_f()) {
      g << "#pragma omp simd\n";
    }
    g << "for (i=0; i<" << (vectorize_f() ? n_padded() : n_) << "; ++i) {\n";
    // Evaluate
    if (str(f_).find("SXFunction")!= std::string::npos) {
      g << g(f_, "arg1", "res1", "iw", "w", "1", local) << ";\n";
    } else {
      g << "if (" << g(f_, "arg1", "res1", "iw", "w", "1", local) << ") return 1;\n";
    }

    // Update input buffers
    for (casadi_int j=0; j<n_in_; ++j) {
      if (f_.nnz_in(j)) {
        casadi_int stride = vectorize_f() ? 1 : f_.nnz_in(j);
        if (inst.arg_null.empty()) {
          g << "if (arg1[" << j << "]) arg1[" << j << "]+=" << stride << ";\n";
        } else {
          if (!inst.arg_null[j]) g << "arg1[" << j << "]+=" << stride << ";\n";
        }
      }
    }
    // Update output buffers
    for (casadi_int j=0; j<n_out_; ++j) {
      if (f_.nnz_out(j)) {
        casadi_int stride = vectorize_f() ? 1 : f_.nnz_out(j);
        if (inst.res_null.empty()) {
          g << "if (res1[" << j << "]) res1[" << j << "]+=" << stride << ";\n";
        } else {
          if (!inst.res_null[j]) g << "res1[" << j << "]+=" << stride << ";\n";
        }
      }
    }
    g << "}\n";
  }

  casadi_int Map::n_padded() const {
    if (vectorize_f()) {
      return n_padded(n_);
    }
    return n_;
  }

  casadi_int Map::n_padded(casadi_int n) {
    return n;
    casadi_int rem = n % GlobalOptions::vector_width_real;
    if (rem==0) return n;
    return n + GlobalOptions::vector_width_real - rem;
  }

  Function Map
  ::get_forward(casadi_int nfwd, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Generate map of derivative
    casadi_assert_dev(!f_orig_.is_null());
    Function df = f_orig_.forward(nfwd);
    Function dm = df.map(n_, parallelization());

    // Strip permuting layer when vectorized
    if (dm.is_a("MXFunction")) {
      dm = dm.get_function(dm.get_function()[0]);
    }

    // Input expressions
    std::vector<MX> arg = dm.mx_in(); // change to layout of reinterpret_layout

    std::vector<MX> res = arg;
    for (casadi_int i=0;i<f_.n_in();++i) {
      MX& x = res[f_.n_in()+f_.n_out()+i];
      casadi_int df_in = df.nnz_in(f_.n_in()+f_.n_out()+i);
      if (false && vectorize_f()) {
        Layout source({n_, df_in}, {n_padded(), df_in});
        Layout target({df_in, n_});
        Relayout ret = Relayout(source, {1, 0}, target);
        //REMOVE uout() << "FOO " << i << std::endl;
        x = permute_layout(x, ret);

        // Layout source({f_.nnz_in(i), n_});
        // Layout target({n_, f_.nnz_in(i)}, {n_padded(), f_.nnz_in(i)});
      }
      if (!vectorize_f()) {
        //REMOVE uout() << "here" << std::endl;
        Layout source({df_in/nfwd,n_,nfwd});
        Layout target({df_in/nfwd, nfwd, n_});
        x = permute_layout(x, Relayout(source, {0, 2, 1}, target, "get_forward_in_"));
        //REMOVE uout() << "baz" << x << "source" << source << "target" << target << std::endl;
      }
    }

    // Get output expressions
    res = dm(res);

    for (casadi_int i=0;i<f_.n_out();++i) {
      casadi_int df_out = df.nnz_out(i);
      MX& x = res[i];
      if (false && vectorize_f()) {

      //Layout source({n_, f_.nnz_out(i)}, {n_padded(), f_.nnz_out(i)});
      //Layout target({f_.nnz_out(i), n_});
      //  Relayout ret = Relayout(source, {1, 0}, target);
        Layout source({df_out, n_});
        Layout target({n_, df_out}, {n_padded(), f_.nnz_out()});
        Relayout ret = Relayout(source, {1, 0}, target);
        //REMOVE uout() << "BAR " << i << std::endl;
        x = permute_layout(x, ret);
      }
      if (!vectorize_f()) {
        //REMOVE uout() << "there" << std::endl;
        Layout source({df_out/nfwd,nfwd,n_});
        Layout target({df_out/nfwd,n_,nfwd});
        x = permute_layout(x,Relayout(source, {0, 2, 1}, target,"get_forward_out"));
       //REMOVE  uout() << "baz" << res[i] << "source" << source << "target" << target << std::endl;
      }
    }

    Dict options = opts;
    options["allow_duplicate_io_names"] = true;

    // Construct return function
    options["always_inline"] = true;
    Function retf(name, arg, res, inames, onames, options);
    static std::vector<Function> leaking;
    leaking.push_back(retf);
    return retf;
  }

  Function Map
  ::get_reverse(casadi_int nadj, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Generate map of derivative
    Function df = f_orig_.reverse(nadj);
    Function dm = df.map(n_, parallelization());

    // Strip permuting layer when vectorized
    if (dm.is_a("MXFunction")) {
      dm = dm.get_function(dm.get_function()[0]);
    }

    // Input expressions
    std::vector<MX> arg = dm.mx_in();

    std::vector<MX> res = arg;
    for (casadi_int i=0;i<f_.n_out();++i) {
      std::vector<casadi_int> dims = {df.nnz_in(f_.n_in()+f_.n_out()+i)/nadj,n_,nadj};
      std::vector<casadi_int> dims_in = {df.nnz_in(f_.n_in()+f_.n_out()+i)/nadj,nadj,n_};
      MX& x = res[f_.n_in()+f_.n_out()+i];
      if (!vectorize_f()) {
        Layout source({df.nnz_in(f_.n_in()+f_.n_out()+i)/nadj,n_,nadj});
        Layout target({df.nnz_in(f_.n_in()+f_.n_out()+i)/nadj,nadj,n_});
        x = permute_layout(x, Relayout(source, {0, 2, 1}, target, "get_forward_in_"));
      }
    }

    // Get output expressions
    res = dm(res);

    for (casadi_int i=0;i<f_.n_in();++i) {
      MX& x = res[i];
      if (!vectorize_f()) {
        Layout source({df.nnz_out(i)/nadj,nadj,n_});
        Layout target({df.nnz_out(i)/nadj,n_,nadj});
        x = permute_layout(x,Relayout(source, {0, 2, 1}, target,"get_forward_out"));
      }
    }

    Dict options = opts;
    options["allow_duplicate_io_names"] = true;
    options["always_inline"] = true;
    // Construct return function
    Function retf(name, arg, res, inames, onames, options);
    static std::vector<Function> leaking;
    leaking.push_back(retf);
    return retf;
  }

  Function Map::get_jacobian(const std::string& name,
                                  const std::vector<std::string>& inames,
                                  const std::vector<std::string>& onames,
                                  const Dict& opts) const {

// Generate map of derivative
    Function Jf = f_.jacobian();

    Function Jmap = Jf.map(n_);

    // Input expressions
    std::vector<MX> arg = Jmap.mx_in();

    std::vector<MX> res = Jmap(arg);

    size_t i=0;
    for (size_t oind = 0; oind < n_out_; ++oind) {
      for (size_t iind = 0; iind < n_in_; ++iind) {
        MX& r = res[i];
        r = sparsity_cast(r, Sparsity::kron(Sparsity::diag(n_), Jf.sparsity_out(i)));
        i++;
      }
    }

    // Construct return function
    Dict custom_opts = opts;
    custom_opts["always_inline"] = true;
    return Function(name, arg, res, inames, onames, custom_opts);
  }

  int Map::eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    // This checkout/release dance is an optimization.
    // Could also use the thread-safe variant f_(arg1, res1, iw, w)
    // in Map::eval_gen
    setup(mem, arg, res, iw, w);
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
    setup(mem, arg, res, iw, w);
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
        << g.arg(j, true) << "+i*" << f_.nnz_in(j) << ": 0;\n";
    }
    g << "res1 = res + " <<  n_out_ << "+i*" <<  sz_res << ";\n";
    for (casadi_int j=0; j<n_out_; ++j) {
      g << "res1[" << j << "] = res[" << j << "] ?"
        << g.res(j,true) << "+i*" << f_.nnz_out(j) << ": 0;\n";
    }
    g << "flag = "
      << g(f_, "arg1", "res1", "iw+i*" + str(sz_iw), "w+i*" + str(sz_w), "1", inst) << " || flag;\n"
      << "}\n"
      << "if (flag) return 1;\n";
  }

  void OmpMap::init(const Dict& opts) {
#ifndef WITH_OPENMP
    casadi_warning("CasADi was not compiled with WITH_OPENMP=ON. "
                   "Falling back to serial evaluation.");
#endif // WITH_OPENMP
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
    setup(mem, arg, res, iw, w);
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
