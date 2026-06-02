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


#include "linsol_qr.hpp"
#include "casadi/core/global_options.hpp"

namespace casadi {

  extern "C"
  int CASADI_LINSOL_QR_EXPORT
  casadi_register_linsol_qr(LinsolInternal::Plugin* plugin) {
    plugin->creator = LinsolQr::creator;
    plugin->name = "qr";
    plugin->doc = LinsolQr::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &LinsolQr::options_;
    plugin->deserialize = &LinsolQr::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_QR_EXPORT casadi_load_linsol_qr() {
    LinsolInternal::registerPlugin(casadi_register_linsol_qr);
  }

  LinsolQr::LinsolQr(const std::string& name, const Sparsity& sp)
    : LinsolInternal(name, sp) {
  }

  LinsolQr::~LinsolQr() {
    clear_mem();
  }

  const Options LinsolQr::options_
  = {{&LinsolInternal::options_},
     {{"eps",
       {OT_DOUBLE,
        "Minimum R entry before singularity is declared [1e-12]"}},
      {"cache",
       {OT_DOUBLE,
        "Amount of factorisations to remember (thread-local) [0]"}}
     }
  };

  void LinsolQr::init(const Dict& opts) {
    // Call the init method of the base class
    LinsolInternal::init(opts);

    // Read options
    eps_ = 1e-12;
    n_cache_ = 0;
    for (auto&& op : opts) {
      if (op.first=="eps") {
        eps_ = op.second;
      } else if (op.first=="cache") {
        n_cache_ = op.second;
      }
    }

    // Symbolic factorization
    sp_.qr_sparse(sp_v_, sp_r_, prinv_, pc_);
  }

  // Sign (+1/-1) of a permutation: (-1)^(n - number of cycles)
  static double permutation_sign(const std::vector<casadi_int>& p) {
    std::vector<bool> seen(p.size(), false);
    casadi_int ncycle = 0;
    for (casadi_int i=0; i<static_cast<casadi_int>(p.size()); ++i) {
      if (!seen[i]) {
        ncycle++;
        for (casadi_int j=i; !seen[j]; j=p[j]) seen[j] = true;
      }
    }
    return (static_cast<casadi_int>(p.size())-ncycle) % 2 ? -1 : 1;
  }

  void LinsolQr::finalize() {
    cache_stride_ = sp_.nnz()+sp_v_.nnz()+sp_r_.nnz()+ncol();
    // Determinant of the (data-independent) row and column pivoting
    pivot_sign_ = permutation_sign(prinv_) * permutation_sign(pc_);
    LinsolInternal::finalize();
  }

  int LinsolQr::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<LinsolQrMemory*>(mem);

    // Memory for numerical solution
    m->v.resize(sp_v_.nnz());
    m->r.resize(sp_r_.nnz());
    m->beta.resize(ncol());
    m->w.resize(nrow() + ncol());

    m->cache.resize(cache_stride_*n_cache_);
    m->cache_loc.resize(n_cache_, -1);

    return 0;
  }

  int LinsolQr::sfact(void* mem, const double* A) const {
    return 0;
  }

  int LinsolQr::nfact(void* mem, const double* A) const {
    auto m = static_cast<LinsolQrMemory*>(mem);

    // Check for a cache hit
    double* cache = nullptr;
    bool cache_hit = casadi_cache_check(A, get_ptr(m->cache), get_ptr(m->cache_loc),
      cache_stride_, n_cache_, sp_.nnz(), &cache);

    if (cache && cache_hit) {
      cache += sp_.nnz();
      // Retrieve from cache and return early
      casadi_copy(cache, sp_v_.nnz(), get_ptr(m->v)); cache+=sp_v_.nnz();
      casadi_copy(cache, sp_r_.nnz(), get_ptr(m->r)); cache+=sp_r_.nnz();
      casadi_copy(cache, ncol(), get_ptr(m->beta)); cache+=ncol();
      return 0;
    }

    // Cache miss -> compute result
    casadi_qr(sp_, A, get_ptr(m->w),
              sp_v_, get_ptr(m->v), sp_r_, get_ptr(m->r),
              get_ptr(m->beta), get_ptr(prinv_), get_ptr(pc_));
    // Check singularity
    double rmin;
    casadi_int irmin, nullity;
    nullity = casadi_qr_singular(&rmin, &irmin, get_ptr(m->r), sp_r_, get_ptr(pc_), eps_);
    if (nullity) {
      if (verbose_) {
        print("Singularity detected: Rank %lld<%lld\n", ncol()-nullity, ncol());
        print("First singular R entry: %g<%g, corresponding to row %lld\n", rmin, eps_, irmin);
        casadi_qr_colcomb(get_ptr(m->w), get_ptr(m->r), sp_r_, get_ptr(pc_), eps_, 0);
        print("Linear combination of columns:\n[");
        for (casadi_int k=0; k<ncol(); ++k) print(k==0 ? "%g" : ", %g", m->w[k]);
        print("]\n");
      }
      return 1;
    }

    if (cache) { // Store result in cache
      casadi_copy(A, sp_.nnz(), cache); cache+=sp_.nnz();
      casadi_copy(get_ptr(m->v), sp_v_.nnz(), cache); cache+=sp_v_.nnz();
      casadi_copy(get_ptr(m->r), sp_r_.nnz(), cache); cache+=sp_r_.nnz();
      casadi_copy(get_ptr(m->beta), ncol(), cache); cache+=ncol();
    }
    return 0;
  }

  int LinsolQr::solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    auto m = static_cast<LinsolQrMemory*>(mem);
    casadi_qr_solve(x, nrhs, tr,
                    sp_v_, get_ptr(m->v), sp_r_, get_ptr(m->r),
                    get_ptr(m->beta), get_ptr(prinv_), get_ptr(pc_), get_ptr(m->w));
    return 0;
  }

  double LinsolQr::det(void* mem, const double* A) const {
    auto m = static_cast<LinsolQrMemory*>(mem);
    return pivot_sign_ * casadi_det(sp_v_, get_ptr(m->v), sp_r_, get_ptr(m->r),
                                    get_ptr(m->beta));
  }

  void LinsolQr::generate_factorize(CodeGenerator& g, const std::string& A) const {
    // Codegen the integer vectors
    std::string prinv = g.constant(prinv_);
    std::string pc = g.constant(pc_);
    std::string sp = g.sparsity(sp_);
    std::string sp_v = g.sparsity(sp_v_);
    std::string sp_r = g.sparsity(sp_r_);

    // Carve the factorization workspace from w (reserved by sz_w_fact())
    g.local("qr_v", "casadi_real", "*");
    g << "qr_v = w;\n";
    g.local("qr_r", "casadi_real", "*");
    g << "qr_r = w+" << sp_v_.nnz() << ";\n";
    g.local("qr_beta", "casadi_real", "*");
    g << "qr_beta = w+" << sp_v_.nnz() + sp_r_.nnz() << ";\n";
    g.local("qr_w", "casadi_real", "*");
    g << "qr_w = w+" << sp_v_.nnz() + sp_r_.nnz() + ncol() << ";\n";

    if (n_cache_) {
      // Place the cache in a block to scope its (untouched) local arrays
      g << "{\n";
      g << "casadi_real *c;\n";
      g << "casadi_real cache[" << cache_stride_*n_cache_ << "];\n";
      g << "int cache_loc[" << n_cache_ << "] = {";
      for (casadi_int i=0;i<n_cache_;++i) {
        g << "-1,";
      }
      g << "};\n";
      g << "if (" << g.cache_check(A, "cache", "cache_loc",
        cache_stride_, n_cache_, sp_.nnz(), "&c") << ") {\n";
      casadi_int offset = sp_.nnz();
      g.comment("Retrieve from cache");
      g << g.copy("c+" + str(offset), sp_v_.nnz(), "qr_v") << "\n"; offset+=sp_v_.nnz();
      g << g.copy("c+" + str(offset), sp_r_.nnz(), "qr_r") << "\n"; offset+=sp_r_.nnz();
      g << g.copy("c+" + str(offset), ncol(), "qr_beta") << "\n"; offset+=ncol();
      g << "} else {\n";
    }

    // Factorize
    g << g.qr(sp, A, "qr_w", sp_v, "qr_v", sp_r, "qr_r", "qr_beta", prinv, pc) << "\n";

    if (n_cache_) {
      casadi_int offset = 0;
      g.comment("Store in cache");
      g << g.copy(A, sp_.nnz(), "c") << "\n";; offset+=sp_.nnz();
      g << g.copy("qr_v", sp_v_.nnz(), "c+"+str(offset)) << "\n"; offset+=sp_v_.nnz();
      g << g.copy("qr_r", sp_r_.nnz(), "c+"+str(offset)) << "\n"; offset+=sp_r_.nnz();
      g << g.copy("qr_beta", ncol(), "c+"+str(offset)) << "\n"; offset+=ncol();
      g << "}\n";
      // End of cache block
      g << "}\n";
    }
  }

  void LinsolQr::generate(CodeGenerator& g, const std::string& A, const std::string& x,
                          casadi_int nrhs, bool tr) const {
    generate_factorize(g, A);

    // Solve
    g << g.qr_solve(x, nrhs, tr, g.sparsity(sp_v_), "qr_v", g.sparsity(sp_r_), "qr_r",
                    "qr_beta", g.constant(prinv_), g.constant(pc_), "qr_w") << "\n";
  }

  void LinsolQr::generate_det(CodeGenerator& g, const std::string& A,
                              const std::string& d) const {
    generate_factorize(g, A);

    // Determinant from the factors, corrected by the pivoting sign
    g << d << " = " << (pivot_sign_ < 0 ? "-" : "")
      << g.det(g.sparsity(sp_v_), "qr_v", g.sparsity(sp_r_), "qr_r", "qr_beta") << ";\n";
  }

  LinsolQr::LinsolQr(DeserializingStream& s) : LinsolInternal(s) {
    int version = s.version("LinsolQr", 1, 2);
    s.unpack("LinsolQr::prinv", prinv_);
    s.unpack("LinsolQr::pc", pc_);
    s.unpack("LinsolQr::sp_v", sp_v_);
    s.unpack("LinsolQr::sp_r", sp_r_);
    s.unpack("LinsolQr::eps", eps_);
    if (version>1) {
      s.unpack("LinsolQr::n_cache", n_cache_);
    } else {
      n_cache_ = 1;
    }
  }

  void LinsolQr::serialize_body(SerializingStream &s) const {
    LinsolInternal::serialize_body(s);
    s.version("LinsolQr", 2);
    s.pack("LinsolQr::prinv", prinv_);
    s.pack("LinsolQr::pc", pc_);
    s.pack("LinsolQr::sp_v", sp_v_);
    s.pack("LinsolQr::sp_r", sp_r_);
    s.pack("LinsolQr::eps", eps_);
    s.pack("LinsolQr::n_cache", n_cache_);
  }

} // namespace casadi
