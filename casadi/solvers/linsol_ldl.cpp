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


#include "linsol_ldl.hpp"
#include "casadi/core/global_options.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_LDL_EXPORT
  casadi_register_linsol_ldl(LinsolInternal::Plugin* plugin) {
    plugin->creator = LinsolLdl::creator;
    plugin->name = "ldl";
    plugin->doc = LinsolLdl::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &LinsolLdl::options_;
    plugin->deserialize = &LinsolLdl::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LDL_EXPORT casadi_load_linsol_ldl() {
    LinsolInternal::registerPlugin(casadi_register_linsol_ldl);
  }

  LinsolLdl::LinsolLdl(const std::string& name, const Sparsity& sp)
    : LinsolInternal(name, sp) {
  }

  LinsolLdl::~LinsolLdl() {
    clear_mem();
  }

  const Options LinsolLdl::options_
  = {{&ProtoFunction::options_},
     {{"incomplete",
      {OT_BOOL,
       "Incomplete factorization, without any fill-in"}},
      {"preordering",
       {OT_BOOL,
       "Approximate minimal degree (AMD) preordering"}}
     }
  };

  void LinsolLdl::init(const Dict& opts) {
    // Call the init method of the base class
    LinsolInternal::init(opts);

    // Default options
    incomplete_ = false;
    amd_ = true;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="incomplete") {
        incomplete_ = op.second;
      } else if (op.first=="amd") {
        amd_ = op.second;
      }
    }

    // Symbolic factorization
    if (incomplete_) {
      if (amd_) {
        // Incomplete LDL^T, AMD permutation
        p_ = sp_.amd();
        std::vector<casadi_int> tmp;
        Sparsity Aperm = sp_.sub(p_, p_, tmp);
        sp_Lt_ = triu(Aperm, false);  // no fill-in
      } else {
        p_ = range(sp_.size1());  // no reordering
        sp_Lt_ = triu(sp_, false);  // no fill-in
      }
    } else {
      // Regular LDL^T
      sp_Lt_ = sp_.ldl(p_, amd_);
    }
  }

  int LinsolLdl::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<LinsolLdlMemory*>(mem);

    // Work vectors
    casadi_int nrow = this->nrow();
    m->d.resize(nrow);
    m->l.resize(sp_Lt_.nnz());
    m->w.resize(nrow);

    return 0;
  }

  int LinsolLdl::sfact(void* mem, const double* A) const {
    return 0;
  }

  int LinsolLdl::nfact(void* mem, const double* A) const {
    auto m = static_cast<LinsolLdlMemory*>(mem);
    casadi_ldl(sp_, A, sp_Lt_, get_ptr(m->l), get_ptr(m->d), get_ptr(p_), get_ptr(m->w));
    for (double d : m->d) {
      if (d==0) casadi_warning("LDL factorization has zeros in D");
    }
    return 0;
  }

  int LinsolLdl::solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    auto m = static_cast<LinsolLdlMemory*>(mem);
    casadi_ldl_solve(x, nrhs, sp_Lt_, get_ptr(m->l), get_ptr(m->d), get_ptr(p_), get_ptr(m->w));
    return 0;
  }

  casadi_int LinsolLdl::neig(void* mem, const double* A) const {
    // Count number of negative eigenvalues
    auto m = static_cast<LinsolLdlMemory*>(mem);
    casadi_int nrow = this->nrow();
    casadi_int ret = 0;
    for (casadi_int i=0; i<nrow; ++i) if (m->d[i]<0) ret++;
    return ret;
  }

  casadi_int LinsolLdl::rank(void* mem, const double* A) const {
    // Count number of nonzero eigenvalues
    auto m = static_cast<LinsolLdlMemory*>(mem);
    casadi_int nrow = this->nrow();
    casadi_int ret = 0;
    for (casadi_int i=0; i<nrow; ++i) if (m->d[i]!=0) ret++;
    return ret;
  }

  void LinsolLdl::generate(CodeGenerator& g, const std::string& A, const std::string& x,
                          casadi_int nrhs, bool tr) const {
    // Codegen the integer vectors
    string sp = g.sparsity(sp_);
    string sp_Lt = g.sparsity(sp_Lt_);
    string p = g.constant(p_);

    // Place in block to avoid conflicts caused by local variables
    g << "{\n";
    g.comment("FIXME(@jaeandersson): Memory allocation can be avoided");
    g << "casadi_real lt[" << sp_Lt_.nnz() << "], "
         "d[" << nrow() << "], "
         "w[" << nrow() << "];\n";

    // Factorize
    g << g.ldl(sp, A, sp_Lt, "lt", "d", p, "w") << "\n";

    // Solve
    g << g.ldl_solve(x, nrhs, sp_Lt, "lt", "d", p, "w") << "\n";

    // End of block
    g << "}\n";
  }

  LinsolLdl::LinsolLdl(DeserializingStream& s) : LinsolInternal(s) {
    s.version("LinsolLdl", 1);
    s.unpack("LinsolLdl::p", p_);
    s.unpack("LinsolLdl::sp_Lt", sp_Lt_);
  }

  void LinsolLdl::serialize_body(SerializingStream &s) const {
    LinsolInternal::serialize_body(s);
    s.version("LinsolLdl", 1);
    s.pack("LinsolLdl::p", p_);
    s.pack("LinsolLdl::sp_Lt", sp_Lt_);
  }

} // namespace casadi
