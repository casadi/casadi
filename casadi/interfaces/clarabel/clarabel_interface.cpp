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

#include "clarabel_interface.hpp"
#include "casadi/core/nlp_tools.hpp"
#include <clarabel_runtime_str.h>  // analogous to highs_runtime_str.h
namespace casadi {



extern "C"
int CASADI_CONIC_CLARABEL_EXPORT
casadi_register_conic_clarabel(Conic::Plugin* plugin) {
  plugin->creator = ClarabelInterface::creator;
  plugin->name = "clarabel";
  plugin->doc = ClarabelInterface::meta_doc.c_str();
  plugin->version = CASADI_VERSION;
  plugin->options = &ClarabelInterface::options_;
  plugin->deserialize = &ClarabelInterface::deserialize;
  return 0;
}

extern "C"
void CASADI_CONIC_CLARABEL_EXPORT casadi_load_conic_clarabel() {
  Conic::registerPlugin(casadi_register_conic_clarabel);
}


ClarabelInterface::ClarabelInterface(const std::string& name,
                                     const std::map<std::string, Sparsity>& st)
  : Conic(name, st) {
}

const Options ClarabelInterface::options_
  = {{&Conic::options_},
     {{"clarabel",
       {OT_DICT,
        "Options to be passed to Clarabel."}},
     }
   };

void ClarabelInterface::init(const Dict& opts) {
  // Call the init method of the base class
  Conic::init(opts);
  // Read user options for Clarabel (if provided)
  for (auto&& op : opts) {
    if (op.first=="clarabel") {
      opts_ = op.second;
    }
  }
  // Initialize dependent (read-only) members
  init_dependent();
  set_clarabel_prob();

  // Allocate memory for arguments and work vectors; the work sizes are determined
  // by a call to our clarabel_work routine.
  casadi_int sz_arg, sz_res, sz_w, sz_iw;
  casadi_clarabel_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);

  alloc_arg(sz_arg, true);
  alloc_res(sz_res, true);
  alloc_iw(sz_iw, true);
  alloc_w(sz_w, true);
}

void ClarabelInterface::init_dependent() {
  // For Clarabel we need two CSC representations:
  // one for the quadratic cost (matrix P, built here from H_) and one for the constraints (A_).
  // (We assume that H_ and A_ have been set up by the base class.)
  colindp_.resize(H_.size2()+1);
  rowp_.resize(H_.nnz());
  colinda_.resize(A_.size2()+1);
  rowa_.resize(A_.nnz());
  copy_vector(H_.colind(), colindp_);
  copy_vector(H_.row(), rowp_);
  copy_vector(A_.colind(), colinda_);
  copy_vector(A_.row(), rowa_);
  // (We omit integrality because in this example Clarabel does not support integer vars.)
}

void codegen_local(CodeGenerator& g, const std::string& name, const std::vector<int>& v) {
  std::string n = name + "[]";
  g.local(n, "static const int");
  std::stringstream init;
  init << "{";
  for (casadi_int i = 0; i < v.size(); ++i) {
    init << v[i];
    if (i < v.size()-1) init << ", ";
  }
  if (v.empty()) init << "0"; // ISO C forbids empty initializer braces
  init << "}";
  g.init_local(n, init.str());
}

void ClarabelInterface::set_clarabel_prob(CodeGenerator& g) const {
  g << "p.qp = &p_qp;\n";
  codegen_local(g, "colindp", colindp_);
  codegen_local(g, "rowp", rowp_);
  codegen_local(g, "colinda", colinda_);
  codegen_local(g, "rowa", rowa_);
  g << "p.colindp = colindp;\n";
  g << "p.rowp = rowp;\n";
  g << "p.colinda = colinda;\n";
  g << "p.rowa = rowa;\n";
  g << "casadi_clarabel_setup(&p);\n";
}

void ClarabelInterface::set_clarabel_prob() {
  p_.qp = &p_qp_;
  p_.colindp = get_ptr(colindp_);
  p_.rowp = get_ptr(rowp_);
  p_.colinda = get_ptr(colinda_);
  p_.rowa = get_ptr(rowa_);

  casadi_clarabel_setup(&p_);
}

void ClarabelInterface::codegen_init_mem(CodeGenerator& g) const {
  g << "clarabel_init_mem(&" + codegen_mem(g) + ");\n";
  g << "return 0;\n";
}

void ClarabelInterface::codegen_free_mem(CodeGenerator& g) const {
  g << "clarabel_free_mem(&" + codegen_mem(g) + ");\n";
}

int ClarabelInterface::init_mem(void* mem) const {
  if (Conic::init_mem(mem)) return 1;
  if (!mem) return 1;
  auto m = static_cast<ClarabelMemory*>(mem);
  clarabel_init_mem(&m->d);
  m->add_stat("preprocessing");
  m->add_stat("solver");
  m->add_stat("postprocessing");
  return 0;
}

void ClarabelInterface::free_mem(void* mem) const {
  auto m = static_cast<ClarabelMemory*>(mem);
  clarabel_free_mem(&m->d);
  delete m;
}

void ClarabelInterface::set_work(void* mem, const double**& arg, double**& res,
                                 casadi_int*& iw, double*& w) const {
  auto m = static_cast<ClarabelMemory*>(mem);
  Conic::set_work(mem, arg, res, iw, w);
  m->d.prob = &p_;
  m->d.qp = &m->d_qp;
  casadi_clarabel_init(&m->d, &arg, &res, &iw, &w);
}

int ClarabelInterface::solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
  auto m = static_cast<ClarabelMemory*>(mem);
  m->fstats.at("solver").tic();
  casadi_clarabel_solve(&m->d, arg, res, iw, w);
  m->fstats.at("solver").toc();
  return 0;
}

ClarabelInterface::~ClarabelInterface() {
  clear_mem();
}

void ClarabelInterface::codegen_body(CodeGenerator& g) const {
  qp_codegen_body(g);
  g.add_auxiliary(CodeGenerator::AUX_PROJECT);
  g.add_auxiliary(CodeGenerator::AUX_SCAL);
  g.add_auxiliary(CodeGenerator::AUX_SPARSIFY);
  g.add_auxiliary(CodeGenerator::AUX_MAX);
  g.add_auxiliary(CodeGenerator::AUX_SPARSITY);
  g.add_auxiliary(CodeGenerator::AUX_SUM);
  g.add_auxiliary(CodeGenerator::AUX_FILL);
  g.add_auxiliary(CodeGenerator::AUX_CLIP_MIN);
  g.add_auxiliary(CodeGenerator::AUX_CLIP_MAX);
  g.add_auxiliary(CodeGenerator::AUX_DOT);
  g.add_auxiliary(CodeGenerator::AUX_BILIN);
  g.add_include("interfaces/clarabel_c_api.h");

  g.auxiliaries << g.sanitize_source(clarabel_runtime_str, {"casadi_real"});

  g.local("d", "struct casadi_clarabel_data*");
  g.init_local("d", "&" + codegen_mem(g));
  g.local("p", "struct casadi_clarabel_prob");
  set_clarabel_prob(g);

  g << "d->prob = &p;\n";
  g << "d->qp = &d_qp;\n";
  g << "casadi_clarabel_init(d, &arg, &res, &iw, &w);\n";

  g << "casadi_clarabel_solve(d, arg, res, iw, w);\n";
  g << "if (!d_qp.success) {\n";
  if (error_on_fail_) {
    g << "return -1000;\n";
  } else {
    g << "return -1;\n";
  }
  g << "}\n";
  g << "return 0;\n";
}

Dict ClarabelInterface::get_stats(void* mem) const {
  Dict stats = Conic::get_stats(mem);
  auto m = static_cast<ClarabelMemory*>(mem);
  stats["return_status"] = m->d.return_status;
  // (Additional clarabel‐specific stats can be added here.)
  return stats;
}

ClarabelInterface::ClarabelInterface(DeserializingStream& s) : Conic(s) {
  s.version("ClarabelInterface", 1);
  s.unpack("ClarabelInterface::opts", opts_);
  init_dependent();
  set_clarabel_prob();
}

void ClarabelInterface::serialize_body(SerializingStream &s) const {
  Conic::serialize_body(s);
  s.version("ClarabelInterface", 1);
  s.pack("ClarabelInterface::opts", opts_);
}

} // end namespace casadi