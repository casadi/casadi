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


#include "superscs_interface.hpp"
#include "casadi/core/casadi_misc.hpp"

#include <linsys/amatrix.h>
#include <linsys/common.h>
#include <cstring>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_SUPERSCS_EXPORT
  casadi_register_conic_superscs(Conic::Plugin* plugin) {
    plugin->creator = SuperscsInterface::creator;
    plugin->name = "superscs";
    plugin->doc = SuperscsInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &SuperscsInterface::options_;
    plugin->deserialize = &SuperscsInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_SUPERSCS_EXPORT casadi_load_conic_superscs() {
    Conic::registerPlugin(casadi_register_conic_superscs);
  }

  SuperscsInterface::SuperscsInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  SuperscsInterface::~SuperscsInterface() {
    clear_mem();
  }

  const Options SuperscsInterface::options_
  = {{&Conic::options_},
     {{"superscs",
       {OT_DICT,
        "Options to be passed to superscs."}}
     }
  };

  void SuperscsInterface::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    ScsData dummy;
    dummy.stgs = &settings_;
    scs_set_default_settings(&dummy);

    settings_.eps = 1e-6;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="superscs") {
        const Dict& opts = op.second;
        for (auto&& op : opts) {
          if (op.first=="normalize") {
            settings_.normalize = op.second;
          } else if (op.first=="scale") {
            settings_.scale = op.second;
          } else if (op.first=="rho_x") {
            settings_.rho_x = op.second;
          } else if (op.first=="max_time_milliseconds") {
            settings_.max_time_milliseconds = op.second;
          } else if (op.first=="max_iters") {
            settings_.max_iters = op.second;
          } else if (op.first=="previous_max_iters") {
            settings_.previous_max_iters = op.second;
          } else if (op.first=="eps") {
            settings_.eps = op.second;
          } else if (op.first=="alpha") {
            settings_.alpha = op.second;
          } else if (op.first=="cg_rate") {
            settings_.cg_rate = op.second;
          } else if (op.first=="verbose") {
            settings_.verbose = op.second;
          } else if (op.first=="warm_start") {
            settings_.warm_start = op.second;
          } else if (op.first=="do_super_scs") {
            settings_.do_super_scs = op.second;
          } else if (op.first=="k0") {
            settings_.k0 = op.second;
          } else if (op.first=="c_bl") {
            settings_.c_bl = op.second;
          } else if (op.first=="k1") {
            settings_.k1 = op.second;
          } else if (op.first=="k2") {
            settings_.k2 = op.second;
          } else if (op.first=="c1") {
            settings_.c1 = op.second;
          } else if (op.first=="sse") {
            settings_.sse = op.second;
          } else if (op.first=="ls") {
            settings_.ls = op.second;
          } else if (op.first=="beta") {
            settings_.beta = op.second;
          } else if (op.first=="sigma") {
            settings_.sigma = op.second;
          } else if (op.first=="direction") {
            if (op.second=="restarted_broyden") {
              settings_.direction = restarted_broyden;
            } else if (op.second=="anderson_acceleration") {
              settings_.direction = anderson_acceleration;
            } else if (op.second=="fixed_point_residual") {
              settings_.direction = fixed_point_residual;
            } else if (op.second=="full_broyden") {
              settings_.direction = full_broyden;
            } else {
              casadi_error("Unknown argument for direction.");
            }
          } else if (op.first=="thetabar") {
            settings_.thetabar = op.second;
          } else if (op.first=="memory") {
            settings_.memory = op.second;
          } else if (op.first=="tRule") {
            settings_.tRule = op.second;
          } else if (op.first=="broyden_init_scaling") {
            settings_.broyden_init_scaling = op.second;
          } else if (op.first=="do_record_progress") {
            settings_.do_record_progress = op.second;
          } else if (op.first=="do_override_streams") {
            settings_.do_override_streams = op.second;
          } else {
            casadi_error("Not recognised");
          }
        }
      }
    }

    // Create a helper function that trasforms the Hessian into
    // a factorised representation
    //
    // 1/2 x' H x ->  1/2 || F x ||_2^2
    HL_sp_ = H_.ldl(Hp_);
    MX P = DM::eye(nx_)(Slice(), Hp_);

    // Arguments for F Function are the LDL parameters
    MX L = MX::sym("L", HL_sp_);
    MX D = MX::sym("D", Sparsity::diag(nx_));
    F_ = Function("F", {L, D}, {sqrt(2)*mtimes(mtimes(sqrt(D), DM::eye(nx_)+ L), P.T())});

    // Note: Quadratic reformulation
    //
    // min       1/2 x'Hx + G'x
    //  x   st    lb <= Ax <= ub
    //
    //   is transformed to
    //
    //  min     y + G'x
    //  x,y st    lb <= Ax <= ub
    //           || F x || <= 1+y (part 1)
    //                   y <= -1  (part 2)

    // Initialize SDP to SOCP memory
    sdp_to_socp_init(sdp_to_socp_mem_);

    // Canonical form of superscs:
    // min  <c,x>
    //  x,s
    //
    //    subject to  Ax+s = b
    //                s in K^l x K^q1 x K^q2 ...
    //
    //    with:
    //            x in K^l <=> x>=0
    //       (t, x) in K^q <=> ||x|| <= t

    // Note: soc reordering
    //  We will get cones from CasADi
    //  in a form [A1;A2] x  <=> || A1 x || <= A2 x
    //  while we need [A2;A1]
    const SDPToSOCPMem& sm = sdp_to_socp_mem_;
    for (casadi_int i=1;i<sm.r.size();++i) {
      perturb_.push_back(sm.r[i]-1);
      for (casadi_int k=sm.r[i-1];k<sm.r[i]-1;++k) {
        perturb_.push_back(k);
      }
    }

    casadi_int offset = 0;

    // A: simple bounds
    // Index matrix for x <= ubx
    IM B(Sparsity::diag(nx_), range(nx_), false);
    offset += nx_;
    // and lbx <= x
    B = vertcat(B, IM(Sparsity::diag(nx_), range(offset, offset+nx_), false));
    offset += nx_;

    // A: linear constraints
    // Index matrix for A x <= ubx
    IM A(A_, range(offset, offset+A_.nnz()), false);
    offset += A_.nnz();
    // Index matrix for lbx <= A x
    A = vertcat(A, IM(A_, range(offset, offset+A_.nnz()), false));
    offset += A_.nnz();

    // A: quadratic reformulation (part 2)
    // Index matrix for y<=-1
    IM Ae(Sparsity::unit(nx_+1, nx_).T(), offset);
    offset += 1;

    // SOCP helper constraints
    const Sparsity& sp = sm.map_Q.sparsity();
    const casadi_int* colind = sp.colind();
    const casadi_int* row = sp.row();
    const casadi_int* ddata = sm.map_Q.ptr();

    // A: conic constraints
    casadi_int n_c = sp.size2();

    std::vector<casadi_int> c_row, c_col, c_data;

    for (casadi_int j=0; j<n_c; ++j) {
      casadi_int i = perturb_[j];
      for (casadi_int k=colind[i]; k<colind[i+1]-1; ++k) {
        c_col.push_back(row[k]);
        c_row.push_back(j);
        c_data.push_back(offset+ddata[k]);
      }
    }
    offset += Q_.nnz();

    IM C = IM::triplet(c_row, c_col, c_data, n_c, nx_);

    // A: conic constraints stemming from quadratic reformulation (part 1)
    IM F = IM(F_.sparsity_out(0), range(offset, offset+F_.nnz_out(0)));
    offset+= F_.nnz_out(0);

    F = blockcat(IM(1, nx_), offset, F, IM(nx_, 1));
    F = vertcat(F, horzcat(IM(1, nx_), offset+1));

    // Augment an empty column, because of the iontroduction of 'y'
    // in the quadratic reformulation
    A = horzcat(A, IM(A.size1(), 1));
    B = horzcat(B, IM(B.size1(), 1));
    C = horzcat(C, IM(C.size1(), 1));

    At_ = IM::vertcat({B, A, Ae, C, F});

    lookup_ = lookupvector(At_.nonzeros(), 2*nx_ + 2*A_.nnz() + 1 + Q_.nnz()+F_.nnz_out(0)+2);

    alloc_w(At_.nnz(), true);
    alloc_w(At_.size1(), true);
  }

  int SuperscsInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<SuperscsMemory*>(mem);

    m->ldl_d.resize(nx_);
    m->ldl_l.resize(HL_sp_.nnz());
    m->ldl_w.resize(nx_);
    m->F_res.resize(F_.sparsity_out(0).nnz());
    m->g.resize(nx_+1);

    m->data.n = At_.size2();
    m->data.m = At_.size1();
    m->data.stgs = &m->settings;

    m->A.m = At_.size1();
    m->A.n = At_.size2();
    m->at_colind = At_.sparsity().get_colind();
    m->A.p = get_ptr(m->at_colind);
    m->at_row = At_.sparsity().get_row();
    m->A.i = get_ptr(m->at_row);
    m->A.x = nullptr;

    m->data.A = &m->A;
    const SDPToSOCPMem& sm = sdp_to_socp_mem_;

    m->cone.f = 0;
    m->cone.qsize = sm.r.size();
    m->q = diff(sm.r);
    m->q.push_back(nx_+2);
    m->cone.q = get_ptr(m->q);

    m->sol = nullptr;
    m->info = nullptr;

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();
    return 0;
  }

  int SuperscsInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<SuperscsMemory*>(mem);
    const SDPToSOCPMem& sm = sdp_to_socp_mem_;

    // Perform LDL factorization
    bool H_all_zero = true;
    for (casadi_int k=0;k<H_.nnz();++k) {
      H_all_zero = H_all_zero && arg[CONIC_H][k]==0;
    }

    // Circumvent bug in ldl?
    if (H_all_zero) {
      casadi_clear(get_ptr(m->ldl_l), HL_sp_.nnz());
      casadi_clear(get_ptr(m->ldl_d), nx_);
    } else {
      casadi_ldl(H_, arg[CONIC_H], HL_sp_, get_ptr(m->ldl_l), get_ptr(m->ldl_d),
        get_ptr(Hp_), get_ptr(m->ldl_w));
    }

    F_({get_ptr(m->ldl_l), get_ptr(m->ldl_d)}, {get_ptr(m->F_res)});

    double* a_ptr = w; w+= At_.nnz();
    m->A.x = a_ptr;

    const casadi_int* lookup = get_ptr(lookup_);

    // A: simple bounds
    for (casadi_int k=0;k<nx_;++k) {
      a_ptr[*lookup++] = 1;
    }
    for (casadi_int k=0;k<nx_;++k) {
      a_ptr[*lookup++] = -1;
    }
    // A: linear constraints
    for (casadi_int k=0;k<A_.nnz();++k) {
      a_ptr[*lookup++] = arg[CONIC_A][k];
    }
    for (casadi_int k=0;k<A_.nnz();++k) {
      a_ptr[*lookup++] = -arg[CONIC_A][k];
    }

    // A: quadratic reformulation (part 2)
    a_ptr[*lookup++] = -1;

    // A: conic constraints
    for (casadi_int k=0;k<Q_.nnz();++k) {
      casadi_int loc = *lookup++;
      if (loc>=0) {
        a_ptr[loc] = -arg[CONIC_Q][k];
      }
    }

    // A: conic constraints stemming from quadratic reformulation (part 1)
    for (casadi_int k=0;k<F_.nnz_out(0);++k) {
      a_ptr[*lookup++] = -m->F_res[k];
    }
    a_ptr[*lookup++] = -1;
    a_ptr[*lookup++] = 1;

    // b: simple bounds
    double* b_ptr = w;

    m->data.b = b_ptr;
    casadi_clear(b_ptr, At_.size1());

    casadi_axpy(nx_, 1.0, arg[CONIC_UBX], b_ptr);
    b_ptr += nx_;
    casadi_axpy(nx_, -1.0, arg[CONIC_LBX], b_ptr);
    b_ptr += nx_;

    // b: linear constraints
    casadi_copy(arg[CONIC_UBA], na_, b_ptr);
    b_ptr += na_;
    casadi_axpy(na_, -1.0, arg[CONIC_LBA], b_ptr);
    b_ptr += na_;

    // b: quadratic reformulation (part 2)
    *b_ptr = -1;
    b_ptr += 1;

    // b: conic constraints
    for (casadi_int j=0; j<sm.map_Q.size2(); ++j) {
      casadi_int i = perturb_[j];
      // Get bound
      b_ptr[j] = sm.map_P[i]==-1 ? 0 : arg[CONIC_P][sm.map_P[i]];
    }
    b_ptr += sm.map_Q.size2();

    // b: conic constraints stemming from quadratic reformulation (part 1)
    *b_ptr  = 1;
    b_ptr  += 1;
    b_ptr  += nx_;
    *b_ptr  = 1;
    b_ptr  += 1;

    // SuperSCS cannot dealing with infinities,
    // and replacing with arbitrarily large bounds
    // is numerically undesirable for the method

    // -> detect inactive rows
    std::vector<casadi_int> s;
    for (casadi_int k=0;k<At_.size1();++k) {
      if (!isinf(m->data.b[k])) s.push_back(k);
    }
    std::vector<casadi_int> sl = lookupvector(s, At_.size1());
    std::vector<casadi_int> mapping;
    Sparsity Anew = At_.sparsity().sub(s, range(At_.size2()), mapping);

    // Trim A
    casadi_copy(Anew.colind(), nx_+2, m->A.p);
    casadi_copy(Anew.row(), Anew.nnz(), m->A.i);
    m->A.m = Anew.size1();
    m->data.m = Anew.size1();
    for (casadi_int k=0;k<Anew.nnz();++k) {
      m->A.x[k] = m->A.x[mapping[k]];
    }

    // Trim b
    for (casadi_int k=0;k<At_.size1();++k) {
      casadi_int e = sl[k];
      if (e>=0) m->data.b[e] = m->data.b[k];
    }

    // How many linear constraints we have
    // depends on the trimming procedure
    m->cone.l = 2*nx_+2*na_+1-(At_.size1()-Anew.size1());

    // Set objective
    casadi_copy(arg[CONIC_G], nx_, get_ptr(m->g));
    // Part related to quadratic reformulation: y
    m->g[nx_] = 1;

    m->data.c = get_ptr(m->g);

    // No other cones
    m->cone.ssize = 0;
    m->cone.ed = 0;
    m->cone.ep = 0;
    m->cone.psize = 0;
    m->cone.p = nullptr;
    m->cone.s = nullptr;

    if (m->sol) scs_free_sol(m->sol);
    if (m->info) scs_free_info(m->info);

    m->sol = scs_init_sol();
    m->info = scs_init_info();

    // Copy settings
    std::memcpy(m->data.stgs, &settings_, sizeof(ScsSettings));

    casadi_int status = scs(&m->data, &m->cone, m->sol, m->info);

    m->success = SCS_SOLVED==status;

    casadi_copy(m->sol->x, nx_, res[CONIC_X]);
    if (res[CONIC_COST])
      *res[CONIC_COST] = casadi_dot(nx_, m->sol->x, get_ptr(m->g));
    if (res[CONIC_COST])
      *res[CONIC_COST] += 0.5*casadi_bilin(arg[CONIC_H], H_, m->sol->x, m->sol->x);

    return 0;
  }

  Dict SuperscsInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    //auto m = static_cast<SuperscsMemory*>(mem);
    //stats["return_status"] = return_status_string(m->return_status);
    return stats;
  }

  SuperscsMemory::SuperscsMemory() {
  }

  SuperscsMemory::~SuperscsMemory() {
    if (sol) scs_free_sol(sol);
    if (info) scs_free_info(info);
  }

  SuperscsInterface::SuperscsInterface(DeserializingStream& s) : Conic(s) {
    s.version("SuperscsInterface", 1);
    ScsData dummy;
    dummy.stgs = &settings_;
    scs_set_default_settings(&dummy);
    s.unpack("SuperscsInterface::settings::normalize", settings_.normalize);
    s.unpack("SuperscsInterface::settings::scale", settings_.scale);
    s.unpack("SuperscsInterface::settings::rho_x", settings_.rho_x);
    s.unpack("SuperscsInterface::settings::max_time_milliseconds", settings_.max_time_milliseconds);
    s.unpack("SuperscsInterface::settings::max_iters", settings_.max_iters);
    s.unpack("SuperscsInterface::settings::previous_max_iters", settings_.previous_max_iters);
    s.unpack("SuperscsInterface::settings::eps", settings_.eps);
    s.unpack("SuperscsInterface::settings::alpha", settings_.alpha);
    s.unpack("SuperscsInterface::settings::cg_rate", settings_.cg_rate);
    s.unpack("SuperscsInterface::settings::verbose", settings_.verbose);
    s.unpack("SuperscsInterface::settings::warm_start", settings_.warm_start);
    s.unpack("SuperscsInterface::settings::do_super_scs", settings_.do_super_scs);
    s.unpack("SuperscsInterface::settings::k0", settings_.k0);
    s.unpack("SuperscsInterface::settings::c_bl", settings_.c_bl);
    s.unpack("SuperscsInterface::settings::k1", settings_.k1);
    s.unpack("SuperscsInterface::settings::k2", settings_.k2);
    s.unpack("SuperscsInterface::settings::c1", settings_.c1);
    s.unpack("SuperscsInterface::settings::sse", settings_.sse);
    s.unpack("SuperscsInterface::settings::ls", settings_.ls);
    s.unpack("SuperscsInterface::settings::beta", settings_.beta);
    s.unpack("SuperscsInterface::settings::sigma", settings_.sigma);
    casadi_int dir;
    s.unpack("SuperscsInterface::settings::direction", dir);
    settings_.direction = static_cast<direction_enum>(dir);
    s.unpack("SuperscsInterface::settings::thetabar", settings_.thetabar);
    s.unpack("SuperscsInterface::settings::memory", settings_.memory);
    s.unpack("SuperscsInterface::settings::tRule", settings_.tRule);
    s.unpack("SuperscsInterface::settings::broyden_init_scaling", settings_.broyden_init_scaling);
    s.unpack("SuperscsInterface::settings::do_record_progress", settings_.do_record_progress);
    s.unpack("SuperscsInterface::settings::do_override_streams", settings_.do_override_streams);
    s.unpack("SuperscsInterface::Hp", Hp_);
    s.unpack("SuperscsInterface::HL_sp", HL_sp_);
    s.unpack("SuperscsInterface::f", F_);
    s.unpack("SuperscsInterface::At", At_);
    s.unpack("SuperscsInterface::lookup", lookup_);
    s.unpack("SuperscsInterface::perturb", perturb_);
    s.unpack("SuperscsInterface::opts", opts_);
    Conic::deserialize(s, sdp_to_socp_mem_);
  }

  void SuperscsInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);
    s.version("SuperscsInterface", 1);
    s.pack("SuperscsInterface::settings::normalize", settings_.normalize);
    s.pack("SuperscsInterface::settings::scale", settings_.scale);
    s.pack("SuperscsInterface::settings::rho_x", settings_.rho_x);
    s.pack("SuperscsInterface::settings::max_time_milliseconds", settings_.max_time_milliseconds);
    s.pack("SuperscsInterface::settings::max_iters", settings_.max_iters);
    s.pack("SuperscsInterface::settings::previous_max_iters", settings_.previous_max_iters);
    s.pack("SuperscsInterface::settings::eps", settings_.eps);
    s.pack("SuperscsInterface::settings::alpha", settings_.alpha);
    s.pack("SuperscsInterface::settings::cg_rate", settings_.cg_rate);
    s.pack("SuperscsInterface::settings::verbose", settings_.verbose);
    s.pack("SuperscsInterface::settings::warm_start", settings_.warm_start);
    s.pack("SuperscsInterface::settings::do_super_scs", settings_.do_super_scs);
    s.pack("SuperscsInterface::settings::k0", settings_.k0);
    s.pack("SuperscsInterface::settings::c_bl", settings_.c_bl);
    s.pack("SuperscsInterface::settings::k1", settings_.k1);
    s.pack("SuperscsInterface::settings::k2", settings_.k2);
    s.pack("SuperscsInterface::settings::c1", settings_.c1);
    s.pack("SuperscsInterface::settings::sse", settings_.sse);
    s.pack("SuperscsInterface::settings::ls", settings_.ls);
    s.pack("SuperscsInterface::settings::beta", settings_.beta);
    s.pack("SuperscsInterface::settings::sigma", settings_.sigma);
    s.pack("SuperscsInterface::settings::direction", static_cast<casadi_int>(settings_.direction));
    s.pack("SuperscsInterface::settings::thetabar", settings_.thetabar);
    s.pack("SuperscsInterface::settings::memory", settings_.memory);
    s.pack("SuperscsInterface::settings::tRule", settings_.tRule);
    s.pack("SuperscsInterface::settings::broyden_init_scaling", settings_.broyden_init_scaling);
    s.pack("SuperscsInterface::settings::do_record_progress", settings_.do_record_progress);
    s.pack("SuperscsInterface::settings::do_override_streams", settings_.do_override_streams);
    s.pack("SuperscsInterface::Hp", Hp_);
    s.pack("SuperscsInterface::HL_sp", HL_sp_);
    s.pack("SuperscsInterface::f", F_);
    s.pack("SuperscsInterface::At", At_);
    s.pack("SuperscsInterface::lookup", lookup_);
    s.pack("SuperscsInterface::perturb", perturb_);
    s.pack("SuperscsInterface::opts", opts_);
    Conic::serialize(s, sdp_to_socp_mem_);
  }

} // namespace casadi
