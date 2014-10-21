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


#include "socp_to_sdp.hpp"

#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"
#include "casadi/core/function/mx_function.hpp"
#include "casadi/core/mx/mx_tools.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_SOCPSOLVER_SDP_EXPORT
  casadi_register_socpsolver_sdp(SocpSolverInternal::Plugin* plugin) {
    plugin->creator = SocpToSdp::creator;
    plugin->name = "sdp";
    plugin->doc = SocpToSdp::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_SOCPSOLVER_SDP_EXPORT casadi_load_socpsolver_sdp() {
    SocpSolverInternal::registerPlugin(casadi_register_socpsolver_sdp);
  }

  SocpToSdp* SocpToSdp::clone() const {
    // Return a deep copy
    SocpToSdp* node = new SocpToSdp(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  SocpToSdp::SocpToSdp(const std::vector<Sparsity> &st) : SocpSolverInternal(st) {
    Adaptor::addOptions();
  }

  SocpToSdp::~SocpToSdp() {
  }

  void SocpToSdp::evaluate() {
    if (print_problem_) printProblem();

    mapping_.setInput(input(SOCP_SOLVER_G), 0);
    mapping_.setInput(input(SOCP_SOLVER_H), 1);
    mapping_.setInput(input(SOCP_SOLVER_E), 2);
    mapping_.setInput(input(SOCP_SOLVER_F), 3);

    mapping_.evaluate();

    solver_.setInput(mapping_.output(0), SDP_SOLVER_F);
    solver_.setInput(mapping_.output(1), SDP_SOLVER_G);
    solver_.setInput(input(SOCP_SOLVER_A), SDP_SOLVER_A);
    solver_.setInput(input(SOCP_SOLVER_C), SDP_SOLVER_C);
    solver_.setInput(input(SOCP_SOLVER_LBX), SDP_SOLVER_LBX);
    solver_.setInput(input(SOCP_SOLVER_UBX), SDP_SOLVER_UBX);
    solver_.setInput(input(SOCP_SOLVER_LBA), SDP_SOLVER_LBA);
    solver_.setInput(input(SOCP_SOLVER_UBA), SDP_SOLVER_UBA);

    solver_.evaluate();

    // Pass the stats
    stats_["sdp_solver_stats"] = solver_.getStats();

    setOutput(solver_.output(SDP_SOLVER_X), SOCP_SOLVER_X);
    setOutput(solver_.output(SDP_SOLVER_COST), SOCP_SOLVER_COST);
    setOutput(solver_.output(SDP_SOLVER_LAM_X), SOCP_SOLVER_LAM_X);
    setOutput(solver_.output(SDP_SOLVER_LAM_A), SOCP_SOLVER_LAM_A);
  }

  void SocpToSdp::init() {
    // Initialize the base classes
    SocpSolverInternal::init();

    /*
     *   || Gi' x + hi ||_2 <=  ei'x + fi
     *
     *        <=>
     *
     *  | (ei' x + fi) I   Gi' x + hi   |   >= 0
     *  | (Gi' x + hi)'     ei' x + fi  |
     *
     *        <=>
     *
     *  | ei[k]I     Gi(:, k)'  |
     *  | Gi(:, k)     ei[k]  |
     *
     *      | k  (n)      \   i  (for each cone)
     *      v              \
     */


    MX G = MX::sym("G", input(SOCP_SOLVER_G).sparsity());
    MX H = MX::sym("H", input(SOCP_SOLVER_H).sparsity());
    MX E = MX::sym("E", input(SOCP_SOLVER_E).sparsity());
    MX F = MX::sym("F", input(SOCP_SOLVER_F).sparsity());


    int i_start;

    // The blocks that will make up Fi of SDP
    std::vector<MX> Fi;

    for (int k=0;k<n_;++k) {

      // The blocks that will end up on the diagonal of Fi of SDP
      std::vector<MX> Fi_d;

      i_start = 0;
      // Loop over all SOCP constraints
      for (int i=0;i<ni_.size();++i) {
        MX Gik = G(Slice(k), Slice(i_start, i_start+ni_[i])).T();
        MX Eik = E[n_*i+k];
        Fi_d.push_back(blockcat(Eik*MX::eye(ni_[i]), Gik, Gik.T(), Eik));
        i_start += ni_[i];
      }
      Fi.push_back(blkdiag(Fi_d));
    }
    // The blocks that will end up on the diagonal of G of SDP
    std::vector<MX> G_d;

    i_start = 0;
    // Loop over all SOCP constraints
    for (int i=0;i<ni_.size();++i) {
      MX Fi  = F[i];
      MX Hi  = H[range(i_start, i_start+ni_[i])];
      G_d.push_back(blockcat(Fi*MX::eye(ni_[i]), Hi, Hi.T(), Fi));
      i_start += ni_[i];
    }

    std::vector<MX> out;
    out.push_back(-horzcat(Fi));
    out.push_back(blkdiag(G_d));

    std::vector<MX> syms;
    syms.push_back(G);
    syms.push_back(H);
    syms.push_back(E);
    syms.push_back(F);

    mapping_ = MXFunction(syms, out);
    mapping_.init();

    log("SocpToSdp::init", "Created mapping function");

    // Create an SdpSolver instance
    solver_ = SdpSolver(getOption(solvername()),
                        sdpStruct("a", input(SOCP_SOLVER_A).sparsity(),
                                  "f", mapping_.output(0).sparsity(),
                                  "g", mapping_.output(1).sparsity()));
    solver_.setSOCPOptions();
    if (hasSetOption(optionsname())) solver_.setOption(getOption(optionsname()));
    solver_.init();

    log("SocpToSdp::init", "Initialized SDP solver");
  }

} // namespace casadi

