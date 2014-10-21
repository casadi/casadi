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


#include "sdqp_to_sdp.hpp"

#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"

#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_SDQPSOLVER_SDP_EXPORT
  casadi_register_sdqpsolver_sdp(SdqpSolverInternal::Plugin* plugin) {
    plugin->creator = SdqpToSdp::creator;
    plugin->name = "sdp";
    plugin->doc = SdqpToSdp::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_SDQPSOLVER_SDP_EXPORT casadi_load_sdqpsolver_sdp() {
    SdqpSolverInternal::registerPlugin(casadi_register_sdqpsolver_sdp);
  }

  SdqpToSdp::SdqpToSdp(const std::vector<Sparsity> &st) : SdqpSolverInternal(st) {
    Adaptor::addOptions();
  }

  SdqpToSdp::~SdqpToSdp() {
  }

  void SdqpToSdp::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    SdqpSolverInternal::deepCopyMembers(already_copied);
    solver_ = deepcopy(solver_, already_copied);
    cholesky_ = deepcopy(cholesky_, already_copied);
    mapping_ = deepcopy(mapping_, already_copied);
  }

  void SdqpToSdp::init() {
    // Initialize the base classes
    SdqpSolverInternal::init();

    cholesky_ = LinearSolver("csparsecholesky", st_[SDQP_STRUCT_H]);
    cholesky_.init();

    MX g_socp = MX::sym("x", cholesky_.getFactorizationSparsity(true));
    MX h_socp = MX::sym("h", n_);

    MX f_socp = sqrt(inner_prod(h_socp, h_socp));
    MX en_socp = 0.5/f_socp;

    MX f_sdqp = MX::sym("f", input(SDQP_SOLVER_F).sparsity());
    MX g_sdqp = MX::sym("g", input(SDQP_SOLVER_G).sparsity());

    std::vector<MX> fi(n_+1);
    MX znp = MX::sparse(n_+1, n_+1);
    for (int k=0;k<n_;++k) {
      MX gk = vertcat(g_socp(ALL, k), DMatrix::sparse(1, 1));
      MX fk = -blockcat(znp, gk, gk.T(), DMatrix::sparse(1, 1));
      // TODO(Joel): replace with ALL
      fi.push_back(blkdiag(f_sdqp(ALL, Slice(f_sdqp.size1()*k, f_sdqp.size1()*(k+1))), fk));
    }
    MX fin = en_socp*DMatrix::eye(n_+2);
    fin(n_, n_+1) = en_socp;
    fin(n_+1, n_) = en_socp;

    fi.push_back(blkdiag(DMatrix::sparse(f_sdqp.size1(), f_sdqp.size1()), -fin));

    MX h0 = vertcat(h_socp, DMatrix::sparse(1, 1));
    MX g = blockcat(f_socp*DMatrix::eye(n_+1), h0, h0.T(), f_socp);

    g = blkdiag(g_sdqp, g);

    IOScheme mappingIn("g_socp", "h_socp", "f_sdqp", "g_sdqp");
    IOScheme mappingOut("f", "g");

    mapping_ = MXFunction(mappingIn("g_socp", g_socp, "h_socp", h_socp,
                                    "f_sdqp", f_sdqp, "g_sdqp", g_sdqp),
                          mappingOut("f", horzcat(fi), "g", g));
    mapping_.init();

    // Create an SdpSolver instance
    solver_ = SdpSolver(getOption(solvername()),
                        sdpStruct("g", mapping_.output("g").sparsity(),
                                  "f", mapping_.output("f").sparsity(),
                                  "a", horzcat(input(SDQP_SOLVER_A).sparsity(),
                                               Sparsity::sparse(nc_, 1))));
    if (hasSetOption(optionsname())) solver_.setOption(getOption(optionsname()));
    solver_.init();

    solver_.input(SDP_SOLVER_C).at(n_)=1;

    // Output arguments
    setNumOutputs(SDQP_SOLVER_NUM_OUT);
    output(SDQP_SOLVER_X) = DMatrix::zeros(n_, 1);

    std::vector<int> r = range(input(SDQP_SOLVER_G).size1());
    output(SDQP_SOLVER_P) = solver_.output(SDP_SOLVER_P).isEmpty() ? DMatrix() :
        solver_.output(SDP_SOLVER_P)(r, r);
    output(SDQP_SOLVER_DUAL) = solver_.output(SDP_SOLVER_DUAL).isEmpty() ? DMatrix() :
        solver_.output(SDP_SOLVER_DUAL)(r, r);
    output(SDQP_SOLVER_COST) = 0.0;
    output(SDQP_SOLVER_DUAL_COST) = 0.0;
    output(SDQP_SOLVER_LAM_X) = DMatrix::zeros(n_, 1);
    output(SDQP_SOLVER_LAM_A) = DMatrix::zeros(nc_, 1);
  }

  void SdqpToSdp::evaluate() {
    cholesky_.setInput(input(SDQP_SOLVER_H));
    for (int k=0;k<cholesky_.input(0).size();++k) {
      cholesky_.input().at(k)*=0.5;
    }
    cholesky_.prepare();
    mapping_.setInput(cholesky_.getFactorization(true), "g_socp");
    std::copy(input(SDQP_SOLVER_C).begin(),
              input(SDQP_SOLVER_C).end(),
              mapping_.input("h_socp").begin());
    cholesky_.solveL(&mapping_.input("h_socp").data().front(), 1, true);
    for (int k=0;k<mapping_.input("h_socp").size();++k) {
      mapping_.input("h_socp").at(k)*= 0.5;
    }

    mapping_.setInput(input(SDQP_SOLVER_F), "f_sdqp");
    mapping_.setInput(input(SDQP_SOLVER_G), "g_sdqp");

    mapping_.evaluate();

    std::copy(input(SDQP_SOLVER_A).begin(),
              input(SDQP_SOLVER_A).end(),
              solver_.input(SDP_SOLVER_A).begin());
    solver_.setInput(mapping_.output("f"), SDP_SOLVER_F);
    solver_.setInput(mapping_.output("g"), SDP_SOLVER_G);
    std::copy(input(SDQP_SOLVER_LBA).begin(),
              input(SDQP_SOLVER_LBA).end(),
              solver_.input(SDP_SOLVER_LBA).begin());
    std::copy(input(SDQP_SOLVER_UBA).begin(),
              input(SDQP_SOLVER_UBA).end(),
              solver_.input(SDP_SOLVER_UBA).begin());
    std::copy(input(SDQP_SOLVER_LBX).begin(),
              input(SDQP_SOLVER_LBX).end(),
              solver_.input(SDP_SOLVER_LBX).begin());
    std::copy(input(SDQP_SOLVER_UBX).begin(),
              input(SDQP_SOLVER_UBX).end(),
              solver_.input(SDP_SOLVER_UBX).begin());

    solver_.evaluate();

    // Pass the stats
    stats_["sdp_solver_stats"] = solver_.getStats();

    std::copy(solver_.output(SDP_SOLVER_X).begin(),
              solver_.output(SDP_SOLVER_X).begin()+n_,
              output(SDQP_SOLVER_X).begin());
    setOutput(solver_.output(SDP_SOLVER_COST),
              SDQP_SOLVER_COST);
    setOutput(solver_.output(SDP_SOLVER_DUAL_COST),
              SDQP_SOLVER_DUAL_COST);
    if (!output(SDQP_SOLVER_DUAL).isEmpty())
        std::copy(solver_.output(SDP_SOLVER_DUAL).begin(),
                  solver_.output(SDP_SOLVER_DUAL).begin()+output(SDQP_SOLVER_DUAL).size(),
                  output(SDQP_SOLVER_DUAL).begin());
    if (!output(SDQP_SOLVER_P).isEmpty())
        std::copy(solver_.output(SDP_SOLVER_P).begin(),
                  solver_.output(SDP_SOLVER_P).begin()+output(SDQP_SOLVER_P).size(),
                  output(SDQP_SOLVER_P).begin());
    std::copy(solver_.output(SDP_SOLVER_X).begin(),
              solver_.output(SDP_SOLVER_X).begin()+n_,
              output(SDQP_SOLVER_X).begin());
    std::copy(solver_.output(SDP_SOLVER_LAM_A).begin(),
              solver_.output(SDP_SOLVER_LAM_A).end(),
              output(SDQP_SOLVER_LAM_A).begin());
    std::copy(solver_.output(SDP_SOLVER_LAM_X).begin(),
              solver_.output(SDP_SOLVER_LAM_X).begin()+n_,
              output(SDQP_SOLVER_LAM_X).begin());
  }

} // namespace casadi
