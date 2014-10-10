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


#include "sdp_solver_internal.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/sparsity_tools.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"

INPUTSCHEME(SDPInput)
OUTPUTSCHEME(SDPOutput)

using namespace std;
namespace casadi {

// Constructor
SdpSolverInternal::SdpSolverInternal(const std::vector<Sparsity> &st) : st_(st) {
  addOption("calc_p", OT_BOOLEAN, true,
            "Indicate if the P-part of primal solution should be allocated and calculated. "
            "You may want to avoid calculating this variable for problems with n large, "
            "as is always dense (m x m).");
  addOption("calc_dual", OT_BOOLEAN, true, "Indicate if dual should be allocated and calculated. "
            "You may want to avoid calculating this variable for problems with n large, "
            "as is always dense (m x m).");
  addOption("print_problem", OT_BOOLEAN, false, "Print out problem statement for debugging.");

  casadi_assert_message(st_.size()==SDP_STRUCT_NUM, "Problem structure mismatch");

  const Sparsity& A = st_[SDP_STRUCT_A];
  const Sparsity& G = st_[SDP_STRUCT_G];
  const Sparsity& F = st_[SDP_STRUCT_F];

  casadi_assert_message(G==G.transpose(), "SdpSolverInternal: Supplied G sparsity must "
                        "symmetric but got " << G.dimString());

  m_ = G.size1();

  nc_ = A.size1();
  n_ = A.size2();

  casadi_assert_message(F.size1()==m_, "SdpSolverInternal: Supplied F sparsity: number of rows ("
                        << F.size1() <<  ")  must match m (" << m_ << ")");

  casadi_assert_message(F.size2()%n_==0, "SdpSolverInternal: Supplied F sparsity: "
                        "number of columns (" << F.size2()
                        <<  ")  must be an integer multiple of n (" << n_
                        << "), but got remainder " << F.size2()%n_);

  // Input arguments
  setNumInputs(SDP_SOLVER_NUM_IN);
  input(SDP_SOLVER_G) = DMatrix(G, 0);
  input(SDP_SOLVER_F) = DMatrix(F, 0);
  input(SDP_SOLVER_A) = DMatrix(A, 0);
  input(SDP_SOLVER_C) = DMatrix::zeros(n_);
  input(SDP_SOLVER_LBX) = -DMatrix::inf(n_);
  input(SDP_SOLVER_UBX) = DMatrix::inf(n_);
  input(SDP_SOLVER_LBA) = -DMatrix::inf(nc_);
  input(SDP_SOLVER_UBA) = DMatrix::inf(nc_);

  for (int i=0;i<n_;i++) {
    Sparsity s = input(SDP_SOLVER_F)(ALL, Slice(i*m_, (i+1)*m_)).sparsity();
    casadi_assert_message(s==s.transpose(),
                          "SdpSolverInternal: Each supplied Fi must be symmetric. "
                          "But got " << s.dimString() <<  " for i = " << i << ".");
  }

  input_.scheme = SCHEME_SDPInput;
  output_.scheme = SCHEME_SDPOutput;

}

void SdpSolverInternal::init() {
  // Call the init method of the base class
  FunctionInternal::init();

  calc_p_ = getOption("calc_p");
  calc_dual_ = getOption("calc_dual");
  print_problem_ = getOption("print_problem");

  // Find aggregate sparsity pattern
  Sparsity aggregate = input(SDP_SOLVER_G).sparsity();
  for (int i=0;i<n_;++i) {
    aggregate = aggregate + input(SDP_SOLVER_F)(ALL, Slice(i*m_, (i+1)*m_)).sparsity();
  }

  // Detect block diagonal structure in this sparsity pattern
  std::vector<int> p;
  std::vector<int> r;
  nb_ = aggregate.stronglyConnectedComponents(p, r);
  block_boundaries_.resize(nb_+1);
  std::copy(r.begin(), r.begin()+nb_+1, block_boundaries_.begin());

  block_sizes_.resize(nb_);
  for (int i=0;i<nb_;++i) {
    block_sizes_[i]=r[i+1]-r[i];
  }

  // Make a mapping function from dense blocks to inversely-permuted block diagonal P
  std::vector< SX > full_blocks;
  for (int i=0;i<nb_;++i) {
    full_blocks.push_back(SX::sym("block", block_sizes_[i], block_sizes_[i]));
  }

  Pmapper_ = SXFunction(full_blocks, blkdiag(full_blocks)(lookupvector(p, p.size()),
                                                         lookupvector(p, p.size())));
  Pmapper_.init();

  if (nb_>0) {
    // Make a mapping function from (G, F) -> (G[p, p]_j, F_i[p, p]j)
    SX G = SX::sym("G", input(SDP_SOLVER_G).sparsity());
    SX F = SX::sym("F", input(SDP_SOLVER_F).sparsity());

    std::vector<SX> in;
    in.push_back(G);
    in.push_back(F);
    std::vector<SX> out((n_+1)*nb_);
    for (int j=0;j<nb_;++j) {
      out[j] = G(p, p)(Slice(r[j], r[j+1]), Slice(r[j], r[j+1]));
    }
    for (int i=0;i<n_;++i) {
      SX Fi = F(ALL, Slice(i*m_, (i+1)*m_))(p, p);
      for (int j=0;j<nb_;++j) {
        out[(i+1)*nb_+j] = Fi(Slice(r[j], r[j+1]), Slice(r[j], r[j+1]));
      }
    }
    mapping_ = SXFunction(in, out);
    mapping_.init();
  }

  // Output arguments
  setNumOutputs(SDP_SOLVER_NUM_OUT);
  output(SDP_SOLVER_X) = DMatrix::zeros(n_, 1);
  output(SDP_SOLVER_P) = calc_p_? DMatrix(Pmapper_.output().sparsity(), 0) : DMatrix();
  output(SDP_SOLVER_DUAL) = calc_dual_? DMatrix(Pmapper_.output().sparsity(), 0) : DMatrix();
  output(SDP_SOLVER_COST) = 0.0;
  output(SDP_SOLVER_DUAL_COST) = 0.0;
  output(SDP_SOLVER_LAM_X) = DMatrix::zeros(n_, 1);
  output(SDP_SOLVER_LAM_A) = DMatrix::zeros(nc_, 1);

}

SdpSolverInternal::~SdpSolverInternal() {
}

void SdpSolverInternal::printProblem(std::ostream &stream) const {
  stream << "SDP Problem statement -- start" << std::endl;

  stream << "f: "<< std::endl;  input(SDP_SOLVER_F).printDense(stream);
  stream << "c: "<< std::endl;  input(SDP_SOLVER_C).printDense(stream);
  stream << "g: "<< std::endl;  input(SDP_SOLVER_G).printDense(stream);
  stream << "a: "<< std::endl;  input(SDP_SOLVER_A).printDense(stream);
  stream << "lba: " << input(SDP_SOLVER_LBA) << std::endl;
  stream << "uba: " << input(SDP_SOLVER_UBA) << std::endl;
  stream << "lbx: " << input(SDP_SOLVER_LBX) << std::endl;
  stream << "ubx: " << input(SDP_SOLVER_UBX) << std::endl;

  stream << "SDP Problem statement -- end" << std::endl;
}

void SdpSolverInternal::evaluate() {
  throw CasadiException("SdpSolverInternal::evaluate: Not implemented");
}

void SdpSolverInternal::solve() {
  throw CasadiException("SdpSolverInternal::solve: Not implemented");
}

void SdpSolverInternal::checkInputs() const {
  for (int i=0;i<input(SDP_SOLVER_LBX).size();++i) {
    casadi_assert_message(input(SDP_SOLVER_LBX).at(i)<=input(SDP_SOLVER_UBX).at(i),
                          "LBX[i] <= UBX[i] was violated for i=" << i
                          << ". Got LBX[i]=" << input(SDP_SOLVER_LBX).at(i)
                          << " and UBX[i]=" << input(SDP_SOLVER_UBX).at(i));
  }
  for (int i=0;i<input(SDP_SOLVER_LBA).size();++i) {
    casadi_assert_message(input(SDP_SOLVER_LBA).at(i)<=input(SDP_SOLVER_UBA).at(i),
                          "LBA[i] <= UBA[i] was violated for i=" << i
                          << ". Got LBA[i]=" << input(SDP_SOLVER_LBA).at(i)
                          << " and UBA[i]=" << input(SDP_SOLVER_UBA).at(i));
  }
}

  std::map<std::string, SdpSolverInternal::Plugin> SdpSolverInternal::solvers_;

  const std::string SdpSolverInternal::infix_ = "sdpsolver";

} // namespace casadi
