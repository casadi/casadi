/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef SOLVE_IMPL_HPP
#define SOLVE_IMPL_HPP

#include "solve.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "../fx/fx_internal.hpp"

using namespace std;

namespace CasADi{

  template<bool Tr>
  Solve<Tr>::Solve(const MX& r, const MX& A){
    casadi_assert_message(r.size2() == A.size1(),"Solve::Solve: dimension mismatch.");
    setDependencies(r,A);
    setSparsity(r.sparsity());
  }

  template<bool Tr>
  void Solve<Tr>::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "(";
    } else if(part==1){
      stream << "/";
    } else {
      if(Tr) stream << "'";
      stream << ")";
    }
  }

  template<bool Tr>
  void Solve<Tr>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<bool Tr>
  void Solve<Tr>::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<bool Tr>
  template<typename T, typename MatV, typename MatVV>
  void Solve<Tr>::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    if(input[0]!=output[0]){
      copy(input[0]->begin(),input[0]->end(),output[0]->begin());
    }
    casadi_error("not implemented");

    // Forward sensitivities
    for(int d=0; d<nfwd; ++d){
      if(fwdSeed[d][0]!=fwdSens[d][0]){
        copy(fwdSeed[d][0]->begin(),fwdSeed[d][0]->end(),fwdSens[d][0]->begin());
      }
      casadi_error("not implemented");
    }

    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      casadi_error("not implemented");
      if(adjSeed[d][0]!=adjSens[d][0]){
        transform(adjSeed[d][0]->begin(),adjSeed[d][0]->end(),adjSens[d][0]->begin(),adjSens[d][0]->begin(),std::plus<T>());
        adjSeed[d][0]->setZero();
      }
    }
  }

  template<bool Tr>
  void Solve<Tr>::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    const MX& r = *input[0];
    const MX& A = *input[1];
    MX& x = *output[0];
    if(!output_given){
      x = A->getSolve(r,Tr);
    }

    // Forward sensitivities, collect the right hand sides
    int nfwd = fwdSens.size();
    vector<MX> rhs;
    vector<int> row_offset(1,0);
    for(int d=0; d<nfwd; ++d){
      const MX& r_dot = *fwdSeed[d][0];
      const MX& A_dot = *fwdSeed[d][1];
      rhs.push_back(r_dot - mul(x,Tr ? trans(A_dot) : A_dot));
      row_offset.push_back(row_offset.back()+r_dot.size1());
    }

    // Solve for all directions at once
    MX lhs = A->getSolve(vertcat(rhs),Tr);
    
    // Slit up the sensitivities
    for(int d=0; d<nfwd; ++d){
      MX& x_dot = *fwdSens[d][0];
      x_dot = lhs(Slice(row_offset[d],row_offset[d+1]),Slice());
      x_dot = MX();
    }
  
    // Adjoint sensitivities, collect right hand sides
    int nadj = adjSeed.size();
    rhs.resize(0);
    row_offset.resize(1);
    for(int d=0; d<nadj; ++d){
      const MX& x_dot = *adjSeed[d][0];      
      rhs.push_back(x_dot);
      row_offset.push_back(row_offset.back()+x_dot.size1());
    }

    // Solve for all directions at once
    lhs = A->getSolve(vertcat(rhs),!Tr);
    
    // Split up the sensitivities
    for(int d=0; d<nadj; ++d){
      MX r_bar = lhs(Slice(row_offset[d],row_offset[d+1]),Slice());      
      *adjSens[d][0] += r_bar;
      *adjSens[d][1] -= mul(trans(x),r_bar);
    }
  }
  
  //  template<bool Tr>
  //  void Solve<Tr>::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
  //    casadi_error("not implemented");
  //  }

  template<bool Tr>
  void Solve<Tr>::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    // Check if inplace
    bool inplace = arg.at(0).compare(res.front())==0;

    // Copy first argument if not inplace
    if(!inplace){      
      stream << "  for(i=0; i<" << this->size() << "; ++i) " << res.front() << "[i]=" << arg.at(0) << "[i];" << endl;
    }

    casadi_error("not implemented");
  }

} // namespace CasADi

#endif // SOLVE_IMPL_HPP

