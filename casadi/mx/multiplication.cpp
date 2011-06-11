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

#include "multiplication.hpp"
#include "../matrix/matrix_tools.hpp"
#include <vector>

using namespace std;

namespace CasADi{

Multiplication::Multiplication(const MX& x, const MX& y){
  casadi_assert_message(x.size2() == y.size1(),"Multiplication::Multiplication: dimension mismatch");
  setDependencies(x,y);

  // Check if the arguments are dense
  x_dense_ = x.size()==x.numel();
  y_dense_ = y.size()==y.numel();
  
  if(x_dense_ && y_dense_){
    // Dense matrix
    setSparsity(CRSSparsity(x.size1(),y.size2(),true));
    
  } else {
    // Find the mapping corresponding to the transpose of y (no need to form the transpose explicitly)
    CRSSparsity y_trans_sparsity = y->sparsity().transpose(y_trans_map_);

    // Create the sparsity pattern for the matrix-matrix product
    CRSSparsity spres = x->sparsity().patternProduct(y_trans_sparsity, prod_map_);
  
    // Save sparsity
    setSparsity(spres);
  }
}

Multiplication* Multiplication::clone() const{
  return new Multiplication(*this);
}

void Multiplication::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "prod(" << args.at(0) << "," << args.at(1) << ")";
}

void Multiplication::evaluateDenseDense(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  // Get dimensions
  int nx1 = dep(0).size1();
  int nx2 = dep(0).size2();
  //int ny1 = dep(1).size1();
  int ny2 = dep(1).size2();
  //int nz1 = nx1;
  int nz2 = ny2;
    
  for(int i=0; i<nx1; ++i){
    for(int j=0; j<ny2; ++j){
      // Add scalar product
      double sum = 0;
      for(int k=0; k<nx2; ++k){
        sum += input[0][k+i*nx2] * input[1][j+k*ny2];
      }
      output[j+i*nz2] = sum;
    }
  }
  
  // Forward sensitivities: dot(Z) = dot(X)*Y + X*dot(Y)
  for(int d=0; d<nfwd; ++d){
    for(int i=0; i<nx1; ++i){
      for(int j=0; j<ny2; ++j){
        // Add scalar product
        double sum = 0;
        for(int k=0; k<nx2; ++k){
          sum += fwdSeed[0][d][k+i*nx2] * input[1][j+k*ny2] + input[0][k+i*nx2] * fwdSeed[1][d][j+k*ny2];
        }
        fwdSens[d][j+i*nz2] = sum;
      }
    }
  }
  
  // Adjoint sensitivities
  for(int d=0; d<nadj; ++d){
    for(int i=0; i<nx1; ++i){
      for(int j=0; j<ny2; ++j){
        for(int k=0; k<nx2; ++k){
          adjSens[0][d][k+i*nx2] += adjSeed[d][j+i*nz2]*input[1][j+k*ny2];
          adjSens[1][d][j+k*ny2] += adjSeed[d][j+i*nz2]*input[0][k+i*nx2];
        }
      }
    }
  }
}

void Multiplication::evaluateSparseSparse(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  // Carry out the actual multiplication
  for(int i=0; i<prod_map_.size(); ++i){
    output[i] = 0;
    for(std::vector< std::pair<int,int> >::const_iterator it=prod_map_[i].begin(); it!=prod_map_[i].end(); ++it)
      output[i] += input[0][it->first] * input[1][y_trans_map_[it->second]];
  }

  // Forward sensitivities: dot(Z) = dot(X)*Y + X*dot(Y)
  for(int d=0; d<nfwd; ++d){
    for(int i=0; i<prod_map_.size(); ++i){
      fwdSens[d][i] = 0;
      for(std::vector< std::pair<int,int> >::const_iterator it=prod_map_[i].begin(); it!=prod_map_[i].end(); ++it){
        fwdSens[d][i] += 
            fwdSeed[0][d][it->first] * input[1][y_trans_map_[it->second]] + 
            input[0][it->first] * fwdSeed[1][d][y_trans_map_[it->second]];
      }
    }
  }

  // Adjoint sensitivities
  for(int d=0; d<nadj; ++d){
    for(int i=0; i<prod_map_.size(); ++i){
      for(std::vector< std::pair<int,int> >::const_iterator it=prod_map_[i].begin(); it!=prod_map_[i].end(); ++it){
        adjSens[0][d][it->first] += adjSeed[d][i]*input[1][y_trans_map_[it->second]];
        adjSens[1][d][y_trans_map_[it->second]] += adjSeed[d][i]*input[0][it->first];
      }
    }
  }
}

void Multiplication::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  if(x_dense_ && y_dense_){
    evaluateDenseDense(input, output, fwdSeed, fwdSens, adjSeed, adjSens, nfwd, nadj);
  } else {
    evaluateSparseSparse(input, output, fwdSeed, fwdSens, adjSeed, adjSens, nfwd, nadj);
  }
}


} // namespace CasADi

