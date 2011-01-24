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
  setDependencies(x,y);
  if(x.size2() != y.size1()) throw CasadiException("Multiplication::dimension mismatch");
  setSparsity(CRSSparsity(x.size1(),y.size2(),true));
}

Multiplication* Multiplication::clone() const{
  return new Multiplication(*this);
}

void Multiplication::print(std::ostream &stream) const{
  stream << "prod(" << dep(0) << "," << dep(1) << ")";
}

void Multiplication::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  if(dep(0).size() == dep(0).numel() && dep(1).size() == dep(1).numel()){
    // Dense
  
    // Get dimensions
    int nx1 = dep(0).size1();
    int nx2 = dep(0).size2();
    int ny1 = dep(1).size1();
    int ny2 = dep(1).size2();
    int nz1 = nx1;
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
            adjSens[0][d][j+k*ny2] += adjSeed[d][j+i*nz2]*input[1][k+i*nx2];
            adjSens[1][d][j+k*ny2] += adjSeed[d][j+i*nz2]*input[0][k+i*nx2];
          }
        }
      }
    }
  } else {
    // Sparse
    throw CasadiException("sparse matrix multiplication not implemented in MX");
  }
}


} // namespace CasADi

