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

#include "vertcat.hpp"
#include "../stl_vector_tools.hpp"
#include <iterator>
#include <algorithm>

using namespace std;

namespace CasADi{

// Constructor
Vertcat::Vertcat(const vector<MX>& d){
  casadi_assert(d.size()>=2);
  setDependencies(d);
  int sz1=0;
  int sz2=d[0].size2();
  for(int i=0; i<ndep(); ++i){
    sz1 += d[0].size1();
    if(sz2!=d[i].size2())
      throw CasadiException("Vertcat: dimension mismatch");
  }
  
/*  CRSSparsity sp = d[0].sparsity();
  for(int i=1; i<d.size(); ++i){
    sp.append(d[i].sparsity());
  }*/
  
  
  // 
  setSparsity(CRSSparsity(sz1,sz2,true));
}

Vertcat* Vertcat::clone() const{
  return new Vertcat(*this);
}

void Vertcat::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "[";
  stream << args.at(0);
  for(int i=1; i<args.size(); ++i)
    stream << " ; " << args.at(i);
  stream << "]";
}

void Vertcat::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  int offset = 0;
  for(int i=0; i<ndep(); ++i){
    // Number of rows
    int ni = dep(i).size();
    
    // Copy values
    for(int k=0; k<ni; ++k)
      output[k+offset] = input[i][k];
    
    // Transmit forward derivatives
    for(int d=0; d<nfwd; ++d)
      for(int k=0; k<ni; ++k)
        fwdSens[d][k+offset] = fwdSeed[i][d][k];
    
    // Transmit adjoint derivatives
    for(int d=0; d<nadj; ++d)
      for(int k=0; k<ni; ++k)
        adjSens[i][d][k] += adjSeed[d][k+offset];
    
    // Update offset
    offset += ni;
  }
}


} // namespace CasADi
