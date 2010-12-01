/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "mx_constant.hpp"
#include <cassert>
#include <vector>
#include <algorithm>
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

MXConstant::MXConstant(const double *x, int n, int m, char order){
  sz.nrow = n;
  sz.ncol = m;
  data.resize(n*m);

  if(order=='R'){
    for(int i=0; i<n; ++i)
      for(int j=0; j<m; ++j)
	data[i + n*j] = x[i + n*j];
  } else {
    assert(order=='C');
    for(int i=0; i<n; ++i)
      for(int j=0; j<m; ++j)
	data[i + n*j] = x[j + m*i];
  }
}

MXConstant::MXConstant(int n, int m){
  sz.nrow = n;
  sz.ncol = m;
  data.resize(n*m);
}

MXConstant* MXConstant::clone() const{
  return new MXConstant(*this);
}

void MXConstant::print(std::ostream &stream) const{
  if(sz.nrow==1 && sz.ncol==1)
    stream << data;
  else if(sz.nrow==1){
    stream << "[";
    stream << data[0];
    for(int i=1; i<sz.ncol; ++i)
      stream << "," << data[i];
     stream << "]";
  } else {
    stream << "[";
    for(int i=0; i<sz.nrow; ++i){
      stream << data[i];
      for(int j=1; j<sz.ncol; ++j){
	stream << "," << data[i + sz.nrow*j];
      }
      if(i!=sz.nrow-1)
	stream << ";";
    }    
    stream << "]";
  }
}

void MXConstant::evaluate(int fsens_order, int asens_order){
  if(fsens_order==0){
    val(0) = data;
  } else {
    vector<double>& res = val(1);
    for(vector<double>::iterator it=res.begin(); it!=res.end(); ++it) 
      *it = 0;
  }
}

// void MXConstant::evaluateAdj(){
// 
// }


double MXConstant::operator()(int i, int j) const{
  assert(i<sz.nrow && j<sz.ncol);
  return data[i+sz.nrow*j];
}

double& MXConstant::operator()(int i, int j){
  assert(i<sz.nrow && j<sz.ncol);
  return data[i+sz.nrow*j];  
}

double MXConstant::operator[](int k) const{
  return data[k];
}

double& MXConstant::operator[](int k){
  return data[k];
}

ostream& operator<<(ostream &stream, const MXConstant& x){
  return stream << x.data;
}

bool MXConstant::isConstant() const{
  return true;
}


} // namespace CasADi

