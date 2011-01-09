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

#include "matrix_matrix_op.hpp"
#include <cassert>
#include <vector>
#include <sstream>

using namespace std;

namespace CasADi{

MatrixMatrixOp::MatrixMatrixOp(OPERATION op_, const MX& x, const MX& y) : op(op_), MXNode(x,y){
  assert(x.size1() == y.size1() || x.size2() == y.size2());
  nrow_ = x.size1();
  ncol_ = x.size2();
}

MatrixMatrixOp* MatrixMatrixOp::clone() const{
  return new MatrixMatrixOp(*this);
}

void MatrixMatrixOp::print(std::ostream &stream) const{
  stringstream sx; sx << dep(0);
  stringstream sy; sy << dep(1);
  print_c[op](stream,sx.str(),sy.str());
}

void MatrixMatrixOp::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  if(fsens_order==0){
  const vector<double>& x = input(0);  // first argument
  const vector<double>& y = input(1);  // second argument
  vector<double>& res = output();
  
  for(int i=0; i<res.size(); ++i)
    nfun0[op](x[i],y[i],&res[i]);
  } else {

    const vector<double>& x = input(0);  // first argument
    const vector<double>& dx = fwdSeed(0); // first argument derivative
    const vector<double>& y = input(1);  // second argument
    const vector<double>& dy = fwdSeed(1); // second argument derivative
    vector<double>& res = fwdSens();

  double tmp[3];

  for(int i=0; i<res.size(); ++i){
    nfun1[op](x[i],y[i],tmp);
    res[i] = tmp[1]*dx[i] + tmp[2]*dy[i]; // chain rule
  }
  }
  
  if(asens_order>0){
    const vector<double>& x = input(0);  // first  argument
    vector<double>& dx = adjSens(0); // first argument derivative
    const vector<double>& y = input(1);  // second argument
    vector<double>& dy = adjSens(1); // second argument derivative
    const vector<double>& res = adjSeed();

    double tmp[3];

    for(int i=0; i<res.size(); ++i){
      nfun1[op](x[i],y[i],tmp);
      dx[i] += res[i]*tmp[1];
      dy[i] += res[i]*tmp[2];
    }
  }
}


} // namespace CasADi

