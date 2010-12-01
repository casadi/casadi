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

#include "unary_op.hpp"
#include <cassert>
#include <vector>
#include <sstream>

using namespace std;

namespace CasADi{

UnaryOp::UnaryOp(OPERATION op_, const MX& x) : op(op_), MXNode(x){
  sz = x.size();
}

UnaryOp* UnaryOp::clone() const{
  return new UnaryOp(*this);
}

void UnaryOp::print(std::ostream &stream) const{
  stringstream sx; sx << dep(0);
  print_c[op](stream,sx.str(),"nan");
}

void UnaryOp::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  if(fsens_order==0){
  const vector<double>& x = dep(0)->val(0);  // first argument
  vector<double>& res = val(0);
  
  for(int i=0; i<res.size(); ++i)
    nfun0[op](x[i],0,&res[i]);
  } else {

    const vector<double>& x = dep(0)->val(0);  // first argument
    const vector<double>& dx = dep(0)->val(1); // first argument derivative
    vector<double>& res = val(1);

    double tmp[3];
  
    for(int i=0; i<res.size(); ++i){
      nfun1[op](x[i],0,tmp);
      res[i] = tmp[1]*dx[i]; // chain rule
    }
  }
  
  if(asens_order>0){
    const vector<double>& x = dep(0)->val(0);  // first  argument
    vector<double>& dx = dep(0)->val(1); // first argument derivative
    const vector<double>& res = val(1);

    double tmp[3];
  
    for(int i=0; i<res.size(); ++i){
      nfun1[op](x[i],0,tmp);
      dx[i] += res[i]*tmp[1];
    }
  }
}

} // namespace CasADi

