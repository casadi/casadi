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

#include "reordering.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

Reordering::Reordering(){
}

void Reordering::init(){
  // Call the base class initialization
  MXNode::init();
  
  // Allocate indices if this has not been done
  if(nzind_.empty()){
    nzind_.resize(output().size());
    fill(nzind_.begin(),nzind_.end(),0);
  }
  
  if(argind_.empty()){
    argind_.resize(output().size());
    fill(argind_.begin(),argind_.end(),0);
  }
  if(nzind_.empty()){
    nzind_.resize(output().size());
    fill(nzind_.begin(),nzind_.end(),0);
  }
  
  if(argind_.size() != output().size() || nzind_.size() != output().size()) 
    throw CasadiException("Reordering::init: dimension mismatch");
    
}


void Reordering::evaluate(int fsens_order, int asens_order){
  if(fsens_order==0){
    vector<double>& res = output();
    for(int k=0; k<res.size(); ++k)
      res[k] = input(k2l(k))[k2k(k)];
  } else {
    vector<double>& fsens = fwdSens();
    for(int k=0; k<fsens.size(); ++k)
      fsens[k] = fwdSeed(k2l(k))[k2k(k)];
  }
  
  if(asens_order>0){
    const vector<double>& aseed = adjSeed();
    for(int k=0; k<aseed.size(); ++k)
      adjSens(k2l(k))[k2k(k)] += aseed[k];
  }
}

int Reordering::k2l(int k) const{
  return argind_[k];
}

int Reordering::k2k(int k) const{
  return nzind_[k];
}

void Reordering::print(std::ostream &stream) const{
  stream << "reordering(" << dep(0) << "," << nzind_;
  if(ndep()>1) stream << "," << argind_;
  stream << ")";
}


} // namespace CasADi
