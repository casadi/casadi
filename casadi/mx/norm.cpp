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

#include "norm.hpp"


#include "mx_tools.hpp"

using namespace std;
namespace CasADi{

Norm::Norm(const MX& x){
  setDependencies(x);
  setSparsity(CRSSparsity(1,1,true));
}

void Norm::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  throw CasadiException("Norm::evaluate not implemented");
}

void Norm::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  throw CasadiException("Norm::evaluateSX not implemented");
}

void Norm::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  throw CasadiException("Norm::evaluateMX not implemented");
}

MX Norm::adFwd(const std::vector< MX > & jx) {
  // Number of derivative directions
  int ndir = jx[0].size2();

  // Return a not a number
  return MX(1,ndir,numeric_limits<double>::quiet_NaN());
}


Norm2::Norm2(const MX& x) : Norm(x){
}

Norm2* Norm2::clone() const{
  return new Norm2(*this);
}

void Norm2::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  vector<double> &outputd = output[0]->data();
  const vector<double> &inputd = input[0]->data();

  double temp=0;
  for (int k=0;k<dep(0).size();k++) {
    temp+=inputd[k]*inputd[k];
  }
  outputd[0]=sqrt(temp);
  // Propagate forward seeds
  for(int d=0; d<nfwd; ++d){
    fwdSens[d][0]->data()[0]=0;
    for(int k=0; k<dep(0).size(); k++){
      fwdSens[d][0]->data()[0] += inputd[k]/outputd[0] * fwdSeed[d][0]->data()[k];
    }
  }

  // Propagate adjoint seeds
  for(int d=0; d<nadj; ++d){
    if (adjSeed[d][0]->data()[0]==0) continue;
    for(int k=0; k<dep(0).size(); k++){
      adjSens[d][0]->data()[k] += inputd[k]/outputd[0] * adjSeed[d][0]->data()[0];
    }
  }
  

}

MX Norm2::adFwd(const std::vector< MX > & jx	) {
  MX ret;
  ret.assignNode(this);
  return trans(prod(jx.at(0),dep(0)))/ret;
}

void Norm2::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "||" << args.at(0) << "||_2"; 
}

Norm1::Norm1(const MX& x) : Norm(x){
}

Norm1* Norm1::clone() const{
  return new Norm1(*this);
}

void Norm1::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "||" << args.at(0) << "||_1"; 
}

void Norm1::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  vector<double> &outputd = output[0]->data();
    const vector<double> &inputd = input[0]->data();

  if(nfwd==0 && nadj==0){
   double temp=0;
   for (int k=0;k<dep(0).size();k++) temp+=std::abs(inputd[k]);
   outputd[0]=temp;
   return; 
  }

  // Propagate forward seeds
  for(int d=0; d<nfwd; ++d){
    fwdSens[d][0]->data()[0]=0;
    for(int k=0; k<dep(0).size(); k++){
      if (fwdSeed[d][0]->data()[k]==0) continue;
      if (inputd[k] < 0) {
        fwdSens[d][0]->data()[0] -= fwdSeed[d][0]->data()[k];
      } else if (inputd[k] > 0) {
        fwdSens[d][0]->data()[0] += fwdSeed[d][0]->data()[k];
      } else {
        fwdSens[d][0]->data()[0] += std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  // Propagate adjoint seeds
  for(int d=0; d<nadj; ++d){
    if (adjSeed[d][0]->data()[0]==0) continue;
    for(int k=0; k<dep(0).size(); k++){
      if (inputd[k] < 0) {
        adjSens[d][0]->data()[k] -=  adjSeed[d][0]->data()[0];
      } else if (inputd[k] > 0) {
        adjSens[d][0]->data()[k] +=  adjSeed[d][0]->data()[0];
      } else {
        adjSens[d][0]->data()[k] += std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

}

NormInf::NormInf(const MX& x) : Norm(x){
}

NormInf* NormInf::clone() const{
  return new NormInf(*this);
}

void NormInf::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "||" << args.at(0) << "||_inf"; 
}

void NormInf::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  vector<double> &outputd = output[0]->data();
  const vector<double> &inputd = input[0]->data();

  double temp=std::numeric_limits<double>::infinity();
  double a;
  for (int k=0;k<dep(0).size();k++) {
    a = std::abs(inputd[k]);
    if (a>temp) temp=a;
  }
  outputd[0]=temp;
  

  if (nadj!=0) throw CasadiException("NormInf::evaluate not implemented");
}

} // namespace CasADi

