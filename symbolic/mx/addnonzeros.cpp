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

#include "addnonzeros.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../fx/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace CasADi{

  AddNonzeros::AddNonzeros(const MX& y, const MX& x, const std::vector<int>& nz) : NonzerosBase(nz){
    setSparsity(y.sparsity());
    setDependencies(y,x);
    casadi_assert(nz.size()==x.size());
  }

  AddNonzeros* AddNonzeros::clone() const{
    return new AddNonzeros(*this);
  }

  void AddNonzeros::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void AddNonzeros::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<typename T, typename MatV, typename MatVV>
  void AddNonzeros::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){

    // Number of sensitivities
    int nadj = adjSeed.size();
    int nfwd = fwdSens.size();
    
    // Nondifferentiated outputs
    const vector<T>& idata0 = input[0]->data();
    typename vector<T>::const_iterator idata_it = input[1]->begin();
    vector<T>& odata = output[0]->data();
    if(&idata0 != &odata){
      copy(idata0.begin(),idata0.end(),odata.begin());
    }
    for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k, ++idata_it){
      if(*k>=0) odata[*k] += *idata_it;
    }
    
    // Forward sensitivities
    for(int d=0; d<nfwd; ++d){
      const vector<T>& fseed0 = fwdSeed[d][0]->data();
      typename vector<T>::const_iterator fseed_it = fwdSeed[d][1]->begin();
      vector<T>& fsens = fwdSens[d][0]->data();
      if(&fseed0 != &fsens){
	copy(fseed0.begin(),fseed0.end(),fsens.begin());
      }
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k, ++fseed_it){
	if(*k>=0) fsens[*k] += *fseed_it;
      }
    }
      
    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      vector<T>& aseed = adjSeed[d][0]->data();
      vector<T>& asens0 = adjSens[d][0]->data();
      typename vector<T>::iterator asens_it = adjSens[d][1]->begin();
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k, ++asens_it){
	if(*k>=0) *asens_it += aseed[*k];
      }
      if(&aseed != &asens0){
	transform(aseed.begin(),aseed.end(),asens0.begin(),asens0.begin(),std::plus<T>());
	fill(aseed.begin(),aseed.end(),0);
      }
    }
  }
  
  void AddNonzeros::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    // Get references to the assignment operations and data
    bvec_t *outputd = get_bvec_t(output[0]->data());
    bvec_t *inputd0 = get_bvec_t(input[0]->data());
    bvec_t *inputd = get_bvec_t(input[1]->data());

    // Propate sparsity
    if(fwd){
      if(outputd != inputd0){
	copy(inputd0,inputd0+input[0]->size(),outputd);
      }
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k, ++inputd){
	if(*k>=0) outputd[*k] |= *inputd;
      }
    } else {
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k, ++inputd){
	if(*k>=0) *inputd |= outputd[*k];
      }
      if(outputd != inputd0){
	for(int k=0; k<input[0]->size(); ++k){
	  inputd0[k] |= outputd[k];
	  outputd[k] = 0;
	}
      }
    }
  }

  void AddNonzeros::printPart(std::ostream &stream, int part) const{
    switch(part){
    case 1: stream << nz_ << " += "; break;
    }
  }

  void AddNonzeros::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    evaluateMXBase(input,output,fwdSeed,fwdSens,adjSeed,adjSens,output_given,true);
  }
  
  Matrix<int> AddNonzeros::mapping(int iind) const {
    return Matrix<int>(sparsity(),nz_);
  }

  bool AddNonzeros::isAddition() const{
    // Check sparsity
    if(!(sparsity() == dep(1).sparsity()))
      return false;
      
    // Check if the nonzeros follow in increasing order
    for(int k=0; k<nz_.size(); ++k){
      if(nz_[k] != k) return false;
    }
    
    // True if reached this point
    return true;
  }

  void AddNonzeros::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    // Check if inplace
    bool inplace = arg.at(0).compare(res.front())==0;

    // Copy first argument if not implace
    if(!inplace){      
      stream << "  for(i=0; i<" << size() << "; ++i) " << res.front() << "[i]=" << arg.at(0) << "[i];" << endl;
    }

    if(Slice::isSlice(nz_)){
      // Compact if a Slice (TODO: Move to separate class)
      Slice s(nz_);
      stream << "  for(rr=" << res.front() << "+" << s.start_ << ", ss=" << arg.at(1) << "; rr!=" << res.front() << "+" << s.stop_ << "; rr+=" << s.step_ << ") *rr += *ss++;" << endl;
    } else {
      // Condegen the indices
      int ind = gen.getConstant(nz_,true);
      
      // Codegen the additions
      stream << "  for(ii=s" << ind << ", rr=" << res.front() << ", ss=" << arg.at(1) << "; ii!=s" << ind << "+" << nz_.size() << "; ++ii, ++ss) if(*ii>=0) rr[*ii] += *ss;" << endl;
    }
  }

  void AddNonzeros::simplifyMe(MX& ex){
    // Simplify if addition
    if(isAddition()){
      MX t = dep(1);
      ex += t;
    }
  }

} // namespace CasADi
