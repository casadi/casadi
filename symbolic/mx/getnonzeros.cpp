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

#include "getnonzeros.hpp"
#include "mapping.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../fx/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace CasADi{

  GetNonzeros::GetNonzeros(const CRSSparsity& sp, const MX& y, const std::vector<int>& nz) : nz_(nz){
    setSparsity(sp);
    setDependencies(y);
  }

  GetNonzeros* GetNonzeros::clone() const{
    return new GetNonzeros(*this);
  }

  void GetNonzeros::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void GetNonzeros::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<typename T, typename MatV, typename MatVV>
  void GetNonzeros::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){

    // Number of sensitivities
    int nadj = adjSeed.size();
    int nfwd = fwdSens.size();
    
    // Nondifferentiated outputs
    const vector<T>& idata = input[0]->data();
    typename vector<T>::iterator odata_it = output[0]->begin();
    for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k){
      *odata_it++ = *k>=0 ? idata[*k] : 0;
    }
    
    // Forward sensitivities
    for(int d=0; d<nfwd; ++d){
      const vector<T>& fseed = fwdSeed[d][0]->data();
      typename vector<T>::iterator fsens_it = fwdSens[d][0]->begin();
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k){
     	*fsens_it++ = *k>=0 ? fseed[*k] : 0;
      }
    }
      
    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      typename vector<T>::iterator aseed_it = adjSeed[d][0]->begin();
      vector<T>& asens = adjSens[d][0]->data();
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k){
	if(*k>=0) asens[*k] += *aseed_it;
	*aseed_it++ = 0;
      }
    }
  }
  
  void GetNonzeros::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    // Get references to the assignment operations and data
    bvec_t *outputd = get_bvec_t(output[0]->data());
    bvec_t *inputd = get_bvec_t(input[0]->data());
    
    // Propate sparsity
    if(fwd){
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k){
	*outputd++ = *k>=0 ? inputd[*k] : 0;
      }
    } else {
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k){
	if(*k>=0) inputd[*k] |= *outputd;
	*outputd++ = 0;
      }
    }
  }

  void GetNonzeros::printPart(std::ostream &stream, int part) const{
    switch(part){
    case 1: stream << nz_; break;
    }
  }

  void GetNonzeros::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){

    // Output sparsity
    const CRSSparsity &osp = sparsity();
    const vector<int>& ocol = osp.col();
    vector<int> orow = osp.getRow();
    
    // Input sparsity
    const CRSSparsity &isp = dep().sparsity();
    const vector<int>& icol = isp.col();
    vector<int> irow = isp.getRow();

    // Number of derivative directions
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Sparsity patterns in sparse triplet format of quantities being calculated
    vector<int> r_row, r_col, r_nz;
  
    // Find out which matrix elements that we are trying to calculate
    vector<int> el_wanted(nz_.size());
    for(int k=0; k<nz_.size(); ++k){
      el_wanted[k] = orow[k] + ocol[k]*osp.size1();
    }
    
    // We next need to resort the assignment vector by inputs instead of outputs
    // Start by counting the number of input nonzeros corresponding to each output nonzero
    vector<int> inz_count(icol.size()+1,0);
    for(vector<int>::const_iterator it=nz_.begin(); it!=nz_.end(); ++it){
      casadi_assert_message(*it>=0,"Not implemented");
      inz_count[*it+1]++;
    }
    
    // Cumsum to get index offset for input nonzero
    for(int i=0; i<icol.size(); ++i){
      inz_count[i+1] += inz_count[i];
    }
    
    // Get the order of assignments
    vector<int> nz_order(nz_.size());
    for(int k=0; k<nz_.size(); ++k){
      // Save the new index
      nz_order[inz_count[nz_[k]]++] = k;
    }
    
    // Find out which elements are given
    vector<int>& el_known = inz_count; // Reuse memory
    el_known.resize(nz_.size());
    for(int k=0; k<nz_.size(); ++k){
      // Get output nonzero
      int inz_k = nz_[nz_order[k]];
      
      // Get element
      el_known[k] = irow[inz_k] + icol[inz_k]*isp.size1();
    }
    
    // Temporary vector
    vector<int> temp(el_known);
    
    // Evaluate the nondifferentiated function
    if(!output_given){
      
      // Get the matching nonzeros (temp==el_known)
      input[0]->sparsity().getNZInplace(temp);
      
      // Add to sparsity pattern
      for(int k=0; k<nz_.size(); ++k){
	if(temp[k]!=-1){
	  r_nz.push_back(temp[k]);
	  r_col.push_back(ocol[nz_order[k]]);
	  r_row.push_back(orow[nz_order[k]]);
	}
      }

      // Create a sparsity pattern from vectors
      CRSSparsity r_sp = sp_triplet(osp.size1(),osp.size2(),r_row,r_col,temp);      
      if(r_nz.size()==0){
	*output[0] = MX::zeros(r_sp);
      } else {
	*output[0] = (*input[0])->getGetNonzeros(r_sp,r_nz);
      }
    }
    
    // Forward sensitivities
    for(int d=0; d<nfwd; ++d){
      
      // Get the matching nonzeros
      temp.resize(el_known.size());
      copy(el_known.begin(),el_known.end(),temp.begin());
      fwdSeed[d][0]->sparsity().getNZInplace(temp);
      
      // Add to sparsity pattern
      r_nz.clear();
      r_col.clear();
      r_row.clear();
      for(int k=0; k<nz_.size(); ++k){
	if(temp[k]!=-1){
	  r_nz.push_back(temp[k]);
	  r_col.push_back(ocol[nz_order[k]]);
	  r_row.push_back(orow[nz_order[k]]);
	}
      }

      // Create a sparsity pattern from vectors
      CRSSparsity f_sp = sp_triplet(osp.size1(),osp.size2(),r_row,r_col,temp);
      if(r_nz.size()==0){
	*fwdSens[d][0] = MX::zeros(f_sp);
      } else {
	*fwdSens[d][0] = (*fwdSeed[d][0])->getGetNonzeros(f_sp,r_nz);
      }
    }
    
    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){

      // Get the matching nonzeros
      r_nz.resize(el_wanted.size());
      copy(el_wanted.begin(),el_wanted.end(),r_nz.begin());
      adjSeed[d][0]->sparsity().getNZInplace(r_nz);
      
      // Check if nothing to add
      bool nothing_to_add = true;
      for(int i=0; i<nz_.size(); ++i){
	if(nz_[i]>=0 && r_nz[i]>=0){
	  nothing_to_add = false;
	  break;
	}
      }

      // Quick return if nothing to add
      if(!nothing_to_add){

	for(int iter=0; iter<2; ++iter){
	  // Use r_row as temporary memory
	  vector<int>& temp2 = r_row;

	  // Get the corresponding output 
	  temp2.resize(el_known.size());
	  copy(el_known.begin(),el_known.end(),temp2.begin());
	  adjSens[d][0]->sparsity().getNZInplace(temp2);

	  // Resort in the order of the inputs
	  temp.resize(temp2.size());
	  for(int i=0; i<nz_order.size(); ++i){
	    temp[nz_order[i]] = temp2[i];
	  }

	  // Check if any additions aren't included in the current value of the sensitivity
	  bool spilled = false;
	  for(int i=0; i<nz_.size(); ++i){
	    if(r_nz[i]>=0 && nz_[i]>=0 && temp[i]<0){
	      spilled = true;
	      break;
	    }
	  }

	  // All additions fit
	  if(!spilled) break;

	  // Densify the sensitivitity and make another loop (never more than two needed)
	  casadi_assert(iter<2);

	  // Create a new pattern which includes both the the previous seed and the addition
	  vector<unsigned char> tmp1;
	  CRSSparsity sp = adjSens[d][0]->sparsity().patternUnion(dep().sparsity(),tmp1);
	  MX t = (*adjSens[d][0])->getDensification(sp);
	  *adjSens[d][0] = t;
	}

	// Compress the nonzeros to the size of the seed
	int n=0;	
	for(int i=0; i<nz_.size(); ++i){
	  if(r_nz[i]>=0){
	    r_nz[n++] = nz_[i]>=0 ? temp[i] : -1;
	  }
	}
	r_nz.resize(n);

	// Add to the element
	*adjSens[d][0] = (*adjSeed[d][0])->getAddNonzeros(*adjSens[d][0],r_nz);
      }
      
      // Clear adjoint seeds
      *adjSeed[d][0] = MX();
    }
  }
  
  Matrix<int> GetNonzeros::mapping(int iind) const {
    return Matrix<int>(sparsity(),nz_);
  }

  bool GetNonzeros::isIdentity() const{
    // Check sparsity
    if(!(sparsity() == dep().sparsity()))
      return false;
      
    // Check if the nonzeros follow in increasing order
    for(int k=0; k<nz_.size(); ++k){
      if(nz_[k] != k) return false;
    }
    
    // True if reached this point
    return true;
  }

  void GetNonzeros::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    if(nz_.size()==1 && nz_.front()>=0){
      // Compact if just a scalar (TODO: extend to slices)
      stream << "  " << res.front() << "[0]=" << arg.front() << "[" << nz_.front() << "];" << endl;
    } else {
      // Condegen the indices
      int ind = gen.getConstant(nz_,true);
      
      // Codegen the assignments
      stream << "  for(ii=s" << ind << ", rr=" << res.front() << ", ss=" << arg.front() << "; ii!=s" << ind << "+" << nz_.size() << "; ++ii) *rr++ = *ii>=0 ? ss[*ii] : 0;" << endl;
    }
  }

  void GetNonzeros::simplifyMe(MX& ex){
    // Simplify if identity
    if(isIdentity()){
      MX t = dep(0);
      ex = t;
    }
  }

  MX GetNonzeros::getGetNonzeros(const CRSSparsity& sp, const std::vector<int>& nz) const{
    // Eliminate recursive calls
    vector<int> nz_new(nz);
    for(vector<int>::iterator i=nz_new.begin(); i!=nz_new.end(); ++i){
      *i = nz_[*i];
    }
    return dep()->getGetNonzeros(sp,nz_new);
  }

} // namespace CasADi
