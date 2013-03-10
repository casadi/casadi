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

#include "setnonzeros.hpp"
#include "mapping.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../fx/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace CasADi{

  SetNonzeros::SetNonzeros(const MX& y, const MX& x, const std::vector<int>& nz) : nz_(nz){
    setSparsity(y.sparsity());
    setDependencies(y,x);
    casadi_assert(nz.size()==x.size());
  }

  SetNonzeros* SetNonzeros::clone() const{
    return new SetNonzeros(*this);
  }

  void SetNonzeros::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void SetNonzeros::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<typename T, typename MatV, typename MatVV>
  void SetNonzeros::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){

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
      if(*k>=0) odata[*k] = *idata_it;
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
	if(*k>=0) fsens[*k] = *fseed_it;
      }
    }
      
    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      vector<T>& aseed = adjSeed[d][0]->data();
      vector<T>& asens0 = adjSens[d][0]->data();
      typename vector<T>::iterator asens_it = adjSens[d][1]->begin();
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k, ++asens_it){
	if(*k>=0){
	  *asens_it += aseed[*k];
	  aseed[*k] = 0;
	}
      }
      if(&aseed != &asens0){
	transform(aseed.begin(),aseed.end(),asens0.begin(),asens0.begin(),std::plus<T>());
	fill(aseed.begin(),aseed.end(),0);
      }
    }
  }
  
  void SetNonzeros::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
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
	if(*k>=0) outputd[*k] = *inputd;
      }
    } else {
      for(vector<int>::const_iterator k=nz_.begin(); k!=nz_.end(); ++k, ++inputd){
	if(*k>=0){
	  *inputd |= outputd[*k];
	  outputd[*k] = 0;
	}
      }
      if(outputd != inputd0){
	for(int k=0; k<input[0]->size(); ++k){
	  inputd0[k] |= outputd[k];
	  outputd[k] = 0;
	}
      }
    }
  }

  void SetNonzeros::printPart(std::ostream &stream, int part) const{
    switch(part){
    case 1: stream << nz_ << " = "; break;
    }
  }

  void SetNonzeros::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    casadi_error("not implemented");

    // Output sparsity
    const CRSSparsity &osp = sparsity();
    const vector<int>& ocol = osp.col();
    vector<int> orow = osp.getRow();

    // Input sparsity (first input same as output)
    const CRSSparsity &isp = dep(1).sparsity();
    const vector<int>& icol = isp.col();
    vector<int> irow = isp.getRow();

    // Number of derivative directions
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Sparsity patterns in sparse triplet format of quantities being calculated
    vector<int> r_row, r_col, r_nz;
  
    // Find out which matrix elements that we are trying to add
    vector<int> el_input(nz_.size());
    for(int k=0; k<nz_.size(); ++k){
      el_input[k] = irow[k] + icol[k]*isp.size1();
    }
    
    // We next need to resort the addition vector by outputs instead of inputs
    // Start by counting the number of outputs nonzeros corresponding to each input nonzero
    vector<int> onz_count(ocol.size()+1,0);
    for(vector<int>::const_iterator it=nz_.begin(); it!=nz_.end(); ++it){
      onz_count[*it+1]++;
    }
    
    // Cumsum to get index offset for output nonzero
    for(int i=0; i<ocol.size(); ++i){
      onz_count[i+1] += onz_count[i];
    }
    
    // Get the order of assignments
    vector<int> nz_order(nz_.size());
    for(int k=0; k<nz_.size(); ++k){
      // Save the new index
      nz_order[onz_count[nz_[k]]++] = k;
    }
    
    // Find out the destination of the elements
    vector<int>& el_output = onz_count; // Reuse memory
    el_output.resize(nz_.size());
    for(int k=0; k<nz_.size(); ++k){
      // Get input nonzero
      int onz_k = nz_[nz_order[k]];
      
      // Get element
      el_output[k] = orow[onz_k] + ocol[onz_k]*osp.size1();
    }
    
    // Temporary vector
    vector<int> temp(el_input);
    
    // Evaluate the nondifferentiated function
    if(!output_given){

      // Get the matching nonzeros in the source
      input[1]->sparsity().getNZInplace(temp);
      
      // Add to sparsity pattern
      for(int k=0; k<nz_.size(); ++k){
	if(temp[k]!=-1){
	  r_nz.push_back(temp[k]);
	  r_col.push_back(icol[k]);
	  r_row.push_back(irow[k]);
	}
      }
      
      // Check if anything needs to be added
      if(r_nz.size()==0){

	// Simple assignment of nothing to be added
	*output[0] = *input[0];
      } else {

	// Get argument to be added with ignored entries removed
	MX arg;
	if(temp.size()==nz_.size()){
	  // No need to drop any entries
	  arg = *input[1];
	} else {
	  // Drop entries that do will not be used
	  CRSSparsity r_sp = sp_triplet(isp.size1(),isp.size2(),r_row,r_col,temp);
	  arg = (*input[1])->getGetNonzeros(r_sp,r_nz);
	}
      }

     

// Get the matching nonzeros in the destination
//       temp.resize(el_destination.size());
//       copy(el_destination.begin(),el_destination.end(),temp.begin());
//       input[0]->sparsity().getNZInplace(temp);




//       if(r_nz.size()==0){
// 	// Nothing to add
// 	*output[0] = *input[0];
//       } else {


	
// 	// Now get the destination in the output vector
// 	for(vector<int>::iterator i=r_nz.begin(); i!=r_nz.end(); ++i){
// 	  *i = el_destination[*i];
// 	}
	//	temp.resize(el_known.size());
	//	copy(el_known.begin(),el_known.end(),temp.begin());
	//	arg->sparsity().getNZInplace(temp);
	

    }

    casadi_error("not ready");

#if 0
    
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
      temp.resize(el_wanted.size());
      copy(el_wanted.begin(),el_wanted.end(),temp.begin());
      adjSeed[d][0]->sparsity().getNZInplace(temp);
      
      // Add to sparsity pattern
      r_nz.clear();
      r_col.clear();
      r_row.clear();
      for(int k=0; k<nz_.size(); ++k){
	if(temp[k]!=-1){
	  r_nz.push_back(temp[k]);
	  r_col.push_back(icol[nz_[k]]);
	  r_row.push_back(irow[nz_[k]]);
	}
      }

      // Create a sparsity pattern from vectors
      CRSSparsity a_sp = sp_triplet(isp.size1(),isp.size2(),r_row,r_col,temp,true);
      if(r_nz.size()>0){
	// Create a mapping matrix
	MX s = MX::create(new Mapping(a_sp));
	s->assign(*adjSeed[d][0],r_nz,temp,true);
	
	// Save to adjoint sensitivities
	*adjSens[d][0] += s;
      }
      
      // Clear adjoint seeds
      *adjSeed[d][0] = MX();
    }

#endif

  }
  
  Matrix<int> SetNonzeros::mapping(int iind) const {
    return Matrix<int>(sparsity(),nz_);
  }

  bool SetNonzeros::isAssignment() const{
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

  void SetNonzeros::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    // Check if inplace
    bool inplace = arg.at(0).compare(res.front())==0;

    // Copy first argument if not implace
    if(!inplace){      
      stream << "  for(i=0; i<" << size() << "; ++i) " << res.front() << "[i]=" << arg.at(0) << "[i];" << endl;
    }

    if(nz_.size()==1){
      // Compact if just a scalar (TODO: extend to slices)
      stream << "  " << res.front() << "[" << nz_.front() << "]=" << arg.at(1) << "[0];" << endl;
    } else {
      // Condegen the indices
      int ind = gen.getConstant(nz_,true);
      
      // Codegen the assignments
      stream << "  for(ii=s" << ind << ", rr=" << res.front() << ", ss=" << arg.at(1) << "; ii!=s" << ind << "+" << nz_.size() << "; ++ii, ++ss) if(*ii>=0) rr[*ii] = *ss;" << endl;
    }
  }

  void SetNonzeros::simplifyMe(MX& ex){
    // Simplify if addition
    if(isAssignment()){
      MX t = dep(1);
      ex = t;
    }
  }

} // namespace CasADi
