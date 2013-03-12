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

#include "nonzeros_base.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../fx/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace CasADi{

  NonzerosBase:: ~NonzerosBase(){
  }

  void NonzerosBase::evaluateMXBase(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given, bool add){

    // Number of derivative directions
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Output sparsity
    const CRSSparsity &osp = sparsity();
    const vector<int>& ocol = osp.col();
    vector<int> orow = osp.getRow();
    
    // Input sparsity (first input same as output)
    const CRSSparsity &isp = dep(1).sparsity();
    const vector<int>& icol = isp.col();
    vector<int> irow = isp.getRow();
      
    // Get all output elements (this time without duplicates)
    vector<int> el_output;
    osp.getElements(el_output,false);
    
    // Temporary for sparsity pattern unions
    vector<unsigned char> tmp1;

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<int> r_rowind, r_col, r_nz;

    // Nondifferentiated function and forward sensitivities
    int first_d = output_given ? 0 : -1;
    for(int d=0; d<nfwd; ++d){

      // Get references to arguments and results
      MX& arg = d<0 ? *input[1] : *fwdSeed[d][1];
      MX arg0 = d<0 ? *input[0] : *fwdSeed[d][0];
      MX& res = d<0 ? *output[0] : *fwdSens[d][0];

      // Get all the arguments
      arg.sparsity().getElements(r_nz,false);

      // Get the corresponding nz locations in the input sparsity pattern
      isp.getNZInplace(r_nz);

      // Filter out ignored entries and check if there is anything to set at all
      bool elements_to_set = false;
      for(vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k){
	if(*k>=0){
	  if(nz_[*k]>=0){
	    elements_to_set = true;
	  } else {
	    *k = -1;
	  }
	}
      }

      // Quick continue of no elements to set
      if(!elements_to_set) continue;
     
      // Get the nz locations in the argument corresponding to the inputs
      vector<int> &r_nz2 = r_col; // Reuse memory
      r_nz2.resize(el_output.size());
      copy(el_output.begin(),el_output.end(),r_nz2.begin());
      arg0.sparsity().getNZInplace(r_nz2);
      
      // Enlarge the sparsity pattern of the arguments if not all assignments fit
      for(vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k){
	if(*k>=0 && r_nz2[nz_[*k]]<0){
	  
	  // Create a new pattern which includes both the the previous seed and the addition
	  CRSSparsity sp = arg0.sparsity().patternUnion(sparsity(),tmp1);
	  arg0 = arg0->getDensification(sp);

	  // Recalculate the nz locations in the arguments corresponding to the inputs
	  copy(el_output.begin(),el_output.end(),r_nz2.begin());
	  arg0.sparsity().getNZInplace(r_nz2);

	  break;
	}
      }

      // Have r_nz point to locations in the result instead of the output
      for(vector<int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k){
	if(*k>=0){
	  *k = r_nz2[nz_[*k]];
	}
      }

      // Add to the element to the sensitivity
      if(add){
	res = arg->getAddNonzeros(arg0,r_nz);
      } else {
	res = arg->getSetNonzeros(arg0,r_nz);
      }
    }

    // Quick return if no adjoints
    if(nadj==0) return;

    // We next need to resort the assignment vector by outputs instead of inputs
    // Start by counting the number of output nonzeros corresponding to each input nonzero
    vector<int>& onz_count = el_output;
    onz_count.resize(ocol.size()+1);
    fill(onz_count.begin(),onz_count.end(),0);
    for(vector<int>::const_iterator it=nz_.begin(); it!=nz_.end(); ++it){
      casadi_assert_message(*it>=0,"Not implemented");
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

    // Find out which elements are being set
    el_output.resize(nz_.size());
    for(int k=0; k<nz_.size(); ++k){
      // Get output nonzero
      int onz_k = nz_[nz_order[k]];
      
      // Get element (note: may contain duplicates)
      el_output[k] = orow[onz_k] + ocol[onz_k]*osp.size1();
    }

    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){

      // Get an owning references to the seeds and sensitivities and clear the seeds for the next run
      MX aseed = *adjSeed[d][0];      
      *adjSeed[d][0] = MX();
      MX& asens0 = *adjSens[d][0];
      MX& asens = *adjSens[d][1];
      
      // Get the matching nonzeros
      r_nz.resize(el_output.size());
      copy(el_output.begin(),el_output.end(),r_nz.begin());
      aseed.sparsity().getNZInplace(r_nz);
      
      // Add to sparsity pattern
      int n=0;
      r_col.clear();
      r_rowind.resize(isp.size1()+1); // Row count
      fill(r_rowind.begin(),r_rowind.end(),0);
      for(int k=0; k<nz_.size(); ++k){
	if(r_nz[k]!=-1){
	  r_nz[n++] = r_nz[k];
	  r_col.push_back(icol[nz_order[k]]);
	  r_rowind[1+irow[nz_order[k]]]++;
	}
      }
      r_nz.resize(n);
      for(int i=1; i<r_rowind.size(); ++i) r_rowind[i] += r_rowind[i-1]; // row count -> row offset

      if(r_nz.size()==0){
	// Nothing to set/add
	asens0 = aseed;
      } else {
	// Create a sparsity pattern from vectors
	CRSSparsity f_sp(isp.size1(),isp.size2(),r_col,r_rowind);
	asens += aseed->getGetNonzeros(f_sp,r_nz);

	if(add){
	  // The corresponding nonzeros remain in the seed
	  asens0 = aseed;
	} else {
	  // The corresponding nonzeros disappear from the seed
	  asens0 = MX::zeros(f_sp)->getSetNonzeros(aseed,r_nz);
	}
      }
    }
  }
  
} // namespace CasADi
