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
    vector<T>& odata = output[0]->data();
    for(int k=0; k<nz_.size(); ++k){
      odata[k] = idata[nz_[k]];
    }
    
    // Forward sensitivities
    for(int d=0; d<nfwd; ++d){
      const vector<T>& fseed = fwdSeed[d][0]->data();
      vector<T>& fsens = fwdSens[d][0]->data();
      for(int k=0; k<nz_.size(); ++k){
	fsens[k] = fseed[nz_[k]];
      }
    }
      
    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      vector<T>& aseed = adjSeed[d][0]->data();
      vector<T>& asens = adjSens[d][0]->data();
      for(int k=0; k<nz_.size(); ++k){
	asens[nz_[k]] += aseed[k];
	aseed[k] = 0;
      }
    }
  }
  
  void GetNonzeros::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    // Get references to the assignment operations and data
    bvec_t *outputd = get_bvec_t(output[0]->data());
    bvec_t *inputd = get_bvec_t(input[0]->data());
    
    // Propate sparsity
    for(int k=0; k<nz_.size(); ++k){
      if(fwd){
	outputd[k] = inputd[nz_[k]];
      } else {
	inputd[nz_[k]] |= outputd[k];
	outputd[k] = 0;
      }
    }
  }

  void GetNonzeros::printPart(std::ostream &stream, int part) const{
    switch(part){
    case 1: stream << nz_; break;
    }
  }

  void GetNonzeros::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    casadi_assert(input.size()==1);

    // Sparsity
    const CRSSparsity &sp = sparsity();
    vector<int> row = sp.getRow();
  
    // Number of derivative directions
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Function evaluation in sparse triplet format
    vector<int> r_row, r_col, r_inz;
  
    // Sensitivity matrices in sparse triplet format for all forward sensitivities
    vector<vector<int> > f_row(nfwd), f_col(nfwd), f_inz(nfwd);
  
    // Sensitivity matrices in sparse triplet format for all adjoint sensitivities
    vector<vector<int> > a_row(nadj), a_col(nadj), a_onz(nadj);
    
    // Output sparsity
    const CRSSparsity &osp = sparsity();
    const vector<int>& ocol = osp.col();
    vector<int> orow = osp.getRow();
    
    // Input sparsity
    const CRSSparsity &isp = dep().sparsity();
    const vector<int>& icol = isp.col();
    vector<int> irow = isp.getRow();
    
    // Find out which matrix elements that we are trying to calculate
    vector<int> el_wanted(nz_.size()); // TODO: Move outside loop
    for(int k=0; k<nz_.size(); ++k){
      el_wanted[k] = orow[k] + ocol[k]*osp.size1();
    }
    
    // We next need to resort the assigns vector with increasing inputs instead of outputs
    // Start by counting the number of inputs corresponding to each nonzero
    vector<int> inz_count(icol.size()+1,0);
    for(vector<int>::const_iterator it=nz_.begin(); it!=nz_.end(); ++it){
      inz_count[*it+1]++;
    }
    
    // Cumsum to get index offset for input nonzero
    for(int i=0; i<icol.size(); ++i){
      inz_count[i+1] += inz_count[i];
    }
    
    // Get the order of assignments
    vector<int> nz_order(nz_.size()); // TODO: Move allocation outside loop
    for(int k=0; k<nz_.size(); ++k){
      // Save the new index
      nz_order[inz_count[nz_[k]]++] = k;
    }
    
    // Find out which matrix elements that we are trying to calculate
    vector<int> el_known(nz_.size()); // TODO: Move allocation outside loop
    for(int k=0; k<nz_.size(); ++k){
      // Get output nonzero
      int inz_k = nz_[nz_order[k]];
      
      // Get element
      el_known[k] = irow[inz_k] + icol[inz_k]*isp.size1();
    }
    
    // Temporary vectors
    vector<int> temp = el_known;
    
    // Evaluate the nondifferentiated function
    if(!output_given){
      
      // Get the matching nonzeros
      input[0]->sparsity().getNZInplace(temp);
      
      // Add to sparsity pattern
      for(int k=0; k<nz_.size(); ++k){
	if(temp[k]!=-1){
	  r_inz.push_back(temp[k]);
	  r_col.push_back(ocol[nz_order[k]]);
	  r_row.push_back(orow[nz_order[k]]);
	}
      }
    }
    
    // Forward sensitivities
    for(int d=0; d<nfwd; ++d){
      
      // Get the matching nonzeros
      copy(el_known.begin(),el_known.end(),temp.begin());
      fwdSeed[d][0]->sparsity().getNZInplace(temp);
      
      // Add to sparsity pattern
      for(int k=0; k<nz_.size(); ++k){
	if(temp[k]!=-1){
	  f_inz[d].push_back(temp[k]);
	  f_col[d].push_back(ocol[nz_order[k]]);
	  f_row[d].push_back(orow[nz_order[k]]);
	}
      }
    }
    
    // Continue of no adjoint sensitivities
    if(nadj>0){
      
      // Resize temp to be able to hold el_wanted
      temp.resize(el_wanted.size());
      
      // Adjoint sensitivities
      for(int d=0; d<nadj; ++d){
	
	// Get the matching nonzeros
	copy(el_wanted.begin(),el_wanted.end(),temp.begin());
	adjSeed[d][0]->sparsity().getNZInplace(temp);
	
	// Add to sparsity pattern
	for(int k=0; k<nz_.size(); ++k){
	  if(temp[k]!=-1){
	    a_onz[d].push_back(temp[k]);
	    a_col[d].push_back(icol[nz_[k]]);
	    a_row[d].push_back(irow[nz_[k]]);
	  }
	}
      }
    }
    
    // Non-differentiated output
    if(!output_given){
      
      // Output sparsity
      const CRSSparsity &osp = sparsity();
      
      // Create a sparsity pattern from vectors
      CRSSparsity r_sp = sp_triplet(osp.size1(),osp.size2(),r_row,r_col);
      
      if(r_inz.size()>0){
	*output[0] = (*input[0])->getGetNonzeros(r_sp,r_inz);
      } else {
	*output[0] = MX::zeros(r_sp);
      }
    }
    
    // Create the forward sensitivity matrices
    for(int d=0; d<nfwd; ++d){
    
      // Output sparsity
      const CRSSparsity &osp = sparsity();
      
      // Create a sparsity pattern from vectors
      CRSSparsity f_sp = sp_triplet(osp.size1(),osp.size2(),f_row[d],f_col[d]);
      if(f_inz[d].size()>0){
	*fwdSens[d][0] = (*fwdSeed[d][0])->getGetNonzeros(f_sp,f_inz[d]);
      } else {
	*fwdSens[d][0] = MX::zeros(f_sp);
      }
    }
  
    // Create the adjoint sensitivity matrices
    vector<int> a_inz;
    for(int d=0; d<nadj; ++d){

      // Input sparsity
      const CRSSparsity &isp = input[0]->sparsity();
      
      // Create a sparsity pattern from vectors
      CRSSparsity a_sp = sp_triplet(isp.size1(),isp.size2(),a_row[d],a_col[d],a_inz,true);

      // Create a mapping matrix
      MX s = MX::create(new Mapping(a_sp));
      s->assign(*adjSeed[d][0],a_onz[d],a_inz,true);
      
      // Save to adjoint sensitivities
      *adjSens[d][0] += s;
      
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
    if(nz_.size()==1){
      // Compact if just a scalar
      stream << "  " << res.front() << "[0]=" << arg.front() << "[" << nz_.front() << "];" << endl;
    } else {
      // Condegen the indices
      int ind = gen.getConstant(nz_,true);
      
      // Codegen the assignments
      stream << "  for(i=0; i<" << nz_.size() << "; ++i) " << res.front() << "[i]=" << arg.front() << "[s" << ind << "[i]];" << endl;
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
