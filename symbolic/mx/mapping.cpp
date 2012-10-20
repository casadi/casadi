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

#include "mapping.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../fx/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"

const bool ELIMINATE_NESTED = true;

using namespace std;

namespace CasADi{

Mapping::Mapping(const CRSSparsity& sp){
  setSparsity(sp);
  output_sorted_.resize(sp.size());
}

Mapping* Mapping::clone() const{
  return new Mapping(*this);
}

template<typename T>
void Mapping::evaluateBlock(int iind, int oind, const vector<T>& idata, vector<T>& odata, bool fwd) const{
  // Get references to the assignment operations
  const IOMap& assigns = index_output_sorted_[oind][iind];
  
  if(fwd){
    // Assignment operations
    for(IOMap::const_iterator it=assigns.begin(); it!=assigns.end(); ++it)
      odata[it->second] += idata[it->first];

  } else {
    // Assignment operations
    for(IOMap::const_iterator it=assigns.begin(); it!=assigns.end(); ++it)
      odata[it->first] += idata[it->second];
  }
}

void Mapping::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
}

void Mapping::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
}

template<typename T, typename MatV, typename MatVV>
void Mapping::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){
  // Number of sensitivities
  int nadj = adjSeed.size();
  int nfwd = fwdSens.size();

  // Loop over outputs
  for(int oind=0; oind<output.size(); ++oind){

    // Clear output and forward sensitivities
    if(output[oind]!=0)
      output[oind]->setZero();
    for(int d=0; d<nfwd; ++d)
      if(fwdSens[d][oind]!=0)
        fwdSens[d][oind]->setZero();
    
    // Loop over inputs
    for(int iind=0; iind<input.size(); ++iind){
    
      // Nondifferentiated outputs
      if(input[iind]!=0 && output[oind]!=0)
        evaluateBlock(iind,oind,input[iind]->data(),output[oind]->data(),true);

      // Forward sensitivities
      for(int d=0; d<nfwd; ++d){
        if(fwdSeed[d][iind]!=0 && fwdSens[d][oind]!=0)
          evaluateBlock(iind,oind,fwdSeed[d][iind]->data(),fwdSens[d][oind]->data(),true);
      }
      
      // Adjoint sensitivities
      for(int d=0; d<nadj; ++d){
        if(adjSeed[d][oind]!=0 && adjSens[d][iind]!=0)
          evaluateBlock(iind,oind,adjSeed[d][oind]->data(),adjSens[d][iind]->data(),false);
      }
    }
  }
}

void Mapping::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
  
  // Loop over outputs
  for(int oind=0; oind<output.size(); ++oind){

    // Clear output
    if(fwd && output[oind]!=0){
      bvec_t *outputd = get_bvec_t(output[oind]->data());
      fill_n(outputd,output[oind]->size(),0);
    }
    
    // Loop over inputs
    for(int iind=0; iind<input.size(); ++iind){
    
      // Nondifferentiated outputs
      if(input[iind]!=0 && output[oind]!=0){

        // Get references to the assignment operations and data
        const IOMap& assigns = index_output_sorted_[oind][iind];
        bvec_t *outputd = get_bvec_t(output[oind]->data());
        bvec_t *inputd = get_bvec_t(input[iind]->data());
        
        // Propate sparsity
        for(IOMap::const_iterator it=assigns.begin(); it!=assigns.end(); ++it){
          if(fwd){
            outputd[it->second] |= inputd[it->first];
          } else {
            inputd[it->first] |= outputd[it->second];
          }
        }
      }
    }
  }
}

void Mapping::printPart(std::ostream &stream, int part) const{
  if(ndep()==0){
    stream << "sparse(" << size1() << "," << size2() << ")";
  } else if(numel()==1 && size()==1 && ndep()==1 && output_sorted_[0].size()==1){
    if(part==1)
      if(dep(0).numel()>1)
        stream << "[" << output_sorted_[0][0].inz << "]";
  } else {
    if(part==0){
      if(isTranspose()){
	stream << "trans(";
      } else {
	stream << "mapping(";
	if(sparsity().dense())            stream << "dense";
	else if(sparsity().diagonal())    stream << "diagonal";
	else                              stream << "sparse";
	stream << " " << size1() << "-by-" << size2() << " matrix, dependencies: [";
      }
    } else if(part==ndep()){
      if(isTranspose()){
	stream << ")";
      } else {
	stream << "], nonzeros: [";
	for(int k=0; k<output_sorted_.size(); ++k){
	  for(int kk=0; kk<output_sorted_[k].size(); ++kk){
	    if(ndep()>1){
	      stream << output_sorted_[k][kk].inz << "(" << output_sorted_[k][kk].iind << ")";
	    } else {
	      stream << output_sorted_[k][kk].inz;
	    }
	    stream << ",";
	  }
	}
	stream << "])";
      }
    } else {
      stream << ",";
    }
  }
}

void Mapping::assign(const MX& d, const std::vector<int>& inz, bool add){
  assign(d,inz,range(inz.size()),add);
}

void Mapping::assign(const MX& d, const std::vector<int>& inz, const std::vector<int>& onz, bool add){
  // Quick return if no elements
  if(inz.empty()) return;
  
  casadi_assert(!d.isNull());
  
  if(ELIMINATE_NESTED && d->getOp()==OP_MAPPING){ // Move this logic to init!
    // Clear the existing element if we are not adding
    if(!add){
      for(int k=0; k<onz.size(); ++k){
        output_sorted_[onz[k]].clear();
      }
    }
    
    // Eliminate if a mapping node
    const Mapping* dnode = static_cast<const Mapping*>(d.get());
    vector<MX> d2 = dnode->dep_;
    
    // Split the vector according to dependency index
    vector<vector<int> > inz2(d2.size()), onz2(d2.size());
    for(int k=0; k<inz.size(); ++k){
      
      // Get the sum
      const std::vector<OutputNZ>& sum = dnode->output_sorted_[inz[k]];
      
      // Add the elements in the sum
      for(std::vector<OutputNZ>::const_iterator it2=sum.begin(); it2!=sum.end(); ++it2){
        inz2[it2->iind].push_back(it2->inz);
        onz2[it2->iind].push_back(onz[k]);
      }
    }
    
    // Call the function recursively
    for(int i=0; i<d2.size(); ++i){
      assign(d2[i],inz2[i],onz2[i],true);
    }
  } else {
    // Add the node if it is not already a dependency
    std::map<const MXNode*, int>::const_iterator it = depmap_.find(static_cast<const MXNode*>(d.get()));
    int depind;
    if(it==depmap_.end()){
      depind = addDependency(d);
      depmap_[static_cast<const MXNode*>(d.get())] = depind;
    } else {
      depind = it->second;
    }
    
    // Save the mapping
    for(int k=0; k<inz.size(); ++k){
      OutputNZ new_el = {inz[k],depind};
      if(!add) output_sorted_[onz[k]].clear();
      output_sorted_[onz[k]].push_back(new_el);
    }
  }
}

void Mapping::init(){
  // Call init of the base class
  MXNode::init();
  
  // Clear the runtime
  index_output_sorted_.resize(1);
  index_output_sorted_[0].resize(ndep());
  for(int iind=0; iind<index_output_sorted_[0].size(); ++iind){
    index_output_sorted_[0][iind].clear();
  }
  
  // For all the outputs
  for(int onz=0; onz<output_sorted_.size(); ++onz){

    // Get the sum
    const std::vector<OutputNZ>& sum = output_sorted_[onz];
    
    // For all elements in the sum
    for(std::vector<OutputNZ>::const_iterator it=sum.begin(); it!=sum.end(); ++it){
      
      // Add the element to the runtime
      index_output_sorted_[0][it->iind].push_back(pair<int,int>(it->inz,onz));
    }
  }
}

void Mapping::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  // Sparsity
  const CRSSparsity &sp = sparsity();
  const vector<int>& rowind = sp.rowind();
  const vector<int>& col = sp.col();
  vector<int> row = sp.getRow();

  // Dimensions
  int d1=sp.size1(), d2=sp.size2();
  
  // Number of derivative directions
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();

  // Number of inputs and outputs
  int niind = input.size();
  int noind = output.size();
  
  // Function evaluation in sparse triplet format
  vector<vector<int> > r_row, r_col, r_inz, r_offset;
  if(!output_given){
    r_row.resize(noind);
    r_col.resize(noind);
    r_inz.resize(noind);
    r_offset.resize(noind,vector<int>(1,0));
  }
  
  // Sensitivity matrices in sparse triplet format for all forward sensitivities
  vector<vector<vector<int> > > f_row(nfwd), f_col(nfwd), f_inz(nfwd), f_offset(nfwd);
  for(int d=0; d<nfwd; ++d){
    f_row[d].resize(noind);
    f_col[d].resize(noind);
    f_inz[d].resize(noind);
    f_offset[d].resize(noind,vector<int>(1,0));
  }
  
  // Sensitivity matrices in sparse triplet format for all adjoint sensitivities
  vector<vector<vector<int> > > a_row(nadj), a_col(nadj), a_onz(nadj), a_offset(nadj);
  for(int d=0; d<nadj; ++d){
    a_row[d].resize(niind);
    a_col[d].resize(niind);
    a_onz[d].resize(niind);
    a_offset[d].resize(niind,vector<int>(1,0));
  }
    
  // For all outputs
  for(int oind=0; oind<output.size(); ++oind){
    if(output[oind]==0) continue; // Skip if output doesn't exist
    
    // Output sparsity
    const CRSSparsity &osp = sparsity(oind);
    const vector<int>& ocol = osp.col();
    vector<int> orow = osp.getRow();
    
    // For all inputs
    for(int iind=0; iind<input.size(); ++iind){
      
      // Skip if input doesn't exist
      if(input[iind]==0) continue;

      // Input sparsity
      const CRSSparsity &isp = dep(iind).sparsity();
      const vector<int>& icol = isp.col();
      vector<int> irow = isp.getRow();

      // Get references to the assignment operations
      const IOMap& assigns = index_output_sorted_[oind][iind];

      // Find out which matrix elements that we are trying to calculate
      vector<int> el_wanted(assigns.size()); // TODO: Move outside loop
      for(int k=0; k<assigns.size(); ++k){
        el_wanted[k] = orow[assigns[k].second] + ocol[assigns[k].second]*osp.size1();
      }

      // We next need to resort the assigns vector with increasing inputs instead of outputs
      // Start by counting the number of inputs corresponding to each nonzero
      vector<int> inz_count(icol.size()+1,0);
      for(IOMap::const_iterator it=assigns.begin(); it!=assigns.end(); ++it){
        inz_count[it->first+1]++;
      }
      
      // Cumsum to get index offset for input nonzero
      for(int i=0; i<icol.size(); ++i){
        inz_count[i+1] += inz_count[i];
      }
      
      // Get the order of assignments
      vector<int> assigns_order(assigns.size()); // TODO: Move allocation outside loop
      for(int k=0; k<assigns.size(); ++k){
        // Save the new index
        assigns_order[inz_count[assigns[k].first]++] = k;
      }

      // Find out which matrix elements that we are trying to calculate
      vector<int> el_known(assigns.size()); // TODO: Move allocation outside loop
      for(int k=0; k<assigns.size(); ++k){
        // Get output nonzero
        int inz_k = assigns[assigns_order[k]].first;
        
        // Get element
        el_known[k] = irow[inz_k] + icol[inz_k]*isp.size1();
      }

      // Temporary vectors
      vector<int> temp = el_known;
      
      // Evaluate the nondifferentiated function
      if(!output_given){
        
        // Get the matching nonzeros
        input[iind]->sparsity().getNZInplace(temp);

        // Add to sparsity pattern
        for(int k=0; k<assigns.size(); ++k){
          if(temp[k]!=-1){
            r_inz[oind].push_back(temp[k]);
            r_col[oind].push_back(ocol[assigns[assigns_order[k]].second]);
            r_row[oind].push_back(orow[assigns[assigns_order[k]].second]);
          }
        }
        r_offset[oind].push_back(r_inz[oind].size());
      }
      
      // Forward sensitivities
      for(int d=0; d<nfwd; ++d){
        
        // Get the matching nonzeros
        copy(el_known.begin(),el_known.end(),temp.begin());
        fwdSeed[d][iind]->sparsity().getNZInplace(temp);

        // Add to sparsity pattern
        for(int k=0; k<assigns.size(); ++k){
          if(temp[k]!=-1){
            f_inz[d][oind].push_back(temp[k]);
            f_col[d][oind].push_back(ocol[assigns[assigns_order[k]].second]);
            f_row[d][oind].push_back(orow[assigns[assigns_order[k]].second]);
          }
        }
        f_offset[d][oind].push_back(f_inz[d][oind].size());
      }
      
      // Continue of no adjoint sensitivities
      if(nadj==0) continue;
      
      // Resize temp to be able to hold el_wanted
      temp.resize(el_wanted.size());
      
      // Adjoint sensitivities
      for(int d=0; d<nadj; ++d){
        
        // Get the matching nonzeros
        copy(el_wanted.begin(),el_wanted.end(),temp.begin());
        adjSeed[d][oind]->sparsity().getNZInplace(temp);

        // Add to sparsity pattern
        for(int k=0; k<assigns.size(); ++k){
          if(temp[k]!=-1){
            a_onz[d][iind].push_back(temp[k]);
            a_col[d][iind].push_back(icol[assigns[k].first]);
            a_row[d][iind].push_back(irow[assigns[k].first]);
          }
        }
        a_offset[d][iind].push_back(a_onz[d][iind].size());
      }
    }
  }
  
  // Non-differentiated output
  vector<int> r_onz, inz_local, onz_local;
  if(!output_given){
    // For all outputs
    for(int oind=0; oind<output.size(); ++oind){
      if(output[oind]==0) continue; // Skip if output doesn't exist

      // Output sparsity
      const CRSSparsity &osp = sparsity(oind);

      // Create a sparsity pattern from vectors
      CRSSparsity r_sp = sp_triplet(osp.size1(),osp.size2(),r_row[oind],r_col[oind],r_onz,true);
      
      // Create a mapping matrix
      *output[oind] = MX::create(new Mapping(r_sp));
      
      // Add all the dependencies
      for(int iind=0; iind<input.size(); ++iind){
        if(input[iind]==0) continue; // Skip if input doesn't exist
        
        // Get the elements corresponding to the input index
        int iind0 = r_offset[oind][iind], iind1 = r_offset[oind][iind+1];

        // Get the input nonzeros
        inz_local.resize(iind1-iind0);
        copy(r_inz[oind].begin()+iind0,r_inz[oind].begin()+iind1,inz_local.begin());
        
        // Get the output nonzeros
        onz_local.resize(iind1-iind0);
        copy(r_onz.begin()+iind0,r_onz.begin()+iind1,onz_local.begin());
        
        // Save to mapping
        (*output[oind])->assign(*input[iind],inz_local,onz_local,true);
      }
    }
  }
    
  // Create the forward sensitivity matrices
  vector<int>& f_onz = r_onz; // reuse
  for(int d=0; d<nfwd; ++d){
    
    // For all outputs
    for(int oind=0; oind<output.size(); ++oind){
      if(output[oind]==0) continue; // Skip if output doesn't exist

      // Output sparsity
      const CRSSparsity &osp = sparsity(oind);

      // Create a sparsity pattern from vectors
      CRSSparsity f_sp = sp_triplet(osp.size1(),osp.size2(),f_row[d][oind],f_col[d][oind],f_onz,true);
      
      // Create a mapping matrix
      *fwdSens[d][oind] = MX::create(new Mapping(f_sp));
      
      // Add all the dependencies
      for(int iind=0; iind<input.size(); ++iind){
        if(input[iind]==0) continue; // Skip if input doesn't exist
        
        // Get the elements corresponding to the input index
        int iind0 = f_offset[d][oind][iind], iind1 = f_offset[d][oind][iind+1];

        // Get the input nonzeros
        inz_local.resize(iind1-iind0);
        copy(f_inz[d][oind].begin()+iind0,f_inz[d][oind].begin()+iind1,inz_local.begin());
        
        // Get the output nonzeros
        onz_local.resize(iind1-iind0);
        copy(f_onz.begin()+iind0,f_onz.begin()+iind1,onz_local.begin());
        
        // Save to mapping
        (*fwdSens[d][oind])->assign(*fwdSeed[d][iind],inz_local,onz_local,true);
      }
    }
  }
  
  // Create the adjoint sensitivity matrices
  vector<int>& a_inz=f_onz; // reuse
  for(int d=0; d<nadj; ++d){

    // For all inputs
    for(int iind=0; iind<input.size(); ++iind){
      if(input[iind]==0) continue; // Skip if input doesn't exist
        
      // Input sparsity
      const CRSSparsity &isp = input[iind]->sparsity();
      
      // Create a sparsity pattern from vectors
      CRSSparsity a_sp = sp_triplet(isp.size1(),isp.size2(),a_row[d][iind],a_col[d][iind],a_inz,true);
      
      // Create a mapping matrix
      MX s = MX::create(new Mapping(a_sp));
     
      // Add all dependencies
      for(int oind=0; oind<output.size(); ++oind){
        if(output[oind]==0) continue; // Skip if output doesn't exist
          
        // Get the elements corresponding to the input index
        int oind0 = a_offset[d][iind][oind], oind1 = a_offset[d][iind][oind+1];

        // Get the input nonzeros
        inz_local.resize(oind1-oind0);
        copy(a_inz.begin()+oind0,a_inz.begin()+oind1,inz_local.begin());
        
        // Get the output nonzeros
        onz_local.resize(oind1-oind0);
        copy(a_onz[d][iind].begin()+oind0,a_onz[d][iind].begin()+oind1,onz_local.begin());
        
        // Save to mapping
        s->assign(*adjSeed[d][oind],onz_local,inz_local,true);
      }
      
      // Save to adjoint sensitivities
      *adjSens[d][iind] += s;
    }
  }
}

Matrix<int> Mapping::mapping(int iind) const {
  casadi_assert_message(iind < ndep(),"Mapping::mapping(int): first argument (" << iind << ") must be smaller than ndep (" << ndep() << ").");
  // TODO: make this efficient
  std::vector< int > row;
  std::vector< int > col;
  sparsity().getSparsity(row,col);
  Matrix<int> ret(size1(),size2());
  for (int k=0;k<output_sorted_.size();++k) { // Loop over output non-zeros
    for (int i=0;i<output_sorted_[k].size(); ++i) { // Loop over elements to be summed
      const OutputNZ &el = output_sorted_[k][i];
      if (el.iind==iind) ret(row[k],col[k]) = el.inz;
    }
  }
  return ret;
}

std::vector<int> Mapping::getDepInd() const {
  // TODO: make this efficient
  std::vector<int> ret(size());
  for (int k=0;k<output_sorted_.size();++k) { // Loop over output non-zeros
    for (int i=0;i<output_sorted_[k].size(); ++i) { // Loop over elements to be summed
      const OutputNZ &el = output_sorted_[k][i];
      ret[k] = el.iind;
    }
  }
  return ret;
}

bool Mapping::isIdentity() const{
  // Make sure that there is at least one dependency
  if(ndep()<1) return false;
  
  // Check sparsity
  if(!(sparsity() == dep(0).sparsity()))
    return false;
      
  // Check if the nonzeros follow in increasing order
  for(int k=0; k<output_sorted_.size(); ++k){
    if(output_sorted_[k].size()!=1) return false;
    const OutputNZ &e = output_sorted_[k].front();
    if(e.inz != k || e.iind !=0) return false;
  }
    
  // True if reached this point
  return true;
}

bool Mapping::isTranspose() const{
  // Make sure that there is at least one dependency
  if(ndep()!=1) return false;
  
  // Get the dependency
  const MX& d = dep(0);
  
  // First check if sparsity patterns are transposes of each other
  if(!d.sparsity().isTranspose(sparsity()))
    return false;
  
  // Get all the elements of the transpose
  vector<int> d_elements = d.sparsity().getElements(false);
  
  // Sparsity pattern of the transpose
  int sz1 = size1();
  int sz2 = size2();
  const vector<int>& rowind = sparsity().rowind();
  const vector<int>& col = sparsity().col();
  
  // Loop over the rows of the matrix
  for(int i=0; i<sz1; ++i){
    
    // Loop over the nonzeros of the row
    for(int k=rowind[i]; k<rowind[i+1]; ++k){
      
      // Get the column
      int j=col[k];
      
      // Make sure that there is a single output
      if(output_sorted_[k].size()!=1) return false;
      const OutputNZ &e = output_sorted_[k].front();
      if(e.iind !=0) return false;
      
      // Get the element
      int el = d_elements[e.inz];
      
      // Make sure that it is the transpose
      if(el!=j+i*sz2) return false;
    }
  }
  
  // True if reached this point
  return true;
}

} // namespace CasADi
