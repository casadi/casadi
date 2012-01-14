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
#include "inverse_mapping.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../fx/sx_function.hpp"

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


void Mapping::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
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

bool Mapping::isReady() const{
  return true;
}
    
void Mapping::printPart(std::ostream &stream, int part) const{
  casadi_assert(isReady());

  if(ndep()==0){
    stream << "sparse(" << size1() << "," << size2() << ")";
  } else if(numel()==1 && size()==1 && ndep()==1 && output_sorted_[0].size()==1){
    if(part==1)
      if(dep(0).numel()>1)
        stream << "[" << output_sorted_[0][0].inz << "]";
  } else {
    if(part==0){
      stream << "mapping(";
      if(sparsity().dense())            stream << "dense";
      else if(sparsity().diagonal())    stream << "diagonal";
      else                              stream << "sparse";
      stream << " " << size1() << "-by-" << size2() << " matrix, dependencies: [";
    } else if(part==ndep()){
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
    } else {
      stream << ",";
    }
  }
}

void Mapping::assign(const MX& d, const IOMap& iomap){
  casadi_assert(!d.isNull());
  //const std::vector<int>& nzind_ = nzmap_.data();
  
  // Quick return if no elements
  if(iomap.empty()) return;
  
  if(ELIMINATE_NESTED && d->isMapping()){ // Move this logic to init!
    // Eliminate if a mapping node
    const Mapping* dnode = static_cast<const Mapping*>(d.get());
    vector<MX> d2 = dnode->dep_;
    vector<IOMap> iomap2(d2.size());
    for(IOMap::const_iterator it=iomap.begin(); it!=iomap.end(); it++){
      // Get the sum
      const std::vector<OutputNZ>& sum = dnode->output_sorted_[it->first];
      
      // Add the elements in the sum
      for(std::vector<OutputNZ>::const_iterator it2=sum.begin(); it2!=sum.end(); ++it2){
        iomap2[it2->iind].push_back(pair<int,int>(it2->inz, it->second));
      }
    }
    
    // Call the function recursively
    for(int i=0; i<d2.size(); ++i){
      assign(d2[i],iomap2[i]);
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
    for(IOMap::const_iterator it=iomap.begin(); it!=iomap.end(); ++it){
      OutputNZ new_el = {it->first,depind};
      output_sorted_[it->second].clear(); // FIXME
      output_sorted_[it->second].push_back(new_el);
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
  casadi_assert(isReady());
  
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

  // For all outputs
  for(int oind=0; oind<output.size(); ++oind){

    // Skip if output doesn't exist
    if(output[oind]==0) continue;
    
    // Output sparsity
    const CRSSparsity &osp = sparsity(oind);
    const vector<int>& ocol = osp.col();
    vector<int> orow = osp.getRow();

    // Zero output
    MX zero_output = MX::sparse(d1,d2);
    
    // Initialize output to zero
    if(!output_given){
      *output[oind] = zero_output;
    }
    
    // Initialize forward derivatives to zero
    for(int d=0; d<nfwd; ++d){
      *fwdSens[d][oind] = zero_output;
    }
    
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
      
      // Step 1: Find out which matrix elements that we are trying to calculate
      vector<int> row_wanted(assigns.size()), col_wanted(assigns.size()), el_wanted(assigns.size());
      for(int k=0; k<assigns.size(); ++k){
        row_wanted[k] = orow[assigns[k].second];
        col_wanted[k] = ocol[assigns[k].second];
        el_wanted[k] = row_wanted[k] + col_wanted[k]*osp.size1();
      }
      
      // Step 2: Find out which matrix elements that we are trying to access
      vector<int> row_known(assigns.size()), col_known(assigns.size()), el_known(assigns.size());
      for(int k=0; k<assigns.size(); ++k){
        row_known[k] = irow[assigns[k].first];
        col_known[k] = icol[assigns[k].first];
        el_known[k] = row_known[k] + col_known[k]*isp.size1();
      }
      
      // Find the order of the inputs in the inputs
      vector<int> iorder(isp.size(),-1);
      for(int k=0; k<assigns.size(); ++k) 
        iorder[assigns[k].first] = k;
      int iorderk=0;
      for(vector<int>::const_iterator it=iorder.begin(); it!=iorder.end(); ++it){
        if(*it>=0) iorder[iorderk++] = *it;
      }
      iorder.resize(iorderk);
      
      // Evaluate the nondifferentiated function
      if(!output_given){
        // Get references to the data
        MX& ip = *output[oind];
        MX& op = *input[iind];
      }
      
      // The elements in the block
      vector<int> elems(assigns.size());
      
      // Forward sensitivities
      for(int d=0; d<nfwd; ++d){
        
      }
      
      // Adjoint sensitivities
      for(int d=0; d<nadj; ++d){
        
      }
    }
  }
  
  casadi_assert_message(output_given,"not implemented");
  casadi_assert(isReady());
  if(nadj>0){
    // Number of inputs
    int n = input.size();
    
    // All the input and output nonzeros and indices
    vector<vector<int> > input_el(n), input_el_sorted(n);
    vector<vector<int> > output_el(n), output_el_sorted(n);
    vector<int> input_nz;
    
    // Temporary vector
    vector<int> temp;
    
    // For all inputs
    for(int iind=0; iind<n; ++iind){
      
      // Find input/output pairs
      input_nz.clear();
      
      // Loop over rows
      for(int i=0; i<d1; ++i){
        // Loop over nonzeros
        for(int k=rowind[i]; k<rowind[i+1]; ++k){
          // Get column
          int j=col[k];
          
          // Add if index matches
          if(output_sorted_[k][0].iind==iind){
            input_nz.push_back(output_sorted_[k][0].inz);
            output_el[iind].push_back(j + i*d2);
          }
        }
      }
      
      // At this point we shall construct input_el and at the same time sort according to it
      input_el[iind].resize(input_nz.size());
      input_el_sorted[iind].resize(input_nz.size());
      output_el_sorted[iind].resize(input_nz.size());
      
      // Get input sparsity
      const CRSSparsity& sp_in = input[iind]->sparsity();
      const vector<int>& c_in = sp_in.col();
      const vector<int>& r_in = sp_in.rowind();
      int d1_in = sp_in.size1();
      int d2_in = sp_in.size2();
      
      // Add extra temporary elements if necessary
      temp.resize(max(temp.size(),c_in.size()),-1);
      
      // Mark inputs
      for(int k=0; k<input_nz.size(); ++k){
        temp[input_nz[k]] = k;
      }
      
      // Loop over rows
      int el = 0;
      for(int i=0; i<d1_in; ++i){
        // Loop over nonzeros
        for(int k=r_in[i]; k<r_in[i+1]; ++k){
          // Get column
          int j=c_in[k];
          
          // Save to vector if nonzero requested
          if(temp[k]>=0){
            input_el[iind][temp[k]] = input_el_sorted[iind][el] = j + i*d2_in;
            output_el_sorted[iind][temp[k]] = output_el[iind][el];
            el++;
          }
        }
      }
      casadi_assert(el==output_el[iind].size());
      
      // Unmark inputs
      for(int k=0; k<input_nz.size(); ++k){
        temp[input_nz[k]] = -1;
      }
      
      // At this point we have two vector pairs, output_el and input_el containing the elements of the mapping matrix
      // sorted according to output and output_el_sorted and input_el_sorted containing the same information but 
      // sorted according to input
  /*    cout << "output_el = " << output_el << endl;
      cout << "input_el = " << input_el << endl;
      cout << "output_el_sorted = " << output_el_sorted << endl;
      cout << "input_el_sorted = " << input_el_sorted << endl;
      cout << endl;*/
    }
    
    // Now for all forward directions
    for(int d=0; d<nfwd; ++d){
      
    }
    
    // Now for all adjoint directions NOTE: This is a quick-hack and need to be fixed for larger systems
    for(int d=0; d<nadj; ++d){
      for(int iind=0; iind<n; ++iind){
        for(int k=0; k<input_el[iind].size(); ++k){
          // Which element to get
          int el = output_el[iind][k];
          int o1 = el / adjSeed[d][0]->size2();
          int o2 = el % adjSeed[d][0]->size2();
          
          // Get the seed
          MX seed = (*adjSeed[d][0])(o1,o2);
          if(!isZero(seed)){
            int el = input_el[iind][k];
            int i1 = el / adjSens[d][iind]->size2();
            int i2 = el % adjSens[d][iind]->size2();
            (*adjSens[d][iind])(i1,i2) += seed;
          }
        }
      }
      
/*      // Sparsity of the adjoint sensitivities
      vector<CRSSparsity> sp_adjsens(n);
      for(int iind=0; iind<n; ++iind){
        sp_adjsens[iind] = CRSSparsity(0,input[iind]->size1());
      }
      
      // Loop over the rows of the adjoint seed
      for(int i=0; i<d1; ++i){
        
        
        
      }*/
      
      
        
// /*        // Sparsity of the adjoint seed
//         const CRSSparsity& sp_a = adjSeed[d][iind]->sparsity();
//         const vector<int>& c_a = sp_a.col();
//         const vector<int>& r_a = sp_a.rowind();
//         int d1_a = sp_a.size1();
//         int d2_a = sp_a.size2();*/

        
        
        
        
//        sp_dep
//      }
      
    }
    return;
  }
  
  
  
  
  
  
  
  
  
  
  
  // Quick return if no inputs
  if(nfwd>0 && input.empty()){
    for(int oind=0; oind<output.size(); ++oind){
      if(fwdSens[0][oind]!=0){
        *fwdSens[0][oind] = MX::sparse(d1,d2);
        for(int d=0; d<nfwd; ++d){
          *fwdSens[d][oind] = *fwdSens[0][oind];
        }
      }
    }
    return;
  }
  
  // Mapping each nonzero to a nonzero of the forward sensitivity matrix, of -1 if none
  vector<int> nzind_sens(output_sorted_.size()); // BUG? Initialize to -1?
  
  // Mapping from input nonzero index to forward sensitivity nonzero index, or -1
  vector<int> fsens_ind;

  // Mapping for a specific input
  vector<int> &nz = fsens_ind; // reuse memory
  vector<int> nzd;

  // For all forward directions
  for(int d=0; d<nfwd; ++d){

    // Number of nonzeros of the output
    int nnz=0;
    
    // For all inputs
    for(int iind=0; iind<input.size(); ++iind){
      // Mapping from input nonzero index to forward sensitivity nonzero index, or -1
      fsens_ind.resize(input[iind]->size());
      fill(fsens_ind.begin(),fsens_ind.end(),-1);
      
      // Get sparsity of the input and forward sensitivity
      int id1=input[iind]->size1(), id2=input[iind]->size2();
      const vector<int>& rowind_i = input[iind]->sparsity().rowind();
      const vector<int>& rowind_f = fwdSeed[d][iind]->sparsity().rowind();
      const vector<int>& col_i = input[iind]->sparsity().col();
      const vector<int>& col_f = fwdSeed[d][iind]->sparsity().col();
      
      // Loop over rows of input and forward sensitivities
      for(int i=0; i<id1; ++i){
        
        // Nonzero of the forward sensitivity
        int el_f = rowind_f[i];

        // Column of the forward sensitivity (-1 if no element)
        int j_f = el_f==rowind_f[i+1] ? -1 : col_f[el_f];
        
        // Loop over nonzeros of the input
        for(int el=rowind_i[i]; el<rowind_i[i+1]; ++el){
          
          // Column of the input
          int j=col_i[el];
          
          // Continue to the same entry in the forward sensitivity matrix
          while(el_f < rowind_f[i+1] && j_f<j){
            el_f++;
            j_f = col_f[el_f];
          }
          
          // Add element to temp vector or -1 of no corresponding entry
          fsens_ind[el] = j==j_f ? el_f : -1;
        }
      }

      // Update nonzero vector
      for(int el=0; el<nzind_sens.size(); ++el){
        // Check if dependency index match
        if(output_sorted_[el][0].iind==iind){
          
          // Point nzind_sens to the nonzero index of the sensitivity
          nzind_sens[el] = fsens_ind[output_sorted_[el][0].inz];
          
          // Count the number of nonzeros
          nnz += int(nzind_sens[el]>=0);
        }
      }
    }
    
    // Sparsity of the return matrix
    CRSSparsity sp_fsens(d1,d2);
    sp_fsens.reserve(nnz,d1);
    
    // Get references to the vectors
    vector<int>& col_fsens = sp_fsens.colRef();
    vector<int>& rowind_fsens = sp_fsens.rowindRef();
    
    // Loop over rows of the resulting matrix
    for(int i=0; i<d1; ++i){
      
      // Loop over the nonzero elements of the resulting matrix
      for(int el=rowind[i]; el<rowind[i+1]; ++el){
      
        // If the corresponding entry exists in the sensitivity matrix
        if(nzind_sens[el]>=0){
          
          // Get column
          int j=col[el];
          
          // Add nonzero
          col_fsens.push_back(j);
        }
      }
      
      // Save upper bound on nonzero index of the row
      rowind_fsens[i+1] = col_fsens.size();
    }
    
    // Return matrix
    *fwdSens[d][0] = MX::create(new Mapping(sp_fsens));
    
    // Add the dependencies
    for(int dp=0; dp<input.size(); ++dp){
      
      // Get the local nonzeros
      nz.clear();
      nzd.clear();
      int el=0;
      for(int k=0; k<nzind_sens.size(); ++k){
        // If a nonzero
        if(nzind_sens[k]>=0){
        
          // If dependency matches
          if(output_sorted_[k][0].iind==dp){
            nz.push_back(el);
            nzd.push_back(nzind_sens[k]);
          }
          
          // Next nonzero
          el++;
        }
      }
      
      // Save to return matrix
      (*fwdSens[d][0])->addDependency(*fwdSeed[d][dp],nzd,nz);
    }
  }
}

void Mapping::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  for(int k=0; k<size(); ++k){
    (*output[0])[k] = (*input[output_sorted_[k][0].iind])[output_sorted_[k][0].inz];
  }
}

} // namespace CasADi
