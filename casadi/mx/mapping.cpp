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

void Mapping::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
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
  casadi_assert_message(output_given,"not implemented");
  
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

  // Sensitivity matrices in sparse triplet format for all forward sensitivities
  vector<vector<vector<int> > > f_row(nfwd), f_col(nfwd), f_inz(nfwd), f_offset(nfwd);
  for(int d=0; d<nfwd; ++d){
    int noind = output.size();
    f_row[d].resize(noind);
    f_col[d].resize(noind);
    f_inz[d].resize(noind);
    f_offset[d].resize(noind,vector<int>(1,0));
  }
  
  // Sensitivity matrices in sparse triplet format for all adjoint sensitivities
  vector<vector<vector<int> > > a_row(nadj), a_col(nadj), a_onz(nadj), a_offset(nadj);
  for(int d=0; d<nadj; ++d){
    int niind = input.size();
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

      // Find out which matrix elements that we are trying to access
      vector<int> el_known(assigns.size());
      for(int k=0; k<assigns.size(); ++k){
        el_known[k] = irow[assigns[k].first] + icol[assigns[k].first]*isp.size1();
      }

      // Temporary vectors
      vector<int> temp = el_known;
      
      // Evaluate the nondifferentiated function
      if(!output_given){
        casadi_assert(0);
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
            f_col[d][oind].push_back(ocol[assigns[k].second]);
            f_row[d][oind].push_back(orow[assigns[k].second]);
          }
        }
        f_offset[d][oind].push_back(f_inz[d][oind].size());
      }
      
      // Continue of no adjoint sensitivities
      if(nfwd==0) continue;
      
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
      vector<int> assigns_order(assigns.size()); // // TODO: Move outside loop
      for(int k=0; k<assigns.size(); ++k){
        // Save the new index
        assigns_order[inz_count[assigns[k].first]++] = k;
      }

      // Find out which matrix elements that we are trying to calculate
      vector<int> el_wanted(assigns.size()); // TODO: Move outside loop
      for(int k=0; k<assigns.size(); ++k){
        // Get output nonzero
        int onz_k = assigns[el_wanted[k]].second;
        
        // Get element
        el_wanted[k] = orow[onz_k] + ocol[onz_k]*osp.size1();
      }
      
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
            a_col[d][iind].push_back(icol[assigns[assigns_order[k]].second]);
            a_row[d][iind].push_back(irow[assigns[assigns_order[k]].second]);
          }
        }
        a_offset[d][iind].push_back(a_onz[d][iind].size());
      }
    }
  }
    
  // Create the forward sensitivity matrices
  vector<int> f_onz, inz_local, onz_local;
  for(int d=0; d<nfwd; ++d){
    
    // For all outputs
    for(int oind=0; oind<output.size(); ++oind){
      if(output[oind]==0) continue; // Skip if output doesn't exist

      // Output sparsity
      const CRSSparsity &osp = sparsity(oind);

      // Create a sparsity pattern from vectors
      f_onz.clear();
      CRSSparsity f_sp = sp_triplet(osp.size1(),osp.size2(),f_row[d][oind],f_col[d][oind],f_onz,false,true);
      
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
        (*fwdSens[d][oind])->addDependency(*fwdSeed[d][iind],inz_local,onz_local);
      }
    }
  }
  
  // Create the adjoint sensitivity matrices
  vector<int>& a_inz=f_onz; // reuse
  for(int d=0; d<nadj; ++d){
    continue;
  
    // For all inputs
    for(int iind=0; iind<input.size(); ++iind){
      if(input[iind]==0) continue; // Skip if input doesn't exist
        
      // Input sparsity
      const CRSSparsity &isp = dep(iind).sparsity();
        
      // Create a sparsity pattern from vectors
      a_inz.clear();
      CRSSparsity a_sp = sp_triplet(isp.size1(),isp.size2(),a_row[d][iind],a_col[d][iind],a_inz,false,true);
        
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
        s->addDependency(*adjSeed[d][iind],onz_local,inz_local);
      }
      
      // Save to adjoint sensitivities
      //*adjSens[d][iind] += s; // FIXME: DOES NOT YET WORK, DEBUGGING NEEDED
    }
  }
  
  // NOTE: Everything below should be deleted
  
//  return;
  //casadi_assert(nadj==0);
  
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
    }
  }
}

} // namespace CasADi
