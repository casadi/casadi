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

#include "evaluation.hpp"
#include "../fx/fx_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../fx/x_function.hpp"
#include "jacobian_reference.hpp"

using namespace std;

namespace CasADi{

Evaluation::Evaluation(const FX& fcn, const vector<MX>& dep) : fcn_(fcn) {
  // Argument checking
  if (dep.size()!=fcn.getNumInputs()) {
    std::stringstream s;
    s << "Evaluation::Evaluation: number of passed-in dependencies (" << dep.size() << ") should match number of inputs of function (" << fcn.getNumInputs() << ").";
    throw CasadiException(s.str());
  }
  // Assumes initialised
  for (int i=0;i<fcn.getNumInputs();i++) {
     if (dep[i].isNull())
       continue;
      if (dep[i].size1()!=fcn.input(i).size1() || dep[i].size2()!=fcn.input(i).size2()) {
        std::stringstream s;
        s << "Evaluation::shapes of passed-in dependencies should match shapes of inputs of function." << std::endl;
        s << "Input argument " << i << " has shape (" << fcn.input(i).size1() << "," << fcn.input(i).size2() << ") while a shape (" << dep[i].size1() << "," << dep[i].size2() << ") was supplied.";
        throw CasadiException(s.str());
      }     
   }
  setDependencies(dep);
  setSparsity(CRSSparsity(1,1,true));
}

Evaluation* Evaluation::clone() const{
  return new Evaluation(*this);
}

void Evaluation::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << fcn_ << ".call(" << args << ")";
}

void Evaluation::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  
  // Pass the input and forward seeds to the function
  for(int i=0; i<input.size(); ++i){
    if(input[i] != 0 && input[i]->size() !=0 ){
      fcn_.setInput(input[i]->data(),i);
      for(int d=0; d<nfwd; ++d){
        fcn_.setFwdSeed(fwdSeed[d][i]->data(),i,d);
      }
    }
  }
  
  // Pass the adjoint seed to the function
  for(int i=0; i<output.size(); ++i){
    for(int d=0; d<nadj; ++d){
      if(adjSeed[d][0]!=0 && adjSeed[d][0]->size() != 0){
        fcn_.setAdjSeed(adjSeed[d][0]->data(),i,d);
      }
    }
  }

  // Evaluate
  fcn_.evaluate(nfwd, nadj);
  
  // Get the outputs and forward sensitivities
  for(int i=0; i<output.size(); ++i){
    if(output[i] != 0 && output[i]->size() !=0 ){
      fcn_.getOutput(output[i]->data(),i);
      for(int d=0; d<nfwd; ++d){
        fcn_.getFwdSens(fwdSens[d][i]->data(),i,d);
      }
    }
  }
  
  // Get the adjoint sensitivities
  for(int i=0; i<input.size(); ++i){
    for(int d=0; d<nadj; ++d){
      if(adjSens[d][i] != 0 && adjSens[d][i]->size() != 0){
        const vector<double>& asens = fcn_.adjSens(i,d).data();
        for(int j=0; j<asens.size(); ++j)
          adjSens[d][i]->data()[j] += asens[j];
      }
    }
  }
}

    
const CRSSparsity& Evaluation::sparsity(int oind){
  return fcn_.output(oind).sparsity();
}

FX& Evaluation::getFunction(){ 
  return fcn_;
}

void Evaluation::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  // Make sure that the function is an X-function
  XFunction fcn = shared_cast<XFunction>(fcn_);
  casadi_assert_message(!fcn.isNull(),"Function not an SXFunction or MXFunction");
  vector<SXMatrix> arg(input.size());
  for(int i=0; i<arg.size(); ++i){
    if(input[i]!=0)
      arg[i] = *input[i];
  }
  
  std::vector<SXMatrix> res = fcn.eval(arg);
  
  for(int i=0; i<res.size(); ++i){
    if(output[i]!=0)
      *output[i] = res[i];
  }
}

void Evaluation::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  // Evaluate function
  if(!output_given){
    // Evaluate the function symbolically
    vector<MX> arg(input.size());
    for(int i=0; i<arg.size(); ++i){
      if(input[i])
        arg[i] = *input[i];
    }
    vector<MX> res = fcn_.call(arg);
    for(int i=0; i<res.size(); ++i){
      if(output[i])
        *output[i] = res[i];
    }
  }

  // Sensitivities
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  if(nfwd>0 || nadj>0){
    // Loop over outputs
    for(int oind=0; oind<output.size(); ++oind){
      // Skip of not used
      if(output[oind]==0) continue;
          
      // Output dimensions
      int od1 = output[oind]->size1();
      int od2 = output[oind]->size2();
      
      // Loop over inputs
      for(int iind=0; iind<input.size(); ++iind){
        // Skip of not used
        if(input[iind]==0) continue;

        // Input dimensions
        int id1 = input[iind]->size1();
        int id2 = input[iind]->size2();

        // Create a Jacobian node
        MX J = MX::create(new JacobianReference(*output[oind],iind));
        if(isZero(J)) continue;
        
        // Forward sensitivities
        for(int d=0; d<nfwd; ++d){
          MX fsens_d = prod(J,vec(*fwdSeed[d][iind]));
          if(!isZero(fsens_d)){
            // Reshape sensitivity contribution if necessary
            if(od2>1) fsens_d = reshape(fsens_d,od1,od2);
            
            // Save or add to vector
            if(fwdSens[d][oind]->isNull()){
              *fwdSens[d][oind] = fsens_d;
            } else {
              *fwdSens[d][oind] += fsens_d;
            }
          }
          
          // If no contribution added, set to zero
          if(fwdSens[d][oind]->isNull()){
            *fwdSens[d][oind] = MX::zeros(od1,od2);
          }
        }
        
        // Adjoint sensitivities
        for(int d=0; d<nadj; ++d){
          MX asens_d = prod(trans(J),vec(*adjSeed[d][oind]));
          if(!isZero(asens_d)){
            // Reshape sensitivity contribution if necessary
            if(id2>1) asens_d = reshape(asens_d,id1,id2);
            
            // Add to vector
            *adjSens[d][iind] += asens_d;
          }
        }
      }
    }
  }
}

void Evaluation::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  MXNode::deepCopyMembers(already_copied);
  fcn_ = deepcopy(fcn_,already_copied);
}

void Evaluation::propagateSparsity(const DMatrixPtrV& input, DMatrixPtrV& output){
  // Clear the outputs
  for(int oind=0; oind<output.size(); ++oind){
    // Skip of not used
    if(output[oind]==0) continue;
        
    // Get data array for output and clear it
    bvec_t *outputd = get_bvec_t(output[oind]->data());
    fill_n(outputd,output[oind]->size(),0);
  }
  
  // Temporary variable to hold the input and output as a dense vectors
  vector<bvec_t> input_dense, output_dense;
  
  // Loop over inputs
  for(int iind=0; iind<input.size(); ++iind){
    // Skip of not used
    if(input[iind]==0) continue;

    // Get the sparsity of the input
    const CRSSparsity& sp_in = input[iind]->sparsity();
    int id1 = sp_in.size1();
    int id2 = sp_in.size2();
    casadi_assert(id1==fcn_.input(iind).size1());
    casadi_assert(id2==fcn_.input(iind).size2());
    const vector<int>& irowind = sp_in.rowind();
    const vector<int>& icol = sp_in.col();

    // Get data array for input
    bvec_t *inputd = get_bvec_t(input[iind]->data());

    // Resize input vector
    input_dense.resize(sp_in.numel());
    
    // Copy the input vector
    fill(input_dense.begin(),input_dense.end(),0);
    
    // Copy nonzeros
    for(int i=0; i<id1; ++i){
      for(int el=irowind[i]; el<irowind[i+1]; ++el){
        // Get column
        int j=icol[el];
        
        // Copy element
        input_dense[j+i*id2] = inputd[el];
      }
    }
    
    // Loop over outputs
    for(int oind=0; oind<output.size(); ++oind){

      // Skip of not used
      if(output[oind]==0) continue;

      // Get the sparsity of the Jacobian block
      CRSSparsity& sp = fcn_.jacSparsity(iind,oind);
      if(sp.isNull() || sp.size()==0) continue; // Skip if zero
      int d1 = sp.size1();
      int d2 = sp.size2();
      const vector<int>& rowind = sp.rowind();
      const vector<int>& col = sp.col();

      // Get the sparsity of the output
      const CRSSparsity& sp_out = output[oind]->sparsity();
      int od1 = sp_out.size1();
      int od2 = sp_out.size2();
      casadi_assert(od1==fcn_.output(oind).size1());
      if(od2!=fcn_.output(oind).size2())
        cout << fcn_.get() << ": iind = " << iind << ", oind = " << oind << endl;
      casadi_assert(od2==fcn_.output(oind).size2());
      const vector<int>& orowind = sp_out.rowind();
      const vector<int>& ocol = sp_out.col();

      // Make sure that the Jacobian dimensions are consistent with the inputs and outputs
      casadi_assert(d1==sp_out.numel());
      casadi_assert(d2==sp_in.numel());

      // Get data array for output
      bvec_t *outputd = get_bvec_t(output[oind]->data());

      // Resize dense output vector
      output_dense.resize(sp_out.numel());
      
      // Clear the output vector (the parts we will use)
      for(int i=0; i<od1; ++i){
        for(int el=orowind[i]; el<orowind[i+1]; ++el){
          // Get column
          int j=ocol[el];
          
          // Clear element
          output_dense[j+i*od2] = 0;
        }
      }
      
      // Carry out the sparse matrix-vector multiplication
      for(int i=0; i<d1; ++i){
        for(int el=rowind[i]; el<rowind[i+1]; ++el){
          // Get column
          int j=col[el];
          
          // Propagate dependencies
          output_dense[i] |= input_dense[j];
        }
      }
      
      // Get the dependencies
      for(int i=0; i<od1; ++i){
        for(int el=orowind[i]; el<orowind[i+1]; ++el){
          // Get column
          int j=ocol[el];
          
          // Clear element
          outputd[el] |= output_dense[j+i*od2];
        }
      }
    }
  }
}



} // namespace CasADi
