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

#include "jacobian_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "../pre_c99_support.hpp"
#include "../matrix/sparsity_tools.hpp"

#include <stack>
#include <typeinfo>
#include <cassert>

using namespace std;

namespace CasADi{
  
JacobianInternal::JacobianInternal(const FX& fcn, const std::vector<std::pair<int,int> >& jblocks) : fcn_(fcn), jblocks_(jblocks){
  addOption("ad_mode",            OT_STRING,  "default");  // "forward", "adjoint" or "default" (default meaning checking both)
}
  
JacobianInternal::JacobianInternal(const FX& fcn, int iind, int oind) : fcn_(fcn), iind_(iind), oind_(oind){
  addOption("finite_differences", OT_BOOLEAN,   false);    // Using finite differences instead of automatic differentiation
  addOption("ad_mode",            OT_STRING,  "default");  // "forward", "adjoint" or "default", i.e. forward if n_<=m_, otherwise adjoint
  setOption("sparse", false);
  
  //casadi_assert_message(fcn_.input(iind_).size() == fcn_.input(iind_).numel(),"Jacobian not yet implemented for sparse input");
  //casadi_assert_message(fcn_.output(oind_).size()== fcn_.output(oind_).numel(),"Jacobian not yet implemented for sparse output");
  
  // make sure that input and output are vectors (not matrices)
  //casadi_assert(fcn_.input(iind_).size2() == 1 or fcn_.input(iind_).size1() == 1);
  //casadi_assert(fcn_.output(oind_).size2() == 1 or fcn_.output(oind_).size1() == 1);

  // get the dimensions
  n_ = fcn_.input(iind_).size();
  m_ = fcn_.output(oind_).size();

  input_ = fcn_->input_;
  
  output_.resize(1);
  
  // Input/Output Sparsity 
  if (fcn_.input(iind_).size() == fcn_.input(iind_).numel() && fcn_.output(oind_).size()== fcn_.output(oind_).numel()) {
    output(0) = DMatrix(m_,n_,0); // dense
  } else {
    // sparse
    output(0) = DMatrix(
                  sp_rowcol(
                    getNZDense(fcn_.output(oind_).sparsity()),
                    getNZDense(fcn_.input(iind_).sparsity()),
                    fcn_.output(oind_).numel(), fcn_.input(iind_).numel()
                  ),
                  0
                );
  }
  
}



JacobianInternal::~JacobianInternal(){
}

JacobianInternal* JacobianInternal::clone() const{
  JacobianInternal* node = new JacobianInternal(*this);
  node->fcn_ = shared_cast<FX>(fcn_.clone());
  return node;
}

void JacobianInternal::init(){
  // Use old or new (sparse) Jacobian calculation algorithm
  if(jblocks_.empty()){     // Old (to be depreciated) algorithm
    // Call the init function of the base class
    FXInternal::init();

    if(!fcn_.isInit()) fcn_.init();
    
    // Number of directions that we can calculate at a time
    nadir_fcn_ = fcn_.getOption("number_of_adj_dir");
    nfdir_fcn_ = fcn_.getOption("number_of_fwd_dir");
    casadi_assert(nadir_fcn_>0 || nfdir_fcn_>0);
  
    // Use finite differences?
    use_fd_ = getOption("finite_differences").toBool();
    
    // Use forward or adjoint ad?
    if(use_fd_){
      use_ad_fwd_ = true;
    } else {
      if(getOption("ad_mode") == "forward")
        use_ad_fwd_ = true;
      else if(getOption("ad_mode") == "adjoint")
        use_ad_fwd_ = false;
      else if(getOption("ad_mode") == "default")
        use_ad_fwd_ = nadir_fcn_==0 || n_<=m_;
      else
        throw CasadiException("unknown ad mode: " +  getOption("ad_mode").toString());
    }
  } else { // New (sparse) algorithm
    
    // Initialize the function if not already initialized
    if(!fcn_.isInit()) fcn_.init();

    // Same input dimensions as function
    setNumInputs(fcn_.getNumInputs());
    for(int i=0; i<getNumInputs(); ++i){
      input(i) = fcn_.input(i);
    }

    // Get the output dimensions (specified by jblocks)
    setNumOutputs(jblocks_.size());
    for(int i=0; i<jblocks_.size(); ++i){
      
      // Get input and output indices
      int oind = jblocks_[i].first;
      int iind = jblocks_[i].second;
      
      // Jacobian block or nondifferentiated function output
      if(iind<0){
        // Nondifferentiated function
        output(i) = fcn_.output(oind);
      } else {
        // Jacobian of input iind with respect to output oind
        output(i) = DMatrix(fcn_.jacSparsity(iind,oind),0);
      }
    }
    
    // Only now, when we know the sparsity of the inputs and outputs, we are anle to call the init function of the base class
    FXInternal::init();

    // Number of directions that we can calculate at a time
    nadir_fcn_ = fcn_.getOption("number_of_adj_dir");
    nfdir_fcn_ = fcn_.getOption("number_of_fwd_dir");
    casadi_assert_message(nadir_fcn_>0 || nfdir_fcn_>0, "Function does not support directional derivatives, neither \"number_of_fwd_dir\" nor \"number_of_adj_dir\" is positive");
    
    
    
    // Quick hack!
    iind_ = -1;
    oind_ = -1;
    js_trans_ = CRSSparsity();
    for(int i=0; i<jblocks_.size(); ++i){
      if(jblocks_[i].second>=0){
        casadi_assert_message(js_trans_.isNull(), "Maximum one Jacobian block currently supported");
        js_ = output(i).sparsity();
        js_trans_ = js_.transpose(js_trans_mapping_);
        
        oind_ = jblocks_[i].first;
        iind_ = jblocks_[i].second;
        vector<pair<int,int> > jblocks_no_f(1,jblocks_[i]);
        D1_.resize(1);
        D2_.resize(1);
        fcn_->getPartition(jblocks_no_f,D1_,D2_);

      
        if(verbose()){
          int nfwd = D1_.front().isNull() ? 0 : D1_.front().size1();
          int nadj = D2_.front().isNull() ? 0 : D2_.front().size1();
          cout << "JacobianInternal::init: " << nfwd << " forward directions and " << nadj << " adjoint directions needed for the jacobian" << endl;
        }
      }
    }
  }
}

void JacobianInternal::evaluate(int nfdir, int nadir){
  
  // Pass the argument to the function
  for(int i=0; i<input_.size(); ++i)
    fcn_.setInput(input(i),i);

    // Use old or new (sparse) Jacobian calculation algorithm
  if(jblocks_.empty()){     // Old (to be depreciated) algorithm

    vector<double>& res2 = output().data();
    
    if(use_ad_fwd_){ // forward AD if less inputs than outputs
      // Clear the forward seeds
      for(int i=0; i<fcn_.getNumInputs(); ++i)
        for(int dir=0; dir<nfdir_fcn_; ++dir)
          fcn_.fwdSeed(i,dir).setZero();

      
      // Calculate the forward sensitivities, nfdir_fcn_ directions at a time
      for(int ofs=0; ofs<n_; ofs += nfdir_fcn_){
          for(int dir=0; dir<nfdir_fcn_ && ofs+dir<n_; ++dir){
            // Pass forward seeds
            int i=ofs+dir;
            vector<double>& fseed = fcn_.fwdSeed(iind_,dir).data();
            fill(fseed.begin(),fseed.end(),0);
            fseed[i] = 1;
          }
        
        // Evaluate the AD forward algorithm
        fcn_.evaluate(nfdir_fcn_,0);
            
        // Get the output seeds
        for(int dir=0; dir<nfdir_fcn_ && ofs+dir<n_; ++dir){
          // Save to the result
          int i=ofs+dir;
          const vector<double>& fsens = fcn_.fwdSens(oind_,dir).data();
          for(int j=0; j<m_; ++j){
            res2[i+j*n_] = fsens[j];
          }
        }
      }
    } else { // adjoint AD
      for(int ofs=0; ofs<m_; ofs += nadir_fcn_){
          for(int dir=0; dir<nadir_fcn_ && ofs+dir<m_; ++dir){
            // Pass forward seeds
            int j=ofs+dir;
            vector<double>& aseed = fcn_.adjSeed(oind_,dir).data();
            fill(aseed.begin(),aseed.end(),0);
            aseed[j] = 1;
          }
        
        // Evaluate the AD forward algorithm
        fcn_.evaluate(0,nadir_fcn_);
            
        // Get the output seeds
        for(int dir=0; dir<nadir_fcn_ && ofs+dir<m_; ++dir){
          // Save to the result
          int j=ofs+dir;
          const vector<double>& asens = fcn_.adjSens(iind_,dir).data();
          for(int i=0; i<n_; ++i){
            res2[i+j*n_] = asens[i];
          }
        }
      }
    }
  } else { // New (sparse) algorithm
    
    // Has the function been called at least once?
    bool called_once = false;
    
    // Loop over outputs
    for(int ind=0; ind<jblocks_.size(); ++ind){
      
      // Get input and output indices
      int oind = jblocks_[ind].first;
      int iind = jblocks_[ind].second;
      
      // Skip if nondifferentiated function output
      if(iind<0) continue;
      
      // Get the number of forward and adjoint sweeps
      int nfwd = D1_.front().isNull() ? 0 : D1_.front().size1();
      int nadj = D2_.front().isNull() ? 0 : D2_.front().size1();
      casadi_assert(nfwd || nadj);
      casadi_assert(!(nfwd && nadj));
      
      // Calculate the forward sensitivities, nfdir_fcn_ directions at a time
      for(int ofs=0; ofs<nfwd; ofs += nfdir_fcn_){
        
        // Number of directions
        int ndir = std::min(nfdir_fcn_,nfwd-ofs);
        
        // Clear the forward seeds
        for(int i=0; i<fcn_.getNumInputs(); ++i)
          for(int dir=0; dir<ndir; ++dir)
            fcn_.fwdSeed(i,dir).setZero();
        
        // Pass forward seeds
        for(int dir=0; dir<ndir; ++dir){
          int d = ofs+dir;
          for(int el=D1_[0].rowind(d); el<D1_[0].rowind(d+1); ++el){
            // Get the direction
            int j=D1_[0].col(el);
            
            // Set the seed
            Matrix<double>& fseed = fcn_.fwdSeed(iind,dir);
            casadi_assert(j<fseed.numel());
            int i_fseed = j/fseed.size2();
            int j_fseed = j - i_fseed*fseed.size2();
            int nzind = fcn_.fwdSeed(iind,dir).sparsity().getNZ(i_fseed,j_fseed);
            if(nzind>=0){
              fcn_.fwdSeed(iind,dir).elem(i_fseed,j_fseed) = 1;
            }
          }
        }
      
        // Evaluate the AD forward algorithm
        fcn_.evaluate(ndir,0);
        called_once = true;

        // Get the output seeds
        for(int dir=0; dir<ndir; ++dir){
          int d = ofs+dir;
          for(int el=D1_[0].rowind(d); el<D1_[0].rowind(d+1); ++el){
            // Get the direction
            int j=D1_[0].col(el);

            // Save to the result
            const Matrix<double>& fsens = fcn_.fwdSens(oind,dir);
            
            // Loop over the rows using this variable
            for(int el2=js_trans_.rowind(j); el2<js_trans_.rowind(j+1); ++el2){
              
              // Get the row
              int i = js_trans_.col(el2);
              
              // Store the value
              casadi_assert(i<fsens.numel());
              int i_fsens = i/fsens.size2();
              int j_fsens = i - i_fsens*fsens.size2();
              int nzind = fsens.sparsity().getNZ(i_fsens,j_fsens);
              //casadi_assert(nzind>=0);
              if(nzind>=0){
                output(ind).data()[js_trans_mapping_[el2]] = fsens.elem(i_fsens,j_fsens); // quick hack!
              }
            }
          }
        }
      }

      // Calculate the forward sensitivities, nfdir_fcn_ directions at a time
      for(int ofs=0; ofs<nadj; ofs += nadir_fcn_){
        // Number of directions
        int ndir = std::min(nadir_fcn_,nadj-ofs);
        
        // Clear the adjoint seeds
        for(int i=0; i<fcn_.getNumOutputs(); ++i)
          for(int dir=0; dir<ndir; ++dir)
            fcn_.adjSeed(i,dir).setZero();
          
        // Set the adjoint seeds
        for(int dir=0; dir<ndir; ++dir){
          int d = ofs+dir;
          for(int el=D2_[0].rowind(d); el<D2_[0].rowind(d+1); ++el){
            // Get the output
            int i=D2_[0].col(el);

            // Save to the result
            Matrix<double>& aseed = fcn_.adjSeed(oind,dir);
            
            // Store the value
            int i_aseed = i/aseed.size2();
            int j_aseed = i - i_aseed*aseed.size2();
            int nzind = aseed.sparsity().getNZ(i_aseed,j_aseed);
            casadi_assert(nzind>=0);
            if(nzind>=0){
              aseed.elem(i_aseed,j_aseed) += 1;
            }
          }
        }
        
        // Evaluate the AD adjoint algorithm
        fcn_.evaluate(0,ndir);
        called_once = true;

        // Set the adjoint sensitivities
        for(int dir=0; dir<ndir; ++dir){
          int d = ofs+dir;
          for(int el=D2_[0].rowind(d); el<D2_[0].rowind(d+1); ++el){
            // Get the direction
            int i=D2_[0].col(el);

            // Save to the result
            const Matrix<double>& asens = fcn_.adjSens(iind,dir);
            
            // Loop over the rows using this variable
            for(int el2=js_.rowind(i); el2<js_.rowind(i+1); ++el2){
              
              // Get the column
              int j = js_.col(el2);
              
              // Store the value
              int i_asens = j/asens.size2();
              int j_asens = j - i_asens*asens.size2();
              int nzind = asens.sparsity().getNZ(i_asens,j_asens);
              casadi_assert(nzind>=0);
              output(ind).data()[el2] = asens.elem(i_asens,j_asens); // quick hack!
            }
          }
        }
      }
    }
    
    // Evaluate if necessary
    if(!called_once) fcn_.evaluate();
    
    // Loop over outputs
    for(int ind=0; ind<jblocks_.size(); ++ind){
      // Get input and output indices
      int oind = jblocks_[ind].first;
      int iind = jblocks_[ind].second;
      if(iind<0){ 
        // Nondifferentiated function
        fcn_.output(oind).get(output(ind));
      }
    }
  }
}

// void JacobianInternal::evaluateNew(int fsens_order, int asens_order){
//   casadi_assert(fsens_order==0);
//   casadi_assert(asens_order==0);
//   
//   // Loop over the seeding directions
// /*  for(int i=0; i<*/
//   
//   
//   
// }



// void JacobianInternal::evaluate(int fsens_order, int asens_order){
//   // Pass the argument to the function
//   for(int i=0; i<input_.size(); ++i)
//     fcn_.input(i).set(input(i).data());
// 
//   vector<double>& res2 = output().data();
//   
//   int el = 0; // running index
//     
//   if(n_<=m_){ // forward AD if less inputs than outputs
//     // Input and output seed vectors
//     std::vector<double>& fseed = fcn_.input(iind_).dataF();
//     std::vector<double>& fsens = fcn_.output(oind_).dataF();
//     
//     // forward AD in all directions
//     for(int i=0; i<n_; ++i){
//       // Give a seed in the k-th direction
//       fill(fseed.begin(),fseed.end(),0);
//       fseed[i] = 1;
// 
//       // Execute forward AD algorithm
//       fcn_.evaluate(1,0);
// 
//       // Save to the result
//       if(sparse_jac_){
//         while(el<cc_.size() && cc_[el]==i){
//           res2[elind_[el]] = fsens[rr_[el]];
//           el++;
//         }
//       } else {
//         for(int j=0; j<m_; ++j){
//           res2[i+j*n_] = fsens[j];
//         }
//       }
//     }
//   } else { // adjoint AD
//     // Input and output seed vectors
//     std::vector<double>& aseed = fcn_.output(oind_).dataA();
//     std::vector<double>& asens = fcn_.input(iind_).dataA();
//     
//     // Adjoint AD in all directions
//     for(int j=0; j<m_; ++j){
//       // Give a backward seed in the j-th direction
//       fill(aseed.begin(),aseed.end(),0);
//       aseed[j] = 1;
// 
//       // Execute the adjoint AD algorithm
//       fcn_.evaluate(0,1);
// 
//       // Save to the result
//       if(sparse_jac_){
//         while(el<rr_.size() && rr_[el]==j){
//           res2[elind_[el]] = asens[cc_[el]];
//           el++;
//         }
//       } else {
//         for(int i=0; i<n_; ++i){
//           res2[i+j*n_] = asens[i];
//         }
//       }
//     }
//   }
// }



// void JacobianInternal::compress(){
//   casadi_assert(isInit());
//   
// }

} // namespace CasADi

