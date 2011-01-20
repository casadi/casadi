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

#include <stack>
#include <typeinfo>
#include <cassert>

using namespace std;

namespace CasADi{
  
JacobianInternal::JacobianInternal(const FX& fcn, int iind, int oind) : fcn_(fcn), iind_(iind), oind_(oind){
  addOption("finite_differences", OT_BOOLEAN,   false);    // Using finite differences instead of automatic differentiation
  addOption("ad_mode",            OT_STRING,  "default");  // "forward", "adjoint" or "default", i.e. forward if n_<=m_, otherwise adjoint
  setOption("sparse", false);

  // make sure that input and output are vectors (not matrices)
  if(fcn_->input(iind_).size2() != 1) throw CasadiException("JacobianInternal::JacobianInternal: input not a vector");
  if(fcn_->result(oind_).size2() != 1) throw CasadiException("JacobianInternal::JacobianInternal: output not a vector");

  // get the dimensions
  
  n_ = fcn_.input(iind_).size1();
  m_ = fcn_.result(oind_).size1();

  input_ = fcn_->input_;
  
  output_.resize(1);
  result(0).resize(m_,n_);
//  output[1] = fcn_->output[oind_];
}



JacobianInternal::~JacobianInternal(){
}


void JacobianInternal::init(){
  // Call the init function of the base class
  FXInternal::init();

  // Number of directions that we can calculate at a time
  nadir_fcn_ = getOption("number_of_adj_dir").toInt();
  nfdir_fcn_ = getOption("number_of_fwd_dir").toInt();
  
  // Use finite differences?
  use_fd_ = getOption("finite_differences").toBool();
  
  // Use forward or adjoint ad?
  if(use_fd_){
    if(getOption("ad_mode") == "forward")
      use_ad_fwd_ = true;
    else if(getOption("ad_mode") == "adjoint")
      use_ad_fwd_ = false;
    else if(getOption("ad_mode") == "default")
      use_ad_fwd_ = n_<=m_;
    else
      throw CasadiException("unknown ad mode: " +  getOption("ad_mode").toString());
  }

  sparse_jac_ = getOption("sparse").toBool();
  
  // Quick return if dense
  if(!sparse_jac_){
    return;
  }
  
  assert(0);
  throw CasadiException("JacobianInternal::init: sparse deactivated");
  
  // Input and output seed vectors
  std::vector<double>& dir = fcn_.fwdSeed(iind_);
  std::vector<double>& bdir = fcn_.fwdSens(oind_); // TODO: fix this
  
  if(use_ad_fwd_){ // forward AD if less inputs than outputs 
      assert(0);

      // forward AD in all directions
      for(int i=0; i<n_; ++i){ // (COLUMN-BY-COLUMN)
	
	// Give a seed in the k-th direction
        fill(dir.begin(),dir.end(),0);
	dir[i] = numeric_limits<double>::quiet_NaN();
      
	// Execute forward AD algorithm
	fcn_.evaluate(1,0);

	// Save to the result
	for(int j=0; j<m_; ++j)
	  if(ISNAN(bdir[j])){
	    rr_.push_back(j);
	    cc_.push_back(i);
	  }
      }
      
    } else { // adjoint AD 
      // Adjoint AD in all directions
      for(int j=0; j<m_; ++j){ // (ROW-BY-ROW)
	// Give a backward seed in the j-th direction
        fill(bdir.begin(),bdir.end(),0);
	bdir[j] = numeric_limits<double>::quiet_NaN();
            
	// Execute the adjoint AD algorithm
	fcn_.evaluate(0,1);

	// Save to the result
	for(int i=0; i<n_; ++i)
	  if(ISNAN(dir[i])){
	    rr_.push_back(j);
	    cc_.push_back(i);
	  }
      }
    }

  // Set the sparsity of the output
//  output_.at(0).setSparsity(rr_,cc_,elind_);
    
}

void JacobianInternal::evaluate(int fsens_order, int asens_order){
  // Pass the argument to the function
  for(int i=0; i<input_.size(); ++i)
    fcn_.setInput(input(i),i);
  vector<double>& res2 = result();
  
  int el = 0; // running index

  if(use_ad_fwd_){ // forward AD if less inputs than outputs
    // Calculate the forward sensitivities, nfdir_fcn_ directions at a time
    for(int ofs=0; ofs<n_; ofs += nfdir_fcn_){
        for(int dir=0; dir<nfdir_fcn_ && ofs+dir<n_; ++dir){
          // Pass forward seeds
          int i=ofs+dir;
          vector<double>& fseed = fcn_.fwdSeed(iind_,dir);
          fill(fseed.begin(),fseed.end(),0);
          fseed[i] = 1;
        }
      
      // Evaluate the AD forward algorithm
      fcn_.evaluate(1,0);
          
      // Get the output seeds
      for(int dir=0; dir<nfdir_fcn_ && ofs+dir<n_; ++dir){
        // Save to the result
        int i=ofs+dir;
        const vector<double>& fsens = fcn_.fwdSens(oind_,dir);
        if(sparse_jac_){
          while(el<cc_.size() && cc_[el]==i){
            res2[elind_[el]] = fsens[rr_[el]];
            el++;
          }
          } else {
            for(int j=0; j<m_; ++j){
              res2[i+j*n_] = fsens[j];
            }
          }
      }
    }
  } else { // adjoint AD
    for(int ofs=0; ofs<m_; ofs += nadir_fcn_){
        for(int dir=0; dir<nadir_fcn_ && ofs+dir<m_; ++dir){
          // Pass forward seeds
          int j=ofs+dir;
          vector<double>& aseed = fcn_.adjSeed(oind_,dir);
          fill(aseed.begin(),aseed.end(),0);
          aseed[j] = 1;
        }
      
      // Evaluate the AD forward algorithm
      fcn_.evaluate(0,1);
          
      // Get the output seeds
      for(int dir=0; dir<nadir_fcn_ && ofs+dir<m_; ++dir){
        // Save to the result
        int j=ofs+dir;
        const vector<double>& asens = fcn_.adjSens(iind_,dir);
        if(sparse_jac_){
          while(el<rr_.size() && rr_[el]==j){
            res2[elind_[el]] = asens[cc_[el]];
            el++;
          }
        } else {
          for(int i=0; i<n_; ++i){
            res2[i+j*n_] = asens[i];
          }
        }
      }
    }
  }
}



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





















} // namespace CasADi

