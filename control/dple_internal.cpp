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

#include "dple_internal.hpp"
#include <cassert>
#include "../symbolic/std_vector_tools.hpp"
#include "../symbolic/matrix/matrix_tools.hpp"
#include "../symbolic/mx/mx_tools.hpp"
#include "../symbolic/sx/sx_tools.hpp"
#include "../symbolic/fx/mx_function.hpp"
#include "../symbolic/fx/sx_function.hpp"

INPUTSCHEME(DPLEInput)
OUTPUTSCHEME(DPLEOutput)

using namespace std;
namespace CasADi{

  DpleInternal::DpleInternal(const std::vector< Sparsity > & A, const std::vector< Sparsity > &V,int nfwd, int nadj) : A_(A), V_(V), nfwd_(nfwd), nadj_(nadj) {

    // set default options
    setOption("name","unnamed_dple_solver"); // name of the function 
    
    addOption("const_dim",OT_BOOLEAN,true,"Assume constant dimension of P");
    addOption("pos_def",OT_BOOLEAN,false,"Assume P positive definite");
    
    addOption("error_unstable",OT_BOOLEAN,false,"Throw an exception when it is detected that Product(A_i,i=N..1) has eigenvalues greater than 1-eps_unstable");
    addOption("eps_unstable",OT_REAL,1e-4,"A margin for unstability detection");
    
    if (nfwd_==0 && nadj_==0) {
      input_.scheme = SCHEME_DPLEInput;
      output_.scheme = SCHEME_DPLEOutput;
    }

  }

  DpleInternal::~DpleInternal(){ 

  }

  void DpleInternal::init(){
    
    const_dim_ = getOption("const_dim");
    pos_def_ = getOption("pos_def");
    error_unstable_ = getOption("error_unstable");
    eps_unstable_ = getOption("eps_unstable");
    
    // Dimension sanity checks
    casadi_assert_message(A_.size()==V_.size(),"A and V arguments must be of same length, but got " << A_.size() << " and " << V_.size() << ".");
    K_ = A_.size();
    for (int k=0;k<K_;++k) {
      casadi_assert_message(V_[k].isSymmetric(),"V_i must be symmetric but got " << V_[k].dimString() << " for i = " << k << ".");

      casadi_assert_message(A_[k].size1()==V_[k].size1(),"First dimension of A (" << A_[k].size1() << ") must match dimension of symmetric V_i (" << V_[k].size1() << ")" << " for i = " << k << ".");
    }
    
    if (const_dim_) {
      int n = A_[0].size1();
       for (int k=1;k<K_;++k) {
         casadi_assert_message(A_[k].size1()==n,"You have set const_dim option, but found an A_i with dimension ( " << A_[k].dimString() << " ) deviating from n = " << n << " at i = " << k << ".");
      }
    }

    // Allocate inputs
    setNumInputs(DPLE_NUM_IN*(1+nfwd_) + DPLE_NUM_OUT*nadj_);
    
    for (int i=0;i<nfwd_+1;++i) {
      if (const_dim_) {
        input(DPLE_NUM_IN*i+DPLE_A)  = DMatrix::zeros(horzcat(A_));
        input(DPLE_NUM_IN*i+DPLE_V)  = DMatrix::zeros(horzcat(V_));
      } else {
        input(DPLE_NUM_IN*i+DPLE_A)  = DMatrix::zeros(blkdiag(A_));
        input(DPLE_NUM_IN*i+DPLE_V)  = DMatrix::zeros(blkdiag(V_));
      }
    }
    for (int i=0;i<nadj_;++i) {
      if (const_dim_) {
        input(DPLE_NUM_IN*(1+nfwd_)+DPLE_NUM_OUT*i+DPLE_P)  = DMatrix::zeros(horzcat(A_));
      } else {
        input(DPLE_NUM_IN*(1+nfwd_)+DPLE_NUM_OUT*i+DPLE_P)  = DMatrix::zeros(blkdiag(A_));
      }
    }
    
    // Allocate outputs
    std::vector<Sparsity> P; 
    for (int k=0;k<K_;++k) {
      P.push_back(Sparsity::dense(V_[k].size1(),V_[k].size1()));
    }
    setNumOutputs(DPLE_NUM_OUT*(1+nfwd_) + DPLE_NUM_IN*nadj_);
    for (int i=0;i<nfwd_+1;++i) {
      if (const_dim_) {
        output(DPLE_NUM_OUT*i+DPLE_P) = DMatrix::zeros(horzcat(P));
      } else {
        output(DPLE_NUM_OUT*i+DPLE_P) = DMatrix::zeros(blkdiag(P));
      }
    }
    for (int i=0;i<nadj_;++i) {
      if (const_dim_) {
        output(DPLE_NUM_OUT*(nfwd_+1)+DPLE_NUM_IN*i+DPLE_A)  = DMatrix::zeros(horzcat(A_));
        output(DPLE_NUM_OUT*(nfwd_+1)+DPLE_NUM_IN*i+DPLE_V)  = DMatrix::zeros(horzcat(V_));
      } else {
        output(DPLE_NUM_OUT*(nfwd_+1)+DPLE_NUM_IN*i+DPLE_A)  = DMatrix::zeros(blkdiag(A_));
        output(DPLE_NUM_OUT*(nfwd_+1)+DPLE_NUM_IN*i+DPLE_V)  = DMatrix::zeros(blkdiag(V_));
      }
    }
  
    FXInternal::init();
      
  }

  void DpleInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    FXInternal::deepCopyMembers(already_copied);
  }


} // namespace CasADi


