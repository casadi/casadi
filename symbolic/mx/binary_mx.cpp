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

#include "binary_mx.hpp"
#include "binary_mx_impl.hpp"
#include "mx_tools.hpp"
#include <vector>
#include <sstream>
#include <cassert>
#include "../matrix/sparsity_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

  MatrixScalarOp::MatrixScalarOp(Operation op, const MX& x, const MX& y) : BinaryMX<false,true>(op,x,y){
    setSparsity(x.sparsity());
  }

  template<typename T, typename MatV, typename MatVV> 
  void MatrixScalarOp::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Get data
    vector<T>& output0 = output[0]->data();
    const vector<T> &input0 = input[0]->data();
    const vector<T> &input1 = input[1]->data();
  
    if(nfwd==0 && nadj==0){
      casadi_math<T>::fun(op_,getPtr(input0),input1[0],getPtr(output0),input0.size());
    } else {

      // Function and partial derivatives
      T f,pd[2];

      for(int el=0; el<input0.size(); ++el){
	casadi_math<T>::fun(op_,input0[el],input1[0],f);
	casadi_math<T>::der(op_,input0[el],input1[0],f,pd);
	output0[el] = f;

	// Propagate forward seeds
	for(int d=0; d<nfwd; ++d){
	  fwdSens[d][0]->data()[el] = pd[0]*fwdSeed[d][0]->data()[el] + pd[1]*fwdSeed[d][1]->data()[0];
	}
    
	// Propagate adjoint seeds
	for(int d=0; d<nadj; ++d){
	  T s = adjSeed[d][0]->data()[el];
	  adjSeed[d][0]->data()[el] = 0;
	  adjSens[d][0]->data()[el] += s*pd[0];
	  adjSens[d][1]->data()[0]  += s*pd[1];
	}
      }
    }
  }

  void MatrixScalarOp::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void MatrixScalarOp::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  ScalarMatrixOp::ScalarMatrixOp(Operation op, const MX& x, const MX& y) : BinaryMX<true,false>(op,x,y){
    setSparsity(y.sparsity());
  }

  template<typename T, typename MatV, typename MatVV> 
  void ScalarMatrixOp::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Get data
    vector<T>& output0 = output[0]->data();
    const vector<T> &input0 = input[0]->data();
    const vector<T> &input1 = input[1]->data();
  
    if(nfwd==0 && nadj==0){
      casadi_math<T>::fun(op_,input0[0],getPtr(input1),getPtr(output0),input1.size());
    } else {

      // Function and partial derivatives
      T f,pd[2];

      for(int el=0; el<input1.size(); ++el){      
	casadi_math<T>::fun(op_,input0[0],input1[el],f);
	casadi_math<T>::der(op_,input0[0],input1[el],f,pd);
	output0[el] = f;

	// Propagate forward seeds
	for(int d=0; d<nfwd; ++d){
	  fwdSens[d][0]->data()[el] = pd[0]*fwdSeed[d][0]->data()[0] + pd[1]*fwdSeed[d][1]->data()[el];
	}
    
	// Propagate adjoint seeds
	for(int d=0; d<nadj; ++d){
	  T s = adjSeed[d][0]->data()[el];
	  adjSeed[d][0]->data()[el] = 0;
	  adjSens[d][0]->data()[0]  += s*pd[0];
	  adjSens[d][1]->data()[el] += s*pd[1];
	}
      }
    }
  }

  void ScalarMatrixOp::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void ScalarMatrixOp::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }


  MatrixMatrixOp::MatrixMatrixOp(Operation op, const MX& x, const MX& y) : BinaryMX<false,false>(op,x,y){
    setSparsity(x.sparsity());
  }

  template<typename T, typename MatV, typename MatVV> 
  void MatrixMatrixOp::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Get data
    vector<T>& output0 = output[0]->data();
    const vector<T> &input0 = input[0]->data();
    const vector<T> &input1 = input[1]->data();
  
    if(nfwd==0 && nadj==0){
      casadi_math<T>::fun(op_,getPtr(input0),getPtr(input1),getPtr(output0),input0.size());
    } else {

      // Function and partial derivatives
      T f,pd[2];

      for(int el=0; el<input0.size(); ++el){
	casadi_math<T>::fun(op_,input0[el],input1[el],f);
	casadi_math<T>::der(op_,input0[el],input1[el],f,pd);
	output0[el] = f;

	// Propagate forward seeds
	for(int d=0; d<nfwd; ++d){
	  fwdSens[d][0]->data()[el] = pd[0]*fwdSeed[d][0]->data()[el] + pd[1]*fwdSeed[d][1]->data()[el];
	}
    
	// Propagate adjoint seeds
	for(int d=0; d<nadj; ++d){
	  T s = adjSeed[d][0]->data()[el];
	  adjSeed[d][0]->data()[el] = 0;
	  adjSens[d][0]->data()[el] += s*pd[0];
	  adjSens[d][1]->data()[el] += s*pd[1];
	}
      }
    }
  }

  void MatrixMatrixOp::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void MatrixMatrixOp::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void MatrixMatrixOp::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    bvec_t *input0 = get_bvec_t(input[0]->data());
    bvec_t *input1 = get_bvec_t(input[1]->data());
    bvec_t *outputd = get_bvec_t(output[0]->data());
    for(int el=0; el<output[0]->size(); ++el){
      if(fwd){
	outputd[el] = input0[el] | input1[el];
      } else {
	bvec_t s = outputd[el];
	outputd[el] = bvec_t(0);
	input0[el] |= s;
	input1[el] |= s;
      }
    }
  }

  void MatrixScalarOp::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    bvec_t *input0 = get_bvec_t(input[0]->data());
    bvec_t *input1 = get_bvec_t(input[1]->data());
    bvec_t *outputd = get_bvec_t(output[0]->data());
    for(int el=0; el<output[0]->size(); ++el){
      if(fwd){
	outputd[el] = input0[el] | input1[0];
      } else {
	bvec_t s = outputd[el];
	outputd[el] = bvec_t(0);
	input0[el] |= s;
	input1[0]  |= s;
      }
    }
  }

  void ScalarMatrixOp::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    bvec_t *input0 = get_bvec_t(input[0]->data());
    bvec_t *input1 = get_bvec_t(input[1]->data());
    bvec_t *outputd = get_bvec_t(output[0]->data());
    for(int el=0; el<output[0]->size(); ++el){
      if(fwd){
	outputd[el] = input0[0] | input1[el];
      } else {
	bvec_t s = outputd[el];
	outputd[el] = bvec_t(0);
	input0[0]  |= s;
	input1[el] |= s;
      }
    }
  }

} // namespace CasADi

