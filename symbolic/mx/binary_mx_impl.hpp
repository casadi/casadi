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

#ifndef BINARY_MX_IMPL_HPP
#define BINARY_MX_IMPL_HPP

#include "binary_mx.hpp"
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
  
  template<bool ScX, bool ScY>
  BinaryMX<ScX,ScY>::BinaryMX(Operation op, const MX& x, const MX& y) : op_(op){
    setDependencies(x,y);
    if(ScX){
      setSparsity(y.sparsity());
    } else {
      setSparsity(x.sparsity());
    }
  }

  template<bool ScX, bool ScY>
  BinaryMX<ScX,ScY>::~BinaryMX(){
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX,ScY>::printPart(std::ostream &stream, int part) const{
    if(part==0){
      casadi_math<double>::printPre(op_,stream);
    } else if(part==1){
      casadi_math<double>::printSep(op_,stream);
    } else {
      casadi_math<double>::printPost(op_,stream);
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX,ScY>::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    // Evaluate function
    MX f; // Function value
    if(output_given){
      f = *output[0];
    } else {
      casadi_math<MX>::fun(op_,*input[0],*input[1],f);
    }

    // Number of forward directions
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();
    if(nfwd>0 || nadj>0){
      // Get partial derivatives
      MX pd[2];
      casadi_math<MX>::der(op_,*input[0],*input[1],f,pd);
    
      // Propagate forward seeds
      for(int d=0; d<nfwd; ++d){
        *fwdSens[d][0] = pd[0]*(*fwdSeed[d][0]) + pd[1]*(*fwdSeed[d][1]);
      }
    
      // Propagate adjoint seeds
      for(int d=0; d<nadj; ++d){
        MX s = *adjSeed[d][0];
        *adjSeed[d][0] = MX();
        for(int c=0; c<2; ++c){
          // Get increment of sensitivity c
          MX t = pd[c]*s;
          
          // If dimension mismatch (i.e. one argument is scalar), then sum all the entries
          if(!t.scalar() && t.shape() != dep(c).shape()){
            if(pd[c].shape()!=s.shape()) pd[c] = MX(s.sparsity(),pd[c]);
            t = inner_prod(pd[c],s);
          }
          
          // Propagate the seeds
          *adjSens[d][c] += t;
        }
      }
    }

    if(!output_given){
      *output[0] = f;
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX,ScY>::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{

    // Print loop and right hand side
    stream << "  for(i=0; i<" << sparsity().size() << "; ++i) ";
    stream << res.at(0) << "[i]";

    // Check if inplace
    bool inplace;
    switch(op_){
    case OP_ADD:
    case OP_SUB:
    case OP_MUL:
    case OP_DIV:
      inplace = res[0].compare(arg[0]) == 0;
    default:
      inplace = false;
    }
    
    if(inplace){
      casadi_math<double>::printSep(op_,stream);
      stream << "=";
      stream << arg.at(1) << (ScY ? "[0]" : "[i]");
    } else {
      stream << "=";
      casadi_math<double>::printPre(op_,stream);
      stream << arg.at(0) << (ScX ? "[0]" : "[i]");
      casadi_math<double>::printSep(op_,stream);
      stream << arg.at(1) << (ScY ? "[0]" : "[i]");
      casadi_math<double>::printPost(op_,stream);
    }

    stream << ";" << endl;
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX,ScY>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX,ScY>::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<bool ScX, bool ScY>
  template<typename T, typename MatV, typename MatVV>
  void BinaryMX<ScX,ScY>::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Get data
    vector<T>& output0 = output[0]->data();
    const vector<T> &input0 = input[0]->data();
    const vector<T> &input1 = input[1]->data();
  
    if(nfwd==0 && nadj==0){
      if(!ScX && !ScY){
        casadi_math<T>::fun(op_, getPtr(input0), getPtr(input1), getPtr(output0), output0.size());
      } else if(ScX){
        casadi_math<T>::fun(op_, input0[0],      getPtr(input1), getPtr(output0), output0.size());
      } else {
        casadi_math<T>::fun(op_, getPtr(input0), input1[0],      getPtr(output0), output0.size());
      }

    } else {

      // Function and partial derivatives
      T f,pd[2];

      for(int el=0; el<output0.size(); ++el){
        casadi_math<T>::fun(op_,input0[ScX ? 0 : el],input1[ScY ? 0 : el],f);
        casadi_math<T>::der(op_,input0[ScX ? 0 : el],input1[ScY ? 0 : el],f,pd);
        output0[el] = f;

        // Propagate forward seeds
        for(int d=0; d<nfwd; ++d){
          fwdSens[d][0]->data()[el] = pd[0]*fwdSeed[d][0]->data()[ScX ? 0 : el] + pd[1]*fwdSeed[d][1]->data()[ScY ? 0 : el];
        }
    
        // Propagate adjoint seeds
        for(int d=0; d<nadj; ++d){
          T s = adjSeed[d][0]->data()[el];
          adjSeed[d][0]->data()[el] = 0;
          adjSens[d][0]->data()[ScX ? 0 : el] += s*pd[0];
          adjSens[d][1]->data()[ScY ? 0 : el] += s*pd[1];
        }
      }
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX,ScY>::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    bvec_t *input0 = get_bvec_t(input[0]->data());
    bvec_t *input1 = get_bvec_t(input[1]->data());
    bvec_t *outputd = get_bvec_t(output[0]->data());
    for(int el=0; el<output[0]->size(); ++el){
      if(fwd){
        outputd[el] = input0[ScX ? 0 : el] | input1[ScY ? 0 : el];
      } else {
        bvec_t s = outputd[el];
        outputd[el] = bvec_t(0);
        input0[ScX ? 0 : el] |= s;
        input1[ScY ? 0 : el] |= s;
      }
    }
  }

  template<bool ScX, bool ScY>
  MX BinaryMX<ScX,ScY>::getUnary(int op) const{
    switch(op_){
    default: break; // no rule
    }

    // Fallback to default implementation
    return MXNode::getUnary(op);
  }

  template<bool ScX, bool ScY>
  MX BinaryMX<ScX,ScY>::getBinary(int op, const MX& y, bool scX, bool scY) const{
    switch(op_){
    case OP_ADD:
      if(op==OP_SUB && y.isEqual(dep(0),maxDepth())) return dep(1);
      if(op==OP_SUB && y.isEqual(dep(1),maxDepth())) return dep(0);
      break;
    case OP_SUB:
      if(op==OP_SUB && y.isEqual(dep(0),maxDepth())) return -dep(1);
      if(op==OP_ADD && y.isEqual(dep(1),maxDepth())) return dep(0);
      break;
    default: break; // no rule
    }

    // Fallback to default implementation
    return MXNode::getBinary(op,y,scX,scY);    
  }


} // namespace CasADi

#endif // BINARY_MX_IMPL_HPP
