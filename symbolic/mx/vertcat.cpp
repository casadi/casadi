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

#include "vertcat.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../fx/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace CasADi{

  Vertcat::Vertcat(const vector<MX>& x){
    setDependencies(x);
    
    // Construct the sparsity
    casadi_assert(!x.empty());
    CRSSparsity sp = x.front().sparsity();
    for(vector<MX>::const_iterator i=x.begin()+1; i!=x.end(); ++i){
      sp.append(i->sparsity());
    }

    setSparsity(sp);
  }

  Vertcat* Vertcat::clone() const{
    return new Vertcat(*this);
  }

  void Vertcat::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void Vertcat::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<typename T, typename MatV, typename MatVV>
  void Vertcat::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){
    // Number of derivatives
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Nondifferentiated outputs and forward sensitivities
    for(int d=-1; d<nfwd; ++d){
      typename vector<T>::iterator res_it = d==-1 ? output[0]->data().begin() : fwdSens[d][0]->data().begin();
      for(int i=0; i<input.size(); ++i){
        const vector<T>& arg_i = d==-1 ? input[i]->data() : fwdSeed[d][i]->data();
        copy(arg_i.begin(),arg_i.end(),res_it);
        res_it += arg_i.size();
      }
    }
    
    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      typename vector<T>::iterator arg_it = adjSeed[d][0]->data().begin();
      for(int i=0; i<input.size(); ++i){
        vector<T>& res_i = adjSens[d][i]->data();
        transform(res_i.begin(),res_i.end(),arg_it,res_i.begin(),std::plus<T>());
        fill_n(arg_it, res_i.size(), 0);
        arg_it += res_i.size();
      }
    }
  }

  void Vertcat::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    bvec_t *res_ptr = get_bvec_t(output[0]->data());
    for(int i=0; i<input.size(); ++i){
      vector<double>& arg_i = input[i]->data();
      bvec_t *arg_i_ptr = get_bvec_t(arg_i);
      if(fwd){
        copy(arg_i_ptr, arg_i_ptr+arg_i.size(), res_ptr);
        res_ptr += arg_i.size();
      } else {
        for(int k=0; k<arg_i.size(); ++k){
          *arg_i_ptr++ |= *res_ptr;
          *res_ptr++ = 0;
        }
      }
    }
  }

  void Vertcat::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "vertcat(";
    } else if(part==ndep()){
      stream << ")";
    } else {
      stream << ",";
    }
  }

  void Vertcat::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    if(!output_given){
      *output[0] = vertcat(getVector(input));
    }
    
    // Forward sensitivities
    int nfwd = fwdSens.size();
    for(int d = 0; d<nfwd; ++d){
      *fwdSens[d][0] = vertcat(getVector(fwdSeed[d]));
    }
    
    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for(int d=0; d<nadj; ++d){
      int row_offset = 0;
      MX& aseed = *adjSeed[d][0];
      for(int i=0; i<input.size(); ++i){
        MX& asens = *adjSens[d][i];
        int nrow = asens.size1();
        asens += aseed(Slice(row_offset,row_offset+nrow),Slice());
        row_offset += nrow;
      }
      casadi_assert(row_offset == aseed.size1());
      aseed = MX();
    }
  }

  void Vertcat::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    int nz_offset = 0;
    for(int i=0; i<arg.size(); ++i){
      int nz = dep(i).size();
      stream << "  for(i=0; i<" << nz << "; ++i) " << res.front() << "[i+" << nz_offset << "] = " << arg.at(i) << "[i];" << endl;
      nz_offset += nz;
    }
    casadi_assert(nz_offset == size());
  }

  MX Vertcat::getGetNonzeros(const CRSSparsity& sp, const std::vector<int>& nz) const{
    // Get the first nonnegative nz
    int nz_test = -1;
    for(vector<int>::const_iterator i=nz.begin(); i!=nz.end(); ++i){
      if(*i>=0){
        nz_test = *i;
        break;
      }
    }

    // Quick return if none
    if(nz_test<0) return MX::zeros(sp);
    
    // Find out to which dependency it might depend
    int begin=0, end=0;
    int i;
    for(i=0; i<ndep(); ++i){
      begin = end;
      end += dep(i).size();
      if(nz_test < end) break;
    }
    
    // Check if any nz refer to a different nonzero
    for(vector<int>::const_iterator j=nz.begin(); j!=nz.end(); ++j){
      if(*j>=0 && (*j < begin || *j >= end)){
        
        // Fallback to the base class
        return MXNode::getGetNonzeros(sp,nz);
      }
    }
    
    // All nz refer to the same dependency, update the nonzero indices
    if(begin==0){
      return dep(i)->getGetNonzeros(sp,nz);
    } else {
      vector<int> nz_new(nz);
      for(vector<int>::iterator j=nz_new.begin(); j!=nz_new.end(); ++j){
        if(*j>=0) *j -= begin;
      }
      return dep(i)->getGetNonzeros(sp,nz_new);
    }
  }

} // namespace CasADi
