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

#include "split.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"
#include "../casadi_options.hpp"

using namespace std;

namespace CasADi{

  Split::Split(const MX& x, const std::vector<int>& offset) : offset_(offset){
    setDependencies(x);
    setSparsity(Sparsity::scalar());
  }

  Split::~Split(){
  }
  
  void Split::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,itmp,rtmp);
  }

  void Split::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp, std::vector<SXElement>& rtmp){
    evaluateGen<SXElement,SXPtrV,SXPtrVV>(input,output,itmp,rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void Split::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp, std::vector<T>& rtmp){
    // Number of derivatives
    int nx = offset_.size()-1;

    const MatV& arg = input;
    MatV& res = output;
    for(int i=0; i<nx; ++i){
      int nz_first = offset_[i];
      int nz_last = offset_[i+1];
      if(res[i]!=0){
        copy(arg[0]->begin()+nz_first, arg[0]->begin()+nz_last, res[i]->begin());
      }
    }
  }

  void Split::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    int nx = offset_.size()-1;
    for(int i=0; i<nx; ++i){
      if(output[i]!=0){
        bvec_t *arg_ptr = get_bvec_t(input[0]->data()) + offset_[i];
        vector<double>& res_i = output[i]->data();
        bvec_t *res_i_ptr = get_bvec_t(res_i);
        for(int k=0; k<res_i.size(); ++k){
          if(fwd){        
            *res_i_ptr++ = *arg_ptr++;
          } else {
            *arg_ptr++ |= *res_i_ptr;          
            *res_i_ptr++ = 0;
          }
        }
      }
    }
  }

  void Split::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    int nx = res.size();
    for(int i=0; i<nx; ++i){
      int nz_first = offset_[i];
      int nz_last = offset_[i+1];
      int nz = nz_last-nz_first;
      if(res.at(i).compare("0")!=0){
        stream << "  for(i=0; i<" << nz << "; ++i) " << res.at(i) << "[i] = " << arg.at(0) << "[i+" << nz_first << "];" << endl;
      }
    }
  }

  Horzsplit::Horzsplit(const MX& x, const std::vector<int>& offset) : Split(x,offset){
    
    // Split up the sparsity pattern
    output_sparsity_ = horzsplit(x.sparsity(),offset_);

    // Have offset_ refer to the nonzero offsets instead of column offsets
    offset_.resize(1);
    for(std::vector<Sparsity>::const_iterator it=output_sparsity_.begin(); it!=output_sparsity_.end(); ++it){
      offset_.push_back(offset_.back() + it->size());
    }
  }

  void Horzsplit::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "horzsplit(";
    } else {
      stream << ")";
    }
  }

  void Horzsplit::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();
    int nx = offset_.size()-1;

    // Get column offsets
    vector<int> col_offset;
    col_offset.reserve(offset_.size());
    col_offset.push_back(0);
    for(std::vector<Sparsity>::const_iterator it=output_sparsity_.begin(); it!=output_sparsity_.end(); ++it){
      col_offset.push_back(col_offset.back() + it->size2());
    }
    
    // Non-differentiated output and forward sensitivities
    int first_d = output_given ? 0 : -1;
    for(int d=first_d; d<nfwd; ++d){
      const MXPtrV& arg = d<0 ? input : fwdSeed[d];
      MXPtrV& res = d<0 ? output : fwdSens[d];
      MX& x = *arg[0];
      vector<MX> y = horzsplit(x,col_offset);
      for(int i=0; i<nx; ++i){
        if(res[i]!=0){
          *res[i] = y[i];
        }
      }
    }

    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      if(adjSens[d][0]!=0){
        vector<MX> v;
        for(int i=0; i<nx; ++i){
          MX* x_i = adjSeed[d][i];          
          if(x_i!=0){
            v.push_back(*x_i);
            *x_i = MX();
          } else {
            v.push_back(MX::sparse(output_sparsity_[i].shape()));
          }
        }
        adjSens[d][0]->addToSum(horzcat(v));
      }
    }
  }

  Vertsplit::Vertsplit(const MX& x, const std::vector<int>& offset) : Split(x,offset){
    
    // Split up the sparsity pattern
    output_sparsity_ = vertsplit(x.sparsity(),offset_);

    // Have offset_ refer to the nonzero offsets instead of column offsets
    offset_.resize(1);
    for(std::vector<Sparsity>::const_iterator it=output_sparsity_.begin(); it!=output_sparsity_.end(); ++it){
      offset_.push_back(offset_.back() + it->size());
    }
  }

  void Vertsplit::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "vertsplit(";
    } else {
      stream << ")";
    }
  }

  void Vertsplit::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();
    int nx = offset_.size()-1;

    // Get row offsets
    vector<int> row_offset;
    row_offset.reserve(offset_.size());
    row_offset.push_back(0);
    for(std::vector<Sparsity>::const_iterator it=output_sparsity_.begin(); it!=output_sparsity_.end(); ++it){
      row_offset.push_back(row_offset.back() + it->size1());
    }
    
    // Non-differentiated output and forward sensitivities
    int first_d = output_given ? 0 : -1;
    for(int d=first_d; d<nfwd; ++d){
      const MXPtrV& arg = d<0 ? input : fwdSeed[d];
      MXPtrV& res = d<0 ? output : fwdSens[d];
      MX& x = *arg[0];
      vector<MX> y = vertsplit(x,row_offset);
      for(int i=0; i<nx; ++i){
        if(res[i]!=0){
          *res[i] = y[i];
        }
      }
    }

    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      if(adjSens[d][0]!=0){
        vector<MX> v;
        for(int i=0; i<nx; ++i){
          MX* x_i = adjSeed[d][i];          
          if(x_i!=0){
            v.push_back(*x_i);
            *x_i = MX();
          } else {
            v.push_back(MX::sparse(output_sparsity_[i].shape()));
          }
        }
        adjSens[d][0]->addToSum(vertcat(v));
      }
    }
  }

  MX Horzsplit::getHorzcat(const std::vector<MX>& x) const{
    // Check x length
    if(x.size()!=getNumOutputs()){
      return MXNode::getHorzcat(x);
    }

    // Check x content
    for(int i=0; i<x.size(); ++i){
      if(!(x[i]->isOutputNode() && x[i]->getFunctionOutput()==i && x[i]->dep().get()==this)){
        return MXNode::getHorzcat(x);
      }
    }

    // OK if reached this point
    return dep();
  }

  MX Vertsplit::getVertcat(const std::vector<MX>& x) const{
    // Check x length
    if(x.size()!=getNumOutputs()){
      return MXNode::getVertcat(x);
    }

    // Check x content
    for(int i=0; i<x.size(); ++i){
      if(!(x[i]->isOutputNode() && x[i]->getFunctionOutput()==i && x[i]->dep().get()==this)){
        return MXNode::getVertcat(x);
      }
    }

    // OK if reached this point
    return dep();
  }

} // namespace CasADi
