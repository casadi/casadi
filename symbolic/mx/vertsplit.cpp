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

#include "vertsplit.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../fx/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace CasADi{

  Vertsplit::Vertsplit(const std::vector<MX>& x, const MX& y){
    std::vector<MX> xy(x);
    xy.push_back(y);
    setDependencies(xy);
    setSparsity(CRSSparsity(1, 1, true));
  }

  Vertsplit* Vertsplit::clone() const{
    return new Vertsplit(*this);
  }

  void Vertsplit::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void Vertsplit::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<typename T, typename MatV, typename MatVV>
  void Vertsplit::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){
    // Number of derivatives
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();

    // Nondifferentiated outputs and forward sensitivities
    for(int d=-1; d<nfwd; ++d){
      casadi_error("not implemented");
      typename vector<T>::iterator res_it = d==-1 ? output[0]->data().begin() : fwdSens[d][0]->data().begin();
      for(int i=0; i<input.size(); ++i){
        const vector<T>& arg_i = d==-1 ? input[i]->data() : fwdSeed[d][i]->data();
        copy(arg_i.begin(),arg_i.end(),res_it);
        res_it += arg_i.size();
      }
    }
    
    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      casadi_error("not implemented");
      typename vector<T>::iterator arg_it = adjSeed[d][0]->data().begin();
      for(int i=0; i<input.size(); ++i){
        vector<T>& res_i = adjSens[d][i]->data();
        transform(res_i.begin(),res_i.end(),arg_it,res_i.begin(),std::plus<T>());
        fill_n(arg_it, res_i.size(), 0);
        arg_it += res_i.size();
      }
    }
  }

  void Vertsplit::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    casadi_error("not implemented");
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

  void Vertsplit::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "vertsplit(";
    } else {
      stream << ")";
    }
  }

  void Vertsplit::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    if(!output_given){
      casadi_error("not implemented");
      //*output[0] = getVertcat(getVector(input));
    }
    
    // Forward sensitivities
    int nfwd = fwdSens.size();
    for(int d = 0; d<nfwd; ++d){
      casadi_error("not implemented");      
      //*fwdSens[d][0] = getVertcat(getVector(fwdSeed[d]));
    }
    
    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for(int d=0; d<nadj; ++d){
      casadi_error("not implemented");
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

} // namespace CasADi
