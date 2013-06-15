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
    int nx = output.size();

    // Nondifferentiated outputs and forward sensitivities
    for(int d=-1; d<nfwd; ++d){      
      const MatV& arg = d<0 ? input : fwdSeed[d];
      MatV& res = d<0 ? output : fwdSens[d];
      typename vector<T>::const_iterator arg_it = arg[nx]->data().begin();
      for(int i=0; i<nx; ++i){
        vector<T>& res_i = res[i]->data();
        copy(arg[i]->begin(),arg[i]->end(),res_i.begin());
        transform(res_i.begin(),res_i.end(),arg_it,res_i.begin(),std::plus<T>());
        arg_it += res_i.size();
      }
    }
    
    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      typename vector<T>::iterator asens_it = adjSens[d][nx]->data().begin();
      for(int i=0; i<nx; ++i){
        vector<T>& aseed_i = adjSeed[d][i]->data();
        vector<T>& asens_i = adjSens[d][i]->data();
        transform(aseed_i.begin(),aseed_i.end(),asens_i.begin(),asens_i.begin(),std::plus<T>());
        transform(aseed_i.begin(),aseed_i.end(),asens_it,asens_it,std::plus<T>());
        fill(aseed_i.begin(), aseed_i.end(), 0);
        asens_it += aseed_i.size();
      }
    }
  }

  void Vertsplit::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    int nx = output.size();
    bvec_t *arg_ptr = get_bvec_t(input[nx]->data());
    for(int i=0; i<nx; ++i){
      vector<double>& arg_i = input[i]->data();
      vector<double>& res_i = output[i]->data();
      bvec_t *arg_i_ptr = get_bvec_t(arg_i);
      bvec_t *res_i_ptr = get_bvec_t(res_i);
      for(int k=0; k<arg_i.size(); ++k){
        if(fwd){        
          *res_i_ptr++ = *arg_i_ptr++ | *arg_ptr++;
        } else {
          *arg_i_ptr++ |= *res_i_ptr;
          *arg_ptr++ |= *res_i_ptr;          
          *res_i_ptr++ = 0;
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
