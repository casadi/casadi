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

#include "transpose.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"

using namespace std;

namespace CasADi{

  Transpose::Transpose(const MX& x){
    setDependencies(x);
    setSparsity(x.sparsity().transpose());
  }

  void Transpose::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,itmp,rtmp);
  }

 void DenseTranspose::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,itmp,rtmp);
  }

  void Transpose::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, std::vector<int>& itmp, std::vector<SX>& rtmp){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,itmp,rtmp);
  }

  void DenseTranspose::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, std::vector<int>& itmp, std::vector<SX>& rtmp){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,itmp,rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void Transpose::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp, std::vector<T>& rtmp){

    // Get sparsity patterns
    //const vector<int>& x_rowind = input[0]->rowind();
    const vector<int>& x_col = input[0]->col();
    const vector<int>& xT_rowind = output[0]->rowind();
    
    const vector<T>& x = input[0]->data();
    vector<T>& xT = output[0]->data();

    // Transpose
    copy(xT_rowind.begin(),xT_rowind.end(),itmp.begin());
    for(int el=0; el<x_col.size(); ++el){
      xT[itmp[x_col[el]]++] = x[el];
    }
  }

  template<typename T, typename MatV, typename MatVV>
  void DenseTranspose::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp, std::vector<T>& rtmp){

    // Get sparsity patterns
    int x_nrow = input[0]->size1();
    int x_ncol = input[0]->size2();
    
    const vector<T>& x = input[0]->data();
    vector<T>& xT = output[0]->data();
    for(int i=0; i<x_nrow; ++i){
      for(int j=0; j<x_ncol; ++j){
        xT[i+j*x_nrow] = x[j+i*x_ncol];
      }
    }
  }
  
  void Transpose::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp, bool fwd){
    // Access the input
    bvec_t *x = get_bvec_t(input[0]->data());
    //const vector<int>& x_rowind = input[0]->rowind();
    const vector<int>& x_col = input[0]->col();

    // Access the output
    bvec_t *xT = get_bvec_t(output[0]->data());
    const vector<int>& xT_rowind = output[0]->rowind();

    // Offset for each row of the result
    copy(xT_rowind.begin(),xT_rowind.end(),itmp.begin());

    // Loop over the nonzeros of the argument
    for(int el=0; el<x_col.size(); ++el){

      // Get the column
      int j = x_col[el];

      // Copy nonzero
      if(fwd){
        xT[itmp[j]++] = x[el];
      } else {
        int elT = itmp[j]++;
        x[el] |= xT[elT];
        xT[elT] = 0;
      }
    }    
  }

  void DenseTranspose::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp, bool fwd){
    // Access the input
    bvec_t *x = get_bvec_t(input[0]->data());
    int x_nrow = input[0]->size1();
    int x_ncol = input[0]->size2();

    // Access the output
    bvec_t *xT = get_bvec_t(output[0]->data());

    // Loop over the elements
    for(int i=0; i<x_nrow; ++i){
      for(int j=0; j<x_ncol; ++j){
        int el = j+i*x_ncol;
        int elT = i+j*x_nrow;
        if(fwd){
          xT[elT] = x[el];
        } else {
          x[el] |= xT[elT];
          xT[elT] = 0;
        }
      }
    }
  }

  void Transpose::printPart(std::ostream &stream, int part) const{
    if(part==0){
    } else {
      stream << "'";
    }
  }

  void Transpose::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    if(!output_given)
      *output[0] = trans(*input[0]);

    // Forward sensitivities
    int nfwd = fwdSens.size();
    for(int d=0; d<nfwd; ++d){
      *fwdSens[d][0] = trans(*fwdSeed[d][0]);
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for(int d=0; d<nadj; ++d){
      *adjSens[d][0] += trans(*adjSeed[d][0]);
      *adjSeed[d][0] = MX();
    }
  }

  void Transpose::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    gen.addAuxiliary(CodeGenerator::AUX_TRANS);

    stream << "  casadi_trans(";
    stream << arg.front() << ",s" << gen.getSparsity(dep().sparsity()) << ",";
    stream << res.front() << ",s" << gen.getSparsity(sparsity()) << ",iii);" << endl;
  }

  void DenseTranspose::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    stream << "  for(i=0; i<" << dep().size1() << "; ++i) ";
    stream << "for(j=0; j<" << dep().size2() << "; ++j) ";
    stream << res.front() << "[i+j*" << dep().size1() << "] = " << arg.front() << "[j+i*" << dep().size2() << "];" << endl;
  }

} // namespace CasADi
