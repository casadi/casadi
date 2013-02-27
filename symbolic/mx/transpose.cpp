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

  Transpose* Transpose::clone() const{
    return new Transpose(*this);
  }

  void Transpose::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens, std::vector<int>& itmp, std::vector<double>& rtmp){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens,itmp,rtmp);
  }

  void Transpose::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens, std::vector<int>& itmp, std::vector<SX>& rtmp){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens,itmp,rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void Transpose::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens, std::vector<int>& itmp, std::vector<T>& rtmp){
    // Access the input
    const vector<T>& x = input[0]->data();
    const vector<int>& x_rowind = input[0]->rowind();
    const vector<int>& x_col = input[0]->col();

    // Access the output
    vector<T>& xT = output[0]->data();
    const vector<int>& xT_rowind = output[0]->rowind();

    // Offset for each row of the result
    copy(xT_rowind.begin(),xT_rowind.end()-1,itmp.begin());

    // Loop over the rows of the input
    for(int i=0; i<x_rowind.size()-1; ++i){

      // Loop over the nonzeros
      for(int el=x_rowind[i]; el<x_rowind[i+1]; ++el){

	// Get the column
	int j = x_col[el];

	// Copy nonzero
	xT[itmp[j]++] = x[el];
      }      
    }
  }

  void Transpose::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp, bool fwd){
    // Access the input
    bvec_t *x = get_bvec_t(input[0]->data());
    const vector<int>& x_rowind = input[0]->rowind();
    const vector<int>& x_col = input[0]->col();

    // Access the output
    bvec_t *xT = get_bvec_t(output[0]->data());
    const vector<int>& xT_rowind = output[0]->rowind();

    // Offset for each row of the result
    copy(xT_rowind.begin(),xT_rowind.end()-1,itmp.begin());

    // Loop over the rows of the input
    for(int i=0; i<x_rowind.size()-1; ++i){

      // Loop over the nonzeros
      for(int el=x_rowind[i]; el<x_rowind[i+1]; ++el){

	// Get the column
	int j = x_col[el];

	// Copy nonzero
	if(fwd){
	  xT[itmp[j]++] = x[el];
	} else {
	  x[el] |= xT[itmp[j]++];
	}
      }      
    }
  }

  void Transpose::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "trans(";
    } else {
      stream << ")";
    }
  }

  void Transpose::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    *output[0] = trans(*input[0]);
  }

  void Transpose::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    gen.addAuxiliary(CodeGenerator::AUX_TRANS);

    stream << "  casadi_trans(";
    stream << arg.front() << ",s" << gen.getSparsity(dep().sparsity()) << ",";
    stream << res.front() << ",s" << gen.getSparsity(sparsity()) << ",itmp);" << endl;
  }

} // namespace CasADi
