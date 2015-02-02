/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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

namespace casadi {

  Transpose::Transpose(const MX& x) {
    setDependencies(x);
    setSparsity(x.sparsity().T());
  }

  void Transpose::evaluateD(const DMatrix** input, DMatrix** output,
                            int* itmp, double* rtmp) {
    evaluateGen<double, DMatrix>(input, output, itmp, rtmp);
  }

 void DenseTranspose::evaluateD(const DMatrix** input, DMatrix** output,
                                int* itmp, double* rtmp) {
    evaluateGen<double, DMatrix>(input, output, itmp, rtmp);
  }

  void Transpose::evaluateSX(const SX** input, SX** output,
                             int* itmp, SXElement* rtmp) {
    evaluateGen<SXElement, SX>(input, output, itmp, rtmp);
  }

  void DenseTranspose::evaluateSX(const SX** input, SX** output,
                                  int* itmp, SXElement* rtmp) {
    evaluateGen<SXElement, SX>(input, output, itmp, rtmp);
  }

  template<typename T, typename Mat>
  void Transpose::evaluateGen(const Mat** input, Mat** output, int* itmp, T* rtmp) {

    // Get sparsity patterns
    //const vector<int>& x_colind = input[0]->colind();
    const int* x_row = input[0]->row();
    int x_sz = input[0]->nnz();
    const int* xT_colind = output[0]->colind();
    int xT_ncol = output[0]->size2();

    const vector<T>& x = input[0]->data();
    vector<T>& xT = output[0]->data();

    // Transpose
    copy(xT_colind, xT_colind+xT_ncol+1, itmp);
    for (int el=0; el<x_sz; ++el) {
      xT[itmp[x_row[el]]++] = x[el];
    }
  }

  template<typename T, typename Mat>
  void DenseTranspose::evaluateGen(const Mat** input, Mat** output, int* itmp, T* rtmp) {

    // Get sparsity patterns
    int x_ncol = input[0]->size2();
    int x_nrow = input[0]->size1();

    const vector<T>& x = input[0]->data();
    vector<T>& xT = output[0]->data();
    for (int i=0; i<x_ncol; ++i) {
      for (int j=0; j<x_nrow; ++j) {
        xT[i+j*x_ncol] = x[j+i*x_nrow];
      }
    }
  }

  void Transpose::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                                    std::vector<double>& rtmp, bool fwd) {
    // Access the input
    bvec_t *x = get_bvec_t(input[0]->data());
    //const int* x_colind = input[0]->colind();
    const int* x_row = input[0]->row();
    int x_sz = input[0]->nnz();

    // Access the output
    bvec_t *xT = get_bvec_t(output[0]->data());
    const int* xT_colind = output[0]->colind();
    int xT_ncol = output[0]->size2();

    // Offset for each col of the result
    copy(xT_colind, xT_colind+xT_ncol+1, itmp.begin());

    // Loop over the nonzeros of the argument
    for (int el=0; el<x_sz; ++el) {

      // Get the row
      int j = x_row[el];

      // Copy nonzero
      if (fwd) {
        xT[itmp[j]++] = x[el];
      } else {
        int elT = itmp[j]++;
        x[el] |= xT[elT];
        xT[elT] = 0;
      }
    }
  }

  void DenseTranspose::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output,
                                         std::vector<int>& itmp,
                                         std::vector<double>& rtmp, bool fwd) {
    // Access the input
    bvec_t *x = get_bvec_t(input[0]->data());
    int x_ncol = input[0]->size2();
    int x_nrow = input[0]->size1();

    // Access the output
    bvec_t *xT = get_bvec_t(output[0]->data());

    // Loop over the elements
    for (int i=0; i<x_ncol; ++i) {
      for (int j=0; j<x_nrow; ++j) {
        int el = j+i*x_nrow;
        int elT = i+j*x_ncol;
        if (fwd) {
          xT[elT] = x[el];
        } else {
          x[el] |= xT[elT];
          xT[elT] = 0;
        }
      }
    }
  }

  void Transpose::printPart(std::ostream &stream, int part) const {
    if (part==0) {
    } else {
      stream << "'";
    }
  }

  void Transpose::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                             MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                             bool output_given) {
    if (!output_given)
      *output[0] = input[0]->T();

    // Forward sensitivities
    int nfwd = fwdSens.size();
    for (int d=0; d<nfwd; ++d) {
      *fwdSens[d][0] = fwdSeed[d][0]->T();
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for (int d=0; d<nadj; ++d) {
      adjSens[d][0]->addToSum(adjSeed[d][0]->T());
      *adjSeed[d][0] = MX();
    }
  }

  void Transpose::generateOperation(std::ostream &stream,
                                    const std::vector<std::string>& arg,
                                    const std::vector<std::string>& res,
                                    CodeGenerator& gen) const {
    gen.addAuxiliary(CodeGenerator::AUX_TRANS);

    stream << "  casadi_trans(";
    stream << arg.front() << ", s" << gen.getSparsity(dep().sparsity()) << ", ";
    stream << res.front() << ", s" << gen.getSparsity(sparsity()) << ", iii);" << endl;
  }

  void DenseTranspose::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                         const std::vector<std::string>& res,
                                         CodeGenerator& gen) const {
    stream << "  for (i=0; i<" << dep().size2() << "; ++i) ";
    stream << "for (j=0; j<" << dep().size1() << "; ++j) ";
    stream << res.front() << "[i+j*" << dep().size2() << "] = " << arg.front()
           << "[j+i*" << dep().size1() << "];" << endl;
  }

} // namespace casadi
