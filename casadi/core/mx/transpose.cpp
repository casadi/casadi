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

  void Transpose::evalD(const cpv_double& input, const pv_double& output,
                            int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

 void DenseTranspose::evalD(const cpv_double& input, const pv_double& output,
                                int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  void Transpose::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                             int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  void DenseTranspose::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                                  int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<typename T>
  void Transpose::evalGen(const std::vector<const T*>& input,
                          const std::vector<T*>& output, int* itmp, T* rtmp) {

    // Get sparsity patterns
    //const vector<int>& x_colind = input[0]->colind();
    const int* x_row = dep(0).row();
    int x_sz = dep(0).nnz();
    const int* xT_colind = sparsity().colind();
    int xT_ncol = sparsity().size2();

    const T* x = input[0];
    T* xT = output[0];

    // Transpose
    copy(xT_colind, xT_colind+xT_ncol+1, itmp);
    for (int el=0; el<x_sz; ++el) {
      xT[itmp[x_row[el]]++] = x[el];
    }
  }

  template<typename T>
  void DenseTranspose::evalGen(const std::vector<const T*>& input,
                               const std::vector<T*>& output, int* itmp, T* rtmp) {
    // Get sparsity patterns
    int x_nrow = dep().size1();
    int x_ncol = dep().size2();

    const T* x = input[0];
    T* xT = output[0];
    for (int i=0; i<x_ncol; ++i) {
      for (int j=0; j<x_nrow; ++j) {
        xT[i+j*x_ncol] = x[j+i*x_nrow];
      }
    }
  }

  void Transpose::spFwd(const cpv_bvec_t& arg,
                        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    // Shortands
    const bvec_t *x = arg[0];
    bvec_t *xT = res[0];

    // Get sparsity
    int nz = nnz();
    const int* x_row = dep().row();
    const int* xT_colind = sparsity().colind();
    int xT_ncol = sparsity().size2();

    // Loop over the nonzeros of the argument
    copy(xT_colind, xT_colind+xT_ncol+1, itmp);
    for (int el=0; el<nz; ++el) {
      xT[itmp[*x_row++]++] = *x++;
    }
  }

  void Transpose::spAdj(const pv_bvec_t& arg,
                        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    // Shortands
    bvec_t *x = arg[0];
    bvec_t *xT = res[0];

    // Get sparsity
    int nz = nnz();
    const int* x_row = dep().row();
    const int* xT_colind = sparsity().colind();
    int xT_ncol = sparsity().size2();

    // Loop over the nonzeros of the argument
    copy(xT_colind, xT_colind+xT_ncol+1, itmp);
    for (int el=0; el<nz; ++el) {
      int elT = itmp[*x_row++]++;
      *x++ |= xT[elT];
      xT[elT] = 0;
    }
  }

  void DenseTranspose::spFwd(const cpv_bvec_t& arg,
                             const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    // Shorthands
    const bvec_t *x = arg[0];
    bvec_t *xT = res[0];
    int x_nrow = dep().size1();
    int x_ncol = dep().size2();

    // Loop over the elements
    for (int rr=0; rr<x_nrow; ++rr) {
      for (int cc=0; cc<x_ncol; ++cc) {
        *xT++ = x[rr+cc*x_nrow];
      }
    }
  }

  void DenseTranspose::spAdj(const pv_bvec_t& arg,
                             const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    // Shorthands
    bvec_t *x = arg[0];
    bvec_t *xT = res[0];
    int x_nrow = dep().size1();
    int x_ncol = dep().size2();

    // Loop over the elements
    for (int rr=0; rr<x_nrow; ++rr) {
      for (int cc=0; cc<x_ncol; ++cc) {
        x[rr+cc*x_nrow] |= *xT;
        *xT++ = 0;
      }
    }
  }

  void Transpose::printPart(std::ostream &stream, int part) const {
    if (part!=0) {
      stream << "'";
    }
  }

  void Transpose::eval(const MXPtrV& input, MXPtrV& output) {
    *output[0] = input[0]->T();
  }

  void Transpose::evalFwd(const MXPtrVV& fwdSeed, MXPtrVV& fwdSens) {
    for (int d=0; d<fwdSens.size(); ++d) {
      *fwdSens[d][0] = fwdSeed[d][0]->T();
    }
  }

  void Transpose::evalAdj(MXPtrVV& adjSeed, MXPtrVV& adjSens) {
    for (int d=0; d<adjSeed.size(); ++d) {
      adjSens[d][0]->addToSum(adjSeed[d][0]->T());
      *adjSeed[d][0] = MX();
    }
  }

  void Transpose::generate(std::ostream &stream,
                                    const std::vector<int>& arg,
                                    const std::vector<int>& res,
                                    CodeGenerator& gen) const {
    gen.addAuxiliary(CodeGenerator::AUX_TRANS);

    stream << "  casadi_trans("
           << gen.work(arg[0]) << ", s" << gen.addSparsity(dep().sparsity()) << ", "
           << gen.work(res[0]) << ", s" << gen.addSparsity(sparsity()) << ", iii);" << endl;
  }

  void DenseTranspose::generate(std::ostream &stream,
                                         const std::vector<int>& arg,
                                         const std::vector<int>& res,
                                         CodeGenerator& gen) const {
    stream << "  for (i=0, rr=" << gen.work(res[0]) << ", "
           << "cs=" << gen.work(arg[0]) << "; i<" << dep().size2() << "; ++i) "
           << "for (j=0; j<" << dep().size1() << "; ++j) "
           << "rr[i+j*" << dep().size2() << "] = *cs++;" << endl;
  }

} // namespace casadi
