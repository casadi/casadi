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


#ifndef CASADI_MULTIPLICATION_CPP
#define CASADI_MULTIPLICATION_CPP

#include "multiplication.hpp"
#include "casadi_misc.hpp"
#include "function_internal.hpp"
#include "serializing_stream.hpp"

using namespace std;

namespace casadi {

  Multiplication::Multiplication(const MX& z, const MX& x, const MX& y, const Dict& opts) {
    casadi_assert(x.size2() == y.size1() && x.size1() == z.size1()
      && y.size2() == z.size2(),
      "Multiplication::Multiplication: dimension mismatch. Attempting to multiply "
      + x.dim() + " with " + y.dim()
      + " and add the result to " + z.dim());

    //if (x->is_zero()) {}
    if (z.is_zero()) {
      set_dep(x, y); // z must be first (cfr n_inplace)
    } else {
      set_dep(z, x, y); // z must be first (cfr n_inplace)
    }
    set_sparsity(z.sparsity());
  }

  std::string Multiplication::disp(const std::vector<std::string>& arg) const {
    if (n_dep()==3) {
      return "mac(" + arg.at(1) + "," + arg.at(2) + "," + arg.at(0) + ")";
    } else {
      return "mtimes(" + arg.at(0) + "," + arg.at(1) + ")";
    }
  }

  int Multiplication::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Multiplication::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }


  template<typename T>
  int Multiplication::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    if (n_dep()==3) {
      if (arg[0]!=res[0]) copy(arg[0], arg[0]+dep(0).nnz(), res[0]);
      casadi_mtimes(arg[1], dep(1).sparsity(),
                arg[2], dep(2).sparsity(),
                res[0], sparsity(), w, false);
    } else {
      casadi_clear(res[0], nnz());
      casadi_mtimes(arg[0], dep(0).sparsity(),
                arg[1], dep(1).sparsity(),
                res[0], sparsity(), w, false);
    }
    return 0;
  }


  void Multiplication::ad_forward(const std::vector<std::vector<MX> >& fseed,
                               std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      if (n_dep()==3) {
        fsens[d][0] = fseed[d][0]
          + mac(dep(1), fseed[d][2], MX::zeros(dep(0).sparsity()))
          + mac(fseed[d][1], dep(2), MX::zeros(dep(0).sparsity()));
      } else {
        fsens[d][0] = mac(dep(0), fseed[d][1], MX::zeros(sparsity()))
          + mac(fseed[d][0], dep(1), MX::zeros(sparsity()));
      }
    }
  }

  void Multiplication::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                               std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      if (n_dep()==3) {
        asens[d][1] += mac(aseed[d][0], dep(2).T(), MX::zeros(dep(1).sparsity()));
        asens[d][2] += mac(dep(1).T(), aseed[d][0], MX::zeros(dep(2).sparsity()));
        asens[d][0] += aseed[d][0];
      } else {
        asens[d][0] += mac(aseed[d][0], dep(1).T(), MX::zeros(dep(0).sparsity()));
        asens[d][1] += mac(dep(0).T(), aseed[d][0], MX::zeros(dep(1).sparsity()));
      }
    }
  }

  void Multiplication::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    if (n_dep()==3) {
      res[0] = mac(arg[1], arg[2], arg[0]);
    } else {
      res[0] = mtimes(arg[0], arg[1]);
    }
  }

  int Multiplication::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    if (n_dep()==3) {
      copy_fwd(arg[0], res[0], nnz());
      Sparsity::mul_sparsityF(arg[1], dep(1).sparsity(),
                             arg[2], dep(2).sparsity(),
                             res[0], sparsity(), w);
    } else {
      fill_n(res[0], nnz(), 0);
      Sparsity::mul_sparsityF(arg[0], dep(0).sparsity(),
                              arg[1], dep(1).sparsity(),
                              res[0], sparsity(), w);
    }
    return 0;
  }

  int Multiplication::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    if (n_dep()==3) {
      Sparsity::mul_sparsityR(arg[1], dep(1).sparsity(),
                             arg[2], dep(2).sparsity(),
                             res[0], sparsity(), w);
      copy_rev(arg[0], res[0], nnz());
    } else {
      Sparsity::mul_sparsityR(arg[0], dep(0).sparsity(),
                             arg[1], dep(1).sparsity(),
                             res[0], sparsity(), w);
      fill_n(res[0], nnz(), 0);
    }
    return 0;
  }

  void Multiplication::generate(CodeGenerator& g,
                                const std::vector<casadi_int>& arg,
                                const std::vector<casadi_int>& res, bool prefer_inline) const {
    if (n_dep()==3) {
      // Copy first argument if not inplace
      if (arg[0]!=res[0]) {
        g << g.copy(g.work(arg[0], nnz()), nnz(), g.work(res[0], nnz())) << '\n';
      }

      // Perform sparse matrix multiplication
     g << g.mtimes(g.work(arg[1], dep(1).nnz()), dep(1).sparsity(),
                          g.work(arg[2], dep(2).nnz()), dep(2).sparsity(),
                          g.work(res[0], nnz()), sparsity(), "w", false) << '\n';
    } else {
      g << g.clear(g.work(res[0], nnz()), nnz()) << '\n';
      // Perform sparse matrix multiplication
      g << g.mtimes(g.work(arg[0], dep(0).nnz()), dep(0).sparsity(),
                            g.work(arg[1], dep(1).nnz()), dep(1).sparsity(),
                            g.work(res[0], nnz()), sparsity(), "w", false) << '\n';
    }

  }

  void DenseMultiplication::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res, bool prefer_inline) const {
    casadi_int nrow_x, nrow_y, ncol_y;
    std::string x, y;
    if (n_dep()==3) {
      // Copy first argument if not inplace
      if (arg[0]!=res[0]) {
        g << g.copy(g.work(arg[0], nnz()), nnz(),
        g.work(res[0], nnz())) << '\n';
      }
      
      nrow_x = dep(1).size1();
      nrow_y = dep(2).size1();
      ncol_y = dep(2).size2();
      x = g.work(arg[1], dep(1).nnz());
      y = g.work(arg[2], dep(2).nnz());
    } else {
      nrow_x = dep(0).size1();
      nrow_y = dep(1).size1();
      ncol_y = dep(1).size2();

      x = g.work(arg[0], dep(0).nnz());
      y = g.work(arg[1], dep(1).nnz());
    }

    casadi_int flops = ncol_y*nrow_x*nrow_y;


    if (g.blasfeo) {//flops>100) {
      /*g << "cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans," << nrow_x << "," << ncol_y << "," << nrow_y << ",1.0,";
      g << g.work(arg[0], dep(0).nnz()) << "," << nrow_x << ",";
      g << g.work(arg[1], dep(1).nnz()) << "," << nrow_y << ",";
      g << "1.0,";
      g << g.work(res[0], dep(2).nnz()) << "," << nrow_x << ");\n";*/
      g.local("ta", "char");
      g.local("tb", "char");
      g.local("blasfeo_m", "int");
      g.local("blasfeo_n", "int");
      g.local("blasfeo_k", "int");
      g.local("blasfeo_lda", "int");
      g.local("blasfeo_ldb", "int");
      g.local("blasfeo_ldc", "int");
      g.local("blasfeo_alpha", "casadi_real");
      g.local("blasfeo_beta", "casadi_real");
      g << "ta = 'N';";
      g << "tb = 'N';";
      g << "blasfeo_m = " << nrow_x << ";";
      g << "blasfeo_n = " << ncol_y << ";";
      g << "blasfeo_k = " << nrow_y << ";";
      g << "blasfeo_lda = " << nrow_x << ";";
      g << "blasfeo_ldb = " << nrow_y << ";";
      g << "blasfeo_ldc = " << nrow_x << ";";
      g << "blasfeo_alpha = 1.0;";
      if (n_dep()==3) {
        g << "blasfeo_beta = 1.0;";
      } else {
        g << "blasfeo_beta = 0.0;";
      }
      g << "blas_" << (g.casadi_real_type=="double"? "d" : "s") << "gemm(&ta,&tb,&blasfeo_m,&blasfeo_n,&blasfeo_k,&blasfeo_alpha,";
      g << "(casadi_real*)" << x << ",&blasfeo_lda,";
      g << "(casadi_real*)" << y << ",&blasfeo_ldb,";
      g << "&blasfeo_beta,";
      g << g.work(res[0], nnz()) << ",&blasfeo_ldc);\n";
    } else {
      g.local("rr", "casadi_real", "*");
      g.local("ss", "const casadi_real", "*");
      g.local("tt", "const casadi_real", "*");
      g.local("i", "casadi_int");
      g.local("j", "casadi_int");
      g.local("k", "casadi_int");
      g << g.clear(g.work(res[0], nnz()), nnz()) << '\n';
      g << "for (i=0, rr=" << g.work(res[0], nnz()) <<"; i<" << ncol_y << "; ++i)"
        << " for (j=0; j<" << nrow_x << "; ++j, ++rr)"
        << " for (k=0, ss=" << x << "+j, tt="
        << y << "+i*" << nrow_y << "; k<" << nrow_y << "; ++k)"
        << " *rr += ss[k*" << nrow_x << "]**tt++;\n";
    }

    /*g << "blasfeo_dgemm_normal('N','N'," << nrow_x << "," << ncol_y << "," << nrow_y << ",1.0,";
    g << g.work(arg[0], dep(0).nnz()) << "," << nrow_x << ",";
    g << g.work(arg[1], dep(1).nnz()) << "," << nrow_y << ",";
    g << "1.0,";
    g << g.work(res[0], dep(2).nnz()) << "," << nrow_x << ");\n";*/
    
    /*
*/
  }

  void Multiplication::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Multiplication::dense", false);
  }

  void DenseMultiplication::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s); // NOLINT
    s.pack("Multiplication::dense", true);
  }

  MXNode* Multiplication::deserialize(DeserializingStream& s) {
    bool dense;
    s.unpack("Multiplication::dense", dense);
    if (dense) {
      return new DenseMultiplication(s);
    } else {
      return new Multiplication(s);
    }

  }

} // namespace casadi

#endif // CASADI_MULTIPLICATION_CPP
