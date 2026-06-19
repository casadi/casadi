/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "blas_blasfeo.hpp"
#include "casadi/core/code_generator.hpp"
#include "casadi/core/exception.hpp"

#include <climits>

namespace casadi {

  static void blasfeo_dgemm(int transa, int transb,
                            casadi_int m, casadi_int n, casadi_int k,
                            double alpha,
                            const double* A, casadi_int lda,
                            const double* B, casadi_int ldb,
                            double beta,
                            double* C, casadi_int ldc) {
    casadi_assert(m   <= INT_MAX && n   <= INT_MAX && k   <= INT_MAX
               && lda <= INT_MAX && ldb <= INT_MAX && ldc <= INT_MAX,
        "BLAS 'blasfeo' plugin: matrix dimension exceeds 32-bit BLAS ABI.");

    char ta = (transa == CASADI_BLAS_TRANS) ? 'T' : 'N';
    char tb = (transb == CASADI_BLAS_TRANS) ? 'T' : 'N';
    int m_ = static_cast<int>(m), n_ = static_cast<int>(n), k_ = static_cast<int>(k);
    int lda_ = static_cast<int>(lda), ldb_ = static_cast<int>(ldb), ldc_ = static_cast<int>(ldc);

    // BLASFEO header takes non-const pointers; semantics are still read-only
    // for A, B, alpha, beta. const_cast is safe here.
    blasfeo_blas_dgemm(&ta, &tb, &m_, &n_, &k_, &alpha,
                       const_cast<double*>(A), &lda_,
                       const_cast<double*>(B), &ldb_,
                       &beta, C, &ldc_);
  }

  static const char* BLASFEO_DECL =
  "/* BLAS \"blasfeo\" plugin: namespaced Fortran ABI, link with -lblasfeo */\n"
  "extern void blasfeo_blas_dgemm(char* transa, char* transb,\n"
  "                               int* m, int* n, int* k,\n"
  "                               double* alpha,\n"
  "                               double* A, int* lda,\n"
  "                               double* B, int* ldb,\n"
  "                               double* beta,\n"
  "                               double* C, int* ldc);";

  static void blasfeo_codegen_mtimes(CodeGenerator& g,
      const std::string& A, casadi_int m, casadi_int k,
      const std::string& B, casadi_int n, const std::string& C) {
    g.add_external(BLASFEO_DECL);
    g.local("blas_tn", "char");      g.init_local("blas_tn", "'N'");
    g.local("blas_one", "double");   g.init_local("blas_one", "1.0");
    g.local("blas_m", "int");
    g.local("blas_n", "int");
    g.local("blas_k", "int");
    g << "blas_m = " << m << "; blas_n = " << n << "; blas_k = " << k << ";\n";
    g << "blasfeo_blas_dgemm(&blas_tn, &blas_tn, &blas_m, &blas_n, &blas_k, "
         "&blas_one, (double*)" << A << ", &blas_m, (double*)" << B
      << ", &blas_k, &blas_one, " << C << ", &blas_m);\n";
  }

  extern "C" int CASADI_BLAS_BLASFEO_EXPORT
  casadi_register_blas_blasfeo(Blas::Plugin* plugin) {
    plugin->name = "blasfeo";
    plugin->doc = BlasBlasfeo::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->exposed.dgemm = &blasfeo_dgemm;
    plugin->exposed.codegen_mtimes = &blasfeo_codegen_mtimes;
    plugin->options = nullptr;
    plugin->deserialize = nullptr;
    plugin->creator = nullptr;
    return 0;
  }

  extern "C" void CASADI_BLAS_BLASFEO_EXPORT casadi_load_blas_blasfeo() {
    Blas::registerPlugin(casadi_register_blas_blasfeo);
  }

} // namespace casadi
