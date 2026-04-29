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


#include "blas_classic.hpp"
#include "casadi/core/code_generator.hpp"
#include "casadi/core/exception.hpp"

#include <climits>

namespace casadi {

  static void classic_dgemm(int transa, int transb,
                            casadi_int m, casadi_int n, casadi_int k,
                            double alpha,
                            const double* A, casadi_int lda,
                            const double* B, casadi_int ldb,
                            double beta,
                            double* C, casadi_int ldc) {
    // Fortran BLAS-32 takes int by reference. Guard against ILP64 / huge
    // problems before silently truncating.
    casadi_assert(m   <= INT_MAX && n   <= INT_MAX && k   <= INT_MAX
               && lda <= INT_MAX && ldb <= INT_MAX && ldc <= INT_MAX,
        "BLAS 'classic' plugin: matrix dimension exceeds 32-bit BLAS ABI. "
        "Rebuild against an ILP64 BLAS or use the \"reference\" plugin.");

    const char ta = (transa == CASADI_BLAS_TRANS) ? 'T' : 'N';
    const char tb = (transb == CASADI_BLAS_TRANS) ? 'T' : 'N';
    int m_ = static_cast<int>(m), n_ = static_cast<int>(n), k_ = static_cast<int>(k);
    int lda_ = static_cast<int>(lda), ldb_ = static_cast<int>(ldb), ldc_ = static_cast<int>(ldc);

    dgemm_(&ta, &tb, &m_, &n_, &k_, &alpha,
           A, &lda_, B, &ldb_, &beta, C, &ldc_);
  }

  static void classic_daxpy(casadi_int n, double alpha,
                            const double* x, double* y) {
    casadi_assert(n <= INT_MAX,
        "BLAS 'classic' plugin: vector length exceeds 32-bit BLAS ABI.");
    int n_ = static_cast<int>(n), inc = 1;
    daxpy_(&n_, &alpha, x, &inc, y, &inc);
  }

  static double classic_ddot(casadi_int n,
                             const double* x, const double* y) {
    casadi_assert(n <= INT_MAX,
        "BLAS 'classic' plugin: vector length exceeds 32-bit BLAS ABI.");
    int n_ = static_cast<int>(n), inc = 1;
    return ddot_(&n_, x, &inc, y, &inc);
  }

  static void classic_dscal(casadi_int n, double alpha, double* x) {
    casadi_assert(n <= INT_MAX,
        "BLAS 'classic' plugin: vector length exceeds 32-bit BLAS ABI.");
    int n_ = static_cast<int>(n), inc = 1;
    dscal_(&n_, &alpha, x, &inc);
  }

  static double classic_dnrm2(casadi_int n, const double* x) {
    casadi_assert(n <= INT_MAX,
        "BLAS 'classic' plugin: vector length exceeds 32-bit BLAS ABI.");
    int n_ = static_cast<int>(n), inc = 1;
    return dnrm2_(&n_, x, &inc);
  }

  static double classic_dasum(casadi_int n, const double* x) {
    casadi_assert(n <= INT_MAX,
        "BLAS 'classic' plugin: vector length exceeds 32-bit BLAS ABI.");
    int n_ = static_cast<int>(n), inc = 1;
    return dasum_(&n_, x, &inc);
  }

  static const char* CLASSIC_DECL =
  "/* BLAS \"classic\" plugin: Fortran ABI dgemm_, link with -lblas/-lopenblas/-lmkl. */\n"
  "/* Override the symbol at compile time with -DCASADI_BLAS_DGEMM=my_dgemm. */\n"
  "#ifndef CASADI_BLAS_DGEMM\n"
  "#define CASADI_BLAS_DGEMM dgemm_\n"
  "#endif\n"
  "extern void CASADI_BLAS_DGEMM(const char* transa, const char* transb,\n"
  "                              const int* m, const int* n, const int* k,\n"
  "                              const double* alpha,\n"
  "                              const double* A, const int* lda,\n"
  "                              const double* B, const int* ldb,\n"
  "                              const double* beta,\n"
  "                              double* C, const int* ldc);";

  static const char* CLASSIC_DAXPY_DECL =
  "#ifndef CASADI_BLAS_DAXPY\n"
  "#define CASADI_BLAS_DAXPY daxpy_\n"
  "#endif\n"
  "extern void CASADI_BLAS_DAXPY(const int* n, const double* alpha,\n"
  "                              const double* x, const int* incx,\n"
  "                              double* y, const int* incy);";

  static const char* CLASSIC_DDOT_DECL =
  "#ifndef CASADI_BLAS_DDOT\n"
  "#define CASADI_BLAS_DDOT ddot_\n"
  "#endif\n"
  "extern double CASADI_BLAS_DDOT(const int* n,\n"
  "                               const double* x, const int* incx,\n"
  "                               const double* y, const int* incy);";

  static const char* CLASSIC_DSCAL_DECL =
  "#ifndef CASADI_BLAS_DSCAL\n"
  "#define CASADI_BLAS_DSCAL dscal_\n"
  "#endif\n"
  "extern void CASADI_BLAS_DSCAL(const int* n, const double* alpha,\n"
  "                              double* x, const int* incx);";

  static const char* CLASSIC_DNRM2_DECL =
  "#ifndef CASADI_BLAS_DNRM2\n"
  "#define CASADI_BLAS_DNRM2 dnrm2_\n"
  "#endif\n"
  "extern double CASADI_BLAS_DNRM2(const int* n, const double* x, const int* incx);";

  static const char* CLASSIC_DASUM_DECL =
  "#ifndef CASADI_BLAS_DASUM\n"
  "#define CASADI_BLAS_DASUM dasum_\n"
  "#endif\n"
  "extern double CASADI_BLAS_DASUM(const int* n, const double* x, const int* incx);";

  static void classic_codegen_mtimes(CodeGenerator& g,
      const std::string& A, casadi_int m, casadi_int k,
      const std::string& B, casadi_int n, const std::string& C) {
    g.add_external(CLASSIC_DECL);
    // Locals shared with other Fortran-ABI BLAS plugins. g.local is
    // idempotent on (name,type) match, so multiple plugins coexisting in
    // the same Function reuse the same locals safely.
    g.local("blas_tn", "char");      g.init_local("blas_tn", "'N'");
    g.local("blas_one", "double");   g.init_local("blas_one", "1.0");
    g.local("blas_m", "int");
    g.local("blas_n", "int");
    g.local("blas_k", "int");
    g << "blas_m = " << m << "; blas_n = " << n << "; blas_k = " << k << ";\n";
    g << "CASADI_BLAS_DGEMM(&blas_tn, &blas_tn, &blas_m, &blas_n, &blas_k, "
         "&blas_one, " << A << ", &blas_m, " << B << ", &blas_k, "
         "&blas_one, " << C << ", &blas_m);\n";
  }


  static void classic_codegen_axpy_aux(CodeGenerator& g,
      const std::vector<std::string>& inst) {
    static const char* SRC =
        "// SYMBOL \"axpy\"\n"
        "void casadi_axpy(casadi_int n, casadi_real alpha,\n"
        "    const casadi_real* x, casadi_real* y) {\n"
        "  int n_ = (int)n, inc = 1;\n"
        "  CASADI_BLAS_DAXPY(&n_, &alpha, x, &inc, y, &inc);\n"
        "}\n";
    g.auxiliaries << CLASSIC_DAXPY_DECL << "\n";
    g.auxiliaries << g.sanitize_source(SRC, inst);
  }

  static void classic_codegen_scal_aux(CodeGenerator& g,
      const std::vector<std::string>& inst) {
    static const char* SRC =
        "// SYMBOL \"scal\"\n"
        "void casadi_scal(casadi_int n, casadi_real alpha, casadi_real* x) {\n"
        "  int n_ = (int)n, inc = 1;\n"
        "  CASADI_BLAS_DSCAL(&n_, &alpha, x, &inc);\n"
        "}\n";
    g.auxiliaries << CLASSIC_DSCAL_DECL << "\n";
    g.auxiliaries << g.sanitize_source(SRC, inst);
  }

  static void classic_codegen_dot_aux(CodeGenerator& g,
      const std::vector<std::string>& inst) {
    static const char* SRC =
        "// SYMBOL \"dot\"\n"
        "casadi_real casadi_dot(casadi_int n,\n"
        "    const casadi_real* x, const casadi_real* y) {\n"
        "  int n_ = (int)n, inc = 1;\n"
        "  return CASADI_BLAS_DDOT(&n_, x, &inc, y, &inc);\n"
        "}\n";
    g.auxiliaries << CLASSIC_DDOT_DECL << "\n";
    g.auxiliaries << g.sanitize_source(SRC, inst);
  }

  static void classic_codegen_nrm2_aux(CodeGenerator& g,
      const std::vector<std::string>& inst) {
    static const char* SRC =
        "// SYMBOL \"norm_2\"\n"
        "casadi_real casadi_norm_2(casadi_int n, const casadi_real* x) {\n"
        "  int n_ = (int)n, inc = 1;\n"
        "  return CASADI_BLAS_DNRM2(&n_, x, &inc);\n"
        "}\n";
    g.auxiliaries << CLASSIC_DNRM2_DECL << "\n";
    g.auxiliaries << g.sanitize_source(SRC, inst);
  }

  static void classic_codegen_asum_aux(CodeGenerator& g,
      const std::vector<std::string>& inst) {
    // Match casadi_norm_1<T>'s null-pointer semantics: returns 0 on x==NULL.
    static const char* SRC =
        "// SYMBOL \"norm_1\"\n"
        "casadi_real casadi_norm_1(casadi_int n, const casadi_real* x) {\n"
        "  int n_ = (int)n, inc = 1;\n"
        "  if (!x) return 0;\n"
        "  return CASADI_BLAS_DASUM(&n_, x, &inc);\n"
        "}\n";
    g.auxiliaries << CLASSIC_DASUM_DECL << "\n";
    g.auxiliaries << g.sanitize_source(SRC, inst);
  }

  extern "C" int CASADI_BLAS_CLASSIC_EXPORT
  casadi_register_blas_classic(Blas::Plugin* plugin) {
    plugin->name = "classic";
    plugin->doc = BlasClassic::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->exposed.dgemm = &classic_dgemm;
    plugin->exposed.codegen_mtimes = &classic_codegen_mtimes;
    plugin->exposed.daxpy = &classic_daxpy;
    plugin->exposed.ddot  = &classic_ddot;
    plugin->exposed.dscal = &classic_dscal;
    plugin->exposed.dnrm2 = &classic_dnrm2;
    plugin->exposed.dasum = &classic_dasum;
    plugin->exposed.codegen_axpy_aux = &classic_codegen_axpy_aux;
    plugin->exposed.codegen_dot_aux  = &classic_codegen_dot_aux;
    plugin->exposed.codegen_scal_aux = &classic_codegen_scal_aux;
    plugin->exposed.codegen_nrm2_aux = &classic_codegen_nrm2_aux;
    plugin->exposed.codegen_asum_aux = &classic_codegen_asum_aux;
    plugin->options = nullptr;
    plugin->deserialize = nullptr;
    plugin->creator = nullptr;
    return 0;
  }

  extern "C" void CASADI_BLAS_CLASSIC_EXPORT casadi_load_blas_classic() {
    Blas::registerPlugin(casadi_register_blas_classic);
  }

} // namespace casadi
