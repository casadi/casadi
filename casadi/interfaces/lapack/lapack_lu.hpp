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


#ifndef CASADI_LAPACK_LU_HPP
#define CASADI_LAPACK_LU_HPP

#include "casadi/core/function/linsol_internal.hpp"
#include <casadi/interfaces/lapack/casadi_linsol_lapacklu_export.h>

extern "C" {
  /// LU-Factorize dense matrix (lapack)
  void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

  /// Solve a system of equation using an LU-factorized matrix (lapack)
  void dgetrs_(char* trans, int *n, int *nrhs, double *a,
               int *lda, int *ipiv, double *b, int *ldb, int *info);

  /// Calculate col and row scaling
  void dgeequ_(int *m, int *n, double *a, int *lda, double *r, double *c,
               double *colcnd, double *rowcnd, double *amax, int *info);

  /// Equilibrate the system
  void dlaqge_(int *m, int *n, double *a, int *lda, double *r, double *c,
               double *colcnd, double *rowcnd, double *amax, char *equed);
}

namespace casadi {

/** \defgroup plugin_Linsol_lapacklu
*
   * This class solves the linear system <tt>A.x=b</tt> by making an LU factorization of A: \n
   * <tt>A = L.U</tt>, with L lower and U upper triangular
   *
*/

/** \pluginsection{Linsol,lapacklu} */
/// \cond INTERNAL
  struct CASADI_LINSOL_LAPACKLU_EXPORT LapackLuMemory : public LinsolMemory {
    // Matrix
    std::vector<double> mat;

    /// Pivoting elements
    std::vector<int> ipiv;

    /// Col and row scaling
    std::vector<double> r, c;

    /// Type of scaling during the last equilibration
    char equed;
  };

  /** \brief \pluginbrief{Linsol,lapacklu}
   *
   * @copydoc Linsol_doc
   * @copydoc plugin_Linsol_lapacklu
   *
   */
  class CASADI_LINSOL_LAPACKLU_EXPORT LapackLu : public LinsolInternal {
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LapackLu(const std::string& name, const Sparsity& sparsity);

    /** \brief  Create a new Linsol */
    static LinsolInternal* creator(const std::string& name, const Sparsity& sp) {
      return new LapackLu(name, sp);
    }

    /// Destructor
    virtual ~LapackLu();

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /// Initialize the solver
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new LapackLuMemory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<LapackLuMemory*>(mem);}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    // Factorize the linear system
    virtual void linsol_factorize(void* mem, const double* A) const;

    // Solve the linear system
    virtual void linsol_solve(void* mem, double* x, int nrhs, bool tr) const;

    /// A documentation string
    static const std::string meta_doc;

  protected:

    /// Equilibrate?
    bool equilibriate_;

    /// Allow the equilibration to fail
    bool allow_equilibration_failure_;

    // Get name of the plugin
    virtual const char* plugin_name() const { return "lapacklu";}
  };

/// \endcond

} // namespace casadi

#endif // CASADI_LAPACK_LU_HPP
