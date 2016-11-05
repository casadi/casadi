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


#ifndef CASADI_SLICOT_DPLE_HPP
#define CASADI_SLICOT_DPLEL_HPP

#include "../../core/function/dple_impl.hpp"
#include "../../core/function/linsol.hpp"
#include <casadi/interfaces/slicot/casadi_dple_slicot_export.h>

/** \defgroup plugin_Dple_slicot
 *
 * An efficient solver for Discrete Periodic Lyapunov Equations using SLICOT

 * Uses Periodic Schur Decomposition ('psd') and does not assume positive definiteness.
 * Based on Periodic Lyapunov equations: some applications and new algorithms.
 * Int. J. Control, vol. 67, pp. 69-87, 1997.
 *
 * Overview of the method:
 *   J. Gillis
 *   Practical Methods for Approximate Robust Periodic Optimal Control ofNonlinear Mechanical Systems,
 *   PhD Thesis, KULeuven, 2015
*/

/** \pluginsection{Dple,slicot} */

/// \cond INTERNAL
namespace casadi {


  // Forward declaration
  class SlicotDple;

  struct CASADI_DPLE_SLICOT_EXPORT SlicotDpleMemory {

    /// T Hessenberg-triangular data
    /// Z Schur form multiplier data
    /// X Schur form multiplier data
    // Xbar Schur form multiplier data
    /// VZ Transformed V data
    /// nnKa Temp data  (n x n) x K
    /// nnKb Temp data  (n x n) x K
    /// eig_real Real parts of eigenvalues
    /// eig_imag Imaginary parts of eigenvalues

    /// Temp data  F
    /// Temp data  FF
    /// dwork Work vector for periodic Schur form

    double *VZ, *T, *Z, *X, *Xbar, *nnKa, *nnKb, *eig_real, *eig_imag, *F, *FF, *A, *B, *dwork;
    int* partition;

    /// Solvers for low-order Discrete Periodic Sylvester Equations
    std::vector< std::vector< Linsol> > dpse_solvers;

    /// Constructor
    SlicotDpleMemory() {}

    /// Destructor
    ~SlicotDpleMemory() {}
  };

  /** \brief \pluginbrief{Dple,slicot}
   *
   * An efficient solver for Discrete Periodic Lyapunov Equations using SLICOT
   *
   * @copydoc Dple_doc
   * @copydoc plugin_Dple_slicot

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DPLE_SLICOT_EXPORT SlicotDple : public Dple {
  public:
    /** \brief  Constructor */
    explicit SlicotDple();

    /** \brief  Constructor
     * \param st \structargument{Dple}
     */
    SlicotDple(const std::string& name, const SpDict & st);

    /** \brief  Create a new QP Solver */
    static Dple* creator(const std::string& name,
                          const SpDict& st) {
      return new SlicotDple(name, st);
    }

    /** \brief  Destructor */
    virtual ~SlicotDple();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "slicot";}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new SlicotDpleMemory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<SlicotDpleMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    /** \brief  Evaluate numerically */
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /// A documentation string
    static const std::string meta_doc;

    SlicotDple(const SpDict & st);

  private:
    /// Dimension of state-space
    int n_;

    int partindex(const SlicotDpleMemory* m, int i, int j, int k, int r, int c) const;

    /// Numerical zero, used in periodic Schur form
    double psd_num_zero_;

    /// Linear solver name
    std::string linear_solver_;

    /// Options to be passed to linear solver constructor
    Dict linear_solver_options_;

  };

  void slicot_mb03vd(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2, double * tau,
                     int ldtau, double * dwork=0);

  void slicot_mb03vy(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2,
                     const double * tau, int ldtau, double * dwork=0, int ldwork=0);

  void slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi, int iloz, int ihiz,
                     double *h, int ldh1, int ldh2, double* z, int ldz1, int ldz2, double* wr,
                     double *wi, double * dwork=0, int ldwork=0);

  void slicot_periodic_schur(int n, int K, const double* a,
                             double* t,  double * z,
                             double* dwork, double* eig_real,
                             double *eig_imag, double num_zero=0);


} // namespace casadi

/// \endcond
#endif // CASADI_SLICOT_DPLE_HPP
