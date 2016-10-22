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


#ifndef CASADI_PSD_INDEF_DPLE_INTERNAL_HPP
#define CASADI_PSD_INDEF_DPLE_INTERNAL_HPP

#include "../../core/function/dple_impl.hpp"
#include <casadi/interfaces/slicot/casadi_dple_slicot_export.h>

/** \defgroup plugin_DpleSolver_slicot
 *
 * An efficient solver for Discrete Periodic Lyapunov Equations using SLICOT

 * Uses Periodic Schur Decomposition ('psd') and does not assume positive definiteness.
 * Based on Periodic Lyapunov equations: some applications and new algorithms.
 * Int. J. Control, vol. 67, pp. 69-87, 1997.
*/

/** \pluginsection{Dple,slicot} */

/// \cond INTERNAL
namespace casadi {


  // Forward declaration
  class SlicotDple;

  struct CASADI_DPLE_SLICOT_EXPORT SlicotDpleMemory {

    /// Constructor
    SlicotDpleMemory();

    /// Destructor
    ~SlicotDpleMemory();
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
  class CASADI_DPLESOLVER_SLICOT_EXPORT SlicotDple : public Dple {
  public:
    /** \brief  Constructor */
    explicit SlicotDple();

    /** \brief  Constructor
     * \param st \structargument{Dple}
     */
    SlicotDple(const std::map<std::string, std::vector<Sparsity> > & st,
                         int nrhs, bool transp);

    /** \brief  Create a new QP Solver */
    static Dple* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new SlicotDple(name, st, 1, false);
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

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    /** \brief  Evaluate numerically */
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /// A documentation string
    static const std::string meta_doc;

    SlicotDple(const std::map<std::string, std::vector<Sparsity> > & st,
                         int nrhs=1, bool transp=false);

    /** \brief  Destructor */
    virtual ~SlicotDple();

    /** \brief  Clone */
    virtual SlicotDple* clone() const;

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function getDerForward(const std::string& name, int nfwd, Dict& opts);
    virtual int numDerForward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function getDerReverse(const std::string& name, int nadj, Dict& opts);
    virtual int numDerReverse() const { return 64;}
    ///@}

  private:
    /// Dimension of state-space
    int n_;

    /// Hessenberg-triangular data
    std::vector<double> T_;

    /// Schur form multiplier data
    std::vector<double> Z_;

    /// Schur form multiplier data
    std::vector<double> X_;

    // Schur form multiplier data
    std::vector<double> Xbar_;

    /// Transformed V data
    std::vector<double> VZ_;

    /// Temp data  (n x n) x K
    std::vector< Matrix<double> > nnKa_;

    /// Temp data  (n x n) x K
    std::vector< Matrix<double> > nnKb_;

    /// Real parts of eigenvalues
    std::vector< double > eig_real_;

    /// Imaginary parts of eigenvalues
    std::vector< double > eig_imag_;

    /// Solvers for low-order Discrete Periodic Sylvester Equations
    std::vector< std::vector< LinearSolver> > dpse_solvers_;

    std::vector<int> partition_;

    int partindex(int i, int j, int k, int r, int c);

    /// Temp data  F
    std::vector< double > F_;

    /// Temp data  FF
    std::vector< double > FF_;

    /// Work vector for periodic Schur form
    std::vector< double > dwork_;

    /// Numerical zero, used in periodic Schur form
    double psd_num_zero_;

  };

  void slicot_mb03vd(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2, double * tau,
                     int ldtau, double * dwork=0);

  void slicot_mb03vy(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2,
                     const double * tau, int ldtau, double * dwork=0, int ldwork=0);

  void slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi, int iloz, int ihiz,
                     double *h, int ldh1, int ldh2, double* z, int ldz1, int ldz2, double* wr,
                     double *wi, double * dwork=0, int ldwork=0);

  void slicot_periodic_schur(int n, int K, const std::vector< double > & a,
                             std::vector< double > & t, std::vector< double > & z,
                             std::vector<double> &eig_real, std::vector<double> &eig_imag,
                             double num_zero=0);

  void slicot_periodic_schur(int n, int K, const std::vector< double > & a,
                             std::vector< double > & t,
                             std::vector< double > & z, std::vector<double> &dwork,
                             std::vector<double> &eig_real, std::vector<double> &eig_imag,
                             double num_zero=0);

  void slicot_periodic_schur(const std::vector< Matrix<double> > & a,
                             std::vector< Matrix<double> > & t,
                             std::vector< Matrix<double> > & z,
                             std::vector< double > & eig_real,
                             std::vector< double > & eig_imag,
                             double num_zero);

} // namespace casadi

/// \endcond
#endif // CASADI_PSD_INDEF_DPLE_INTERNAL_HPP
