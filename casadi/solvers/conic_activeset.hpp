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


#ifndef CASADI_CONIC_ACTIVESET_HPP
#define CASADI_CONIC_ACTIVESET_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/solvers/casadi_conic_activeset_export.h>


/** \defgroup plugin_Conic_activeset
 Solve QPs using an active-set method
*/

/** \pluginsection{Conic,activeset} */

/// \cond INTERNAL
namespace casadi {

  template<typename T1>
  struct casadi_qp_prob {
    // Dimensions
    casadi_int nx, na, nz;
    // Smallest nonzero number
    T1 dmin;
    // Infinity
    T1 inf;
    // Dual to primal error
    T1 du_to_pr;
    // Print iterations
    int print_iter;
    // Sparsity patterns
    const casadi_int *sp_a, *sp_h, *sp_at, *sp_kkt;
    // Symbolic QR factorization
    const casadi_int *prinv, *pc, *sp_v, *sp_r;
  };

  template<typename T1>
  struct casadi_qp_data {
    // Problem structure
    const casadi_qp_prob<T1>* prob;
    // Cost
    T1 f;
    // QP data
    const T1 *nz_a, *nz_h, *g;
    // Vectors
    T1 *z, *lbz, *ubz, *infeas, *tinfeas, *lam, *w, *dz, *dlam;
    casadi_int *iw, *neverzero, *neverlower, *neverupper;
    // Numeric QR factorization
    T1 *nz_at, *nz_kkt, *beta, *nz_v, *nz_r;
    // Message buffer
    char msg[40];
    // Stepsize
    T1 tau;
    // Singularity
    casadi_int sing;
    // Smallest diagonal value for the QR factorization
    T1 mina;
    casadi_int imina;
    // Primal and dual error, corresponding index
    T1 pr, du;
    casadi_int ipr, idu;
  };

  struct CASADI_CONIC_ACTIVESET_EXPORT ConicActiveSetMemory : public ConicMemory {
  };

  /** \brief \pluginbrief{Conic,activeset}

      @copydoc Conic_doc
      @copydoc plugin_Conic_activeset

      \author Joel Andersson
      \date 2018
  */
  class CASADI_CONIC_ACTIVESET_EXPORT ConicActiveSet : public Conic {
  public:
    /** \brief  Create a new Solver */
    explicit ConicActiveSet(const std::string& name,
                     const std::map<std::string, Sparsity> &st);

    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new ConicActiveSet(name, st);
    }

    /** \brief  Destructor */
    ~ConicActiveSet() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "as";}

    // Get name of the class
    std::string class_name() const override { return "ConicActiveSet";}

    /** \brief Create memory block */
    void* alloc_mem() const override { return new ConicActiveSetMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<ConicActiveSetMemory*>(mem);}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Initialize */
    void init(const Dict& opts) override;

    /** \brief Solve the QP */
    int eval(const double** arg, double** res,
             casadi_int* iw, double* w, void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;
    // Memory structure
    casadi_qp_prob<double> p_;
    // A transpose sparsity
    Sparsity AT_;

    // KKT matrix sparsity
    Sparsity kkt_;

    // QR factorization
    std::vector<casadi_int> prinv_, pc_;
    Sparsity sp_v_, sp_r_;

    ///@{
    // Options
    casadi_int max_iter_;
    double tol_;
    bool print_iter_, print_header_;
    double du_to_pr_;
    ///@}

  };

} // namespace casadi
/// \endcond
#endif // CASADI_CONIC_ACTIVESET_HPP
