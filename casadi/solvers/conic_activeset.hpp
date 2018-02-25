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

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    int eval(const double** arg, double** res,
             casadi_int* iw, double* w, void* mem) const override;

    /** Print a vector */
    void print_vector(const char* id, const double* x, casadi_int n) const;

    /** Print an integer vector */
    void print_ivector(const char* id, const casadi_int* x, casadi_int n) const;

    /** Print signs for a vector */
    void print_signs(const char* id, const double* x, casadi_int n) const;

    /// A documentation string
    static const std::string meta_doc;

    // KKT system sparsity
    Sparsity kkt_, AT_;

    // KKT with diagonal
    Sparsity kktd_;

    // QR factorization
    std::vector<casadi_int> prinv_, pc_;
    Sparsity sp_v_, sp_r_;

    ///@{
    // Options
    casadi_int max_iter_;
    double tol_;
    ///@}

  };

} // namespace casadi
/// \endcond
#endif // CASADI_CONIC_ACTIVESET_HPP
