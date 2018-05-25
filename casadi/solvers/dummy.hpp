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


#ifndef CASADI_DUMMY_HPP
#define CASADI_DUMMY_HPP

#include "casadi/core/nlpsol_impl.hpp"
#include <casadi/solvers/casadi_nlpsol_dummy_export.h>

/** \defgroup plugin_Nlpsol_dummy
   Dummy NLP solver, just for showing how to right your own CasADi NLP solver
*/

/** \pluginsection{Nlpsol,dummy} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_DUMMY_EXPORT DummyMemory : public NlpsolMemory {
    /// Memory for the solver: variables that are changed during solve()

    /// Gradient of the objective function
    double *gf;

    /// Current Jacobian
    double *Jk;

    // Statistics
    int iter_count;
  };

  /** \brief  \pluginbrief{Nlpsol,dummy}
  *  @copydoc NLPSolver_doc
  *  @copydoc plugin_Nlpsol_dummy
  */
  class CASADI_NLPSOL_DUMMY_EXPORT Dummy : public Nlpsol {
  public:
    explicit Dummy(const std::string& name, const Function& nlp);
    ~Dummy() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "dummy";}

    // Name of the class
    std::string class_name() const override { return "Dummy";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new Dummy(name, nlp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new DummyMemory();}

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<DummyMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    // Option values
    std::string my_string_;
    int my_int_;
    double my_double_;

    // Variables that are not changed during solve()
    Sparsity J_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_DUMMY_HPP
