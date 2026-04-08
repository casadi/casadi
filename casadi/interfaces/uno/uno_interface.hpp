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


#ifndef CASADI_UNO_INTERFACE_HPP
#define CASADI_UNO_INTERFACE_HPP

#include "uno_nlp.hpp"
#include <casadi/interfaces/uno/casadi_nlpsol_uno_export.h>

#include "casadi/core/nlpsol_impl.hpp"

/** \defgroup plugin_Nlpsol_uno Title
    \par

    David Kiessling

  Uno interface

    \identifier{22c} */

/** \pluginsection{Nlpsol,uno} */

/// \cond INTERNAL

namespace casadi {
  // Forward declaration
  class UnoInterface;

  /*------------------------
  Definition of UnoMemory
  ------------------------*/

  struct CASADI_NLPSOL_UNO_EXPORT UnoMemory : public NlpsolMemory {
    const UnoInterface& self;

    void* model;
    void* solver;
    void* uno_nlp;
    const char* return_status;
    /// Constructor
    UnoMemory(const UnoInterface& uno_interface);

    /// Destructor
    ~UnoMemory();
  };

  /*------------------------------
  Definition of class UnoInterface
  -------------------------------*/

  /** \brief \pluginbrief{Nlpsol,uno}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_uno
  */
  class CASADI_NLPSOL_UNO_EXPORT UnoInterface : public Nlpsol {
    friend class UnoNlp;
  public:

    explicit UnoInterface(const std::string& name, const Function& nlp);
    ~UnoInterface() override;

    // Hessian Sparsity
   Sparsity hesslag_sp_;

   // Jacobian sparsity
   Sparsity jacg_sp_;

    // Get name of the plugin
    const char* plugin_name() const override { return "uno";}

    // Get name of the class
    std::string class_name() const override { return "UnoInterface";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new UnoInterface(name, nlp);
    }

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new UnoMemory(*this);}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<UnoMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // UNO options
    Dict opts_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_UNO_INTERFACE_HPP