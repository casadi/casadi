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


#ifndef CASADI_QP_TO_NLP_HPP
#define CASADI_QP_TO_NLP_HPP

#include "casadi/core/function/qpsol.hpp"
#include "casadi/core/function/nlpsol.hpp"

#include <casadi/solvers/casadi_qpsol_nlpsol_export.h>


/** \defgroup plugin_Qpsol_nlp
   Solve QPs using an Nlpsol
*/

/** \pluginsection{Qpsol,nlp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Qpsol,nlp}

      @copydoc Qpsol_doc
      @copydoc plugin_Qpsol_nlp

      \author Joris Gillis
      \date 2011
  */
  class CASADI_QPSOL_NLPSOL_EXPORT QpToNlp : public Qpsol {
  public:
    /** \brief  Create a new Solver */
    explicit QpToNlp(const std::string& name,
                     const std::map<std::string, Sparsity> &st);

    /** \brief  Create a new QP Solver */
    static Qpsol* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new QpToNlp(name, st);
    }

    /** \brief  Destructor */
    virtual ~QpToNlp();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "nlpsol";}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief  Initialize */
    virtual void init(const Dict& opts);

    virtual void eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const;

    /// A documentation string
    static const std::string meta_doc;

    /// Solve with
    Function solver_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_QP_TO_NLP_HPP
