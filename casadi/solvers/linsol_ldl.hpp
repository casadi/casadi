/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#ifndef CASADI_LINSOL_LDL_HPP
#define CASADI_LINSOL_LDL_HPP

/** \defgroup plugin_Linsol_ldl Title
    \par

  * Linear solver using sparse direct LDL factorization

    \identifier{233} */

/** \pluginsection{Linsol,ldl} */

/// \cond INTERNAL
#include "casadi/core/linsol_internal.hpp"
#include <casadi/solvers/casadi_linsol_ldl_export.h>

namespace casadi {
  struct CASADI_LINSOL_LDL_EXPORT LinsolLdlMemory : public LinsolMemory {
    std::vector<double> l, d, w;
  };

  /** \brief \pluginbrief{LinsolInternal,ldl}
   * @copydoc LinsolInternal_doc
   * @copydoc plugin_LinsolInternal_ldl
   */
  class CASADI_LINSOL_LDL_EXPORT LinsolLdl : public LinsolInternal {
  public:

    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LinsolLdl(const std::string& name, const Sparsity& sp);

    /** \brief  Create a new LinsolInternal */
    static LinsolInternal* creator(const std::string& name, const Sparsity& sp) {
      return new LinsolLdl(name, sp);
    }

    // Destructor
    ~LinsolLdl() override;

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new LinsolLdlMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<LinsolLdlMemory*>(mem);}

    // Symbolic factorization
    int sfact(void* mem, const double* A) const override;

    // Factorize the linear system
    int nfact(void* mem, const double* A) const override;

    // Solve the linear system
    int solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const override;

    /// Generate C code
    void generate(CodeGenerator& g, const std::string& A, const std::string& x,
                  casadi_int nrhs, bool tr) const override;

    /// Number of negative eigenvalues
    casadi_int neig(void* mem, const double* A) const override;

    /// Matrix rank
    casadi_int rank(void* mem, const double* A) const override;

    /// A documentation string
    static const std::string meta_doc;

    // Get name of the plugin
    const char* plugin_name() const override { return "ldl";}

    // Get name of the class
    std::string class_name() const override { return "LinsolLdl";}

    // Symbolic factorization
    std::vector<casadi_int> p_;
    Sparsity sp_Lt_;

    ///@{
    // Options
    bool incomplete_, amd_;
    ///@}

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new LinsolLdl(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit LinsolLdl(DeserializingStream& s);
  };

} // namespace casadi

/// \endcond

#endif // CASADI_LINSOL_LDL_HPP
