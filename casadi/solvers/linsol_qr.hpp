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


#ifndef CASADI_LINSOL_QR_HPP
#define CASADI_LINSOL_QR_HPP

/** \defgroup plugin_Linsol_qr Title
    \par

  * Linear solver using sparse direct QR factorization

    \identifier{22z} */

/** \pluginsection{Linsol,qr} */

/// \cond INTERNAL
#include "casadi/core/linsol_internal.hpp"
#include <casadi/solvers/casadi_linsol_qr_export.h>

namespace casadi {
  struct CASADI_LINSOL_QR_EXPORT LinsolQrMemory : public LinsolMemory {
    std::vector<double> v, r, beta, w;
    std::vector<double> cache;

    // Cache locations sorted by access time
    std::vector<int> cache_loc;
  };

  /** \brief \pluginbrief{LinsolInternal,qr}
   * @copydoc LinsolInternal_doc
   * @copydoc plugin_LinsolInternal_qr
   */
  class CASADI_LINSOL_QR_EXPORT LinsolQr : public LinsolInternal {
  public:

    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LinsolQr(const std::string& name, const Sparsity& sp);

    /** \brief  Create a new LinsolInternal */
    static LinsolInternal* creator(const std::string& name, const Sparsity& sp) {
      return new LinsolQr(name, sp);
    }

    // Destructor
    ~LinsolQr() override;

    // Initialize the solver
    void init(const Dict& opts) override;

    /// Finalize the object creation
    void finalize() override;

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Create memory block */
    void* alloc_mem() const override { return new LinsolQrMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<LinsolQrMemory*>(mem);}

    // Symbolic factorization
    int nfact(void* mem, const double* A) const override;

    // Factorize the linear system
    int sfact(void* mem, const double* A) const override;

    // Solve the linear system
    int solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const override;

    /// Generate C code
    void generate(CodeGenerator& g, const std::string& A, const std::string& x,
                  casadi_int nrhs, bool tr) const override;

    // Get name of the plugin
    const char* plugin_name() const override { return "qr";}

    // Get name of the class
    std::string class_name() const override { return "LinsolQr";}

    /// A documentation string
    static const std::string meta_doc;

    /// Symbolic factorization
    std::vector<casadi_int> prinv_, pc_;
    Sparsity sp_v_, sp_r_;
    double eps_;

    /// Cache size
    casadi_int n_cache_;
    casadi_int cache_stride_;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new LinsolQr(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit LinsolQr(DeserializingStream& s);
  };

} // namespace casadi

/// \endcond

#endif // CASADI_LINSOL_QR_HPP
