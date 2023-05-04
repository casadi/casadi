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


#ifndef CASADI_LSQR_HPP
#define CASADI_LSQR_HPP

#include "casadi/core/linsol_internal.hpp"
#include <casadi/solvers/casadi_linsol_lsqr_export.h>

/** \defgroup plugin_Linsol_symbolicqr Title
    \par

    Linear solver for sparse least-squares problems
    Inspired from https://github.com/scipy/scipy/blob/v0.14.0/scipy/sparse/linalg/isolve/lsqr.py#L96

    \identifier{230} */

/** \pluginsection{Linsol,symbolicqr} */

/// \cond INTERNAL

namespace casadi {

  /** \brief Memory for SymbolicQR  */
  struct CASADI_LINSOL_LSQR_EXPORT LsqrMemory : public LinsolMemory {
    // Work vectors
    std::vector<const double*> arg;
    std::vector<double> w;

    std::vector<double> A;
  };

  /** \brief \pluginbrief{Linsol,symbolicqr}

      @copydoc Linsol_doc
      @copydoc plugin_Linsol_symbolicqr
      \author Joel Andersson
      \date 2013
  */
  class CASADI_LINSOL_LSQR_EXPORT Lsqr : public LinsolInternal {
  public:
    // Constructor
    Lsqr(const std::string& name, const Sparsity& sp);

    // Destructor
    ~Lsqr() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "lsqr";}

    // Name of the class
    std::string class_name() const override { return "Lsqr";}

    /** \brief  Create a new Linsol */
    static LinsolInternal* creator(const std::string& name, const Sparsity& sp) {
      return new Lsqr(name, sp);
    }

    /** \brief Create memory block */
    void* alloc_mem() const override { return new LsqrMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<LsqrMemory*>(mem);}

    // Factorize the linear system
    int nfact(void* mem, const double* A) const override;

    // Solve the linear system
    int solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const override;

    /// Generate C code
    void generate(CodeGenerator& g, const std::string& A, const std::string& x,
                  casadi_int nrhs, bool tr) const override;

    /// A documentation string
    static const std::string meta_doc;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new Lsqr(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit Lsqr(DeserializingStream& s) : LinsolInternal(s) {}
  };

} // namespace casadi

/// \endcond
#endif // CASADI_LSQR_HPP
