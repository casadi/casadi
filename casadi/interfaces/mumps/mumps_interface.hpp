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


#ifndef CASADI_MUMPS_INTERFACE_HPP
#define CASADI_MUMPS_INTERFACE_HPP

#include "casadi/core/linsol_internal.hpp"
#include <casadi/interfaces/mumps/casadi_linsol_mumps_export.h>

#include <mumps_seq/mpi.h>
#include <dmumps_c.h>

/** \defgroup plugin_Linsol_mumps
 * Interface to the sparse direct linear solver MUMPS
 * Works for symmetric indefinite systems
 * \author Joel Andersson
 * \date 2019
 */

/** \pluginsection{Linsol,mumps} */
/// \cond INTERNAL
namespace casadi {
  struct CASADI_LINSOL_MUMPS_EXPORT MumpsMemory : public LinsolMemory {
    // Constructor
    MumpsMemory();

    // Destructor
    ~MumpsMemory();

    // Memory block
    DMUMPS_STRUC_C* id;

    // Sparsity
    std::vector<int> irn, jcn;

    // Nonzeros
    std::vector<double> nz;
  };

  /** \brief \pluginbrief{Linsol,mumps}
   * @copydoc Linsol_doc
   * @copydoc plugin_Linsol_mumps
   */
  class CASADI_LINSOL_MUMPS_EXPORT MumpsInterface : public LinsolInternal {
  public:

    // Create a linear solver given a sparsity pattern and a number of right hand sides
    MumpsInterface(const std::string& name, const Sparsity& sp);

    /** \brief  Create a new Linsol */
    static LinsolInternal* creator(const std::string& name, const Sparsity& sp) {
      return new MumpsInterface(name, sp);
    }

    // Destructor
    ~MumpsInterface() override;

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new MumpsMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<MumpsMemory*>(mem);}

    // Factorize the linear system
    int nfact(void* mem, const double* A) const override;

    // Solve the linear system
    int solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const override;

    /// A documentation string
    static const std::string meta_doc;

    // Get name of the plugin
    const char* plugin_name() const override { return "mumps";}

    // Get name of the class
    std::string class_name() const override { return "MumpsInterface";}

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new MumpsInterface(s); }

    ///@{
    // Options
    bool symmetric_, posdef_;
    ///@}

  protected:
    /** \brief Deserializing constructor */
    explicit MumpsInterface(DeserializingStream& s) : LinsolInternal(s) {}
  };

} // namespace casadi

/// \endcond

#endif // CASADI_MUMPS_INTERFACE_HPP
