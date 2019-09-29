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


#ifndef CASADI_OSQP_INTERFACE_HPP
#define CASADI_OSQP_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/osqp/casadi_conic_osqp_export.h>

// OSQP header
extern "C" {
#include "osqp.h" // NOLINT(build/include)
}

/** \defgroup plugin_Conic_osqp
    Interface to the OSQP Solver for quadratic programming
*/

/** \pluginsection{Conic,osqp} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_CONIC_OSQP_EXPORT OsqpMemory : public ConicMemory {
    // Structures
    OSQPWorkspace* work;

    /// Constructor
    OsqpMemory();

    /// Destructor
    ~OsqpMemory();
  };

  /** \brief \pluginbrief{Conic,osqp}

      @copydoc Conic_doc
      @copydoc plugin_Conic_osqp

  */
  class CASADI_CONIC_OSQP_EXPORT OsqpInterface : public Conic {
  public:
    /** \brief  Create a new Solver */
    explicit OsqpInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st);

    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                                     const std::map<std::string, Sparsity>& st) {
      return new OsqpInterface(name, st);
    }

    /** \brief  Destructor */
    ~OsqpInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "osqp";}

    // Get name of the class
    std::string class_name() const override { return "OsqpInterface";}

    ///@{
    /** \brief const Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new OsqpMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<OsqpMemory*>(mem);}

    /// Solve the QP
    int solve(const double** arg, double** res,
        casadi_int* iw, double* w, void* mem) const override;

    /// Can discrete variables be treated
    bool integer_support() const override { return true;}

    /// Can psd constraints be treated
    bool psd_support() const override { return true;}

    /// A documentation string
    static const std::string meta_doc;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Number of nonzeros in upper part of Hessian
    casadi_int nnzHupp_, nnzA_;

    OSQPSettings settings_;

    bool warm_start_primal_, warm_start_dual_;

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Codegen alloc_mem */
    void codegen_init_mem(CodeGenerator& g) const override;

    /** \brief Codegen free_mem */
    void codegen_free_mem(CodeGenerator& g) const override;

    /** \brief Thread-local memory object type */
    std::string codegen_mem_type() const override { return "OSQPWorkspace"; }

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new OsqpInterface(s); }

  protected:
     /** \brief Deserializing constructor */
    explicit OsqpInterface(DeserializingStream& e);
  };

} // namespace casadi

/// \endcond
#endif // CASADI_OSQP_INTERFACE_HPP
