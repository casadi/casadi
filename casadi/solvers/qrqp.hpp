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


#ifndef CASADI_QRQP_HPP
#define CASADI_QRQP_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/solvers/casadi_conic_qrqp_export.h>

/** \defgroup plugin_Conic_qrqp
 Solve QPs using an active-set method
*/

/** \pluginsection{Conic,qrqp} */

/// \cond INTERNAL
namespace casadi {
  struct CASADI_CONIC_QRQP_EXPORT QrqpMemory : public ConicMemory {
    const char* return_status;
  };

  /** \brief \pluginbrief{Conic,qrqp}

      @copydoc Conic_doc
      @copydoc plugin_Conic_qrqp

      \author Joel Andersson
      \date 2018
  */
  class CASADI_CONIC_QRQP_EXPORT Qrqp : public Conic {
  public:
    /** \brief  Create a new Solver */
    explicit Qrqp(const std::string& name,
                  const std::map<std::string, Sparsity> &st);

    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new Qrqp(name, st);
    }

    /** \brief  Destructor */
    ~Qrqp() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "qrqp";}

    // Get name of the class
    std::string class_name() const override { return "Qrqp";}

    /** \brief Create memory block */
    void* alloc_mem() const override { return new QrqpMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<QrqpMemory*>(mem);}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Initialize */
    void init(const Dict& opts) override;

    /** \brief Solve the QP */
    int solve(const double** arg, double** res,
             casadi_int* iw, double* w, void* mem) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /// A documentation string
    static const std::string meta_doc;
    // Memory structure
    casadi_qp_prob<double> p_;
    // KKT system and its QR factorization
    Sparsity AT_, kkt_, sp_v_, sp_r_;
    // KKT system permutation
    std::vector<casadi_int> prinv_, pc_;
    ///@{
    // Options
    bool print_iter_, print_header_, print_info_, print_lincomb_;
    ///@}

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new Qrqp(s); }

  protected:
     /** \brief Deserializing constructor */
    explicit Qrqp(DeserializingStream& s);

  private:
    void set_qp_prob();
  };

} // namespace casadi
/// \endcond
#endif // CASADI_QRQP_HPP
