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


#ifndef CASADI_MADNLP_INTERFACE_HPP
#define CASADI_MADNLP_INTERFACE_HPP

//#include <casadi/interfaces/madnlp/casadi_nlpsol_madnlp_export.h>
#include <iostream>
#include "casadi/core/nlpsol_impl.hpp"
#include "casadi/core/timing.hpp"
#include "MadnlpCInterface.h"

namespace casadi {
  #include "madnlp_runtime.hpp"
}

/** */

#define CASADI_NLPSOL_MADNLP_EXPORT __attribute__((visibility("default")))
/** \pluginsection{Nlpsol,madnlp} **/

/// \cond INTERNAL
namespace casadi {

struct CASADI_NLPSOL_MADNLP_EXPORT MadnlpMemory : public NlpsolMemory {
  // Problem data structure
  casadi_madnlp_data<double> d;
};

/** \brief \pluginbrief{Nlpsol,madnlp}

    @copydoc Nlpsol_doc
    @copydoc plugin_Nlpsol_madnlp
*/
class CASADI_NLPSOL_MADNLP_EXPORT MadnlpInterface : public Nlpsol {
 public:
  Sparsity gradf_sp_;
  Sparsity jacg_sp_;
  Sparsity hesslag_sp_;

  explicit MadnlpInterface(const std::string& name, const Function& nlp);
  ~MadnlpInterface() override;

  // Get name of the plugin
  const char* plugin_name() const override { return "madnlp";}

  // Get name of the class
  std::string class_name() const override { return "MadnlpInterface";}

  /** \brief  Create a new NLP Solver */
  static Nlpsol* creator(const std::string& name, const Function& nlp) {
    return new MadnlpInterface(name, nlp);
  }

  ///@{
  /** \brief Options */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  // Initialize the solver
  void init(const Dict& opts) override;

  /** \brief Create memory block */
  void* alloc_mem() const override { return new MadnlpMemory();}

  /** \brief Initalize memory block */
  int init_mem(void* mem) const override;

  /** \brief Free memory block */
  void free_mem(void *mem) const override;

  /// Get all statistics
  Dict get_stats(void* mem) const override;

  /** \brief Set the (persistent) work vectors */
  void set_work(void* mem, const double**& arg, double**& res,
                        casadi_int*& iw, double*& w) const override;

  // Solve the NLP
  int solve(void* mem) const override;

  /// Exact Hessian?
  bool exact_hessian_;

  /// All MADNLP options
  Dict opts_;

  /// A documentation string
  static const std::string meta_doc;

  // Options

  /// Data for convexification
  ConvexifyData convexify_data_;

  /// convexify?
  bool convexify_;

  void set_madnlp_prob();
  void set_madnlp_prob(CodeGenerator& g) const;

  /** \brief Generate code for the function body */
  void codegen_body(CodeGenerator& g) const override;

  /** \brief Generate code for the declarations of the C function */
  void codegen_declarations(CodeGenerator& g) const override;

  /** \brief Codegen alloc_mem */
  void codegen_init_mem(CodeGenerator& g) const override;

  /** \brief Codegen free_mem */
  void codegen_free_mem(CodeGenerator& g) const override;

  /** \brief Thread-local memory object type */
  std::string codegen_mem_type() const override { return "struct casadi_madnlp_data"; }

  /** \brief Serialize an object without type information */
  void serialize_body(SerializingStream &s) const override;

  /** \brief Deserialize into MX */
  static ProtoFunction* deserialize(DeserializingStream& s) { return new MadnlpInterface(s); }

 protected:
  /** \brief Deserializing constructor */
  explicit MadnlpInterface(DeserializingStream& s);

 private:
  // Memory structure
  casadi_madnlp_prob<double> p_;

  std::vector<casadi_int> nws_;
  std::vector<casadi_int> ngs_;
};

} // namespace casadi
/// \endcond

#endif // CASADI_MADNLP_INTERFACE_HPP
