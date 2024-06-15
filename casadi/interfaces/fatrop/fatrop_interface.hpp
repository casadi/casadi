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


#ifndef CASADI_FATROP_INTERFACE_HPP
#define CASADI_FATROP_INTERFACE_HPP

#include <casadi/interfaces/fatrop/casadi_nlpsol_fatrop_export.h>
#include "casadi/core/nlpsol_impl.hpp"
#include "casadi/core/timing.hpp"
#include <ocp/OCPAbstract.hpp>
#include <ocp/StageOCPApplication.hpp>
#include <ocp/OCPCInterface.h>

namespace casadi {
  #include "fatrop_runtime.hpp"
}

/** */

/** \pluginsection{Nlpsol,fatrop} **/

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_FATROP_EXPORT FatropMemory : public NlpsolMemory {
    // Problem data structure
    casadi_fatrop_data<double> d;


  };

  /** \brief \pluginbrief{Nlpsol,fatrop}

      @copydoc Nlpsol_doc
      @copydoc plugin_Nlpsol_fatrop
  */
  class CASADI_NLPSOL_FATROP_EXPORT FatropInterface : public Nlpsol {
  public:
    Sparsity jacg_sp_;
    Sparsity hesslag_sp_;

    explicit FatropInterface(const std::string& name, const Function& nlp);
    ~FatropInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "fatrop";}

    // Get name of the class
    std::string class_name() const override { return "FatropInterface";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new FatropInterface(name, nlp);
    }

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new FatropMemory();}

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

    /// All FATROP options
    Dict opts_;

    /// A documentation string
    static const std::string meta_doc;

    // Options

    /// Data for convexification
    ConvexifyData convexify_data_;

    /// convexify?
    bool convexify_;

    void set_fatrop_prob();
    void set_fatrop_prob(CodeGenerator& g) const;

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Codegen alloc_mem */
    void codegen_init_mem(CodeGenerator& g) const override;

    /** \brief Codegen free_mem */
    void codegen_free_mem(CodeGenerator& g) const override;

    /** \brief Thread-local memory object type */
    std::string codegen_mem_type() const override { return "struct casadi_fatrop_data"; }

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new FatropInterface(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit FatropInterface(DeserializingStream& s);

  private:
    // Memory structure
    casadi_fatrop_prob<double> p_;

    static Sparsity blocksparsity(casadi_int rows, casadi_int cols,
                                   const std::vector<casadi_ocp_block>& blocks, bool eye=false);
    static void blockptr(std::vector<double *>& vs, std::vector<double>& v,
      const std::vector<casadi_ocp_block>& blocks, bool eye=false);
    Sparsity Isp_, ABsp_, CDsp_, RSQsp_;

    std::vector< casadi_ocp_block > AB_blocks_, CD_blocks_, RSQ_blocks_, I_blocks_;

    std::vector<casadi_int> nxs_;
    std::vector<casadi_int> nus_;
    std::vector<casadi_int> ngs_;
    casadi_int N_;

    // An enum field for the structure detection
    enum StructureDetection {
      STRUCTURE_NONE,
      STRUCTURE_AUTO,
      STRUCTURE_MANUAL
    };
    StructureDetection structure_detection_;

    std::vector<casadi_int> AB_offsets_, CD_offsets_, RSQ_offsets_, I_offsets_;
    bool debug_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_FATROP_INTERFACE_HPP
