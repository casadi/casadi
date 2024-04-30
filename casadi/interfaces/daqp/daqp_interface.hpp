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

#ifndef CASADI_DAQP_INTERFACE_HPP
#define CASADI_DAQP_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/daqp/casadi_conic_daqp_export.h>

#include <api.h>
#include <string>

namespace casadi {
  #include "daqp_runtime.hpp"
}

/** \defgroup plugin_Conic_daqp Title
    \par

      Interface to Daqp solver for sparse Quadratic Programs,
      see daqp.dev for more information and
      https://www.maths.ed.ac.uk/hall/Daqp/DaqpOptions.html
      for a list of options.

    \identifier{27p} */

/** \pluginsection{Conic,daqp} */

/// \cond INTERNAL

namespace casadi {

  struct CASADI_CONIC_DAQP_EXPORT DaqpMemory : public ConicMemory {
    // Problem data structure
    casadi_daqp_data<double> d;

  };

  /** \brief \pluginbrief{Conic,daqp}

      @copydoc Daqp_doc
      @copydoc plugin_Conic_daqp


      \author Felix Lenders, Ruud Kassing, Joris Gillis
      \date 2021
  */
  class CASADI_CONIC_DAQP_EXPORT DaqpInterface : public Conic {
  public:
    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new DaqpInterface(name, st);
    }

    /// Constructor using sparsity patterns
    explicit DaqpInterface(const std::string& name,
                            const std::map<std::string, Sparsity>& st);

    /// Destructor
    ~DaqpInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "daqp";}

    // Get name of the class
    std::string class_name() const override { return "DaqpInterface";}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    void set_daqp_prob();
    void set_daqp_prob(CodeGenerator& g) const;

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Codegen decref for init_mem*/
    void codegen_init_mem(CodeGenerator& g) const override;

    /** \brief Codegen for free_mem */
    void codegen_free_mem(CodeGenerator& g) const override;

    /** \brief Thread-local memory object type */
    std::string codegen_mem_type() const override { return "struct casadi_daqp_data"; }

    // Initialize the solver
    void init(const Dict& opts) override;

    // Initialize dependant read-only members of class
    void init_dependent();

    /** \brief Create memory block */
    void* alloc_mem() const override { return new DaqpMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override;

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Solve the QP
    int solve(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    /// All Daqp options
    Dict opts_;

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new DaqpInterface(s); }

    /// Can discrete variables be treated
    bool integer_support() const override { return true; }

  protected:
     /** \brief Deserializing constructor */
    explicit DaqpInterface(DeserializingStream& s);

  private:

    // Memory structure
    casadi_daqp_prob<double> p_;

  };
} // end namespace casadi
/// \endcond
#endif // CASADI_DAQP_INTERFACE_HPP
