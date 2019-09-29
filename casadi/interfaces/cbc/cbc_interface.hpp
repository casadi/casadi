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

#ifndef CASADI_CBC_INTERFACE_HPP
#define CASADI_CBC_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/cbc/casadi_conic_cbc_export.h>

#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcEventHandler.hpp"
#include <string>

/** \defgroup plugin_Conic_cbc

      Interface to Cbc solver for sparse Quadratic Programs
*/

/** \pluginsection{Conic,cbc} */

/// \cond INTERNAL

namespace casadi {

  struct CASADI_CONIC_CBC_EXPORT CbcMemory : public ConicMemory {
    /// Constructor
    CbcMemory();

    /// Destructor
    ~CbcMemory();

    std::vector<int> colind, row;

    int return_status;
    int secondary_return_status;

    casadi_int iter_count;
    casadi_int node_count;

  };

  /** \brief \pluginbrief{Conic,cbc}

      @copydoc Conic_doc
      @copydoc plugin_Conic_cbc


      \author Attila Kozma, Joel Andersson
      \date 2012
  */
  class CASADI_CONIC_CBC_EXPORT CbcInterface : public Conic {
  public:
    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new CbcInterface(name, st);
    }

    /// Constructor using sparsity patterns
    explicit CbcInterface(const std::string& name,
                            const std::map<std::string, Sparsity>& st);

    /// Destructor
    ~CbcInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "cbc";}

    // Get name of the class
    std::string class_name() const override { return "CbcInterface";}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new CbcMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<CbcMemory*>(mem);}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Solve the QP
    int solve(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    /// All CBC options
    Dict opts_;

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new CbcInterface(s); }

    /// Can discrete variables be treated
    bool integer_support() const override { return true;}

  protected:
     /** \brief Deserializing constructor */
    explicit CbcInterface(DeserializingStream& s);

  private:
    // Conversion of string to enum for options
    static std::map<std::string, CbcModel::CbcIntParam> param_map_int;
    static std::map<std::string, CbcModel::CbcDblParam> param_map_double;
    static std::map<std::string, OsiIntParam> osi_param_map_int;
    static std::map<std::string, OsiDblParam> osi_param_map_double;

    void copy_cbc_results(const CbcModel& model, double** res) const;

    // SOS structure
    std::vector< std::vector<int> > sos_groups_;
    std::vector< std::vector<double> > sos_weights_;
    std::vector<casadi_int> sos_types_;

    bool hot_start_;

  };
} // end namespace casadi
/// \endcond
#endif // CASADI_CBC_INTERFACE_HPP
