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


#ifndef CASADI_SUPERSCS_INTERFACE_HPP
#define CASADI_SUPERSCS_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/superscs/casadi_conic_superscs_export.h>

#include <scs_parser.h>

/** \defgroup plugin_Conic_superscs
    Interface to the SuperSCS solver for conic programming


    Joris Gillis, 2019
*/

/** \pluginsection{Conic,superscs} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_CONIC_SUPERSCS_EXPORT SuperscsMemory : public ConicMemory {

    ScsSolution* sol;
    ScsData data;;

    ScsInfo * info;
    ScsCone cone;
    ScsAMatrix A;
    ScsSettings settings;

    int return_status;

    std::vector<casadi_int> at_colind, at_row, q;
    std::vector<double> ldl_d, ldl_l, ldl_w, F_res, g;

    /// Constructor
    SuperscsMemory();

    /// Destructor
    ~SuperscsMemory();
  };

  /** \brief \pluginbrief{Conic,superscs}

      @copydoc Conic_doc
      @copydoc plugin_Conic_superscs

  */
  class CASADI_CONIC_SUPERSCS_EXPORT SuperscsInterface : public Conic {
  public:
    /** \brief  Create a new Solver */
    explicit SuperscsInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st);

    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                                     const std::map<std::string, Sparsity>& st) {
      return new SuperscsInterface(name, st);
    }

    /** \brief  Destructor */
    ~SuperscsInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "superscs";}

    // Get name of the class
    std::string class_name() const override { return "SuperscsInterface";}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new SuperscsMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<SuperscsMemory*>(mem);}

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

    /// Superscs options
    Dict opts_;

    /// SDP to SOCP conversion memory
    SDPToSOCPMem sdp_to_socp_mem_;

    IM At_;

    std::vector<casadi_int> lookup_;

    // Symbolic of H
    std::vector<casadi_int> Hp_;
    Sparsity HL_sp_;

    Function F_;

    ScsSettings settings_;

    // Socp perturbation
    std::vector<casadi_int> perturb_;

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new SuperscsInterface(s); }

  protected:
     /** \brief Deserializing constructor */
    explicit SuperscsInterface(DeserializingStream& s);
  };

} // namespace casadi

/// \endcond
#endif // CASADI_SUPERSCS_INTERFACE_HPP
