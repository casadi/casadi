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


#ifndef CASADI_CONIC_IMPL_HPP
#define CASADI_CONIC_IMPL_HPP

#include "conic.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"
#include "im.hpp"

/// \cond INTERNAL
namespace casadi {

  struct CASADI_EXPORT ConicMemory : public FunctionMemory {
    // Problem data structure
    casadi_qp_data<double> d_qp;

  };

  /// Internal class
  class CASADI_EXPORT Conic : public FunctionInternal, public PluginInterface<Conic> {
  public:
    // Memory structure
    casadi_qp_prob<double> p_qp_;

    // Constructor
    Conic(const std::string& name, const std::map<std::string, Sparsity> &st);

    // Destructor
    ~Conic() override = 0;

    ///@{
    /** \brief Number of function inputs and outputs

        \identifier{1j6} */
    size_t get_n_in() override { return CONIC_NUM_IN;}
    size_t get_n_out() override { return CONIC_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs

        \identifier{1j7} */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs

        \identifier{1j8} */
    std::string get_name_in(casadi_int i) override { return conic_in(i);}
    std::string get_name_out(casadi_int i) override { return conic_out(i);}
    /// @}

    ///@{
    /** \brief Options

        \identifier{1j9} */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Solve the QP
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const final;

    /// Solve the QP
    virtual int solve(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const = 0;

    // Initialize
    void init(const Dict& opts) override;

    /** \brief Initalize memory block

        \identifier{1ja} */
    int init_mem(void* mem) const override;

    /** \brief Set the (persistent) work vectors

        \identifier{1jb} */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    /// \brief Check if the numerical values of the supplied bounds make sense
    virtual void check_inputs(const double* lbx, const double* ubx,
                             const double* lba, const double* uba) const;

    /** Generate native code in the interfaced language for debugging */
    virtual void generateNativeCode(std::ostream& file) const;

    // Creator function for internal class
    typedef Conic* (*Creator)(const std::string& name,
                              const std::map<std::string, Sparsity>& st);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "conic";}

    /** \brief Check if the function is of a particular type

        \identifier{1jc} */
    bool is_a(const std::string& type, bool recursive) const override;

    /** \brief Get default input value

        \identifier{1jd} */
    double get_default_in(casadi_int ind) const override;

    /// Can discrete variables be treated
    virtual bool integer_support() const { return false;}

    /// Can psd constraints be treated
    virtual bool psd_support() const { return false;}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief Generate code for the function body

        \identifier{24y} */
    void qp_codegen_body(CodeGenerator& g) const;

  protected:
    /// Options
    std::vector<bool> discrete_;
    bool print_problem_;

    /// Problem structure
    Sparsity H_, A_, Q_, P_;

    /// Number of decision variables
    casadi_int nx_;

    /// The number of constraints (counting both equality and inequality) == A.size1()
    casadi_int na_;

    /// The shape of psd constraint matrix
    casadi_int np_;

    /// SDP to SOCP conversion memory
    struct SDPToSOCPMem {
      // Block partition vector for SOCP (block i runs from r[i] to r[i+1])
      std::vector<casadi_int> r;

      // Tranpose of A, and corresponding mapping
      Sparsity AT;
      std::vector<casadi_int> A_mapping;

      // Aggregate SOCP helper constraints (lhs)
      IM map_Q;

      // Aggregate SOCP helper constraints (rhs)
      std::vector<casadi_int> map_P;

      // Maximum size of ind/val vectors
      casadi_int indval_size;
    };

    /// SDP to SOCP conversion initialization
    void sdp_to_socp_init(SDPToSOCPMem& mem) const;

    void serialize(SerializingStream &s, const SDPToSOCPMem& m) const;
    void deserialize(DeserializingStream &s, SDPToSOCPMem& m);

  public:
      /** \brief Serialize an object without type information

          \identifier{1je} */
    void serialize_body(SerializingStream &s) const override;
    /** \brief Serialize type information

        \identifier{1jf} */
    void serialize_type(SerializingStream &s) const override;

    /** \brief String used to identify the immediate FunctionInternal subclass

        \identifier{1jg} */
    std::string serialize_base_function() const override { return "Conic"; }
    /** \brief Deserialize with type disambiguation

        \identifier{1jh} */
    static ProtoFunction* deserialize(DeserializingStream& s);

  protected:

    /** \brief Deserializing constructor

        \identifier{1ji} */
    explicit Conic(DeserializingStream& s);
  private:
    void set_qp_prob();
  };


} // namespace casadi
/// \endcond
#endif // CASADI_CONIC_IMPL_HPP
