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


#ifndef CASADI_INTERPOLANT_IMPL_HPP
#define CASADI_INTERPOLANT_IMPL_HPP

#include "interpolant.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Internal class
      @copydoc Interpolant_doc
  */
  class CASADI_EXPORT Interpolant
  : public FunctionInternal, public PluginInterface<Interpolant> {
  public:
    /// Constructor
    Interpolant(const std::string& name,
                const std::vector<double>& grid,
                const std::vector<casadi_int>& offset,
                const std::vector<double>& values,
                casadi_int m);

    /// Destructor
    ~Interpolant() override;

    /** \brief Get type name */
    std::string class_name() const override {return "Interpolant";}

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override { return 1+has_parametric_values()+has_parametric_grid();}
    size_t get_n_out() override { return 1;}
    ///@}

    bool is_diff_in(casadi_int i) override { return i==0; }

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(casadi_int i) override;
    std::string get_name_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Initialize
    void init(const Dict& opts) override;

    /// Convert from (optional) lookup modes labels to enum
    static std::vector<casadi_int> interpret_lookup_mode(const std::vector<std::string>& modes,
        const std::vector<double>& grid, const std::vector<casadi_int>& offset,
        const std::vector<casadi_int>& margin_left=std::vector<casadi_int>(),
        const std::vector<casadi_int>& margin_right=std::vector<casadi_int>());

    static void stack_grid(const std::vector< std::vector<double> >& grid,
      std::vector<casadi_int>& offset, std::vector<double>& stacked);

    static void check_grid(const std::vector< std::vector<double> >& grid);
    static void check_grid(const std::vector<casadi_int>& grid);

    static std::vector<double> meshgrid(const std::vector< std::vector<double> >& grid);

    // Creator function for internal class
    typedef Interpolant* (*Creator)(const std::string& name,
                                    const std::vector<double>& grid,
                                    const std::vector<casadi_int>& offset,
                                    const std::vector<double>& values,
                                    casadi_int m);

    /** \brief  Comstruct a new Interpolant */
    static Function construct(const std::string& solver, const std::string& name,
                      const std::vector<double>& grid,
                      const std::vector<casadi_int>& offset,
                      const std::vector<double>& values,
                      casadi_int m,
                      const Dict& opts);

    typedef Function (* DoInline)(const std::string& name,
                    const std::vector<double>& grid,
                    const std::vector<casadi_int>& offset,
                    const std::vector<double>& values,
                    casadi_int m,
                    const Dict& opts);

    // No static functions exposed
    struct Exposed{ DoInline do_inline; };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    // Number of dimensions
    casadi_int ndim_;

    // Number of outputs
    casadi_int m_;

    // Number of inputs
    casadi_int batch_x_;

    // Input grid
    std::vector<double> grid_;

    // Offset for each dimension
    std::vector<casadi_int> offset_;

    // Values at gridpoints
    std::vector<double> values_;

    // Lookup modes
    std::vector<std::string> lookup_modes_;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;
    /** \brief Serialize type information */
    void serialize_type(SerializingStream &s) const override;

    /** \brief String used to identify the immediate FunctionInternal subclass */
    std::string serialize_base_function() const override { return "Interpolant"; }
    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s);

    /** \brief Is parametric? */
    bool has_parametric_values() const { return values_.empty(); }

    /** \brief Is parametric? */
    bool has_parametric_grid() const { return grid_.empty(); }

    casadi_int arg_values() const;
    casadi_int arg_grid() const;

    /** \brief Size of the flattened coefficients vector */
    casadi_int coeff_size() const;
    static casadi_int coeff_size(const std::vector<casadi_int>& offset, casadi_int m);

  protected:

    bool arg_values(casadi_int i) const;
    bool arg_grid(casadi_int i) const;

    /** \brief Deserializing constructor */
    explicit Interpolant(DeserializingStream& s);


  };

} // namespace casadi
/// \endcond

#endif // CASADI_INTERPOLANT_IMPL_HPP
