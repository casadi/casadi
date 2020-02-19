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
                const MX& grid,
                const std::vector<casadi_int>& offset,
                const MX& values,
                casadi_int m);

    /// Destructor
    ~Interpolant() override;

    /** \brief Get type name */
    std::string class_name() const override {return "Interpolant";}

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

    /// Finalize initialization
    void finalize() override;

    /// Convert from (optional) lookup modes labels to enum
    static std::vector<casadi_int> interpret_lookup_mode(const std::vector<std::string>& modes,
        const MX& grid, const std::vector<casadi_int>& offset,
        const std::vector<casadi_int>& margin_left=std::vector<casadi_int>(),
        const std::vector<casadi_int>& margin_right=std::vector<casadi_int>());

    static void stack_grid(const std::vector< std::vector<double> >& grid,
      std::vector<casadi_int>& offset, std::vector<double>& stacked);

    static void stack_grid(const std::vector< MX >& grid,
      std::vector<casadi_int>& offset, MX& stacked);

    static void check_grid(const std::vector<MX>& grid);

    static std::vector< std::vector<double> > parse_grid(const std::vector< DM >& grid);
    static std::vector<MX> parse_grid(const std::vector< MX >& grid) { return grid; }
    static std::vector<double> parse_grid(const DM & grid);
    static MX parse_grid(const MX & grid) { return grid; }

    static std::vector<double> meshgrid(const std::vector< std::vector<double> >& grid);
    template<typename T>
    static T meshgrid(const std::vector< T >& grid);

    // Creator function for internal class
    typedef Interpolant* (*Creator)(const std::string& name,
                                    const MX& grid,
                                    const std::vector<casadi_int>& offset,
                                    const MX& values,
                                    casadi_int m);

    /** \brief  Comstruct a new Interpolant */
    static Function construct(const std::string& solver, const std::string& name,
                      const MX& grid,
                      const std::vector<casadi_int>& offset,
                      const MX& values,
                      casadi_int m,
                      const Dict& opts);

    typedef Function (* DoInline)(const std::string& name,
                    const MX& grid,
                    const std::vector<casadi_int>& offset,
                    const MX& values,
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
    const std::vector<double>* grid_ptr_;

    // Offset for each dimension
    std::vector<casadi_int> offset_;

    // Values at gridpoints
    const std::vector<double>* values_ptr_;

    // Lookup modes
    std::vector<std::string> lookup_modes_;

    MX grid_;
    std::vector<double> grid_vec_;
    MX values_;
    std::vector<double> values_vec_;

    /** \brief Generate code for the declarations of the C function */
    virtual void codegen_declarations(CodeGenerator& g) const;

    std::string codegen_values(CodeGenerator& g) const;
    std::string codegen_grid(CodeGenerator& g) const;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;
    /** \brief Serialize type information */
    void serialize_type(SerializingStream &s) const override;

    /** \brief String used to identify the immediate FunctionInternal subclass */
    std::string serialize_base_function() const override { return "Interpolant"; }
    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s);

    /** \brief Size of the flattened coefficients vector */
    casadi_int coeff_size() const;
    static casadi_int coeff_size(const std::vector<casadi_int>& offset, casadi_int m);

    char type() const override;
    
  protected:

    /** \brief Deserializing constructor */
    explicit Interpolant(DeserializingStream& s);


  };


  template<typename M>
  M meshgrid_fund(const M& next, const std::vector<M>& grid) {
    if (grid.empty()) return next;
    M n = horzcat(repmat(next, grid.front().size1()), vec(repmat(grid.front().T(), next.size1(), 1)));
    return meshgrid_fund(n, std::vector<M>(grid.begin()+1, grid.end()));
  }

  template<typename M>
  M Interpolant::meshgrid(const std::vector< M >& grid) {
    M ret = meshgrid_fund(grid.front(), std::vector<M>(grid.begin()+1, grid.end()));
    return vec(ret.T());
  }

} // namespace casadi
/// \endcond

#endif // CASADI_INTERPOLANT_IMPL_HPP
