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


#ifndef CASADI_BLAZING_SPLINE_IMPL_HPP
#define CASADI_BLAZING_SPLINE_IMPL_HPP

#include "blazing_spline.hpp"
#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {
  class CASADI_EXPORT BlazingSplineFunction : public FunctionInternal {
  public:
    /** \brief Constructor */
    BlazingSplineFunction(
        const std::string& name,
        const std::vector< std::vector<double> >& knots,
        casadi_int diff_order);

    /** \brief Get type name */
    std::string class_name() const override { return "BlazingSplineFunction";}

    /** \brief Destructor */
    ~BlazingSplineFunction() override;

    /** \brief List merge opportunitities */
    void merge(const std::vector<MX>& arg,
        std::vector<MX>& subs_from, std::vector<MX>& subs_to) const override;

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Initialize */
    void init(const Dict& opts) override;

    void init_derived_members();

    /** \brief Is codegen supported?*/
    bool has_codegen() const override { return true;}

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    ///@{
    /** \brief Jacobian of all outputs with respect to all inputs */
    bool has_jacobian() const override;
    Function get_jacobian(const std::string& name,
                          const std::vector<std::string>& inames,
                          const std::vector<std::string>& onames,
                          const Dict& opts) const override;
    ///@}

    casadi_int diff_order_;
    std::vector< std::vector<double> > knots_;

    // Derived fiels
    std::vector<casadi_int> knots_offset_;
    std::vector<double> knots_stacked_;

    // Coefficient tensor size
    casadi_int nc_, ndc_, nddc_;

    ///@{
    /** \brief Number of function inputs and outputs*/
    size_t get_n_in() override;
    size_t get_n_out() override;
    ///@}

    /** \brief Which inputs are differentiable? */
    bool get_diff_in(casadi_int i) override;

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

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s);

    /** \brief String used to identify the immediate FunctionInternal subclass */
    std::string serialize_base_function() const override { return "BlazingSplineFunction"; }

    protected:
        /** \brief Deserializing constructor*/
        explicit BlazingSplineFunction(DeserializingStream& s);
  };


} // namespace casadi
/// \endcond

#endif // CASADI_BLAZING_SPLINE_IMPL_HPP
