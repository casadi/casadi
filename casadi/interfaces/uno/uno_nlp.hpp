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

#ifndef CASADI_UNO_NLP_HPP
#define CASADI_UNO_NLP_HPP

#include "Uno_C_API.h"

namespace casadi {

// Forward declarations
struct UnoMemory;

class UnoNlp{

public:
   UnoNlp(UnoMemory* mem);
   ~UnoNlp() {}

   // objective
   uno_int objective_function(const double* x, double* objective_value);
   static uno_int objective_function_wrapper(uno_int /*number_variables*/, const double* x, double* objective_value, void* user_data);
   
   uno_int objective_gradient(const double* x, double* gradient);
   static uno_int objective_gradient_wrapper(uno_int /*number_variables*/, const double* x, double* gradient, void* /*user_data*/);

   // constraints
   uno_int constraint_functions(const double* x, double* constraint_values);
   static uno_int constraint_functions_wrapper(uno_int /*number_variables*/, uno_int /*number_constraints*/, const double* x, double* constraint_values, void* user_data);

   uno_int jacobian(const double* x, double* jacobian_values);
   static uno_int jacobian_wrapper(uno_int /*number_variables*/, uno_int /*number_jacobian_nonzeros*/, const double* x, double* jacobian_values, void* user_data);
   // Hessian

//    void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
//          SymmetricMatrix<double>& hessian) const override;

//    void get_initial_primal_point(std::vector<double>& x) const override;
//    void get_initial_dual_point(std::vector<double>& multipliers) const override;
//    void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override;

//    [[nodiscard]] const std::vector<size_t>& get_linear_constraints() const override;

private:
   UnoMemory* mem_;

//    mutable std::vector<double> casadi_tmp_gradient{};
//    mutable std::vector<double> casadi_tmp_multipliers{};
//    mutable std::vector<double> casadi_tmp_constraint_jacobian{};
//    mutable std::vector<double> casadi_tmp_hessian{};

//    std::vector<Interval> variables_bounds;
//    std::vector<Interval> constraint_bounds;
//    std::vector<BoundType> variable_status; /*!< Status of the variables (EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES) */
//    std::vector<FunctionType> constraint_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */
//    std::vector<BoundType> constraint_status; /*!< Status of the constraints (EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES,
//  * UNBOUNDED) */

//    std::vector<size_t> linear_constraints;

//    void generate_variables();
//    void generate_constraints();
//   //  void set_function_types(std::string file_name);

//    void set_number_hessian_nonzeros();
//    [[nodiscard]] size_t compute_hessian_number_nonzeros(double objective_multiplier, const std::vector<double>& multipliers) const;
};

} //namespace casadi

#endif // CASADI_UNO_NLP_HPP
