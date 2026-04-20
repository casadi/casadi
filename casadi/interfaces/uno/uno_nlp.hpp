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
#include "casadi/core/generic_type.hpp"
#include <string>

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

   uno_int lagrangian_hessian(const double* x, double objective_multiplier, const double* multipliers, double* hessian_values);
   static uno_int lagrangian_hessian_wrapper(uno_int /*number_variables*/, uno_int /*number_constraints*/, uno_int /*number_hessian_nonzeros*/,
            const double* x, double objective_multiplier, const double* multipliers, double* hessian_values, void* user_data);

   static void set_uno_option(void* solver, const std::string& name, const GenericType& value);
   static void insert_casadi_options(void* solver, Dict opts);

private:
   UnoMemory* mem_;
};

} //namespace casadi

#endif // CASADI_UNO_NLP_HPP
