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


#ifndef CASADI_CORE_HPP
#define CASADI_CORE_HPP

// Scalar expressions (why do I need to put it up here?)
#include "sx/sx_element.hpp"

// Generic tools
#include "polynomial.hpp"
#include "matrix/generic_matrix_tools.hpp"
#include "matrix/generic_expression_tools.hpp"
#include "std_vector_tools.hpp"
#include "functor.hpp"
#include "casadi_options.hpp"
#include "casadi_meta.hpp"

// Matrices
#include "matrix/matrix.hpp"
#include "matrix/matrix_tools.hpp"
#include "matrix/sparsity_tools.hpp"

// Scalar expressions
#include "sx/sx_tools.hpp"

// Matrix expressions
#include "mx/mx.hpp"
#include "mx/mx_tools.hpp"

// Functions
#include "function/sx_function.hpp"
#include "function/mx_function.hpp"
#include "function/external_function.hpp"
#include "function/linear_solver.hpp"
#include "function/nlp_solver.hpp"
#include "function/integrator.hpp"
#include "function/implicit_function.hpp"
#include "function/custom_function.hpp"
#include "function/simulator.hpp"
#include "function/parallelizer.hpp"
#include "function/control_simulator.hpp"
#include "function/qp_solver.hpp"
#include "function/homotopy_nlp_solver.hpp"
#include "function/stabilized_qp_solver.hpp"
#include "function/lp_solver.hpp"
#include "function/sdp_solver.hpp"
#include "function/socp_solver.hpp"
#include "function/qcqp_solver.hpp"
#include "function/sdqp_solver.hpp"
#include "function/nullspace.hpp"
#include "function/lr_dle_solver.hpp"
#include "function/lr_dple_solver.hpp"
#include "function/dple_solver.hpp"
#include "function/cle_solver.hpp"
#include "function/dle_solver.hpp"

// Misc
#include "misc/integration_tools.hpp"
#include "misc/symbolic_nlp.hpp"
#include "misc/variable.hpp"
#include "misc/symbolic_ocp.hpp"
#include "misc/xml_file.hpp"

#endif // CASADI_CORE_HPP
