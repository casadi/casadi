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


      #include "simple_indef_dle_internal.hpp"
      #include <string>

      const std::string casadi::SimpleIndefDleInternal::meta_doc=
      "\n"
"Solving the Discrete Lyapunov Equations with a regular LinearSolver\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| compressed_solv | OT_BOOLEAN      | true            | When a system   |\n"
"| e               |                 |                 | with sparse rhs |\n"
"|                 |                 |                 | arises,         |\n"
"|                 |                 |                 | compress toa    |\n"
"|                 |                 |                 | smaller system  |\n"
"|                 |                 |                 | with dense rhs. |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| linear_solver   | OT_STRING       | GenericType()   | User-defined    |\n"
"|                 |                 |                 | linear solver   |\n"
"|                 |                 |                 | class. Needed   |\n"
"|                 |                 |                 | for             |\n"
"|                 |                 |                 | sensitivities.  |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| linear_solver_o | OT_DICTIONARY   | GenericType()   | Options to be   |\n"
"| ptions          |                 |                 | passed to the   |\n"
"|                 |                 |                 | linear solver.  |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
