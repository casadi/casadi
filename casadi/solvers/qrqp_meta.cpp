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


      #include "qrqp.hpp"
      #include <string>

      const std::string casadi::Qrqp::meta_doc=
      "\n"
"\n"
"\n"
"Solve QPs using an active-set method\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_22y \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------+--------------------------------------------+\n"
"|       Id        |   Type    |                Description                 |\n"
"+=================+===========+============================================+\n"
"| constr_viol_tol | OT_DOUBLE | Constraint violation tolerance [1e-8].     |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| dual_inf_tol    | OT_DOUBLE | Dual feasibility violation tolerance       |\n"
"|                 |           | [1e-8]                                     |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| max_iter        | OT_INT    | Maximum number of iterations [1000].       |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| min_lam         | OT_DOUBLE | Smallest multiplier treated as inactive    |\n"
"|                 |           | for the initial active set [0].            |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| print_header    | OT_BOOL   | Print header [true].                       |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| print_info      | OT_BOOL   | Print info [true].                         |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| print_iter      | OT_BOOL   | Print iterations [true].                   |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| print_lincomb   | OT_BOOL   | Print dependant linear combinations of     |\n"
"|                 |           | constraints [false]. Printed numbers are   |\n"
"|                 |           | 0-based indices into the vector of [simple |\n"
"|                 |           | bounds;linear bounds]                      |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
