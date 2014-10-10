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


      #include "smith_lr_dle_internal.hpp"
      #include <string>

      const std::string casadi::SmithLrDleInternal::meta_doc=
      "\n"
"Solving the Low-rank Discrete Lyapunov Equations with Smith iterations\n"
"\n"
"DleSolversmith  LrDleSolversmith\n"
"\n"
"Implementation details:\n"
"We avoid ever holding P in memory as it might be large norm_inf_mul_tt\n"
"was used to obtain a stopping criteria\n"
"\n"
"We avoid memory allocation in evaluate. All sparsity pattern\n"
"calculations have been done at init\n"
"\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| max_iter        | OT_INTEGER      | 100             | Maximum number  |\n"
"|                 |                 |                 | of iterations   |\n"
"|                 |                 |                 | for the         |\n"
"|                 |                 |                 | algorithm       |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| print_iteration | OT_BOOLEAN      | false           | Print           |\n"
"|                 |                 |                 | information     |\n"
"|                 |                 |                 | about each      |\n"
"|                 |                 |                 | iteration       |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| tol             | OT_REAL         | 0.000           | Tolerance for   |\n"
"|                 |                 |                 | satisfying the  |\n"
"|                 |                 |                 | Lyapunov        |\n"
"|                 |                 |                 | equation.       |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
">List of available stats\n"
"\n"
"+------------+\n"
"|     Id     |\n"
"+============+\n"
"| iter_count |\n"
"+------------+\n"
"\n"
"\n"
"\n"
"\n"
;
