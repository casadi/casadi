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


      #include "qrsqp.hpp"
      #include <string>

      const std::string casadi::Qrsqp::meta_doc=
      "\n"
"\n"
"\n"
"A textbook SQPMethod\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_22u \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------------+-----------+--------------------------------------+\n"
"|          Id           |   Type    |             Description              |\n"
"+=======================+===========+======================================+\n"
"| beta                  | OT_DOUBLE | Line-search parameter, restoration   |\n"
"|                       |           | factor of stepsize                   |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| c1                    | OT_DOUBLE | Armijo condition, coefficient of     |\n"
"|                       |           | decrease in merit                    |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| hessian_approximation | OT_STRING | limited-memory|exact                 |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| lbfgs_memory          | OT_INT    | Size of L-BFGS memory.               |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| max_iter              | OT_INT    | Maximum number of SQP iterations     |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| max_iter_ls           | OT_INT    | Maximum number of linesearch         |\n"
"|                       |           | iterations                           |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| merit_memory          | OT_INT    | Size of memory to store history of   |\n"
"|                       |           | merit function values                |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| min_iter              | OT_INT    | Minimum number of SQP iterations     |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| min_step_size         | OT_DOUBLE | The size (inf-norm) of the step size |\n"
"|                       |           | should not become smaller than this. |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| print_header          | OT_BOOL   | Print the header with problem        |\n"
"|                       |           | statistics                           |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| print_iteration       | OT_BOOL   | Print the iterations                 |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| qpsol                 | OT_STRING | The QP solver to be used by the SQP  |\n"
"|                       |           | method [qrqp]                        |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| qpsol_options         | OT_DICT   | Options to be passed to the QP       |\n"
"|                       |           | solver                               |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| regularize            | OT_BOOL   | Automatic regularization of Lagrange |\n"
"|                       |           | Hessian.                             |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| tol_du                | OT_DOUBLE | Stopping criterion for dual          |\n"
"|                       |           | infeasability                        |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| tol_pr                | OT_DOUBLE | Stopping criterion for primal        |\n"
"|                       |           | infeasibility                        |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
