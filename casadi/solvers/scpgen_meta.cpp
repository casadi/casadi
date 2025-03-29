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


      #include "scpgen.hpp"
      #include <string>

      const std::string casadi::Scpgen::meta_doc=
      "\n"
"\n"
"\n"
"A structure-exploiting sequential quadratic programming (to be come \n"
"sequential convex programming) method for nonlinear programming.\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_232 \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------------+-----------------+--------------------------------+\n"
"|          Id           |      Type       |          Description           |\n"
"+=======================+=================+================================+\n"
"| beta                  | OT_DOUBLE       | Line-search parameter,         |\n"
"|                       |                 | restoration factor of stepsize |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| c1                    | OT_DOUBLE       | Armijo condition, coefficient  |\n"
"|                       |                 | of decrease in merit           |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| codegen               | OT_BOOL         | C-code generation              |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| hessian_approximation | OT_STRING       | gauss-newton|exact             |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| lbfgs_memory          | OT_INT          | Size of L-BFGS memory.         |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| max_iter              | OT_INT          | Maximum number of SQP          |\n"
"|                       |                 | iterations                     |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| max_iter_ls           | OT_INT          | Maximum number of linesearch   |\n"
"|                       |                 | iterations                     |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| merit_memsize         | OT_INT          | Size of memory to store        |\n"
"|                       |                 | history of merit function      |\n"
"|                       |                 | values                         |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| merit_start           | OT_DOUBLE       | Lower bound for the merit      |\n"
"|                       |                 | function parameter             |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| name_x                | OT_STRINGVECTOR | Names of the variables.        |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| print_header          | OT_BOOL         | Print the header with problem  |\n"
"|                       |                 | statistics                     |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| print_x               | OT_INTVECTOR    | Which variables to print.      |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| qpsol                 | OT_STRING       | The QP solver to be used by    |\n"
"|                       |                 | the SQP method                 |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| qpsol_options         | OT_DICT         | Options to be passed to the QP |\n"
"|                       |                 | solver                         |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| reg_threshold         | OT_DOUBLE       | Threshold for the              |\n"
"|                       |                 | regularization.                |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| regularize            | OT_BOOL         | Automatic regularization of    |\n"
"|                       |                 | Lagrange Hessian.              |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| tol_du                | OT_DOUBLE       | Stopping criterion for dual    |\n"
"|                       |                 | infeasability                  |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| tol_pr                | OT_DOUBLE       | Stopping criterion for primal  |\n"
"|                       |                 | infeasibility                  |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| tol_pr_step           | OT_DOUBLE       | Stopping criterion for the     |\n"
"|                       |                 | step size                      |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"| tol_reg               | OT_DOUBLE       | Stopping criterion for         |\n"
"|                       |                 | regularization                 |\n"
"+-----------------------+-----------------+--------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
