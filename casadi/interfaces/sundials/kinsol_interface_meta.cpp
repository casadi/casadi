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


      #include "kinsol_interface.hpp"
      #include <string>

      const std::string casadi::KinsolInterface::meta_doc=
      "\n"
"\n"
"\n"
"KINSOL interface from the Sundials suite\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_226 \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+---------------------------+-----------------+----------------------------+\n"
"|            Id             |      Type       |        Description         |\n"
"+===========================+=================+============================+\n"
"| abstol                    | OT_DOUBLE       | Stopping criterion         |\n"
"|                           |                 | tolerance                  |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| disable_internal_warnings | OT_BOOL         | Disable KINSOL internal    |\n"
"|                           |                 | warning messages           |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| exact_jacobian            | OT_BOOL         | Use exact Jacobian         |\n"
"|                           |                 | information                |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| f_scale                   | OT_DOUBLEVECTOR | Equation scaling factors   |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| iterative_solver          | OT_STRING       | gmres|bcgstab|tfqmr        |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| linear_solver_type        | OT_STRING       | dense|banded|iterative|use |\n"
"|                           |                 | r_defined                  |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| lower_bandwidth           | OT_INT          | Lower bandwidth for banded |\n"
"|                           |                 | linear solvers             |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| max_iter                  | OT_INT          | Maximum number of Newton   |\n"
"|                           |                 | iterations. Putting 0 sets |\n"
"|                           |                 | the default value of       |\n"
"|                           |                 | KinSol.                    |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| max_krylov                | OT_INT          | Maximum Krylov space       |\n"
"|                           |                 | dimension                  |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| pretype                   | OT_STRING       | Type of preconditioner     |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| print_level               | OT_INT          | Verbosity level            |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| strategy                  | OT_STRING       | Globalization strategy     |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| u_scale                   | OT_DOUBLEVECTOR | Variable scaling factors   |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| upper_bandwidth           | OT_INT          | Upper bandwidth for banded |\n"
"|                           |                 | linear solvers             |\n"
"+---------------------------+-----------------+----------------------------+\n"
"| use_preconditioner        | OT_BOOL         | Precondition an iterative  |\n"
"|                           |                 | solver                     |\n"
"+---------------------------+-----------------+----------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
