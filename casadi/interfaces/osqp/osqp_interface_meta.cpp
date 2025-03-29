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


      #include "osqp_interface.hpp"
      #include <string>

      const std::string casadi::OsqpInterface::meta_doc=
      "\n"
"\n"
"\n"
"Interface to the OSQP Solver for quadratic programming\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_220 \n"
"\n"
"\n"
"Interface to the PROXQP Solver for quadratic programming\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_243 \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-------------------+---------+--------------------------------------------+\n"
"|        Id         |  Type   |                Description                 |\n"
"+===================+=========+============================================+\n"
"| osqp              | OT_DICT | const Options to be passed to osqp.        |\n"
"+-------------------+---------+--------------------------------------------+\n"
"| warm_start_dual   | OT_BOOL | Use lam_a0 and lam_x0 input to warmstart   |\n"
"|                   |         | [Default: truw].                           |\n"
"+-------------------+---------+--------------------------------------------+\n"
"| warm_start_primal | OT_BOOL | Use x0 input to warmstart [Default: true]. |\n"
"+-------------------+---------+--------------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
