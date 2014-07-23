/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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


      #include "ooqp_internal.hpp"
      #include <string>

      const std::string casadi::OOQPInternal::meta_doc=
      "\n"
"Interface to the OOQP Solver for quadratic programming The current\n"
"implementation assumes that OOQP is configured with the MA27 sparse\n"
"linear solver.\n"
"\n"
"NOTE: when doing multiple calls to evaluate(), check if you need to\n"
"reInit();\n"
"\n"
"\n"
">List of available options\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| artol           | OT_REAL         | 0.000           | tolerance as    |\n"
"|                 |                 |                 | provided with   |\n"
"|                 |                 |                 | setArTol to     |\n"
"|                 |                 |                 | OOQP            |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| mutol           | OT_REAL         | 0.000           | tolerance as    |\n"
"|                 |                 |                 | provided with   |\n"
"|                 |                 |                 | setMuTol to     |\n"
"|                 |                 |                 | OOQP            |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| print_level     | OT_INTEGER      | 0               | Print level.    |\n"
"|                 |                 |                 | OOQP listens to |\n"
"|                 |                 |                 | print_level 0,  |\n"
"|                 |                 |                 | 10 and 100      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
