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


      #include "newton.hpp"
      #include <string>

      const std::string casadi::Newton::meta_doc=
      "\n"
"\n"
"\n"
"Implements simple newton iterations to solve an implicit function.\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_236 \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------+--------------------------------------------+\n"
"|       Id        |   Type    |                Description                 |\n"
"+=================+===========+============================================+\n"
"| abstol          | OT_DOUBLE | Stopping criterion tolerance on max(|F|)   |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| abstolStep      | OT_DOUBLE | Stopping criterion tolerance on step size  |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| line_search     | OT_BOOL   | Enable line-search (default: true)         |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| max_iter        | OT_INT    | Maximum number of Newton iterations to     |\n"
"|                 |           | perform before returning.                  |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"| print_iteration | OT_BOOL   | Print information about each iteration     |\n"
"+-----------------+-----------+--------------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
