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


      #include "hpipm_interface.hpp"
      #include <string>

      const std::string casadi::HpipmInterface::meta_doc=
      "\n"
"\n"
"Interface to HPIPM Solver\n"
"\n"
"In order to use this interface, you must:\n"
"\n"
"\n"
"Decision variables must only by state and control, and the variable \n"
"ordering must be [x0 u0 x1 u1 ...]\n"
"\n"
"The constraints must be in order: [ gap0 lincon0 gap1 lincon1 ]\n"
"\n"
"gap: Ak+1 = Ak xk + Bk uk lincon: yk= Ck xk + Dk uk\n"
"\n"
"\n"
"\n"
"\n"
"::\n"
"\n"
"         A0 B0 -I\n"
"         C0 D0\n"
"                A1 B1 -I\n"
"                C1 D1\n"
"\n"
"\n"
"\n"
"where I must be a diagonal sparse matrix\n"
"Either supply all of N, nx, ng, nu options or rely on automatic \n"
"detection\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_242 \n"
"\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-------+--------------+--------------------------------------------------+\n"
"|  Id   |     Type     |                   Description                    |\n"
"+=======+==============+==================================================+\n"
"| N     | OT_INT       | OCP horizon                                      |\n"
"+-------+--------------+--------------------------------------------------+\n"
"| hpipm | OT_DICT      | Options to be passed to hpipm                    |\n"
"+-------+--------------+--------------------------------------------------+\n"
"| inf   | OT_DOUBLE    | Replace infinities by this amount [default: 1e8] |\n"
"+-------+--------------+--------------------------------------------------+\n"
"| ng    | OT_INTVECTOR | Number of non-dynamic constraints, length N+1    |\n"
"+-------+--------------+--------------------------------------------------+\n"
"| nu    | OT_INTVECTOR | Number of controls, length N                     |\n"
"+-------+--------------+--------------------------------------------------+\n"
"| nx    | OT_INTVECTOR | Number of states, length N+1                     |\n"
"+-------+--------------+--------------------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
