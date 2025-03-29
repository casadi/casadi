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


      #include "fatrop_interface.hpp"
      #include <string>

      const std::string casadi::FatropInterface::meta_doc=
      "\n"
"\n"
">List of available options\n"
"\n"
"+---------------------+--------------+-------------------------------------+\n"
"|         Id          |     Type     |             Description             |\n"
"+=====================+==============+=====================================+\n"
"| N                   | OT_INT       | OCP horizon                         |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| convexify_margin    | OT_DOUBLE    | When using a convexification        |\n"
"|                     |              | strategy, make sure that the        |\n"
"|                     |              | smallest eigenvalue is at least     |\n"
"|                     |              | this (default: 1e-7).               |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| convexify_strategy  | OT_STRING    | NONE|regularize|eigen-              |\n"
"|                     |              | reflect|eigen-clip. Strategy to     |\n"
"|                     |              | convexify the Lagrange Hessian      |\n"
"|                     |              | before passing it to the solver.    |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| debug               | OT_BOOL      | Produce debug information (default: |\n"
"|                     |              | false)                              |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| fatrop              | OT_DICT      | Options to be passed to fatrop      |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| ng                  | OT_INTVECTOR | Number of non-dynamic constraints,  |\n"
"|                     |              | length N+1                          |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| nu                  | OT_INTVECTOR | Number of controls, length N+1      |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| nx                  | OT_INTVECTOR | Number of states, length N+1        |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| structure_detection | OT_STRING    | NONE | auto | manual                |\n"
"+---------------------+--------------+-------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
