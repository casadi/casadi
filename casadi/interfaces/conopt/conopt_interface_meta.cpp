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


      #include "conopt_interface.hpp"
      #include <string>

      const std::string casadi::ConoptInterface::meta_doc=
      "\n"
"\n"
"\n"
"CONOPT interface\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+---------------+-----------+-----------------------------------+\n"
"|      Id       |   Type    |            Description            |\n"
"+===============+===========+===================================+\n"
"| conopt        | OT_DICT   | Options to be passed to CONOPT    |\n"
"+---------------+-----------+-----------------------------------+\n"
"| debug         | OT_BOOL   | Print debug output: constraint    |\n"
"|               |           | values at each FDEval, solution   |\n"
"|               |           | vector, and option echo           |\n"
"+---------------+-----------+-----------------------------------+\n"
"| exact_hessian | OT_BOOL   | Provide exact Hessian to CONOPT   |\n"
"+---------------+-----------+-----------------------------------+\n"
"| optfile       | OT_STRING | Path to a CONOPT option file (for |\n"
"|               |           | string-valued CR-cells such as    |\n"
"|               |           | Algorithm)                        |\n"
"+---------------+-----------+-----------------------------------+\n"
"| warm_start    | OT_BOOL   | Warm-start CONOPT using           |\n"
"|               |           | multipliers from a prior solve to |\n"
"|               |           | infer basis status (IniStat=2)    |\n"
"+---------------+-----------+-----------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
