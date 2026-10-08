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


      #include "daqp_interface.hpp"
      #include <string>

      const std::string casadi::DaqpInterface::meta_doc=
      "\n"
"\n"
"\n"
"Interface to DAQP for convex quadratic programs with optional binary variables.\n"
"See https://darnstrom.github.io/daqp/parameters for solver settings.\n"
"Set warm_start=true to use x0 as a primal start or MIQP incumbent.\n"
"lam_x0 and lam_a0 initialize the active set when supplied; otherwise x0 does.\n"
"Workspace allocations and unchanged matrix factors are reused between calls.\n"
"Each solve starts cold unless warm_start_previous=true or an explicit start is used.\n"
"Explicit warm starts override the retained active set.\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_29l \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+---------------------+---------+------------------------------------------+\n"
"| Id                  | Type    | Description                              |\n"
"+=====================+=========+==========================================+\n"
"| daqp                | OT_DICT | Settings passed to DAQP.                 |\n"
"+---------------------+---------+------------------------------------------+\n"
"| warm_start          | OT_BOOL | Use primal/dual starts (default: false). |\n"
"+---------------------+---------+------------------------------------------+\n"
"| warm_start_previous | OT_BOOL | Reuse previous solve state (false).      |\n"
"+---------------------+---------+------------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
