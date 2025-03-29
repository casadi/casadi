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


      #include "cplex_interface.hpp"
      #include <string>

      const std::string casadi::CplexInterface::meta_doc=
      "\n"
"\n"
"\n"
"Interface to Cplex solver for sparse Quadratic Programs\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_22a \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+----------------+-----------------------+---------------------------------+\n"
"|       Id       |         Type          |           Description           |\n"
"+================+=======================+=================================+\n"
"| cplex          | OT_DICT               | Options to be passed to CPLEX   |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| dep_check      | OT_INT                | Detect redundant constraints.   |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| dump_filename  | OT_STRING             | The filename to dump to.        |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| dump_to_file   | OT_BOOL               | Dumps QP to file in CPLEX       |\n"
"|                |                       | format.                         |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| mip_start      | OT_BOOL               | Hot start integers with x0      |\n"
"|                |                       | [Default false].                |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| qp_method      | OT_INT                | Determines which CPLEX          |\n"
"|                |                       | algorithm to use.               |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| sos_groups     | OT_INTVECTORVECTOR    | Definition of SOS groups by     |\n"
"|                |                       | indices.                        |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| sos_types      | OT_INTVECTOR          | Specify 1 or 2 for each SOS     |\n"
"|                |                       | group.                          |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| sos_weights    | OT_DOUBLEVECTORVECTOR | Weights corresponding to SOS    |\n"
"|                |                       | entries.                        |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| tol            | OT_DOUBLE             | Tolerance of solver             |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| version_suffix | OT_STRING             | Specify version of cplex to     |\n"
"|                |                       | load. We will attempt to load l |\n"
"|                |                       | ibcplex<version_suffix>.[so|dll |\n"
"|                |                       | |dylib]. Default value is taken |\n"
"|                |                       | from CPLEX_VERSION env          |\n"
"|                |                       | variable.                       |\n"
"+----------------+-----------------------+---------------------------------+\n"
"| warm_start     | OT_BOOL               | Use warm start with simplex     |\n"
"|                |                       | methods (affects only the       |\n"
"|                |                       | simplex methods).               |\n"
"+----------------+-----------------------+---------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
