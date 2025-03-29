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


      #include "symbolic_qr.hpp"
      #include <string>

      const std::string casadi::SymbolicQr::meta_doc=
      "\n"
"\n"
"\n"
"Linear solver for sparse least-squares problems Inspired from https://github.com/scipy/scipy/blob/v0.14.0/scipy/sparse/linalg/isolve/lsqr.py#L96 \n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_230 \n"
"\n"
"\n"
"\n"
"Linsol based on QR factorization with sparsity pattern based \n"
"reordering  without partial pivoting\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_231 \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-------+---------+----------------------------------------------------+\n"
"|  Id   |  Type   |                    Description                     |\n"
"+=======+=========+====================================================+\n"
"| fopts | OT_DICT | Options to be passed to generated function objects |\n"
"+-------+---------+----------------------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
