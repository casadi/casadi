/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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


      #include "fixed_smith_dle_internal.hpp"
      #include <string>

      const std::string casadi::FixedSmithDleInternal::meta_doc=
      "\n"
"Solving the Discrete Lyapunov Equations with a fixed number of smith\n"
"iterations.\n"
"\n"
"This plugin uses Smith iterations.\n"
"\n"
"Th basic idea is to exploit the fact that the discrete algebraic\n"
"Lyapunov operator f(X) = AXA^T + V has a fixed point when A is stable.\n"
"\n"
"The pure Smith iterations are:\n"
"\n"
"\n"
"\n"
"::\n"
"\n"
"  X_{-1} = 0\n"
"  X_0 = V\n"
"  k = 0\n"
"  while ||X_k - X_{k-1} || < do\n"
"    X_{k+1} = A X_k A^T + V\n"
"    k += 1\n"
"  end\n"
"  \n"
"  P = X_k\n"
"  \n"
"\n"
"\n"
"\n"
"With frequency doubling, we have:\n"
"\n"
"\n"
"\n"
"::\n"
"\n"
"  X_{-1} = 0\n"
"  X_0 = V\n"
"  V_0 = V\n"
"  A_0 = A\n"
"  k = 0\n"
"  while ||X_k - X_{k-1} || < do\n"
"    X_{k+1} = A_k X_k A_k^T + V_k\n"
"    V_{k+1} = A_k V_k A_k^T + V_k\n"
"    A_{k+1} = A_k A_k\n"
"    k += 1\n"
"  end\n"
"  \n"
"  P = X_k\n"
"  \n"
"\n"
"\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+---------------+------------+---------+----------------------------+\n"
"|      Id       |    Type    | Default |        Description         |\n"
"+===============+============+=========+============================+\n"
"| freq_doubling | OT_BOOLEAN | false   | Use frequency doubling     |\n"
"+---------------+------------+---------+----------------------------+\n"
"| iter          | OT_INTEGER | 100     | Number of Smith iterations |\n"
"+---------------+------------+---------+----------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
