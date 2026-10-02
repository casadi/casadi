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



      #include "acados_integrator.hpp"
      #include <string>

      const std::string casadi::AcadosInterface::meta_doc=
      "\n"
"Experimental: acados sim (ERK/IRK) as an Integrator plugin.\n"
"p is NOT differentiable by default: derivatives w.r.t. p are structurally zero.\n"
"Pass is_diff_in with true for p to get exact sensitivities w.r.t. p (p enters\n"
"acados as extra controls, in separate inlined calls used only by expressions\n"
"that need them; the hessian then costs O((nx+nu+np)^2) memory).\n"
"Semi-explicit DAEs with scheme irk (z0: guess for the algebraic states).\n"
"zf is not supported: nan (acados only reports z at the start of an interval).\n"
"\n"
;
