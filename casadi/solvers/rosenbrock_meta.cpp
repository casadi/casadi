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


      #include "rosenbrock.hpp"
      #include <string>

      const std::string casadi::Rosenbrock::meta_doc=
      "\n"
"\n"
"\n"
"Rosenbrock-Wanner integrator for stiff ODEs and semi-explicit index-1\n"
" DAEs, with adaptive or fixed steps.\n"
"\n"
"An 8-stage, stiffly accurate, L-stable method of order 5 with an \n"
"embedded method of order 4. Each step needs one Jacobian and one \n"
"sparse QR factorization; the method is linearly implicit, without \n"
"Newton iterations. Steps cross the output times within a control \n"
"group, which are interpolated by a dense output of order 4. Step sizes\n"
" are chosen from the nominal solution only, so forward sensitivities \n"
"are the exact derivatives of the discrete solution (internal numerical\n"
" differentiation) and adjoint sensitivities are the exact discrete \n"
"adjoint of the taped steps. With step_size, each run of steps over a \n"
"control group is divided into equal steps instead, and the tape has a \n"
"known bound.\n"
"\n"
"Zero-crossing events (DAE output zero, option transition) are located \n"
"on the dense output of the step in which a function becomes negative; \n"
"integration restarts from the state after the transition. \n"
"Sensitivities, forward, adjoint and forward-over-adjoint, include \n"
"those of the event times.\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_2ka \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+--------------------+--------------+--------------------------------------+\n"
"|         Id         |     Type     |             Description              |\n"
"+====================+==============+======================================+\n"
"| abstol             | OT_DOUBLE    | Absolute tolerance [1e-8]            |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| control_group_size | OT_INTVECTOR | Number of output intervals in each   |\n"
"|                    |              | control group, over which the        |\n"
"|                    |              | controls, one column each, are       |\n"
"|                    |              | constant: a run of steps [default: 1 |\n"
"|                    |              | for each output interval]            |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| init_max_iter      | OT_INT       | Maximum number of Newton iterations  |\n"
"|                    |              | for consistent algebraic states [20] |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| init_tol           | OT_DOUBLE    | Tolerance on the scaled Newton step  |\n"
"|                    |              | when making the algebraic states     |\n"
"|                    |              | consistent [1e-3]                    |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| max_num_steps      | OT_INT       | Adaptive steps: maximum number of    |\n"
"|                    |              | steps per output interval [10000]    |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| max_step_size      | OT_DOUBLE    | Adaptive steps: maximum step size    |\n"
"|                    |              | [default: unbounded]                 |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| max_tape_steps     | OT_INT       | Adaptive steps: capacity, in steps,  |\n"
"|                    |              | of a tape in the work vector, for    |\n"
"|                    |              | adjoint sensitivities; needed for    |\n"
"|                    |              | code generation [0: tape of          |\n"
"|                    |              | unbounded size]. With fixed steps,   |\n"
"|                    |              | the tape is sized exactly.           |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| min_step_size      | OT_DOUBLE    | Adaptive steps: minimum step size    |\n"
"|                    |              | [1e-14]                              |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| quad_err_con       | OT_BOOL      | Adaptive steps: include the          |\n"
"|                    |              | quadratures in the error control     |\n"
"|                    |              | [true]                               |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| reltol             | OT_DOUBLE    | Relative tolerance [1e-6]            |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| step0              | OT_DOUBLE    | Adaptive steps: initial step size    |\n"
"|                    |              | [default: estimated]                 |\n"
"+--------------------+--------------+--------------------------------------+\n"
"| step_size          | OT_DOUBLE    | Fixed steps: each run of steps, over |\n"
"|                    |              | a control group or up to the end     |\n"
"|                    |              | without controls, is divided into    |\n"
"|                    |              | equal steps of at most step_size [0: |\n"
"|                    |              | adaptive steps]                      |\n"
"+--------------------+--------------+--------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
