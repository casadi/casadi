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


      #include "cvodes_interface.hpp"
      #include <string>

      const std::string casadi::CvodesInterface::meta_doc=
      "\n"
"\n"
"\n"
"Interface to CVodes from the Sundials suite.\n"
"\n"
"A call to evaluate will integrate to the end.\n"
"\n"
"You can retrieve the entire state trajectory as follows, after the \n"
"evaluate call: Call reset. Then call integrate(t_i) and getOuput for a\n"
" series of times t_i.\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_228 \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------------------+-----------+--------------------------------+\n"
"|             Id              |   Type    |          Description           |\n"
"+=============================+===========+================================+\n"
"| abstol                      | OT_DOUBLE | Absolute tolerence for the IVP |\n"
"|                             |           | solution                       |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| always_recalculate_jacobian | OT_BOOL   | Recalculate Jacobian before    |\n"
"|                             |           | factorizations, even if        |\n"
"|                             |           | Jacobian is current [default:  |\n"
"|                             |           | true]                          |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| disable_internal_warnings   | OT_BOOL   | Disable SUNDIALS internal      |\n"
"|                             |           | warning messages               |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| fsens_all_at_once           | OT_BOOL   | Calculate all right hand sides |\n"
"|                             |           | of the sensitivity equations   |\n"
"|                             |           | at once                        |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| fsens_err_con               | OT_BOOL   | include the forward            |\n"
"|                             |           | sensitivities in all error     |\n"
"|                             |           | controls                       |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| interpolation_type          | OT_STRING | Type of interpolation for the  |\n"
"|                             |           | adjoint sensitivities          |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| linear_multistep_method     | OT_STRING | Integrator scheme: BDF|adams   |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| linear_solver               | OT_STRING | A custom linear solver creator |\n"
"|                             |           | function [default: qr]         |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| linear_solver_options       | OT_DICT   | Options to be passed to the    |\n"
"|                             |           | linear solver                  |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| max_krylov                  | OT_INT    | Maximum Krylov subspace size   |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| max_multistep_order         | OT_INT    | Maximum order for the          |\n"
"|                             |           | (variable-order) multistep     |\n"
"|                             |           | method                         |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| max_num_steps               | OT_INT    | Maximum number of integrator   |\n"
"|                             |           | steps                          |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| max_order                   | OT_DOUBLE | Maximum order                  |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| max_step_size               | OT_DOUBLE | Max step size [default: 0/inf] |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| min_step_size               | OT_DOUBLE | Min step size [default: 0/0.0] |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| newton_scheme               | OT_STRING | Linear solver scheme in the    |\n"
"|                             |           | Newton method:                 |\n"
"|                             |           | DIRECT|gmres|bcgstab|tfqmr     |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| nonlin_conv_coeff           | OT_DOUBLE | Coefficient in the nonlinear   |\n"
"|                             |           | convergence test               |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| nonlinear_solver_iteration  | OT_STRING | Nonlinear solver type:         |\n"
"|                             |           | NEWTON|functional              |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| quad_err_con                | OT_BOOL   | Should the quadratures affect  |\n"
"|                             |           | the step size control          |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| reltol                      | OT_DOUBLE | Relative tolerence for the IVP |\n"
"|                             |           | solution                       |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| scale_abstol                | OT_BOOL   | Scale absolute tolerance by    |\n"
"|                             |           | nominal value                  |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| second_order_correction     | OT_BOOL   | Second order correction in the |\n"
"|                             |           | augmented system Jacobian      |\n"
"|                             |           | [true]                         |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| sensitivity_method          | OT_STRING | Sensitivity method:            |\n"
"|                             |           | SIMULTANEOUS|staggered         |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| step0                       | OT_DOUBLE | initial step size [default:    |\n"
"|                             |           | 0/estimated]                   |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| steps_per_checkpoint        | OT_INT    | Number of steps between two    |\n"
"|                             |           | consecutive checkpoints        |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| stop_at_end                 | OT_BOOL   | [DEPRECATED] Stop the          |\n"
"|                             |           | integrator at the end of the   |\n"
"|                             |           | interval                       |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"| use_preconditioner          | OT_BOOL   | Precondition the iterative     |\n"
"|                             |           | solver [default: true]         |\n"
"+-----------------------------+-----------+--------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
