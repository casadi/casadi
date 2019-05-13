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

// C++ (and CasADi) from here on
#include <casadi/casadi.hpp>
using namespace casadi;
using namespace std;

int main(){


            Dict options;

            options["inf"] = 1e6;
            Dict hpipm;


/*/
            hpipm["mu0"] = 1e4;            // initial value for duality measure (initial value for complementarity slackness, maximum element in cost functions)
            hpipm["alpha_min"] = 1e-16;     // minimum step size,  exit cond on step length
            hpipm["res_g_max"] = 1e-6;     // exit cond on inf norm of residuals
            hpipm["res_b_max"] = 1e-8;     // exit cond on inf norm of residuals
            hpipm["res_d_max"] = 1e-8;     // exit cond on inf norm of residuals
            hpipm["res_m_max"] = 1e-0;     // exit cond on inf norm of residuals
            hpipm["iter_max"] = 1000;        // exit cond in iter number
            hpipm["stat_max"] = 1000;        // iterations saved in stat
            hpipm["pred_corr"] = 1;        // use Mehrotra"s predictor-corrector IPM algorithm
            hpipm["cond_pred_corr"] = 1;   // conditional Mehrotra"s predictor-corrector
            hpipm["itref_pred_max"] = 0;   // max number of iterative refinement steps for predictor step
            hpipm["itref_corr_max"] = 4;   // max number of iterative refinement steps for corrector step
            hpipm["reg_prim"] = 1e-7;      // reg of primal hessian
            hpipm["lq_fact"] = 0;         // 0 syrk+potrf, 1 mix, 2 lq
            hpipm["lam_min"] = 1e-1;       // min value in lam vector
            hpipm["t_min"] = 1e-1;         // min value in t vector
            hpipm["warm_start"] = 0;       // 0 no warm start, 1 warm start primal sol
            hpipm["abs_form"] = 1;         // absolute IPM formulation
            hpipm["comp_dual_sol"] = 0;    // dual solution (only for abs_form==1)
            hpipm["comp_res_exit"] = 0;    // compute residuals on exit (only for abs_form==1 and comp_dual_sol==1)
            */
            
            options["hpipm"] = hpipm;

            DMDict args;
            std::vector<std::string> keys = {"a","h","g","p","q","x0","lam_x0","lam_a0","lba","uba","lbx","ubx"};
            for (auto k : keys) {
              args[k] = DM::from_file("/home/jgillis/Downloads/dumped_files/solver.000000.in." + k + ".mtx");
            }

            SpDict qp_struct;
            qp_struct["h"] = args["h"].sparsity();
            qp_struct["a"] = args["a"].sparsity();
            
            auto solver = conic("solver","hpipm",qp_struct, options);
            solver(args);

            


  return 0;
}
