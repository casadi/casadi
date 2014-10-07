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

// Lqr.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include "Ode.hpp"

#include <string>
#include <map>
#include <fstream>

#include <core/sx/sx_tools.hpp>
#include <core/function/sx_function.hpp>

// #define EVAL_DYNAMICS_FUNCTION
// #define EVAL_Q_FUNCTION

#define IDX_BACKWARD_SWEEP_INPUTS_X_K       0
#define IDX_BACKWARD_SWEEP_INPUTS_U_K       1
#define IDX_BACKWARD_SWEEP_INPUTS_COST_0_K  2
#define IDX_BACKWARD_SWEEP_INPUTS_COST_X_K  3
#define IDX_BACKWARD_SWEEP_INPUTS_COST_U_K  4
#define IDX_BACKWARD_SWEEP_INPUTS_COST_XX_K 5
#define IDX_BACKWARD_SWEEP_INPUTS_COST_XU_K 6
#define IDX_BACKWARD_SWEEP_INPUTS_COST_UU_K 7
#define IDX_BACKWARD_SWEEP_INPUTS_V_0_KP1   8
#define IDX_BACKWARD_SWEEP_INPUTS_V_X_KP1   9
#define IDX_BACKWARD_SWEEP_INPUTS_V_XX_KP1  10
#define NUM_BACKWARD_SWEEP_INPUTS 11

#define IDX_BACKWARD_SWEEP_OUTPUTS_U_OPEN_LOOP_K   0
#define IDX_BACKWARD_SWEEP_OUTPUTS_FEEDBACK_GAIN_K 1
#define IDX_BACKWARD_SWEEP_OUTPUTS_V_0_K           2
#define IDX_BACKWARD_SWEEP_OUTPUTS_V_X_K           3
#define IDX_BACKWARD_SWEEP_OUTPUTS_V_XX_K          4
#define NUM_BACKWARD_SWEEP_OUTPUTS 5

#define IDX_COST_INPUTS_X_K 0
#define IDX_COST_INPUTS_U_K 1
#define NUM_COST_INPUTS 2

#define IDX_COST_OUTPUTS_COST_0_K  0
#define IDX_COST_OUTPUTS_COST_X_K  1
#define IDX_COST_OUTPUTS_COST_U_K  2
#define IDX_COST_OUTPUTS_COST_XX_K 3
#define IDX_COST_OUTPUTS_COST_XU_K 4
#define IDX_COST_OUTPUTS_COST_UU_K 5
#define NUM_COST_OUTPUTS 6

class Lqr
{
public:
     Lqr(Ode & _ode, double t0_, double tf_, int N_,
	 casadi::SX (*cost_)(std::map<std::string,casadi::SX> state,
			     std::map<std::string,casadi::SX> action,
			     int timestep,
			     int N_));
     ~Lqr();

     void runBackwardSweep(void);
     void runForwardSweep(void);

     // trajectory
     std::vector<casadi::DMatrix> V_0;
     std::vector<casadi::DMatrix> V_x;
     std::vector<casadi::DMatrix> V_xx;
     std::vector<casadi::DMatrix> cost_0;
     std::vector<casadi::DMatrix> cost_x;
     std::vector<casadi::DMatrix> cost_u;
     std::vector<casadi::DMatrix> cost_xx;
     std::vector<casadi::DMatrix> cost_xu;
     std::vector<casadi::DMatrix> cost_uu;
     std::vector<casadi::DMatrix> xTrajectory;
     std::vector<casadi::DMatrix> xNominalTrajectory;
     std::vector<casadi::DMatrix> uTrajectory;
     std::vector<casadi::DMatrix> uOpenLoop;
     std::vector<casadi::DMatrix> feedbackGain;

     // regularization
     casadi::DMatrix stateRegularization;
     casadi::DMatrix actionRegularization;

     // for debugging:
     std::vector<casadi::DMatrix> Q0Trajectory;
     std::vector<casadi::DMatrix> QxTrajectory;
     std::vector<casadi::DMatrix> QuTrajectory;
     std::vector<casadi::DMatrix> QxxTrajectory;
     std::vector<casadi::DMatrix> QxuTrajectory;
     std::vector<casadi::DMatrix> QuuTrajectory;

     std::vector<casadi::DMatrix> f0Trajectory;
     std::vector<casadi::DMatrix> functionTrajectory;
     std::vector<casadi::DMatrix> fuTrajectory;

     // set actions bounds
     void boundAction( std::vector<double> lb_, std::vector<double> ub_ );
     void boundAction( std::vector<double> lb_, std::vector<double> ub_, int timestep );

     int N;

private:
     double t0;
     double tf;

     std::vector< std::vector<double> > ubAction;
     std::vector< std::vector<double> > lbAction;

     Ode & ode;

     // SXFunctions
     std::vector<casadi::SXFunction> costFunctions;
     casadi::SXFunction backwardSweepFcn;
     casadi::SXFunction dynamicsFcn;
     casadi::SXFunction qFcn;

     void setupBackwardSweepFunction(void);
     void setupCostFunctions(void);
     void takeBackwardStep(int timestep);
     void takeForwardStep(int timestep);

     casadi::SX (*costFcnExt)(std::map<std::string,casadi::SX> state,
			      std::map<std::string,casadi::SX> action,
			      int timestep,
			      int N_);

};
