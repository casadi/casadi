// Ddp.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include "Ode.hpp"

#include <string>
#include <map>
#include <fstream>

#include <symbolic/sx/sx_tools.hpp>
#include <symbolic/function/sx_function.hpp>

#define IDX_BACKWARD_SWEEP_INPUTS_U_K    0
#define IDX_BACKWARD_SWEEP_INPUTS_Q_0_K  1
#define IDX_BACKWARD_SWEEP_INPUTS_Q_X_K  2
#define IDX_BACKWARD_SWEEP_INPUTS_Q_U_K  3
#define IDX_BACKWARD_SWEEP_INPUTS_Q_XX_K 4
#define IDX_BACKWARD_SWEEP_INPUTS_Q_XU_K 5
#define IDX_BACKWARD_SWEEP_INPUTS_Q_UU_K 6
#define NUM_BACKWARD_SWEEP_INPUTS 7

#define IDX_BACKWARD_SWEEP_OUTPUTS_U_OPEN_LOOP_K   0
#define IDX_BACKWARD_SWEEP_OUTPUTS_FEEDBACK_GAIN_K 1
#define IDX_BACKWARD_SWEEP_OUTPUTS_V_0_K           2
#define IDX_BACKWARD_SWEEP_OUTPUTS_V_X_K           3
#define IDX_BACKWARD_SWEEP_OUTPUTS_V_XX_K          4
#define NUM_BACKWARD_SWEEP_OUTPUTS 5


#define IDX_Q_INPUTS_X_K      0
#define IDX_Q_INPUTS_U_K      1
#define IDX_Q_INPUTS_V_0_KP1  2
#define IDX_Q_INPUTS_V_X_KP1  3
#define IDX_Q_INPUTS_V_XX_KP1 4
#define NUM_Q_INPUTS 5

#define IDX_Q_OUTPUTS_Q_0_K  0
#define IDX_Q_OUTPUTS_Q_X_K  1
#define IDX_Q_OUTPUTS_Q_U_K  2
#define IDX_Q_OUTPUTS_Q_XX_K 3
#define IDX_Q_OUTPUTS_Q_XU_K 4
#define IDX_Q_OUTPUTS_Q_UU_K 5
#define NUM_Q_OUTPUTS 6

class Ddp
{
public:
     Ddp(Ode & _ode, double t0_, double tf_, int N_,
	 CasADi::SX (*cost_)(std::map<std::string,CasADi::SX> state,
			     std::map<std::string,CasADi::SX> action,
			     int timestep,
			     int N_));
     ~Ddp();

     void runBackwardSweep(void);
     void runForwardSweep(void);

     // trajectory
     std::vector<CasADi::DMatrix> V_0;
     std::vector<CasADi::DMatrix> V_x;
     std::vector<CasADi::DMatrix> V_xx;
     std::vector<CasADi::DMatrix> xTrajectory;
     std::vector<CasADi::DMatrix> xNominalTrajectory;
     std::vector<CasADi::DMatrix> uTrajectory;
     std::vector<CasADi::DMatrix> uOpenLoop;
     std::vector<CasADi::DMatrix> feedbackGain;

     // regularization
     CasADi::DMatrix stateRegularization;
     CasADi::DMatrix actionRegularization;

     // for debugging:
     std::vector<CasADi::DMatrix> Q0Trajectory;
     std::vector<CasADi::DMatrix> QxTrajectory;
     std::vector<CasADi::DMatrix> QuTrajectory;
     std::vector<CasADi::DMatrix> QxxTrajectory;
     std::vector<CasADi::DMatrix> QxuTrajectory;
     std::vector<CasADi::DMatrix> QuuTrajectory;

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
     std::vector<CasADi::SXFunction> qFunctions;
     CasADi::SXFunction backwardSweepFcn;
     CasADi::SXFunction dynamicsFcn;
     CasADi::SXFunction qFcn;

     void setupBackwardSweepFunction(void);
     void setupQFunctions(void);
     void takeBackwardStep(int timestep);
     void takeForwardStep(int timestep);

     CasADi::SX (*costFcnExt)(std::map<std::string,CasADi::SX> state,
			      std::map<std::string,CasADi::SX> action,
			      int timestep,
			      int N_);

};
