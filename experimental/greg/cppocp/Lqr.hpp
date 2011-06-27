// Lqr.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include "Ode.hpp"

#include <string>
#include <map>
#include <fstream>

#include <casadi/sx/sx_tools.hpp>
#include <casadi/fx/sx_function.hpp>

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
#define NUM_BACKWARD_SWEEP_INPUTS           11

#define IDX_BACKWARD_SWEEP_OUTPUTS_U_FEEDFORWARD_K   0
#define IDX_BACKWARD_SWEEP_OUTPUTS_U_FEEDBACK_GAIN_K 1
#define IDX_BACKWARD_SWEEP_OUTPUTS_V_0_K             2
#define IDX_BACKWARD_SWEEP_OUTPUTS_V_X_K             3
#define IDX_BACKWARD_SWEEP_OUTPUTS_V_XX_K            4
#define NUM_BACKWARD_SWEEP_OUTPUTS                   5


#define IDX_COST_INPUTS_X_K 0
#define IDX_COST_INPUTS_U_K 1
#define NUM_COST_INPUTS     2

#define IDX_COST_OUTPUTS_COST_0_K  0
#define IDX_COST_OUTPUTS_COST_X_K  1
#define IDX_COST_OUTPUTS_COST_U_K  2
#define IDX_COST_OUTPUTS_COST_XX_K 3
#define IDX_COST_OUTPUTS_COST_XU_K 4
#define IDX_COST_OUTPUTS_COST_UU_K 5
#define NUM_COST_OUTPUTS           6

class Lqr
{
public:
	Lqr(Ode & _ode, double t0_, double tf_, int N_,
		CasADi::SX (*cost_)(std::map<std::string,CasADi::SX> state,
							std::map<std::string,CasADi::SX> action,
							int timestep,
							int N_));
	~Lqr();

	void runBackwardSweep(void);
	void runForwardSweep(void);

	// trajectory
	std::vector<CasADi::DMatrix> V_0;
	std::vector<CasADi::DMatrix> V_x;
	std::vector<CasADi::DMatrix> V_xx;
	std::vector<CasADi::DMatrix> u_feedforward;
	std::vector<CasADi::DMatrix> u_feedback_gain;
	std::vector<CasADi::DMatrix> x_trajectory;
	std::vector<CasADi::DMatrix> u_trajectory;
	std::vector<CasADi::DMatrix> cost_0;
	std::vector<CasADi::DMatrix> cost_x;
	std::vector<CasADi::DMatrix> cost_u;
	std::vector<CasADi::DMatrix> cost_xx;
	std::vector<CasADi::DMatrix> cost_xu;
	std::vector<CasADi::DMatrix> cost_uu;
	
private:
	double t0;
	double tf;
	int N;

	Ode & ode;

	// SXFunctions
	std::vector<CasADi::SXFunction> costFunctions;
	CasADi::SXFunction backwardSweepFcn;
	CasADi::SXFunction dynamicsFcn;
	CasADi::SXFunction qFcn;

	void setupBackwardSweepFunction(void);
	void setupCostFunctions(void);
	void takeBackwardStep(int timestep);
	void takeForwardStep(int timestep);

	CasADi::SX (*costFcnExt)(std::map<std::string,CasADi::SX> state,
							 std::map<std::string,CasADi::SX> action,
							 int timestep,
							 int N_);

};
