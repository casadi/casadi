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

#define IDX_INPUTS_X_K      0
#define IDX_INPUTS_U_K      1
#define IDX_INPUTS_V_0_KP1  2
#define IDX_INPUTS_V_X_KP1  3
#define IDX_INPUTS_V_XX_KP1 4
#define LQR_NUM_INPUTS      5

#define IDX_LQR_OUTPUTS_U_FEEDFORWARD_K   0
#define IDX_LQR_OUTPUTS_U_FEEDBACK_GAIN_K 1
#define IDX_LQR_OUTPUTS_V_0_K             2
#define IDX_LQR_OUTPUTS_V_X_K             3
#define IDX_LQR_OUTPUTS_V_XX_K            4
#define LQR_NUM_OUTPUTS                   5

class Lqr
{
public:
	Lqr(Ode & _ode, double t0_, double tf_, int N_,
		CasADi::SX (*cost_)(std::map<std::string,CasADi::SX> state,
							std::map<std::string,CasADi::SX> action));
	~Lqr();

	void runBackwardSweep(void);
	void runForwardSweep(void);
	
private:
	double t0;
	double tf;
	int N;

	Ode & ode;

	// trajectory
	std::vector<CasADi::DMatrix> V_0;
	std::vector<CasADi::DMatrix> V_x;
	std::vector<CasADi::DMatrix> V_xx;
	std::vector<CasADi::DMatrix> u_feedforward;
	std::vector<CasADi::DMatrix> u_feedback_gain;
	std::vector<CasADi::DMatrix> x_trajectory;
	std::vector<CasADi::DMatrix> u_trajectory;

	// SXFunctions
	CasADi::SXFunction ilqr_fcn;
	// CasADi::SXFunction cost_fcn;
	// CasADi::SXFunction dynamics_fcn;
	// CasADi::SXFunction Q_fcn;

	void setupFunctions(void);
	void takeBackwardStep(int timestep);
	void takeForwardStep(int timestep);

	CasADi::SX (*costFcnExt)(std::map<std::string,CasADi::SX> state,
							 std::map<std::string,CasADi::SX> action);

};
