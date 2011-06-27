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


class Lqr
{
public:
	Lqr(Ode & _ode, double t0_, double tf_, int N_,
		CasADi::SX (*cost_)(std::map<std::string,CasADi::SX> state,
							std::map<std::string,CasADi::SX> action));
	~Lqr();

	CasADi::SXMatrix takeBackwardsStep(CasADi::SXMatrix Vxx_kp1, CasADi::SXMatrix Vx_kp1, CasADi::SXMatrix V0_kp1);
	
	// void writeOctaveOutput( std::ofstream & f, double * xOpt);
	
private:
	double t0;
	double tf;
	int N;

	Ode & ode;
	
	// cost function/gradients/hessians
	CasADi::SXFunction cost;
	CasADi::SXFunction cost_x;
	CasADi::SXFunction cost_u;
	CasADi::SXFunction cost_xx;
	CasADi::SXFunction cost_xu;
	CasADi::SXFunction cost_uu;

	// dynamics function/gradients
	CasADi::SXFunction f;
	CasADi::SXFunction f_x; // A matrix
	CasADi::SXFunction f_u; // B matrix
	// CasADi::SXFunction f_xx;
	// CasADi::SXFunction f_xu;
	// CasADi::SXFunction f_uu;

	// Q functions, gradients, hessians
	CasADi::SXFunction Q_0_fcn;
	CasADi::SXFunction Q_x_fcn;
	CasADi::SXFunction Q_u_fcn;
	CasADi::SXFunction Q_xx_fcn;
	CasADi::SXFunction Q_xu_fcn;
	CasADi::SXFunction Q_uu_fcn;

	void setupFunctions(void);



	CasADi::SX (*costFcnExt)(std::map<std::string,CasADi::SX> state,
							 std::map<std::string,CasADi::SX> action);

};
