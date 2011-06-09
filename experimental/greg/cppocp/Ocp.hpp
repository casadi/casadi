// Ocp.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include "Ode.hpp"

#include <casadi/sx/sx_tools.hpp>
#include <casadi/stl_vector_tools.hpp>

#include <interfaces/ipopt/ipopt_solver.hpp>

class Ocp
{
public:
	Ocp(Ode * _ode);
	~Ocp(void);

	void addNonlconIneq(CasADi::SXMatrix gNew);
	void addNonlconEq(CasADi::SXMatrix gNew);

	CasADi::SX objFun;
	CasADi::SXMatrix g;
	std::vector<double> gMin;
	std::vector<double> gMax;

	void solve(void);

private:

protected:
	Ode * ode;

};
