// OcpMultipleShooting.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include "Ocp.hpp"

#include <string>

class OcpMultipleShooting : public Ocp
{
public:
	OcpMultipleShooting(Ode * _ode);
	~OcpMultipleShooting();
	void setTimeInterval(CasADi::SX _t0, CasADi::SX _tf);

	int N;
	CasADi::SX t0;
	CasADi::SX tf;
	int getBigN(void);
	void discretize(int _N);

	std::vector<CasADi::SX>designVariables;
	std::vector<double>lb;
	std::vector<double>ub;
	std::vector<double>guess;

	CasADi::SXMatrix getStateMat(int timeStep);
	CasADi::SXMatrix getActionMat(int timeStep);
	CasADi::SXMatrix getParamMat(void);

	CasADi::SX getState(std::string x, int timeStep);
	CasADi::SX getAction(std::string u, int timeStep);
	CasADi::SX getParam(std::string p);

	void boundStateAction(std::string xu, double _lb, double _ub, int timeStep);
	void boundParam(std::string p, double _lb, double _ub);

private:

	int getStateActionIdx(std::string xu, int timeStep);
	int getParamIdx(std::string p);

	
};
