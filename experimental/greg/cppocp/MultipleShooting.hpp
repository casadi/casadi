// MultipleShooting.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include "Ode.hpp"

#include <string>
#include <map>

class MultipleShooting
{
public:
	MultipleShooting(std::string name_,
					 Ode & _ode,
					 CasADi::SX t0_,
					 CasADi::SX tf_,
					 int N_,
					 std::vector<double>&lb_,
					 std::vector<double>&ub_,
					 std::vector<double>&guess_,
					 CasADi::SXMatrix & dv_,
					 int idx0_);
	~MultipleShooting();

	int getBigN(void);

	void boundStateAction(std::string xu, double lb_, double ub_, int timeStep);
	void boundStateAction(std::string xu, double lb_, double ub_);
    void setStateActionGuess(std::string xu, double guess_, int timeStep);
    void setStateActionGuess(std::string xu, double guess_);

	void writeMatlabOutput( const char * filename, double * xOpt);

	CasADi::SX getState(std::string x, int timeStep);
	CasADi::SX getAction(std::string u, int timeStep);

	CasADi::SXMatrix getStateMat(int timeStep);
	CasADi::SXMatrix getActionMat(int timeStep);

	int N;
	
	CasADi::SXMatrix getDynamicsConstraintError(int timeStep, std::map<std::string,CasADi::SX> params);

private:
	std::string name;

	Ode & ode;
	CasADi::SXMatrix & dv;
	std::vector<double>&lb;
	std::vector<double>&ub;
	std::vector<double>&guess;
	int idx0;

	CasADi::SX t0;
	CasADi::SX tf;

	int getIdx(std::string xu, int timeStep);
	
};
