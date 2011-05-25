// Ode.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include <string>
#include <map>

#include <casadi/sx/sx_tools.hpp>

class Ode
{

public:
	~Ode(void);
	Ode(std::string _name);

	std::string name;
	int locked;

	// states/actions/params variables/functions
	std::map<std::string,int> states;
	std::map<std::string,int> actions;
	std::map<std::string,int> params;
	
	void addState(std::string _newState);
	void addAction(std::string _newAction);
	void addParam(std::string _newParam);

	int nx(void);
	int nu(void);
	int np(void);
	int nxu(void);

	void (*dxdt)(std::map<std::string,CasADi::SX> &xDot,
				 std::map<std::string,CasADi::SX> state,
				 std::map<std::string,CasADi::SX> action,
				 std::map<std::string,CasADi::SX> param,
				 CasADi::SX t);


	void assertUnlocked(void);


	int isState(std::string stateName);
	int isAction(std::string actionName);
	int isParam(std::string paramName);

	CasADi::SXMatrix dxVectorDt( CasADi::SXMatrix x, CasADi::SXMatrix u, CasADi::SXMatrix p, CasADi::SX t );
	CasADi::SXMatrix rk4Step( CasADi::SXMatrix x0Vec,
							  CasADi::SXMatrix u0Vec,
							  CasADi::SXMatrix u1Vec,
							  CasADi::SXMatrix pVec,
							  CasADi::SX t0,
							  CasADi::SX t1);

private:
	void assertUniqueName(std::string newName);

	
};
