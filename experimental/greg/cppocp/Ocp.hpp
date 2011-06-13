// Ocp.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include <string>
#include <map>

#include "Ode.hpp"
#include "MultipleShooting.hpp"

#include <casadi/sx/sx_tools.hpp>
#include <casadi/stl_vector_tools.hpp>

class Ocp
{
public:
	Ocp(void);
	~Ocp(void);

	void addNonlconIneq(CasADi::SXMatrix gNew);
	void addNonlconEq(CasADi::SXMatrix gNew);

	CasADi::SXMatrix designVariables;
	std::vector<double>lb;
	std::vector<double>ub;
	std::vector<double>guess;

	CasADi::SX objFun;
	CasADi::SXMatrix g;
	std::vector<double> gMin;
	std::vector<double> gMax;

	MultipleShooting & addMultipleShooting(std::string name, Ode & ode, CasADi::SX t0, CasADi::SX tf, int N);

	CasADi::SX & addParam(std::string _newParam);
	CasADi::SX & getParam(std::string p);
	void boundParam(std::string p, double lb_, double ub_);
	void setParamGuess(std::string p, double guess_);


	void writeMatlabOutput( const char * filename, double * xOpt);

private:

	// multiple shooting instances
	std::map<std::string,MultipleShooting*> ms;

	// params
	std::map<std::string,int> paramsIdx;
	std::map<std::string,CasADi::SX> params;
	int np(void);
	int isParam(std::string paramName);


	CasADi::SXMatrix getParamMat(void);

	int getParamIdx(std::string p);



protected:
};
