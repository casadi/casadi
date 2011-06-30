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
	void addNonlconIneq(CasADi::SXMatrix gNew, std::string name);	
	void addNonlconEq(CasADi::SXMatrix gNew);
	void addNonlconEq(CasADi::SXMatrix gNew, std::string name);

	CasADi::SXMatrix designVariables;
	std::vector<double>lb;
	std::vector<double>ub;
	std::vector<double>guess;

	CasADi::SX objFun;
	CasADi::SXMatrix g;
	std::vector<int> gSizes;
	std::vector<std::string> gLabels;
	std::vector<double> gMin;
	std::vector<double> gMax;

	MultipleShooting & addMultipleShooting(std::string name, Ode & ode, CasADi::SX t0, CasADi::SX tf, int N);

	CasADi::SX & addParam(std::string _newParam);
	CasADi::SX & getParam(std::string p);
	double getParamSolution(std::string p);
	void boundParam(std::string p, double lb_, double ub_);
	void setParamGuess(std::string p, double guess_);
	
	CasADi::DMatrix getStateSolution(int timestep);
	CasADi::DMatrix getActionSolution(int timestep);

	void setState(  CasADi::DMatrix x, int timestep );
	void setAction( CasADi::DMatrix u, int timestep );
	void setStates(  std::vector<CasADi::DMatrix> & x);
	void setActions( std::vector<CasADi::DMatrix> & u);

	void writeOctaveOutput( std::string name );
	void writeSolution( const char * filename );
	void loadGuess( const char * filename );

	// multiple shooting instances
	std::map<std::string,MultipleShooting*> ms;

	std::vector<double> xopt;

private:

	void assertUniqueName(std::string s);
	int isMultipleShooting(std::string msName);



	// params
	std::map<std::string,int> paramsIdx;
	std::map<std::string,CasADi::SX> params;
	int np(void);
	int isParam(std::string paramName);


	CasADi::SXMatrix getParamMat(void);

	int getParamIdx(std::string p);



protected:
};
