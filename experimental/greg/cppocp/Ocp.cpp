// Ocp.cpp
// Greg Horn
// Casadi 2011

#include "Ocp.hpp"
#include <cstdlib>

#include <casadi/fx/sx_function.hpp>

using namespace std;
using namespace CasADi;

Ocp::Ocp(){}

Ocp::~Ocp()
{
	map<string,MultipleShooting*>::const_iterator msIter;
	for (msIter = ms.begin(); msIter != ms.end(); msIter++){
//		cout << "destroying " << ms[msIter->first] << endl;
		delete msIter->second;
	}
}

void Ocp::addNonlconIneq( SXMatrix gNew )
{
	if (gNew.size2() != 1){
		cerr << "gNew.size2() != 1" << endl;
		throw 1;
	}

	g = vertcat(g, gNew);

	for (int k=0; k<gNew.size1(); k++){
		gMin.push_back(-1.0e50);
		gMax.push_back(0.0);
		//	-numeric_limits<double>::infinity()
	}
}

void Ocp::addNonlconEq( SXMatrix gNew )
{
	if (gNew.size2() != 1){
		cerr << "gNew.size2() != 1" << endl;
		throw 1;
	}

	g = vertcat(g, gNew);

	for (int k=0; k<gNew.size1(); k++){
		gMin.push_back(0);
		gMax.push_back(0);
	}
}

/****************** multiple shooting stuff ****************/
MultipleShooting & Ocp::addMultipleShooting(string name, Ode & ode, SX t0, SX tf, int N)
{
	cout << "Ocp::addMultipleShooting needs assertUniqueName()\n";

	int numNew = ode.nxu()*N;
	designVariables = vertcat( designVariables, SXMatrix(create_symbolic(name, numNew)) );
	for (int k=0; k<numNew; k++){
		lb.push_back(-1e50);
		ub.push_back(-1e50);
		guess.push_back(0);
	}

	ms[name] = new MultipleShooting(name, ode, t0, tf, N, lb, ub, guess, designVariables, designVariables.size1() - ode.nxu()*N, params);

	// dynamics constraint
	for (int k=0; k<N-1; k++)
		addNonlconEq( ms[name]->getDynamicsConstraintError(k, params) );

	return *(ms[name]);
}


/****************** param stuff ******************/
// get param
SX & Ocp::getParam(string p)
{
	return params[p];
}

// add a param
SX & Ocp::addParam(string _newParam)
{
	cout << "Ocp::addParam needs assertUniqueName()\n";
	//	assertUniqueName(_newParam);

//	designVariables = vertcat( designVariables, SXMatrix(create_symbolic(_newParam, 1)) );
	int idx = designVariables.size1();

	SXMatrix newDv = create_symbolic(_newParam, 1);
	designVariables = vertcat( designVariables, newDv );
	guess.push_back(0);
	lb.push_back(-1e50);
	ub.push_back(1e50);

	paramsIdx[_newParam] = idx;
	params[_newParam] = designVariables.at(idx);

	return getParam(_newParam);

}

void Ocp::boundParam(string p, double lb_, double ub_)
{
	int idx = paramsIdx[p];
	lb[idx] = lb_;
	ub[idx] = ub_;
	if (guess[idx] < lb[idx]) guess[idx] = lb[idx];
	if (guess[idx] > ub[idx]) guess[idx] = ub[idx];
}



void Ocp::setParamGuess(string p, double guess_)
{
	int idx = paramsIdx[p];
	guess[idx] = guess_;
	if (guess[idx] < lb[idx] ){
		cerr << "Requested initial guess " << p << " == " << guess[idx];
		cerr << " is less than lb[" << idx << "] == " << lb[idx] << endl;
		cerr << "Setting guess " << p << " = " << lb[idx] << endl;
		guess[idx] = lb[idx];
	}
	if (guess[idx] > ub[idx] ){
		cerr << "Requested initial guess " << p << " == " << guess[idx];
		cerr << " is greater than ub[" << idx << "] == " << ub[idx] << endl;
		cerr << "Setting guess " << p << " = " << ub[idx] << endl;
		guess[idx] = ub[idx];
	}
}



void Ocp::writeMatlabOutput( const char * filename, double * xOpt)
{
	char filename2[200];
	sprintf(filename2, "%s.m", filename);
	ofstream f(filename2);

	if (!f){
		cerr << "error opening " << filename2 << endl;
		exit(1);
	}
	f.precision(10);

	f << "function p = " << filename << "()" << endl;

	map<string, int>::const_iterator iter;

	// params
	for (iter = paramsIdx.begin(); iter != paramsIdx.end(); iter++){
		int idx = paramsIdx[ iter->first ];
		f << "p." << iter->first << " = " << xOpt[idx] << ";" << endl;
	}
	f << endl;

	f.close();
}
