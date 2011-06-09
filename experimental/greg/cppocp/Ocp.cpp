// Ocp.cpp
// Greg Horn
// Casadi 2011

#include "Ocp.hpp"

#include <casadi/fx/sx_function.hpp>

using namespace std;
using namespace CasADi;

Ocp::Ocp(Ode * _ode)
{
	ode = _ode;
}

Ocp::~Ocp(){}

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


//    def addNonlcon(self, lhs, rel, rhs):
//        if rel == "==":
//            self._addNonlconEq(lhs - rhs)
//        elif rel == "<=":
//            self._addNonlconIneq(lhs - rhs)
//        elif rel == ">=":
//            self._addNonlconIneq(rhs - lhs)
//        else:
//            errStr = "invalid relation \""+rel+"\""
//            raise ValueError(errStr)

//class MultiStageOcp():
//    def __init__(self, ssOcps):
//        if not isinstance(ssOcps, list):
//            ssOcps = [ssOcps]
//
//        self.ssOcps = {}
//        for ssOcp in ssOcps:
//            self.ssOcps[ssOcp.ode.name] = ssOcps
//
//        self.bigBigN = sum([ssOcp._getBigN() for ssOcp in self.ssOcps.values()])
//        self.designVariables = C.symbolic('designVariables', self.bigBigN)
//        #self.designVariables = C.MX('designVariables', self.bigBigN)
//
//
