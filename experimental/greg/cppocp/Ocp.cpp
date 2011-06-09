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

void Ocp::addNonlconIneq( vector<SX> gNew )
{
	// loop over new constraint g and append to internal constraints G
	vector<SX>::const_iterator iter;
	for (iter = gNew.begin(); iter != gNew.end(); ++iter){
		g.push_back(*iter);
		gMin.push_back(-1.0e50);
		gMax.push_back(0.0);
		//	-numeric_limits<double>::infinity()
	}
}

//void Ocp::addNonlconEq( vector<SX> gNew )
//{
//	// loop over new constraint g and append to internal constraints G
//	vector<SX>::const_iterator iter;
//	for (iter = gNew.begin(); iter != gNew.end(); ++iter){
//		G.push_back(*iter);
//		Gmin.push_back(0);
//		Gmax.push_back(0);
//	}
//
//}
//

void Ocp::addNonlconEq( SXMatrix gNew )
{
	if (gNew.size2() != 1){
		cerr << "gNew.size2() != 1" << endl;
		throw 1;
	}

	// loop over new constraint g and append to internal constraints G
	SXMatrix::const_iterator iter;
	for (iter = gNew.begin(); iter != gNew.end(); ++iter){
		g.push_back(*iter);
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
