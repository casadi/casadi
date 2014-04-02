#include <iostream>

#include <symbolic/std_vector_tools.hpp>
#include <symbolic/sx/sx_tools.hpp>
#include <symbolic/function/sx_function.hpp>
//#include <symbolic/function/jacobian.hpp>

#include <SnoptInterface.hpp>

#include <string>
#include <map>

using namespace CasADi;
using namespace std;

int
main()
{
	//Matrix<SX>x = ssym("x", 2);
	vector<SX>x = ssym("x", 2).data();
	//vector<SX>x = ssym("x",2,1,1);

	SX f;
	vector<SX>g(2);

//   F[0] = x[1];
//   F[1] = x[0]*x[0] + 4*x[1]*x[1];
//   F[2] = (x[0] - 2)*(x[0] - 2) + x[1]*x[1];

	f = x[1];

	g[0] = x[0]*x[0] + 4*x[1]*x[1];
	g[1] = (x[0] - 2)*(x[0] - 2) + x[1]*x[1];


	vector<SX>fg(3);
	//vector<SX>fg(3);
	fg[0] = f;
	fg[1] = g[0];
	fg[2] = g[1];

	SXFunction fgfcn(x, fg); // objective function

	//	fgfcn.setOption("ad_mode","reverse");
	//	fgfcn.setOption("symbolic_jacobian",false);

	SnoptInterface si(fgfcn);

	si.xlow[0] = 0;
	si.Fupp[1] = 4.0;
	si.Fupp[2] = 5.0;

	si.x[0] = 1;
	si.x[1] = 1;

	si.run();
}
