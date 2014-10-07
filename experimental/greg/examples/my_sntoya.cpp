/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <iostream>

#include <core/std_vector_tools.hpp>
#include <core/sx/sx_tools.hpp>
#include <core/function/sx_function.hpp>
//#include <core/function/jacobian.hpp>

#include <SnoptInterface.hpp>

#include <string>
#include <map>

using namespace casadi;
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
