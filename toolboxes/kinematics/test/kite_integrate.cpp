/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
#include <fstream>
#include <string>

#include <kinetics.hpp>
#include <casadi/stl_vector_tools.hpp>


#include "kiteode.hpp"


using namespace std;
//using namespace KinVec;




int main( ){ 

   
    int nStates=17;
    Matrix q("q",17,1);
    
    SX t=0.0;

    Matrix p = 0.0;
    Matrix d = 0.0;
    Matrix u = 0.0;
    Matrix h = zeros(9,1);
    Matrix F=ode(t,q,p,d,u,h);
        


    
    Matrix arg;arg << q;
  
    SXFunction temp(arg,F);

	double x0test[17]={0,1,0.1,0.2,0.3,0.4,0.5,10,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,0,0.1,0.05};
	std::vector<double> res(17);
	temp.setArgument(x0test);
	temp.evaluate();
	temp.getValue(&res[0]);

    temp.generateCode("ode.c");
    

  
  
   cout << "Num eval "<< res << endl;


   cout << "Take a nap, taking symbolic derivative" << endl;

    F.jacobian(q);
}


