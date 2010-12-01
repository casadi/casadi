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


