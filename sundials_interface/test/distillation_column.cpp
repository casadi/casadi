#include <integrator/sundials_interface/cvodes_integrator.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <casadi/expression_tools.hpp>

#include <fstream>

using namespace OPTICON;
using namespace std;

static double u0[15] = {1, 1, 1, 2, 2, 3, 3, 5, 6, 7, 9, 10, 12, 15, 15};

int main(){
  try{

  // Parameters
  double V = 100;
  double m = 0.1;
  double mc = 0.1;
  int N = 5;

  // Create a single-stage optimal control problem
  OCP_old ocp;

  // Time variable
  Matrix t = ocp.getTime();

  // Specify the time horizon
  Matrix t0 = 0,  tf = 1.4;
  ocp.setInitialTime(t0);
  ocp.setFinalTime(tf);

  // Differential states
  Matrix M0("M0");
  Matrix x("x",N+2);
  Matrix MD("MD");
  Matrix xD("xD");
  Matrix alpha("alpha");
  
  // Augmented state vector
  Matrix X;
  X << M0 << x << MD << xD << alpha;
 
  // Derivative of X
  Matrix XDOT("XDOT",X.size());
  ocp.addState(X,XDOT);

  // Control
  Matrix R("R");
  ocp.addControl(R);

  Matrix R0(&u0[0], 15,1);
  Matrix T0 = linspace(0,1,15);
  R0 = pw_lin(t, T0, R0);

  ocp.guessSolution(R,R0);


//  ocp.guessSolution(R,15*t*t);

  // Differential equation
  Matrix L = R/(1+R)*V;

  Matrix y(N+1,1);
  for(int i=0; i<N+1; ++i){
    y[i] = x[i]*(alpha+1)/(x[i]+alpha);
  }

  Matrix G;
  G << -V+L;
  G << (L*x[1]-V*y[0] + (V-L)*x[0])/M0; 
 for(int i=1; i<N+1; ++i)
    G << (L*x[i+1]-V*y[i] + V*y[i-1] - L*x[i])/m;
  G << V/mc*(y[N] - x[N+1]);
  G << V-L;
  G << (V-L)/MD*(x[N+1]-xD);
  G << 0;

  ocp.addEquation(XDOT, G);

  // Augment with sensitivity equations
  Matrix D1("D1",X.size()), D2("D2",X.size());
  Matrix D1DOT("D1DOT",D1.size()), D2DOT("D2DOT",D2.size());

  ocp.addState(D1,D1DOT);
  ocp.addState(D2,D2DOT);

  Matrix JG = jacobian(G,X);
  ocp.addEquation(D1DOT, JG*D1);
  ocp.addEquation(D2DOT, JG*D2);


//   ocp.addEquation(M0.der(t), -V+L);
//   ocp.addEquation(x[0].der(t), (L*x[1]-V*y[0] + (V-L)*x[0])/M0); 
//   for(int i=1; i<N+1; ++i)
//     ocp.addEquation(x[i].der(t), (L*x[i+1]-V*y[i] + V*y[i-1] - L*x[i])/m);
//   ocp.addEquation(x[N+1].der(t), V/mc*(y[N] - x[N+1]));
//   ocp.addEquation(MD.der(t), V-L);
//   ocp.addEquation(xD.der(t), (V-L)/MD*(x[N+1]-xD));
//   ocp.addEquation(alpha.der(t), 0);

  // Initial conditions
  ocp.addInitialCondition(M0, 100);
  ocp.addInitialCondition(x[0], 0.5);
  for(int i=1; i<N+2; ++i)
    ocp.addInitialCondition(x[i], 1);
  ocp.addInitialCondition(MD, 0.1);
  ocp.addInitialCondition(xD, 1);
  ocp.addInitialCondition(alpha, 0.2);

  Matrix D1_0(D1.size());
  D1_0[1] = 1;
  ocp.addInitialCondition(D1, D1_0);

  Matrix D2_0(D2.size());
  D2_0[D2.size()-1] = 1;
  ocp.addInitialCondition(D2, D2_0);

  // LAGRANGE TERM
 //  ocp.addOutput(u*u, 10, true); // just quadrature
 //  ocp.addOutput(u*u, 10, true, 1, u); // quadrature and sensitivities, first order

  // forward sensitivities of v with respect to u at time tf
//  ocp.addOutput(v, 10, false, 1, u);

  // The desired output times
  int nt_out = 100; // number of outputs
  Matrix t_out = linspace(0,tf,nt_out);
  
  // Set output function
  int oind_M0 = ocp.addOutput(M0, t_out);
  int oind_x0 = ocp.addOutput(x[0], t_out);
  int oind_MD = ocp.addOutput(MD, t_out);
  int oind_xD = ocp.addOutput(xD, t_out);


  int oind_D1 = ocp.addOutput(D1[1], t_out);
  int oind_D2 = ocp.addOutput(D2[D2.size()-1], t_out);

  // Forward sensitivities
  Matrix Y;
  Y << X << D1 << D2;
  ocp.addOutput(Y, t_out, false, 1, R);


  // Print to screen
  std::cout << ocp;

  // REFORMULATE THE PROBLEM

  // Discretize controls into 20 uniform intervals
  Matrix R_disc = parametrizeControls(ocp,R,20);
  std::cout << R_disc << std::endl;

  // Make the DE explicit
//   ocp.makeExplicit();

  // Reformulate switches as integer controls and zero-crossing functions
  makeSmooth(ocp);
 
 // Allocate an integrator
  CvodesIntegrator integrator(ocp);

  // Set the linear solver
   integrator->setOption("linear_solver","dense");
//  integrator->setOption("linear_solver","band");
//  integrator->setOption("linear_solver","sparse");

  // Use exact jacobian
  integrator->setOption("exact_jacobian",true);

  // Upper and lower band-widths (only relevant when using a band linear solver)
//   integrator->setOption("mupper",1);
//   integrator->setOption("mlower",1);

  // set tolerances 
  integrator->setOption("reltol",1e-6);
  integrator->setOption("abstol",1e-8);

  // Initialize the integrator
  integrator->init();

  // Integrate once for visualization
  integrator->evaluateFwd();

  // Print stats
  integrator->printStats();

  // Create a file for saving the results
  std::ofstream resfile;
  resfile.open ("results_dest.txt");

  // Save results to file
  resfile << "t_out " << t_out << std::endl;
  resfile << "M0 " << integrator->output[oind_M0].data() << endl;
  resfile << "x0 " << integrator->output[oind_x0].data() << endl;
  resfile << "MD " << integrator->output[oind_MD].data() << endl;
  resfile << "xD " << integrator->output[oind_xD].data() << endl;

  resfile << "D1 " << integrator->output[oind_D1].data() << endl;
  resfile << "D2 " << integrator->output[oind_D2].data() << endl;
//  resfile << "R0 " << R0.evaluate(t,t_out) << std::endl;
  
  // Close the results file
  resfile.close();


  return 0;
} catch (const char * str){
  std::cerr << str << std::endl;
  return 1;
}

}
