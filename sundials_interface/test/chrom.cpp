#include <integrator/sundials_interface/cvodes_integrator.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <casadi/expression_tools.hpp>

#include <fstream>

using namespace std;
using namespace OPTICON;

int main(){
  try{

  // Create a single-stage optimal control problem
  OCP_old ocp;

  // Time variable
  Matrix t = ocp.getTime();

  // Specify the time horizon
  Matrix t0 = 0,  tf = 100;
  ocp.setInitialTime(t0);
  ocp.setFinalTime(tf);

  // Discretization
  int n = 200;

  // Constants
  double c_in = 1.43e-4;   // Concentration at inlet
  double L = 1.3;
  double Dc = 1e-3;
  double lambda = 1.2;
//   double nu[] = {0,4.7}
  double u = 0.03; 

  double dz = L/n;
  double k_a = 0.01, k_d = 0.01, q_m = 10;

  // Differential states
  Matrix c("c",n), cdot("cdot",n);
  Matrix q("q",n), qdot("qdot",n);
  ocp.addState(c,cdot);
  ocp.addState(q,qdot);

  // Control
/*  Matrix u("u");
  ocp.addControl(u);*/
//   ocp.guessSolution(u,3*sin(t));


  // Differential equation
  Matrix fc(n,1), fq(n,1); // right-hand side of the ode

  fc[0] = -u*(c[0] - c_in)/dz + Dc*(-c[0]+c[1])/dz/dz;
  for(int i=1; i<n-1; ++i)
    fc[i] = -u*(c[i] - c[i-1])/dz + Dc*(c[i-1]-2*c[i]+c[i+1])/dz/dz;
  fc[n-1] = -u*(c[n-1] - c[n-2])/dz + Dc*(-c[n-1]+c[n-2])/dz/dz;


  for(int i=0; i<n; ++i)
    fq[i] = k_a*c[i] * (q_m - q[i]) - k_d*q[i];
   
  ocp.addEquation(cdot, fc);
  ocp.addEquation(qdot, fq);

  // Initial conditions
  Matrix c_init = Matrix(c.size());
  ocp.addInitialCondition(c, c_init);

  Matrix q_init = Matrix(q.size());
  ocp.addInitialCondition(q, q_init);

  // The desired output times
  int nt_out = 100; // number of outputs
  Matrix t_out = linspace(t0,tf,nt_out);

  // Set output function
  int oind_c = ocp.addOutput(t_out, c);
  int oind_q = ocp.addOutput(t_out, q);

  // Print to screen
  cout << ocp;

 // Allocate an integrator
  CvodesIntegrator integrator(ocp);

  // Set the linear solver
//   integrator->setOption("linear_solver","dense");
  integrator->setOption("linear_solver","band");
//  integrator->setOption("linear_solver","sparse");

  // Use exact jacobian
  integrator->setOption("exact_jacobian",false);

  // Upper and lower band-widths (only relevant when using a band linear solver)
  integrator->setOption("mupper",1);
  integrator->setOption("mlower",1);

  // set tolerances 
  integrator->setOption("reltol",1e-6);
  integrator->setOption("abstol",1e-10);

  // Initialize the integrator
  integrator->init();
  cout << "initialized" << endl;

  // Integrate once for visualization
  integrator->evaluate();
  cout << "solved" << endl;

  // Print integrator statistics
  integrator->printStats();

  // Create a file for saving the results
  ofstream resfile;
  resfile.open ("results_chrom.txt");

  // Save results to file
  resfile << "t_out " << t_out << endl;
  resfile << "c " << integrator->output[oind_c].data() << endl;
  resfile << "q " << integrator->output[oind_q].data() << endl;

  // Close the results file
  resfile.close();


  return 0;
} catch (const char * str){
  cerr << str << endl;
  return 1;
}

}
