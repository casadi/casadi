#include <casadi/fx/fx_tools.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/matrix/matrix_tools.hpp>
#include <casadi/stl_vector_tools.hpp>

#include <interfaces/ipopt/ipopt_solver.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>

#include <optimal_control/multiple_shooting.hpp>
#include <optimal_control/multiple_shooting_internal.hpp>

using namespace CasADi;
using namespace CasADi::Sundials;
using namespace CasADi::OptimalControl;
using namespace std;


int main(){
  //Final time (fixed)
  double tf = 20.0;

  // Infinity
  double inf = numeric_limits<double>::infinity();

  // Declare variables (use simple, efficient DAG)
  SX t("t"); //time
  SX x("x"), y("y"), u("u"), L("cost");
  
  // All states
  vector<SX> xx(3);  xx[0] = x;  xx[1] = y;  xx[2] = L;

  //ODE right hand side
  vector<SX> f(3);
  f[0] = (1 - y*y)*x - y + u;
  f[1] = x;
  f[2] = x*x + y*y + u*u;
  
  //ODE right hand side function
  vector<SXMatrix> rhs_in(ODE_NUM_IN);
  rhs_in[ODE_T] = t;
  rhs_in[ODE_Y] = xx;
  rhs_in[ODE_P] = u;
  SXFunction rhs(rhs_in,f);

  //Create an integrator (CVodes)
  CVodesIntegrator I(rhs);
  I.setOption("abstol",1e-8); //abs. tolerance
  I.setOption("reltol",1e-8); //rel. tolerance
  I.setOption("steps_per_checkpoint",500);
  I.setOption("stop_at_end",true);
  I.init();
  
  //Numboer of shooting nodes
  int ns = 50;

  // Number of differential states
  int nx = 3;
  
  // Number of controls
  int nu = 1;
  
  // Mayer objective function
  Matrix<SX> xf = symbolic("xf",nx,1);
  SXFunction mterm(xf, xf[nx-1]);

  // Create a multiple shooting discretization
  MultipleShooting ms(I,mterm);
  ms.setOption("number_of_grid_points",ns);
  ms.setOption("final_time",tf);
  
  ms.init();

  //Control bounds
  double u_min[] = {-0.75};
  double u_max[] = {1.0};
  double u_init[] = {0.0};
  
  for(int k=0; k<ns; ++k){
    copy(u_min,u_min+nu,ms.input(OCP_LBU).begin()+k*nu);
    copy(u_max,u_max+nu,ms.input(OCP_UBU).begin()+k*nu);
    copy(u_init,u_init+nu,ms.input(OCP_U_INIT).begin()+k*nu);
  }
  
  fill(ms.input(OCP_LBX).begin(),ms.input(OCP_LBX).end(),-inf);
  fill(ms.input(OCP_UBX).begin(),ms.input(OCP_UBX).end(),inf);
  fill(ms.input(OCP_X_INIT).begin(),ms.input(OCP_X_INIT).end(),0);

  // Initial condition
  ms.input(OCP_LBX)[0] = ms.input(OCP_UBX)[0] = 0;
  ms.input(OCP_LBX)[1] = ms.input(OCP_UBX)[1] = 1;
  ms.input(OCP_LBX)[2] = ms.input(OCP_UBX)[2] = 0;

  // Final condition
  ms.input(OCP_LBX)[ns*nx+0] = ms.input(OCP_UBX)[ns*nx+0] = 0;
  ms.input(OCP_LBX)[ns*nx+1] = ms.input(OCP_UBX)[ns*nx+1] = 0;

  IpoptSolver solver(ms.getF(),ms.getG(),FX(),ms.getJ());
  solver.setOption("tol",1e-5);
  solver.setOption("hessian_approximation", "limited-memory");
  solver.setOption("max_iter",100);
  solver.setOption("linear_solver","ma57");
  //  solver.setOption("derivative_test","first-order");

  //solver.setOption("verbose",true);
  solver.init();

  // Pass NLP solver
  ms.setNLPSolver(solver);
  
  // Solve the problem
  ms.evaluate(0,0);

  return 0;
}



