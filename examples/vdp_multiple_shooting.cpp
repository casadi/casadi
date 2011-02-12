#include <casadi/fx/fx_tools.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/matrix/matrix_tools.hpp>
#include <casadi/stl_vector_tools.hpp>

#include <interfaces/ipopt/ipopt_solver.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>

#include <optimal_control/multiple_shooting.hpp>

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

  //Control bounds
  double u_min[] = {-0.75};
  double u_max[] = {1.0};
  double u_init[] = {0.0};
  
  //State bounds and initial guess
  double x_min[] = {-inf,-inf,-inf};
  double x_max[] = {inf, inf, inf};
  double x_init[] = {0,0,0};

  //State bounds at the initial time
  double x0_min[] = {0,1,0};
  double x0_max[] = {0,1,0};

  //State bounds at the final time
  double xf_min[] = {0,0,-inf};
  double xf_max[] = {0,0, inf};
  
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
  MultipleShooting ms(I,mterm,ns,nx,nu);
  
  // Copy data
  ms.tf_ = tf;
  
  copy(u_min,u_min+nu,ms.u_min_.begin());
  copy(u_max,u_max+nu,ms.u_max_.begin());
  copy(u_init,u_init+nu,ms.u_init_.begin());
  
  copy(x_min,x_min+nx,ms.x_min_.begin());
  copy(x_max,x_max+nx,ms.x_max_.begin());
  copy(x_init,x_init+nx,ms.x_init_.begin());
  
  copy(x0_min,x0_min+nx,ms.x0_min_.begin());
  copy(x0_max,x0_max+nx,ms.x0_max_.begin());
  
  copy(xf_min,xf_min+nx,ms.xf_min_.begin());
  copy(xf_max,xf_max+nx,ms.xf_max_.begin());

  ms.init();
  
  IpoptSolver solver(ms.F_,ms.G_,FX(),ms.J_);
  solver.setOption("tol",1e-5);
  solver.setOption("hessian_approximation", "limited-memory");
  solver.setOption("max_iter",100);
  solver.setOption("linear_solver","ma57");
//  solver.setOption("derivative_test","first-order");

  //solver.setOption("verbose",true);
  solver.init();

  //Set bounds and initial guess
  solver.setInput(ms.V_min_,  NLP_LBX);
  solver.setInput(ms.V_max_,  NLP_UBX);
  solver.setInput(ms.V_init_,  NLP_X_INIT);
  
  solver.setInput(ms.G_min_,NLP_LBG);
  solver.setInput(ms.G_max_,NLP_UBG);

  //Solve the problem
  solver.solve();

  return 0;
}



