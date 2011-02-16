#include <casadi/stl_vector_tools.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <interfaces/sundials/idas_integrator.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>
#include <interfaces/sundials/kinsol_solver.hpp>

#include <casadi/fx/fx_tools.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/matrix/matrix_tools.hpp>
#include <casadi/fx/jacobian.hpp>

#include <optimal_control/fmi_parser.hpp>
#include <optimal_control/ocp_tools.hpp>
#include <optimal_control/variable_tools.hpp>
#include <optimal_control/multiple_shooting.hpp>
#include <optimal_control/multiple_shooting_internal.hpp>

using namespace CasADi;
using namespace CasADi::Sundials;
using namespace CasADi::OptimalControl;
using namespace std;

double Tc_ref = 280;
double c_ref = 338.775766;
double T_ref = 280.099198;

double c_init = 956.271065;
double T_init = 250.051971;

int main(){

  // Allocate a parser and load the xml
  FMIParser parser("../examples/python/cstr/modelDescription.xml");

  // Dump representation to screen
//   cout << "XML representation:" << endl;
//   parser.print();

  // Obtain the symbolic representation of the OCP
  OCP ocp = parser.parse();

  // Print the ocp to screen
  ocp.print();

  // Scale the OCP
  OCP scaled_ocp = ocp.scale();
  
  // Sort the variables according to type
  OCPVariables var(ocp.variables);
  
  // Variables
  SX t = var.t.sx();
  Matrix<SX> x = sx(var.x);
  Matrix<SX> xdot = der(var.x);
  Matrix<SX> z = sx(var.z);
  Matrix<SX> p = sx(var.p);
  Matrix<SX> u = sx(var.u);

  Matrix<double> x_sca = nominal(var.x);
  Matrix<double> z_sca = nominal(var.z);
  Matrix<double> u_sca = nominal(var.u);
  Matrix<double> p_sca = nominal(var.p);

  // Initial guess for the state
  DMatrix x0(3,1,0);
  x0[0] = 0.0;
  x0[1] = T_init;
  x0[2] = c_init;
  
  // Initial guess for the control
  DMatrix u0 = 280/u_sca;
  
  // Integrator instance
  Integrator integrator;

  // Create an implicit function residual
  vector<Matrix<SX> > impres_in(ODE_NUM_IN+1);
  impres_in[0] = xdot;
  impres_in[1+ODE_T] = t;
  impres_in[1+ODE_Y] = x;
  impres_in[1+ODE_P] = u;
  SXFunction impres(impres_in,scaled_ocp.dae);
  
  // Create an implicit function (KINSOL)
  KinsolSolver ode(impres);
  ode.init();
  
  // DAE residual
  vector<Matrix<SX> > dae_in(DAE_NUM_IN);
  dae_in[DAE_T] = t;
  dae_in[DAE_Y] = x;
  dae_in[DAE_YDOT] = xdot;
  dae_in[DAE_Z] = z;
  dae_in[DAE_P] = u;
  SXFunction dae(dae_in,scaled_ocp.dae);

  bool use_kinsol = true;
  if(use_kinsol){
    // Create an ODE integrator (CVodes)
    integrator = CVodesIntegrator(ode);
    
  } else {
    // Create DAE integrator (IDAS)
    integrator = IdasIntegrator(dae);
    
  }

  // Set integrator options
  integrator.setOption("number_of_fwd_dir",1);
  integrator.setOption("number_of_adj_dir",0);
  integrator.setOption("exact_jacobian",true);
  integrator.setOption("fsens_err_con",true);
  integrator.setOption("quad_err_con",true);
  integrator.setOption("abstol",1e-8);
  integrator.setOption("reltol",1e-8);
  integrator.init();

  // Mayer objective function
  Matrix<SX> xf = symbolic("xf",x.size(),1);
  SXFunction mterm(xf, xf[0]);
  
  // Number of shooting nodes
  int num_nodes = 100;

  // Create a multiple shooting discretization
  MultipleShooting ms(integrator,mterm);
  ms.setOption("number_of_grid_points",num_nodes);
  ms.setOption("final_time",ocp.tf);
//  ms.setOption("parallelization","openmp");
  ms.init();

  for(int k=0; k<num_nodes; ++k){
    ms.input(OCP_U_INIT)(0,k) = 280/u_sca[0];
    ms.input(OCP_LBU)(0,k) = 230/u_sca[0];
    ms.input(OCP_UBU)(0,k) = 370/u_sca[0];
  }

  for(int k=0; k<=num_nodes; ++k){
    for(int i=0; i<x.size(); ++i){
      ms.input(OCP_X_INIT)(i,k) = x0[i]/x_sca[i];
    }
  }
  
  // Bounds on state
  double inf = numeric_limits<double>::infinity();
  ms.input(OCP_LBX).setAll(-inf);
  ms.input(OCP_UBX).setAll(inf);
  for(int k=1; k<=num_nodes; ++k)
    ms.input(OCP_UBX)(1,k) = 350/x_sca[1];

  // Initial condition
  ms.input(OCP_LBX)(0,0) = ms.input(OCP_UBX)(0,0) = 0/x_sca[0];
  ms.input(OCP_LBX)(1,0) = ms.input(OCP_UBX)(1,0) = 250.052/x_sca[1];
  ms.input(OCP_LBX)(2,0) = ms.input(OCP_UBX)(2,0) = 956.271/x_sca[2];
  
  IpoptSolver solver(ms.getF(),ms.getG(),FX(),ms.getJ());
  solver.setOption("tol",1e-5);
  solver.setOption("hessian_approximation", "limited-memory");
  solver.setOption("max_iter",100);
  solver.setOption("linear_solver","ma57");
  //  solver.setOption("derivative_test","first-order");
  //  solver.setOption("verbose",true);
  
  solver.init();

  // Pass to ocp solver
  ms.setNLPSolver(solver);
  
  // Solve the problem
  ms.solve();
  
  cout << solver.output() << endl;
  
  
  return 0;
}
