#include <casadi/stl_vector_tools.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <interfaces/sundials/idas_integrator.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>
#include <interfaces/sundials/kinsol_solver.hpp>
#include <interfaces/csparse/csparse.hpp>

#include <casadi/fx/fx_tools.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/matrix/matrix_tools.hpp>
#include <casadi/fx/jacobian.hpp>

#include <optimal_control/fmi_parser.hpp>
#include <optimal_control/ocp_tools.hpp>
#include <optimal_control/variable_tools.hpp>
#include <optimal_control/multiple_shooting.hpp>

using namespace CasADi;
using namespace CasADi::Interfaces;
using namespace CasADi::Sundials;
using namespace CasADi::OptimalControl;
using namespace std;

int main(){

  // Allocate a parser and load the xml
  FMIParser parser("../examples/python/cstr/modelDescription.xml");

  // Obtain the symbolic representation of the OCP
  OCP ocp = parser.parse();

  // Print the ocp to screen
  ocp.print();

  // Scale the OCP
  OCP scaled_ocp = ocp.scale();
    
  // Correct the inital guess and bounds on variables
  ocp.u_[0].setStart(280);
  ocp.u_[0].setMin(230);
  ocp.u_[0].setMax(370);
  
  // Correct bound on state
  ocp.x_[1].setMax(350);
  
  // Variables
  SX t = ocp.t_.sx();
  Matrix<SX> x = sx(ocp.x_);
  Matrix<SX> xdot = der(ocp.x_);
  Matrix<SX> z = sx(ocp.z_);
  Matrix<SX> p = sx(ocp.p_);
  Matrix<SX> u = sx(ocp.u_);

  // Initial guess and bounds for the state
  vector<double> x0 = getStart(ocp.x_,true);
  vector<double> xmin = getMin(ocp.x_,true);
  vector<double> xmax = getMax(ocp.x_,true);
  
  // Initial guess and bounds for the control
  vector<double> u0 = getStart(ocp.u_,true);
  vector<double> umin = getMin(ocp.u_,true);
  vector<double> umax = getMax(ocp.u_,true);
  
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
  ode.setLinearSolver(CSparse(CRSSparsity()));
  ode.setOption("linear_solver","user_defined");
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

  // Initial condition
  for(int i=0; i<x.size(); ++i){
    ms.input(OCP_X_INIT)(i,0) = ms.input(OCP_LBX)(i,0) = ms.input(OCP_UBX)(i,0) = x0[i];
  }

  // State bounds
  for(int k=1; k<=num_nodes; ++k){
    for(int i=0; i<x.size(); ++i){
      ms.input(OCP_X_INIT)(i,k) = x0[i];
      ms.input(OCP_LBX)(i,k) = xmin[i];
      ms.input(OCP_UBX)(i,k) = xmax[i];
    }
  }

  // Control bounds
  for(int k=0; k<num_nodes; ++k){
    for(int i=0; i<u.size(); ++i){
      ms.input(OCP_U_INIT)(i,k) = u0[i];
      ms.input(OCP_LBU)(i,k) = umin[i];
      ms.input(OCP_UBU)(i,k) = umax[i];
    }
  }

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
