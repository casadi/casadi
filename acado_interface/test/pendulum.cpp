#include <acado_interface/acado_integrator.hpp>
#include <casadi/stl_vector_tools.hpp>

using namespace CasADi;
using namespace std;

int main( ){
    // Variables
    SX t("t");       // time
    SX phi("phi");   // the angle phi
    SX dphi("dphi"); // the first derivative of phi w.r.t time
    SX F("F");       // a force acting on the pendulum
    SX l("l");       // the length of the pendulum

    // Constants
    double m     = 1.0  ;    // the mass of the pendulum
    double g     = 9.81 ;    // the gravitational constant
    double alpha = 2.0  ;    // frictional constant

    // Differential equation
    SX der_phi = dphi;
    SX der_dphi = -(m*g/l)*sin(phi) - alpha*dphi + F/(m*l);

    // state vector
    vector<SX> x;
    x.push_back(phi);
    x.push_back(dphi);
    
    // Right hand side of the ode
    vector<SX> xdot;
    xdot.push_back(der_phi);
    xdot.push_back(der_dphi);
  
    // parameters
    vector<SX> p;
    p.push_back(F);
    p.push_back(l);

    // the arguments of the ode
    vector< vector<SX> > ffcn_arg(ODE_NUM_IN);
    ffcn_arg[ODE_T] = vector<SX>(1,t);
    ffcn_arg[ODE_Y] = x;
    ffcn_arg[ODE_P] = p;
  
    // create an ode rhs function
    SXFunction ffcn(ffcn_arg,xdot);
    
    // Create an integrator instance
    ACADOIntegrator integrator(ffcn);
    integrator.setOption("print_level","medium");
    integrator.setOption("ad_order",2);
  
    // Initialize
    integrator.init();

    // Set time horizon
    integrator.input(INTEGRATOR_T0).set(0.0);
    integrator.input(INTEGRATOR_TF).set(1.0);

    // Set initial values
    vector<double> x_start(2, 0.0);
    integrator.input(INTEGRATOR_X0).set(x_start);
    
    // Set parameters and controls
    vector<double> up0(2, 1.0);
    integrator.input(INTEGRATOR_P).set(up0);
    
    // Set forward seeds
    vector<double> x0_seed1(2,0.0);
    x0_seed1[0] = 1;
    integrator.input(INTEGRATOR_X0).setF(x0_seed1);

    vector<double> x0_seed2(2,0.0);
    x0_seed2[0] = 1;
    integrator.input(INTEGRATOR_X0).setF(x0_seed2); // TODO: second order!!!

    vector<double> p_seed1(2,0.0);
    integrator.input(INTEGRATOR_P).setF(p_seed1);

    vector<double> p_seed2(2,0.0);
    integrator.input(INTEGRATOR_P).setF(p_seed2); // TODO: second order!!!
    
    // Solve
    integrator.evaluate(2,0);
    
    // Display
    cout << "x = " << integrator.output(INTEGRATOR_XF).data() << endl;
    cout << "Dx = " << integrator.output(INTEGRATOR_XF).dataF() << endl;
//    cout << "DDx = " << integrator.output(INTEGRATOR_XF,2) << endl;

    
    return 0;
}



