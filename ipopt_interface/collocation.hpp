#ifndef COLLOCATION_HPP
#define COLLOCATION_HPP

#include <vector>
#include "ipopt_interface.hpp"
#include <casadi/fx/sx_function.hpp>
#include <acado_modelling/ocp_tools.hpp>

namespace IpoptInterface{
    using namespace OPTICON;


class Collocation{
  public:

    // default constructor
    Collocation();

    void optimize();
  
    // Degree of interpolating polynomial
    int K;
  
    // Number of finite elements
    int N;

    // Derivative of the state vector
    SXMatrix xdot;
    SXFunction ffcn;

    SXMatrix c;
    SXFunction cfcn;

    SXMatrix m;
    SXFunction mfcn;

        CollocationPoints cp;               

    // Final time
    double tf;
    
    SXMatrix x;
    vector<double> x_lb, x_ub, x_init;

    SXMatrix u;
    vector<double> u_lb, u_ub, u_init;
                
    SX t;
    
  // Collocated times
    vector<vector<double> > T;

    vector<double> x_opt;
};

} // namespace IpoptInterface

#endif // COLLOCATION_HPP
