#include "ipopt_solver_c.h"
#include "ipopt_solver.hpp"
#include "casadi/c_interface/fx_c.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>

using namespace CasADi;
using namespace std;


extern "C"
int casadi_ipopt_solver(fx_ref fcn, fx_ref f, fx_ref g, fx_ref h, fx_ref j){
 try{
    get_fx(fcn) = IpoptSolver(get_fx(f),get_fx(g),get_fx(h),get_fx(j));
    return 0;
  } catch(...){
    return 1;
  }
}

