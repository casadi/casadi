#include "sundials_interface_c.h"
#include "cvodes_integrator.hpp"
#include "idas_integrator.hpp"
#include "casadi/c_interface/fx_c.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>

using namespace CasADi::Sundials;
using namespace std;

extern "C"
int casadi_cvodes_integrator(fx_ref fcn, fx_ref ffcn){
 try{
    get_fx(fcn) = CVodesIntegrator(get_fx(ffcn));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_idas_integrator(fx_ref fcn, fx_ref ffcn){
 try{
    get_fx(fcn) = IdasIntegrator(get_fx(ffcn));
    return 0;
  } catch(...){
    return 1;
  }
}

