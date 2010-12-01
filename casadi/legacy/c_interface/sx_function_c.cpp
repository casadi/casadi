#include "sx_function_c.h"
#include "../fx/sx_function.hpp"
#include "sx_matrix_c.hpp"
#include "fx_c.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>

using namespace CasADi;
using namespace std;


extern "C"
int casadi_sx_function(fx_ref fcn, sx_matrix_vec iv, sx_matrix_vec ov){
 try{
    get_fx(fcn) = SXFunction(get_sx_matrix_vec(iv),get_sx_matrix_vec(ov));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_function_evaluate_symbolically(fx_ref fcn, sx_matrix_vec iv, sx_matrix_ref res){
 try{
    SXFunction fun(get_fx(fcn));
    get_sx_matrix(res) = fun(get_sx_matrix_vec(iv));
    return 0;
  } catch(...){
    return 1;
  }
}
