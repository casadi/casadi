#include "mx_function_c.h"
#include "../fx/mx_function.hpp"
#include "mx_c.hpp"
#include "fx_c.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>

using namespace CasADi;
using namespace std;


extern "C"
int casadi_mx_function(fx_ref fcn, mx_vec iv, mx_vec ov){
 try{
    get_fx(fcn) = MXFunction(get_mx_vec(iv),get_mx_vec(ov));
    return 0;
 } catch(exception &e){
    cerr << e.what();
    return 1;
 } catch(...){
    return 1;
 }
}
