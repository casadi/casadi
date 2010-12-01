#include "fx_c.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "stl_string_c.hpp"
#include "stl_vector_c.hpp"
#include "../stl_vector_tools.hpp"

using namespace CasADi;
using namespace std;

// Get a reference
FX& get_fx(fx_ref ref){
  if(ref==0) throw CasadiException("get_fx failed: null pointer");
  FX *r = (FX*)(ref);
  return *r;
}

extern "C"
int casadi_ptr_to_int(void* ref, int *ref_int){
#ifdef __i386__
  *ref_int = (int)ref;
  return 0;
#else
  return 1;
#endif
}

extern "C"
int casadi_ptr_to_long(void* ref, long *ref_int){
#ifdef __x86_64__
  *ref_int = (long)ref;
  return 0;
#else
  return 1;
#endif
}

extern "C"
fx_ref casadi_fx_new(){
  return new FX();
}


extern "C"
int casadi_fx_delete(fx_ref ref){
  FX* ptr = &get_fx(ref);
  delete ptr;
  return 0;
}

extern "C"
int casadi_fx_print_cout(fx_ref ref){
  try{
    cout << get_fx(ref) << endl;
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_print_cerr(fx_ref ref){
  try{
    cout << get_fx(ref) << endl;
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_init(fx_ref ref){
  try{
    get_fx(ref).init();
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_setoption_string(fx_ref ref, const char* opname, const char* opval){
  try{
    get_fx(ref).setOption(opname,opval);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_getoption_double(fx_ref ref, const char* name, vector_ptr str){
  try{
    get_stl_vector(ref) = get_fx(ref).getOption(name).toDoubleVector();
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_getoption_string(fx_ref ref, const char* name, string_ptr str){
  try{
    get_stl_string(str) = get_fx(ref).getOption(name).toString();
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_option_is_string(fx_ref ref, const char* name, int *is_string){
  try{
    *is_string = get_fx(ref).getOption(name)->is_string;
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_print_options(fx_ref ref){
  try{
    get_fx(ref).printOptions();
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_setoption_double(fx_ref ref, const char* opname, const double* opval, int n){
  try{
    vector<double> temp(opval,opval+n);
    get_fx(ref).setOption(opname,temp);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_input_size(fx_ref ref, int ind, int *sz){
  try{
    *sz = get_fx(ref).getInputSize(ind);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_getinput(fx_ref ref, int ind, int ord, double* val){
  try{
    get_fx(ref).getInput(val,ind,ord);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_setinput(fx_ref ref, int ind, int ord, const double* val){
 try{
    get_fx(ref).setInput(val,ind,ord);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_output_size(fx_ref ref, int ind, int *sz){
  try{
    *sz = get_fx(ref).getOutputSize(ind);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_getoutput(fx_ref ref, int ind, int ord, double* val){
  try{
    get_fx(ref).getOutput(val,ind,ord);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_setoutput(fx_ref ref, int ind, int ord, const double* val){
  try{    
    get_fx(ref).setOutput(val,ind,ord);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_evaluate(fx_ref ref){
  try{
    get_fx(ref).evaluate();
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

// extern "C"
// int casadi_fx_evaluate_fwd(fx_ref ref){
//   try{
//       get_fx(ref).evaluate(1,0);
//     return 0;
//   } catch(exception &e){
//     cerr << e.what();
//     return 1;
//   } catch(...){
//     return 1;
//   }  
// }

// extern "C"
// int casadi_fx_evaluate_adj(fx_ref ref){
//   try{
//     get_fx(ref).evaluateAdj();
//     return 0;
//   } catch(exception &e){
//     cerr << e.what();
//     return 1;
//   } catch(...){
//     return 1;
//   }
// }

extern "C"
int casadi_fx_print_string(fx_ref ref,string_ptr str){
  try{
    // Print to a temporary variable
    stringstream ss;
    ss << get_fx(ref);
    get_stl_string(str) = ss.str();
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}


extern "C"
int casadi_fx_jacobian(fx_ref ref, fx_ref res, int iind, int oind){
  try{
    get_fx(res) = get_fx(ref).jacobian(iind,oind);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_fx_hessian(fx_ref ref, fx_ref res, int iind, int oind){
  try{
    get_fx(res) = get_fx(ref).hessian(iind,oind);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}
