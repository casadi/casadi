#include "stl_string_c.hpp"

using namespace std;

std::string& get_stl_string(string_ptr str){
  if(str==0) throw "getmxvec failed: null pointer";
  string *r = (string*)(str);
  return *r;  
}

extern "C"
string_ptr casadi_string_new(void){
  try{
    return new string();
  } catch(...){
    return 0;
  }
}

extern "C"
int casadi_string_delete(string_ptr str){
  try{
    string *ptr = &get_stl_string(str);
    delete ptr;
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_string_assign(string_ptr str, const char* s){
  try{
    get_stl_string(str) = string(s);
  } catch(...){
    return 1;
  }
}

extern "C"
const char* casadi_string_get(string_ptr str){
  try{
    return get_stl_string(str).c_str();
  } catch(...){
    return 0;
  }
}


