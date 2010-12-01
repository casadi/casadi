#include "stl_vector_c.hpp"

using namespace std;

std::vector<double>& get_stl_vector(vector_ptr v){
  if(v==0) throw "get_stl_vector failed: null pointer";
  vector<double> *r = (vector<double>*)(v);
  return *r;
}

extern "C"
vector_ptr casadi_vector_new(void){
  try{
    return new vector<double>();
  } catch(...){
    return 0;
  }
}

extern "C"
int casadi_vector_delete(vector_ptr v){
  try{
    vector<double> *ptr = &get_stl_vector(v);
    delete ptr;
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_vector_size(vector_ptr v){
  try{
    return get_stl_vector(v).size();
  } catch(...){
    return -1;
  }
}

extern "C"
const double* casadi_vector_get_ptr(vector_ptr v){
  try{
    return &get_stl_vector(v)[0];
  } catch(...){
    return 0;
  }
}

extern "C"
int casadi_vector_resize(vector_ptr v, int len){
  try{
    get_stl_vector(v).resize(len);
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_vector_set(vector_ptr v, const double* val){
  try{
    vector<double>& vec = get_stl_vector(v);
    copy(val,val+vec.size(),vec.begin());
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_vector_get(vector_ptr v, double* val){
  try{
    const vector<double>& vec = get_stl_vector(v);
    copy(vec.begin(),vec.end(),val);
  } catch(...){
    return 1;
  }  
}

