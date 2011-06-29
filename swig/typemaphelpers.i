%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/matrix.hpp"
#include <sstream>
#include "casadi/casadi_exception.hpp"

// to allow for typechecking
#include "casadi/sx/sx.hpp"

// to typecheck for MX
#include "casadi/mx/mx.hpp"
%}


%include "typemaps.i"

/// Generic typemap structure
%inline %{

#ifdef SWIGPYTHON
#define GUESTOBJECT PyObject * p
#endif // SWIGPYTHON

#ifdef  SWIGOCTAVE
#define GUESTOBJECT const octave_value& p
#endif // SWIGOCTAVE

/** Check if Guest object is of a particular SWIG type */
bool istype(GUESTOBJECT, swig_type_info *type) {
  return SWIG_IsOK(SWIG_ConvertPtr(p, 0, type, 0));
}

template<class T>
class meta {
  public:
    /// Check if Python object is of type T
    static bool isa(GUESTOBJECT) {
      return istype(p,*meta<T>::name);
    };
    /// Convert Python object to pointer of type T
    static bool get_ptr(GUESTOBJECT,T*& m) {
      void *pd = 0 ;
      int res = SWIG_ConvertPtr(p, &pd,*meta<T>::name, 0 );
      if (!SWIG_IsOK(res)) {
        return false;
      }
      m = reinterpret_cast< T*  >(pd);
      return true;
    };
    /// Convert Guest object to type T
    /// This function must work when isa(GUESTOBJECT) too
    static int as(GUESTOBJECT,T& m) {
        T *t = (T *)(0);
        int res = SWIG_CheckState(swig::asptr(p, &t));
        if (res) m=*t;
        return res;
    }
    /// Check if Guest object could ultimately be converted to type T
    /// may return true when isa(GUESTOBJECT), but this is not required.
    static bool couldbe(GUESTOBJECT) { 
        int res = swig::asptr(p, (T**)(0));
        return SWIG_CheckState(res);
    }
    static swig_type_info** name;
    static char expected_message[];
    
    // Vector specific stuff
    
    #ifdef SWIGPYTHON
    static bool couldbe_sequence(PyObject * p) {
      if(PySequence_Check(p) && !meta< CasADi::Matrix<CasADi::SX> >::isa(p) && !meta< CasADi::MX >::isa(p)) {
        PyObject *it = PyObject_GetIter(p);
        PyObject *pe;
        int i=0;
        while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
          if (!meta< T >::couldbe(pe)) {
            Py_DECREF(pe);Py_DECREF(it);return false;
          }
          Py_DECREF(pe);
        }
        Py_DECREF(it);
        return true;
      } else {
        return false;
      }
    }
    #endif // SWIGPYTHON
    
    // Assumes that p is a PYTHON sequence
    #ifdef SWIGPYTHON
    static int as_vector(PyObject * p, std::vector<T> &m) {
      PyObject *it = PyObject_GetIter(p);
      PyObject *pe;
      m.resize(PySequence_Size(p));
      int i=0;
      while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
        bool result=meta< T >::as(pe,m[i++]);
        if (!result) {
          Py_DECREF(pe);Py_DECREF(it);
          return false;
        }
        Py_DECREF(pe);
      }
      Py_DECREF(it);
      return true;
    }
    #endif // SWIGPYTHON
    
    #ifdef SWIGPYTHON
    // Would love to make this const T&, but looks like not allowed
    static bool toPython(T &, PyObject *&p);
    #endif //SWIGPYTHON
};

%}

%inline %{
#define NATIVERETURN(Type, m) if (meta<Type>::isa(p)) { Type *mp; int result = meta<Type>::get_ptr(p,mp); if (!result) return false; m=*mp; return true;}
%}


%define %my_generic_const_typemap(Precedence,Type...) 
%typemap(in) const Type & (Type m) {
  if (meta< Type >::isa($input)) { // Type object get passed on as-is, and fast.
    int result = meta< Type >::get_ptr($input,$1);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,"Type cast failed");
  } else {
    bool result=meta< Type >::as($input,m);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,meta< Type >::expected_message);
    $1 = &m;
  }
}

%typemap(typecheck,precedence=Precedence) const Type & { $1 = meta< Type >::isa($input) || meta< Type >::couldbe($input); }
%typemap(freearg) const Type  & {}

%enddef


/// std::vector< Type >
#ifdef SWIGPYTHON
%define %meta_vector(Type) 
%inline %{
template<> char meta< std::vector< Type > >::expected_message[] = "Expecting sequence(Type)";

template <>
int meta< std::vector< Type > >::as(PyObject * p,std::vector< Type > &m) {
  NATIVERETURN(std::vector< Type >,m)
  return meta< Type >::as_vector(p,m);
}

template <>
bool meta< std::vector< Type > >::couldbe(PyObject * p) {
  return meta< std::vector< Type > >::isa(p) ||  meta< Type >::couldbe_sequence(p);
}
%}
%enddef
#endif //SWIGPYTHON

/// std::pair< TypeA, TypeB >
#ifdef SWIGPYTHON
%define %meta_pair(TypeA,TypeB) 
%inline %{
template <>
int meta< std::pair< TypeA, TypeB > >::as(PyObject * p,std::pair< TypeA, TypeB > &m) {
  if(!PySequence_Check(p)) return false;
  if(PySequence_Size(p)!=2) return false;
  PyObject * first =  PySequence_GetItem(p,0);
  PyObject * second = PySequence_GetItem(p,1);
  bool result = meta< TypeA  >::as(p,m.first) && meta< TypeB  >::as(p,m.second);
  
  Py_DECREF(first);Py_DECREF(second);
  return result;   
}

template <>
bool meta< std::pair< TypeA, TypeB > >::couldbe(PyObject * p) {
  if (meta< std::pair< TypeA, TypeB > >::isa(p)) return true; 
  if(!PySequence_Check(p)) return false;
  if(PySequence_Size(p)!=2) return false;
  PyObject * first =  PySequence_GetItem(p,0);
  PyObject * second = PySequence_GetItem(p,1);
  
  bool success = (meta< TypeA >::isa(first) || meta< TypeA >::couldbe(first)) &&
                 (meta< TypeB >::isa(second) || meta< TypeB >::couldbe(second));
  Py_DECREF(first);Py_DECREF(second);              
  return success;
}

template <>
bool meta < std::pair< TypeA, TypeB > >::toPython(std::pair< TypeA, TypeB > &m, PyObject *&p) {
  p = PyTuple_New(2);
  PyObject *first = 0;
  first  = SWIG_NewPointerObj((new TypeA(static_cast< TypeA& >(m.first ))), *meta< TypeA >::name , SWIG_POINTER_OWN |  0 );
  PyObject *second = 0;
  second = SWIG_NewPointerObj((new TypeB(static_cast< TypeB& >(m.second))), *meta< TypeB >::name , SWIG_POINTER_OWN |  0 );
  
  if (first==0 || second==0) return false;
  
  PyTuple_SetItem(p, 0, first);
  PyTuple_SetItem(p, 1, second);
  
  return true;
}
%}
%enddef
#endif //SWIGPYTHON

#ifdef SWIGPYTHON
%inline%{
/** Check PyObjects by class name */
bool PyObjectHasClassName(PyObject* p, const char * name) {
  PyObject * classo = PyObject_GetAttrString( p, "__class__");
  PyObject * classname = PyObject_GetAttrString( classo, "__name__");
  
  bool ret = strcmp(PyString_AsString(classname),name)==0;
  Py_DECREF(classo);Py_DECREF(classname);
	return ret;
}

bool PyIsSequence(PyObject* p) {
  return PySequence_Check(p) && !meta< CasADi::Matrix<CasADi::SX> >::isa(p) && !meta< CasADi::MX >::isa(p);
}

%}
#endif // SWIGPYTHON

%{
template<> swig_type_info** meta< double >::name = &SWIGTYPE_p_double;
template<> swig_type_info** meta< int >::name = &SWIGTYPE_p_int;
// template<> swig_type_info** meta< std::vector<double> >::name = &SWIGTYPE_p_std__vectorT_double_p_std__allocatorT_double_p_t_t; BUG: Since it was no longer used anywhere, it caused an undeclared type error!
template<> swig_type_info** meta< std::vector<int> >::name = &SWIGTYPE_p_std__vectorT_int_std__allocatorT_int_t_t;
%}

/// int
#ifdef SWIGPYTHON
%inline %{
template<> char meta< int >::expected_message[] = "Expecting number";

template <>
int meta< int >::as(PyObject * p, int &m) {
  NATIVERETURN(int,m)
  NATIVERETURN(double,m)
  PyObject *r = PyNumber_Float(p);
  int ret = !(r==NULL);
  if (ret)
    m = PyFloat_AsDouble(r);
  Py_DECREF(r);
  return ret;
}

template <> bool meta< int >::couldbe(PyObject * p) { return PyInt_Check(p) || PyBool_Check(p) || PyFloat_Check(p) ;}

%}
#endif //SWIGPYTHON

/// double
#ifdef SWIGPYTHON
%inline %{
template<> char meta< double >::expected_message[] = "Expecting number";

template <>
int meta< double >::as(PyObject * p, double &m) {
  NATIVERETURN(double,m)
  PyObject *r = PyNumber_Float(p);
  int ret = !(r==NULL);
  if (ret)
    m = PyFloat_AsDouble(r);
  Py_DECREF(r);
  return ret;
}

template <> bool meta< double >::couldbe(PyObject * p) { return PyInt_Check(p) || PyBool_Check(p) || PyFloat_Check(p) ;}

%}
#endif //SWIGPYTHON

/// std::vector<double>
#ifdef SWIGOCTAVE
%inline %{
template<> char meta< std::vector<double> >::expected_message[] = "Expecting (1xn) array(number)";

template <>
int meta< std::vector<double> >::as(const octave_value& p, std::vector<double> &m) {
  NATIVERETURN(std::vector<double>, m);
  if(p.is_real_matrix()){
    const Matrix &mat = p.matrix_value();
    if (!(mat.rows()==1)) return false;
    m.resize(mat.cols());
    for(int j=0; j<mat.cols(); ++j) m[j] = mat(0,j);
    return true;
  }
}

template <> bool meta< std::vector<double> >::couldbe(const octave_value& p) { 
  if(p.is_real_matrix()){
    const Matrix &mat = p.matrix_value();
    return (mat.rows()==1 );
  } else {
    return false;
  }
}

%}
#endif //SWIGOCTAVE

/// std::vector<int>
#ifdef SWIGOCTAVE
%inline %{
template<> char meta< std::vector<int> >::expected_message[] = "Expecting (1xn) array(number)";

template <>
int meta< std::vector<int> >::as(const octave_value& p, std::vector<int> &m) {
  NATIVERETURN(std::vector<int>, m);
  if(p.is_real_matrix()){
    const Matrix &mat = p.matrix_value();
    if (!(mat.rows()==1)) return false;
    m.resize(mat.cols());
    for(int j=0; j<mat.cols(); ++j) m[j] = mat(0,j);
    return true;
  }
}

template <> bool meta< std::vector<int> >::couldbe(const octave_value& p) { 
  if(p.is_real_matrix()){
    const Matrix &mat = p.matrix_value();
    return (mat.rows()==1 );
  } else {
    return false;
  }
}

%}
#endif //SWIGOCTAVE
