%{
#include "casadi/mx/mx.hpp"
#include "casadi/mx/mx_tools.hpp"


%}

%include "casadi/mx/mx.hpp"

%inline %{
template<> swig_type_info** meta< CasADi::MX >::name = &SWIGTYPE_p_CasADi__MX;
template<> swig_type_info** meta< std::vector< CasADi::MX> >::name = &SWIGTYPE_p_std__vectorT_CasADi__MX_std__allocatorT_CasADi__MX_t_t;
template<> swig_type_info** meta< std::vector< CasADi::Matrix<double> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_double_t_std__allocatorT_CasADi__MatrixT_double_t_t_t;
%}


/// std::vector< CasADi::Matrix<double> >
#ifdef SWIGPYTHON
%inline %{
template<> char meta< std::vector< CasADi::Matrix<double> > >::expected_message[] = "Expecting sequence(DMatrix)";

template <>
bool meta< std::vector< CasADi::Matrix<double> > >::couldbe(PyObject * p) {
  return meta< std::vector<  CasADi::Matrix<double> > >::isa(p) || meta< CasADi::Matrix<double> >::couldbe_sequence(p);
}

template <>
int meta< std::vector<  CasADi::Matrix<double> > >::as(PyObject * p,std::vector< CasADi::Matrix<double> > &m) {
  return meta< CasADi::Matrix<double> >::as_vector(p,m);
}

%}
#endif

/// CasADi::MX
#ifdef SWIGPYTHON
%inline %{
template<> char meta< CasADi::MX >::expected_message[] = "Expecting (MX, number)";

template <>
bool meta< CasADi::MX >::couldbe(PyObject * p) {
  return (meta< CasADi::MX >::isa(p) || couldbePyNumber(p) || meta< CasADi::Matrix<double> >::couldbe(p) );
}

template <>
int meta< CasADi::MX >::as(PyObject * p,CasADi::MX &m) {
  NATIVERETURN(CasADi::MX,m)
  if(meta< CasADi::Matrix<double> >::couldbe(p)) {
    CasADi::DMatrix mt;
    bool result=meta< CasADi::Matrix<double> >::as(p,mt);
    if (!result)
      return false;
    m = CasADi::MX(mt);
  } else if (couldbePyNumber(p)) {
    double res;
    int result = getPyNumber(p,&res);
    if (!result)
      return false;
    m=CasADi::MX(res);
  } else {
    return false;
  }
  return true;
}
%}
#endif //SWIGPYTHON


/// std::vector< CasADi::MX >
#ifdef SWIGPYTHON
%inline %{
template<> char meta< std::vector< CasADi::MX > >::expected_message[] = "Expecting sequence(MX, number)";

template <>
int meta< std::vector< CasADi::MX > >::as(PyObject * p,std::vector< CasADi::MX > &m) {
  NATIVERETURN(std::vector< CasADi::MX >,m)
  if( PyIsSequence(p) ) {
    return meta< CasADi::MX >::as_vector(p,m);
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    int num=0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (PyInt_Check(key)) {
        if (PyInt_AsLong(key)+1>num)
          num=PyInt_AsLong(key)+1;
      } else if (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0) {
        if (PyInt_Check(value))
          num=PyInt_AsLong(value);
      } else {
        return false;
      }
    }
    m.resize(num);
    pos=0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (PyInt_Check(key)) {
        bool result=meta< CasADi::MX >::as(value,m[PyInt_AsLong(key)]);
        if (!result)
          return false;
      }
    }
    return true;
  } else {
    return false;
  }
return true;
}

template <>
bool meta< std::vector< CasADi::MX > >::couldbe(PyObject * p) {
  if(meta< CasADi::MX >::couldbe_sequence(p)) {
    return true;
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (!((PyInt_Check(key) || (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0)) && meta< CasADi::MX >::couldbe(value)))
        return false;
    }
    return true;
  }
  return meta< std::vector<CasADi::MX> >::isa(p);
}

%}
#endif //SWIGPYTHON


#ifdef SWIGPYTHON
%my_generic_const_typemap(std::vector< CasADi::MX >,PRECEDENCE_MXVector);
%my_generic_const_typemap(std::vector< CasADi::Matrix<double> >,PRECEDENCE_DMatrixVector);
#endif //SWIGPYTHON


%include "casadi/mx/mx_tools.hpp"

#ifdef SWIGPYTHON
namespace CasADi{
  %extend MX{
    %pythoncode %{
    def __getitem__(self,s):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,list):
        return self.getNZ(s)
      elif isinstance(s,slice):
        s = (s,[0])
      if isinstance(s,tuple):
        if len(s)!=2:
          raise Exception("get/setitem can only do 1D or 2D indexing")
        s = list(s)
        if isinstance(s[0],int) and isinstance(s[1],int):
          for k in range(2):
            if s[k]<0:
              s[k]=s[k]+self.shape[k]
        else:
          for k in range(2):
            if isinstance(s[k],slice):
              J = s[k].indices(self.shape[k])
              s[k] = range(J[0],J[1],J[2])
            elif isinstance(s[k],int):
              if s[k]<0:
                s[k]=s[k]+self.shape[k]
              s[k] = [s[k]]
      return self.getitem(s)
    %}

    %pythoncode %{
    def __setitem__(self,s,val):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,list):
        self.setNZ(s,val)
        return
      elif isinstance(s,slice):
        s = (s,[0])
      if isinstance(s,tuple):
        if len(s)!=2:
          raise Exception("get/setitem can only do 1D or 2D indexing")
        s = list(s)
        if isinstance(s[0],int) and isinstance(s[1],int):
          for k in range(2):
            if s[k]<0:
              s[k]=s[k]+self.shape[k]
        else:
          for k in range(2):
            if isinstance(s[k],slice):
              J = s[k].indices(self.shape[k])
              s[k] = range(J[0],J[1],J[2])
            elif isinstance(s[k],int):
              if s[k]<0:
                s[k]=s[k]+self.shape[k]
              s[k] = [s[k]]
      self.setitem(s,val)
    %}


    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
    %}
    
    %pythoncode %{
        @property
        def T(self):
            return trans(self)
    %}

  };
} // namespace CasADi
#endif // SWIGPYTHON

