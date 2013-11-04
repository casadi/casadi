

namespace CasADi {
#ifdef SWIG
%callback("%s_cb");
#endif
void py_c_wrapper(CFunction &f, int nfwd, int nadj, void* user_data){
  PyObject *p = static_cast<PyObject *>(user_data);
  if (PyObject_GetAttrString( p, "__call__")) {
     PyObject * nfwd_py = PyInt_FromLong(nfwd);
     PyObject * nadj_py = PyInt_FromLong(nadj);
     PyObject *r = PyObject_CallFunctionObjArgs(p, nfwd_py, nadj_py, NULL);
     
     if (!r) {
       PyErr_Print();
       throw CasADi::CasadiException("py_c_wrapper: python method execution raised an Error.");
     }
     
     Py_DECREF(nfwd_py);
     Py_DECREF(nadj_py);
     
     if (r) Py_DECREF(r);
  } else {
    throw CasadiException("py_c_wrapper:Supplied python object must be callable.");
  }
  
}
#ifdef SWIG
%nocallback;
#endif


CFunction PyFunction_helper(PyObject * fun) {
  if (!PyObject_GetAttrString( fun, "__call__")) {
    throw CasADi::CasadiException("Supplied python method must be callable.");
  }
  CFunctionWrapper cfunctionwrapper = py_c_wrapper;
  CFunction cfunction(cfunctionwrapper);
  cfunction.setOption("user_data",static_cast<void*>(fun));
  
  return cfunction;
}

}
