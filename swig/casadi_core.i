/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


%module(package="casadi") casadi_core

//  The things needed to make each casadi_*.i  compilable by itself: typemaps
%include "common.i"

%include <casadi/core/function/schemes_metadata.hpp>

// Init hooks
#ifdef SWIGPYTHON
#ifdef WITH_PYTHON_INTERRUPTS
%{
#include <pythonrun.h>
void SigIntHandler(int) {
  std::cerr << "Keyboard Interrupt" << std::endl;
  signal(SIGINT, SIG_DFL);
  kill(getpid(), SIGINT);
}
%}

%init %{
PyOS_setsig(SIGINT, SigIntHandler);
%}
#endif // WITH_PYTHON_INTERRUPTS

%pythoncode%{
try:
  from numpy import pi, inf
except:
  pass

try:
  from numpy import sin, cos, tan, sqrt, log, exp, floor, ceil, fmod, fmin, fmax, sinh, cosh, tanh, arcsin, arccos, arctan, arctan2, fabs, sign, arctanh, arcsinh, arccosh, copysign
except:
  sin = lambda x: x.sin()
  cos = lambda x: x.cos()
  tan = lambda x: x.tan()
  arcsin = lambda x: x.arcsin()
  arccos = lambda x: x.arccos()
  arctan = lambda x: x.arctan()
  sqrt = lambda x: x.sqrt()
  log = lambda x: x.log()
  exp = lambda x: x.exp()
  floor = lambda x: x.floor()
  ceil = lambda x: x.ceil()
  fmin = lambda x,y: x.fmin(y)
  fmax = lambda x,y: x.fmax(y)
  sinh = lambda x: x.sinh()
  cosh = lambda x: x.cosh()
  tanh = lambda x: x.tanh()
  fabs = lambda x: x.fabs()
  sign = lambda x: x.sign()
  arctan2 = lambda x,y: x.arctan2(y)
  arctanh = lambda x: x.arctanh()
  arcsinh = lambda x: x.arcsinh()
  arccosh = lambda x: x.arccosh()
  copysign = lambda x,y: x.copysign(y)
%}
#endif // SWIGPYTHON

#ifdef SWIGMATLAB
%rename(plus) __add__;
%rename(minus) __sub__;
%rename(uminus) operator-;
%rename(uplus) operator+;
%rename(times) __mul__;
%rename(mtimes) mul;
%rename(rdivide) __div__;
%rename(ldivide) __rdiv__;
%rename(mrdivide) __mrdivide__;
%rename(mldivide) __mldivide__;
%rename(power) __pow__;
%rename(mpower) __mpower__;
%rename(lt) __lt__;
%rename(gt) __gt__;
%rename(le) __le__;
%rename(ge) __ge__;
%rename(ne) __ne__;
%rename(eq) __eq__;
//%rename(and) logic_and;
//%rename(or) logic_or;
//%rename(not) logic_not;
%rename(trans) transpose;

// Workarounds, pending proper fix
%rename(truediv) __truediv__;
%rename(nonzero) __nonzero__;
%rename(constpow) __constpow__;
%rename(copysign) __copysign__;
%rename(rpow) __rpow__;
%rename(radd) __radd__;
%rename(rsub) __rsub__;
%rename(rmul) __rmul__;
%rename(rtruediv) __rtruediv__;
%rename(rmldivide) __rmldivide__;
%rename(rmrdivide) __rmrdivide__;
%rename(rmpower) __rmpower__;
%rename(rconstpow) __rconstpow__;
%rename(rge) __rge__;
%rename(rgt) __rgt__;
%rename(rle) __rle__;
%rename(rlt) __rlt__;
%rename(req) __req__;
%rename(rne) __rne__;
%rename(rfmin) __rfmin__;
%rename(rfmax) __rfmax__;
%rename(rarctan2) __rarctan2__;
%rename(rcopysign) __rcopysign__;
%rename(hash) __hash__;
#endif // SWIGMATLAB

#ifdef SWIGPYTHON
%pythoncode %{
def prod(self,*args):
    raise Exception("'prod' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
def dot(self,*args):
    raise Exception("'dot' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
    
class NZproxy:
  def __init__(self,matrix):
    self.matrix = matrix
    
  def __getitem__(self,s):
    return self.matrix.__NZgetitem__(s)

  def __setitem__(self,s,val):
    return self.matrix.__NZsetitem__(s,val)

  def __len__(self):
    return self.matrix.size()
    
  def __iter__(self):
    for k in range(len(self)):
      yield self[k]
%}

%define %matrix_convertors
%pythoncode %{
        
    def toMatrix(self):
        import numpy as n
        return n.matrix(self.toArray())

    def __iter__(self):
      for k in self.nz:
        yield k
        
%}
%enddef 
%define %matrix_helpers(Type)
%pythoncode %{
    @property
    def shape(self):
        return (self.size1(),self.size2())
        
    def reshape(self,arg):
        return _casadi_core.reshape(self,arg)
        
    @property
    def T(self):
        return _casadi_core.transpose(self)
        
    def __getitem__(self,s):
        if isinstance(s,tuple) and len(s)==2:
          return self.__Cgetitem__(s[0],s[1])  
        return self.__Cgetitem__(s)

    def __setitem__(self,s,val):
        if isinstance(s,tuple) and len(s)==2:
          return self.__Csetitem__(s[0],s[1],val)  
        return self.__Csetitem__(s,val)
        
    @property
    def nz(self):
      return NZproxy(self)
        
    def prod(self,*args):
        raise Exception("'prod' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
     
%}
%enddef 
    
%define %python_array_wrappers(arraypriority)
%pythoncode %{

  __array_priority__ = arraypriority

  def __array_wrap__(self,out_arr,context=None):
    if context is None:
      return out_arr
    name = context[0].__name__
    args = list(context[1])
    
    if len(context[1])==3:
      raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b' where 'a' is a numpy type. This is not supported, and cannot be supported without changing numpy." % name)

    if "vectorized" in name:
        name = name[:-len(" (vectorized)")]
    
    conversion = {"multiply": "mul", "divide": "div", "true_divide": "div", "subtract":"sub","power":"pow","greater_equal":"ge","less_equal": "le", "less": "lt", "greater": "gt"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self and not(context[1][0] is self):
      name = 'r' + name
      args.reverse()
    if not(hasattr(self,name)) or ('mul' in name):
      name = '__' + name + '__'
    fun=getattr(self, name)
    return fun(*args[1:])
     
     
  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc) and isinstance(args[1][0],n.ufunc) and len(args[1])>1 and args[1][0].nin==len(args[1][1]):
      if len(args[1][1])==3:
        raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b'. This is not supported when 'a' is a numpy type, and cannot be supported without changing numpy itself. Either upgrade a to a CasADi type first, or use 'a = a + b'. " % args[1][0].__name__)
      return n.array([n.nan])
    else:
      if hasattr(self,'__array_custom__'):
        return self.__array_custom__(*args,**kwargs)
      else:
        return self.toArray()
      
%}
%enddef
#endif // SWIGPYTHON

#ifdef SWIGXML
%define %matrix_helpers(Type)
%enddef
#endif

#ifdef SWIGMATLAB
%define %matrix_helpers(Type)
%enddef
#endif

#ifndef SWIGPYTHON
%define %matrix_convertors
%enddef
#endif

#ifdef SWIGPYTHON
%rename(__Cgetitem__) indexed_zero_based;
%rename(__Cgetitem__) indexed;
%rename(__Csetitem__) indexed_zero_based_assignment;
%rename(__Csetitem__) indexed_assignment;
%rename(__NZgetitem__) nz_indexed_zero_based;
%rename(__NZgetitem__) nz_indexed;
%rename(__NZsetitem__) nz_indexed_zero_based_assignment;
%rename(__NZsetitem__) nz_indexed_assignment;
#endif

#ifndef SWIGXML
%include "typemaps.i"
#endif

%include <casadi/core/std_vector_tools.i>
%include <casadi/core/weak_ref.i>
%include <casadi/core/options_functionality.i>
%include <casadi/core/casadi_calculus.i>
%include <casadi/core/matrix/sparsity.i>
%include <casadi/core/matrix/slice.i>
%include <casadi/core/matrix/generic_expression.i>
%include <casadi/core/matrix/generic_matrix.i>
%include <casadi/core/matrix/matrix.i>
%include <casadi/core/sx/sx_element.i>
%include <casadi/core/mx/mx.i>
%include <casadi/core/matrix/matrix_tools.i>
%include <casadi/core/matrix/generic_matrix_tools.i>
%include <casadi/core/matrix/generic_expression_tools.i>
%include <casadi/core/matrix/sparsity_tools.i>
%include <casadi/core/sx/sx_tools.i>
%include <casadi/core/mx/mx_tools.i>


#ifdef SWIGOCTAVE
%rename(__paren__) indexed_one_based;
#endif

#ifdef SWIGPYTHON
%rename(__getitem__) indexed_zero_based;
#endif

#ifdef SWIGPYTHON
%pythoncode %{
def attach_return_type(f,t):
  if not(hasattr(f,'func_annotations')):
    f.func_annotations = {}
  if not(isinstance(getattr(f,'func_annotations'),dict)):
    raise Exception("Cannot annotate this python Method to be a sparsitygenerator. Method has func_annotations attribute with unknown type.")
  f.func_annotations["return"] = t
  return f

def pyderivativegenerator(f):
  return attach_return_type(f,Function)

def pyevaluate(f):
  return attach_return_type(f,None)
  
def pycallback(f):
  return attach_return_type(f,int)
  

def pyfunction(inputs,outputs):
  def wrap(f):
    
    @pyevaluate
    def fcustom(f2):
      res = f([f2.getInput(i) for i in range(f2.getNumInputs())])
      if not isinstance(res,list):
        res = [res]
      for i in range(f2.getNumOutputs()):
        f2.setOutput(res[i],i)
    Fun = CustomFunction(fcustom,inputs,outputs)
    Fun.setOption("name","CustomFunction")
    return Fun
  return wrap
  
def PyFunction(obj,inputs,outputs):
    @pyevaluate
    def fcustom(f):
      obj.evaluate([f.input(i) for i in range(f.getNumInputs())],[f.output(i) for i in range(f.getNumOutputs())])
      
    Fun = CustomFunction(fcustom,inputs,outputs)
    Fun.setOption("name","CustomFunction")
    if hasattr(obj,'getDerivative'):
      @pyderivativegenerator
      def derivativewrap(f,nfwd,nadj):
        return obj.getDerivative(f,nfwd,nadj)
      Fun.setOption("derivative_generator",derivativewrap)
      
    elif hasattr(obj,'fwd') or hasattr(obj,'adj'):
      @pyderivativegenerator
      def derivativewrap(f,nfwd,nadj):
        num_in = f.getNumInputs()
        num_out = f.getNumOutputs()
        
        @pyevaluate
        def der(f2):
          all_inputs = [f2.input(i) for i in range(f2.getNumInputs())]
          all_outputs = [f2.output(i) for i in range(f2.getNumOutputs())]
          inputs=all_inputs[:num_in]
          outputs=all_outputs[:num_out]
          fwd_seeds=zip(*[iter(all_inputs[num_in:num_in*(nfwd+1)])]*num_in)
          fwd_sens=zip(*[iter(all_outputs[num_out:num_out*(nfwd+1)])]*num_out)
          adj_seeds=zip(*[iter(all_inputs[num_in*(nfwd+1):])]*num_out)
          adj_sens=zip(*[iter(all_outputs[num_out*(nfwd+1):])]*num_in)
          if hasattr(obj,'fwd') and nfwd>0:
            obj.fwd(inputs,outputs,fwd_seeds,fwd_sens)
          if hasattr(obj,'adj') and nadj>0:
            obj.adj(inputs,outputs,adj_seeds,adj_sens)
          
        DerFun = CustomFunction(der,inputs+nfwd*inputs+nadj*outputs,outputs+nfwd*outputs+nadj*inputs)
        DerFun.setOption("name","CustomFunction_derivative")
        DerFun.init()
        return DerFun
 
      Fun.setOption("derivative_generator",derivativewrap)
    
    if not(hasattr(obj,'getDerivative')) and hasattr(obj,'fwd') and not hasattr(obj,'adj'):
      Fun.setOption("ad_mode","forward")
    if not(hasattr(obj,'getDerivative')) and not hasattr(obj,'fwd') and hasattr(obj,'adj'):
      Fun.setOption("ad_mode","reverse")
    return Fun
  
%}
#endif

%include <casadi/core/function/io_interface.i>
%include <casadi/core/function/io_scheme.i>
%include <casadi/core/function/io_scheme_vector.i>

%include <casadi/core/function/function.hpp>
%template(Pair_Function_Function) std::pair<casadi::Function,casadi::Function>;
%include <casadi/core/function/sx_function.hpp>
%include <casadi/core/function/mx_function.hpp>
%include <casadi/core/function/linear_solver.hpp>
%include <casadi/core/function/implicit_function.hpp>
%include <casadi/core/function/integrator.hpp>
%template(IntegratorVector) std::vector<casadi::Integrator>;
%include <casadi/core/function/simulator.hpp>
%include <casadi/core/function/control_simulator.hpp>
%include <casadi/core/function/nlp_solver.hpp>
%include <casadi/core/function/homotopy_nlp_solver.hpp>
%include <casadi/core/function/qp_solver.hpp>
%include <casadi/core/function/stabilized_qp_solver.hpp>
%include <casadi/core/function/lp_solver.hpp>
%include <casadi/core/function/sdp_solver.hpp>
%include <casadi/core/function/socp_solver.hpp>
%include <casadi/core/function/qcqp_solver.hpp>
%include <casadi/core/function/sdqp_solver.hpp>
%include <casadi/core/function/external_function.hpp>
%include <casadi/core/function/parallelizer.hpp>
%include <casadi/core/function/custom_function.hpp>
%include <casadi/core/functor.hpp>
%include <casadi/core/function/nullspace.hpp>
%include <casadi/core/function/dple_solver.hpp>
%include <casadi/core/function/dle_solver.hpp>
%include <casadi/core/function/lr_dple_solver.hpp>
%include <casadi/core/function/lr_dle_solver.hpp>
%include <casadi/core/function/cle_solver.hpp>

%include "autogenerated.i"

%include <casadi/core/casadi_options.hpp>
%include <casadi/core/casadi_meta.hpp>
%include <casadi/core/misc/integration_tools.hpp>
%include <casadi/core/misc/symbolic_nlp.hpp>
%include <casadi/core/misc/variable.hpp>
%include <casadi/core/misc/symbolic_ocp.hpp>
%include <casadi/core/misc/xml_file.hpp>
