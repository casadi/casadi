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

%include "casadi/core/function/schemes_metadata.hpp"

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
#endif // SWIGPYTHON

%include "casadi/core/std_vector_tools.i"
%include "casadi/core/weak_ref.i"
%include "casadi/core/options_functionality.i"
%include "casadi/core/casadi_calculus.i"

#ifdef SWIGPYTHON
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

%include "matrix.i"

#ifndef SWIGXML
%include "typemaps.i"
#endif

%include <casadi/core/matrix/sparsity.i>
%include <casadi/core/matrix/slice.i>
%include <casadi/core/matrix/generic_expression.i>
%include <casadi/core/matrix/generic_matrix.i>
%include <casadi/core/matrix/matrix.i>
%include <casadi/core/sx/sx_element.i>

%include "mx.i"

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

%include "casadi/core/function/function.hpp"
%include "casadi/core/function/sx_function.hpp"
%include "casadi/core/function/mx_function.hpp"
%include "casadi/core/function/linear_solver.hpp"
%include "casadi/core/function/implicit_function.hpp"
%include "casadi/core/function/integrator.hpp"
%include "casadi/core/function/simulator.hpp"
%include "casadi/core/function/control_simulator.hpp"
%include "casadi/core/function/nlp_solver.hpp"
%include "casadi/core/function/homotopy_nlp_solver.hpp"
%include "casadi/core/function/qp_solver.hpp"
%include "casadi/core/function/stabilized_qp_solver.hpp"
%include "casadi/core/function/lp_solver.hpp"
%include "casadi/core/function/sdp_solver.hpp"
%include "casadi/core/function/socp_solver.hpp"
%include "casadi/core/function/qcqp_solver.hpp"
%include "casadi/core/function/sdqp_solver.hpp"
%include "casadi/core/function/external_function.hpp"
%include "casadi/core/function/parallelizer.hpp"
%include "casadi/core/function/custom_function.hpp"
%include "casadi/core/functor.hpp"
%include "casadi/core/function/nullspace.hpp"
%include "casadi/core/function/dple_solver.hpp"
%include "casadi/core/function/dle_solver.hpp"
%include "casadi/core/function/lr_dple_solver.hpp"
%include "casadi/core/function/lr_dle_solver.hpp"
%include "casadi/core/function/cle_solver.hpp"

%template(IntegratorVector) std::vector<casadi::Integrator>;
%template(Pair_Function_Function) std::pair<casadi::Function,casadi::Function>;

// autogenerated
%include "autogenerated.i"

%include "casadi/core/casadi_options.hpp"
%include "casadi/core/casadi_meta.hpp"
%include "casadi/core/misc/integration_tools.hpp"

%template(PrintSymbolicNLP)        casadi::PrintableObject<casadi::SymbolicNLP>;
%include "casadi/core/misc/symbolic_nlp.hpp"

%template(PrintVariable)        casadi::PrintableObject<casadi::Variable>;
%include "casadi/core/misc/variable.hpp"

%template(PrintSymbolicOCP)        casadi::PrintableObject<casadi::SymbolicOCP>;
%include "casadi/core/misc/symbolic_ocp.hpp"

%include "casadi/core/misc/xml_file.hpp"
