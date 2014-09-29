/*
 *    $self file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

%{
#include "casadi/core/function/io_interface.hpp"
#include "casadi/core/function/function.hpp"
#include "casadi/core/function/sx_function.hpp"
#include "casadi/core/function/mx_function.hpp"
#include "casadi/core/function/linear_solver.hpp"
#include "casadi/core/function/implicit_function.hpp"
#include "casadi/core/function/integrator.hpp"
#include "casadi/core/function/simulator.hpp"
#include "casadi/core/function/control_simulator.hpp"
#include "casadi/core/function/nlp_solver.hpp"
#include "casadi/core/function/homotopy_nlp_solver.hpp"
#include "casadi/core/function/qp_solver.hpp"
#include "casadi/core/function/stabilized_qp_solver.hpp"
#include "casadi/core/function/lp_solver.hpp"
#include "casadi/core/function/sdp_solver.hpp"
#include "casadi/core/function/socp_solver.hpp"
#include "casadi/core/function/qcqp_solver.hpp"
#include "casadi/core/function/sdqp_solver.hpp"
#include "casadi/core/function/external_function.hpp"
#include "casadi/core/function/parallelizer.hpp"
#include "casadi/core/function/custom_function.hpp"
#include "casadi/core/function/nullspace.hpp"
#include "casadi/core/function/dple_solver.hpp"
#include "casadi/core/function/dle_solver.hpp"
#include "casadi/core/function/lr_dple_solver.hpp"
#include "casadi/core/function/lr_dle_solver.hpp"
#include "casadi/core/function/cle_solver.hpp"
%}

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

%include "casadi/core/function/io_interface.hpp"
%template(IOInterfaceFunction) casadi::IOInterface<casadi::Function>;

%rename(__call__original__) casadi::IOScheme::operator();
%include "casadi/core/function/io_scheme.hpp"
#ifdef SWIGPYTHON
%extend casadi::IOScheme {
%template(__call__original__) operator()< casadi::Sparsity >;
%template(__call__original__) operator()< casadi::MX> ;
%template(__call__original__) operator()< casadi::Matrix<casadi::SXElement> >;
%template(__call__original__) operator()< casadi::Matrix<double> >;

%pythoncode %{
  def __call__(self,*dummy,**kwargs):
    if len(dummy)>1: return self.__call__original__(dummy[0],dummy[1:])
    return self.__call__original__(kwargs.keys(),kwargs.values())

%}


}
#endif

%template(PrintIOSchemeVectorMX)        casadi::PrintableObject<casadi::IOSchemeVector< casadi::MX> >;
%template(PrintIOSchemeVectorSX)        casadi::PrintableObject<casadi::IOSchemeVector< casadi::Matrix<casadi::SXElement> > >;
%template(PrintIOSchemeVectorD)        casadi::PrintableObject<casadi::IOSchemeVector< casadi::Matrix<double> > >;
%template(PrintIOSchemeVectorSparsity)        casadi::PrintableObject<casadi::IOSchemeVector< casadi::Sparsity> >;
%include "casadi/core/function/io_scheme_vector.hpp"

%template(IOSchemeVectorMX) casadi::IOSchemeVector< casadi::MX >;
%template(IOSchemeVectorSX) casadi::IOSchemeVector< casadi::Matrix<casadi::SXElement> >;
%template(IOSchemeVectorD) casadi::IOSchemeVector< casadi::Matrix<double> >;
%template(IOSchemeVectorSparsity) casadi::IOSchemeVector< casadi::Sparsity >;
#ifdef SWIGPYTHON
%extend casadi::IOSchemeVector< casadi::MX > {
%pythoncode %{
  def __iter__(self):
    return iter(self.data)
%}
}
%extend casadi::IOSchemeVector< casadi::Matrix<casadi::SXElement> > {
%pythoncode %{
  def __iter__(self):
    return iter(self.data)
%}
}
%extend casadi::IOSchemeVector< casadi::Matrix<double> > {
%pythoncode %{
  def __iter__(self):
    return iter(self.data)
%}
}
%extend casadi::IOSchemeVector< casadi::Sparsity > {
%pythoncode %{
  def __iter__(self):
    return iter(self.data)
%}
}

#endif

%extend casadi::IOInterface<casadi::Function> {
  casadi::Matrix<double> getInput(int iind=0) const             { static_cast<const casadi::Function*>($self)->assertInit(); return $self->input(iind);}
  casadi::Matrix<double> getInput(const std::string &iname) const             { return $self->input($self->inputSchemeEntry(iname)); }
  casadi::Matrix<double> getOutput(int oind=0) const            { static_cast<const casadi::Function*>($self)->assertInit(); return $self->output(oind);}
}

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
