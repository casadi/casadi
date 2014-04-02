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
#include "symbolic/function/io_interface.hpp"
#include "symbolic/function/function.hpp"
#include "symbolic/function/sx_function.hpp"
#include "symbolic/function/mx_function.hpp"
#include "symbolic/function/linear_solver.hpp"
#include "symbolic/function/symbolic_qr.hpp"
#include "symbolic/function/implicit_function.hpp"
#include "symbolic/function/integrator.hpp"
#include "symbolic/function/simulator.hpp"
#include "symbolic/function/control_simulator.hpp"
#include "symbolic/function/nlp_solver.hpp"
#include "symbolic/function/homotopy_nlp_solver.hpp"
#include "symbolic/function/qp_solver.hpp"
#include "symbolic/function/stabilized_qp_solver.hpp"
#include "symbolic/function/lp_solver.hpp"
#include "symbolic/function/ocp_solver.hpp"
#include "symbolic/function/sdp_solver.hpp"
#include "symbolic/function/socp_solver.hpp"
#include "symbolic/function/qcqp_solver.hpp"
#include "symbolic/function/sdqp_solver.hpp"
#include "symbolic/function/external_function.hpp"
#include "symbolic/function/parallelizer.hpp"
#include "symbolic/function/custom_function.hpp"
#include "symbolic/function/nullspace.hpp"
%}

#ifdef SWIGOCTAVE
%rename(__paren__) indexed_one_based;
#endif

#ifdef SWIGPYTHON
%rename(__getitem__) indexed_zero_based;
#endif

%include "symbolic/function/io_interface.hpp"
%template(IOInterfaceFunction) CasADi::IOInterface<CasADi::Function>;

%rename(__call__original__) CasADi::IOScheme::operator();
%include "symbolic/function/io_scheme.hpp"
#ifdef SWIGPYTHON
%extend CasADi::IOScheme {
%template(__call__original__) operator()< CasADi::Sparsity >;
%template(__call__original__) operator()< CasADi::MX> ;
%template(__call__original__) operator()< CasADi::Matrix<CasADi::SXElement> >;

%pythoncode %{
  def __call__(self,*dummy,**kwargs):
    if len(dummy)>1: return self.__call__original__(dummy[0],dummy[1:])
    return self.__call__original__(kwargs.keys(),kwargs.values())

%}


}
#endif

%include "symbolic/function/io_scheme_vector.hpp"
%template(IOSchemeVectorMX) CasADi::IOSchemeVector< CasADi::MX >;
%template(IOSchemeVectorSX) CasADi::IOSchemeVector< CasADi::Matrix<CasADi::SXElement> >;
%template(IOSchemeVectorSparsity) CasADi::IOSchemeVector< CasADi::Sparsity >;
#ifdef SWIGPYTHON
%extend CasADi::IOSchemeVector< CasADi::MX > {
%pythoncode %{
  def __iter__(self):
    return iter(self.data)
%}
}
%extend CasADi::IOSchemeVector< CasADi::Matrix<CasADi::SXElement> > {
%pythoncode %{
  def __iter__(self):
    return iter(self.data)
%}
}
%extend CasADi::IOSchemeVector< CasADi::Sparsity > {
%pythoncode %{
  def __iter__(self):
    return iter(self.data)
%}
}

#endif
%include "symbolic/function/function.hpp"

%include "symbolic/function/sx_function.hpp"
%include "symbolic/function/mx_function.hpp"
%include "symbolic/function/linear_solver.hpp"
%include "symbolic/function/symbolic_qr.hpp"
%include "symbolic/function/implicit_function.hpp"
%include "symbolic/function/integrator.hpp"
%include "symbolic/function/simulator.hpp"
%include "symbolic/function/control_simulator.hpp"
%include "symbolic/function/nlp_solver.hpp"
%include "symbolic/function/homotopy_nlp_solver.hpp"
%include "symbolic/function/qp_solver.hpp"
%include "symbolic/function/stabilized_qp_solver.hpp"
%include "symbolic/function/lp_solver.hpp"
%include "symbolic/function/ocp_solver.hpp"
%include "symbolic/function/sdp_solver.hpp"
%include "symbolic/function/socp_solver.hpp"
%include "symbolic/function/qcqp_solver.hpp"
%include "symbolic/function/sdqp_solver.hpp"
%include "symbolic/function/external_function.hpp"
%include "symbolic/function/parallelizer.hpp"
%include "symbolic/function/custom_function.hpp"
%include "symbolic/functor.hpp"
%include "symbolic/function/nullspace.hpp"

%template(IntegratorVector) std::vector<CasADi::Integrator>;
%template(Pair_Function_Function) std::pair<CasADi::Function,CasADi::Function>;
