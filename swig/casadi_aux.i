/*
 *    This file is part of CasADi.
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
#include <sstream>
#include "casadi/symbolic/std_vector_tools.hpp"
#include "casadi/symbolic/printable_object.hpp"
#include "casadi/symbolic/shared_object.hpp"
#include "casadi/symbolic/weak_ref.hpp"
#include "casadi/symbolic/generic_type.hpp"
#include "casadi/symbolic/options_functionality.hpp"
#include "casadi/symbolic/casadi_calculus.hpp"
%}

%include "casadi/symbolic/std_vector_tools.hpp"
VECTOR_TOOLS_TEMPLATES(int)
VECTOR_TOOLS_TEMPLATES(double)

%define VECTOR_REPR(type)
%extend std::vector< type >{
  std::string __repr__(){ return casadi::getRepresentation(*$self); }
  std::string __str__(){ return casadi::getDescription(*$self); }
};
%enddef

%include "casadi/symbolic/printable_object.hpp"
%include "casadi/symbolic/shared_object.hpp"
%include "casadi/symbolic/weak_ref.hpp"
%include "casadi/symbolic/casadi_types.hpp"
%include "casadi/symbolic/generic_type.hpp"
%include "casadi/symbolic/options_functionality.hpp"
%include "casadi/symbolic/casadi_calculus.hpp"

namespace casadi {
  %extend OptionsFunctionality {
    void setOption(const std::string &name, const std::string& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, const std::vector<int>& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, const std::vector<double>& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, double val){$self->setOption(name,val);}
    void setOption(const std::string &name, int val){$self->setOption(name,val);} 
    void setOption(const std::string &name, bool val){$self->setOption(name,val);}  
  }
} // namespace casadi

