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
#ifndef CASADI_IO_INTERFACE_I
#define CASADI_IO_INTERFACE_I

%include <casadi/core/sx/sx_element.i>
 //%include <casadi/core/mx/mx.i>
%include <casadi/core/options_functionality.i>
%include <casadi/core/function/io_scheme.i>

%include <casadi/core/function/io_interface.hpp>

%template(IOInterfaceFunction) casadi::IOInterface<casadi::Function>;

%extend casadi::IOInterface<casadi::Function> {
  casadi::Matrix<double> getInput(int iind=0) const             { static_cast<const casadi::Function*>($self)->assertInit(); return $self->input(iind);}
  casadi::Matrix<double> getInput(const std::string &iname) const             { return $self->input($self->inputSchemeEntry(iname)); }
  casadi::Matrix<double> getOutput(int oind=0) const            { static_cast<const casadi::Function*>($self)->assertInit(); return $self->output(oind);}
}

#endif // CASADI_IO_INTERFACE_I
