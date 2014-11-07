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
#ifndef CASADI_PRINTABLE_OBJECT_I
#define CASADI_PRINTABLE_OBJECT_I

%include <casadi/core/printable_object.hpp>

%template(PrintSharedObject) casadi::PrintableObject<casadi::SharedObject>;
%template(PrintSlice)        casadi::PrintableObject<casadi::Slice>;
%template(PrintIMatrix)      casadi::PrintableObject<casadi::Matrix<int> >;
%template(PrintDMatrix)      casadi::PrintableObject<casadi::Matrix<double> >;
%template(PrintSXElement)    casadi::PrintableObject<casadi::SXElement>;
//%template(PrintSX)           casadi::PrintableObject<casadi::Matrix<casadi::SXElement> >;
%template(PrintSymbolicNLP)     casadi::PrintableObject<casadi::SymbolicNLP>;
%template(PrintVariable)        casadi::PrintableObject<casadi::Variable>;
%template(PrintSymbolicOCP)     casadi::PrintableObject<casadi::SymbolicOCP>;

%template(PrintIOSchemeVectorMX)     casadi::PrintableObject<casadi::IOSchemeVector< casadi::MX> >;
%template(PrintIOSchemeVectorSX)     casadi::PrintableObject<casadi::IOSchemeVector< casadi::Matrix<casadi::SXElement> > >;
%template(PrintIOSchemeVectorD)      casadi::PrintableObject<casadi::IOSchemeVector< casadi::Matrix<double> > >;
%template(PrintIOSchemeVectorSparsity)        casadi::PrintableObject<casadi::IOSchemeVector< casadi::Sparsity> >;
%template(PrintIOSchemeVectorSparsityVector)  casadi::PrintableObject<casadi::IOSchemeVector< std::vector< casadi::Sparsity> > >;

#endif // CASADI_PRINTABLE_OBJECT_I
