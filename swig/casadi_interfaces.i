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


#ifdef WITH_OOQP
%{
#include "interfaces/ooqp/ooqp_solver.hpp"
%}
%include "interfaces/ooqp/ooqp_solver.hpp"
#endif

#ifdef WITH_SQIC
%{
#include "interfaces/sqic/sqic_solver.hpp"
#include "interfaces/sqic/stabilized_sqic_solver.hpp"
%}
%include "interfaces/sqic/sqic_solver.hpp"
%include "interfaces/sqic/stabilized_sqic_solver.hpp"
#endif

// IPOPT
#ifdef WITH_IPOPT
%include "ipopt_interface.i"
#endif

// QPOASES
#ifdef WITH_QPOASES
%include "qpoases_interface.i"
%{
#include "interfaces/qpoases/qpoases_solver.hpp"
%}
%include "interfaces/qpoases/qpoases_solver.hpp"
#endif

// Sundials
#ifdef WITH_SUNDIALS
%include "sundials_interface.i"
#endif

// Slicot
#ifdef WITH_SLICOT
%include "slicot_interface.i"
#endif

#ifdef WITH_LAPACK
%include "lapack_interface.i"
#endif

#ifdef WITH_KNITRO
%{ 
  #include "interfaces/knitro/knitro_solver.hpp"
%}
%include "interfaces/knitro/knitro_solver.hpp"
#endif

#ifdef WITH_CPLEX
%{ 
  #include "interfaces/cplex/cplex_solver.hpp"
%}
%include "interfaces/cplex/cplex_solver.hpp"
#endif

#ifdef WITH_LIFTOPT
%{
  #include "interfaces/liftopt/liftopt_solver.hpp"
%}
%include "interfaces/liftopt/liftopt_solver.hpp"
#endif

#ifdef WITH_CSPARSE
%{
#include "interfaces/csparse/csparse.hpp"
#include "interfaces/csparse/csparse_cholesky.hpp"
%}
%include "interfaces/csparse/csparse.hpp"
%include "interfaces/csparse/csparse_cholesky.hpp"
#endif

#ifdef WITH_DSDP
%{
#include "interfaces/dsdp/dsdp_solver.hpp"
%}
%include "interfaces/dsdp/dsdp_solver.hpp"
#endif

#ifdef WITH_GSL
%{
#include "interfaces/gsl/gsl_integrator.hpp"
%}
%include "interfaces/gsl/gsl_integrator.hpp"
#endif

#ifdef WITH_WORHP
%{
#include "interfaces/worhp/worhp_solver.hpp"
%}
%include "interfaces/worhp/worhp_solver.hpp"
#endif


#ifdef WITH_SNOPT
%{
#include "interfaces/snopt/snopt_solver.hpp"
%}
%include "interfaces/snopt/snopt_solver.hpp"
#endif
