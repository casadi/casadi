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

#include "ipopt_nlp.hpp"
#include "ipopt_internal.hpp"

namespace CasADi{

IpoptUserClass::IpoptUserClass(IpoptInternal* solver){
  this->solver = solver;
}

IpoptUserClass::~IpoptUserClass(){
}

// returns the size of the problem
bool IpoptUserClass::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  solver->get_nlp_info(n,m,nnz_jac_g,nnz_h_lag);

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool IpoptUserClass::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  return solver->get_bounds_info(n,x_l,x_u,m,g_l,g_u);
}

// returns the initial point for the problem
bool IpoptUserClass::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  return solver->get_starting_point(n,init_x,x,init_z,z_L,z_U,m,init_lambda,lambda);
}

// returns the value of the objective function
bool IpoptUserClass::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  return solver->eval_f(n,x,new_x,obj_value);
}

// return the gradient of the objective function grad_{x} f(x)
bool IpoptUserClass::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  return solver->eval_grad_f(n,x,new_x,grad_f);
}

// return the value of the constraints: g(x)
bool IpoptUserClass::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  return solver->eval_g(n,x,new_x,m,g);
}

// return the structure or values of the jacobian
bool IpoptUserClass::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  return solver->eval_jac_g(n,x,new_x,m,nele_jac,iRow,jCol,values);
}


bool IpoptUserClass::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{  
  return solver->eval_h(x,new_x,obj_factor,lambda,new_lambda,nele_hess,iRow,jCol,values);
}


void IpoptUserClass::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
				  const IpoptData* ip_data,
				  IpoptCalculatedQuantities* ip_cq)
{
  solver->finalize_solution(x,z_L,z_U,g,lambda,obj_value);
}

Index IpoptUserClass::get_number_of_nonlinear_variables(){
  return solver->get_number_of_nonlinear_variables();
}

bool IpoptUserClass::get_list_of_nonlinear_variables(Index num_nonlin_vars, Index* pos_nonlin_vars){
  return solver->get_list_of_nonlinear_variables(num_nonlin_vars,pos_nonlin_vars);
}




} // namespace CasADi
