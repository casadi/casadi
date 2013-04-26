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
#include <ctime>

namespace CasADi{

IpoptUserClass::IpoptUserClass(IpoptInternal* solver){
  this->solver = solver;
  n_ = solver->nx_;
  m_ = solver->ng_;
  
#ifdef WITH_IPOPT_CALLBACK 
  x_ = new double[n_];
  g_ = new double[m_];
  z_L_ = new double[n_];
  z_U_ = new double[n_];
  lambda_ = new double[m_];
#endif // WITH_IPOPT_CALLBACK 

}

IpoptUserClass::~IpoptUserClass(){
#ifdef WITH_IPOPT_CALLBACK 
  if (x_) delete [] x_;
  if (g_) delete [] g_;
  if (z_U_) delete [] z_U_;
  if (z_L_) delete [] z_L_;
  if (lambda_) delete [] lambda_;
#endif // WITH_IPOPT_CALLBACK 

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
  solver->finalize_solution(x,z_L,z_U,g,lambda,obj_value, ip_data->iter_count());
}

  
bool IpoptUserClass::intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value,
                                       Number inf_pr, Number inf_du,
                                       Number mu, Number d_norm,
                                       Number regularization_size,
                                       Number alpha_du, Number alpha_pr,
                                       Index ls_trials,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq) {

  // Only do the callback every few iterations
  if (iter % solver->callback_step_!=0) return true;

  /// Code copied from TNLPAdapter::FinalizeSolution
  /// See also: http://list.coin-or.org/pipermail/ipopt/2010-July/002078.html
  // http://list.coin-or.org/pipermail/ipopt/2010-April/001965.html
  
  bool full_callback = false;
    
#ifdef WITH_IPOPT_CALLBACK 
  double time1 = clock();
  OrigIpoptNLP* orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
  if (!orignlp) return true;
  TNLPAdapter* tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
  if (!tnlp_adapter) return true;
  
  const Vector& x = *ip_data->curr()->x();
  const Vector& z_L = *ip_data->curr()->z_L();
  const Vector& z_U = *ip_data->curr()->z_U();
  const Vector& c = *ip_cq->curr_c();
  const Vector& d = *ip_cq->curr_d();
  const Vector& y_c = *ip_data->curr()->y_c();
  const Vector& y_d = *ip_data->curr()->y_d();
  
  std::fill_n(x_, n_, 0);
  std::fill_n(g_, m_, 0);
  std::fill_n(z_L_, n_, 0);
  std::fill_n(z_U_, n_, 0);
  std::fill_n(lambda_, m_, 0);
  
  tnlp_adapter->ResortX(x, x_);             // no further steps needed
  tnlp_adapter->ResortG(y_c, y_d, lambda_); // no further steps needed
  tnlp_adapter->ResortG(c, d, g_);
  // Copied from Ipopt source: To Ipopt, the equality constraints are presented with right
  // hand side zero, so we correct for the original right hand side.
  const Index* c_pos = tnlp_adapter->P_c_g_->ExpandedPosIndices();
  Index n_c_no_fixed = tnlp_adapter->P_c_g_->NCols();
  for (Index i=0; i<n_c_no_fixed; i++) {
    g_[c_pos[i]] += tnlp_adapter->c_rhs_[i];
  }
    
  tnlp_adapter->ResortBnds(z_L, z_L_, z_U, z_U_);
  // Copied from Ipopt source: Hopefully the following is correct to recover the bound
  // multipliers for fixed variables (sign ok?)
  if (tnlp_adapter->fixed_variable_treatment_==TNLPAdapter::MAKE_CONSTRAINT && tnlp_adapter->n_x_fixed_>0) {
    const DenseVector* dy_c = static_cast<const DenseVector*>(&y_c);
    const Number* values = dy_c->Values();
    Index n_c_no_fixed = y_c.Dim() - tnlp_adapter->n_x_fixed_;
    for (Index i=0; i<tnlp_adapter->n_x_fixed_; i++) {
      z_L_[tnlp_adapter->x_fixed_map_[i]] = Max(0., -values[n_c_no_fixed+i]);
      z_U_[tnlp_adapter->x_fixed_map_[i]] = Max(0., values[n_c_no_fixed+i]);
    }
  }
  double time2 = clock();
  solver->t_callback_prepare_ += double(time2-time1)/CLOCKS_PER_SEC;
  full_callback = true;
#endif // WITH_IPOPT_CALLBACK 

  return solver->intermediate_callback(x_,z_L_,z_U_,g_,lambda_,obj_value,iter,inf_pr,inf_du,mu,d_norm,regularization_size,alpha_du,alpha_pr,ls_trials,full_callback);

                 
}


  
Index IpoptUserClass::get_number_of_nonlinear_variables(){
  return solver->get_number_of_nonlinear_variables();
}

bool IpoptUserClass::get_list_of_nonlinear_variables(Index num_nonlin_vars, Index* pos_nonlin_vars){
  return solver->get_list_of_nonlinear_variables(num_nonlin_vars,pos_nonlin_vars);
}

bool IpoptUserClass::get_var_con_metadata(Index n,StringMetaDataMapType& var_string_md, IntegerMetaDataMapType& var_integer_md, NumericMetaDataMapType& var_numeric_md,
                                          Index m,StringMetaDataMapType& con_string_md, IntegerMetaDataMapType& con_integer_md, NumericMetaDataMapType& con_numeric_md){
  
  return solver->get_var_con_metadata(n, var_string_md, var_integer_md, var_numeric_md, m, con_string_md, con_integer_md, con_numeric_md);
}

void IpoptUserClass::finalize_metadata(Index n, const StringMetaDataMapType& var_string_md,
                                       const IntegerMetaDataMapType& var_integer_md, const NumericMetaDataMapType& var_numeric_md,
                                       Index m, const StringMetaDataMapType& con_string_md,
                                       const IntegerMetaDataMapType& con_integer_md, const NumericMetaDataMapType& con_numeric_md){
  solver->finalize_metadata(n,var_string_md,var_integer_md,var_numeric_md,m,con_string_md,con_integer_md,con_numeric_md);
}




} // namespace CasADi
