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


#ifndef CASADI_BONMIN_NLP_HPP
#define CASADI_BONMIN_NLP_HPP

#include <BonTMINLP.hpp>
#include <IpIpoptCalculatedQuantities.hpp>
#include <IpIpoptData.hpp>

#include <iostream>

#include <casadi/interfaces/bonmin/casadi_nlpsol_bonmin_export.h>

/// \cond INTERNAL
using namespace Ipopt;
using namespace Bonmin;

namespace casadi {
  // Forward declarations
  class BonminInterface;
  struct BonminMemory;

  class CASADI_NLPSOL_BONMIN_EXPORT BonminUserClass : public TMINLP {

  public:
    BonminUserClass(const BonminInterface& bonminInterface, BonminMemory* mem);
    virtual ~BonminUserClass();

    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u);

    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda);

    /** Method to return the objective value */
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

    /** Method to return the constraint residuals */
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

    /** Method to return:
     *   1) The structure of the Jacobian (if "values" is NULL)
     *   2) The values of the Jacobian (if "values" is not NULL)
     */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index *jCol,
                            Number* values);

    /** Method to return:
     *   1) The structure of the hessian of the Lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the Lagrangian (if "values" is not NULL)
     */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values);

    /** This method is called when the algorithm is complete
     * so the TNLP can store/write the solution */
    virtual void finalize_solution(TMINLP::SolverReturn status,
                                   Index n, const Number* x,
                                   Number obj_value);

    /** Specify the number of variables that appear in the Hessian */
    virtual Index get_number_of_nonlinear_variables();

    /** Specify which variables that appear in the Hessian */
    virtual bool get_list_of_nonlinear_variables(Index num_nonlin_vars, Index* pos_nonlin_vars);


    /** This method is called at every iteration */
    virtual bool intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value,
                                       Number inf_pr, Number inf_du,
                                       Number mu, Number d_norm,
                                       Number regularization_size,
                                       Number alpha_du, Number alpha_pr,
                                       Index ls_trials,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq);

  virtual bool get_variables_types(Index n, VariableType* var_types);
  virtual bool get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types);
  virtual bool get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types);
  virtual const SosInfo * sosConstraints() const {return 0;}
  virtual const BranchingInfo* branchingInfo() const {return 0;}

  private:
    BonminUserClass(const BonminUserClass&);
    BonminUserClass& operator=(const BonminUserClass&);
    const BonminInterface& solver_;
    BonminMemory* mem_;

    double * x_;
    double * z_L_;
    double * z_U_;
    double * g_;
    double * lambda_;
    int n_;
    int m_;
    double obj_value_;
  };

} // namespace casadi
/// \endcond

#endif //BONMIN_NLP
