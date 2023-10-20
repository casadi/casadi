/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include <ocp/OCPAbstract.hpp>
#include <ocp/StageOCPApplication.hpp>
#include <casadi/core/casadi_misc.hpp>

/**
  @requires_conic("hpipm")
  @requires_conic("qpoases")
  def test_hpipm_timevarying(self):
    def mat(a):
      def fl(a):
        return float(a) if len(a)>0 else 0
      return sparsify(DM([list(map(fl,i.split("\t"))) for i in a.split("\n") if len(i)>0]))
    def vec(a):
      return DM(list(map(float,a.split("\n"))))
    N = 2
    A = """
1	0.2	1	-1	0	0	0	0	0	0	0	0
-0.1	0.4	0	0	-1	0	0	0	0	0	0	0
0.3	0.2	0	0	0	-1	0	0	0	0	0	0
2	0	0.3	0	0	0	0	0	0	0	0	0
1	1	0.4	0	0	0	0	0	0	0	0	0
	0	0	1	4	2	1	0.3	-1	0	0	0
	0	0	3	1	0	1	0.2	0	-1	0	0
	0	0	1	1	1	1	1	0	0	0	0
	0	0	0	0	0	0	0	2	4	0	-1
	0	0	0	0	0	0	0	2	3	1	0
	0	0	0	0	0	0	0	0	0	0	3"""
    A = mat(A)
    nx = [2,3,2,1]
    nu = [1, 2,1]
    ng = [2, 1, 1, 1]
    N = 3
    H = """7	0	0.2	0	0	0	0	0	0	0	0	0
	7	0.3	0	0	0	0	0	0	0	0	0
0.2	0.3	1	0	0	0	0	0	0	0	0	0
	0	0	3	0	0	0	1	0	0	0	0
	0	0	0	2	0.1	0	0.7	0	0	0	0
	0	0	0	0.1	1	0	1	0	0	0	0
	0	0	0	0	0	1	0.1	0	0	0	0
	0	0	1	0.7	1	0.1	2	0	0	0	0
	0	0	0	0	0	0	0	6	0	1	0
	0	0	0	0	0	0	0	0	6	0	0
	0	0	0	0	0	0	0	1	0	4	0
	0	0	0	0	0	0	0	0	0	0	9
"""
    H = mat(H)
    options = {"hpipm":{"iter_max":100,"res_g_max":1e-10,"res_b_max":1e-10,"res_d_max":1e-10,"res_m_max":1e-10}}
    #solver = conic('solver', 'hpipm', {"a": A.sparsity(), "h": H.sparsity()},{"N":N,"nx":nx,"nu":nu,"ng":ng,"tol":1e-12,"mu0":2,"max_iter":20})
    solver = conic('solver', 'hpipm', {"a": A.sparsity(), "h": H.sparsity()},options)
    solver_ref = conic('solver', 'qpoases', {"a": A.sparsity(), "h": H.sparsity()})

    g = vec("""1
1
0.2
0.4
1
0.5
0.3
1
0.6
1
1
0.7""")
    lbg = vec("""0
    0
    0
    -2
    -2
    0
    0
    -2
    0
    -2
    -2""")

    ubg = vec("""0
    0
    0
    2
    2
    0
    0
    2
    0
    2
    2""")

    lbx = vec("""0.5
    0.2
    -1
    -1
    -1
    -1
    -1
    -1
    -1
    -1
    -1
    -1""")
    ubx = vec("""0.5
    0.2
    1
    1
    1
    1
    1
    1
    1
    1
    1
    1""")

*/

namespace casadi {

class CasadiStructuredQP : public fatrop::OCPAbstract {
  /// @brief number of states for time step k
  /// @param k: time step
  fatrop_int get_nxk(const fatrop_int k) const override {
    std::vector<fatrop_int> res = {2,3,2,1};
    return res[k];
  }
  /// @brief number of inputs for time step k
  /// @param k: time step
  fatrop_int get_nuk(const fatrop_int k) const override {
    std::vector<fatrop_int> res = {1, 2,1};
    return res[k];
  };
  /// @brief number of equality constraints for time step k
  /// @param k: time step
  fatrop_int get_ngk(const fatrop_int k) const override {
    std::vector<fatrop_int> res = {2, 0,0,0};
    return res[k]; // 2 from lbx

  };
  /// @brief  number of stage parameters for time step k
  /// @param k: time step
  fatrop_int get_n_stage_params_k(const fatrop_int k) const override { return 0;}
  /// @brief  number of global parameters
  fatrop_int get_n_global_params() const override { return 0;}
  /// @brief default stage parameters for time step k
  /// @param stage_params: pointer to array of size n_stage_params_k
  /// @param k: time step
  fatrop_int get_default_stage_paramsk(double *stage_params, const fatrop_int k) const override { return 0;}
  /// @brief default global parameters
  /// @param global_params: pointer to array of size n_global_params
  fatrop_int get_default_global_params(double *global_params) const override{ return 0; }
  /// @brief number of inequality constraints for time step k
  /// @param k: time step
  virtual fatrop_int get_ng_ineq_k(const fatrop_int k) const {
    std::vector<fatrop_int> res_lbg = {2, 1, 0, 0};
    std::vector<fatrop_int> res_lbx = {1, 5, 3, 1};
    return res_lbg[k]+res_lbx[k];
  }
  /// @brief horizon length
  fatrop_int get_horizon_length() const override { return 4; }
  /// @brief  discretized dynamics
  /// it evaluates the vertical concatenation of A_k^T, B_k^T, and b_k^T from the linearized dynamics x_{k+1} = A_k x_k + B_k u_k + b_k. 
  /// The matrix is in column major format.
  /// @param states_kp1: pointer to nx_{k+1}-array states of time step k+1
  /// @param inputs_k: pointer to array inputs of time step k
  /// @param states_k: pointer to array states of time step k
  /// @param stage_params_k: pointer to array stage parameters of time step k
  /// @param global_params: pointer to array global parameters
  /// @param res: pointer to (nu+nx+1 x nu+nx)-matrix 
  /// @param k: time step
  fatrop_int eval_BAbtk(
      const double *states_kp1,
      const double *inputs_k,
      const double *states_k,
      const double *stage_params_k,
      const double *global_params,
      MAT *res,
      const fatrop_int k) override {
        printf("eval_BAbtk k=%d\n", k);
        if (k==0) {
          std::vector<double> r = {1, 0.2, 1,   0,
                                      -0.1,	0.4,0 ,  0,
                                      0.3,	0.2,	0, 0};
          int out_m = 4; // rows
          int out_n = 3; // cols
          PACKMAT(out_m, out_n, get_ptr(r), out_m, res, 0, 0);
          blasfeo_print_dmat(out_m, out_n,  res, 0, 0);

          std::vector<double> r2 = {1, -0.1, 0.3, 0.2, 0.4, 0.2, 1, 0, 0, 0, 0, 0};

          blasfeo_pack_tran_dmat(out_n, out_m, get_ptr(r2), out_n, res, 0, 0);

          blasfeo_print_dmat(out_m, out_n,  res, 0, 0);

          //casadi_error("foo");

        } else if (k==1) {
          std::vector<double> r = {1,	4,	2,	1,	0.3, 0,
                                    3,	1,	00,	1,	0.2, 0};
          int out_m = 6; // rows
          int out_n = 2; // cols
          PACKMAT(out_m, out_n, get_ptr(r), out_m, res, 0, 0);
          blasfeo_print_dmat(out_m, out_n,  res, 0, 0);


          std::vector<double> r2 = {1,3,4,1,2,0,1,1,0.3,0.2,0,0,0,0,0};

          //blasfeo_pack_tran_dmat(out_n, out_m, get_ptr(r2), out_n, res, 0, 0);

          //blasfeo_print_dmat(out_m, out_n,  res, 0, 0);

        } else if (k==2) {
          std::vector<double> r = {2, 4, 0, 0};
          int out_m = 4; // rows
          int out_n = 1; // cols
          PACKMAT(out_m, out_n, get_ptr(r), out_m, res, 0, 0);
          blasfeo_print_dmat(out_m, out_n,  res, 0, 0);
        }

        
      };
  /// @brief  stagewise Lagrangian Hessian
  /// It evaluates is the vertical concatenation of (1) the Hessian of the Lagrangian to the concatenation of (u_k, x_k) (2) the first order derivative of the Lagrangian Hessian to the concatenation of (u_k, x_k). 
  /// The matrix is in column major format.
  /// @param objective_scale: scale factor for objective function (usually 1.0)
  /// @param inputs_k: pointer to array inputs of time step k
  /// @param states_k: pointer to array states of time step k
  /// @param lam_dyn_k: pointer to array dual variables for dynamics of time step k
  /// @param lam_eq_k: pointer to array dual variables for equality constraints of time step k
  /// @param lam_eq_ineq_k: pointer to array dual variables for inequality constraints of time step k
  /// @param stage_params_k: pointer to array stage parameters of time step k
  /// @param global_params: pointer to array global parameters
  /// @param res: pointer to (nu+nx+1 x nu+nx)-matrix. 
  /// @param k
  /// @return
  fatrop_int eval_RSQrqtk(
      const double *objective_scale,
      const double *inputs_k,
      const double *states_k,
      const double *lam_dyn_k,
      const double *lam_eq_k,
      const double *lam_eq_ineq_k,
      const double *stage_params_k,
      const double *global_params,
      MAT *res,
      const fatrop_int k) override {

      }
  /// @brief stagewise equality constraints Jacobian. 
  /// It evaluates the vertical concatenation of (1) the Jacobian of the equality constraints to the concatenation of (u_k, x_k) (2) the equality constraints evaluated at u_k, x_k.
  /// The matrix is in column major format.
  /// @param inputs_k: pointer to array inputs of time step k
  /// @param states_k: pointer to array states of time step k
  /// @param stage_params_k: pointer to array stage parameters of time step k
  /// @param global_params: pointer to array global parameters
  /// @param res: pointer to (nu+nx+1 x ng)-matrix.
  /// @param k: time step
  /// @return
  fatrop_int eval_Ggtk(
      const double *inputs_k,
      const double *states_k,
      const double *stage_params_k,
      const double *global_params,
      MAT *res,
      const fatrop_int k) override {
    printf("eval_Ggtk k=%d\n", k);
    if (k==0) {
      std::vector<double> r = {2, 0, 0.3, 0,
                              1, 1, 0.4, 0};
      int out_m = 4; // rows
      int out_n = 2; // cols
      PACKMAT(out_m, out_n, get_ptr(r), out_m, res, 0, 0);
      blasfeo_print_dmat(out_m, out_n,  res, 0, 0);
    }
  }
  /// @brief stagewise inequality constraints Jacobian. 
  /// It evaluates the vertical concatenation of (1) the Jacobian of the inequality constraints to the concatenation of (u_k, x_k) (2) the inequality constraints evaluated at u_k, x_k. 
  /// The matrix is in column major format.
  /// @param inputs_k: pointer to array inputs of time step k
  /// @param states_k: pointer to array states of time step k
  /// @param stage_params_k: pointer to array stage parameters of time step k
  /// @param global_params_ko: pointer to array global parameters
  /// @param res: pointer to (nu+nx+1 x ng_ineq)-matrix, column major format
  /// @param k : time step
  /// @return
  fatrop_int eval_Ggt_ineqk(
      const double *inputs_k,
      const double *states_k,
      const double *stage_params_k,
      const double *global_params,
      MAT *res,
      const fatrop_int k) {
    printf("eval_Ggt_ineqk k=%d\n", k);
  }
  /// @brief the dynamics constraint violation (b_k = -x_{k+1} + f_k(u_k, x_k, p_k, p))
  /// @param states_kp1: pointer to array states of time step k+1
  /// @param inputs_k: pointer to array inputs of time step k
  /// @param states_k: pointer to array states of time step k
  /// @param stage_params_k: pointer to array stage parameters of time step k
  /// @param global_params: pointer to array global parameters
  /// @param res: pointer to array nx_{k+1}-vector
  /// @param k: time step
  /// @return
  fatrop_int eval_bk(
      const double *states_kp1,
      const double *inputs_k,
      const double *states_k,
      const double *stage_params_k,
      const double *global_params,
      double *res,
      const fatrop_int k) override {
  printf("eval_bk k=%d\n", k);
      }
  /// @brief the equality constraint violation (g_k = g_k(u_k, x_k, p_k, p))
  /// @param inputs_k: pointer to array inputs of time step k
  /// @param states_k: pointer to array states of time step k
  /// @param stage_params_k: pointer to array stage parameters of time step k
  /// @param global_params: pointer to array global parameters
  /// @param res: pointer to array ng-vector
  /// @param k: time step
  fatrop_int eval_gk(
      const double *states_k,
      const double *inputs_k,
      const double *stage_params_k,
      const double *global_params,
      double *res,
      const fatrop_int k) override {
    printf("eval_gk k=%d\n", k);
  }
  /// @brief the inequality constraint violation (g_ineq_k = g_ineq_k(u_k, x_k, p_k, p))
  /// @param inputs_k: pointer to array inputs of time step k
  /// @param states_k: pointer to array states of time step k
  /// @param stage_params_k: pointer to array stage parameters of time step k
  /// @param global_params: pointer to array global parameters
  /// @param res: pointer to array ng_ineq-vector
  /// @param k: time step
  fatrop_int eval_gineqk(
      const double *states_k,
      const double *inputs_k,
      const double *stage_params_k,
      const double *global_params,
      double *res,
      const fatrop_int k)  override {
      printf("eval_gineqk k=%d\n", k);
  }
  /// @brief gradient of the objective function (not the Lagrangian!) to the concatenation of (u_k, x_k)
  /// @param objective_scale: pointer to objective scale
  /// @param inputs_k: pointer to array inputs of time step k
  /// @param states_k: pointer to array states of time step k
  /// @param stage_params_k: pointer to array stage parameters of time step k
  /// @param global_params: pointer to array global parameters
  /// @param res: pointer to (nu+nx)-array
  /// @param k: time step
  fatrop_int eval_rqk(
      const double *objective_scale,
      const double *inputs_k,
      const double *states_k,
      const double *stage_params_k,
      const double *global_params,
      double *res,
      const fatrop_int k) override {
      printf("eval_rqk k=%d\n", k);
  }
  /// @brief objective function value 
  /// @param objective_scale: pointer to array objective scale
  /// @param inputs_k: pointer to array inputs of time step k
  /// @param states_k: pointer to array states of time step k
  /// @param stage_params_k: pointer to array stage parameters of time step k
  /// @param global_params: pointer to array global parameters
  /// @param res: pointer to double
  /// @param k: time step
  fatrop_int eval_Lk(
      const double *objective_scale,
      const double *inputs_k,
      const double *states_k,
      const double *stage_params_k,
      const double *global_params,
      double *res,
      const fatrop_int k) override {
      printf("eval_Lk k=%d\n", k);
  }
  /// @brief the bounds of the inequalites at stage k
  /// @param lower: pointer to ng_ineq-vector
  /// @param upper: pointer to ng_ineq-vector
  /// @param k: time step
  fatrop_int get_boundsk(double *lower, double *upper, const fatrop_int k) const override {
      printf("get_boundsk k=%d\n", k);
  }
  /// @brief default initial guess for the states of stage k
  /// @param xk: pointer to states of time step k 
  /// @param k: time step
  fatrop_int get_initial_xk(double *xk, const fatrop_int k) const override {
    printf("get_initial_xk k=%d\n", k);
  }
  /// @brief default initial guess for the inputs of stage k
  /// @param uk: pointer to inputs of time step k
  /// @param k: time step
  fatrop_int get_initial_uk(double *uk, const fatrop_int k) const override {
    printf("get_initial_uk k=%d\n", k);
  }

};

} // namespace casadi

int main() {

  int m = 5;
  int n = 7;


  size_t s = blasfeo_memsize_dmat(m, n);
  std::vector<double> mem(s, 0);

  blasfeo_dmat res;
  blasfeo_create_dmat(m, n, &res, &mem.front());

  blasfeo_print_dmat(m, n, &res, 0, 0);


  std::vector<double> a(100);
  for (int i=0;i<100;++i) a[i] = i+1;


  //blasfeo_pack_dmat(2, 3, &a.front(), 10, &res, 1, 3);

  //blasfeo_pack_tran_dmat(2, 3, &a.front(), 10, &res, 1, 3);
  blasfeo_pack_tran_dmat(1, 3, &a.front(), 2, &res, 1, 3);

  blasfeo_print_dmat(m, n, &res, 0, 0);
  return 0;

  /*blasfeo_pack_dmat(int m, int n, double *A, int lda, struct blasfeo_dmat *sB, int bi, int bj);

  blasfeo_pack_tran_dmat(1, n, d->CD+p->CD_offsets[k]+(i-start), ng_ineq, res, 0, column++);

  blasfeo_print_dmat(m, n,  res, 0, 0);*/



  casadi::CasadiStructuredQP qp;

  fatrop::OCPApplication app(std::make_shared<casadi::CasadiStructuredQP>(qp));
  app.build();

  app.optimize();

}