#include "sleqp_func.hpp"

#include <sleqp.h>

#include "sleqp_interface.hpp"

namespace casadi {

  static SLEQP_RETCODE casadi_func_set(SleqpFunc* func,
                                       SleqpVec* value,
                                       SLEQP_VALUE_REASON reason,
                                       bool* reject,
                                       void* func_data)
  {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(func_data);

    SLEQP_CALL(sleqp_vec_to_raw(value, m->x));

    return SLEQP_OKAY;
  }

  static SLEQP_RETCODE casadi_func_nonzeros(SleqpFunc* func,
                                            int* obj_grad_nnz,
                                            int* cons_val_nnz,
                                            int* cons_jac_nnz,
                                            int* hess_prod_nnz,
                                            void* func_data)
  {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(func_data);
    return SLEQP_OKAY;
  }

  static SLEQP_RETCODE casadi_func_obj_val(SleqpFunc* func,
                                           double* obj_val,
                                           void* func_data)
  {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(func_data);

    m->arg[0] = m->x;
    m->arg[1] = m->d_nlp.p;
    m->res[0] = obj_val;

    try {

      if(m->interface->calc_function(m, "nlp_f") == 0) {
        return SLEQP_OKAY;
      }

    } catch (std::exception& ex) {
      std::cout << "Error: " << ex.what() << std::endl;
    }

    return SLEQP_ERROR;
  }

  static SLEQP_RETCODE casadi_func_obj_grad(SleqpFunc* func,
                                            SleqpVec* obj_grad,
                                            void* func_data)
  {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(func_data);

    m->arg[0] = m->x;
    m->arg[1] = m->d_nlp.p;
    m->res[0] = nullptr;
    m->res[1] = m->grad_fk;

    try {

      if(m->interface->calc_function(m, "nlp_grad_f") != 0) {
        return SLEQP_ERROR;
      }

      const int num_vars = sleqp_func_num_vars(func);

      SLEQP_CALL(sleqp_vec_set_from_raw(obj_grad,
                                        m->grad_fk,
                                        num_vars,
                                        0.));

    } catch (std::exception& ex) {
      std::cout << "Error: " << ex.what() << std::endl;
    }

    return SLEQP_OKAY;
  }

  static SLEQP_RETCODE casadi_func_cons_val(SleqpFunc* func,
                                            SleqpVec* cons_val,
                                            void* func_data)
  {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(func_data);

    m->arg[0] = m->x;
    m->arg[1] = m->d_nlp.p;
    m->res[0] = m->gk;

    try {

      if(m->interface->calc_function(m, "nlp_g") != 0) {
        return SLEQP_ERROR;
      }

      const int num_cons = sleqp_func_num_cons(func);

      SLEQP_CALL(sleqp_vec_set_from_raw(cons_val,
                                        m->gk,
                                        num_cons,
                                        0.));

    } catch (std::exception& ex) {
      std::cout << "Error: " << ex.what() << std::endl;
    }

    return SLEQP_OKAY;
  }

  static SLEQP_RETCODE casadi_func_cons_jac(SleqpFunc* func,
                                            SleqpMat* cons_jac,
                                            void* func_data)
  {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(func_data);

    m->arg[0] = m->x;
    m->arg[1] = m->d_nlp.p;
    m->res[0] = nullptr;
    m->res[1] = m->jac_gk;

    try {

      if(m->interface->calc_function(m, "nlp_jac_g") != 0) {
        return SLEQP_ERROR;
      }

    } catch (std::exception& ex) {
      std::cout << "Error: " << ex.what() << std::endl;
    }

    casadi_int ncol = m->interface->jacg_sp_.size2();
    const casadi_int* colind = m->interface->jacg_sp_.colind();
    const casadi_int* row = m->interface->jacg_sp_.row();

    const casadi_int nnz = colind[ncol];

    SLEQP_CALL(sleqp_mat_reserve(cons_jac, nnz));

    double* cons_jac_data = sleqp_mat_data(cons_jac);
    int* cons_jac_rows = sleqp_mat_rows(cons_jac);
    int* cons_jac_cols = sleqp_mat_cols(cons_jac);

    for(int k = 0; k < nnz; ++k) {
      cons_jac_data[k] = m->jac_gk[k];
      cons_jac_rows[k] = row[k];
    }

    for(int k = 0; k <= ncol; ++k) {
      cons_jac_cols[k] = colind[k];
    }

    SLEQP_CALL(sleqp_mat_set_nnz(cons_jac, nnz));

    return SLEQP_OKAY;
  }

  static SLEQP_RETCODE casadi_func_hess_prod(SleqpFunc* func,
                                             const SleqpVec* direction,
                                             const SleqpVec* cons_duals,
                                             SleqpVec* product,
                                             void* func_data)
  {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(func_data);
    return SLEQP_OKAY;
  }

  static SLEQP_RETCODE casadi_func_free(void* func_data)
  {
    SLEQPMemory* m = static_cast<SLEQPMemory*>(func_data);
    return SLEQP_OKAY;
  }

  void casadi_sleqp_func_create(SleqpFunc** star,
                                int num_vars,
                                int num_cons,
                                SLEQPMemory* m)
  {
    SleqpFuncCallbacks callbacks{};

    callbacks.set_value = casadi_func_set;
    callbacks.nonzeros  = casadi_func_nonzeros;
    callbacks.obj_val   = casadi_func_obj_val;
    callbacks.obj_grad  = casadi_func_obj_grad;
    callbacks.cons_val  = casadi_func_cons_val;
    callbacks.cons_jac  = casadi_func_cons_jac;
    callbacks.hess_prod = casadi_func_hess_prod;
    callbacks.func_free = casadi_func_free;


    SLEQP_CALL_EXC(sleqp_func_create(star,
                                     &callbacks,
                                     num_vars,
                                     num_cons,
                                     static_cast<void*>(m)));
  }

} // namespace casadi
