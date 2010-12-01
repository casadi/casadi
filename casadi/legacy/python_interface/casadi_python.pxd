cdef extern from "casadi/c_interface/stl_string_c.h":
    ctypedef void* string_ptr
    ctypedef char* const_char_ptr "const char*"

    string_ptr casadi_string_new()
    int casadi_string_delete(string_ptr str)
    int casadi_string_assign(string_ptr str, const_char_ptr s)
    const_char_ptr casadi_string_get(string_ptr str)

cdef extern from "casadi/c_interface/stl_vector_c.h":
    ctypedef void* vector_ptr
    ctypedef double* const_double_ptr "const double*"

    vector_ptr casadi_vector_new()
    int casadi_vector_delete(vector_ptr str)
    int casadi_vector_size(vector_ptr v)
    const_double_ptr casadi_vector_get_ptr(vector_ptr v)
    int casadi_vector_resize(vector_ptr v, int len)
    int casadi_vector_set(vector_ptr v, const_double_ptr val)
    int casadi_vector_get(vector_ptr v, double* val)

cdef extern from "casadi/c_interface/sx_c.h":
    ctypedef void* sx_ptr
    ctypedef void* sx_vec_ptr
    sx_ptr casadi_sx_new()
    int casadi_sx_symbol(sx_ptr ptr, const_char_ptr name)
    int casadi_sx_constant(sx_ptr ptr, double value)
    int casadi_sx_delete(sx_ptr ptr)
    int casadi_sx_print(sx_ptr ptr, string_ptr str)
    int casadi_sx_binary(sx_ptr r, int op, sx_ptr x, sx_ptr y)
    int casadi_sx_unary(sx_ptr r, int op, sx_ptr x)
    sx_vec_ptr casadi_sx_vec_new()
    int casadi_sx_vec_delete(sx_vec_ptr v)
    int casadi_sx_vec_push_back(sx_vec_ptr v, sx_ptr ptr)
    int casadi_sx_vec_size(sx_vec_ptr v)
    int casadi_sx_vec_print(sx_vec_ptr ptr, string_ptr str)

cdef extern from "casadi/c_interface/sx_matrix_c.h":
    ctypedef void* sx_matrix_ref
    ctypedef void* sx_matrix_vec

    sx_matrix_ref casadi_sx_matrix_new()
    int casadi_sx_matrix_symbol(sx_matrix_ref ref, const_char_ptr name, int nrow, int ncol)
    int casadi_sx_matrix_constant(sx_matrix_ref ref, const_double_ptr data, int nrow, int ncol, char order)
    int casadi_sx_matrix_sx(sx_matrix_ref ref, sx_ptr scalar)
    int casadi_sx_matrix_sx_vec(sx_matrix_ref ref, sx_vec_ptr v)
    int casadi_sx_matrix_delete(sx_matrix_ref ref)
    int casadi_sx_matrix_print(sx_matrix_ref ref, string_ptr str)
    int casadi_sx_matrix_binary(sx_matrix_ref r, int op, sx_matrix_ref x, sx_matrix_ref y)
    int casadi_sx_matrix_unary(sx_matrix_ref r, int op, sx_matrix_ref x)
    int casadi_sx_matrix_prod(sx_matrix_ref r, sx_matrix_ref x, sx_matrix_ref y)
    int casadi_sx_matrix_vertcat(sx_matrix_ref r, sx_matrix_vec v)
    int casadi_sx_matrix_horzcat(sx_matrix_ref r, sx_matrix_vec v)
    sx_matrix_vec casadi_sx_matrix_vec_new()
    int casadi_sx_matrix_vec_delete(sx_matrix_vec v)
    int casadi_sx_matrix_vec_push_back(sx_matrix_vec v, sx_matrix_ref ref)
    int casadi_sx_matrix_vec_size(sx_matrix_vec v)
    int casadi_sx_matrix_vec_print(sx_matrix_vec ref, string_ptr str)
    int casadi_sx_matrix_size(sx_matrix_ref ref, int *sz)
    int casadi_sx_matrix_size1(sx_matrix_ref ref, int *sz)
    int casadi_sx_matrix_size2(sx_matrix_ref ref, int *sz)
    int casadi_sx_matrix_transpose(sx_matrix_ref res, sx_matrix_ref ref)

cdef extern from "casadi/c_interface/mx_c.h":
    ctypedef void* mx_ref
    ctypedef void* mx_vec

    mx_ref casadi_mx_new()
    int casadi_mx_symbol(mx_ref ref, const_char_ptr name, int nrow, int ncol)
    int casadi_mx_constant(mx_ref ref, const_double_ptr  data, int nrow, int ncol, char order)
    int casadi_mx_delete(mx_ref ref)
    int casadi_mx_print_string(mx_ref ref, string_ptr str)
    int casadi_mx_binary(mx_ref r, int op, mx_ref x, mx_ref y)
    int casadi_mx_unary(mx_ref r, int op, mx_ref x)

    mx_vec casadi_mx_vec_new()
    int casadi_mx_vec_delete(mx_vec v)
    int casadi_mx_vec_push_back(mx_vec v, mx_ref ref)
    int casadi_mx_vec_size(mx_vec v)

cdef extern from "casadi/c_interface/fx_c.h":
    ctypedef void* fx_ref
    ctypedef void* mxarray "const mx_ref*"
    
    fx_ref casadi_fx_new()
    int casadi_fx_delete(fx_ref ref)
    int casadi_fx_print_cout(fx_ref ref)
    int casadi_fx_print_cerr(fx_ref ref)
    int casadi_fx_print_string(fx_ref ref, string_ptr str)
    int casadi_fx_setoption_string(fx_ref ref, const_char_ptr, const_char_ptr)
    int casadi_fx_setoption_double(fx_ref ref, const_char_ptr, const_double_ptr, int n)
    int casadi_fx_getoption_string(fx_ref ref, const_char_ptr name, string_ptr str)
    int casadi_fx_getoption_double(fx_ref ref, const_char_ptr name, vector_ptr str)
    int casadi_fx_option_is_string(fx_ref ref, const_char_ptr name, int *is_string)
    int casadi_fx_print_option(fx_ref ref, const_char_ptr name, string_ptr str)
    int casadi_fx_print_options(fx_ref ref)
    int casadi_fx_input_size(fx_ref ref, int ind, int *sz)
    int casadi_fx_getinput(fx_ref ref, int ind, int ord, double* val)
    int casadi_fx_setinput(fx_ref ref, int ind, int ord, const_double_ptr val)
    int casadi_fx_output_size(fx_ref ref, int ind, int *sz)
    int casadi_fx_getoutput(fx_ref ref, int ind, int ord, double* val)
    int casadi_fx_setoutput(fx_ref ref, int ind, int ord, const_double_ptr val)
    int casadi_fx_init(fx_ref ref)
    int casadi_fx_evaluate(fx_ref ref)
    int casadi_fx_evaluate_fwd(fx_ref ref)
    int casadi_fx_evaluate_adj(fx_ref ref)

cdef extern from "casadi/c_interface/mx_function_c.h":
    int casadi_mx_function(fx_ref fcn, mx_vec iv, mx_ref ov)

cdef extern from "casadi/c_interface/sx_function_c.h":
    int casadi_sx_function(fx_ref fcn, sx_matrix_vec iv, sx_matrix_vec ov)

cdef extern from "casadi/c_interface/integrator_c.h":
    ctypedef void* integrator_ref
    int casadi_integrator_integrate(fx_ref ref, double t_out)
    int casadi_integrator_reset(fx_ref ref, int with_sens)

cdef extern from "sundials_interface/sundials_interface_c.h":
    int casadi_cvodes_integrator(fx_ref fcn, fx_ref ffcn)
    int casadi_idas_integrator(fx_ref fcn, fx_ref ffcn)

cdef extern from "ipopt_interface/ipopt_solver_c.h":
    int casadi_ipopt_solver(fx_ref fcn, fx_ref f, fx_ref g, fx_ref h, fx_ref j)
