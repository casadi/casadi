// NOLINTBEGIN(modernize-use-using,modernize-deprecated-headers)

#pragma once

#include <stddef.h>
#include <stdint.h>

#define ALPAQA_DL_ABI_VERSION 0xA1A000000005

#ifdef ALPAQA_DL_PROBLEM_EXPORT
#elif defined(_WIN32)
#define ALPAQA_DL_PROBLEM_EXPORT
#elif defined(__GNUC__)
#define ALPAQA_DL_PROBLEM_EXPORT __attribute__((visibility("default")))
#else
#define ALPAQA_DL_PROBLEM_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#define ALPAQA_DEFAULT(...) {__VA_ARGS__}
#else
#define ALPAQA_DEFAULT(...)
#endif

typedef double alpaqa_real_t;
typedef ptrdiff_t alpaqa_length_t;
typedef alpaqa_length_t alpaqa_index_t;
typedef uint64_t alpaqa_dl_abi_version_t;

typedef enum {
    /// No type was specified (discouraged).
    alpaqa_register_arg_unspecified = 0,
    /// The @ref alpaqa_register_arg_t::data pointer points to a C++ `std::any`
    /// object. Use `reinterpret_cast<std::any *>(data)` to convert back.
    /// Then use `std::any_cast` to get the actual value out.
    /// @warning `std::any` relies on RTTI that is included in both the
    ///          dynamically loaded problem module and in the code that
    ///          actually loads that module. Types that need to be passed
    ///          between the loading and loaded code need to be exported as
    ///          public classes to ensure that only a single symbol ends up
    ///          in the final application. If the RTTI symbols are private,
    ///          each library will end up with its own local copy, and the
    ///          type information won't compare equal when using libc++,
    ///          causing `std::any_cast` to fail.
    alpaqa_register_arg_std_any = 1,
    /// The @ref alpaqa_register_arg_t::data pointer points to a dynamic C++
    /// `std::span` of `std::string_view`.
    /// Use `reinterpret_cast<std::span<std::string_view> *>(data)` to
    /// convert back.
    alpaqa_register_arg_strings = 2,
    /// The @ref alpaqa_register_arg_t::data pointer points to a C++ tuple of
    /// `pybind11::args` and `pybind11::kwargs`.
    /// Use `reinterpret_cast<std::tuple<pybind11::args, pybind11::kwargs> *>(data)`
    /// to convert back.
    alpaqa_register_arg_py_args = 3,
} alpaqa_register_arg_type_t;

/// User-provided argument that is passed to the problem registration functions.
struct alpaqa_register_arg_t {
    /// Pointer to the user-provided argument passed to the problem registration
    /// functions. Use the @ref type member to determine its type.
    void *data ALPAQA_DEFAULT(nullptr);
    /// Specifies the type of the data pointed to by @ref data.
    alpaqa_register_arg_type_t
        type ALPAQA_DEFAULT(alpaqa_register_arg_unspecified);
};
typedef struct alpaqa_register_arg_t alpaqa_register_arg_t;

/// @see @ref alpaqa::sparsity::Symmetry
typedef enum {
    alpaqa_unsymmetric = 0,
    alpaqa_upper       = 1,
    alpaqa_lower       = 2,
} alpaqa_symmetry;

/// @see @ref alpaqa::sparsity::Dense
struct alpaqa_dense_t {
    alpaqa_length_t rows ALPAQA_DEFAULT(0), cols ALPAQA_DEFAULT(0);
    alpaqa_symmetry symmetry ALPAQA_DEFAULT(alpaqa_unsymmetric);
};
typedef struct alpaqa_dense_t alpaqa_dense_t;

/// @see @ref alpaqa::sparsity::SparseCSC
typedef struct {
    alpaqa_length_t rows, cols;
    alpaqa_symmetry symmetry;
    alpaqa_length_t nnz;
    const int *inner_idx;
    const int *outer_ptr;
    /// @see @ref alpaqa::sparsity::SparseCSC::Order
    enum {
        alpaqa_sparse_csc_unsorted    = 0,
        alpaqa_sparse_csc_sorted_rows = 1,
    } order;
} alpaqa_sparse_csc_t;

/// @see @ref alpaqa::sparsity::SparseCSC
typedef struct {
    alpaqa_length_t rows, cols;
    alpaqa_symmetry symmetry;
    alpaqa_length_t nnz;
    const long *inner_idx;
    const long *outer_ptr;
    /// @see @ref alpaqa::sparsity::SparseCSC::Order
    enum {
        alpaqa_sparse_csc_l_unsorted    = 0,
        alpaqa_sparse_csc_l_sorted_rows = 1,
    } order;
} alpaqa_sparse_csc_l_t;

/// @see @ref alpaqa::sparsity::SparseCSC
typedef struct {
    alpaqa_length_t rows, cols;
    alpaqa_symmetry symmetry;
    alpaqa_length_t nnz;
    const long long *inner_idx;
    const long long *outer_ptr;
    /// @see @ref alpaqa::sparsity::SparseCSC::Order
    enum {
        alpaqa_sparse_csc_ll_unsorted    = 0,
        alpaqa_sparse_csc_ll_sorted_rows = 1,
    } order;
} alpaqa_sparse_csc_ll_t;

/// @see @ref alpaqa::sparsity::SparseCOO
typedef struct {
    alpaqa_length_t rows, cols;
    alpaqa_symmetry symmetry;
    alpaqa_length_t nnz;
    const int *row_indices;
    const int *col_indices;
    /// @see @ref alpaqa::sparsity::SparseCOO::Order
    enum {
        alpaqa_sparse_coo_unsorted                = 0,
        alpaqa_sparse_coo_sorted_by_cols_and_rows = 1,
        alpaqa_sparse_coo_sorted_by_cols_only     = 2,
        alpaqa_sparse_coo_sorted_by_rows_and_cols = 3,
        alpaqa_sparse_coo_sorted_by_rows_only     = 4,
    } order;
    int first_index;
} alpaqa_sparse_coo_t;

/// @see @ref alpaqa::sparsity::SparseCOO
typedef struct {
    alpaqa_length_t rows, cols;
    alpaqa_symmetry symmetry;
    alpaqa_length_t nnz;
    const long *row_indices;
    const long *col_indices;
    /// @see @ref alpaqa::sparsity::SparseCOO::Order
    enum {
        alpaqa_sparse_coo_l_unsorted                = 0,
        alpaqa_sparse_coo_l_sorted_by_cols_and_rows = 1,
        alpaqa_sparse_coo_l_sorted_by_cols_only     = 2,
        alpaqa_sparse_coo_l_sorted_by_rows_and_cols = 3,
        alpaqa_sparse_coo_l_sorted_by_rows_only     = 4,
    } order;
    long first_index;
} alpaqa_sparse_coo_l_t;

/// @see @ref alpaqa::sparsity::SparseCOO
typedef struct {
    alpaqa_length_t rows, cols;
    alpaqa_symmetry symmetry;
    alpaqa_length_t nnz;
    const long long *row_indices;
    const long long *col_indices;
    /// @see @ref alpaqa::sparsity::SparseCOO::Order
    enum {
        alpaqa_sparse_coo_ll_unsorted                = 0,
        alpaqa_sparse_coo_ll_sorted_by_cols_and_rows = 1,
        alpaqa_sparse_coo_ll_sorted_by_cols_only     = 2,
        alpaqa_sparse_coo_ll_sorted_by_rows_and_cols = 3,
        alpaqa_sparse_coo_ll_sorted_by_rows_only     = 4,
    } order;
    long long first_index;
} alpaqa_sparse_coo_ll_t;

/// Sparsity of matrices.
/// @see @ref alpaqa::sparsity::Sparsity
typedef struct {
    union {
        alpaqa_dense_t dense;
        alpaqa_sparse_csc_t sparse_csc;
        alpaqa_sparse_csc_l_t sparse_csc_l;
        alpaqa_sparse_csc_ll_t sparse_csc_ll;
        alpaqa_sparse_coo_t sparse_coo;
        alpaqa_sparse_coo_l_t sparse_coo_l;
        alpaqa_sparse_coo_ll_t sparse_coo_ll;
    };
    enum {
        alpaqa_sparsity_dense,
        alpaqa_sparsity_sparse_csc,
        alpaqa_sparsity_sparse_csc_l,
        alpaqa_sparsity_sparse_csc_ll,
        alpaqa_sparsity_sparse_coo,
        alpaqa_sparsity_sparse_coo_l,
        alpaqa_sparsity_sparse_coo_ll,
    } kind;
} alpaqa_sparsity_t;

/// C API providing function pointers to problem functions.
/// Used by @ref alpaqa::dl::DLProblem.
struct alpaqa_problem_functions_t {
    /// Number of decision variables.
    /// @see @ref alpaqa::TypeErasedProblem::get_num_variables()
    alpaqa_length_t num_variables ALPAQA_DEFAULT(0);
    /// Number of constraints.
    /// @see @ref alpaqa::TypeErasedProblem::get_num_constraints()
    alpaqa_length_t num_constraints ALPAQA_DEFAULT(0);
    /// Name of the problem.
    /// @see @ref alpaqa::TypeErasedProblem::get_name()
    const char *name ALPAQA_DEFAULT(nullptr);

    // clang-format off
    /// Cost function.
    /// @see @ref alpaqa::TypeErasedProblem::eval_objective()
    alpaqa_real_t (*eval_objective)(
        void *instance,
        const alpaqa_real_t *x) ALPAQA_DEFAULT(nullptr);
    /// Gradient of the cost function.
    /// @see @ref alpaqa::TypeErasedProblem::eval_objective_gradient()
    void (*eval_objective_gradient)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *grad_fx) ALPAQA_DEFAULT(nullptr);
    /// Constraints function.
    /// @see @ref alpaqa::TypeErasedProblem::eval_constraints()
    void (*eval_constraints)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *gx) ALPAQA_DEFAULT(nullptr);
    /// Gradient-vector product of the constraints function.
    /// @see @ref alpaqa::TypeErasedProblem::eval_constraints_gradient_product()
    void (*eval_constraints_gradient_product)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t *grad_gxy) ALPAQA_DEFAULT(nullptr);

    /// Difference between point and its projection onto the general constraint
    /// set D.
    /// @see @ref alpaqa::TypeErasedProblem::eval_projecting_difference_constraints()
    void (*eval_projecting_difference_constraints)(
        void *instance,
        const alpaqa_real_t *z,
        alpaqa_real_t *e) ALPAQA_DEFAULT(nullptr);
    /// Project the Lagrange multipliers.
    /// @see @ref alpaqa::TypeErasedProblem::eval_projection_multipliers()
    void (*eval_projection_multipliers)(
        void *instance,
        alpaqa_real_t *y, alpaqa_real_t M);
    /// Proximal gradient step.
    /// @see @ref alpaqa::TypeErasedProblem::eval_proximal_gradient_step()
    /// If not set, the default implementation from
    /// @ref alpaqa::BoxConstrProblem is used.
    alpaqa_real_t (*eval_proximal_gradient_step)(
        void *instance,
        alpaqa_real_t γ,
        const alpaqa_real_t *x,
        const alpaqa_real_t *grad_ψ,
        alpaqa_real_t *x̂,
        alpaqa_real_t *p) ALPAQA_DEFAULT(nullptr);
    /// Active indices for proximal operator.
    /// @see @ref alpaqa::TypeErasedProblem::eval_inactive_indices_res_lna()
    /// If not set, the default implementation from
    /// @ref alpaqa::BoxConstrProblem is used.
    alpaqa_index_t (*eval_inactive_indices_res_lna)(
        void *instance,
        alpaqa_real_t γ,
        const alpaqa_real_t *x,
        const alpaqa_real_t *grad_ψ,
        alpaqa_index_t *J) ALPAQA_DEFAULT(nullptr);

    /// Jacobian of the constraints function.
    /// @see @ref alpaqa::TypeErasedProblem::eval_constraints_jacobian()
    void (*eval_constraints_jacobian)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *J_values) ALPAQA_DEFAULT(nullptr);
    /// Get the sparsity pattern of the Jacobian of the constraints function.
    /// @see @ref alpaqa::TypeErasedProblem::get_constraints_jacobian_sparsity()
    alpaqa_sparsity_t (*get_constraints_jacobian_sparsity)(
        void *instance) ALPAQA_DEFAULT(nullptr);
    /// Gradient of specific constraint function.
    /// @see @ref alpaqa::TypeErasedProblem::eval_grad_gi()
    void (*eval_grad_gi)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_index_t i,
        alpaqa_real_t *grad_gi) ALPAQA_DEFAULT(nullptr);
    /// Hessian-vector product of the Lagrangian.
    /// @see @ref alpaqa::TypeErasedProblem::eval_lagrangian_hessian_product()
    void (*eval_lagrangian_hessian_product)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t scale,
        const alpaqa_real_t *v,
        alpaqa_real_t *Hv) ALPAQA_DEFAULT(nullptr);
    /// Hessian of the Lagrangian.
    /// @see @ref alpaqa::TypeErasedProblem::eval_lagrangian_hessian()
    void (*eval_lagrangian_hessian)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t scale,
        alpaqa_real_t *H_values) ALPAQA_DEFAULT(nullptr);
    /// Get the sparsity pattern of the Hessian of the Lagrangian.
    /// @see @ref alpaqa::TypeErasedProblem::get_lagrangian_hessian_sparsity()
    alpaqa_sparsity_t (*get_lagrangian_hessian_sparsity)(
        void *instance) ALPAQA_DEFAULT(nullptr);
    /// Hessian-vector product of the augmented Lagrangian.
    /// @see @ref alpaqa::TypeErasedProblem::eval_augmented_lagrangian_hessian_product()
    void (*eval_augmented_lagrangian_hessian_product)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        alpaqa_real_t scale,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        const alpaqa_real_t *v,
        alpaqa_real_t *Hv) ALPAQA_DEFAULT(nullptr);
    /// Hessian of the augmented Lagrangian.
    /// @see @ref alpaqa::TypeErasedProblem::eval_augmented_lagrangian_hessian()
    void (*eval_augmented_lagrangian_hessian)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        alpaqa_real_t scale,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        alpaqa_real_t *H_values) ALPAQA_DEFAULT(nullptr);
    /// Get the sparsity pattern of the Hessian of the augmented Lagrangian.
    /// @see @ref alpaqa::TypeErasedProblem::get_augmented_lagrangian_hessian_sparsity()
    alpaqa_sparsity_t (*get_augmented_lagrangian_hessian_sparsity)(
        void *instance) ALPAQA_DEFAULT(nullptr);

    /// Cost and its gradient.
    /// @see @ref alpaqa::TypeErasedProblem::eval_objective_and_gradient()
    alpaqa_real_t (*eval_objective_and_gradient)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *grad_fx) ALPAQA_DEFAULT(nullptr);
    /// Cost and constraints.
    /// @see @ref alpaqa::TypeErasedProblem::eval_objective_and_constraints()
    alpaqa_real_t (*eval_objective_and_constraints)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *g) ALPAQA_DEFAULT(nullptr);
    /// Gradient of the cost and gradient-vector product of the constraints.
    /// @see @ref alpaqa::TypeErasedProblem::eval_objective_gradient_and_constraints_gradient_product()
    void (*eval_objective_gradient_and_constraints_gradient_product)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t *grad_f,
        alpaqa_real_t *grad_gxy) ALPAQA_DEFAULT(nullptr);
    /// Gradient of the Lagrangian.
    /// @see @ref alpaqa::TypeErasedProblem::eval_lagrangian_gradient()
    void (*eval_lagrangian_gradient)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t *grad_L,
        alpaqa_real_t *work_n) ALPAQA_DEFAULT(nullptr);

    /// Augmented Lagrangian.
    /// @see @ref alpaqa::TypeErasedProblem::eval_augmented_lagrangian()
    alpaqa_real_t (*eval_augmented_lagrangian)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        alpaqa_real_t *ŷ) ALPAQA_DEFAULT(nullptr);
    /// Gradient of the augmented Lagrangian.
    /// @see @ref alpaqa::TypeErasedProblem::eval_augmented_lagrangian_gradient()
    void (*eval_augmented_lagrangian_gradient)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        alpaqa_real_t *grad_ψ,
        alpaqa_real_t *work_n,
        alpaqa_real_t *work_m) ALPAQA_DEFAULT(nullptr);
    /// Augmented Lagrangian and its gradient.
    /// @see @ref alpaqa::TypeErasedProblem::eval_augmented_lagrangian_and_gradient()
    alpaqa_real_t (*eval_augmented_lagrangian_and_gradient)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        alpaqa_real_t *grad_ψ,
        alpaqa_real_t *work_n,
        alpaqa_real_t *work_m) ALPAQA_DEFAULT(nullptr);

    /// Provide the initial values for the bounds of
    /// @ref alpaqa::BoxConstrProblem::variable_bounds, i.e. the box constraints
    /// on the decision variables.
    void (*initialize_variable_bounds)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub) ALPAQA_DEFAULT(nullptr);
    /// Provide the initial values for the bounds of
    /// @ref alpaqa::BoxConstrProblem::general_bounds, i.e. the general
    /// constraints.
    void (*initialize_general_bounds)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub) ALPAQA_DEFAULT(nullptr);
    /// Provide the initial values for @ref alpaqa::BoxConstrProblem::l1_reg,
    /// the ℓ₁-regularization factor.
    /// This function is called twice:
    ///  1. with @p lambda set to `nullptr`, and the user should set the size.
    ///  2. with @p lambda pointing to an array of that size, and the user
    ///     should initialize this array.
    void (*initialize_l1_reg)(
        void *instance,
        alpaqa_real_t *lambda,
        alpaqa_length_t *size) ALPAQA_DEFAULT(nullptr);
    // clang-format on
};
typedef struct alpaqa_problem_functions_t alpaqa_problem_functions_t;

/// Opaque type for a C++-only map of extra functions.
typedef struct alpaqa_function_dict_s alpaqa_function_dict_t;
/// Opaque type for a C++-only exception pointer.
typedef struct alpaqa_exception_ptr_s alpaqa_exception_ptr_t;

/// @note When used in C, you should initialize this struct by passing a pointer
///       to your instance to the @ref ALPAQA_PROBLEM_REGISTER_INIT macro.
///       In C++, this is not necessary, because all members have default
///       initializers.
struct alpaqa_problem_register_t {
    /// To check whether the loaded problem is compatible with the version of
    /// the solver.
    alpaqa_dl_abi_version_t abi_version ALPAQA_DEFAULT(ALPAQA_DL_ABI_VERSION);
    /// Owning pointer.
    void *instance ALPAQA_DEFAULT(nullptr);
    /// Non-owning pointer, lifetime at least as long as @ref instance.
    alpaqa_problem_functions_t *functions ALPAQA_DEFAULT(nullptr);
    /// Pointer to the function to clean up @ref instance.
    void (*cleanup)(void *) ALPAQA_DEFAULT(nullptr);
    /// Pointer to a map of extra functions (C++ only).
    /// @see @ref alpaqa::register_function
    /// @see @ref alpaqa::register_member_function
    alpaqa_function_dict_t *extra_functions ALPAQA_DEFAULT(nullptr);
    /// Pointer to an exception that ocurred during problem creation.
    alpaqa_exception_ptr_t *exception ALPAQA_DEFAULT(nullptr);
};
typedef struct alpaqa_problem_register_t alpaqa_problem_register_t;

struct alpaqa_control_problem_functions_t {
    alpaqa_length_t N ALPAQA_DEFAULT(0), nx ALPAQA_DEFAULT(0),
        nu ALPAQA_DEFAULT(0), nh ALPAQA_DEFAULT(0), nh_N ALPAQA_DEFAULT(0),
        nc ALPAQA_DEFAULT(0), nc_N ALPAQA_DEFAULT(0);

    // clang-format off
    void (*get_U)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub) ALPAQA_DEFAULT(nullptr);
    void (*get_D)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub) ALPAQA_DEFAULT(nullptr);
    void (*get_D_N)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub) ALPAQA_DEFAULT(nullptr);
    void (*get_x_init)(
        void *instance,
        alpaqa_real_t *x_init) ALPAQA_DEFAULT(nullptr);
    void (*eval_f)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *u,
        alpaqa_real_t *fxu) ALPAQA_DEFAULT(nullptr);
    void (*eval_jac_f)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *u,
        alpaqa_real_t *J_fxu) ALPAQA_DEFAULT(nullptr);
    void (*eval_grad_f_prod)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *u,
        const alpaqa_real_t *p,
        alpaqa_real_t *grad_fxu_p) ALPAQA_DEFAULT(nullptr);
    void (*eval_h)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *u,
        alpaqa_real_t *h) ALPAQA_DEFAULT(nullptr);
    void (*eval_h_N)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *h) ALPAQA_DEFAULT(nullptr);
    alpaqa_real_t (*eval_l)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *h) ALPAQA_DEFAULT(nullptr);
    alpaqa_real_t (*eval_l_N)(
        void *instance,
        const alpaqa_real_t *h) ALPAQA_DEFAULT(nullptr);
    void (*eval_qr)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        alpaqa_real_t *qr) ALPAQA_DEFAULT(nullptr);
    void (*eval_q_N)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *h,
        alpaqa_real_t *q) ALPAQA_DEFAULT(nullptr);
    void (*eval_add_Q)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        alpaqa_real_t *Q) ALPAQA_DEFAULT(nullptr);
    void (*eval_add_Q_N)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *h,
        alpaqa_real_t *Q) ALPAQA_DEFAULT(nullptr);
    void (*eval_add_R_masked)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        const alpaqa_index_t *mask,
        alpaqa_real_t *R,
        alpaqa_real_t *work) ALPAQA_DEFAULT(nullptr);
    void (*eval_add_S_masked)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        const alpaqa_index_t *mask,
        alpaqa_real_t *S,
        alpaqa_real_t *work) ALPAQA_DEFAULT(nullptr);
    void (*eval_add_R_prod_masked)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        const alpaqa_index_t *mask_J,
        const alpaqa_index_t *mask_K,
        const alpaqa_real_t *v,
        alpaqa_real_t *out,
        alpaqa_real_t *work) ALPAQA_DEFAULT(nullptr);
    void (*eval_add_S_prod_masked)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        const alpaqa_index_t *mask_K,
        const alpaqa_real_t *v,
        alpaqa_real_t *out,
        alpaqa_real_t *work) ALPAQA_DEFAULT(nullptr);
    alpaqa_length_t (*get_R_work_size)(
        void *instance) ALPAQA_DEFAULT(nullptr);
    alpaqa_length_t (*get_S_work_size)(
        void *instance) ALPAQA_DEFAULT(nullptr);
    void (*eval_constr)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        alpaqa_real_t *c) ALPAQA_DEFAULT(nullptr);
    void (*eval_constr_N)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *c) ALPAQA_DEFAULT(nullptr);
    void (*eval_grad_constr_prod)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *p,
        alpaqa_real_t *grad_cx_p) ALPAQA_DEFAULT(nullptr);
    void (*eval_grad_constr_prod_N)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *p,
        alpaqa_real_t *grad_cx_p) ALPAQA_DEFAULT(nullptr);
    void (*eval_add_gn_hess_constr)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *M,
        alpaqa_real_t *out) ALPAQA_DEFAULT(nullptr);
    void (*eval_add_gn_hess_constr_N)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *M,
        alpaqa_real_t *out) ALPAQA_DEFAULT(nullptr);
    // clang-format on
};
typedef struct alpaqa_control_problem_functions_t
    alpaqa_control_problem_functions_t;

/// @note When used in C, you should initialize this struct by passing a pointer
///       to your instance to the @ref ALPAQA_PROBLEM_REGISTER_INIT macro.
///       In C++, this is not necessary, because all members have default
///       initializers.
struct alpaqa_control_problem_register_t {
    /// To check whether the loaded problem is compatible with the version of
    /// the solver.
    alpaqa_dl_abi_version_t abi_version ALPAQA_DEFAULT(ALPAQA_DL_ABI_VERSION);
    /// Owning pointer.
    void *instance ALPAQA_DEFAULT(nullptr);
    /// Non-owning pointer, lifetime at least as long as @ref instance.
    alpaqa_control_problem_functions_t *functions ALPAQA_DEFAULT(nullptr);
    /// Pointer to the function to clean up @ref instance.
    void (*cleanup)(void *) ALPAQA_DEFAULT(nullptr);
    /// Pointer to a map of extra functions (C++ only).
    alpaqa_function_dict_t *extra_functions ALPAQA_DEFAULT(nullptr);
    /// Pointer to an exception that ocurred during problem creation.
    alpaqa_exception_ptr_t *exception ALPAQA_DEFAULT(nullptr);
};
typedef struct alpaqa_control_problem_register_t
    alpaqa_control_problem_register_t;

#ifdef __cplusplus
}
#endif

#if !defined(__cplusplus) || defined(DOXYGEN)
inline static void
alpaqa_problem_register_init(alpaqa_problem_register_t *self) {
    *self = (alpaqa_problem_register_t){
        .abi_version = ALPAQA_DL_ABI_VERSION,
    };
}
inline static void
alpaqa_control_problem_register_init(alpaqa_control_problem_register_t *self) {
    *self = (alpaqa_control_problem_register_t){
        .abi_version = ALPAQA_DL_ABI_VERSION,
    };
}
/// Initialize an instance of @ref alpaqa_problem_register_t or
/// @ref alpaqa_control_problem_register_t. It initializes all members to zero,
/// except for the ABI version, which is initialized to the current ABI version.
/// Available in C only (unnecessary in C++).
/// @param  self
///         A pointer to the instance to initialize.
#define ALPAQA_PROBLEM_REGISTER_INIT(self)                                          \
    _Generic((self),                                                                \
        alpaqa_problem_register_t *: alpaqa_problem_register_init,                  \
        alpaqa_control_problem_register_t *: alpaqa_control_problem_register_init)( \
        self)
#endif

#if defined(__cplusplus) &&                                                    \
    (__cplusplus > 201703L || (defined(_MSVC_LANG) && _MSVC_LANG > 201703L))

#include <any>
#include <exception>
#include <functional>
#include <map>
#include <string>

struct alpaqa_function_dict_s {
    std::map<std::string, std::any> dict{};
};

struct alpaqa_exception_ptr_s {
    std::exception_ptr exc{};
};

namespace alpaqa {

using function_dict_t             = alpaqa_function_dict_t;
using problem_register_t          = alpaqa_problem_register_t;
using control_problem_register_t  = alpaqa_control_problem_register_t;
using problem_functions_t         = alpaqa_problem_functions_t;
using control_problem_functions_t = alpaqa_control_problem_functions_t;

namespace detail {
/// Custom type for which we can export the RTTI to support std::any across
/// shared library boundaries when using libc++.
template <class Signature>
struct ALPAQA_DL_PROBLEM_EXPORT function_wrapper_t {
    std::function<Signature> function;
};
template <class Signature>
function_wrapper_t(std::function<Signature>) -> function_wrapper_t<Signature>;
} // namespace detail

/// Make the given function available to alpaqa.
/// @see @ref alpaqa::dl::DLProblem::call_extra_func
/// @see @ref alpaqa::dl::DLControlProblem::call_extra_func
template <class Func>
void register_function(function_dict_t *&extra_functions, std::string name,
                       Func &&func) {
    if (extra_functions == nullptr)
        extra_functions = new function_dict_t{};
    extra_functions->dict.insert_or_assign(
        std::move(name),
        detail::function_wrapper_t{std::function{std::forward<Func>(func)}});
}

template <class Func>
void register_function(problem_register_t &result, std::string name,
                       Func &&func) {
    register_function(result.extra_functions, std::move(name),
                      std::forward<Func>(func));
}

template <class Func>
void register_function(control_problem_register_t &result, std::string name,
                       Func &&func) {
    register_function(result.extra_functions, std::move(name),
                      std::forward<Func>(func));
}

template <class Result, class T, class Ret, class... Args>
void register_member_function(Result &result, std::string name,
                              Ret (T::*member)(Args...)) {
    register_function(result, std::move(name),
                      [member](void *self_, Args... args) -> Ret {
                          auto *self = reinterpret_cast<T *>(self_);
                          return (self->*member)(std::forward<Args>(args)...);
                      });
}

template <class Result, class T, class Ret, class... Args>
void register_member_function(Result &result, std::string name,
                              Ret (T::*member)(Args...) const) {
    register_function(result, std::move(name),
                      [member](const void *self_, Args... args) -> Ret {
                          const auto *self = reinterpret_cast<const T *>(self_);
                          return (self->*member)(std::forward<Args>(args)...);
                      });
}

namespace detail {
/// Overload for non-const-qualified member functions.
/// @see @ref alpaqa::member_caller
template <auto Member, class Class, class Ret, class... Args>
static auto member_caller(Ret (Class::*)(Args...)) {
    return [](void *self_, Args... args) -> Ret {
        auto *self = reinterpret_cast<Class *>(self_);
        return (self->*Member)(std::forward<Args>(args)...);
    };
}
/// Overload for const-qualified member functions.
/// @see @ref alpaqa::member_caller
template <auto Member, class Class, class Ret, class... Args>
static auto member_caller(Ret (Class::*)(Args...) const) {
    return []<class Self>(Self *self_, Args... args) -> Ret
               requires std::is_void_v<Self>
    {
        const auto *self = reinterpret_cast<const Class *>(self_);
        return (self->*Member)(std::forward<Args>(args)...);
    };
}
/// Overload for member variables.
/// @see @ref alpaqa::member_caller
template <auto Member, class Class, class Ret>
static auto member_caller(Ret Class::*) {
    return []<class Self>(Self *self_) -> decltype(auto)
               requires std::is_void_v<Self>
    {
        using CClass = std::conditional_t<std::is_const_v<Self>,
                                          std::add_const_t<Class>, Class>;
        auto *self   = reinterpret_cast<CClass *>(self_);
        return self->*Member;
    };
}
} // namespace detail

/// Wrap the given member function or variable into a (possibly generic) lambda
/// function that accepts the instance as a type-erased void pointer.
///
/// The signature of the resulting function depends on the signature of the
/// member function:
///
/// - `Ret Class::member(args...)` → `Ret(void *self, args...)`
/// - `Ret Class::member(args...) const` → `Ret(void *self, args...)`
/// - `Ret Class::member(args...) const` → `Ret(const void *self, args...)`
///
/// If the given member is a variable:
///
/// - `Type Class::member` → `Type &(void *self)`
/// - `Type Class::member` → `const Type &(const void *self)`
template <auto Member>
static auto member_caller() {
    return detail::member_caller<Member>(Member);
}

/// Cleans up the extra functions registered by @ref register_function.
/// @note   This does not need to be called for the functions returned by the
///         registration function, those functions will be cleaned up by alpaqa.
/// @note   The @ref alpaqa_problem_register_t and
///         @ref alpaqa_control_problem_register_t structs are part of the C API
///         and do not automatically clean up their resources when destroyed,
///         you have to do it manually by calling this function.
inline void unregister_functions(function_dict_t *&extra_functions) {
    delete extra_functions;
    extra_functions = nullptr;
}

} // namespace alpaqa

#endif

#undef ALPAQA_DEFAULT

// NOLINTEND(modernize-use-using,modernize-deprecated-headers)
