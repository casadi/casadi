#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/dl/dl-problem.h>
#include <alpaqa/problem/box-constr-problem.hpp>
#include <alpaqa/util/demangled-typename.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace alpaqa::dl {

class DLLoader {
  protected:
    /// Load a shared library.
    DLLoader(
        /// Filename of the shared library to load.
        std::string so_filename,
        /// Prefix of the symbols in the library.
        std::string symbol_prefix);

  protected:
    std::string so_filename;
    std::string symbol_prefix;

    using dl_handle_t = std::shared_ptr<void>;
    /// Handle to the shared library (returned by `dlopen`).
    dl_handle_t handle;

    /// An associative array of additional functions exposed by the problem.
    std::shared_ptr<function_dict_t> extra_functions;

    /// Open the shared library using `dlopen`
    [[nodiscard]] std::shared_ptr<void> load_lib() const;
    /// Load a function with signature @p F from the library using `dlsym`.
    template <class F>
    [[nodiscard]] F *load_func(std::string_view name) const;

    template <class Signature>
        requires std::is_function_v<Signature>
    const std::function<Signature> &extra_func(const std::string &name) const {
        if (!extra_functions)
            throw std::out_of_range("DLProblem: no extra functions");
        auto it = extra_functions->dict.find(name);
        if (it == extra_functions->dict.end())
            throw std::out_of_range("DLProblem: no extra function named \"" +
                                    name + '"');
        try {
            return std::any_cast<const std::function<Signature> &>(it->second);
        } catch (const std::bad_any_cast &e) {
            throw std::logic_error(
                "DLProblem: incorrect type for extra function \"" + name +
                "\" (stored type: " + demangled_typename(it->second.type()) +
                ')');
        }
    }

  public:
    /// Unique type for calling an extra function that is a member function.
    struct instance_t;

  protected:
    template <class Func>
    struct FuncTag {};

    template <class Ret, class... FArgs, class... Args>
    decltype(auto)
    call_extra_func_helper(const void *instance,
                           FuncTag<Ret(const instance_t *, FArgs...)>,
                           const std::string &name, Args &&...args) const {
        return extra_func<Ret(const void *, FArgs...)>(name)(
            instance, std::forward<Args>(args)...);
    }

    template <class Ret, class... FArgs, class... Args>
    decltype(auto)
    call_extra_func_helper(void *instance, FuncTag<Ret(instance_t *, FArgs...)>,
                           const std::string &name, Args &&...args) {
        return extra_func<Ret(void *, FArgs...)>(name)(
            instance, std::forward<Args>(args)...);
    }

    template <class Ret, class... FArgs, class... Args>
    decltype(auto) call_extra_func_helper(const void *, FuncTag<Ret(FArgs...)>,
                                          const std::string &name,
                                          Args &&...args) const {
        return extra_func<Ret(FArgs...)>(name)(std::forward<Args>(args)...);
    }
};

/// Class that loads a problem using `dlopen`.
///
/// The shared library should export a C function with the name
/// `<symbol_prefix>_register` that accepts a void pointer with user data, and
/// returns a struct of type @ref alpaqa_problem_register_t that contains all
/// data to represent the problem, as well as function pointers for all
/// required operations. See @ref C++/DLProblem/main.cpp
///
/// @note   Copies are shallow, they all share the same problem instance, take
///         that into account when using multiple threads.
///
/// @ingroup    grp_Problems
/// @see @ref   TypeErasedProblem
class DLProblem : private DLLoader, public BoxConstrProblem<DefaultConfig> {
  public:
    USING_ALPAQA_CONFIG(DefaultConfig);

    /// Load a problem from a shared library.
    DLProblem(
        /// Filename of the shared library to load.
        std::string so_filename,
        /// Prefix of the registration function.
        std::string symbol_prefix = "alpaqa_problem",
        /// Pointer to custom user data to pass to the registration function.
        void *user_param = nullptr);

  private:
    /// Problem instance created by the registration function, including the
    /// deleter to destroy it.
    std::shared_ptr<void> instance;
    /// Pointer to the struct of function pointers for evaluating the objective,
    /// constraints, their gradients, etc.
    problem_functions_t *functions = nullptr;

  public:
    // clang-format off
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const;
    real_t eval_f(crvec x) const;
    void eval_grad_f(crvec x, rvec grad_fx) const;
    void eval_g(crvec x, rvec gx) const;
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const;
    void eval_jac_g(crvec x, rindexvec inner_idx, rindexvec outer_ptr, rvec J_values) const;
    length_t get_jac_g_num_nonzeros() const;
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const;
    void eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v, rvec Hv) const;
    void eval_hess_L(crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const;
    length_t get_hess_L_num_nonzeros() const;
    void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv) const;
    void eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const;
    length_t get_hess_ψ_num_nonzeros() const;
    real_t eval_f_grad_f(crvec x, rvec grad_fx) const;
    real_t eval_f_g(crvec x, rvec g) const;
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const;
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const;
    real_t eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const;
    void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const;
    real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const;

    [[nodiscard]] bool provides_eval_f() const;
    [[nodiscard]] bool provides_eval_grad_f() const;
    [[nodiscard]] bool provides_eval_g() const;
    [[nodiscard]] bool provides_eval_grad_g_prod() const;
    [[nodiscard]] bool provides_eval_jac_g() const;
    [[nodiscard]] bool provides_get_jac_g_num_nonzeros() const;
    [[nodiscard]] bool provides_eval_grad_gi() const;
    [[nodiscard]] bool provides_eval_hess_L_prod() const;
    [[nodiscard]] bool provides_eval_hess_L() const;
    [[nodiscard]] bool provides_get_hess_L_num_nonzeros() const;
    [[nodiscard]] bool provides_eval_hess_ψ_prod() const;
    [[nodiscard]] bool provides_eval_hess_ψ() const;
    [[nodiscard]] bool provides_get_hess_ψ_num_nonzeros() const;
    [[nodiscard]] bool provides_eval_f_grad_f() const;
    [[nodiscard]] bool provides_eval_f_g() const;
    [[nodiscard]] bool provides_eval_grad_f_grad_g_prod() const;
    [[nodiscard]] bool provides_eval_grad_L() const;
    [[nodiscard]] bool provides_eval_ψ() const;
    [[nodiscard]] bool provides_eval_grad_ψ() const;
    [[nodiscard]] bool provides_eval_ψ_grad_ψ() const;
    [[nodiscard]] bool provides_get_box_C() const;
    // clang-format on

    using instance_t = DLLoader::instance_t;

    template <class Signature, class... Args>
    decltype(auto) call_extra_func(const std::string &name,
                                   Args &&...args) const {
        return extra_func_helper(instance.get(), FuncTag<Signature>{}, name,
                                 std::forward<Args>(args)...);
    }

    template <class Signature, class... Args>
    decltype(auto) call_extra_func(const std::string &name, Args &&...args) {
        return call_extra_func_helper(instance.get(), FuncTag<Signature>{},
                                      name, std::forward<Args>(args)...);
    }
};

#if ALPAQA_WITH_OCP

/// Class that loads an optimal control problem using `dlopen`.
///
/// The shared library should export a C function with the name
/// `<symbol_prefix>_register` that accepts a void pointer with user data, and
/// returns a struct of type @ref alpaqa_control_problem_register_t that
/// contains all data to represent the problem, as well as function pointers for
/// all required operations. See @ref C++/DLProblem/main.cpp
///
/// @note   Copies are shallow, they all share the same problem instance, take
///         that into account when using multiple threads.
///
/// @ingroup    grp_Problems
/// @see @ref   TypeErasedControlProblem
class DLControlProblem : private DLLoader {
  public:
    USING_ALPAQA_CONFIG(DefaultConfig);
    using Box = alpaqa::Box<config_t>;

    /// Load a problem from a shared library.
    DLControlProblem(
        /// Filename of the shared library to load.
        std::string so_filename,
        /// Prefix of the registration function.
        std::string symbol_prefix = "alpaqa_control_problem",
        /// Pointer to custom user data to pass to the registration function.
        void *user_param = nullptr);

  private:
    /// Problem instance created by the registration function, including the
    /// deleter to destroy it.
    std::shared_ptr<void> instance;
    /// Pointer to the struct of function pointers for evaluating the objective,
    /// constraints, their gradients, etc.
    control_problem_functions_t *functions = nullptr;

  public:
    length_t get_N() const { return functions->N; }
    length_t get_nx() const { return functions->nx; }
    length_t get_nu() const { return functions->nu; }
    length_t get_nh() const { return functions->nh; }
    length_t get_nh_N() const { return functions->nh_N; }
    length_t get_nc() const { return functions->nc; }
    length_t get_nc_N() const { return functions->nc_N; }

    void check() const {} // TODO

    // clang-format off
    void get_U(Box &U) const;
    void get_D(Box &D) const;
    void get_D_N(Box &D) const;
    void get_x_init(rvec x_init) const;
    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const;
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const;
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const;
    void eval_h(index_t timestep, crvec x, crvec u, rvec h) const;
    void eval_h_N(crvec x, rvec h) const;
    [[nodiscard]] real_t eval_l(index_t timestep, crvec h) const;
    [[nodiscard]] real_t eval_l_N(crvec h) const;
    void eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const;
    void eval_q_N(crvec x, crvec h, rvec q) const;
    void eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const;
    void eval_add_Q_N(crvec x, crvec h, rmat Q) const;
    void eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R, rvec work) const;
    void eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S, rvec work) const;
    void eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J, crindexvec mask_K, crvec v, rvec out, rvec work) const;
    void eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v, rvec out, rvec work) const;
    [[nodiscard]] length_t get_R_work_size() const;
    [[nodiscard]] length_t get_S_work_size() const;
    void eval_constr(index_t timestep, crvec x, rvec c) const;
    void eval_constr_N(crvec x, rvec c) const;
    void eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const;
    void eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const;
    void eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const;
    void eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const;

    [[nodiscard]] bool provides_get_D() const;
    [[nodiscard]] bool provides_get_D_N() const;
    [[nodiscard]] bool provides_eval_add_Q_N() const;
    [[nodiscard]] bool provides_eval_add_R_prod_masked() const;
    [[nodiscard]] bool provides_eval_add_S_prod_masked() const;
    [[nodiscard]] bool provides_get_R_work_size() const;
    [[nodiscard]] bool provides_get_S_work_size() const;
    [[nodiscard]] bool provides_eval_constr() const;
    [[nodiscard]] bool provides_eval_constr_N() const;
    [[nodiscard]] bool provides_eval_grad_constr_prod() const;
    [[nodiscard]] bool provides_eval_grad_constr_prod_N() const;
    [[nodiscard]] bool provides_eval_add_gn_hess_constr() const;
    [[nodiscard]] bool provides_eval_add_gn_hess_constr_N() const;
    // clang-format on

    using instance_t = DLLoader::instance_t;

    template <class Signature, class... Args>
    decltype(auto) call_extra_func(const std::string &name,
                                   Args &&...args) const {
        return extra_func_helper(instance.get(), FuncTag<Signature>{}, name,
                                 std::forward<Args>(args)...);
    }

    template <class Signature, class... Args>
    decltype(auto) call_extra_func(const std::string &name, Args &&...args) {
        return call_extra_func_helper(instance.get(), FuncTag<Signature>{},
                                      name, std::forward<Args>(args)...);
    }
};

#endif

} // namespace alpaqa::dl