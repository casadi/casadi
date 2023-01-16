#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/interop/dl/dl-problem.h>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/demangled-typename.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace alpaqa::dl {

/// Class that loads a problem using dlopen.
///
/// The shared library should export a C function with the name
/// `<symbol_prefix>_register` that accepts a void pointer with user data, and
/// returns a struct of type @ref alpaqa_problem_register_t that contains all
/// data to represent the problem, as well as function pointers for all
/// required operations. See @ref C++/DLProblem/main.cpp
///
/// @note   Copies are shallow, they all share the same problem instance, take
///         that into account when using multiple threads.
class DLProblem : public BoxConstrProblem<EigenConfigd> {
  public:
    USING_ALPAQA_CONFIG(EigenConfigd);

    /// Load a problem from a shared library.
    DLProblem(
        /// Filename of the shared library to load.
        std::string so_filename,
        /// Prefix of the registration function.
        std::string symbol_prefix = "alpaqa_problem",
        /// Pointer to custom user data to pass to the registration function.
        void *user_param = nullptr);

    /// Unique type
    struct instance_t;

  private:
    std::string so_filename;
    std::string symbol_prefix;

    using dl_handle_t = std::shared_ptr<void>;
    /// Handle to the shared library (returned by `dlopen`).
    dl_handle_t handle;
    /// Problem instance created by the registration function, including the
    /// deleter to destroy it.
    std::shared_ptr<void> instance;
    /// Pointer to the struct of function pointers for evaluating the objective,
    /// constraints, their gradients, etc.
    problem_functions_t *functions = nullptr;
    /// An associative array of additional functions exposed by the problem.
    std::shared_ptr<function_dict_t> extra_functions;

    /// Open the shared library using `dlopen`
    [[nodiscard]] std::shared_ptr<void> load_lib() const;
    /// Load a function with signature @p F from the library using `dlsym`.
    template <class F>
    [[nodiscard]] F *load_func(std::string_view name) const;

  public:
    // clang-format off
    real_t eval_f(crvec x) const;
    void eval_grad_f(crvec x, rvec grad_fx) const;
    void eval_g(crvec x, rvec gx) const;
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const;
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const;
    void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const;
    void eval_hess_L(crvec x, crvec y, rmat H) const;
    real_t eval_f_grad_f(crvec x, rvec grad_fx) const;
    real_t eval_f_g(crvec x, rvec g) const;
    real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const;
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const;
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const;
    real_t eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const;
    void eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n) const;
    void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const;
    real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const;

    [[nodiscard]] bool provides_eval_f() const;
    [[nodiscard]] bool provides_eval_grad_f() const;
    [[nodiscard]] bool provides_eval_g() const;
    [[nodiscard]] bool provides_eval_grad_g_prod() const;
    [[nodiscard]] bool provides_eval_grad_gi() const;
    [[nodiscard]] bool provides_eval_hess_L_prod() const;
    [[nodiscard]] bool provides_eval_hess_L() const;
    [[nodiscard]] bool provides_eval_f_grad_f() const;
    [[nodiscard]] bool provides_eval_f_g() const;
    [[nodiscard]] bool provides_eval_f_grad_f_g() const;
    [[nodiscard]] bool provides_eval_grad_f_grad_g_prod() const;
    [[nodiscard]] bool provides_eval_grad_L() const;
    [[nodiscard]] bool provides_eval_ψ() const;
    [[nodiscard]] bool provides_eval_grad_ψ_from_ŷ() const;
    [[nodiscard]] bool provides_eval_grad_ψ() const;
    [[nodiscard]] bool provides_eval_ψ_grad_ψ() const;
    // clang-format on

    template <class Signature>
        requires std::is_function_v<Signature>
    const std::function<Signature> &extra_func(const std::string &name) const {
        if (!extra_functions)
            throw std::logic_error("DLProblem: no extra functions");
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

    template <class Func>
    struct FuncTag {};

    template <class Ret, class... FArgs, class... Args>
    decltype(auto)
    call_extra_func_helper(FuncTag<Ret(const instance_t *, FArgs...)>,
                           const std::string &name, Args &&...args) const {
        return extra_func<Ret(const void *, FArgs...)>(name)(
            instance.get(), std::forward<Args>(args)...);
    }

    template <class Ret, class... FArgs, class... Args>
    decltype(auto) call_extra_func_helper(FuncTag<Ret(instance_t *, FArgs...)>,
                                          const std::string &name,
                                          Args &&...args) {
        return extra_func<Ret(const void *, FArgs...)>(name)(
            instance.get(), std::forward<Args>(args)...);
    }

    template <class Ret, class... FArgs, class... Args>
    decltype(auto) call_extra_func_helper(FuncTag<Ret(FArgs...)>,
                                          const std::string &name,
                                          Args &&...args) {
        return extra_func<Ret(FArgs...)>(name)(std::forward<Args>(args)...);
    }

    template <class Signature, class... Args>
    decltype(auto) call_extra_func(const std::string &name, Args &&...args) {
        return call_extra_func_helper(FuncTag<Signature>{}, name,
                                      std::forward<Args>(args)...);
    }
};

} // namespace alpaqa::dl