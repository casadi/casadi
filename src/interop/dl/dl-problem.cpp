#include <alpaqa/interop/dl/dl-problem.hpp>

#include <dlfcn.h>
#include <memory>
#include <stdexcept>

namespace alpaqa::dl {

std::shared_ptr<void> DLProblem::load_lib() const {
    assert(!so_filename.empty());
    ::dlerror();
    void *h = ::dlopen(so_filename.c_str(), RTLD_LOCAL | RTLD_NOW);
    if (auto *err = ::dlerror())
        throw std::runtime_error(err);
    return std::shared_ptr<void>{h, &::dlclose};
}

template <class F>
F *DLProblem::load_func(std::string_view name) const {
    assert(handle);
    auto full_name = symbol_prefix + "_" + std::string(name);
    ::dlerror();
    auto *h = ::dlsym(handle.get(), full_name.c_str());
    if (auto *err = ::dlerror())
        throw std::runtime_error(err);
    // We can only hope that the user got the signature right ...
    return reinterpret_cast<F *>(h);
}

DLProblem::DLProblem(std::string so_filename, std::string symbol_prefix,
                     void *user_param)
    : BoxConstrProblem{0, 0}, so_filename(std::move(so_filename)),
      symbol_prefix(std::move(symbol_prefix)) {
    handle              = load_lib();
    auto *register_func = load_func<problem_register_t(void *)>("register");
    auto r              = register_func(user_param);
    // Avoid leaking if std::shared_ptr constructor throws
    std::unique_ptr<void, void (*)(void *)> unique_inst{r.instance, r.cleanup};
    std::unique_ptr<alpaqa_function_dict_t> unique_extra{r.extra_functions};
    // Store data returned by plugin
    instance        = std::shared_ptr<void>{std::move(unique_inst)};
    functions       = r.functions;
    this->n         = r.functions->n;
    this->m         = r.functions->m;
    this->C         = Box{this->n};
    this->D         = Box{this->m};
    extra_functions = std::shared_ptr<function_dict_t>{std::move(unique_extra)};
}

// clang-format off
auto DLProblem::eval_f(crvec x) const -> real_t { return functions->eval_f(instance.get(), x.data()); }
auto DLProblem::eval_grad_f(crvec x, rvec grad_fx) const -> void { return functions->eval_grad_f(instance.get(), x.data(), grad_fx.data()); }
auto DLProblem::eval_g(crvec x, rvec gx) const -> void { return functions->eval_g(instance.get(), x.data(), gx.data()); }
auto DLProblem::eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const -> void { return functions->eval_grad_g_prod(instance.get(), x.data(), y.data(), grad_gxy.data()); }
auto DLProblem::eval_grad_gi(crvec x, index_t i, rvec grad_gi) const -> void { return functions->eval_grad_gi(instance.get(), x.data(), i, grad_gi.data()); }
auto DLProblem::eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const -> void { return functions->eval_hess_L_prod(instance.get(), x.data(), y.data(), v.data(), Hv.data()); }
auto DLProblem::eval_hess_L(crvec x, crvec y, rmat H) const -> void { return functions->eval_hess_L(instance.get(), x.data(), y.data(), H.data()); }
auto DLProblem::eval_f_grad_f(crvec x, rvec grad_fx) const -> real_t { return functions->eval_f_grad_f(instance.get(), x.data(), grad_fx.data()); }
auto DLProblem::eval_f_g(crvec x, rvec g) const -> real_t { return functions->eval_f_g(instance.get(), x.data(), g.data()); }
auto DLProblem::eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const -> real_t { return functions->eval_f_grad_f_g(instance.get(), x.data(), grad_fx.data(), g.data()); }
auto DLProblem::eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const -> void { return functions->eval_grad_f_grad_g_prod(instance.get(), x.data(), y.data(), grad_f.data(), grad_gxy.data()); }
auto DLProblem::eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const -> void { return functions->eval_grad_L(instance.get(), x.data(), y.data(), grad_L.data(), work_n.data()); }
auto DLProblem::eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const -> real_t { return functions->eval_ψ(instance.get(), x.data(), y.data(), Σ.data(), ŷ.data()); }
auto DLProblem::eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n) const -> void { return functions->eval_grad_ψ_from_ŷ(instance.get(), x.data(), ŷ.data(), grad_ψ.data(), work_n.data()); }
auto DLProblem::eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const -> void { return functions->eval_grad_ψ(instance.get(), x.data(), y.data(), Σ.data(), grad_ψ.data(), work_n.data(), work_m.data()); }
auto DLProblem::eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const -> real_t { return functions->eval_ψ_grad_ψ(instance.get(), x.data(), y.data(), Σ.data(), grad_ψ.data(), work_n.data(), work_m.data()); }

bool DLProblem::provides_eval_f() const { return functions->eval_f != nullptr; }
bool DLProblem::provides_eval_grad_f() const { return functions->eval_grad_f != nullptr; }
bool DLProblem::provides_eval_g() const { return functions->eval_g != nullptr; }
bool DLProblem::provides_eval_grad_g_prod() const { return functions->eval_grad_g_prod != nullptr; }
bool DLProblem::provides_eval_grad_gi() const { return functions->eval_grad_gi != nullptr; }
bool DLProblem::provides_eval_hess_L_prod() const { return functions->eval_hess_L_prod != nullptr; }
bool DLProblem::provides_eval_hess_L() const { return functions->eval_hess_L != nullptr; }
bool DLProblem::provides_eval_f_grad_f() const { return functions->eval_f_grad_f != nullptr; }
bool DLProblem::provides_eval_f_g() const { return functions->eval_f_g != nullptr; }
bool DLProblem::provides_eval_f_grad_f_g() const { return functions->eval_f_grad_f_g != nullptr; }
bool DLProblem::provides_eval_grad_f_grad_g_prod() const { return functions->eval_grad_f_grad_g_prod != nullptr; }
bool DLProblem::provides_eval_grad_L() const { return functions->eval_grad_L != nullptr; }
bool DLProblem::provides_eval_ψ() const { return functions->eval_ψ != nullptr; }
bool DLProblem::provides_eval_grad_ψ_from_ŷ() const { return functions->eval_grad_ψ_from_ŷ != nullptr; }
bool DLProblem::provides_eval_grad_ψ() const { return functions->eval_grad_ψ != nullptr; }
bool DLProblem::provides_eval_ψ_grad_ψ() const { return functions->eval_ψ_grad_ψ != nullptr; }
// clang-format on

} // namespace alpaqa::dl