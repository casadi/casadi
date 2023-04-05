#include <alpaqa/dl/dl-problem.hpp>

#include <dlfcn.h>
#include <cassert>
#include <memory>
#include <stdexcept>

namespace alpaqa::dl {

std::shared_ptr<void> DLLoader::load_lib() const {
    assert(!so_filename.empty());
    ::dlerror();
    void *h = ::dlopen(so_filename.c_str(), RTLD_LOCAL | RTLD_NOW);
    if (auto *err = ::dlerror())
        throw std::runtime_error(err);
    return std::shared_ptr<void>{h, &::dlclose};
}

template <class F>
F *DLLoader::load_func(std::string_view name) const {
    assert(handle);
    auto full_name = symbol_prefix + "_" + std::string(name);
    ::dlerror();
    auto *h = ::dlsym(handle.get(), full_name.c_str());
    if (auto *err = ::dlerror())
        throw std::runtime_error(err);
    // We can only hope that the user got the signature right ...
    return reinterpret_cast<F *>(h);
}

DLLoader::DLLoader(std::string so_filename, std::string symbol_prefix)
    : so_filename(std::move(so_filename)),
      symbol_prefix(std::move(symbol_prefix)), handle(load_lib()) {}

DLProblem::DLProblem(std::string so_filename, std::string symbol_prefix,
                     void *user_param)
    : DLLoader{std::move(so_filename), std::move(symbol_prefix)},
      BoxConstrProblem{0, 0} {
    auto *register_func = load_func<problem_register_t(void *)>("register");
    auto r              = register_func(user_param);
    // Avoid leaking if std::shared_ptr constructor throws
    std::unique_ptr<void, void (*)(void *)> unique_inst{r.instance, r.cleanup};
    std::unique_ptr<alpaqa_function_dict_t> unique_extra{r.extra_functions};
    // Store data returned by plugin
    instance  = std::shared_ptr<void>{std::move(unique_inst)};
    functions = r.functions;
    this->n   = r.functions->n;
    this->m   = r.functions->m;
    this->C   = Box{this->n};
    this->D   = Box{this->m};
    if (functions->get_C)
        functions->get_C(instance.get(), this->C.lowerbound.data(),
                         this->C.upperbound.data());
    if (functions->get_D)
        functions->get_D(instance.get(), this->D.lowerbound.data(),
                         this->D.upperbound.data());
    extra_functions = std::shared_ptr<function_dict_t>{std::move(unique_extra)};
}

auto DLProblem::eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂,
                                    rvec p) const -> real_t {
    if (functions->eval_prox_grad_step)
        return functions->eval_prox_grad_step(
            instance.get(), γ, x.data(), grad_ψ.data(), x̂.data(), p.data());
    return BoxConstrProblem<config_t>::eval_prox_grad_step(γ, x, grad_ψ, x̂, p);
}

// clang-format off
auto DLProblem::eval_f(crvec x) const -> real_t { return functions->eval_f(instance.get(), x.data()); }
auto DLProblem::eval_grad_f(crvec x, rvec grad_fx) const -> void { return functions->eval_grad_f(instance.get(), x.data(), grad_fx.data()); }
auto DLProblem::eval_g(crvec x, rvec gx) const -> void { return functions->eval_g(instance.get(), x.data(), gx.data()); }
auto DLProblem::eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const -> void { return functions->eval_grad_g_prod(instance.get(), x.data(), y.data(), grad_gxy.data()); }
auto DLProblem::eval_grad_gi(crvec x, index_t i, rvec grad_gi) const -> void { return functions->eval_grad_gi(instance.get(), x.data(), i, grad_gi.data()); }
auto DLProblem::eval_jac_g(crvec x, rindexvec inner_idx, rindexvec outer_ptr, rvec J_values) const -> void { return functions->eval_jac_g(instance.get(), x.data(), inner_idx.data(), outer_ptr.data(), J_values.size() == 0 ? nullptr : J_values.data()); }
auto DLProblem::get_jac_g_num_nonzeros() const -> length_t { return functions->get_jac_g_num_nonzeros(instance.get()); }
auto DLProblem::eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v, rvec Hv) const -> void { return functions->eval_hess_L_prod(instance.get(), x.data(), y.data(), scale, v.data(), Hv.data()); }
auto DLProblem::eval_hess_L(crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const -> void { return functions->eval_hess_L(instance.get(), x.data(), y.data(), scale, inner_idx.data(), outer_ptr.data(), H_values.size() == 0 ? nullptr : H_values.data()); }
auto DLProblem::get_hess_L_num_nonzeros() const -> length_t { return functions->get_hess_L_num_nonzeros(instance.get()); }
auto DLProblem::eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv) const -> void { return functions->eval_hess_ψ_prod(instance.get(), x.data(), y.data(), Σ.data(), scale, D.lowerbound.data(), D.upperbound.data(), v.data(), Hv.data()); }
auto DLProblem::eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const -> void { return functions->eval_hess_ψ(instance.get(), x.data(), y.data(), Σ.data(), scale, D.lowerbound.data(), D.upperbound.data(), inner_idx.data(), outer_ptr.data(), H_values.size() == 0 ? nullptr : H_values.data()); }
auto DLProblem::get_hess_ψ_num_nonzeros() const -> length_t { return functions->get_hess_ψ_num_nonzeros(instance.get()); }
auto DLProblem::eval_f_grad_f(crvec x, rvec grad_fx) const -> real_t { return functions->eval_f_grad_f(instance.get(), x.data(), grad_fx.data()); }
auto DLProblem::eval_f_g(crvec x, rvec g) const -> real_t { return functions->eval_f_g(instance.get(), x.data(), g.data()); }
auto DLProblem::eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const -> void { return functions->eval_grad_f_grad_g_prod(instance.get(), x.data(), y.data(), grad_f.data(), grad_gxy.data()); }
auto DLProblem::eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const -> void { return functions->eval_grad_L(instance.get(), x.data(), y.data(), grad_L.data(), work_n.data()); }
auto DLProblem::eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const -> real_t { return functions->eval_ψ(instance.get(), x.data(), y.data(), Σ.data(), D.lowerbound.data(), D.upperbound.data(), ŷ.data()); }
auto DLProblem::eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const -> void { return functions->eval_grad_ψ(instance.get(), x.data(), y.data(), Σ.data(), D.lowerbound.data(), D.upperbound.data(), grad_ψ.data(), work_n.data(), work_m.data()); }
auto DLProblem::eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const -> real_t { return functions->eval_ψ_grad_ψ(instance.get(), x.data(), y.data(), Σ.data(), D.lowerbound.data(), D.upperbound.data(), grad_ψ.data(), work_n.data(), work_m.data()); }

bool DLProblem::provides_eval_f() const { return functions->eval_f != nullptr; }
bool DLProblem::provides_eval_grad_f() const { return functions->eval_grad_f != nullptr; }
bool DLProblem::provides_eval_g() const { return functions->eval_g != nullptr; }
bool DLProblem::provides_eval_grad_g_prod() const { return functions->eval_grad_g_prod != nullptr; }
bool DLProblem::provides_eval_jac_g() const { return functions->eval_jac_g != nullptr; }
bool DLProblem::provides_get_jac_g_num_nonzeros() const { return functions->get_jac_g_num_nonzeros != nullptr; }
bool DLProblem::provides_eval_grad_gi() const { return functions->eval_grad_gi != nullptr; }
bool DLProblem::provides_eval_hess_L_prod() const { return functions->eval_hess_L_prod != nullptr; }
bool DLProblem::provides_eval_hess_L() const { return functions->eval_hess_L != nullptr; }
bool DLProblem::provides_get_hess_L_num_nonzeros() const { return functions->get_hess_L_num_nonzeros != nullptr; }
bool DLProblem::provides_eval_hess_ψ_prod() const { return functions->eval_hess_ψ_prod != nullptr; }
bool DLProblem::provides_eval_hess_ψ() const { return functions->eval_hess_ψ != nullptr; }
bool DLProblem::provides_get_hess_ψ_num_nonzeros() const { return functions->get_hess_ψ_num_nonzeros != nullptr; }
bool DLProblem::provides_eval_f_grad_f() const { return functions->eval_f_grad_f != nullptr; }
bool DLProblem::provides_eval_f_g() const { return functions->eval_f_g != nullptr; }
bool DLProblem::provides_eval_grad_f_grad_g_prod() const { return functions->eval_grad_f_grad_g_prod != nullptr; }
bool DLProblem::provides_eval_grad_L() const { return functions->eval_grad_L != nullptr; }
bool DLProblem::provides_eval_ψ() const { return functions->eval_ψ != nullptr; }
bool DLProblem::provides_eval_grad_ψ() const { return functions->eval_grad_ψ != nullptr; }
bool DLProblem::provides_eval_ψ_grad_ψ() const { return functions->eval_ψ_grad_ψ != nullptr; }
bool DLProblem::provides_get_box_C() const { return functions->eval_prox_grad_step == nullptr; }
// clang-format on

#if ALPAQA_WITH_OCP

DLControlProblem::DLControlProblem(std::string so_filename,
                                   std::string symbol_prefix, void *user_param)
    : DLLoader{std::move(so_filename), std::move(symbol_prefix)} {
    auto *register_func =
        load_func<control_problem_register_t(void *)>("register");
    auto r = register_func(user_param);
    // Avoid leaking if std::shared_ptr constructor throws
    std::unique_ptr<void, void (*)(void *)> unique_inst{r.instance, r.cleanup};
    std::unique_ptr<alpaqa_function_dict_t> unique_extra{r.extra_functions};
    // Store data returned by plugin
    instance        = std::shared_ptr<void>{std::move(unique_inst)};
    functions       = r.functions;
    extra_functions = std::shared_ptr<function_dict_t>{std::move(unique_extra)};
}

// clang-format off
auto DLControlProblem::get_U(Box &U) const -> void { return functions->get_U(instance.get(), U.lowerbound.data(), U.upperbound.data()); }
auto DLControlProblem::get_D(Box &D) const -> void { return functions->get_D(instance.get(), D.lowerbound.data(), D.upperbound.data()); }
auto DLControlProblem::get_D_N(Box &D) const -> void { return functions->get_D_N(instance.get(), D.lowerbound.data(), D.upperbound.data()); }
auto DLControlProblem::get_x_init(rvec x_init) const -> void { return functions->get_x_init(instance.get(), x_init.data()); }
auto DLControlProblem::eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const -> void { return functions->eval_f(instance.get(), timestep, x.data(), u.data(), fxu.data()); }
auto DLControlProblem::eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const -> void { return functions->eval_jac_f(instance.get(), timestep, x.data(), u.data(), J_fxu.data()); }
auto DLControlProblem::eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const -> void { return functions->eval_grad_f_prod(instance.get(), timestep, x.data(), u.data(), p.data(), grad_fxu_p.data()); }
auto DLControlProblem::eval_h(index_t timestep, crvec x, crvec u, rvec h) const -> void { return functions->eval_h(instance.get(), timestep, x.data(), u.data(), h.data()); }
auto DLControlProblem::eval_h_N(crvec x, rvec h) const -> void { return functions->eval_h_N(instance.get(), x.data(), h.data()); }
auto DLControlProblem::eval_l(index_t timestep, crvec h) const -> real_t { return functions->eval_l(instance.get(), timestep, h.data()); }
auto DLControlProblem::eval_l_N(crvec h) const -> real_t { return functions->eval_l_N(instance.get(), h.data()); }
auto DLControlProblem::eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const -> void { return functions->eval_qr(instance.get(), timestep, xu.data(), h.data(), qr.data()); }
auto DLControlProblem::eval_q_N(crvec x, crvec h, rvec q) const -> void { return functions->eval_q_N(instance.get(), x.data(), h.data(), q.data()); }
auto DLControlProblem::eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const -> void { return functions->eval_add_Q(instance.get(), timestep, xu.data(), h.data(), Q.data()); }
auto DLControlProblem::eval_add_Q_N(crvec x, crvec h, rmat Q) const -> void { return functions->eval_add_Q_N(instance.get(), x.data(), h.data(), Q.data()); }
auto DLControlProblem::eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R, rvec work) const -> void { return functions->eval_add_R_masked(instance.get(), timestep, xu.data(), h.data(), mask.data(), R.data(), work.data()); }
auto DLControlProblem::eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S, rvec work) const -> void { return functions->eval_add_S_masked(instance.get(), timestep, xu.data(), h.data(), mask.data(), S.data(), work.data()); }
auto DLControlProblem::eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J, crindexvec mask_K, crvec v, rvec out, rvec work) const -> void { return functions->eval_add_R_prod_masked(instance.get(), timestep, xu.data(), h.data(), mask_J.data(), mask_K.data(), v.data(), out.data(), work.data()); }
auto DLControlProblem::eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v, rvec out, rvec work) const -> void { return functions->eval_add_S_prod_masked(instance.get(), timestep, xu.data(), h.data(), mask_K.data(), v.data(), out.data(), work.data()); }
auto DLControlProblem::get_R_work_size() const -> length_t { return functions->get_R_work_size(instance.get()); }
auto DLControlProblem::get_S_work_size() const -> length_t { return functions->get_S_work_size(instance.get()); }
auto DLControlProblem::eval_constr(index_t timestep, crvec x, rvec c) const -> void { return functions->eval_constr(instance.get(), timestep, x.data(), c.data()); }
auto DLControlProblem::eval_constr_N(crvec x, rvec c) const -> void { return functions->eval_constr_N(instance.get(), x.data(), c.data()); }
auto DLControlProblem::eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const -> void { return functions->eval_grad_constr_prod(instance.get(), timestep, x.data(), p.data(), grad_cx_p.data()); }
auto DLControlProblem::eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const -> void { return functions->eval_grad_constr_prod_N(instance.get(), x.data(), p.data(), grad_cx_p.data()); }
auto DLControlProblem::eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const -> void { return functions->eval_add_gn_hess_constr(instance.get(), timestep, x.data(), M.data(), out.data()); }
auto DLControlProblem::eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const -> void { return functions->eval_add_gn_hess_constr_N(instance.get(), x.data(), M.data(), out.data()); }

bool DLControlProblem::provides_get_D() const { return functions->get_D != nullptr; }
bool DLControlProblem::provides_get_D_N() const { return functions->get_D_N != nullptr; }
bool DLControlProblem::provides_eval_add_Q_N() const { return functions->eval_add_Q_N != nullptr; }
bool DLControlProblem::provides_eval_add_R_prod_masked() const { return functions->eval_add_R_prod_masked != nullptr; }
bool DLControlProblem::provides_eval_add_S_prod_masked() const { return functions->eval_add_S_prod_masked != nullptr; }
bool DLControlProblem::provides_get_R_work_size() const { return functions->get_R_work_size != nullptr; }
bool DLControlProblem::provides_get_S_work_size() const { return functions->get_S_work_size != nullptr; }
bool DLControlProblem::provides_eval_constr() const { return functions->eval_constr != nullptr; }
bool DLControlProblem::provides_eval_constr_N() const { return functions->eval_constr_N != nullptr; }
bool DLControlProblem::provides_eval_grad_constr_prod() const { return functions->eval_grad_constr_prod != nullptr; }
bool DLControlProblem::provides_eval_grad_constr_prod_N() const { return functions->eval_grad_constr_prod_N != nullptr; }
bool DLControlProblem::provides_eval_add_gn_hess_constr() const { return functions->eval_add_gn_hess_constr != nullptr; }
bool DLControlProblem::provides_eval_add_gn_hess_constr_N() const { return functions->eval_add_gn_hess_constr_N != nullptr; }
// clang-format on

#endif

} // namespace alpaqa::dl