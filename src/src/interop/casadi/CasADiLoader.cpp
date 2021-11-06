#include <alpaqa/interop/casadi/CasADiFunctionWrapper.hpp>
#include <alpaqa/interop/casadi/CasADiLoader.hpp>

#include <casadi/core/external.hpp>

#include <memory>
#include <optional>
#include <stdexcept>

namespace alpaqa {

#if 0
std::function<alpaqa::Problem::f_sig>
load_CasADi_objective(const std::string &so_name, const std::string &fun_name) {
    return CasADiFun_1Vi1So(casadi::external(fun_name, so_name));
}
std::function<alpaqa::Problem::grad_f_sig>
load_CasADi_gradient_objective(const std::string &so_name,
                               const std::string &fun_name) {
    return CasADiFun_1Vi1Vo(casadi::external(fun_name, so_name));
}
std::function<alpaqa::Problem::g_sig>
load_CasADi_constraints(const std::string &so_name,
                        const std::string &fun_name) {
    return CasADiFun_1Vi1Vo(casadi::external(fun_name, so_name));
}
std::function<alpaqa::Problem::grad_g_prod_sig>
load_CasADi_gradient_constraints_prod(const std::string &so_name,
                                      const std::string &fun_name) {
    return [csf{CasADiFun_2Vi1Vo(casadi::external(fun_name, so_name))}] //
        (alpaqa::crvec x, alpaqa::crvec y, alpaqa::rvec gradprod) {                 //
            if (y.size() == 0)
                gradprod.setZero();
            else
                csf(x, y, gradprod);
        };
}
std::function<alpaqa::Problem::hess_L_sig>
load_CasADi_hessian_lagrangian(const std::string &so_name,
                               const std::string &fun_name) {
    return [csf{CasADiFun_2Vi1Mo(casadi::external(fun_name, so_name))}] //
        (alpaqa::crvec x, alpaqa::crvec y, alpaqa::rmat H) {                        //
            // Fix the stride if the matrix is larger than n
            if (x.rows() != H.rows()) { // TODO: this is probably unnecessary
                for (auto c = x.rows(); c-- > 1;)
                    for (auto r = x.rows(); r-- > 0;)
                        std::swap(H(r, c), H.data()[r + x.rows() * c]);
            }
            csf(x, y, H);
            // Fix the stride if the matrix is larger than n
            if (x.rows() != H.rows()) {
                for (auto c = x.rows(); c-- > 1;)
                    for (auto r = x.rows(); r-- > 0;)
                        std::swap(H(r, c), H.data()[r + x.rows() * c]);
            }
        };
}
std::function<alpaqa::Problem::hess_L_prod_sig>
load_CasADi_hessian_lagrangian_prod(const std::string &so_name,
                                    const std::string &fun_name) {
    return CasADiFun_3Vi1Vo(casadi::external(fun_name, so_name));
}
#endif

template <class F>
auto wrap_load(const std::string &so_name, const char *name, F f) {
    try {
        return f();
    } catch (const std::invalid_argument &e) {
        throw std::invalid_argument("Unable to load function '" + so_name +
                                    ":" + name + "': " + e.what());
    }
}

template <class T, class... Args>
auto wrapped_load(const std::string &so_name, const char *name,
                  Args &&...args) {
    return wrap_load(so_name, name, [&] {
        return T(casadi::external(name, so_name), std::forward<Args>(args)...);
    });
}

constexpr static auto dims = [](auto... a) {
    return std::array<casadi_int, sizeof...(a)>{a...};
};
using dim = std::pair<casadi_int, casadi_int>;

alpaqa::Problem load_CasADi_problem(const std::string &so_name, unsigned n,
                                unsigned m, bool second_order) {

    auto load_g_unknown_dims = [&] {
        CasADiFunctionEvaluator<1, 1> g{casadi::external("g", so_name)};
        if (g.fun.size2_in(0) != 1)
            throw std::invalid_argument(
                "First input argument should be a column vector.");
        if (g.fun.size2_out(0) != 1)
            throw std::invalid_argument(
                "First output argument should be a column vector.");
        if (n == 0)
            n = g.fun.size1_in(0);
        if (m == 0)
            m = g.fun.size1_out(0);
        g.validate_dimensions({dim(n, 1)}, {dim(m, 1)});
        return g;
    };

    auto load_g_known_dims = [&] {
        CasADiFunctionEvaluator<1, 1> g{
            casadi::external("g", so_name), {dim(n, 1)}, {dim(m, 1)}};
        return g;
    };

    CasADiFunctionEvaluator<1, 1> g =
        (n == 0 || m == 0)
            // If not all dimensions are specified, load the function "g" to
            //determine the missing dimensions.
            ? wrap_load(so_name, "g", load_g_unknown_dims)
            // Otherwise, load the function "g" and compare its dimensions to
            // the dimensions specified by the user.
            : wrap_load(so_name, "g", load_g_known_dims);

    auto prob = alpaqa::Problem(n, m);

    prob.f      = wrapped_load<CasADiFun_1Vi1So>(so_name, "f", n);
    prob.grad_f = wrapped_load<CasADiFun_1Vi1Vo>(so_name, "grad_f", n, n);
    prob.g      = CasADiFun_1Vi1Vo(std::move(g));
    auto grad_g =
        wrapped_load<CasADiFun_2Vi1Vo>(so_name, "grad_g", dims(n, m), n);
    prob.grad_g_prod = grad_g;
    if (second_order) {
        alpaqa::vec w    = alpaqa::vec::Zero(m);
        prob.grad_gi = //
            [grad_g, w](alpaqa::crvec x, unsigned i, alpaqa::rvec g) mutable {
                w(i) = 1;
                grad_g(x, w, g);
                w(i) = 0;
            };
        prob.hess_L =                                              //
            wrapped_load<CasADiFun_2Vi1Mo>(so_name, "hess_L",      //
                                           dims(n, m), dim(n, n)); //
        prob.hess_L_prod =                                         //
            wrapped_load<CasADiFun_3Vi1Vo>(so_name, "hess_L_prod", //
                                           dims(n, m, n), n);      //
    }
    return prob;
}

class CasADiParamWrapper
    : public ParamWrapper,
      public std::enable_shared_from_this<CasADiParamWrapper> {

  public:
    struct Functions {
        CasADiFun_2Vi1So f;
        CasADiFun_2Vi1Vo grad_f;
        CasADiFun_2Vi1Vo g;
        CasADiFun_3Vi1Vo grad_g_prod;
        std::optional<CasADiFun_3Vi1Mo> hess_L;
        std::optional<CasADiFun_4Vi1Vo> hess_L_prod;
    } cs;

  private:
    CasADiParamWrapper(unsigned p, Functions &&functions)
        : ParamWrapper(p), cs(std::move(functions)) {}

  public:
    static std::shared_ptr<CasADiParamWrapper> create(unsigned p,
                                                      Functions &&functions) {
        return std::make_shared<CasADiParamWrapper>(CasADiParamWrapper{
            p,
            std::move(functions),
        });
    }

    void wrap(Problem &prob) override {
        auto param = this->shared_from_this();
        prob.f     = [param](crvec x) -> real_t {
            return param->cs.f(x, param->param);
        };
        prob.grad_f = [param](crvec x, rvec gr) {
            param->cs.grad_f(x, param->param, gr);
        };
        prob.g = [param](crvec x, rvec g) -> void {
            param->cs.g(x, param->param, g);
        };
        prob.grad_g_prod = [param](crvec x, crvec y, rvec g) {
            param->cs.grad_g_prod(x, param->param, y, g);
        };
        alpaqa::vec w    = alpaqa::vec::Zero(prob.m);
        prob.grad_gi = [param, w](crvec x, unsigned i, rvec g) mutable {
            w(i) = 1;
            param->cs.grad_g_prod(x, param->param, w, g);
            w(i) = 0;
        };
        if (param->cs.hess_L) {
            prob.hess_L = [param](crvec x, crvec y, rvec g) {
                (*param->cs.hess_L)(x, param->param, y, g);
            };
        }
        if (param->cs.hess_L_prod) {
            prob.hess_L_prod = [param](crvec x, crvec y, crvec v, rvec g) {
                (*param->cs.hess_L_prod)(x, param->param, y, v, g);
            };
        }
    }

    std::shared_ptr<ParamWrapper> clone() const override {
        return std::make_shared<CasADiParamWrapper>(*this);
    }
};

ProblemWithParam load_CasADi_problem_with_param(const std::string &so_name,
                                                unsigned n, unsigned m,
                                                unsigned p, bool second_order) {

    auto load_g_unknown_dims = [&] {
        CasADiFunctionEvaluator<2, 1> g{casadi::external("g", so_name)};
        if (g.fun.size2_in(0) != 1)
            throw std::invalid_argument(
                "First input argument should be a column vector.");
        if (g.fun.size2_in(1) != 1)
            throw std::invalid_argument(
                "Second input argument should be a column vector.");
        if (g.fun.size2_out(0) != 1)
            throw std::invalid_argument(
                "First output argument should be a column vector.");
        if (n == 0)
            n = g.fun.size1_in(0);
        if (m == 0)
            m = g.fun.size1_out(0);
        if (p == 0)
            p = g.fun.size1_in(1);
        g.validate_dimensions({dim(n, 1), dim(p, 1)}, {dim(m, 1)});
        return g;
    };

    auto load_g_known_dims = [&] {
        CasADiFunctionEvaluator<2, 1> g{casadi::external("g", so_name),
                                        {dim(n, 1), dim(p, 1)},
                                        {dim(m, 1)}};
        return g;
    };

    CasADiFunctionEvaluator<2, 1> g =
        (n == 0 || m == 0 || p == 0)
            // If not all dimensions are specified, load the function "g" to
            // determine the missing dimensions.
            ? wrap_load(so_name, "g", load_g_unknown_dims)
            // Otherwise, load the function "g" and compare its dimensions to
            // the dimensions specified by the user.
            : wrap_load(so_name, "g", load_g_known_dims);

    auto prob = ProblemWithParam(n, m);

    prob.wrapper = CasADiParamWrapper::create(
        p,
        {
            wrapped_load<CasADiFun_2Vi1So>(so_name, "f", dims(n, p)),
            wrapped_load<CasADiFun_2Vi1Vo>(so_name, "grad_f", dims(n, p), n),
            CasADiFun_2Vi1Vo(std::move(g)),
            wrapped_load<CasADiFun_3Vi1Vo>(so_name, "grad_g", dims(n, p, m), n),
            second_order ? std::make_optional(wrapped_load<CasADiFun_3Vi1Mo>(
                               so_name, "hess_L", dims(n, p, m), dim(n, n)))
                         : std::nullopt,
            second_order ? std::make_optional(wrapped_load<CasADiFun_4Vi1Vo>(
                               so_name, "hess_L_prod", dims(n, p, m, n), n))
                         : std::nullopt,
        });
    prob.wrapper->wrap(prob);
    return prob;
}

} // namespace alpaqa