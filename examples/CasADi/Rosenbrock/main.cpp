#include <casadi/core/function.hpp>
#include <panoc-alm/alm.hpp>

/// Class for evaluating CasADi functions, allocating the necessary workspace
/// storage in advance for allocation-free evaluations.
template <size_t N_in, size_t N_out>
class CasADiFunctionEvaluator {
  public:
    CasADiFunctionEvaluator(casadi::Function &&f)
        : fun(std::move(f)), iwork(fun.sz_iw()), dwork(fun.sz_w()) {
        assert(N_in == fun.n_in());
        assert(N_out == fun.n_out());
    }

  protected:
    void operator()(const double *const *in, double *const *out) const {
        fun(const_cast<const double **>(in), const_cast<double **>(out),
            iwork.data(), dwork.data(), 0);
    }

  public:
    void operator()(const double *const (&in)[N_in],
                    double *const (&out)[N_out]) const {
        this->operator()(&in[0], &out[0]);
    }

  private:
    casadi::Function fun;
    mutable std::vector<casadi_int> iwork;
    mutable std::vector<double> dwork;
};

/// Wrapper for CasADiFunctionEvaluator with 1 vector input, scalar output.
class CasADiFun_1iso {
  public:
    CasADiFun_1iso(casadi::Function &&f) : fun(std::move(f)) {}

    double operator()(const pa::vec &x) const {
        double out;
        fun({x.data()}, {&out});
        return out;
    }

  private:
    CasADiFunctionEvaluator<1, 1> fun;
};

/// Wrapper for CasADiFunctionEvaluator with 1 vector input, 1 vector output.
class CasADiFun_1i1o {
  public:
    CasADiFun_1i1o(casadi::Function &&f) : fun(std::move(f)) {}

    void operator()(const pa::vec &in, pa::vec &out) const {
        fun({in.data()}, {out.data()});
    }

  private:
    CasADiFunctionEvaluator<1, 1> fun;
};

/// Wrapper for CasADiFunctionEvaluator with 2 vector inputs, 1 vector output.
class CasADiFun_2i1o {
  public:
    CasADiFun_2i1o(casadi::Function &&f) : fun(std::move(f)) {}

    void operator()(const pa::vec &in1, const pa::vec &in2,
                    pa::vec &out) const {
        fun({in1.data(), in2.data()}, {out.data()});
    }

  private:
    CasADiFunctionEvaluator<2, 1> fun;
};

#include <casadi/core/external.hpp>
namespace cs = casadi;

int main(int argc, char *argv[]) {
    const char *so_name =
        "./examples/CasADi/Rosenbrock/librosenbrock_functions.so";
    if (argc > 1)
        so_name = argv[1];

    constexpr auto inf = std::numeric_limits<double>::infinity();
    using pa::real_t;
    using pa::vec;

    pa::Problem p;
    p.n = 3;
    p.m = 1;
    p.C = pa::Box{
        vec::Constant(3, inf),
        vec::Constant(3, -inf),
    };
    p.D = pa::Box{
        vec::Constant(1, 0.),
        vec::Constant(1, 0.),
    };
    p.f      = CasADiFun_1iso(cs::external("f", so_name));
    p.grad_f = CasADiFun_1i1o(cs::external("grad_f", so_name));
    p.g      = CasADiFun_1i1o(cs::external("g", so_name));
    p.grad_g = CasADiFun_2i1o(cs::external("grad_g_v", so_name));

    pa::ALMParams almparam;
    almparam.ε        = 1e-8;
    almparam.δ        = 1e-8;
    almparam.Δ        = 10;
    almparam.Σ₀       = 2;
    almparam.ε₀       = 1;
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.σₘₐₓ     = 1e9;
    almparam.max_iter = 10;

    pa::PANOCParams panocparam;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.lbfgs_mem   = 10;
    panocparam.max_iter    = 50;

    pa::ALMSolver solver{almparam, panocparam};

    vec x(3);
    x << 2.5, 3.0, 0.75;
    vec y(1);
    y << 1;

    auto stats = solver(p, y, x);

    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;

    std::cout << "inner: " << stats.inner_iterations << std::endl;
    std::cout << "outer: " << stats.outer_iterations << std::endl;
}