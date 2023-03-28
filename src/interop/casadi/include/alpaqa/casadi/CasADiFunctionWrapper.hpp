#pragma once

#include <alpaqa/config/config.hpp>

#include <stdexcept>
#include <string>
#include <vector>

#include <casadi/core/function.hpp>
#include <casadi/mem.h>

namespace alpaqa::casadi_loader {

/// Class for evaluating CasADi functions, allocating the necessary workspace
/// storage in advance for allocation-free evaluations.
template <Config Conf, size_t N_in, size_t N_out>
class CasADiFunctionEvaluator {
  public:
    USING_ALPAQA_CONFIG(Conf);
    static_assert(std::is_same_v<real_t, casadi_real>);

    using casadi_dim = std::pair<casadi_int, casadi_int>;

    /// @throws std::invalid_argument
    CasADiFunctionEvaluator(casadi::Function &&f)
        : fun(std::move(f)), iwork(fun.sz_iw()), dwork(fun.sz_w()),
          arg_work(fun.sz_arg()), res_work(fun.sz_res()) {
        using namespace std::literals::string_literals;
        if (N_in != fun.n_in())
            throw std::invalid_argument(
                "Invalid number of input arguments: got "s +
                std::to_string(fun.n_in()) + ", should be " +
                std::to_string(N_in) + ".");
        if (N_out != fun.n_out())
            throw std::invalid_argument(
                "Invalid number of output arguments: got "s +
                std::to_string(fun.n_out()) + ", should be " +
                std::to_string(N_out) + ".");
    }

    /// @throws std::invalid_argument
    CasADiFunctionEvaluator(casadi::Function &&f,
                            const std::array<casadi_dim, N_in> &dim_in,
                            const std::array<casadi_dim, N_out> &dim_out)
        : CasADiFunctionEvaluator{std::move(f)} {
        validate_dimensions(dim_in, dim_out);
    }

    /// @throws std::invalid_argument
    void
    validate_dimensions(const std::array<casadi_dim, N_in> &dim_in   = {},
                        const std::array<casadi_dim, N_out> &dim_out = {}) {
        using namespace std::literals::string_literals;
        static constexpr std::array count{"first",   "second", "third",
                                          "fourth",  "fifth",  "sixth",
                                          "seventh", "eighth"};
        static_assert(N_in <= count.size());
        static_assert(N_out <= count.size());
        auto to_string = [](casadi_dim d) {
            return "(" + std::to_string(d.first) + ", " +
                   std::to_string(d.second) + ")";
        };
        for (size_t n = 0; n < N_in; ++n) {
            auto cs_n = static_cast<casadi_int>(n);
            if (dim_in[n].first != 0 && dim_in[n] != fun.size_in(cs_n))
                throw std::invalid_argument(
                    "Invalid dimension of "s + count[n] +
                    " input argument: got " + to_string(fun.size_in(cs_n)) +
                    ", should be " + to_string(dim_in[n]) + ".");
        }
        for (size_t n = 0; n < N_out; ++n) {
            auto cs_n = static_cast<casadi_int>(n);
            if (dim_out[n].first != 0 && dim_out[n] != fun.size_out(cs_n))
                throw std::invalid_argument(
                    "Invalid dimension of "s + count[n] +
                    " output argument: got " + to_string(fun.size_out(cs_n)) +
                    ", should be " + to_string(dim_out[n]) + ".");
        }
    }

  protected:
    void operator()(const double *const *in, double *const *out) const {
        std::copy_n(in, N_in, arg_work.begin());
        std::copy_n(out, N_out, res_work.begin());
        fun(arg_work.data(), res_work.data(), iwork.data(), dwork.data(), 0);
    }

  public:
    void operator()(const double *const (&in)[N_in],
                    double *const (&out)[N_out]) const {
        this->operator()(&in[0], &out[0]);
    }

  public:
    casadi::Function fun;

  private:
    mutable std::vector<casadi_int> iwork;
    mutable std::vector<double> dwork;
    mutable std::vector<const double *> arg_work;
    mutable std::vector<double *> res_work;
};

} // namespace alpaqa::casadi_loader