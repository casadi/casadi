#pragma once

#include <casadi/core/function.hpp>

#include <alpaqa/util/vec.hpp>

#include <stdexcept>
#include <string>
#include <vector>

/// @addtogroup grp_ExternalProblemLoaders
/// @{

/// Class for evaluating CasADi functions, allocating the necessary workspace
/// storage in advance for allocation-free evaluations.
template <size_t N_in, size_t N_out>
class CasADiFunctionEvaluator {
  public:
    using casadi_dim = std::pair<casadi_int, casadi_int>;

    /// @throws std::invalid_argument
    CasADiFunctionEvaluator(casadi::Function &&f,
                            const casadi_dim (&dim_in)[N_in]   = {},
                            const casadi_dim (&dim_out)[N_out] = {})
        : fun(std::move(f)), iwork(fun.sz_iw()), dwork(fun.sz_w()) {
        if (N_in != fun.n_in())
            throw std::invalid_argument("Invalid number of input arguments.");
        if (N_out != fun.n_out())
            throw std::invalid_argument("Invalid number of output arguments.");
        validate_dimensions(dim_in, dim_out);
    }

    /// @throws std::invalid_argument
    void validate_dimensions(const casadi_dim (&dim_in)[N_in],
                             const casadi_dim (&dim_out)[N_out]) {
        using std::operator""s;
        constexpr static const char *count[]{"first", "second", "third",
                                             "fourth"};
        static_assert(N_in <= 4);
        static_assert(N_out <= 4);
        auto to_string = [](casadi_dim d) {
            return "(" + std::to_string(d.first) + ", " +
                   std::to_string(d.second) + ")";
        };
        for (size_t n = 0; n < N_in; ++n)
            if (dim_in[n].first != 0 && dim_in[n] != fun.size_in(n))
                throw std::invalid_argument(
                    "Invalid dimension of "s + count[n] +
                    " input argument: got " + to_string(fun.size_in(n)) +
                    ", should be " + to_string(dim_in[n]) + ".");
        for (size_t n = 0; n < N_out; ++n)
            if (dim_out[n].first != 0 && dim_out[n] != fun.size_out(n))
                throw std::invalid_argument(
                    "Invalid dimension of "s + count[n] +
                    " output argument: got " + to_string(fun.size_out(n)) +
                    ", should be " + to_string(dim_out[n]) + ".");
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

  public:
    casadi::Function fun;

  private:
    mutable std::vector<casadi_int> iwork;
    mutable std::vector<double> dwork;
};

/// Wrapper for CasADiFunctionEvaluator with 1 vector input, scalar output.
class CasADiFun_1Vi1So {
  public:
    CasADiFun_1Vi1So(casadi::Function &&f, casadi_int dim_in = 0)
        : fun(std::move(f), {{dim_in, 1}}, {{1, 1}}) {}

    double operator()(alpaqa::crvec x) const {
        double out;
        fun({x.data()}, {&out});
        return out;
    }

  private:
    CasADiFunctionEvaluator<1, 1> fun;
};

/// Wrapper for CasADiFunctionEvaluator with 2 vector inputs, scalar output.
class CasADiFun_2Vi1So {
  public:
    CasADiFun_2Vi1So(casadi::Function &&f,
                     const std::array<casadi_int, 2> &dim_in = {})
        : fun(std::move(f), {{dim_in[0], 1}, {dim_in[1], 1}}, {{1, 1}}) {}

    double operator()(alpaqa::crvec x, alpaqa::crvec p) const {
        double out;
        fun({x.data(), p.data()}, {&out});
        return out;
    }

  private:
    CasADiFunctionEvaluator<2, 1> fun;
};

/// Wrapper for CasADiFunctionEvaluator with 1 vector input, 1 vector output.
class CasADiFun_1Vi1Vo {
  public:
    CasADiFun_1Vi1Vo(CasADiFunctionEvaluator<1, 1> &&fun)
        : fun(std::move(fun)) {}
    CasADiFun_1Vi1Vo(casadi::Function &&f, casadi_int dim_in = 0,
                     casadi_int dim_out = 0)
        : fun(std::move(f), {{dim_in, 1}}, {{dim_out, 1}}) {}

    void operator()(alpaqa::crvec in, alpaqa::rvec out) const {
        fun({in.data()}, {out.data()});
    }

  private:
    CasADiFunctionEvaluator<1, 1> fun;
};

/// Wrapper for CasADiFunctionEvaluator with 2 vector inputs, 1 vector output.
class CasADiFun_2Vi1Vo {
  public:
    CasADiFun_2Vi1Vo(CasADiFunctionEvaluator<2, 1> &&fun)
        : fun(std::move(fun)) {}
    CasADiFun_2Vi1Vo(casadi::Function &&f,
                     const std::array<casadi_int, 2> &dim_in = {},
                     casadi_int dim_out                      = 0)
        : fun(std::move(f), {{dim_in[0], 1}, {dim_in[1], 1}}, {{dim_out, 1}}) {}

    void operator()(alpaqa::crvec in1, alpaqa::crvec in2, alpaqa::rvec out) const {
        fun({in1.data(), in2.data()}, {out.data()});
    }

  private:
    CasADiFunctionEvaluator<2, 1> fun;
};

/// Wrapper for CasADiFunctionEvaluator with 2 vector inputs, 1 matrix output.
class CasADiFun_2Vi1Mo {
  public:
    CasADiFun_2Vi1Mo(casadi::Function &&f,
                     const std::array<casadi_int, 2> &dim_in           = {},
                     CasADiFunctionEvaluator<2, 1>::casadi_dim dim_out = {0, 0})
        : fun(std::move(f), {{dim_in[0], 1}, {dim_in[1], 1}}, {dim_out}) {}

    void operator()(alpaqa::crvec in1, alpaqa::crvec in2, alpaqa::rmat out) const {
        fun({in1.data(), in2.data()}, {out.data()});
    }

  private:
    CasADiFunctionEvaluator<2, 1> fun;
};

/// Wrapper for CasADiFunctionEvaluator with 3 vector inputs, 1 matrix output.
class CasADiFun_3Vi1Mo {
  public:
    CasADiFun_3Vi1Mo(casadi::Function &&f,
                     const std::array<casadi_int, 3> &dim_in           = {},
                     CasADiFunctionEvaluator<3, 1>::casadi_dim dim_out = {0, 0})
        : fun(std::move(f),
              {
                  {dim_in[0], 1},
                  {dim_in[1], 1},
                  {dim_in[2], 1},
              },
              {dim_out}) {}

    void operator()(alpaqa::crvec in1, alpaqa::crvec in2, alpaqa::crvec in3,
                    alpaqa::rmat out) const {
        fun({in1.data(), in2.data(), in3.data()}, {out.data()});
    }

  private:
    CasADiFunctionEvaluator<3, 1> fun;
};

/// Wrapper for CasADiFunctionEvaluator with 3 vector inputs, 1 vector output.
class CasADiFun_3Vi1Vo {
  public:
    CasADiFun_3Vi1Vo(casadi::Function &&f,
                     const std::array<casadi_int, 3> &dim_in = {},
                     casadi_int dim_out                      = 0)
        : fun(std::move(f),
              {
                  {dim_in[0], 1},
                  {dim_in[1], 1},
                  {dim_in[2], 1},
              },
              {{dim_out, 1}}) {}

    void operator()(alpaqa::crvec in1, alpaqa::crvec in2, alpaqa::crvec in3,
                    alpaqa::rvec out) const {
        fun({in1.data(), in2.data(), in3.data()}, {out.data()});
    }

  private:
    CasADiFunctionEvaluator<3, 1> fun;
};

/// Wrapper for CasADiFunctionEvaluator with 4 vector inputs, 1 vector output.
class CasADiFun_4Vi1Vo {
  public:
    CasADiFun_4Vi1Vo(casadi::Function &&f,
                     const std::array<casadi_int, 4> &dim_in = {},
                     casadi_int dim_out                      = 0)
        : fun(std::move(f),
              {
                  {dim_in[0], 1},
                  {dim_in[1], 1},
                  {dim_in[2], 1},
                  {dim_in[3], 1},
              },
              {{dim_out, 1}}) {}

    void operator()(alpaqa::crvec in1, alpaqa::crvec in2, alpaqa::crvec in3, alpaqa::crvec in4,
                    alpaqa::rvec out) const {
        fun({in1.data(), in2.data(), in3.data(), in4.data()}, {out.data()});
    }

  private:
    CasADiFunctionEvaluator<4, 1> fun;
};

/// @}
