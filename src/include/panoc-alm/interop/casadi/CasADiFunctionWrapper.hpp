#pragma once

#include <casadi/core/function.hpp>

#include <panoc-alm/util/vec.hpp>

#include <cassert>
#include <vector>

/// @addtogroup grp_ExternalProblemLoaders
/// @{

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

/// @}
