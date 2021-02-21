#pragma once

#include <memory>
#include <panoc-alm/problem.hpp>

/// Wrapper for CUTEst problems loaded from an external shared library.
///
/// @warning    The lifetime of the wrapper should be at least as long as the
///             lifetime of the @ref CUTEstProblem::problem member. Do not make
///             a copy of the problem that could outlive the wrapper.
class CUTEstProblem {
  private:
    CUTEstProblem();

  public:
    CUTEstProblem(CUTEstProblem &&);
    CUTEstProblem &operator=(CUTEstProblem &&);
    ~CUTEstProblem();

  public:
    pa::Problem problem; ///< Problem statement (bounds, objective, constraints)
    pa::vec x0;          ///< Initial value of decision variables
    pa::vec y0;          ///< Initial value of Lagrange multipliers

  private:
    std::unique_ptr<class CUTEstLoader> implementation;

  public:
    friend CUTEstProblem load_CUTEst_problem(const char *so_fname,
                                             const char *outsdif_fname);
};

/// Load a CUTEst problem from the given shared library and OUTSDIF.d file.
/// @related    CUTEstProblem
CUTEstProblem load_CUTEst_problem(const char *so_fname,
                                  const char *outsdif_fname);
