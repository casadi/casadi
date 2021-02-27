#pragma once

#include <iosfwd>
#include <memory>
#include <panoc-alm/problem.hpp>
#include <string>

/// Wrapper for CUTEst problems loaded from an external shared library.
///
/// @warning    The lifetime of the wrapper should be at least as long as the
///             lifetime of the @ref CUTEstProblem::problem member. Do not make
///             a copy of the problem that could outlive the wrapper.
class CUTEstProblem {

  public:
    /// Load a CUTEst problem from the given shared library and OUTSDIF.d file.
    CUTEstProblem(const char *so_fname, const char *outsdif_fname);
    /// @copydoc CUTEstProblem::CUTEstProblem(const char*, const char *)
    CUTEstProblem(const std::string &so_fname,
                  const std::string &outsdif_fname);
    CUTEstProblem(CUTEstProblem &&);
    CUTEstProblem &operator=(CUTEstProblem &&);
    ~CUTEstProblem();

  public:
    /// The report generated by CUTEst.
    ///
    /// @see `man CUTEST_creport` and `man CUTEST_ureport`
    struct Report {
        /// Name of the problem.
        std::string name;

        /// Number of independent variables.
        unsigned nvar = 0;
        /// Number of constraints.
        unsigned ncon = 0;

        enum Status {
            Success         = 0,    ///< Successful call.
            AllocationError = 1,    ///< Array allocation/deallocation error.
            ArrayBoundError = 2,    ///< Array bound error.
            EvaluationError = 3,    ///< Evaluation error.
        } status = Status::Success; ///< Exit status.

        /// Function call counters.
        ///
        /// @note   Note that hessian_times_vector, constraints and constraints_grad
        ///         may account for codes which allow the evaluation of a
        ///         selection of constraints only and may thus be much smaller
        ///         than the number of constraints times the number of
        ///         iterations.
        struct {
            /// Number of calls to the objective function.
            unsigned objective = 0;
            /// Number of calls to the objective gradient.
            unsigned objective_grad = 0;
            /// Number of calls to the objective Hessian.
            unsigned objective_hess = 0;
            /// Number of Hessian times vector products.
            unsigned hessian_times_vector = 0;
            /// Number of calls to the constraint functions.
            unsigned constraints = 0;
            /// Number of calls to the constraint gradients.
            unsigned constraints_grad = 0;
            /// Number of calls to the constraint Hessians.
            unsigned constraints_hess = 0;
        } calls; ///< Function call counters.

        /// CPU time (in seconds) for CUTEST_csetup.
        double time_setup = 0;
        /// CPU time (in seconds) since the end of CUTEST_csetup.
        double time = 0;
    };

    Report get_report() const;

  public:
    pa::Problem problem; ///< Problem statement (bounds, objective, constraints)
    pa::vec x0;          ///< Initial value of decision variables
    pa::vec y0;          ///< Initial value of Lagrange multipliers

  private:
    std::unique_ptr<class CUTEstLoader> implementation;
};

/// @related    CUTEstProblem::Report
std::ostream &operator<<(std::ostream &, const CUTEstProblem::Report &);
/// @related    CUTEstProblem::Report::Status
std::ostream &operator<<(std::ostream &, CUTEstProblem::Report::Status);
