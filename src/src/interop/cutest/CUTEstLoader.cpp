#include <alpaqa/interop/cutest/CUTEstLoader.hpp>

#include <cutest.h>
#include <dlfcn.h>

#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

using namespace std::string_literals;

class CUTEstLoader {
  private:
    static void throw_error(std::string s, int code) {
        throw std::runtime_error(s + " (" + std::to_string(code) + ")");
    }
    static void throw_if_error(std::string s, int code) {
        if (code)
            throw_error(s, code);
    }

  public:
    CUTEstLoader(const char *so_fname, const char *outsdif_fname) {
        // Open the shared library
        so_handle = dlopen(so_fname, RTLD_LAZY);
        if (!so_handle)
            throw std::runtime_error("Failed to open "s + so_fname);

        // Open the OUTSDIF.d file
        integer ierr;
        auto fptr_open = dlfun<decltype(FORTRAN_open)>("fortran_open_");
        fptr_open(&funit, outsdif_fname, &ierr);
        if (ierr)
            throw std::runtime_error("Failed to open "s + outsdif_fname);

        // Get the dimensions of the problem
        integer status;
        auto fptr_cdimen = dlfun<decltype(CUTEST_cdimen)>("cutest_cdimen_");
        fptr_cdimen(&status, &funit, &nvar, &ncon);
        throw_if_error("Failed to call cutest_cdimen", status);

        // Set up the datastructures
        x.resize(nvar);
        x_l.resize(nvar);
        x_u.resize(nvar);
        y.resize(ncon);
        c_l.resize(ncon);
        c_u.resize(ncon);
        equatn.resize(ncon);
        linear.resize(ncon);
        integer e_order = 0;
        integer l_order = 0;
        integer v_order = 0;

        // Constrained problem
        if (ncon > 0) {
            auto fptr_csetup =
                dlfun<decltype(CUTEST_csetup)>("cutest_cint_csetup_");
            fptr_csetup(&status, &funit, &iout, &io_buffer, &nvar, &ncon,
                        x.data(), x_l.data(), x_u.data(), y.data(), c_l.data(),
                        c_u.data(), (logical *)equatn.data(),
                        (logical *)linear.data(), &e_order, &l_order, &v_order);
            throw_if_error("Failed to call cutest_csetup", status);
        }
        // Unconstrained problem
        else {
            auto fptr_usetup = dlfun<decltype(CUTEST_usetup)>("cutest_usetup_");
            fptr_usetup(&status, &funit, &iout, &io_buffer, &nvar, x.data(),
                        x_l.data(), x_u.data());
            throw_if_error("Failed to call cutest_usetup", status);
        }
        if (ncon > 0)
            work.resize(std::max(nvar, ncon));

        eval_obj_p              = ncon > 0                         //
                                      ? dlfun<void>("cutest_cfn_") //
                                      : dlfun<void>("cutest_ufn_");
        eval_obj_grad_p         = ncon > 0                               //
                                      ? dlfun<void>("cutest_cint_cofg_") //
                                      : dlfun<void>("cutest_ugr_");
        eval_constr_p           = ncon > 0         //
                                      ? eval_obj_p //
                                      : nullptr;
        eval_constr_grad_prod_p = ncon > 0                                 //
                                      ? dlfun<void>("cutest_cint_cjprod_") //
                                      : nullptr;
        eval_constr_i_grad_p    = ncon > 0                                //
                                      ? dlfun<void>("cutest_cint_ccifg_") //
                                      : nullptr;
        eval_lagr_hess_prod_p   = ncon > 0 //
                                      ? dlfun<void>("cutest_cint_chprod_")
                                      : dlfun<void>("cutest_cint_uhprod_");
        eval_lagr_hess_p        = ncon > 0                         //
                                      ? dlfun<void>("cutest_cdh_") //
                                      : dlfun<void>("cutest_udh_");
    }

    template <class F>
    static inline constexpr auto call_as(void *funp) {
        return (F *)funp;
    }

    doublereal eval_objective_constrained(alpaqa::crvec x) const {
        assert(x.size() == nvar);
        assert(ncon > 0);
        integer status;
        doublereal f;
        call_as<decltype(CUTEST_cfn)>(eval_obj_p)(&status, &nvar, &ncon,
                                                  x.data(), &f, work.data());
        throw_if_error("Failed to call cutest_cfn", status);
        return f;
    }

    doublereal eval_objective_unconstrained(alpaqa::crvec x) const {
        assert(x.size() == nvar);
        assert(ncon == 0);
        integer status;
        doublereal f;
        call_as<decltype(CUTEST_ufn)>(eval_obj_p)(&status, &nvar, x.data(), &f);
        throw_if_error("Failed to call cutest_ufn", status);
        return f;
    }

    void eval_objective_grad_constrained(alpaqa::crvec x, alpaqa::rvec grad_f) const {
        assert(x.size() == nvar);
        assert(grad_f.size() == nvar);
        assert(ncon > 0);
        integer status;
        logical grad = true;
        call_as<decltype(CUTEST_cofg)>(eval_obj_grad_p)(
            &status, &nvar, x.data(), work.data(), grad_f.data(), &grad);
        throw_if_error("Failed to call cutest_cfn", status);
    }

    void eval_objective_grad_unconstrained(alpaqa::crvec x, alpaqa::rvec grad_f) const {
        assert(x.size() == nvar);
        assert(grad_f.size() == nvar);
        assert(ncon == 0);
        integer status;
        call_as<decltype(CUTEST_ugr)>(eval_obj_grad_p)(&status, &nvar, x.data(),
                                                       grad_f.data());
        throw_if_error("Failed to call cutest_ugr", status);
    }

    void eval_constraints(alpaqa::crvec x, alpaqa::rvec g) const {
        assert(x.size() == nvar);
        assert(g.size() == ncon);
        if (ncon == 0)
            return;
        integer status;
        call_as<decltype(CUTEST_cfn)>(eval_constr_p)(
            &status, &nvar, &ncon, x.data(), work.data(), g.data());
        throw_if_error("Failed to call cutest_cfn", status);
    }

    void eval_constraints_grad_prod(alpaqa::crvec x, alpaqa::crvec v,
                                    alpaqa::rvec grad_g_v) const {
        assert(x.size() == nvar);
        assert(v.size() == ncon);
        assert(grad_g_v.size() == nvar);
        if (ncon == 0) {
            grad_g_v.setZero();
            return;
        }
        integer status;
        logical gotJ   = false;
        logical jtrans = true;
        call_as<decltype(CUTEST_cjprod)>(eval_constr_grad_prod_p)(
            &status, &nvar, &ncon, &gotJ, &jtrans, x.data(), v.data(), &ncon,
            grad_g_v.data(), &nvar);
        throw_if_error("Failed to call cutest_cjprod", status);
    }

    void eval_constraint_i_grad(alpaqa::crvec x, unsigned i,
                                alpaqa::rvec grad_gi) const {
        assert(x.size() == nvar);
        assert(grad_gi.size() == nvar);
        if (ncon == 0) {
            grad_gi.setZero();
            return;
        }
        integer status;
        integer icon = i + 1;
        assert(icon >= 1);
        assert(icon <= ncon);
        alpaqa::real_t ci;
        logical grad = true;
        call_as<decltype(CUTEST_ccifg)>(eval_constr_i_grad_p)(
            &status, &nvar, &icon, x.data(), &ci, grad_gi.data(), &grad);
        throw_if_error("Failed to call cutest_ccifg", status);
    }

    void eval_lagr_hess_prod(alpaqa::crvec x, alpaqa::crvec y, alpaqa::crvec v,
                             alpaqa::rvec Hv) const {
        assert(x.size() == nvar);
        assert(y.size() == ncon);
        assert(v.rows() == nvar);
        assert(Hv.rows() == nvar);
        integer status;
        logical gotH = false;
        if (ncon == 0) {
            call_as<decltype(CUTEST_uhprod)>(eval_lagr_hess_prod_p)(
                &status, &nvar, &gotH, x.data(), v.data(), Hv.data());
            throw_if_error("Failed to call cutest_uhprod", status);
        } else {
            call_as<decltype(CUTEST_chprod)>(eval_lagr_hess_prod_p)(
                &status, &nvar, &ncon, &gotH, x.data(), y.data(),
                const_cast<alpaqa::real_t *>(v.data()), Hv.data());
            // TODO: why is the VECTOR argument not const?
            throw_if_error("Failed to call cutest_chprod", status);
        }
    }

    void eval_lagr_hess(alpaqa::crvec x, alpaqa::crvec y, alpaqa::rmat H) const {
        assert(x.size() == nvar);
        assert(y.size() == ncon);
        assert(H.rows() >= nvar);
        assert(H.cols() >= nvar);
        integer status;
        integer ldH = H.rows();
        if (ncon == 0) {
            call_as<decltype(CUTEST_udh)>(eval_lagr_hess_p)(
                &status, &nvar, x.data(), &ldH, H.data());
            throw_if_error("Failed to call cutest_udh", status);
        } else {
            call_as<decltype(CUTEST_cdh)>(eval_lagr_hess_p)(
                &status, &nvar, &ncon, x.data(), y.data(), &ldH, H.data());
            throw_if_error("Failed to call cutest_cdh", status);
        }
    }

    unsigned count_box_constraints() const {
        return std::count_if(x_l.data(), x_l.data() + nvar,
                             [](alpaqa::real_t x) { return x > -CUTE_INF; }) +
               std::count_if(x_u.data(), x_u.data() + nvar,
                             [](alpaqa::real_t x) { return x < +CUTE_INF; });
    }

    std::string get_name() {
        std::string name;
        name.resize(10);
        integer status;
        dlfun<decltype(CUTEST_probname)>("cutest_probname_")(&status,
                                                             name.data());
        throw_if_error("Failed to call cutest_probname", status);
        auto nspace = std::find_if(name.rbegin(), name.rend(),
                                   [](char c) { return c != ' '; });
        name.resize(name.size() - (nspace - name.rbegin()));
        return name;
    }

    integer get_report(doublereal *calls, doublereal *time) {
        integer status;
        auto fptr_report =
            ncon > 0 ? dlfun<decltype(CUTEST_creport)>("cutest_creport_")
                     : dlfun<decltype(CUTEST_ureport)>("cutest_ureport_");
        fptr_report(&status, calls, time);
        return status;
    }

    ~CUTEstLoader() {
        // Terminate CUTEst
        if (ncon > 0) {
            integer status;
            auto fptr_cterminate =
                dlfun<decltype(CUTEST_cterminate)>("cutest_cterminate_");
            fptr_cterminate(&status);
            throw_if_error("Failed to call cutest_cterminate", status);
        } else {
            integer status;
            auto fptr_uterminate =
                dlfun<decltype(CUTEST_uterminate)>("cutest_uterminate_");
            fptr_uterminate(&status);
            throw_if_error("Failed to call cutest_uterminate", status);
        }

        // Close the OUTSDIF.d file
        integer ierr;
        auto fptr_close = dlfun<decltype(FORTRAN_close)>("fortran_close_");
        fptr_close(&funit, &ierr);
        throw_if_error("Failed to close OUTSDIF.d file", ierr);

        // Close the shared library
        if (so_handle)
            dlclose(so_handle);
    }

    template <class T>
    T *dlfun(const char *name) {
        (void)dlerror();
        auto res = (T *)dlsym(so_handle, name);
        if (const char *error = dlerror(); error)
            throw std::runtime_error(error);
        return res;
    }

    void *so_handle = nullptr; ///< dlopen handle to shared library

    integer funit     = 42; ///< Fortran Unit Number for OUTSDIF.d file
    integer iout      = 6;  ///< Fortran Unit Number for standard output
    integer io_buffer = 11; ///< Fortran Unit Number for internal IO

    integer nvar; ///< Number of decision variabls
    integer ncon; ///< Number of constraints

    using logical_vec = Eigen::Matrix<logical, Eigen::Dynamic, 1>;
    alpaqa::vec x;            ///< decision variable
    alpaqa::vec x_l;          ///< lower bound on x
    alpaqa::vec x_u;          ///< upper bound on x
    alpaqa::vec y;            ///< lagrange multipliers
    alpaqa::vec c_l;          ///< lower bounds on constraints
    alpaqa::vec c_u;          ///< upper bounds on constraints
    logical_vec equatn;   ///< whether the constraint is an equality
    logical_vec linear;   ///< whether the constraint is linear
    mutable alpaqa::vec work; ///< work vector

    void *eval_obj_p              = nullptr;
    void *eval_obj_grad_p         = nullptr;
    void *eval_constr_p           = nullptr;
    void *eval_constr_grad_prod_p = nullptr;
    void *eval_constr_i_grad_p    = nullptr;
    void *eval_lagr_hess_prod_p   = nullptr;
    void *eval_lagr_hess_p        = nullptr;
};

CUTEstProblem::CUTEstProblem(const char *so_fname, const char *outsdif_fname) {
    implementation = std::make_unique<CUTEstLoader>(so_fname, outsdif_fname);
    auto *l        = implementation.get();
    name           = l->get_name();
    number_box_constraints = l->count_box_constraints();
    problem.n              = l->nvar;
    problem.m              = l->ncon;
    problem.C.lowerbound   = std::move(l->x_l);
    problem.C.upperbound   = std::move(l->x_u);
    problem.D.lowerbound   = std::move(l->c_l);
    problem.D.upperbound   = std::move(l->c_u);
    using namespace std::placeholders;
    if (problem.m > 0) {
        problem.f = std::bind(&CUTEstLoader::eval_objective_constrained, l, _1);
        problem.grad_f = std::bind(
            &CUTEstLoader::eval_objective_grad_constrained, l, _1, _2);
    } else {
        problem.f =
            std::bind(&CUTEstLoader::eval_objective_unconstrained, l, _1);
        problem.grad_f = std::bind(
            &CUTEstLoader::eval_objective_grad_unconstrained, l, _1, _2);
    }
    problem.g = std::bind(&CUTEstLoader::eval_constraints, l, _1, _2);
    problem.grad_g_prod =
        std::bind(&CUTEstLoader::eval_constraints_grad_prod, l, _1, _2, _3);
    problem.grad_gi =
        std::bind(&CUTEstLoader::eval_constraint_i_grad, l, _1, _2, _3);
    problem.hess_L_prod =
        std::bind(&CUTEstLoader::eval_lagr_hess_prod, l, _1, _2, _3, _4);
    problem.hess_L = std::bind(&CUTEstLoader::eval_lagr_hess, l, _1, _2, _3);
    x0             = std::move(l->x);
    y0             = std::move(l->y);
}

CUTEstProblem::CUTEstProblem(const std::string &so_fname,
                             const std::string &outsdif_fname)
    : CUTEstProblem(so_fname.c_str(), outsdif_fname.c_str()) {}

CUTEstProblem::CUTEstProblem(CUTEstProblem &&) = default;
CUTEstProblem &CUTEstProblem::operator=(CUTEstProblem &&) = default;
CUTEstProblem::~CUTEstProblem()                           = default;

CUTEstProblem::Report CUTEstProblem::get_report() const {
    doublereal calls[7];
    doublereal time[2];
    Report r;
    using stat_t = decltype(r.status);
    r.status     = static_cast<stat_t>(implementation->get_report(calls, time));
    r.name       = implementation->get_name();
    r.nvar       = implementation->nvar;
    r.ncon       = implementation->ncon;
    r.calls.objective            = calls[0];
    r.calls.objective_grad       = calls[1];
    r.calls.objective_hess       = calls[2];
    r.calls.hessian_times_vector = calls[3];
    r.calls.constraints          = implementation->ncon > 0 ? calls[4] : 0;
    r.calls.constraints_grad     = implementation->ncon > 0 ? calls[5] : 0;
    r.calls.constraints_hess     = implementation->ncon > 0 ? calls[6] : 0;
    r.time_setup                 = time[0];
    r.time                       = time[1];
    return r;
}

const char *enum_name(CUTEstProblem::Report::Status s) {
    using Status = CUTEstProblem::Report::Status;
    switch (s) {
        case Status::Success: return "Success";
        case Status::AllocationError: return "AllocationError";
        case Status::ArrayBoundError: return "ArrayBoundError";
        case Status::EvaluationError: return "EvaluationError";
    }
    throw std::out_of_range(
        "invalid value for alpaqa::CUTEstProblem::Report::Status");
}

std::ostream &operator<<(std::ostream &os, CUTEstProblem::Report::Status s) {
    return os << enum_name(s);
}

std::ostream &operator<<(std::ostream &os, const CUTEstProblem::Report &r) {
    os << "CUTEst problem: " << r.name << "\r\n\n"                 //
       << "Number of variables:   " << r.nvar << "\r\n"            //
       << "Number of constraints: " << r.ncon << "\r\n\n"          //
       << "Status: " << r.status << " (" << +r.status << ")\r\n\n" //
       << "Objective function evaluations:            " << r.calls.objective
       << "\r\n"
       << "Objective function gradient evaluations:   "
       << r.calls.objective_grad << "\r\n"
       << "Objective function Hessian evaluations:    "
       << r.calls.objective_hess << "\r\n"
       << "Hessian times vector products:             "
       << r.calls.objective_hess << "\r\n\n";
    if (r.ncon > 0) {
        os << "Constraint function evaluations:           "
           << r.calls.constraints << "\r\n"
           << "Constraint function gradients evaluations: "
           << r.calls.constraints_grad << "\r\n"
           << "Constraint function Hessian evaluations:   "
           << r.calls.constraints_hess << "\r\n\n";
    }
    return os << "Setup time:       " << r.time_setup << "s\r\n"
              << "Time since setup: " << r.time << "s";
}