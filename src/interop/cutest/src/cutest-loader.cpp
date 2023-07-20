#include <alpaqa/cutest/cutest-loader.hpp>

#include <cutest.h>
#include <dlfcn.h>

#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

/*
CUTEST_cfn        : function and constraints values
CUTEST_cofg       : function value and possibly gradient
CUTEST_cofsg      : function value and possibly gradient
CUTEST_ccfg       : constraint functions values and possibly gradients
CUTEST_clfg       : Lagrangian function value and possibly gradient
CUTEST_cgr        : constraints gradients and gradient of objective/Lagrangian function
CUTEST_csgr       : constraints gradients and gradient of objective/Lagrangian function
CUTEST_ccfsg      : constraint functions values and possibly their gradients in sparse format
CUTEST_ccifg      : single constraint function value and possibly its gradient
CUTEST_ccifsg     : single constraint function value and possibly gradient in sparse format
CUTEST_cgrdh      : constraints gradients, Hessian of Lagrangian function and gradient of objective/Lagrangian function
CUTEST_cdh        : Hessian of the Lagrangian
CUTEST_cdhc       : Hessian of the constraint part of the Lagrangian
CUTEST_cshp       : sparsity pattern of the Hessian of the Lagrangian function
CUTEST_csh        : Hessian of the Lagrangian, in sparse format
CUTEST_cshc       : Hessian of the constraint part of the Lagrangian, in sparse format
CUTEST_ceh        : sparse Lagrangian Hessian matrix in finite element format
CUTEST_cifn       : problem function value
CUTEST_cigr       : gradient of a problem function
CUTEST_cisgr      : gradient of a problem function in sparse format
CUTEST_cidh       : Hessian of a problem function
CUTEST_cish       : Hessian of an individual problem function, in sparse format
CUTEST_csgrsh     : constraints gradients, sparse Lagrangian Hessian and the gradient of either the objective/Lagrangian in sparse format
CUTEST_csgreh     : constraint gradients, the Lagrangian Hessian in finite element format and the gradient of either the objective/Lagrangian in sparse format
CUTEST_chprod     : matrix-vector product of a vector with the Hessian matrix of the Lagrangian
CUTEST_cshprod    : matrix-vector product of a sparse vector with the Hessian matrix of the Lagrangian
CUTEST_chcprod    : matrix-vector product of a vector with the Hessian matrix of the constraint part of the Lagrangian
CUTEST_cshcprod   : matrix-vector product of a spaarse vector with the Hessian matrix of the constraint part of the Lagrangian
CUTEST_cjprod     : matrix-vector product of a vector with the Jacobian of the constraints, or its transpose
CUTEST_csjprod    : matrix-vector product of a sparse vector with the Jacobian of the constraints, or its transpose
CUTEST_cchprods   : matrix-vector products of a vector with each of the Hessian matrices of the constraint functions
*/

#define LOAD_DL_FUNC(f) dlfun<decltype(f)>(#f)

using namespace std::string_literals;

namespace {
void throw_error(std::string s, int code) {
    throw std::runtime_error(s + " (" + std::to_string(code) + ")");
}
void throw_if_error(std::string s, int code) {
    if (code)
        throw_error(s, code);
}
void log_if_error(std::string s, int code) {
    if (code)
        std::cerr << s << " (" << code << ")\n";
}

std::shared_ptr<void> load_lib(const char *so_filename) {
    assert(so_filename);
    ::dlerror();
    void *h = ::dlopen(so_filename, RTLD_LOCAL | RTLD_NOW);
    if (auto *err = ::dlerror())
        throw std::runtime_error(err);
    return std::shared_ptr<void>{h, &::dlclose};
}
} // namespace

namespace alpaqa {

class CUTEstLoader {
  public:
    USING_ALPAQA_CONFIG(CUTEstProblem::config_t);
    using Box = alpaqa::Box<config_t>;

  private:
    struct Cleanup {
        using fun_t = std::function<void()>;
        fun_t func;
        Cleanup() = default;
        Cleanup(fun_t func) : func{std::move(func)} {}
        Cleanup(Cleanup &&) noexcept = default;
        ~Cleanup() { func(); }
    };

    Cleanup load_outsdif(const char *outsdif_fname) {
        integer ierr;
        auto fptr_open  = LOAD_DL_FUNC(FORTRAN_open);
        auto fptr_close = LOAD_DL_FUNC(FORTRAN_close);
        fptr_open(&funit, outsdif_fname, &ierr);
        throw_if_error("Failed to open "s + outsdif_fname, ierr);
        return {[this, fptr_close] {
            integer ierr;
            fptr_close(&funit, &ierr);
            log_if_error("Failed to close OUTSDIF.d file", ierr);
        }};
    }

  public:
    CUTEstLoader(const char *so_fname, const char *outsdif_fname) {
        // Open the shared library
        so_handle = load_lib(so_fname);

        // Open the OUTSDIF.d file
        cleanup_outsdif = load_outsdif(outsdif_fname);

        // Get the dimensions of the problem
        integer status;
        auto fptr_cdimen = LOAD_DL_FUNC(CUTEST_cdimen);
        fptr_cdimen(&status, &funit, &nvar, &ncon);
        throw_if_error("Failed to call cutest_cdimen", status);
    }

    struct ConstrFuncs {
        decltype(CUTEST_cfn) *obj_constr;
        decltype(CUTEST_cint_cofg) *obj_constr_grad;
        decltype(CUTEST_cjprod) *constr_grad_prod;
        decltype(CUTEST_ccifg) *constr_grad_i;
        decltype(CUTEST_cdh) *lagr_hess;
        decltype(CUTEST_csh) *lagr_hess_sparse;
        decltype(CUTEST_cdh) *lagr_hess_prod;
    };

    void setup_problem(rvec x0, rvec y0, Box &C, Box &D) {
        assert(x.size() == nvar);
        assert(C.lowerbound.size() == nvar);
        assert(C.upperbound.size() == nvar);
        assert(y.size() == ncon);
        assert(D.lowerbound.size() == ncon);
        assert(D.upperbound.size() == ncon);
        equatn.resize(ncon);
        linear.resize(ncon);
        integer e_order = 0; // no specific order of equality constraints
        integer l_order = 0; // no specific order of linear constraints
        integer v_order = 0; // no specific order of linear variables
        integer status;

        // Constrained problem
        if (ncon > 0) {
            LOAD_DL_FUNC(CUTEST_csetup)
            (&status, &funit, &iout, &io_buffer, &nvar, &ncon, x.data(),
             C.lowerbound.data(), C.upperbound.data(), y.data(),
             D.lowerbound.data(), D.upperbound.data(), equatn.data(),
             linear.data(), &e_order, &l_order, &v_order);
            throw_if_error("Failed to call cutest_csetup", status);
        }
        // Unconstrained problem
        else {
            LOAD_DL_FUNC(CUTEST_usetup)
            (&status, &funit, &iout, &io_buffer, &nvar, x.data(),
             C.lowerbound.data(), C.upperbound.data());
            throw_if_error("Failed to call cutest_usetup", status);
        }

        if (ncon > 0) {
            work.resize(std::max(nvar, ncon));
            return ConstrFuncs {
                .obj_constr
            };
        }

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
    inline static constexpr auto call_as(void *funp) {
        return (F *)funp;
    }

    doublereal eval_objective_constrained(pa::crvec x) const {
        assert(x.size() == nvar);
        assert(ncon > 0);
        integer status;
        doublereal f;
        call_as<decltype(CUTEST_cfn)>(eval_obj_p)(&status, &nvar, &ncon,
                                                  x.data(), &f, work.data());
        throw_if_error("Failed to call cutest_cfn", status);
        return f;
    }

    doublereal eval_objective_unconstrained(pa::crvec x) const {
        assert(x.size() == nvar);
        assert(ncon == 0);
        integer status;
        doublereal f;
        call_as<decltype(CUTEST_ufn)>(eval_obj_p)(&status, &nvar, x.data(), &f);
        throw_if_error("Failed to call cutest_ufn", status);
        return f;
    }

    void eval_objective_grad_constrained(pa::crvec x, pa::rvec grad_f) const {
        assert(x.size() == nvar);
        assert(grad_f.size() == nvar);
        assert(ncon > 0);
        integer status;
        logical grad = true;
        call_as<decltype(CUTEST_cofg)>(eval_obj_grad_p)(
            &status, &nvar, x.data(), work.data(), grad_f.data(), &grad);
        throw_if_error("Failed to call cutest_cfn", status);
    }

    void eval_objective_grad_unconstrained(pa::crvec x, pa::rvec grad_f) const {
        assert(x.size() == nvar);
        assert(grad_f.size() == nvar);
        assert(ncon == 0);
        integer status;
        call_as<decltype(CUTEST_ugr)>(eval_obj_grad_p)(&status, &nvar, x.data(),
                                                       grad_f.data());
        throw_if_error("Failed to call cutest_ugr", status);
    }

    void eval_constraints(pa::crvec x, pa::rvec g) const {
        assert(x.size() == nvar);
        assert(g.size() == ncon);
        if (ncon == 0)
            return;
        integer status;
        call_as<decltype(CUTEST_cfn)>(eval_constr_p)(
            &status, &nvar, &ncon, x.data(), work.data(), g.data());
        throw_if_error("Failed to call cutest_cfn", status);
    }

    void eval_constraints_grad_prod(pa::crvec x, pa::crvec v,
                                    pa::rvec grad_g_v) const {
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

    void eval_constraint_i_grad(pa::crvec x, unsigned i,
                                pa::rvec grad_gi) const {
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
        pa::real_t ci;
        logical grad = true;
        call_as<decltype(CUTEST_ccifg)>(eval_constr_i_grad_p)(
            &status, &nvar, &icon, x.data(), &ci, grad_gi.data(), &grad);
        throw_if_error("Failed to call cutest_ccifg", status);
    }

    void eval_lagr_hess_prod(pa::crvec x, pa::crvec y, pa::crvec v,
                             pa::rvec Hv) const {
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
                const_cast<pa::real_t *>(v.data()), Hv.data());
            // TODO: why is the VECTOR argument not const?
            throw_if_error("Failed to call cutest_chprod", status);
        }
    }

    void eval_lagr_hess(pa::crvec x, pa::crvec y, pa::rmat H) const {
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
                             [](pa::real_t x) { return x > -CUTE_INF; }) +
               std::count_if(x_u.data(), x_u.data() + nvar,
                             [](pa::real_t x) { return x < +CUTE_INF; });
    }

    std::string get_name() {
        std::string name;
        name.resize(10);
        integer status;
        LOAD_DL_FUNC(CUTEST_probname)(&status, name.data());
        throw_if_error("Failed to call cutest_probname", status);
        auto nspace = std::find_if(name.rbegin(), name.rend(),
                                   [](char c) { return c != ' '; });
        name.resize(name.size() - (nspace - name.rbegin()));
        return name;
    }

    integer get_report(doublereal *calls, doublereal *time) {
        integer status;
        auto fptr_report = ncon > 0 ? LOAD_DL_FUNC(CUTEST_creport)
                                    : LOAD_DL_FUNC(CUTEST_ureport);
        fptr_report(&status, calls, time);
        return status;
    }

    ~CUTEstLoader() {
        // Terminate CUTEst
        if (ncon > 0) {
            integer status;
            auto fptr_cterminate = LOAD_DL_FUNC(CUTEST_cterminate);
            fptr_cterminate(&status);
            log_if_error("Failed to call cutest_cterminate", status);
        } else {
            integer status;
            auto fptr_uterminate = LOAD_DL_FUNC(CUTEST_uterminate);
            fptr_uterminate(&status);
            log_if_error("Failed to call cutest_uterminate", status);
        }
    }

    template <class T>
    T *dlfun(const char *name) {
        (void)dlerror();
        auto res = reinterpret_cast<T *>(dlsym(so_handle, name));
        if (const char *error = dlerror())
            throw std::runtime_error(error);
        return res;
    }

    std::shared_ptr<void> so_handle; ///< dlopen handle to shared library
    Cleanup cleanup_outsdif;  ///< Responsible for closing the OUTSDIF.d file
    Cleanup cutest_terminate; ///< Responsible for calling CUTEST_xterminate

    integer funit     = 42; ///< Fortran Unit Number for OUTSDIF.d file
    integer iout      = 6;  ///< Fortran Unit Number for standard output
    integer io_buffer = 11; ///< Fortran Unit Number for internal IO

    integer nvar; ///< Number of decision variabls
    integer ncon; ///< Number of constraints

    using logical_vec = Eigen::VectorX<logical>;
    vec x;              ///< decision variable
    vec x_l;            ///< lower bound on x
    vec x_u;            ///< upper bound on x
    vec y;              ///< lagrange multipliers
    vec c_l;            ///< lower bounds on constraints
    vec c_u;            ///< upper bounds on constraints
    logical_vec equatn; ///< whether the constraint is an equality
    logical_vec linear; ///< whether the constraint is linear
    mutable vec work;   ///< work vector

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

CUTEstProblem::CUTEstProblem(CUTEstProblem &&)            = default;
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
        "invalid value for pa::CUTEstProblem::Report::Status");
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

} // namespace alpaqa