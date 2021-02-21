#include <panoc-alm/interop/cutest/CUTEstLoader.hpp>

#include <cutest.h>
#include <dlfcn.h>
#include <stdexcept>
#include <string>

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
        auto rosenbr_open = dlfun<decltype(FORTRAN_open)>("fortran_open_");
        rosenbr_open(&funit, outsdif_fname, &ierr);
        if (ierr)
            throw std::runtime_error("Failed to open "s + outsdif_fname);

        // Get the dimensions of the problem
        integer status;
        auto rosenbr_cdimen = dlfun<decltype(CUTEST_cdimen)>("cutest_cdimen_");
        rosenbr_cdimen(&status, &funit, &nvar, &ncon);
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
            auto rosenbr_csetup =
                dlfun<decltype(CUTEST_csetup)>("cutest_cint_csetup_");
            rosenbr_csetup(&status, &funit, &iout, &io_buffer, &nvar, &ncon,
                           x.data(), x_l.data(), x_u.data(), y.data(),
                           c_l.data(), c_u.data(), (logical *)equatn.data(),
                           (logical *)linear.data(), &e_order, &l_order,
                           &v_order);
            throw_if_error("Failed to call cutest_csetup", status);
        }
        // Unconstrained problem
        else {
            auto rosenbr_usetup =
                dlfun<decltype(CUTEST_usetup)>("cutest_usetup_");
            rosenbr_usetup(&status, &funit, &iout, &io_buffer, &nvar, x.data(),
                           x_l.data(), x_u.data());
            throw_if_error("Failed to call cutest_usetup", status);
        }
        if (ncon > 0)
            work.resize(std::max(nvar, ncon));

        eval_objective_p        = ncon > 0                         //
                                      ? dlfun<void>("cutest_cfn_") //
                                      : dlfun<void>("cutest_ufn_");
        eval_objective_grad_p   = ncon > 0                               //
                                      ? dlfun<void>("cutest_cint_cofg_") //
                                      : dlfun<void>("cutest_ugr_");
        eval_constraints_p      = ncon > 0               //
                                      ? eval_objective_p //
                                      : nullptr;
        eval_constraints_grad_p = ncon > 0                                 //
                                      ? dlfun<void>("cutest_cint_cjprod_") //
                                      : nullptr;
    }

    // template <typename Ret, typename... Args>
    // static inline constexpr __attribute__((always_inline)) auto
    // call_as(void *funp, Ret (*)(Args...)) {
    //     // return [funp](Args... args) -> Ret {
    //     //     return ((Ret(*)(Args...))funp)(args...);
    //     // };
    //     return (Ret(*)(Args...))funp;
    // }
    template <class F>
    static inline constexpr auto call_as(void *funp) {
        return (F *)funp;
    }

    doublereal eval_objective_constrained(const pa::vec &x) const {
        assert(ncon > 0);
        integer status;
        doublereal f;
        call_as<decltype(CUTEST_cfn)>(eval_objective_p)(
            &status, &nvar, &ncon, x.data(), &f, work.data());
        throw_if_error("Failed to call cutest_cfn", status);
        return f;
    }

    doublereal eval_objective_unconstrained(const pa::vec &x) const {
        assert(ncon == 0);
        integer status;
        doublereal f;
        call_as<decltype(CUTEST_ufn)>(eval_objective_p)(&status, &nvar,
                                                        x.data(), &f);
        throw_if_error("Failed to call cutest_ufn", status);
        return f;
    }

    void eval_objective_grad_constrained(const pa::vec &x,
                                         pa::vec &grad_f) const {
        assert(ncon > 0);
        integer status;
        logical grad = true;
        call_as<decltype(CUTEST_cofg)>(eval_objective_grad_p)(
            &status, &nvar, x.data(), work.data(), grad_f.data(), &grad);
        throw_if_error("Failed to call cutest_cfn", status);
    }

    void eval_objective_grad_unconstrained(const pa::vec &x,
                                           pa::vec &grad_f) const {
        assert(ncon == 0);
        integer status;
        call_as<decltype(CUTEST_ugr)>(eval_objective_grad_p)(
            &status, &nvar, x.data(), grad_f.data());
        throw_if_error("Failed to call cutest_ugr", status);
    }

    void eval_constraints(const pa::vec &x, pa::vec &g) const {
        if (ncon == 0)
            return;
        integer status;
        call_as<decltype(CUTEST_cfn)>(eval_constraints_p)(
            &status, &nvar, &ncon, x.data(), work.data(), g.data());
        throw_if_error("Failed to call cutest_cfn", status);
    }

    void eval_constraints_grad(const pa::vec &x, const pa::vec &v,
                               pa::vec &grad_g_v) const {
        if (ncon == 0) {
            grad_g_v.setZero();
            return;
        }
        integer status;
        logical gotJ   = false;
        logical jtrans = true;
        call_as<decltype(CUTEST_cjprod)>(eval_constraints_p)(
            &status, &nvar, &ncon, &gotJ, &jtrans, x.data(), v.data(), &ncon,
            grad_g_v.data(), &nvar);
        throw_if_error("Failed to call cutest_cjprod", status);
    }

    ~CUTEstLoader() {
        // Terminate CUTEst
        if (ncon > 0) {
            integer status;
            auto rosenbr_cterminate =
                dlfun<decltype(CUTEST_cterminate)>("cutest_cterminate_");
            rosenbr_cterminate(&status);
            throw_if_error("Failed to call cutest_cterminate", status);
        } else {
            integer status;
            auto rosenbr_uterminate =
                dlfun<decltype(CUTEST_uterminate)>("cutest_uterminate_");
            rosenbr_uterminate(&status);
            throw_if_error("Failed to call cutest_uterminate", status);
        }

        // Close the OUTSDIF.d file
        integer ierr;
        auto rosenbr_close = dlfun<decltype(FORTRAN_close)>("fortran_close_");
        rosenbr_close(&funit, &ierr);
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
    pa::vec x;            ///< decision variable
    pa::vec x_l;          ///< lower bound on x
    pa::vec x_u;          ///< upper bound on x
    pa::vec y;            ///< lagrange multipliers
    pa::vec c_l;          ///< lower bounds on constraints
    pa::vec c_u;          ///< upper bounds on constraints
    logical_vec equatn;   ///< whether the constraint is an equality
    logical_vec linear;   ///< whether the constraint is linear
    mutable pa::vec work; ///< work vector

    void *eval_objective_p        = nullptr;
    void *eval_objective_grad_p   = nullptr;
    void *eval_constraints_p      = nullptr;
    void *eval_constraints_grad_p = nullptr;
};

CUTEstProblem::CUTEstProblem()                 = default;
CUTEstProblem::CUTEstProblem(CUTEstProblem &&) = default;
CUTEstProblem &CUTEstProblem::operator=(CUTEstProblem &&) = default;
CUTEstProblem::~CUTEstProblem()                           = default;

CUTEstProblem load_CUTEst_problem(const char *so_fname,
                                  const char *outsdif_fname) {
    CUTEstProblem p;
    p.implementation = std::make_unique<CUTEstLoader>(so_fname, outsdif_fname);
    auto *l          = p.implementation.get();
    p.problem.n      = l->nvar;
    p.problem.m      = l->ncon;
    p.problem.C.lowerbound = std::move(l->x_l);
    p.problem.C.upperbound = std::move(l->x_u);
    p.problem.D.lowerbound = std::move(l->c_l);
    p.problem.D.upperbound = std::move(l->c_u);
    using namespace std::placeholders;
    if (p.problem.m > 0) {
        p.problem.f =
            std::bind(&CUTEstLoader::eval_objective_constrained, l, _1);
        p.problem.grad_f = std::bind(
            &CUTEstLoader::eval_objective_grad_constrained, l, _1, _2);
    } else {
        p.problem.f =
            std::bind(&CUTEstLoader::eval_objective_unconstrained, l, _1);
        p.problem.grad_f = std::bind(
            &CUTEstLoader::eval_objective_grad_unconstrained, l, _1, _2);
    }
    p.problem.g = std::bind(&CUTEstLoader::eval_constraints, l, _1, _2);
    p.problem.grad_g =
        std::bind(&CUTEstLoader::eval_constraints_grad, l, _1, _2, _3);
    p.x0 = std::move(l->x);
    p.y0 = std::move(l->y);
    return p;
}
