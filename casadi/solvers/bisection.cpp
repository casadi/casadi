

#include "bisection.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace casadi
{
    extern "C" int CASADI_ROOTFINDER_BISECTION_EXPORT casadi_register_rootfinder_bisection(Rootfinder::Plugin *plugin)
    {
        plugin->creator = Bisection::creator;
        plugin->name = "bisection";
        plugin->doc = Bisection::meta_doc.c_str();
        plugin->version = CASADI_VERSION;
        plugin->options = &Bisection::options_;
        plugin->deserialize = &Bisection::deserialize;
        return 0;
    }

    extern "C" void CASADI_ROOTFINDER_BISECTION_EXPORT casadi_load_rootfinder_bisection()
    {
        Rootfinder::registerPlugin(casadi_register_rootfinder_bisection);
    }

    Bisection::Bisection(const std::string &name, const Function &f)
        : Rootfinder(name, f)
    {
    }

    Bisection::~Bisection()
    {
        clear_mem();
    }

    const Options Bisection::options_ = {{&Rootfinder::options_},
                                         {
                                             {"abstol", {OT_DOUBLE, "Stopping criterion tolerance on ||g||__inf)"}},
                                             {"abstol_step", {OT_DOUBLE, "Stopping tolerance on bracket width"}},
                                             {"max_iter", {OT_INT, "Maximum number of Newton iterations to perform before returning."}},
                                             {"lb", {OT_DOUBLE, "lower bound"}},
                                             {"ub", {OT_DOUBLE, "upper bound"}},
                                             {"search_step", {OT_DOUBLE, "Step size for bracket searching"}},
                                             {"max_search", {OT_INT, "Maximum bracket search iterations"}},
                                         }};

    void Bisection::init(const Dict &opts)
    {
        Rootfinder::init(opts);

        max_iter_ = 100;
        abstol_ = 1e-9;
        abstol_step_ = 1e-9;
        lb_ = -1e-12;
        ub_ = 1e12;
        search_step_ = 1.0;
        max_search_ = 100;

        for (auto &op : opts)
        {
            if (op.first == "max_iter")
            {
                max_iter_ = op.second;
            }
            else if (op.first == "abstol")
            {
                abstol_ = op.second;
            }
            else if (op.first == "abstol_step")
            {
                abstol_step_ = op.second;
            }
            else if (op.first == "lb")
            {
                lb_ = op.second;
            }
            else if (op.first == "ub")
            {
                ub_ = op.second;
            }
            else if (op.first == "search_step")
            {
                search_step_ = op.second;
            }
            else if (op.first == "max_search")
            {
                max_search_ = op.second;
            }
        }

        casadi_assert(oracle_.n_in() > 0, "Bisection: the supplied f must have at least one input.");
        casadi_assert(n_ == 1, "Bisection only supports scalar equations (n=1).");
        casadi_assert(lb_ < ub_, "lb must be strictly less than ub.");
    }

    int Bisection::init_mem(void *mem) const
    {
        if (Rootfinder::init_mem(mem)) return 1;
        auto m = static_cast<BisectionMemory *>(mem);
        m->return_status = 0;
        m->iter = 0;
        m->search_iter = 0;
        return 0;
    }

    void Bisection::set_work(void *mem, const double **&arg, double **&res, casadi_int *&iw, double *&w) const
    {
        Rootfinder::set_work(mem, arg, res, iw, w);
    }

    int Bisection::solve(void *mem) const
    {
        auto m = static_cast<BisectionMemory *>(mem);

        double f_val = 0.0;

        auto eval_f = [&](double x) -> double
        {
            for (casadi_int i = 0; i < n_in_; ++i)
                m->arg[i] = m->iarg[i];

            m->arg[iin_] = &x;

            for (casadi_int i = 0; i < n_out_; ++i)
                m->res[i] = nullptr;
            m->res[iout_] = &f_val;

            if (oracle_(m->arg, m->res, m->iw, m->w, 0))
            {
                f_val = std::numeric_limits<double>::quiet_NaN();
            }
            return f_val;
        };

        double x0 = m->iarg[iin_][0];
        x0 = std::max(lb_, std::min(ub_, x0));

        double f0 = eval_f(x0);

        if (std::isnan(f0)) return finish(m, x0, f0, 0.0, -1, false, SOLVER_RET_UNKNOWN);
        if (std::fabs(f0) < abstol_) return finish(m, x0, f0, 0, 2, true, SOLVER_RET_SUCCESS);

        bool bracketed = false;
        double a = x0, b = x0;
        double fa = f0, fb = f0;

        for (m->search_iter = 1; m->search_iter <= max_search_; ++m->search_iter)
        {
            if (a > lb_)
            {
                a = std::max(lb_, a - search_step_);
                fa = eval_f(a);
            }
            if (b < ub_)
            {
                b = std::min(ub_, b + search_step_);
                fb = eval_f(b);
            }

            if (std::isnan(fa) || std::isnan(fb)) return finish(m, a, fa, b - a, -1, false, SOLVER_RET_UNKNOWN);

            if (fa * fb <= 0.0)
            {
                bracketed = true;
                break;
            }

            if (a == lb_ && b == ub_ && fa * fb > 0.0)
            {
                break;
            }
        }

        if (!bracketed)
        {
            return finish(m, x0, f0, b - a, -2, false, SOLVER_RET_UNKNOWN);
        }

        if (fa == 0.0) return finish(m, a, 0.0, b - a, 2, true, SOLVER_RET_SUCCESS);
        if (fb == 0.0) return finish(m, b, 0.0, b - a, 2, true, SOLVER_RET_SUCCESS);

        double mid = a, f_mid_val = fa;

        for (m->iter = 0; m->iter < max_iter_; ++m->iter)
        {
            mid = 0.5 * (a + b);
            f_mid_val = eval_f(mid);

            if (std::isnan(f_mid_val))
                return finish(m, mid, f_mid_val, b - a, -1, false, SOLVER_RET_UNKNOWN);

            if (std::fabs(f_mid_val) < abstol_)
                return finish(m, mid, f_mid_val, b - a, 2, true, SOLVER_RET_SUCCESS);

            if ((b - a) < abstol_step_)
                return finish(m, mid, f_mid_val, b - a, 1, true, SOLVER_RET_SUCCESS);

            if (fa * f_mid_val < 0.0)
            {
                b = mid;
                fb = f_mid_val;
            }
            else
            {
                a = mid;
                fa = f_mid_val;
            }
        }

        return finish(m, mid, f_mid_val, b - a, 0, false, SOLVER_RET_LIMITED);
    }

    int Bisection::finish(BisectionMemory *m, double x_sol, double f_sol, double width, int status, bool success, UnifiedReturnStatus urs) const
    {
        casadi_copy(&x_sol, 1, m->ires[iout_]);
        m->return_status = status;
        m->f_mid = f_sol;
        m->bracket_width = width;
        m->success = success;
        m->unified_return_status = urs;
        return 0;
    }

    std::string Bisection::status_str(int status)
    {
        switch (status)
        {
        case 0:
            return "max_iteration_reached";
        case 1:
            return "converged_bracket";
        case 2:
            return "converged_abstol";
        case -1:
            return "nan_encountered";
        case -2:
            return "failed_to_bracket_root";
        default:
            return "unknown";
        }
    }

    Dict Bisection::get_stats(void *mem) const
    {
        Dict stats = Rootfinder::get_stats(mem);
        auto m = static_cast<BisectionMemory *>(mem);
        stats["return_status"] = status_str(m->return_status);
        stats["search_iter"] = m->search_iter; // 输出给 Python 端查看探索了多少步
        stats["iter_count"] = m->iter;
        stats["f_mid"] = m->f_mid;
        stats["bracket_width"] = m->bracket_width;
        return stats;
    }

    void Bisection::serialize_body(SerializingStream &s) const
    {
        Rootfinder::serialize_body(s);
        s.version("Bisection", 1); // 升级版本号
        s.pack("Bisection::max_iter", max_iter_);
        s.pack("Bisection::max_search", max_search_);
        s.pack("Bisection::search_step", search_step_);
        s.pack("Bisection::abstol", abstol_);
        s.pack("Bisection::abstol_step", abstol_step_);
        s.pack("Bisection::lb", lb_);
        s.pack("Bisection::ub", ub_);
    }

    Bisection::Bisection(DeserializingStream &s) : Rootfinder(s)
    {
        s.version("Bisection", 1);
        s.unpack("Bisection::max_iter", max_iter_);
        s.unpack("Bisection::max_search", max_search_);
        s.unpack("Bisection::search_step", search_step_);
        s.unpack("Bisection::abstol", abstol_);
        s.unpack("Bisection::abstol_step", abstol_step_);
        s.unpack("Bisection::lb", lb_);
        s.unpack("Bisection::ub", ub_);
    }

}