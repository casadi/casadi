

#ifndef CASADI_BISECTION_HPP
#define CASADI_BISECTION_HPP

#include "casadi/core/rootfinder_impl.hpp"
#include <casadi/solvers/casadi_rootfinder_bisection_export.h>

/// \cond
namespace casadi
{
    struct CASADI_ROOTFINDER_BISECTION_EXPORT BisectionMemory : public RootfinderMemory
    {
        int return_status;
        casadi_int iter;
        casadi_int search_iter;
        double f_mid;
        double bracket_width;
    };

    class CASADI_ROOTFINDER_BISECTION_EXPORT Bisection : public Rootfinder
    {
    public:
        explicit Bisection(const std::string &name, const Function &f);

        ~Bisection() override;

        const char *plugin_name() const override { return "bisection"; }

        std::string class_name() const override { return "Bisection"; }

        static Rootfinder *creator(const std::string &name, const Function &f)
        {
            return new Bisection(name, f);
        }

        static const Options options_;

        const Options &get_options() const override { return options_; }

        Dict get_stats(void *mem) const override;

        void init(const Dict &opts) override;

        void *alloc_mem() const override { return new BisectionMemory(); }

        int init_mem(void *mem) const override;

        void free_mem(void *mem) const override { delete static_cast<BisectionMemory *>(mem); }

        void set_work(void *mem, const double **&arg, double **&res, casadi_int *&iw, double *&w) const override;

        int solve(void *mem) const override;

        static const std::string meta_doc;

        // void codegen_body(CodeGenerator &g) const override;

        // void codegen_declarations(CodeGenerator &g) const override;

        void serialize_body(SerializingStream &s) const override;

        static ProtoFunction *deserialize(DeserializingStream &s) { return new Bisection(s); }

    protected:
        explicit Bisection(DeserializingStream &s);

        casadi_int max_iter_;
        casadi_int max_search_;

        double abstol_;
        double abstol_step_;

        double search_step_;

        double lb_;
        double ub_;

        int finish(BisectionMemory *m, double x_sol, double f_sol, double width, int status, bool success, UnifiedReturnStatus urs) const;
        static std::string status_str(int status);
    };

}

/// \endcond
#endif // CASADI_BISECTION_HPP