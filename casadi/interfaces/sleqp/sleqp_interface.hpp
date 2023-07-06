#include <sleqp.h>

#include <casadi/interfaces/sleqp/casadi_nlpsol_sleqp_export.h>
#include "casadi/core/nlpsol_impl.hpp"
#include "casadi/core/timing.hpp"
#include "sleqp/pub_problem.h"
#include "sleqp/pub_settings.h"
#include "sleqp/pub_solver.h"

namespace casadi {

  struct CASADI_NLPSOL_SLEQP_EXPORT SLEQPMemory : public NlpsolMemory {
    SleqpProblem* problem;
    SleqpSettings* settings;

    SleqpVec* primal;
    SleqpSolver* solver;
  };

  class CASADI_NLPSOL_SLEQP_EXPORT SLEQPInterface : public Nlpsol {
  private:
    void clear_mem(SLEQPMemory* m) const;

  public:
    explicit SLEQPInterface(const std::string& name, const Function& nlp);
    ~SLEQPInterface() override;

    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new SLEQPInterface(name, nlp);
    }

    const char* plugin_name() const override { return "SLEQP";}

    // Get name of the class
    std::string class_name() const override { return "SLEQPInterface";}

    static const Options options_;
    static const std::string meta_doc;

    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new SLEQPMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                  casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

  };
};
