#include <sleqp.h>

#include <casadi/interfaces/sleqp/casadi_nlpsol_sleqp_export.h>
#include "casadi/core/nlpsol_impl.hpp"
#include "casadi/core/timing.hpp"

#include "sleqp.h"

// TODO: Use casadi exceptions / error reporting??
#define SLEQP_CALL_EXC(x)                       \
  do {                                          \
    const SLEQP_RETCODE _status = (x);          \
    if(_status != SLEQP_OKAY) {                 \
      throw std::runtime_error("SLEQP error");  \
    }                                           \
  } while(false)

namespace casadi {

  class SLEQPInterface;

  struct CASADI_NLPSOL_SLEQP_EXPORT SLEQPMemory : public NlpsolMemory {

    struct {
      SleqpProblem* problem;
      SleqpSettings* settings;

      SleqpVec* primal;
      SleqpSolver* solver;
    } internal;

    // Current calculated quantities
    double* xk;
    double *gk, *grad_fk, *jac_gk, *hess_lk, *grad_lk;

    // Additional callback data
    double *cb_xk, *cb_lam_xk, *cb_lam_gk;

    bool iteration_callback_ignore_errors;

    const SLEQPInterface* interface;
  };

  class CASADI_NLPSOL_SLEQP_EXPORT SLEQPInterface : public Nlpsol {
  private:
    void clear_mem_at(SLEQPMemory* m) const;

  public:
    explicit SLEQPInterface(const std::string& name, const Function& nlp);
    ~SLEQPInterface() override;

    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new SLEQPInterface(name, nlp);
    }

    const char* plugin_name() const override { return "sleqp";}

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

    Sparsity jacg_sp_;

    /// All SLEQP options
    Dict opts_;

    int max_iter_;
    double max_wall_time_;
    int print_level_;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new SLEQPInterface(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit SLEQPInterface(DeserializingStream& s);
  };
};
