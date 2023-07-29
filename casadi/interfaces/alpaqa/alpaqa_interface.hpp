#include <alpaqa/panoc-alm.hpp>
#include <alpaqa/problem/box-constr-problem.hpp>
#include <alpaqa/problem/problem-with-counters.hpp>

#include <casadi/interfaces/alpaqa/casadi_nlpsol_alpaqa_export.h>
#include "casadi/core/nlpsol_impl.hpp"
#include "casadi/core/timing.hpp"
#include "alpaqa_problem.hpp"

using Direction   = alpaqa::LBFGSDirection<alpaqa::DefaultConfig>;
using InnerSolver = alpaqa::PANOCSolver<Direction>;
using OuterSolver = alpaqa::ALMSolver<InnerSolver>;

namespace casadi {

  class AlpaqaInterface;

  struct CASADI_NLPSOL_ALPAQA_EXPORT AlpaqaMemory : public NlpsolMemory {

    alpaqa::ALMSolver<InnerSolver>* solver;

    const AlpaqaInterface* interface;
  };

  class CASADI_NLPSOL_ALPAQA_EXPORT AlpaqaInterface : public Nlpsol {
    friend class AlpaqaProblem;
  private:
    void clear_mem_at(AlpaqaMemory* m) const;

  public:
    explicit AlpaqaInterface(const std::string& name, const Function& nlp);
    ~AlpaqaInterface() override;

    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new AlpaqaInterface(name, nlp);
    }

    const char* plugin_name() const override { return "alpaqa";}

    // Get name of the class
    std::string class_name() const override { return "AlpaqaInterface";}

    static const Options options_;
    static const std::string meta_doc;

    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new AlpaqaMemory();}

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

    /// All Alpaqa options
    Dict opts_;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new AlpaqaInterface(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit AlpaqaInterface(DeserializingStream& s);

    
  };
};
