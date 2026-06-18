#ifndef CASADI_CONOPT_INTERFACE_HPP
#define CASADI_CONOPT_INTERFACE_HPP

#include "casadi/core/nlpsol_impl.hpp"
#include <conopt.h>
#include <vector>
#include <string>
#include <utility>

namespace casadi {
  class ConoptInterface;

  struct ConoptMemory : public NlpsolMemory {
    const ConoptInterface& self;
    coiHandle_t cntvect;

    int modsta;
    int solsta;
    int iter;
    std::string return_status;
    int success;

    // Caching state for the evaluation block
    std::vector<double> cached_x;
    double cached_f;
    std::vector<double> cached_grad_f;
    std::vector<double> cached_g;
    std::vector<double> cached_jac_g;
    bool cache_valid;

    // Options handling
    std::vector<std::pair<std::string, GenericType>> custom_options;
    size_t current_option_idx;

    // Range-constraint expansion state (recomputed each solve)
    int ng_expanded;
    int numnz_expanded;
    std::vector<int> conopt_to_casadi;   // CONOPT constraint row (0-indexed) → CasADi index
    std::vector<int> casadi_to_conopt_lb_row;  // CasADi index → CONOPT row for lb (or only) side
    std::vector<int> casadi_to_conopt_ub_row;  // CasADi index → CONOPT row for ub side (range only); -1 otherwise
    std::vector<int> conopt_type;        // CONOPT TYPE for each expanded row
    std::vector<double> conopt_rhs;      // CONOPT RHS for each expanded row

    ConoptMemory(const ConoptInterface& interface);
    ~ConoptMemory();
  };

  class ConoptInterface : public Nlpsol {
  public:
    explicit ConoptInterface(const std::string& name, const Function& nlp);
    ~ConoptInterface() override;

    const char* plugin_name() const override { return "conopt"; }
    std::string class_name() const override { return "ConoptInterface"; }

    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new ConoptInterface(name, nlp);
    }

    static const Options options_;
    const Options& get_options() const override { return options_; }

    void init(const Dict& opts) override;
    void* alloc_mem() const override { return new ConoptMemory(*this); }
    int init_mem(void* mem) const override;
    void free_mem(void* mem) const override;
    void set_work(void* mem, const double**& arg, double**& res,
                  casadi_int*& iw, double*& w) const override;
    int solve(void* mem) const override;
    Dict get_stats(void* mem) const override;

    // Sparsities for the problem components
    Sparsity gradf_sp_;
    Sparsity jacg_sp_;
    Sparsity hesslag_sp_;
    bool exact_hessian_;
    Dict opts_; // CONOPT specific options

    // Per-column flag: true if the objective gradient has a nonzero in that column
    std::vector<bool> gradf_col_flag_;

    // Row-indexed (CSR) structure for fast Jacobian scatter in cb_fd_eval
    std::vector<int> jacg_rowstart_;
    std::vector<int> jacg_col_;
    std::vector<int> jacg_nzidx_;

    // Serialization and Deserialization
    void serialize_body(SerializingStream &s) const override;
    static ProtoFunction* deserialize(DeserializingStream& s) { return new ConoptInterface(s); }

    // CONOPT Mandatory Callbacks
    static int COI_CALLCONV cb_read_matrix(double LOWER[], double CURR[], double UPPER[], int VSTA[], int TYPEX[], double RHS[], int ESTA[], int COLSTA[], int ROWNO[], double VALUE[], int NLFLAG[], int NUMVAR, int NUMCON, int NUMNZ, void* USRMEM);
    static int COI_CALLCONV cb_fdevalini(const double X[], const int ROWLIST[], int MODE, int LISTSIZE, int NUMTHREAD, int IGNERR, int* ERRCNT, int NUMVAR, void* USRMEM);
    static int COI_CALLCONV cb_fd_eval(const double X[], double* G, double JAC[], int ROWNO, const int JACNUM[], int MODE, int IGNERR, int* ERRCNT, int NUMVAR, int NUMJAC, int THREAD, void* USRMEM);
    static int COI_CALLCONV cb_fdevalend(int IGNERR, int* ERRCNT, void* USRMEM);
    static int COI_CALLCONV cb_status(int MODSTA, int SOLSTA, int ITER, double OBJVAL, void* USRMEM);
    static int COI_CALLCONV cb_solution(const double XVAL[], const double XMAR[], const int XBAS[], const int XSTA[], const double YVAL[], const double YMAR[], const int YBAS[], const int YSTA[], int NUMVAR, int NUMCON, void* USRMEM);

    // Logging and Messages
    static int COI_CALLCONV cb_message(int SMSG, int DMSG, int NMSG, char* MSGV[], void* USRMEM);
    static int COI_CALLCONV cb_errmsg(int ROWNO, int COLNO, int POSNO, const char* MSG, void* USRMEM);

    // Options and Progress
    static int COI_CALLCONV cb_option(int NCALL, double* RVAL, int* IVAL, int* LVAL, char* NAME, void* USRMEM);
    static int COI_CALLCONV cb_progress(int LEN_INT, const int INTX[], int LEN_RL, const double RL[], const double X[], void* USRMEM);

    // CONOPT 2nd Order Callbacks
    static int COI_CALLCONV cb_2dlagrsize(int* NODRV, int NUMVAR, int NUMCON, int* NHESS, int MAXHESS, void* USRMEM);
    static int COI_CALLCONV cb_2dlagrstr(int HSRW[], int HSCL[], int* NODRV, int NUMVAR, int NUMCON, int NHESS, void* USRMEM);
    static int COI_CALLCONV cb_2dlagrval(const double X[], const double U[], const int HSRW[], const int HSCL[], double HSVL[], int* NODRV, int NUMVAR, int NUMCON, int NHESS, void* USRMEM);

  protected:
    explicit ConoptInterface(DeserializingStream& s);
  };
}
#endif // CASADI_CONOPT_INTERFACE_HPP
