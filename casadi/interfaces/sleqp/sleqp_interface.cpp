#include "sleqp_interface.hpp"


namespace casadi {
  extern "C"
  int CASADI_NLPSOL_SLEQP_EXPORT
  casadi_register_nlpsol_SLEQP(Nlpsol::Plugin* plugin) {
    plugin->creator = SLEQPInterface::creator;
    plugin->name = "SLEQP";
    plugin->doc = SLEQPInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &SLEQPInterface::options_;
    plugin->deserialize = &SLEQPInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_SLEQP_EXPORT casadi_load_nlpsol_SLEQP() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_SLEQP);
  }

  SLEQPInterface::SLEQPInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }

  SLEQPInterface::~SLEQPInterface() {
    clear_mem();
  }

  const Options SLEQPInterface::options_
  = {

  };

  const std::string SLEQPInterface::meta_doc = "";

  void SLEQPInterface::init(const Dict& opts) {
  }

  /** \brief Initalize memory block */
  int SLEQPInterface::init_mem(void* mem) const {
    return 0;
  }

  /// Get all statistics
  Dict SLEQPInterface::get_stats(void* mem) const {
  }

  /** \brief Set the (persistent) work vectors */
  void SLEQPInterface::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
  }

  // Solve the NLP
  int SLEQPInterface::solve(void* mem) const {
    //return calc_function(NULL, "nlp_f")==0;
    return 0;
  }
}
