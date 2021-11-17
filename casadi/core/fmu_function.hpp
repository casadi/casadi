/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_FMU_FUNCTION_IMPL_HPP
#define CASADI_FMU_FUNCTION_IMPL_HPP

#ifdef WITH_FMU

#include "function_internal.hpp"
#include "importer.hpp"
#include "casadi_enum.hpp"

//#include <fmi2FunctionTypes.h>
#include <fmi2Functions.h>
//#include <fmi2TypesPlatform.h>

/// \cond INTERNAL

namespace casadi {

// SYMBOL "jac_prob"
template<typename T1>
struct casadi_jac_prob {
  // Number of outputs, i.e. rows of the Jacobian
  casadi_int n_out;
  // Number of inputs, i.e. columns of the Jacobian
  casadi_int n_in;
  // Number of colors
  casadi_int n_color;
  // Extended Jacobian sparsity
  const casadi_int* sp_ext;
  // Jacobian coloring
  const casadi_int* coloring;
  // Nominal values for inputs, if any
  const T1* nom_in;
  // Index mapping for outputs (i.e. Jacobian rows), if any
  const size_t* map_out;
  // Index mapping for inputs (i.e.  Jacobian columns), if any
  const size_t* map_in;
};

// SYMBOL "jac_setup"
template<typename T1>
void casadi_jac_setup(casadi_jac_prob<T1>* p, const casadi_int* sp_ext,
    const casadi_int* coloring) {
  // Set pointers
  p->sp_ext = sp_ext;
  p->coloring = coloring;
  // Dimensions are given by the sparsity patterns
  p->n_out = sp_ext[0];
  p->n_in = sp_ext[1];
  p->n_color = coloring[1];
  // The following defaults to null
  p->nom_in = 0;
  p->map_out = 0;
  p->map_in = 0;
}

// SYMBOL "jac_work"
template<typename T1>
void casadi_jac_work(const casadi_jac_prob<T1>* p, casadi_int* sz_iw, casadi_int* sz_w) {
  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  // Work vectors in data struct
  *sz_iw += p->n_in;  // iseed
  *sz_w += p->n_in;  // seed
  *sz_iw += p->n_out;  // isens
  *sz_w += p->n_out;  // sens
  *sz_w += p->n_out;  // scal
  *sz_iw += p->n_out;  // wrt
  *sz_iw += p->n_out;  // nzind
}

// SYMBOL "jac_data"
template<typename T1>
struct casadi_jac_data {
  // Number of seeds, sensitivities for the current color
  casadi_int nseed, nsens;
  // Inputs that are being seeded
  casadi_int *iseed;
  // Set of seeds for the seeded inputs
  T1 *seed;
  // Set of outputs for which sensitivities are calculated
  casadi_int *isens;
  // Set of values for the calculated sensitivities
  T1 *sens;
  // Scaling factors for calculated sensitivities
  T1 *scal;
  // Input corresponding to calculated sensitivities
  casadi_int *wrt;
  // Jacobian nonzero corresponding to calculated sensitivities
  casadi_int *nzind;
};

// SYMBOL "jac_init"
template<typename T1>
void casadi_jac_init(const casadi_jac_prob<T1>* p, casadi_jac_data<T1>* d,
    casadi_int** iw, T1** w) {
  // Set work vectors
  d->iseed = *iw; *iw += p->n_in;
  d->seed = *w; *w += p->n_in;
  d->isens = *iw; *iw += p->n_out;
  d->sens = *w; *w += p->n_out;
  d->scal = *w; *w += p->n_out;
  d->wrt = *iw; *iw += p->n_out;
  d->nzind = *iw; *iw += p->n_out;
}

// SYMBOL "jac_pre"
template<typename T1>
void casadi_jac_pre(const casadi_jac_prob<T1>* p, casadi_jac_data<T1>* d, casadi_int c) {
  // Local variables
  casadi_int i, kc, vin, vout, Jk;
  double nom, inv_nom;
  const casadi_int *color_colind, *color_row, *jac_colind, *jac_row;
  // Extract sparsities
  color_colind = p->coloring + 2;
  color_row = color_colind + p->n_color + 1;
  jac_colind = p->sp_ext + 2;
  jac_row = jac_colind + p->n_in + 1;
  // Loop over input indices for color
  d->nseed = d->nsens = 0;
  for (kc = color_colind[c]; kc < color_colind[c + 1]; ++kc) {
    vin = color_row[kc];
    // Nominal value, used as a seed for the column
    nom = p->nom_in ? p->nom_in[vin] : 1;
    inv_nom = 1. / nom;
    // Collect seeds for column
    d->seed[d->nseed] = nom;
    d->iseed[d->nseed] = vin;
    d->nseed++;
    // Request corresponding outputs
    for (Jk = jac_colind[vin]; Jk < jac_colind[vin + 1]; ++Jk) {
      vout = jac_row[Jk];
      d->scal[d->nsens] = inv_nom;
      d->isens[d->nsens] = vout;
      d->wrt[d->nsens] = vin;
      d->nzind[d->nsens] = Jk;
      d->nsens++;
    }
  }
  // Map indices
  if (p->map_in) {
    for (i = 0; i < d->nseed; ++i) d->iseed[i] = p->map_in[d->iseed[i]];
    for (i = 0; i < d->nsens; ++i) d->wrt[i] = p->map_in[d->wrt[i]];
  }
  if (p->map_out) {
    for (i = 0; i < d->nsens; ++i) d->isens[i] = p->map_out[d->isens[i]];
  }
}

// SYMBOL "jac_scale"
template<typename T1>
void casadi_jac_scale(const casadi_jac_prob<T1>* p, casadi_jac_data<T1>* d) {
  // Local variables
  casadi_int i;
  // Scale derivatives
  for (i = 0; i < d->nsens; ++i) d->sens[i] *= d->scal[i];
}

// SYMBOL "get_sub"
template<typename T1>
void casadi_get_sub(T1* sub, const casadi_int* sp_a, const T1* nz_a,
    casadi_int rbegin, casadi_int rend, casadi_int cbegin, casadi_int cend) {
  // Local variables
  casadi_int nr, nc, r, c, k;
  const casadi_int *colind, *row;
  // Quick return if null
  if (sub == 0) return;
  // Extract sparsity
  nr = sp_a[0];
  nc = sp_a[1];
  colind = sp_a + 2;
  row = colind + nc + 1;
  // Loop over columns
  for (c = cbegin; c < cend; ++c) {
    // Loop over nonzeros
    for (k = colind[c]; k < colind[c + 1]; ++k) {
      // Get row
      r = row[k];
      // Save, if in range
      if (r >= rbegin && r < rend) *sub++ = nz_a[k];
    }
  }
}

// Forward declarations
class DaeBuilderInternal;
class FmuFunction;

// Memory object
struct CASADI_EXPORT FmuMemory : public FunctionMemory {
  // Function object
  const FmuFunction& self;
  // Evaluation arguments, work vectors
  const double** arg;
  double** res;
  casadi_int* iw;
  double* w;
  // Extended Jacobian
  double *jac_nz;
  // Adjoint seeds, sensitivities being calculated
  double *aseed, *asens;
  // Component memory
  fmi2Component c;
  // Additional (slave) memory objects
  std::vector<FmuMemory*> slaves;
  // Input and output buffers
  std::vector<double> ibuf_, obuf_;
  // Seeds, sensitivities
  std::vector<double> seed_, sens_;
  // Which entries have been changed
  std::vector<bool> changed_;
  // Which entries are being requested
  std::vector<bool> requested_;
  // Derivative with respect to
  std::vector<size_t> wrt_;
  // Current known/unknown variables
  std::vector<size_t> id_in_, id_out_;
  // Which perturbations are permitted
  std::vector<bool> in_bounds_;
  // Value references
  std::vector<fmi2ValueReference> vr_in_, vr_out_;
  // Work vector (reals)
  std::vector<fmi2Real> v_in_, v_out_, d_in_, d_out_, fd_out_, v_pert_;
  // Constructor
  explicit FmuMemory(const FmuFunction& self) : self(self), c(nullptr) {}
};

/// Variable type
enum class FdMode {FORWARD, BACKWARD, CENTRAL, SMOOTHING, NUMEL};

/// Convert to string
CASADI_EXPORT std::string to_string(FdMode v);

/// Type of parallelization
enum class Parallelization {SERIAL, OPENMP, THREAD, NUMEL};

/// Convert to string
CASADI_EXPORT std::string to_string(Parallelization v);

// Interface to binary FMU (shared between derivatives)
struct CASADI_EXPORT Fmu {
  // Constructor
  Fmu(const std::vector<std::string>& name_in, const std::vector<std::string>& name_out,
    const std::map<std::string, std::vector<size_t>>& scheme,
    const std::map<std::string, std::vector<size_t>>& lc);

  // Initialize
  void init(const DaeBuilderInternal* dae);

  // Reference counter
  int counter_;

  // DLL
  Importer li_;

  // IO scheme, linear combinations
  std::vector<std::string> name_in_, name_out_;
  std::map<std::string, std::vector<size_t>> scheme_, lc_;

  // Mapping from scheme variable to and from FMU variable indices
  std::vector<size_t> iind_, iind_map_, oind_, oind_map_;

  // Meta information about the input/output variable subsets
  std::vector<double> nominal_in_, nominal_out_;
  std::vector<double> min_in_, min_out_;
  std::vector<double> max_in_, max_out_;
  std::vector<std::string> vn_in_, vn_out_;
  std::vector<fmi2ValueReference> vr_in_, vr_out_;

  // Reduced space indices for all inputs and outputs
  std::vector<std::vector<size_t>> ired_, ored_;

  // FMU C API function prototypes. Cf. FMI specification 2.0.2
  fmi2InstantiateTYPE* instantiate_;
  fmi2FreeInstanceTYPE* free_instance_;
  fmi2ResetTYPE* reset_;
  fmi2SetupExperimentTYPE* setup_experiment_;
  fmi2EnterInitializationModeTYPE* enter_initialization_mode_;
  fmi2ExitInitializationModeTYPE* exit_initialization_mode_;
  fmi2EnterContinuousTimeModeTYPE* enter_continuous_time_mode_;
  fmi2GetRealTYPE* get_real_;
  fmi2SetRealTYPE* set_real_;
  fmi2GetBooleanTYPE* get_boolean_;
  fmi2SetBooleanTYPE* set_boolean_;
  fmi2GetIntegerTYPE* get_integer_;
  fmi2SetIntegerTYPE* set_integer_;
  fmi2GetStringTYPE* get_string_;
  fmi2SetStringTYPE* set_string_;
  fmi2GetDirectionalDerivativeTYPE* get_directional_derivative_;

  // Callback functions
  fmi2CallbackFunctions functions_;

  // Path to the FMU resource directory
  std::string resource_loc_;

  // Tolerance
  double fmutol_;

  // Instance name
  std::string instance_name_;

  // GUID
  std::string guid_;

  // Logging?
  bool logging_on_;

  // Collection of variable values, all types
  struct Value {
    std::vector<fmi2Real> v_real;
    std::vector<fmi2Integer> v_integer;
    std::vector<fmi2Boolean> v_boolean;
    std::vector<std::string> v_string;
  };

  // Variables used for initialization, by type
  std::vector<fmi2ValueReference> vr_real_, vr_integer_, vr_boolean_, vr_string_;
  std::vector<fmi2Real> init_real_;
  std::vector<fmi2Integer> init_integer_;
  std::vector<fmi2Boolean> init_boolean_;
  std::vector<std::string> init_string_;


  // Auxilliary variables, by type
  std::vector<std::string> vn_aux_real_, vn_aux_integer_, vn_aux_boolean_, vn_aux_string_;
  std::vector<fmi2ValueReference> vr_aux_real_, vr_aux_integer_, vr_aux_boolean_, vr_aux_string_;
  Value aux_;

  // Sparsity pattern for Jacobian of all outputs w.r.t. all inputs
  Sparsity sp_jac_;

  // Graph coloring for sp_jac_
  Sparsity coloring_jac_;

  // Index lookup for input
  size_t index_in(const std::string& n) const;

  // Index lookup for output
  size_t index_out(const std::string& n) const;

  // Load an FMI function
  signal_t load_function(const std::string& symname);

  // New memory object
  fmi2Component instantiate() const;

  // Free FMU instance
  void free_instance(fmi2Component c) const;

  // Reset solver
  int reset(fmi2Component c);

  // Setup experiment
  void setup_experiment(fmi2Component c) const;

  // Enter initialization mode
  int enter_initialization_mode(fmi2Component c) const;

  // Exit initialization mode
  int exit_initialization_mode(fmi2Component c) const;

  // Copy values set in DaeBuilder to FMU
  int set_values(fmi2Component c) const;

  // Retrieve auxilliary variables from FMU
  int get_aux(fmi2Component c, Value* v) const;

  /** \brief Get stats */
  void get_stats(FmuMemory* m, Dict* stats) const;

  /** \brief Initalize memory block */
  int init_mem(FmuMemory* m) const;

  // Set value
  void set(FmuMemory* m, size_t ind, const double* value) const;

  // Request the calculation of a variable
  void request(FmuMemory* m, size_t ind) const;

  // Gather user inputs and outputs
  void gather_io(FmuMemory* m) const;

  // Calculate all requested variables
  int eval(FmuMemory* m) const;

  // Get a calculated variable
  void get(FmuMemory* m, size_t id, double* value) const;

  // Set seed
  void set_seed(FmuMemory* m, casadi_int nseed, const casadi_int* id, const double* v) const;

  // Request the calculation of a sensitivity
  void request_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const;

  // Get calculated derivatives
  void get_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id, double* v) const;

  // Gather user sensitivities
  void gather_sens(FmuMemory* m) const;

  // Calculate directional derivatives using AD
  int eval_ad(FmuMemory* m) const;

  // Calculate directional derivatives using FD
  int eval_fd(FmuMemory* m) const;

  // Calculate directional derivatives
  int eval_derivative(FmuMemory* m) const;

  // Get Jacobian sparsity for a subset of inputs and outputs
  Sparsity jac_sparsity(const std::vector<size_t>& osub, const std::vector<size_t>& isub) const;

  // Get Jacobian sparsity for an output/input pair
  Sparsity jac_sparsity(size_t oind, size_t iind) const {
    return jac_sparsity(ored_.at(oind), ired_.at(iind));
  }

  // Name of system, per the FMI specification
  static std::string system_infix();

  // DLL suffix, per the FMI specification
  static std::string dll_suffix();

  // Process message
  static void logger(fmi2ComponentEnvironment componentEnvironment,
    fmi2String instanceName,
    fmi2Status status,
    fmi2String category,
    fmi2String message, ...);
};

class CASADI_EXPORT FmuFunction : public FunctionInternal {
 public:
  // FMU (shared between derivative expressions
  Fmu* fmu_;

  // Types of inputs
  enum InputType {REG_INPUT, ADJ_SEED, DUMMY_OUTPUT};

  // Input structure
  struct InputStruct {
    // Type of input
    InputType type;
    // Corresponding index in Fmu
    size_t ind;
  };

  // Information about function inputs
  std::vector<InputStruct> in_;

  // Types of inputs
  enum OutputType {REG_OUTPUT, ADJ_SENS, JAC_OUTPUT};

  // Output structure
  struct OutputStruct {
    // Type of input
    OutputType type;
    // Output index in Fmu
    size_t ind;
    // With-respect-to index in Fmu
    size_t wrt;
    // Selection
    size_t rbegin, rend, cbegin, cend;
    // Constructor
    OutputStruct() : ind(-1), wrt(-1), rbegin(-1), rend(-1), cbegin(-1), cend(-1) {}
  };

  // Information about function outputs
  std::vector<OutputStruct> out_;

  // All Jacobian inputs and outputs
  std::vector<size_t> jac_in_, jac_out_;

  // Nominal values for Jacobian inputs
  std::vector<double> jac_nom_in_;

  // What blocks exist?
  bool has_jac_, has_adj_;

  // User-set options
  bool enable_ad_, validate_ad_;
  double step_, abstol_, reltol_, u_aim_, h_min_, h_max_;
  casadi_int h_iter_;

  // FD method as an enum
  FdMode fd_;

  // Types of parallelization
  Parallelization parallelization_;

  // Stats from initialization
  Dict init_stats_;

  /** \brief Constructor */
  FmuFunction(const std::string& name, Fmu* fmu,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out);

  /** \brief Destructor */
  ~FmuFunction() override;

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunction";}

  ///@{
  /** \brief Options */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /// Initialize
  void init(const Dict& opts) override;

  // Parse input name
  void parse_input(InputStruct* s, const std::string& n) const;

  // Parse output name
  void parse_output(OutputStruct* s, const std::string& n) const;

  // Get sparsity pattern for extended Jacobian
  Sparsity sp_ext_;

  // Graph coloring
  Sparsity coloring_;

  // Jacobian memory
  casadi_jac_prob<double> p_;

  // Work vector size for Jacobian calculation
  casadi_int jac_iw_, jac_w_;

  // Number of parallel tasks
  casadi_int max_n_task_;

  ///@{
  /** \brief Number of function inputs and outputs */
  size_t get_n_in() override { return in_.size();}
  size_t get_n_out() override {return out_.size();}
  ///@}

  /// @{
  /** \brief Retreive sparsities */
  Sparsity get_sparsity_in(casadi_int i) override;
  Sparsity get_sparsity_out(casadi_int i) override;
  /// @}

  // Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

  // Evaluate numerically, single thread
  int eval_thread(FmuMemory* m, casadi_int thread, casadi_int n_thread,
    bool need_jac, bool need_adj) const;

  ///@{
  /** \brief Return sparsity of Jacobian of an output respect to an input */
  bool has_jac_sparsity(casadi_int oind, casadi_int iind) const override;
  Sparsity get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const override;
  ///@}

  ///@{
  /** \brief Full Jacobian */
  bool has_jacobian() const override {return true;}
  Function get_jacobian(const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  ///@{
  /** \brief Reverse mode AD */
  bool has_reverse(casadi_int nadj) const override { return nadj == 1;}
  Function get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  // Helper function
  static bool has_prefix(const std::string& s);

  // Split prefix
  static std::string pop_prefix(const std::string& s, std::string* rem = 0);

  /** \brief Create memory block */
  void* alloc_mem() const override;

  /** \brief Initalize memory block */
  int init_mem(void* mem) const override;

  /** \brief Free memory block */
  void free_mem(void *mem) const override;

  /// Get all statistics
  Dict get_stats(void* mem) const override;
};

} // namespace casadi
/// \endcond

#endif  // WITH_FMU

#endif // CASADI_FMU_FUNCTION_IMPL_HPP
