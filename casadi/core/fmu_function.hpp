/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

#include <fmi2Functions.h>

/// \cond INTERNAL

namespace casadi {

template<typename T1>
void casadi_forward_diff_new(const T1* yk, T1* J, T1 h, casadi_int n_y) {
  // Local variables
  casadi_int i;
  const double *yf, *yc;
  T1 hinv;
  // Inverse of step size
  hinv = 1. / h;
  // Get stencil
  yc = yk;
  yf = yk + n_y;
  // Calculate FD approximation
  for (i = 0; i < n_y; ++i) J[i] = hinv * (yf[i] - yc[i]);
}

template<typename T1>
void casadi_central_diff_new(const T1* yk, T1* J, T1 h, casadi_int n_y) {
  // Local variables
  casadi_int i;
  const T1 *yf, *yc, *yb;
  T1 hinv;
  // Inverse of step size
  hinv = 1. / h;
  // Get stencil
  yb = yk;
  yc = yk + n_y;
  yf = yk + 2 * n_y;
  // Set u and stencils to zero (also supresses warnings)
  for (i = 0; i < n_y; ++i) {
    if (isfinite(yb[i])) {
      if (isfinite(yf[i])) {
        // Both forward and backward allowed
        J[i] = 0.5 * hinv * (yf[i] - yb[i]);
      } else {
        // Backward but not forward allowed
        J[i] = hinv * (yc[i] - yb[i]);
      }
    } else {
      if (isfinite(yf[i])) {
        // Forward but not backward allowed
        J[i] = hinv * (yf[i] - yc[i]);
      } else {
        // Neither forward nor backward possible
        J[i] = std::numeric_limits<T1>::quiet_NaN();
      }
    }
  }
}

template<typename T1>
T1 casadi_central_diff_err(const T1* yk, T1 h, casadi_int n_y, casadi_int i,
    T1 abstol, T1 reltol) {
  // Local variables
  const T1 *yf, *yc, *yb;
  T1 err_trunc, err_round;
  // Get stencil
  yb = yk;
  yc = yk + n_y;
  yf = yk + 2 * n_y;
  // Only consider points where both forward and backward allowed
  if (isfinite(yb[i]) && isfinite(yf[i])) {
    // Truncation error
    err_trunc = yf[i] - 2*yc[i] + yb[i];
    // Roundoff error
    err_round = reltol / h * fmax(fabs(yf[i] - yc[i]), fabs(yc[i] - yb[i])) + abstol;
    // Error quotient estimate
    return err_trunc / err_round;
  } else {
    // Cannot be calculated
    return std::numeric_limits<T1>::quiet_NaN();;
  }
}

template<typename T1>
T1 casadi_smoothing_diff_weights(casadi_int k, T1 yb, T1 yc, T1 yf, T1 *J) {
  // Calculate shifted finite difference approximation, weights
  if (k == 0) {
    // Backward shifted
    // 7.10 in Conte & Carl de Boor: Elementary Numerical Analysis (1972)
    // and 25.3.4 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
    if (J) *J = 3 * yf - 4 * yc + yb;
    // Relative weight is 1
    return 1;
  } else if (k == 1) {
    // Central
    // We give this the relative weight 4 since if all weights are equal,
    // this would amount to a five point formula for the derivative
    // (yb2 - 8*yb + 8*yf - y_f2)/(12*h)
    // cf. 25.3.6 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
    if (J) *J = yf - yb;
    // Relative weight is 4
    return 4;
  } else {
    // Forward shifted, cf. backward shifted above
    if (J) *J = -3 * yb + 4 * yc - yf;
    // Relative weight is 1
    return 1;
  }
}

template<typename T1>
void casadi_smoothing_diff_new(const T1* yk, T1* J, T1 h, casadi_int n_y, T1 smoothing) {
  // Stencil
  T1 yb, yc, yf;
  // Local variables
  T1 Jk, wk, sw, ui, sm;
  casadi_int i, k;
  // Set stencils to zero (also supresses warnings)
  yf = yc = yb = 0;
  for (i=0; i<n_y; ++i) {
    // Reset derivative estimate, sum of weights, error estimate
    J[i] = sw = ui = 0;
    // For backward shifted, central and forward shifted
    for (k = 0; k < 3; ++k) {
      // Get stencil
      yb = yk[i + k * n_y];
      yc = yk[i + (k + 1) * n_y];
      yf = yk[i + (k + 2) * n_y];
      // No contribuation if any value is infinite
      if (!isfinite(yb) || !isfinite(yc) || !isfinite(yf)) continue;
      // Calculate weights
      wk = casadi_smoothing_diff_weights(k, yb, yc, yf, &Jk);
      // Smoothness measure (second order derivative)
      sm = yf - 2*yc + yb;
      sm /= h*h;
      // Modify the weight according to smoothness
      wk /= sm*sm + smoothing;
      sw += wk;
      // Added weighted contribution to weight and error
      J[i] += wk * Jk;
    }
    // If sw is 0, no stencil worked
    if (sw == 0) {
      // Set component to 0, return -1
      J[i] = std::numeric_limits<T1>::quiet_NaN();
    } else {
      // Finalize estimate using the sum of weights and the step length
      J[i] /= 2*h*sw;
    }
  }
}

template<typename T1>
T1 casadi_smoothing_diff_err(const T1* yk, T1 h, casadi_int n_y, casadi_int i,
    T1 abstol, T1 reltol, T1 smoothing) {
  // Stencil
  T1 yb, yc, yf;
  // Local variables
  T1 wk, sw, ui, err_trunc, err_round, sm;
  casadi_int k;
  // Set stencils to zero (also supresses warnings)
  yf = yc = yb = 0;
  // Reset derivative estimate, sum of weights, error estimate
  sw = ui = 0;
  // For backward shifted, central and forward shifted
  for (k = 0; k < 3; ++k) {
    // Get stencil
    yb = yk[i + k * n_y];
    yc = yk[i + (k + 1) * n_y];
    yf = yk[i + (k + 2) * n_y];
    // No contribuation if any value is infinite
    if (!isfinite(yb) || !isfinite(yc) || !isfinite(yf)) continue;
    // Calculate weights
    wk = casadi_smoothing_diff_weights(k, yb, yc, yf, static_cast<T1*>(0));
    // Truncation error
    err_trunc = yf - 2*yc + yb;
    // Roundoff error
    err_round = reltol/h*fmax(fabs(yf - yc), fabs(yc - yb)) + abstol;
    // We use the second order derivative as a smoothness measure
    sm = err_trunc/(h*h);
    // Modify the weight according to smoothness
    wk /= sm*sm + smoothing;
    sw += wk;
    // Added weighted contribution to weight and error
    ui += wk * fabs(err_trunc / err_round);
  }
  // If sw is 0, no stencil worked
  if (sw == 0) {
    // Cannot be calculated
    return std::numeric_limits<T1>::quiet_NaN();;
  } else {
    // Finalize estimate using the sum of weights and the step length
    return ui / sw;
  }
}

// Forward declarations
class DaeBuilderInternal;
class FmuFunction;
struct InputStruct;

// Memory object
struct CASADI_EXPORT FmuMemory : public FunctionMemory {
  // Function object
  const FmuFunction& self;
  // Evaluation inputs
  const double** arg;
  // Evaluation outputs
  double** res;
  // Work vector for star coloring
  casadi_int* star_iw;
  // Extended Jacobian
  double *jac_nz;
  // Extended Hessian
  double *hess_nz;
  // Adjoint seeds, sensitivities being calculated
  double *aseed, *asens, *pert_asens;
  // Memory for Jacobian calculation
  casadi_jac_data<double> d;
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
  // Flip sign?
  std::vector<bool> flip_;
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

/// Length of FD stencil, including unperturbed input
CASADI_EXPORT casadi_int n_fd_points(FdMode v);

/// Offset for FD stencil, i.e. index of unperturbed input
CASADI_EXPORT casadi_int fd_offset(FdMode v);

/// Calculate FD estimate
template<typename T1>
CASADI_EXPORT void finite_diff(FdMode v, const T1* yk, T1* J, T1 h, casadi_int n_y,
    T1 smoothing) {
  switch (v) {
    case FdMode::FORWARD:
    case FdMode::BACKWARD:
      return casadi_forward_diff_new(yk, J, h, n_y);
    case FdMode::CENTRAL:
      return casadi_central_diff_new(yk, J, h, n_y);
    case FdMode::SMOOTHING:
      return casadi_smoothing_diff_new(yk, J, h, n_y, eps);
    default:
      casadi_error("FD mode " + to_string(v) + " not implemented");
  }
}

/// Is an error estimate available
CASADI_EXPORT bool fd_has_err(FdMode v);

/// Calculate FD error estimate
template<typename T1>
CASADI_EXPORT T1 finite_diff_err(FdMode v, const T1* yk, T1 h, casadi_int n_y, casadi_int i,
    T1 abstol, T1 reltol, T1 smoothing) {
  switch (v) {
    case FdMode::CENTRAL:
      return casadi_central_diff_err(yk, h, n_y, i, abstol, reltol);
    case FdMode::SMOOTHING:
      return casadi_smoothing_diff_err(yk, h, n_y, i, abstol, reltol, smoothing);
    default:
      casadi_error("Error estimate cannot be calculated for FD mode " + to_string(v));
  }
  return -1;
}

/// Type of parallelization
enum class Parallelization {SERIAL, OPENMP, THREAD, NUMEL};

/// Convert to string
CASADI_EXPORT std::string to_string(Parallelization v);

// Interface to binary FMU (shared between derivatives)
struct CASADI_EXPORT Fmu {
  // Constructor
  Fmu(const std::vector<std::string>& scheme_in, const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme, const std::vector<std::string>& aux);

  // Initialize
  void init(const DaeBuilderInternal* dae);

  // Reference counter
  int counter_;

  // DLL
  Importer li_;

  // IO scheme
  std::vector<std::string> scheme_in_, scheme_out_;
  std::map<std::string, std::vector<size_t>> scheme_;

  // Auxilliary outputs
  std::vector<std::string> aux_;

  // Mapping from scheme variable to and from FMU variable indices
  std::vector<size_t> iind_, iind_map_, oind_, oind_map_;

  // Meta information about the input/output variable subsets
  std::vector<double> nominal_in_, nominal_out_;
  std::vector<double> min_in_, min_out_;
  std::vector<double> max_in_, max_out_;
  std::vector<std::string> vn_in_, vn_out_;
  std::vector<fmi2ValueReference> vr_in_, vr_out_;

  // Numerical values for inputs
  std::vector<fmi2Real> value_in_;

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
  Value aux_value_;

  // Sparsity pattern for extended Jacobian, Hessian
  Sparsity jac_sp_, hess_sp_;

  // Index lookup for input
  size_t index_in(const std::string& n) const;

  // Index lookup for output
  size_t index_out(const std::string& n) const;

  // Load an FMI function
  template<typename T>
  T* load_function(const std::string& symname);

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

  // Retrieve input variable values from FMU
  int get_in(fmi2Component c, std::vector<fmi2Real>* v) const;

  // Retrieve auxilliary variables from FMU
  int get_aux(fmi2Component c, Value* v) const;

  /** \brief Get stats

      \identifier{ye} */
  void get_stats(FmuMemory* m, Dict* stats,
    const std::vector<std::string>& name_in, const InputStruct* in) const;

  /** \brief Initalize memory block

      \identifier{yf} */
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
  int eval_fd(FmuMemory* m, bool independent_seeds) const;

  // Calculate directional derivatives
  int eval_derivative(FmuMemory* m, bool independent_seeds) const;

  // Get Jacobian sparsity for a subset of inputs and outputs
  Sparsity jac_sparsity(const std::vector<size_t>& osub, const std::vector<size_t>& isub) const;

  // Get Jacobian sparsity for an output/input pair
  Sparsity jac_sparsity(size_t oind, size_t iind) const {
    return jac_sparsity(ored_.at(oind), ired_.at(iind));
  }

  // Get Hessian sparsity for a subset of inputs
  Sparsity hess_sparsity(const std::vector<size_t>& r, const std::vector<size_t>& c) const;

  // Get Jacobian sparsity for an input/input pair
  Sparsity hess_sparsity(size_t r, size_t c) const {
    return hess_sparsity(ired_.at(r), ired_.at(c));
  }

  /// @{
  /** \brief Retreive nominal values

      \identifier{yg} */
  std::vector<double> get_nominal_in(casadi_int i) const;
  std::vector<double> get_nominal_out(casadi_int i) const;
  /// @}

  // Print description of an input
  std::string desc_in(FmuMemory* m, size_t id) const;

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

// Types of inputs
enum class InputType {REG, FWD, ADJ, OUT, ADJ_OUT};

// Input structure
struct CASADI_EXPORT InputStruct {
  // Type of input
  InputType type;
  // Corresponding index in Fmu
  size_t ind;
  // Parse an input string
  static InputStruct parse(const std::string& n, const Fmu* fmu,
    std::vector<std::string>* name_in = 0,
    std::vector<std::string>* name_out = 0);
};

// Types of inputs
enum class OutputType {REG, FWD, ADJ, JAC, JAC_TRANS, JAC_ADJ_OUT, JAC_REG_ADJ, HESS};

// Output structure
struct CASADI_EXPORT OutputStruct {
  // Type of input
  OutputType type;
  // Output index in Fmu
  size_t ind;
  // With-respect-to index in Fmu
  size_t wrt;
  // Selection
  size_t rbegin, rend, cbegin, cend;
  // Parse an output string
  static OutputStruct parse(const std::string& n, const Fmu* fmu,
    std::vector<std::string>* name_in = 0,
    std::vector<std::string>* name_out = 0);
  // Constructor
  OutputStruct() : ind(-1), wrt(-1), rbegin(-1), rend(-1), cbegin(-1), cend(-1) {}
};

// Helper function
CASADI_EXPORT bool has_prefix(const std::string& s);

// Split prefix
CASADI_EXPORT std::string pop_prefix(const std::string& s, std::string* rem = 0);

class CASADI_EXPORT FmuFunction : public FunctionInternal {
 public:
  // FMU (shared between derivative expressions
  Fmu* fmu_;

  // Information about function inputs
  std::vector<InputStruct> in_;

  // Information about function outputs
  std::vector<OutputStruct> out_;

  // All Jacobian inputs and outputs
  std::vector<size_t> jac_in_, jac_out_;

  // Nominal values for Jacobian inputs
  std::vector<double> jac_nom_in_;

  // Sparsity of transpose (if needed)
  std::vector<Sparsity> sp_trans_;
  std::vector<casadi_int> sp_trans_map_;

  // What blocks exist?
  bool has_jac_, has_fwd_, has_adj_, has_hess_;

  // User-set options
  bool enable_ad_, validate_ad_, make_symmetric_, check_hessian_;
  double step_, abstol_, reltol_;
  bool print_progress_, new_jacobian_, new_hessian_, hessian_coloring_;
  std::string validate_ad_file_;

  // FD method as an enum
  FdMode fd_;

  // Types of parallelization
  Parallelization parallelization_;

  // Stats from initialization
  Dict init_stats_;

  /** \brief Constructor

      \identifier{yh} */
  FmuFunction(const std::string& name, Fmu* fmu,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out);

  /** \brief Destructor

      \identifier{yi} */
  ~FmuFunction() override;

  /** \brief Get type name

      \identifier{yj} */
  std::string class_name() const override { return "FmuFunction";}

  ///@{
  /** \brief Options

      \identifier{yk} */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /// Initialize
  void init(const Dict& opts) override;

  // Identify input and output schemes from FmuFunction inputs and outputs
  static void identify_io(
    std::vector<std::string>* scheme_in,
    std::vector<std::string>* scheme_out,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out);

  // Get sparsity pattern for extended Jacobian, Hessian
  Sparsity jac_sp_, hess_sp_;

  // Graph coloring
  Sparsity jac_colors_, hess_colors_;

  // Nonlinearly entering variables
  std::vector<casadi_int> nonlin_;

  // Jacobian memory
  casadi_jac_prob<double> p_;

  // Number of parallel tasks
  casadi_int max_jac_tasks_, max_hess_tasks_, max_n_tasks_;

  ///@{
  /** \brief Number of function inputs and outputs

      \identifier{yl} */
  size_t get_n_in() override { return in_.size();}
  size_t get_n_out() override {return out_.size();}
  ///@}

  /// @{
  /** \brief Retreive sparsities

      \identifier{ym} */
  Sparsity get_sparsity_in(casadi_int i) override;
  Sparsity get_sparsity_out(casadi_int i) override;
  /// @}

  /// @{
  /** \brief Retreive nominal values

      \identifier{yn} */
  std::vector<double> get_nominal_in(casadi_int i) const override;
  std::vector<double> get_nominal_out(casadi_int i) const override;
  /// @}

  // Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

  // Evaluate all tasks numerically, serially or in parallel
  int eval_all(FmuMemory* m, casadi_int n_task,
    bool need_nondiff, bool need_jac, bool need_adj, bool need_hess) const;

  // Evaluate numerically, single thread
  int eval_task(FmuMemory* m, casadi_int task, casadi_int n_task,
    bool need_nondiff, bool need_jac, bool need_adj, bool need_hess) const;

  // Remove NaNs from Hessian (necessary for star coloring approach)
  void remove_nans(double *hess_nz, casadi_int* iw) const;

  // Check extended Hessian
  void check_hessian(FmuMemory* m, const double *hess_nz, casadi_int* iw) const;

  // Make extended Hessian symmetric
  void make_symmetric(double *hess_nz, casadi_int* iw) const;

  ///@{
  /** \brief Return sparsity of Jacobian of an output respect to an input

      \identifier{yo} */
  bool has_jac_sparsity(casadi_int oind, casadi_int iind) const override;
  Sparsity get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const override;
  ///@}

  // Are all inputs/outputs regular?
  bool all_regular() const;

  // Are all inputs/outputs vectors (i.e. not Jacobian or Hessian blocks)?
  bool all_vectors() const;

  ///@{
  /** \brief Full Jacobian

      \identifier{yp} */
  bool has_jacobian() const override;
  Function get_jacobian(const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  ///@{
  /** \brief Reverse mode AD

      \identifier{yq} */
  bool has_reverse(casadi_int nadj) const override;
  Function get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  /** \brief Create memory block

      \identifier{yr} */
  void* alloc_mem() const override;

  /** \brief Initalize memory block

      \identifier{ys} */
  int init_mem(void* mem) const override;

  /** \brief Free memory block

      \identifier{yt} */
  void free_mem(void *mem) const override;

  /// Get all statistics
  Dict get_stats(void* mem) const override;
};

} // namespace casadi
/// \endcond

#endif  // WITH_FMU

#endif // CASADI_FMU_FUNCTION_IMPL_HPP
