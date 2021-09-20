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

// Forward declarations
class DaeBuilderInternal;
class FmuFunction;

// Input structure
class CASADI_EXPORT FmuInput {
 public:
  // Destructor
  virtual ~FmuInput() = 0;
  // Number of elements
  virtual size_t size() const;
  // Access an index
  virtual size_t ind(size_t k) const;
  // Access all indices
  virtual const std::vector<size_t>& ind() const;
  // Get sparsity pattern
  virtual Sparsity sparsity() const = 0;
  // Class name
  virtual std::string class_name() const = 0;
  // It it a regular input?
  virtual bool is_reg() const {return false;}
};

// Regular input
class CASADI_EXPORT RegInput : public FmuInput {
 private:
  // Input indices
  std::vector<size_t> ind_;
 public:
  // Constructor
  RegInput(const std::vector<size_t>& ind) : ind_(ind) {}
  // Destructor
  ~RegInput() override;
  // Number of elements
  size_t size() const override { return ind_.size(); }
  // Access an index
  size_t ind(size_t k) const override { return ind_.at(k);}
  // Access all indices
  const std::vector<size_t>& ind() const override { return ind_;}
  // Get sparsity pattern
  Sparsity sparsity() const override { return Sparsity::dense(size(), 1);}
  // Class name
  std::string class_name() const override { return "RegInput";}
  // It it a regular input?
  bool is_reg() const override {return true;}
};

// Regular input
class CASADI_EXPORT DummyInput : public FmuInput {
 private:
  // Input indices
  size_t dim_;
 public:
  // Constructor
  DummyInput(size_t dim) : dim_(dim) {}
  // Destructor
  ~DummyInput() override;
  // Get sparsity pattern
  Sparsity sparsity() const override { return Sparsity(dim_, 1);}
  // Class name
  std::string class_name() const override { return "DummyInput";}
};

// Output structure
class CASADI_EXPORT FmuOutput {
 public:
  // Destructor
  virtual ~FmuOutput() = 0;
  // Number of elements
  virtual size_t size() const;
  // Access an index
  virtual size_t ind(size_t k) const;
  // Access all indices
  virtual const std::vector<size_t>& ind() const;
  // Access all indices
  virtual const std::vector<size_t>& ind1() const;
  // Access all indices
  virtual const std::vector<size_t>& ind2() const;
  // Get sparsity pattern
  virtual Sparsity sparsity(const FmuFunction& f) const = 0;
  // Class name
  virtual std::string class_name() const = 0;
  // It it a regular output?
  virtual bool is_reg() const {return false;}
  // It it a Jacobian output?
  virtual bool is_jac() const {return false;}
};

// Output structure
class CASADI_EXPORT RegOutput : public FmuOutput {
 public:
  // Output indices
  std::vector<size_t> ind_;
 public:
  // Constructor
  RegOutput(const std::vector<size_t>& ind) : ind_(ind) {}
  // Destructor
  ~RegOutput() override;
  // Number of elements
  size_t size() const override { return ind_.size(); }
  // Access an index
  size_t ind(size_t k) const override { return ind_.at(k);}
  // Access all indices
  const std::vector<size_t>& ind() const override { return ind_;}
  // Get sparsity pattern
  Sparsity sparsity(const FmuFunction& f) const override { return Sparsity::dense(size(), 1);}
  // Class name
  std::string class_name() const override { return "RegOutput";}
  // It it a regular output?
  bool is_reg() const override {return true;}
};

// Output structure
class CASADI_EXPORT JacOutput : public FmuOutput {
 private:
  // Output indices
  std::vector<size_t> ind1_, ind2_;
 public:
  // Constructor
  JacOutput(const std::vector<size_t>& ind1, const std::vector<size_t>& ind2)
    : ind1_(ind1), ind2_(ind2) {}
  // Destructor
  ~JacOutput() override;
  // Get sparsity pattern
  Sparsity sparsity(const FmuFunction& f) const override;
  // Class name
  std::string class_name() const override { return "JacOutput";}
  // It it a Jacobian block?
  bool is_jac() const override {return true;}
  // Access all indices
  const std::vector<size_t>& ind1() const override { return ind1_;}
  // Access all indices
  const std::vector<size_t>& ind2() const override { return ind2_;}
};

// Memory object
struct CASADI_EXPORT FmuMemory : public FunctionMemory {
  // Function object
  const FmuFunction& self;
  // Component memory
  fmi2Component c;
  // Value buffer
  std::vector<double> buffer_;
  // Sensitivities
  std::vector<double> sens_;
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
  // Nominal values
  std::vector<fmi2Real> nominal_out_;
  // Constructor
  explicit FmuMemory(const FmuFunction& self) : self(self), c(nullptr) {}
};

// Interface to binary FMU (shared between derivatives)
struct CASADI_EXPORT Fmu {
  // Constructor
  Fmu(const DaeBuilder& dae);

  // Initialize
  void init();

  // Reference counter
  int counter_;

  // DaeBuilder instance (non-owning reference to avoid circular dependency)
  WeakRef dae_;

  // DLL
  Importer li_;

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

  // Get DaeBuilder instance
  DaeBuilderInternal* dae() const;

  // Load an FMI function
  signal_t load_function(const std::string& symname);

  // Reset solver
  int reset(FmuMemory* m);

  // Enter initialization mode
  int enter_initialization_mode(FmuMemory* m) const;

  // Exit initialization mode
  int exit_initialization_mode(FmuMemory* m) const;

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

  // IO scheme, linear combinations
  std::map<std::string, std::vector<size_t>> scheme_, lc_;

  // Memory
  void* m_;

  // Information about function inputs
  std::vector<FmuInput*> in_;

  // Information about function outputs
  std::vector<FmuOutput*> out_;

  // All Jacobian inputs and outputs
  std::vector<size_t> jac_in_, jac_out_;

  // User-set options
  bool enable_ad_, validate_ad_;
  double step_, abstol_, reltol_, u_aim_, h_min_, h_max_;
  casadi_int h_iter_;

  /// Variable type
  enum FdMode {FORWARD, BACKWARD, CENTRAL, SMOOTHING, N_FDMODE};

  // FD method as an enum
  FdMode fd_;

  // Variables used for initialization, by type
  std::vector<fmi2ValueReference> vr_real_, vr_integer_, vr_boolean_, vr_string_;
  std::vector<fmi2Real> v_real_;
  std::vector<fmi2Integer> v_integer_;
  std::vector<fmi2Boolean> v_boolean_;
  std::vector<std::string> v_string_;

  // Number of perturbations
  casadi_int n_pert() const;

  /** \brief Constructor */
  FmuFunction(const std::string& name, Fmu* fmu,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const std::map<std::string, std::vector<size_t>>& scheme,
      const std::map<std::string, std::vector<size_t>>& lc);

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

  // Get sparsity pattern for extended Jacobian
  Sparsity sp_ext_;

  // Graph coloring
  Sparsity coloring_;

  ///@{
  /** \brief Number of function inputs and outputs */
  size_t get_n_in() override { return in_.size();}
  size_t get_n_out() override {return out_.size();}
  ///@}

  /// @{
  /** \brief Retreive sparsities */
  Sparsity get_sparsity_in(casadi_int i) override {return in_.at(i)->sparsity();}
  Sparsity get_sparsity_out(casadi_int i) override {return out_.at(i)->sparsity(*this);}
  /// @}

  // Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

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

  // Helper function
  static bool has_prefix(const std::string& s);

  // Split prefix
  static std::string pop_prefix(const std::string& s, std::string* rem = 0);

  /** \brief Create memory block */
  void* alloc_mem2() const { return new FmuMemory(*this); }

  /** \brief Initalize memory block */
  int init_mem2(void* mem) const;

  /** \brief Free memory block */
  void free_mem2(void *mem) const;

  // Setup experiment
  void setup_experiment(FmuMemory* m) const;

  // New memory object
  fmi2Component instantiate() const;

  // Set value
  void set(FmuMemory* m, size_t id, double value) const;

  // Request the calculation of a variable
  void request(FmuMemory* m, size_t id, size_t wrt_id = -1) const;

  // Calculate all requested variables
  int eval(FmuMemory* m) const;

  // Get a calculated variable
  void get(FmuMemory* m, size_t id, double* value) const;

  // Set seed
  void set_seed(FmuMemory* m, size_t id, double value) const;

  // Gather user inputs and outputs
  void gather_io(FmuMemory* m) const;

  // Gather user sensitivities
  void gather_sens(FmuMemory* m) const;

  // Calculate directional derivatives using AD
  int eval_ad(FmuMemory* m) const;

  // Calculate directional derivatives using FD
  int eval_fd(FmuMemory* m) const;

  // Calculate directional derivatives
  int eval_derivative(FmuMemory* m) const;

  // Get a derivative
  void get_sens(FmuMemory* m, size_t id, double* value) const;
};

/// Number of entries in enums
template<> struct enum_traits<FmuFunction::FdMode> {
  static const size_t n_enum = FmuFunction::N_FDMODE;
};

/// Convert to string
CASADI_EXPORT std::string to_string(FmuFunction::FdMode v);

} // namespace casadi
/// \endcond

#endif  // WITH_FMU

#endif // CASADI_FMU_FUNCTION_IMPL_HPP
