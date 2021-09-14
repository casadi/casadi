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

struct CASADI_EXPORT Fmu {
  /// Variable type
  enum FdMode {FORWARD, BACKWARD, CENTRAL, SMOOTHING, N_FDMODE};

  // Constructor
  Fmu(const DaeBuilderInternal& self);

  // Destructor
  ~Fmu();

  // Initialize
  void init();

  // Load an FMI function
  signal_t get_function(const std::string& symname);

  // Process message
  static void logger(fmi2ComponentEnvironment componentEnvironment,
    fmi2String instanceName,
    fmi2Status status,
    fmi2String category,
    fmi2String message, ...);

  // New memory object
  int checkout();

  // Free memory object
  void release(int mem);

  // New memory object
  fmi2Component instantiate();

  // Setup experiment
  int setup_experiment(int mem, const FmuFunction& f);

  // Reset solver
  int reset(int mem);

  // Enter initialization mode
  int enter_initialization_mode(int mem);

  // Exit initialization mode
  int exit_initialization_mode(int mem);

  // Set value
  void set(int mem, size_t id, double value);

  // Request the calculation of a variable
  void request(int mem, size_t id, size_t wrt_id = -1);

  // Calculate all requested variables
  int eval(int mem, const FmuFunction& f);

  // Get a calculated variable
  void get(int mem, size_t id, double* value);

  // Set seed
  void set_seed(int mem, size_t id, double value);

  // Gather user inputs and outputs
  void gather_io(int mem);

  // Gather user sensitivities
  void gather_sens(int mem);

  // Calculate directional derivatives using AD
  int eval_ad(int mem, const FmuFunction& f);

  // Calculate directional derivatives using FD
  int eval_fd(int mem, const FmuFunction& f);

  // Calculate directional derivatives
  int eval_derivative(int mem, const FmuFunction& f);

  // Get a derivative
  void get_sens(int mem, size_t id, double* value);

  // Evaluate, non-differentated
  int eval(int mem, const double** arg, double** res, const FmuFunction& f);

  // Evaluate Jacobian numerically, optionally multiply by vector from left
  int eval_jac(int mem, const double** arg, double** res, const FmuFunction& f, bool adj);

  // Get memory object
  fmi2Component memory(int mem);

  // Get memory object and remove from pool
  fmi2Component pop_memory(int mem);

  // DaeBuilder instance
  const DaeBuilderInternal& self_;

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
  fmi2SetRealTYPE* set_real_;
  fmi2SetBooleanTYPE* set_boolean_;
  fmi2GetRealTYPE* get_real_;
  fmi2GetDirectionalDerivativeTYPE* get_directional_derivative_;

  // Callback functions
  fmi2CallbackFunctions functions_;

  // Path to the FMU resource directory
  std::string resource_loc_;

  // Memory object
  struct Memory {
    // Component memory
    fmi2Component c;
    // Currently in use
    bool in_use;
    // Need to initialize
    bool need_init;
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
    // Value references
    std::vector<fmi2ValueReference> vr_in_, vr_out_;
    // Work vector (reals)
    std::vector<fmi2Real> v_in_, v_out_, d_in_, d_out_, fd_out_;
    // Nominal values
    std::vector<fmi2Real> nominal_out_;
    // Constructor
    explicit Memory() : c(0), in_use(false) {}
  };

  // Memory objects
  std::vector<Memory> mem_;
};

// Input structure
class CASADI_EXPORT FmuInput {
 public:
  // Destructor
  virtual ~FmuInput() = 0;
  // Number of elements
  virtual size_t size() const;
  // Access an index
  virtual size_t ind(size_t k) const;
  // Get sparsity pattern
  virtual Sparsity sparsity() const = 0;
  // Class name
  virtual std::string class_name() const = 0;
  // It it a regular input?
  virtual bool is_reg() const {return false;}
};

// Regular input
class CASADI_EXPORT RegInput : public FmuInput {
  friend class FmuFunction;  // workaround for FmuFunction::init
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
  friend class FmuFunction;  // workaround for FmuFunction::init
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
};

class CASADI_EXPORT FmuFunction : public FunctionInternal {
 public:
  // DaeBuilder instance (non-owning reference to avoid circular dependency)
  WeakRef dae_;

  // Information about function inputs
  std::vector<FmuInput*> in_;

  // Information about function outputs
  std::vector<FmuOutput*> out_;

  // User-set options
  bool enable_ad_, validate_ad_;
  double step_, abstol_, reltol_, fmutol_, u_aim_, h_min_, h_max_;
  casadi_int h_iter_;

  // FD method as an enum
  Fmu::FdMode fd_;

  // Number of perturbations
  casadi_int n_pert() const;

  /** \brief Constructor */
  FmuFunction(const std::string& name, const DaeBuilder& dae,
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

  // Jacobian sparsity patterns
  std::vector<std::vector<Sparsity>> sp_jac_;

  // Graph coloring
  Sparsity coloring_;

  // Index offset for concatenated inputs
  std::vector<casadi_int> offset_;

  // Input corresponding to a concatenated input index
  std::vector<casadi_int> offset_map_;

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
  bool has_jac_sparsity(casadi_int oind, casadi_int iind) const override {return true;}
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
  /** \brief Adjoint */
  bool has_reverse(casadi_int nadj) const override {return nadj == 1;}
  Function get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  // Helper function
  static bool has_prefix(const std::string& s);

  // Split prefix
  static std::string pop_prefix(const std::string& s, std::string* rem = 0);
};

/** Jacobian */
class CASADI_EXPORT FmuFunctionJac : public FunctionInternal {
 public:
  /// Constructor
  FmuFunctionJac(const std::string& name) : FunctionInternal(name) {}

  /// Destructor
  ~FmuFunctionJac() override;

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunctionJac";}

  /// Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;
};

/** Adjoint */
class CASADI_EXPORT FmuFunctionAdj : public FunctionInternal {
 public:
  /// Constructor
  FmuFunctionAdj(const std::string& name) : FunctionInternal(name) {}

  /// Destructor
  ~FmuFunctionAdj() override;

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunctionAdj";}

  /// Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;
};

/// Number of entries in enums
template<> struct enum_traits<Fmu::FdMode> {
  static const size_t n_enum = Fmu::N_FDMODE;
};

/// Convert to string
CASADI_EXPORT std::string to_string(Fmu::FdMode v);

} // namespace casadi
/// \endcond

#endif  // WITH_FMU

#endif // CASADI_FMU_FUNCTION_IMPL_HPP
