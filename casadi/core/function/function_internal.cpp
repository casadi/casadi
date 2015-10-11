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


#include "function_internal.hpp"
#include "../mx/casadi_call.hpp"
#include <typeinfo>
#include "../std_vector_tools.hpp"
#include "mx_function.hpp"
#include "external_function.hpp"

#include "../casadi_options.hpp"
#include "../profiling.hpp"

#include <cctype>
#ifdef WITH_DL
#include <cstdlib>
#include <ctime>
#endif // WITH_DL

using namespace std;

namespace casadi {
  FunctionInternal::FunctionInternal(const std::string& name) : name_(name) {
    addOption("verbose",                  OT_BOOLEAN,             false,
              "Verbose evaluation -- for debugging");
    addOption("ad_weight",                OT_REAL,                GenericType(),
              "Weighting factor for derivative calculation."
              "When there is an option of either using forward or reverse mode "
              "directional derivatives, the condition ad_weight*nf<=(1-ad_weight)*na "
              "is used where nf and na are estimates of the number of forward/reverse "
              "mode directional derivatives needed. By default, ad_weight is calculated "
              "automatically, but this can be overridden by setting this option. "
              "In particular, 0 means forcing forward mode and 1 forcing reverse mode. "
              "Leave unset for (class specific) heuristics.");
    addOption("ad_weight_sp",             OT_REAL,                GenericType(),
              "Weighting factor for sparsity pattern calculation calculation."
              "Overrides default behavior. Set to 0 and 1 to force forward and "
              "reverse mode respectively. Cf. option \"ad_weight\".");
    addOption("jac_penalty",             OT_REAL,                 2,
              "When requested for a number of forward/reverse directions,   "
              "it may be cheaper to compute first the full jacobian and then "
              "multiply with seeds, rather than obtain the requested directions "
              "in a straightforward manner. "
              "Casadi uses a heuristic to decide which is cheaper. "
              "A high value of 'jac_penalty' makes it less likely for the heurstic "
              "to chose the full Jacobian strategy. "
              "The special value -1 indicates never to use the full Jacobian strategy");
    addOption("user_data",                OT_VOIDPTR,             GenericType(),
              "A user-defined field that can be used to identify "
              "the function or pass additional information");
    addOption("monitor",                  OT_STRINGVECTOR,        GenericType(),
              "Monitors to be activated", "inputs|outputs");
    addOption("regularity_check",         OT_BOOLEAN,             true,
              "Throw exceptions when NaN or Inf appears during evaluation");
    addOption("inputs_check",             OT_BOOLEAN,             true,
              "Throw exceptions when the numerical values of the inputs don't make sense");
    addOption("gather_stats",             OT_BOOLEAN,             false,
              "Flag to indicate whether statistics must be gathered");
    addOption("input_scheme", OT_STRINGVECTOR, GenericType(), "Custom input scheme");
    addOption("output_scheme", OT_STRINGVECTOR, GenericType(), "Custom output scheme");
    addOption("jit", OT_BOOLEAN, false, "Use just-in-time compiler to speed up the evaluation");
    addOption("compiler", OT_STRING, "clang", "Just-in-time compiler plugin to be used.");
    addOption("jit_options", OT_DICT, GenericType(), "Options to be passed to the jit compiler.");

    verbose_ = false;
    jit_ = false;
    evalD_ = 0;
    user_data_ = 0;
    monitor_inputs_ = false;
    monitor_outputs_ = false;
    sz_arg_ = 0;
    sz_res_ = 0;
    sz_iw_ = 0;
    sz_w_ = 0;
  }

  FunctionInternal::~FunctionInternal() {
  }

  void FunctionInternal::init() {
    setDefaultOptions();

    verbose_ = getOption("verbose");
    jit_ = getOption("jit");
    compilerplugin_ = getOption("compiler").toString();
    if (hasSetOption("jit_options")) jit_options_ = getOption("jit_options");
    regularity_check_ = getOption("regularity_check");

    // Warn for functions with too many inputs or outputs
    casadi_assert_warning(n_in()<10000, "Function " << name_
                          << " has a large number of inputs (" << n_in() << "). "
                          "Changing the problem formulation is strongly encouraged.");
    casadi_assert_warning(n_out()<10000, "Function " << name_
                          << " has a large number of outputs (" << n_out() << "). "
                          "Changing the problem formulation is strongly encouraged.");

    // Resize the matrix that holds the sparsity of the Jacobian blocks
    jac_sparsity_ = jac_sparsity_compact_ =
        SparseStorage<Sparsity>(Sparsity(n_out(), n_in()));
    jac_ = jac_compact_ = SparseStorage<WeakRef>(Sparsity(n_out(), n_in()));

    if (hasSetOption("user_data")) {
      user_data_ = getOption("user_data").toVoidPointer();
    }

    // Pass monitors
    if (hasSetOption("monitor")) {
      const std::vector<std::string>& monitors = getOption("monitor");
      for (std::vector<std::string>::const_iterator it=monitors.begin();it!=monitors.end();it++) {
        monitors_.insert(*it);
      }
    }

    // Custom input scheme
    if (hasSetOption("input_scheme")) {
      ischeme_ = getOption("input_scheme");
    }

    // If input scheme empty, provide default names
    if (ischeme_.empty()) {
      ischeme_.resize(n_in());
      for (size_t i=0; i!=ischeme_.size(); ++i) {
        ischeme_[i] = "i" + CodeGenerator::to_string(i);
      }
    }

    // Custom output scheme
    if (hasSetOption("output_scheme")) {
      oscheme_ = getOption("output_scheme");
    }

    // If output scheme null, provide default names
    if (oscheme_.empty()) {
      oscheme_.resize(n_out());
      for (size_t i=0; i!=oscheme_.size(); ++i) {
        oscheme_[i] = "o" + CodeGenerator::to_string(i);
      }
    }

    monitor_inputs_ = monitored("inputs");
    monitor_outputs_ = monitored("outputs");
    gather_stats_ = getOption("gather_stats");
    inputs_check_ = getOption("inputs_check");
    alloc_arg(0);
    alloc_res(0);

    // Mark the function as initialized
    is_init_ = true;
  }

  void FunctionInternal::finalize() {
    if (jit_) {
      CodeGenerator gen;
      gen.add(shared_from_this<Function>(), "jit_tmp");
      gen.generate("jit_tmp");
      compiler_ = Compiler("jit_tmp.c", compilerplugin_, jit_options_);
      evalD_ = (evalPtr)compiler_.getFunction("jit_tmp");
      casadi_assert_message(evalD_!=0, "Cannot load JIT'ed function.");
    }
  }

  void FunctionInternal::eval(const double** arg, double** res, int* iw, double* w) {
    if (evalD_) {
      evalD_(arg, res, iw, w);
    } else {
      evalD(arg, res, iw, w);
    }
  }

  void FunctionInternal::printDimensions(ostream &stream) const {
    casadi_assert(isInit());
    stream << " Number of inputs: " << n_in() << endl;
    for (int i=0; i<n_in(); ++i) {
      stream << "  Input " << i  << ", a.k.a. \"" << name_in(i) << "\", "
             << input(i).dimString() << ", " << description_in(i) << endl;
    }
    stream << " Number of outputs: " << n_out() << endl;
    for (int i=0; i<n_out(); ++i) {
      stream << "  Output " << i  << ", a.k.a. \"" << name_out(i) << "\", "
             << output(i).dimString() << ", " << description_out(i) << endl;
    }
  }

  void FunctionInternal::print(ostream &stream) const {
    printDimensions(stream);
  }

  void FunctionInternal::repr(ostream &stream) const {
    stream << name_;
  }

  Function FunctionInternal::gradient(int iind, int oind) {
    // Assert scalar
    casadi_assert_message(output(oind).isscalar(),
                          "Only gradients of scalar functions allowed. Use jacobian instead.");

    // Give it a suitable name
    stringstream ss;
    ss << "gradient_" << name_ << "_" << iind << "_" << oind;

    // Output names
    std::vector<std::string> ionames;
    ionames.reserve(1 + n_out());
    ionames.push_back("grad");
    for (int i=0; i<n_out(); ++i) {
      ionames.push_back(oscheme_.at(i));
    }

    // Generate gradient function
    Dict opts = make_dict("input_scheme", ischeme_,
                          "output_scheme", ionames,
                          "jit", jit_, "compiler", compilerplugin_, "jit_options", jit_options_);
    return getGradient(ss.str(), iind, oind, opts);
  }

  Function FunctionInternal::tangent(int iind, int oind) {
    // Assert scalar
    casadi_assert_message(input(iind).isscalar(),
                          "Only tangent of scalar input functions allowed. Use jacobian instead.");

    // Give it a suitable name
    stringstream ss;
    ss << "tangent_" << name_ << "_" << iind << "_" << oind;

    // Output names
    std::vector<std::string> ionames;
    ionames.reserve(1 + n_out());
    ionames.push_back("tangent");
    for (int i=0; i<n_out(); ++i) {
      ionames.push_back(oscheme_.at(i));
    }

    // Generate gradient function
    Dict opts = make_dict("input_scheme", ischeme_,
                          "output_scheme", ionames,
                          "jit", jit_, "compiler", compilerplugin_, "jit_options", jit_options_);
    return getTangent(ss.str(), iind, oind, opts);
  }

  Function FunctionInternal::hessian(int iind, int oind) {
    log("FunctionInternal::hessian");

    // Assert scalar
    casadi_assert_message(output(oind).isscalar(), "Only hessians of scalar functions allowed.");

    // Generate gradient function
    return getHessian(iind, oind);
  }

  Function FunctionInternal::getGradient(const std::string& name, int iind, int oind,
                                         const Dict& opts) {
    Function f = wrapMXFunction();
    f.init();
    return f.gradient(iind, oind);
  }

  Function FunctionInternal::getTangent(const std::string& name, int iind, int oind,
                                        const Dict& opts) {
    Function f = wrapMXFunction();
    f.init();
    return f.tangent(iind, oind);
  }

  MXFunction FunctionInternal::wrapMXFunction() {
    // Construct options of the wrapping MXFunction
    Dict opts;

    // Propagate naming of inputs and outputs
    opts["input_scheme"] = ischeme_;
    opts["output_scheme"] = oscheme_;
    opts["output_scheme"] = oscheme_;

    // Propagate AD parameters
    opts["ad_weight"] = adWeight();
    opts["ad_weight_sp"] = adWeightSp();

    // Propagate JIT
    opts["jit"] = jit_;
    opts["compiler"] = compilerplugin_;
    opts["jit_options"] = jit_options_;

    // Wrap the function
    vector<MX> arg = symbolicInput();
    vector<MX> res = shared_from_this<Function>()(arg);
    return MXFunction("wrap_" + name_, arg, res, opts);
  }

  Function FunctionInternal::getHessian(int iind, int oind) {
    log("FunctionInternal::getHessian");

    // Create gradient function
    log("FunctionInternal::getHessian generating gradient");
    Function g = gradient(iind, oind);

    // Return the Jacobian of the gradient, exploiting symmetry (the gradient has output index 0)
    log("FunctionInternal::getHessian generating Jacobian of gradient");
    return g.jacobian(iind, 0, false, true);
  }

  void FunctionInternal::log(const string& msg) const {
    if (verbose()) {
      userOut() << "CasADi log message: " << msg << endl;
    }
  }

  void FunctionInternal::log(const string& fcn, const string& msg) const {
    if (verbose()) {
      userOut() << "CasADi log message: In \"" << fcn << "\" --- " << msg << endl;
    }
  }

  bool FunctionInternal::verbose() const {
    return verbose_;
  }

  bool FunctionInternal::monitored(const string& mod) const {
    return monitors_.count(mod)>0;
  }

  const Dict & FunctionInternal::getStats() const {
    return stats_;
  }

  GenericType FunctionInternal::getStat(const string & name) const {
    // Locate the statistic
    Dict::const_iterator it = stats_.find(name);

    // Check if found
    if (it == stats_.end()) {
      casadi_error("Statistic: " << name << " has not been set." << endl <<
                   "Note: statistcs are only set after an evaluate call");
    }

    return GenericType(it->second);
  }

  std::vector<MX> FunctionInternal::symbolicInput() const {
    assertInit();
    vector<MX> ret(n_in());
    for (int i=0; i<ret.size(); ++i) {
      stringstream name;
      name << "x_" << i;
      ret[i] = MX::sym(name.str(), input(i).sparsity());
    }
    return ret;
  }

  std::vector<MX> FunctionInternal::symbolicOutput() const {
    assertInit();
    vector<MX> ret(n_out());
    for (int i=0; i<ret.size(); ++i) {
      stringstream name;
      name << "r_" << i;
      ret[i] = MX::sym(name.str(), output(i).sparsity());
    }
    return ret;
  }

  std::vector<MX> FunctionInternal::symbolicOutput(const std::vector<MX>& arg) {
    return shared_from_this<Function>()(arg);
  }

  std::vector<SX> FunctionInternal::symbolicInputSX() const {
    vector<SX> ret(n_in());
    assertInit();
    for (int i=0; i<ret.size(); ++i) {
      stringstream name;
      name << "x_" << i;
      ret[i] = SX::sym(name.str(), input(i).sparsity());
    }
    return ret;
  }
  /// \cond INTERNAL

  void bvec_toggle(bvec_t* s, int begin, int end, int j) {
    for (int i=begin; i<end; ++i) {
      s[i] ^= (bvec_t(1) << j);
    }
  }

  void bvec_clear(bvec_t* s, int begin, int end) {
    for (int i=begin; i<end; ++i) {
      s[i] = 0;
    }
  }


  void bvec_or(bvec_t* s, bvec_t & r, int begin, int end) {
    r = 0;
    for (int i=begin; i<end; ++i) r |= s[i];
  }
  /// \endcond

  Sparsity FunctionInternal::getJacSparsityPlain(int iind, int oind) {
    // Number of nonzero inputs
    int nz_in = input(iind).nnz();

    // Number of nonzero outputs
    int nz_out = output(oind).nnz();

    // Number of forward sweeps we must make
    int nsweep_fwd = nz_in/bvec_size;
    if (nz_in%bvec_size>0) nsweep_fwd++;

    // Number of adjoint sweeps we must make
    int nsweep_adj = nz_out/bvec_size;
    if (nz_out%bvec_size>0) nsweep_adj++;

    // Get weighting factor
    double w = adWeightSp();

    // Use forward mode?
    bool use_fwd = w*nsweep_fwd <= (1-w)*nsweep_adj;

    // Reset the virtual machine
    spInit(use_fwd);

    // Clear the forward seeds/adjoint sensitivities
    for (int ind=0; ind<n_in(); ++ind) {
      vector<double> &v = ibuf_[ind].data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Clear the adjoint seeds/forward sensitivities
    for (int ind=0; ind<n_out(); ++ind) {
      vector<double> &v = obuf_[ind].data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Get seeds and sensitivities
    bvec_t* input_v = get_bvec_t(ibuf_[iind].data());
    bvec_t* output_v = get_bvec_t(obuf_[oind].data());
    bvec_t* seed_v = use_fwd ? input_v : output_v;
    bvec_t* sens_v = use_fwd ? output_v : input_v;

    // Number of sweeps needed
    int nsweep = use_fwd ? nsweep_fwd : nsweep_adj;

    // The number of zeros in the seed and sensitivity directions
    int nz_seed = use_fwd ? nz_in  : nz_out;
    int nz_sens = use_fwd ? nz_out : nz_in;

    // Print
    if (verbose()) {
      userOut() << "FunctionInternal::getJacSparsity: using "
                << (use_fwd ? "forward" : "adjoint") << " mode: ";
      userOut() << nsweep << " sweeps needed for " << nz_seed << " directions" << endl;
    }

    // Progress
    int progress = -10;

    // Temporary vectors
    std::vector<int> jcol, jrow;

    // Loop over the variables, ndir variables at a time
    for (int s=0; s<nsweep; ++s) {

      // Print progress
      if (verbose()) {
        int progress_new = (s*100)/nsweep;
        // Print when entering a new decade
        if (progress_new / 10 > progress / 10) {
          progress = progress_new;
          userOut() << progress << " %"  << endl;
        }
      }

      // Nonzero offset
      int offset = s*bvec_size;

      // Number of local seed directions
      int ndir_local = std::min(bvec_size, nz_seed-offset);

      for (int i=0; i<ndir_local; ++i) {
        seed_v[offset+i] |= bvec_t(1)<<i;
      }

      // Propagate the dependencies
      spEvaluate(use_fwd);

      // Loop over the nonzeros of the output
      for (int el=0; el<nz_sens; ++el) {

        // Get the sparsity sensitivity
        bvec_t spsens = sens_v[el];

        // Clear the sensitivities for the next sweep
        if (!use_fwd) {
          sens_v[el] = 0;
        }

        // If there is a dependency in any of the directions
        if (0!=spsens) {

          // Loop over seed directions
          for (int i=0; i<ndir_local; ++i) {

            // If dependents on the variable
            if ((bvec_t(1) << i) & spsens) {
              // Add to pattern
              jcol.push_back(el);
              jrow.push_back(i+offset);
            }
          }
        }
      }

      // Remove the seeds
      for (int i=0; i<bvec_size && offset+i<nz_seed; ++i) {
        seed_v[offset+i] = 0;
      }
    }

    // Set inputs and outputs to zero
    for (int ind=0; ind<n_in(); ++ind) input(ind).setAll(0);
    for (int ind=0; ind<n_out(); ++ind) output(ind).setAll(0);

    // Construct sparsity pattern
    Sparsity ret = Sparsity::triplet(nz_out, nz_in, use_fwd ? jcol : jrow, use_fwd ? jrow : jcol);

    casadi_msg("Formed Jacobian sparsity pattern (dimension " << ret.size() << ", "
               << ret.nnz() << " nonzeros, " << (100.0*ret.nnz())/ret.numel() << " % nonzeros).");
    casadi_msg("FunctionInternal::getJacSparsity end ");

    // Return sparsity pattern
    return ret;
  }

  Sparsity FunctionInternal::getJacSparsityHierarchicalSymm(int iind, int oind) {
    casadi_assert(spCanEvaluate(true));

    // Number of nonzero inputs
    int nz = input(iind).nnz();

    // Clear the forward seeds/adjoint sensitivities
    for (int ind=0; ind<n_in(); ++ind) {
      vector<double> &v = ibuf_[ind].data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Clear the adjoint seeds/forward sensitivities
    for (int ind=0; ind<n_out(); ++ind) {
      vector<double> &v = obuf_[ind].data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Sparsity triplet accumulator
    std::vector<int> jcol, jrow;

    // Cols/rows of the coarse blocks
    std::vector<int> coarse(2, 0); coarse[1] = nz;

    // Cols/rows of the fine blocks
    std::vector<int> fine;

    // In each iteration, subdivide each coarse block in this many fine blocks
    int subdivision = bvec_size;

    Sparsity r = Sparsity::dense(1, 1);

    // The size of a block
    int granularity = nz;

    int nsweeps = 0;

    bool hasrun = false;

    while (!hasrun || coarse.size()!=nz+1) {
      casadi_msg("Block size: " << granularity);

      // Clear the sparsity triplet acccumulator
      jcol.clear();
      jrow.clear();

      // Clear the fine block structure
      fine.clear();

      Sparsity D = r.starColoring();

      casadi_msg("Star coloring on " << r.dimString() << ": " << D.size2() << " <-> " << D.size1());

      // Reset the virtual machine
      spInit(true);

      // Get seeds and sensitivities
      bvec_t* input_v = get_bvec_t(ibuf_[iind].data());
      bvec_t* output_v = get_bvec_t(obuf_[oind].data());
      bvec_t* seed_v = input_v;
      bvec_t* sens_v = output_v;

      // Clear the seeds
      for (int i=0; i<nz; ++i) seed_v[i]=0;

      // Subdivide the coarse block
      for (int k=0;k<coarse.size()-1;++k) {
        int diff = coarse[k+1]-coarse[k];
        int new_diff = diff/subdivision;
        if (diff%subdivision>0) new_diff++;
        std::vector<int> temp = range(coarse[k], coarse[k+1], new_diff);
        fine.insert(fine.end(), temp.begin(), temp.end());
      }
      if (fine.back()!=coarse.back()) fine.push_back(coarse.back());

      granularity = fine[1] - fine[0];

      // The index into the bvec bit vector
      int bvec_i = 0;

      // Create lookup tables for the fine blocks
      std::vector<int> fine_lookup = lookupvector(fine, nz+1);

      // Triplet data used as a lookup table
      std::vector<int> lookup_col;
      std::vector<int> lookup_row;
      std::vector<int> lookup_value;

      // Loop over all coarse seed directions from the coloring
      for (int csd=0; csd<D.size2(); ++csd) {
        // The maximum number of fine blocks contained in one coarse block
        int n_fine_blocks_max = fine_lookup[coarse[1]]-fine_lookup[coarse[0]];

        int fci_offset = 0;
        int fci_cap = bvec_size-bvec_i;

        // Flag to indicate if all fine blocks have been handled
        bool f_finished = false;

        // Loop while not finished
        while (!f_finished) {

          // Loop over all coarse rows that are found in the coloring for this coarse seed direction
          for (int k=D.colind(csd); k<D.colind(csd+1); ++k) {
            int cci = D.row(k);

            // The first and last rows of the fine block
            int fci_start = fine_lookup[coarse[cci]];
            int fci_end   = fine_lookup[coarse[cci+1]];

            // Local counter that modifies index into bvec
            int bvec_i_mod = 0;

            int value = -bvec_i + fci_offset + fci_start;

            //casadi_assert(value>=0);

            // Loop over the rows of the fine block
            for (int fci = fci_offset;fci<min(fci_end-fci_start, fci_cap);++fci) {

              // Loop over the coarse block cols that appear in the
              // coloring for the current coarse seed direction
              for (int cri=r.colind(cci);cri<r.colind(cci+1);++cri) {
                lookup_col.push_back(r.row(cri));
                lookup_row.push_back(bvec_i+bvec_i_mod);
                lookup_value.push_back(value);
              }

              // Toggle on seeds
              bvec_toggle(seed_v, fine[fci+fci_start], fine[fci+fci_start+1], bvec_i+bvec_i_mod);
              bvec_i_mod++;
            }
          }

          // Bump bvec_i for next major coarse direction
          bvec_i+= min(n_fine_blocks_max, fci_cap);

          // Check if bvec buffer is full
          if (bvec_i==bvec_size || csd==D.size2()-1) {
            // Calculate sparsity for bvec_size directions at once

            // Statistics
            nsweeps+=1;

            // Construct lookup table
            IMatrix lookup = IMatrix::triplet(lookup_row, lookup_col, lookup_value,
                                              bvec_size, coarse.size());

            std::reverse(lookup_col.begin(), lookup_col.end());
            std::reverse(lookup_row.begin(), lookup_row.end());
            std::reverse(lookup_value.begin(), lookup_value.end());
            IMatrix duplicates =
              IMatrix::triplet(lookup_row, lookup_col, lookup_value, bvec_size, coarse.size())
              - lookup;
            duplicates = sparsify(duplicates);
            lookup(duplicates.sparsity()) = -bvec_size;

            // Propagate the dependencies
            spEvaluate(true);

            // Temporary bit work vector
            bvec_t spsens;

            // Loop over the cols of coarse blocks
            for (int cri=0;cri<coarse.size()-1;++cri) {

              // Loop over the cols of fine blocks within the current coarse block
              for (int fri=fine_lookup[coarse[cri]];fri<fine_lookup[coarse[cri+1]];++fri) {
                // Lump individual sensitivities together into fine block
                bvec_or(sens_v, spsens, fine[fri], fine[fri+1]);

                // Loop over all bvec_bits
                for (int bvec_i=0;bvec_i<bvec_size;++bvec_i) {
                  if (spsens & (bvec_t(1) << bvec_i)) {
                    // if dependency is found, add it to the new sparsity pattern
                    int lk = lookup.elem(bvec_i, cri);
                    if (lk>-bvec_size) {
                      jrow.push_back(bvec_i+lk);
                      jcol.push_back(fri);
                      jrow.push_back(fri);
                      jcol.push_back(bvec_i+lk);
                    }
                  }
                }
              }
            }

            // Clear the forward seeds/adjoint sensitivities, ready for next bvec sweep
            for (int ind=0; ind<n_in(); ++ind) {
              vector<double> &v = ibuf_[ind].data();
              if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
            }

            // Clear the adjoint seeds/forward sensitivities, ready for next bvec sweep
            for (int ind=0; ind<n_out(); ++ind) {
              vector<double> &v = obuf_[ind].data();
              if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
            }

            // Clean lookup table
            lookup_col.clear();
            lookup_row.clear();
            lookup_value.clear();
          }

          if (n_fine_blocks_max>fci_cap) {
            fci_offset += min(n_fine_blocks_max, fci_cap);
            bvec_i = 0;
            fci_cap = bvec_size;
          } else {
            f_finished = true;
          }

        }

      }

      // Construct fine sparsity pattern
      r = Sparsity::triplet(fine.size()-1, fine.size()-1, jrow, jcol);

      // There may be false positives here that are not present
      // in the reverse mode that precedes it.
      // This can lead to an assymetrical result
      //  cf. #1522
      r=r*r.T();

      coarse = fine;
      hasrun = true;
    }

    casadi_msg("Number of sweeps: " << nsweeps);
    casadi_msg("Formed Jacobian sparsity pattern (dimension " << r.size() <<
               ", " << r.nnz() << " nonzeros, " << (100.0*r.nnz())/r.numel() << " % nonzeros).");

    return r.T();
  }

  Sparsity FunctionInternal::getJacSparsityHierarchical(int iind, int oind) {
    // Number of nonzero inputs
    int nz_in = input(iind).nnz();

    // Number of nonzero outputs
    int nz_out = output(oind).nnz();

    // Clear the forward seeds/adjoint sensitivities
    for (int ind=0; ind<n_in(); ++ind) {
      vector<double> &v = ibuf_[ind].data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Clear the adjoint seeds/forward sensitivities
    for (int ind=0; ind<n_out(); ++ind) {
      vector<double> &v = obuf_[ind].data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Sparsity triplet accumulator
    std::vector<int> jcol, jrow;

    // Cols of the coarse blocks
    std::vector<int> coarse_col(2, 0); coarse_col[1] = nz_out;
    // Rows of the coarse blocks
    std::vector<int> coarse_row(2, 0); coarse_row[1] = nz_in;

    // Cols of the fine blocks
    std::vector<int> fine_col;

    // Rows of the fine blocks
    std::vector<int> fine_row;

    // In each iteration, subdivide each coarse block in this many fine blocks
    int subdivision = bvec_size;

    Sparsity r = Sparsity::dense(1, 1);

    // The size of a block
    int granularity_row = nz_in;
    int granularity_col = nz_out;

    bool use_fwd = true;

    int nsweeps = 0;

    bool hasrun = false;

    // Get weighting factor
    double w = adWeightSp();

    // Lookup table for bvec_t
    std::vector<bvec_t> bvec_lookup;
    bvec_lookup.reserve(bvec_size);
    for (int i=0;i<bvec_size;++i) {
      bvec_lookup.push_back(bvec_t(1) << i);
    }

    while (!hasrun || coarse_col.size()!=nz_out+1 || coarse_row.size()!=nz_in+1) {
      casadi_msg("Block size: " << granularity_col << " x " << granularity_row);

      // Clear the sparsity triplet acccumulator
      jcol.clear();
      jrow.clear();

      // Clear the fine block structure
      fine_row.clear();
      fine_col.clear();

      // r transpose will be needed in the algorithm
      Sparsity rT = r.T();

      /**       Decide which ad_mode to take           */

      // Forward mode
      Sparsity D1 = rT.unidirectionalColoring(r);
      // Adjoint mode
      Sparsity D2 = r.unidirectionalColoring(rT);

      casadi_msg("Coloring on " << r.dimString() << " (fwd seeps: " << D1.size2() <<
                 " , adj sweeps: " << D2.size1() << ")");

      // Use forward mode?
      int fwd_cost = use_fwd ? granularity_row: granularity_col;
      int adj_cost = use_fwd ? granularity_col: granularity_row;

      // Use whatever required less colors if we tried both (with preference to forward mode)
      if ((w*D1.size2()*fwd_cost <= (1-w)*D2.size2()*adj_cost)) {
        use_fwd = true;
        casadi_msg("Forward mode chosen (fwd cost: " << w*D1.size2()*fwd_cost << ", adj cost: "
                   << (1-w)*D2.size2()*adj_cost << ")");
      } else {
        use_fwd = false;
        casadi_msg("Adjoint mode chosen (adj cost: " << w*D1.size2()*fwd_cost << ", adj cost: "
                   << (1-w)*D2.size2()*adj_cost << ")");
      }

      // Reset the virtual machine
      spInit(use_fwd);

      // Get seeds and sensitivities
      bvec_t* input_v = get_bvec_t(ibuf_[iind].data());
      bvec_t* output_v = get_bvec_t(obuf_[oind].data());
      bvec_t* seed_v = use_fwd ? input_v : output_v;
      bvec_t* sens_v = use_fwd ? output_v : input_v;

      // The number of zeros in the seed and sensitivity directions
      int nz_seed = use_fwd ? nz_in  : nz_out;
      int nz_sens = use_fwd ? nz_out : nz_in;

      // Clear the seeds
      for (int i=0; i<nz_seed; ++i) seed_v[i]=0;

      // Choose the active jacobian coloring scheme
      Sparsity D = use_fwd ? D1 : D2;

      // Adjoint mode amounts to swapping
      if (!use_fwd) {
        std::swap(coarse_col, coarse_row);
        std::swap(granularity_col, granularity_row);
        std::swap(r, rT);
      }

      // Subdivide the coarse block cols
      for (int k=0;k<coarse_col.size()-1;++k) {
        int diff = coarse_col[k+1]-coarse_col[k];
        int new_diff = diff/subdivision;
        if (diff%subdivision>0) new_diff++;
        std::vector<int> temp = range(coarse_col[k], coarse_col[k+1], new_diff);
        fine_col.insert(fine_col.end(), temp.begin(), temp.end());
      }
      // Subdivide the coarse block rows
      for (int k=0;k<coarse_row.size()-1;++k) {
        int diff = coarse_row[k+1]-coarse_row[k];
        int new_diff = diff/subdivision;
        if (diff%subdivision>0) new_diff++;
        std::vector<int> temp = range(coarse_row[k], coarse_row[k+1], new_diff);
        fine_row.insert(fine_row.end(), temp.begin(), temp.end());
      }
      if (fine_row.back()!=coarse_row.back()) fine_row.push_back(coarse_row.back());
      if (fine_col.back()!=coarse_col.back()) fine_col.push_back(coarse_col.back());

      granularity_col = fine_col[1] - fine_col[0];
      granularity_row = fine_row[1] - fine_row[0];

      // The index into the bvec bit vector
      int bvec_i = 0;

      // Create lookup tables for the fine blocks
      std::vector<int> fine_col_lookup = lookupvector(fine_col, nz_sens+1);
      std::vector<int> fine_row_lookup = lookupvector(fine_row, nz_seed+1);

      // Triplet data used as a lookup table
      std::vector<int> lookup_col;
      std::vector<int> lookup_row;
      std::vector<int> lookup_value;

      // Loop over all coarse seed directions from the coloring
      for (int csd=0; csd<D.size2(); ++csd) {

        // The maximum number of fine blocks contained in one coarse block
        int n_fine_blocks_max = fine_row_lookup[coarse_row[1]]-fine_row_lookup[coarse_row[0]];

        int fci_offset = 0;
        int fci_cap = bvec_size-bvec_i;

        // Flag to indicate if all fine blocks have been handled
        bool f_finished = false;

        // Loop while not finished
        while (!f_finished) {

          // Loop over all coarse rows that are found in the coloring for this coarse seed direction
          for (int k=D.colind(csd); k<D.colind(csd+1); ++k) {
            int cci = D.row(k);

            // The first and last rows of the fine block
            int fci_start = fine_row_lookup[coarse_row[cci]];
            int fci_end   = fine_row_lookup[coarse_row[cci+1]];

            // Local counter that modifies index into bvec
            int bvec_i_mod = 0;

            int value = -bvec_i + fci_offset + fci_start;

            // Loop over the rows of the fine block
            for (int fci = fci_offset;fci<min(fci_end-fci_start, fci_cap);++fci) {

              // Loop over the coarse block cols that appear in the coloring
              // for the current coarse seed direction
              for (int cri=rT.colind(cci);cri<rT.colind(cci+1);++cri) {
                lookup_col.push_back(rT.row(cri));
                lookup_row.push_back(bvec_i+bvec_i_mod);
                lookup_value.push_back(value);
              }

              // Toggle on seeds
              bvec_toggle(seed_v, fine_row[fci+fci_start], fine_row[fci+fci_start+1],
                          bvec_i+bvec_i_mod);
              bvec_i_mod++;
            }
          }

          // Bump bvec_i for next major coarse direction
          bvec_i+= min(n_fine_blocks_max, fci_cap);

          // Check if bvec buffer is full
          if (bvec_i==bvec_size || csd==D.size2()-1) {
            // Calculate sparsity for bvec_size directions at once

            // Statistics
            nsweeps+=1;

            // Construct lookup table
            IMatrix lookup = IMatrix::triplet(lookup_row, lookup_col, lookup_value, bvec_size,
                                              coarse_col.size());

            // Propagate the dependencies
            spEvaluate(use_fwd);

            // Temporary bit work vector
            bvec_t spsens;

            // Loop over the cols of coarse blocks
            for (int cri=0;cri<coarse_col.size()-1;++cri) {

              // Loop over the cols of fine blocks within the current coarse block
              for (int fri=fine_col_lookup[coarse_col[cri]];
                   fri<fine_col_lookup[coarse_col[cri+1]];++fri) {
                // Lump individual sensitivities together into fine block
                bvec_or(sens_v, spsens, fine_col[fri], fine_col[fri+1]);

                // Next iteration if no sparsity
                if (!spsens) continue;

                // Loop over all bvec_bits
                for (int bvec_i=0;bvec_i<bvec_size;++bvec_i) {
                  if (spsens & bvec_lookup[bvec_i]) {
                    // if dependency is found, add it to the new sparsity pattern
                    jrow.push_back(bvec_i+lookup.elem(bvec_i, cri));
                    jcol.push_back(fri);
                  }
                }
              }
            }

            // Clear the forward seeds/adjoint sensitivities, ready for next bvec sweep
            for (int ind=0; ind<n_in(); ++ind) {
              vector<double> &v = ibuf_[ind].data();
              if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
            }

            // Clear the adjoint seeds/forward sensitivities, ready for next bvec sweep
            for (int ind=0; ind<n_out(); ++ind) {
              vector<double> &v = obuf_[ind].data();
              if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
            }

            // Clean lookup table
            lookup_col.clear();
            lookup_row.clear();
            lookup_value.clear();
          }

          if (n_fine_blocks_max>fci_cap) {
            fci_offset += min(n_fine_blocks_max, fci_cap);
            bvec_i = 0;
            fci_cap = bvec_size;
          } else {
            f_finished = true;
          }

        }

      }

      // Swap results if adjoint mode was used
      if (use_fwd) {
        // Construct fine sparsity pattern
        r = Sparsity::triplet(fine_row.size()-1, fine_col.size()-1, jrow, jcol);
        coarse_col = fine_col;
        coarse_row = fine_row;
      } else {
        // Construct fine sparsity pattern
        r = Sparsity::triplet(fine_col.size()-1, fine_row.size()-1, jcol, jrow);
        coarse_col = fine_row;
        coarse_row = fine_col;
      }
      hasrun = true;
    }
    casadi_msg("Number of sweeps: " << nsweeps);
    casadi_msg("Formed Jacobian sparsity pattern (dimension " << r.size() <<
               ", " << r.nnz() << " nonzeros, " << (100.0*r.nnz())/r.numel() << " % nonzeros).");

    return r.T();
  }

  Sparsity FunctionInternal::getJacSparsity(int iind, int oind, bool symmetric) {
    // Check if we are able to propagate dependencies through the function
    if (spCanEvaluate(true) || spCanEvaluate(false)) {
      Sparsity sp;
      if (input(iind).nnz()>3*bvec_size && output(oind).nnz()>3*bvec_size) {
        if (symmetric) {
          sp = getJacSparsityHierarchicalSymm(iind, oind);
        } else {
          sp = getJacSparsityHierarchical(iind, oind);
        }
      } else {
        sp = getJacSparsityPlain(iind, oind);
      }
      // There may be false positives here that are not present
      // in the reverse mode that precedes it.
      // This can lead to an assymetrical result
      //  cf. #1522
      if (symmetric) sp=sp*sp.T();
      return sp;
    } else {
      // Dense sparsity by default
      return Sparsity::dense(output(oind).nnz(), input(iind).nnz());
    }
  }

  void FunctionInternal::setJacSparsity(const Sparsity& sp, int iind, int oind, bool compact) {
    if (compact) {
      jac_sparsity_compact_.elem(oind, iind) = sp;
    } else {
      jac_sparsity_.elem(oind, iind) = sp;
    }
  }

  Sparsity& FunctionInternal::jacSparsity(int iind, int oind, bool compact, bool symmetric) {
    casadi_assert_message(isInit(), "Function not initialized.");

    // Get an owning reference to the block
    Sparsity jsp = compact ? jac_sparsity_compact_.elem(oind, iind)
        : jac_sparsity_.elem(oind, iind);

    // Generate, if null
    if (jsp.isNull()) {
      if (compact) {

        // Use internal routine to determine sparsity
        jsp = getJacSparsity(iind, oind, symmetric);

      } else {

        // Get the compact sparsity pattern
        Sparsity sp = jacSparsity(iind, oind, true, symmetric);

        // Enlarge if sparse output
        if (output(oind).numel()!=sp.size1()) {
          casadi_assert(sp.size1()==output(oind).nnz());

          // New row for each old row
          vector<int> row_map = output(oind).sparsity().find();

          // Insert rows
          sp.enlargeRows(output(oind).numel(), row_map);
        }

        // Enlarge if sparse input
        if (input(iind).numel()!=sp.size2()) {
          casadi_assert(sp.size2()==input(iind).nnz());

          // New column for each old column
          vector<int> col_map = input(iind).sparsity().find();

          // Insert columns
          sp.enlargeColumns(input(iind).numel(), col_map);
        }

        // Save
        jsp = sp;
      }
    }

    // If still null, not dependent
    if (jsp.isNull()) {
      jsp = Sparsity(output(oind).nnz(), input(iind).nnz());
    }

    // Return a reference to the block
    Sparsity& jsp_ref = compact ? jac_sparsity_compact_.elem(oind, iind) :
        jac_sparsity_.elem(oind, iind);
    jsp_ref = jsp;
    return jsp_ref;
  }

  void FunctionInternal::getPartition(int iind, int oind, Sparsity& D1, Sparsity& D2,
                                      bool compact, bool symmetric) {
    log("FunctionInternal::getPartition begin");

    // Sparsity pattern with transpose
    Sparsity &AT = jacSparsity(iind, oind, compact, symmetric);
    Sparsity A = symmetric ? AT : AT.T();

    // Get seed matrices by graph coloring
    if (symmetric) {
      casadi_assert(numDerForward()>0);

      // Star coloring if symmetric
      log("FunctionInternal::getPartition starColoring");
      D1 = A.starColoring();
      casadi_msg("Star coloring completed: " << D1.size2() << " directional derivatives needed ("
                 << A.size1() << " without coloring).");

    } else {
      casadi_assert(numDerForward()>0 || numDerReverse()>0);
      // Get weighting factor
      double w = adWeight();

      // Which AD mode?
      bool test_ad_fwd=w<1, test_ad_adj=w>0;

      // Best coloring encountered so far (relatively tight upper bound)
      double best_coloring = numeric_limits<double>::infinity();

      // Test forward mode first?
      bool test_fwd_first = test_ad_fwd && w*A.size1() <= (1-w)*A.size2();
      int mode_fwd = test_fwd_first ? 0 : 1;

      // Test both coloring modes
      for (int mode=0; mode<2; ++mode) {
        // Is this the forward mode?
        bool fwd = mode==mode_fwd;

        // Skip?
        if (!test_ad_fwd && fwd) continue;
        if (!test_ad_adj && !fwd) continue;

        // Perform the coloring
        if (fwd) {
          log("FunctionInternal::getPartition unidirectional coloring (forward mode)");
          int max_colorings_to_test = best_coloring>=w*A.size1() ? A.size1() :
            floor(best_coloring/w);
          D1 = AT.unidirectionalColoring(A, max_colorings_to_test);
          if (D1.isNull()) {
            if (verbose()) userOut() << "Forward mode coloring interrupted (more than "
                               << max_colorings_to_test << " needed)." << endl;
          } else {
            if (verbose()) userOut() << "Forward mode coloring completed: "
                               << D1.size2() << " directional derivatives needed ("
                               << A.size1() << " without coloring)." << endl;
            D2 = Sparsity();
            best_coloring = w*D1.size2();
          }
        } else {
          log("FunctionInternal::getPartition unidirectional coloring (adjoint mode)");
          int max_colorings_to_test = best_coloring>=(1-w)*A.size2() ? A.size2() :
            floor(best_coloring/(1-w));

          D2 = A.unidirectionalColoring(AT, max_colorings_to_test);
          if (D2.isNull()) {
            if (verbose()) userOut() << "Adjoint mode coloring interrupted (more than "
                               << max_colorings_to_test << " needed)." << endl;
          } else {
            if (verbose()) userOut() << "Adjoint mode coloring completed: "
                               << D2.size2() << " directional derivatives needed ("
                               << A.size2() << " without coloring)." << endl;
            D1 = Sparsity();
            best_coloring = (1-w)*D2.size2();
          }
        }
      }

    }
    log("FunctionInternal::getPartition end");
  }

  void FunctionInternal::evaluate() {
    // Allocate temporary memory if needed
    alloc();

    // Get pointers to input arguments
    int n_in = this->n_in();
    vector<const double*> arg(sz_arg());
    for (int i=0; i<n_in; ++i) arg[i]=input(i).ptr();

    // Get pointers to output arguments
    int n_out = this->n_out();
    vector<double*> res(sz_res());
    for (int i=0; i<n_out; ++i) res[i]=output(i).ptr();

    // Call memory-less
    eval(getPtr(arg), getPtr(res), getPtr(iw_tmp_), getPtr(w_tmp_));
  }

  void FunctionInternal::evalD(const double** arg,
                               double** res, int* iw, double* w) {
    // Number of inputs and outputs
    int num_in = n_in();
    int num_out = n_out();

    // Pass the inputs to the function
    for (int i=0; i<num_in; ++i) {
      if (arg[i] != 0) {
        setInputNZ(arg[i], i);
      } else {
        setInput(0., i);
      }
    }

    // Evaluate
    evaluate();

    // Get the outputs
    for (int i=0; i<num_out; ++i) {
      if (res[i] != 0) getOutputNZ(res[i], i);
    }
  }

  void FunctionInternal::evalSX(const SXElement** arg, SXElement** res,
                                int* iw, SXElement* w) {
    casadi_assert(canEvalSX());

    // Number of inputs and outputs
    int num_in = n_in();
    int num_out = n_out();

    // Create input arguments
    vector<SX> argv(num_in);
    for (int i=0; i<num_in; ++i) {
      argv[i] = SX::zeros(input(i).sparsity());
      if (arg[i] != 0) {
        std::copy(arg[i], arg[i]+argv[i].nnz(), argv[i].begin());
      }
    }

    // Evaluate symbolically
    vector<SX> resv;
    evalSX(argv, resv);

    // Collect the result
    for (int i = 0; i < num_out; ++i) {
      if (res[i] != 0) {
        std::copy(resv[i].begin(), resv[i].end(), res[i]);
      }
    }
  }

  void FunctionInternal::evalSX(const std::vector<SX>& arg, std::vector<SX>& res) {
    casadi_assert(canEvalSX());

    // Get the number of inputs and outputs
    int num_in = n_in();
    int num_out = n_out();

    // Check if matching input sparsity
    bool matching_sparsity = true;
    casadi_assert(arg.size()==num_in);
    for (int i=0; matching_sparsity && i<num_in; ++i)
      matching_sparsity = arg[i].sparsity()==input(i).sparsity();

    // Correct input sparsity if needed
    if (!matching_sparsity) {
      vector<SX> arg2(arg);
      for (int i=0; i<num_in; ++i)
        if (arg2[i].sparsity()!=input(i).sparsity())
          arg2[i] = project(arg2[i], input(i).sparsity());
      return evalSX(arg2, res);
    }

    // Allocate results
    res.resize(num_out);
    for (int i=0; i<num_out; ++i)
      if (res[i].sparsity()!=output(i).sparsity())
        res[i] = SX::zeros(output(i).sparsity());

    // Allocate temporary memory if needed
    iw_tmp_.resize(sz_iw());
    vector<SXElement> w_tmp(sz_w());

    // Get pointers to input arguments
    vector<const SXElement*> argp(sz_arg());
    for (int i=0; i<arg.size(); ++i) argp[i]=getPtr(arg[i]);

    // Get pointers to output arguments
    vector<SXElement*> resp(sz_res());
    for (int i=0; i<n_out(); ++i) resp[i]=getPtr(res[i]);

    // Call memory-less
    evalSX(getPtr(argp), getPtr(resp), getPtr(iw_tmp_), getPtr(w_tmp));
  }

  void FunctionInternal::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    MXFunction f = wrapMXFunction();
    f.init();
    f.call(arg, res, true);
  }

  void FunctionInternal::spEvaluate(bool fwd) {
    // Allocate temporary memory if needed
    iw_tmp_.resize(sz_iw());
    w_tmp_.resize(sz_w());
    int *iw = getPtr(iw_tmp_);
    bvec_t *w = reinterpret_cast<bvec_t*>(getPtr(w_tmp_));

    // Get pointers to output arguments
    vector<bvec_t*> res(sz_res());
    for (int i=0; i<n_out(); ++i) res[i]=reinterpret_cast<bvec_t*>(output(i).ptr());

    if (fwd) {
      // Get pointers to input arguments
      vector<const bvec_t*> arg(sz_arg());
      for (int i=0; i<n_in(); ++i) arg[i]=reinterpret_cast<const bvec_t*>(input(i).ptr());

      // Call memory-less
      spFwd(getPtr(arg), getPtr(res), iw, w);
    } else {
      // Get pointers to input arguments
      vector<bvec_t*> arg(sz_arg());
      for (int i=0; i<n_in(); ++i) arg[i]=reinterpret_cast<bvec_t*>(input(i).ptr());

      // Call memory-less
      spAdj(getPtr(arg), getPtr(res), iw, w);
    }
  }

  void FunctionInternal::spEvaluateViaJacSparsity(bool fwd) {
    if (fwd) {
      // Clear the outputs
      for (int oind = 0; oind < n_out(); ++oind) {
        // Get data array for output and clear it
        bvec_t *outputd = get_bvec_t(output(oind).data());
        fill_n(outputd, output(oind).nnz(), 0);
      }
    }

    // Loop over inputs
    for (int iind = 0; iind < n_in(); ++iind) {
      // Skip if no seeds
      if (fwd && input(iind).isempty())
        continue;

      // Get data array for input
      bvec_t *inputd = get_bvec_t(input(iind).data());

      // Loop over outputs
      for (int oind = 0; oind < n_out(); ++oind) {

        // Skip if no seeds
        if (!fwd && output(oind).isempty())
          continue;


        // Save the seeds/sens (#1532)
        // This code should stay in place until refactored away
        std::vector<double> store_in(nnz_in());
        std::vector<double> store_out(nnz_out());

        int offset = 0;
        for (int i=0;i<n_in();++i) {
          std::copy(input(i).data().begin(), input(i).data().begin()+input(i).nnz(),
            store_in.begin()+offset);
          offset+=input(i).nnz();
        }
        offset = 0;
        for (int i=0;i<n_out();++i) {
          std::copy(output(i).data().begin(), output(i).data().begin()+output(i).nnz(),
            store_out.begin()+offset);
          offset+=output(i).nnz();
        }


        // Get the sparsity of the Jacobian block
        Sparsity sp = jacSparsity(iind, oind, true, false);
        if (sp.isNull() || sp.nnz() == 0)
          continue; // Skip if zero

        // recover the seeds/sens
        offset = 0;
        for (int i=0;i<n_in();++i) {
          std::copy(store_in.begin()+offset, store_in.begin()+offset+input(i).nnz(),
            input(i).data().begin());
          offset+=input(i).nnz();
        }
        offset = 0;
        for (int i=0;i<n_out();++i) {
          std::copy(store_out.begin()+offset, store_out.begin()+offset+output(i).nnz(),
            output(i).data().begin());
          offset+=output(i).nnz();
        }

        const int d1 = sp.size2();
        //const int d2 = sp.size1();
        const int* colind = sp.colind();
        const int* row = sp.row();

        // Get data array for output
        bvec_t *outputd = get_bvec_t(output(oind).data());

        // Carry out the sparse matrix-vector multiplication
        for (int cc = 0; cc < d1; ++cc) {
          for (int el = colind[cc]; el < colind[cc + 1]; ++el) {
            // Get row
            int rr = row[el];

            // Propagate dependencies
            if (fwd) {
              outputd[rr] |= inputd[cc];
            } else {
              inputd[cc] |= outputd[rr];
            }
          }
        }
      }
    }
    if (!fwd) {
      for (int oind=0; oind < n_out(); ++oind) {
        vector<double> &w = output(oind).data();
        fill_n(get_bvec_t(w), w.size(), bvec_t(0));
      }
    }
  }

  Function FunctionInternal::jacobian(int iind, int oind, bool compact, bool symmetric) {

    // Return value
    WeakRef cached = compact ? jac_compact_.elem(oind, iind) : jac_.elem(oind, iind);

    // Check if cached
    if (cached.alive()) {
      // Return an owning reference
      return shared_cast<Function>(cached.shared());

    } else {
      // Give it a suitable name
      stringstream ss;
      ss << "jacobian_" << name_ << "_" << iind << "_" << oind;

      // Output names
      std::vector<std::string> ionames;
      ionames.reserve(1 + n_out());
      ionames.push_back("jac");
      for (int i=0; i<n_out(); ++i) {
        ionames.push_back(oscheme_.at(i));
      }

      // Generate a Jacobian
      Dict opts = make_dict("verbose", verbose_,
                            "input_scheme", ischeme_, "output_scheme", ionames,
                            "jit", jit_, "compiler", compilerplugin_,
                            "jit_options", jit_options_);
      Function ret = getJacobian(ss.str(), iind, oind, compact, symmetric, opts);

      // Save in cache
      compact ? jac_compact_.elem(oind, iind) : jac_.elem(oind, iind) = ret;
      return ret;
    }
  }

  void FunctionInternal::setJacobian(const Function& jac, int iind, int oind, bool compact) {
    if (compact) {
      jac_compact_.elem(oind, iind) = jac;
    } else {
      jac_.elem(oind, iind) = jac;
    }
  }

  Function FunctionInternal::
  getJacobian(const std::string& name, int iind, int oind, bool compact, bool symmetric,
              const Dict& opts) {
    return getNumericJacobian(name, iind, oind, compact, symmetric, opts);
  }

  Function FunctionInternal::derForward(int nfwd) {
    casadi_assert(nfwd>=0);

    // Check if there are enough forward directions allocated
    if (nfwd>=derivative_fwd_.size()) {
      derivative_fwd_.resize(nfwd+1);
    }

    // Quick return if already cached
    if (derivative_fwd_[nfwd].alive()) {
      return shared_cast<Function>(derivative_fwd_[nfwd].shared());
    }

    // Give it a suitable name
    stringstream ss;
    ss << "fwd" << nfwd << "_" << name_;
    string name = ss.str();

    // Get the number of inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Names of inputs
    std::vector<std::string> i_names;
    i_names.reserve(n_in + n_out + n_in*nfwd);

    // Nondifferentiated inputs
    for (int i=0; i<n_in; ++i) {
      i_names.push_back("der_" + ischeme_.at(i));
    }

    // Nondifferentiated outputs (given)
    for (int i=0; i<n_out; ++i) {
      i_names.push_back("der_" + oscheme_.at(i));
    }

    // Forward seeds
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<n_in; ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << ischeme_.at(i);
        i_names.push_back(ss.str());
      }
    }

    // Names of outputs
    std::vector<std::string> o_names;
    o_names.reserve(n_out*nfwd);

    // Forward sensitivities
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<n_out; ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << oscheme_.at(i);
        o_names.push_back(ss.str());
      }
    }

    // Options
    Dict opts = make_dict("input_scheme", i_names, "output_scheme", o_names,
                          "jit", jit_, "compiler", compilerplugin_, "jit_options", jit_options_);

    // Return value
    casadi_assert(numDerForward()>0);
    Function ret = getDerForward(name, nfwd, opts);

    // Consistency check for inputs
    for (int i=0; i<ret.n_in(); ++i) {
      const Sparsity& sp = i<n_in ? input(i).sparsity() :
        i<n_in+n_out ? output(i-n_in).sparsity() :
        input((i-n_in-n_out) % n_in).sparsity();
      casadi_assert_message(ret.input(i).size()==sp.size(),
                            "Incorrect shape for " << ret << " input " << i << " \""
                            << i_names.at(i) << "\". Expected " << sp.size()
                            << " but got " << ret.input(i).size());
    }

    // Consistency check for outputs
    for (int i=0; i<ret.n_out(); ++i) {
      const Sparsity& sp = output(i % n_out).sparsity();
      casadi_assert_message(ret.output(i).size()==sp.size(),
                            "Incorrect shape for " << ret << " output " << i << " \""
                            << o_names.at(i) << "\". Expected " << sp.size()
                            << " but got " << ret.output(i).size());
    }

    // Save to cache
    derivative_fwd_[nfwd] = ret;

    // Return generated function
    return ret;
  }

  Function FunctionInternal::derReverse(int nadj) {
    casadi_assert(nadj>=0);

    // Check if there are enough adjoint directions allocated
    if (nadj>=derivative_adj_.size()) {
      derivative_adj_.resize(nadj+1);
    }

    // Quick return if already cached
    if (derivative_adj_[nadj].alive()) {
      return shared_cast<Function>(derivative_adj_[nadj].shared());
    }

    // Give it a suitable name
    stringstream ss;
    ss << "adj" << nadj << "_" << name_;
    string name = ss.str();

    // Get the number of inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Names of inputs
    std::vector<std::string> i_names;
    i_names.reserve(n_in + n_out + n_out*nadj);

    // Nondifferentiated inputs
    for (int i=0; i<n_in; ++i) {
      i_names.push_back("der_" + ischeme_.at(i));
    }

    // Nondifferentiated outputs (given)
    for (int i=0; i<n_out; ++i) {
      i_names.push_back("der_" + oscheme_.at(i));
    }

    // Adjoint seeds
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_out; ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << oscheme_.at(i);
        i_names.push_back(ss.str());
      }
    }

    // Names of outputs
    std::vector<std::string> o_names;
    o_names.reserve(n_in*nadj);

    // Adjoint sensitivities
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_in; ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << ischeme_.at(i);
        o_names.push_back(ss.str());
      }
    }

    // Options
    Dict opts = make_dict("input_scheme", i_names,
                          "output_scheme", o_names,
                          "jit", jit_,
                          "compiler", compilerplugin_,
                          "jit_options", jit_options_);

    // Return value
    casadi_assert(numDerReverse()>0);
    Function ret = getDerReverse(name, nadj, opts);

    // Consistency check for inputs
    for (int i=0; i<ret.n_in(); ++i) {
      const Sparsity& sp = i<n_in ? input(i).sparsity() :
        i<n_in+n_out ? output(i-n_in).sparsity() :
        output((i-n_in-n_out) % n_out).sparsity();
      casadi_assert_message(ret.input(i).size()==sp.size(),
                            "Incorrect shape for " << ret << " input " << i << " \""
                            << i_names.at(i) << "\". Expected " << sp.size()
                            << " but got " << ret.input(i).size());
    }

    // Consistency check for outputs
    for (int i=0; i<ret.n_out(); ++i) {
      const Sparsity& sp = input(i % n_in).sparsity();
      casadi_assert_message(ret.output(i).size()==sp.size(),
                            "Incorrect shape for " << ret << " output " << i << " \""
                            << o_names.at(i) << "\". Expected " << sp.size()
                            << " but got " << ret.output(i).size());
    }

    // Save to cache
    derivative_adj_[nadj] = ret;

    // Return generated function
    return ret;
  }

  void FunctionInternal::setDerForward(const Function& fcn, int nfwd) {

    // Check if there are enough forward directions allocated
    if (nfwd>=derivative_fwd_.size()) {
      derivative_fwd_.resize(nfwd+1);
    }

    // Save to cache
    derivative_fwd_[nfwd] = fcn;
  }

  void FunctionInternal::setDerReverse(const Function& fcn, int nadj) {

    // Check if there are enough adjoint directions allocated
    if (nadj>=derivative_adj_.size()) {
      derivative_adj_.resize(nadj+1);
    }

    // Save to cache
    derivative_adj_[nadj] = fcn;
  }

  Function FunctionInternal::getDerForward(const std::string& name, int nfwd, Dict& opts) {
    // TODO(@jaeandersson): Fallback on finite differences
    casadi_error("FunctionInternal::getDerForward not defined for class "
                 << typeid(*this).name());
  }

  Function FunctionInternal::getDerReverse(const std::string& name, int nadj, Dict& opts) {
    casadi_error("FunctionInternal::getDerReverse not defined for class "
                 << typeid(*this).name());
  }

  int FunctionInternal::nnz_in() const {
    int ret=0;
    for (int iind=0; iind<n_in(); ++iind) {
      ret += input(iind).nnz();
    }
    return ret;
  }

  int FunctionInternal::nnz_out() const {
    int ret=0;
    for (int oind=0; oind<n_out(); ++oind) {
      ret += output(oind).nnz();
    }
    return ret;
  }

  int FunctionInternal::numel_in() const {
    int ret=0;
    for (int iind=0; iind<n_in(); ++iind) {
      ret += input(iind).numel();
    }
    return ret;
  }

  int FunctionInternal::numel_out() const {
    int ret=0;
    for (int oind=0; oind<n_out(); ++oind) {
      ret += output(oind).numel();
    }
    return ret;
  }

  void FunctionInternal::call(const MXVector& arg, MXVector& res,
                              bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");

    // Lo logic for inlining yet
    bool inline_function = always_inline;

    if (inline_function) {
      // Evaluate the function symbolically
      evalMX(arg, res);

    } else {
      // Create a call-node
      assertInit();

      // Argument checking
      casadi_assert_message(
        arg.size()<=n_in(), "Function::call: number of passed-in dependencies ("
        << arg.size() << ") should not exceed the number of inputs of the function ("
        << n_in() << ").");

      // Assumes initialized
      for (int i=0; i<arg.size(); ++i) {
        if (arg[i].isNull() || arg[i].isempty() || input(i).isempty()) continue;
        casadi_assert_message(
          arg[i].size2()==input(i).size2() && arg[i].size1()==input(i).size1(),
          "Evaluation::shapes of passed-in dependencies should match shapes of inputs of function."
          << endl << ischeme_.at(i) <<  " has shape (" << input(i).size2()
          << ", " << input(i).size1() << ") while a shape (" << arg[i].size2() << ", "
          << arg[i].size1() << ") was supplied.");
      }
      res = Call::create(shared_from_this<Function>(), arg);
    }
  }

  void FunctionInternal::call(const SXVector& arg, SXVector& res,
                        bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    casadi_assert_message(!never_inline, "SX expressions do not support call-nodes");
    evalSX(arg, res);
  }

  void FunctionInternal::call(const DMatrixVector& arg, DMatrixVector& res,
                        bool always_inline, bool never_inline) {
    for (int i=0;i<arg.size();++i) {
      setInput(arg[i], i);
    }
    evaluate();
    res.resize(n_out());
    for (int i=0;i<res.size();++i) {
      res[i]=output(i);
    }
  }

  Function FunctionInternal::
  getNumericJacobian(const std::string& name, int iind, int oind, bool compact, bool symmetric,
                     const Dict& opts) {
    Function f = wrapMXFunction();
    f.init();
    return f->getNumericJacobian(name, iind, oind, compact, symmetric, opts);
  }

  Function FunctionInternal::fullJacobian() {
    if (full_jacobian_.alive()) {
      // Return cached Jacobian
      return shared_cast<Function>(full_jacobian_.shared());
    } else {
      // Options
      string name = name_ + "_jac";
      Dict opts;
      opts["input_scheme"] = ischeme_;
      opts["output_scheme"] = std::vector<std::string>(1, "jac");

      Function ret = getFullJacobian(name, opts);

      // Consistency check
      casadi_assert(ret.n_in()==n_in());
      casadi_assert(ret.n_out()==1);

      // Cache it for reuse and return
      full_jacobian_ = ret;
      return ret;
    }
  }

  Function FunctionInternal::getFullJacobian(const std::string& name, const Dict& opts) {
    casadi_assert(numDerForward()>0 || numDerReverse()>0);

    // Number inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Symbolic inputs of the full Jacobian function under construction
    vector<MX> ret_argv = symbolicInput(), argv, resv;

    // Symbolic input of the SISO function formed for generating the Jacobian
    MX arg;

    // Reuse the same symbolic primitives, if single input
    if (n_in==1) {
      argv = ret_argv;
      arg = argv.front();
      resv = symbolicOutput(argv);
    } else {
      // Collect Sparsity patterns of inputs
      vector<Sparsity> sp_argv(n_in);
      vector<int> row_offset(n_in+1, 0);
      for (int i=0; i<n_in; ++i) {
        sp_argv[i] = vec(input(i).sparsity());
        row_offset[i+1] = row_offset[i]+sp_argv[i].numel();
      }
      Sparsity sp_arg = vertcat(sp_argv);

      // Create symbolic primitive
      arg = MX::sym("x", sp_arg);

      // Split up and reshape to correct shape
      argv = vertsplit(arg, row_offset);
      for (int i=0; i<n_in; ++i) {
        argv[i] = reshape(argv[i], input(i).sparsity());
      }

      // Evaluate symbolically
      resv = shared_from_this<Function>()(argv);
    }

    // Reuse the same output, if possible
    MX res = n_out==1 ? resv.front() : veccat(resv);

    // Form Jacobian
    MX J;
    {
      MXFunction tmp("tmp", make_vector(arg), make_vector(res),
                     make_dict("ad_weight", adWeight()));
      J = tmp.jac();
    }

    // Make sure argv is the input of J
    if (n_in!=1) {
      J = substitute(J, arg, veccat(ret_argv));
    }

    // Form an expression for the full Jacobian
    return MXFunction(name, ret_argv, make_vector(J), opts);
  }

  void FunctionInternal::generateFunction(CodeGenerator& g,
                                          const std::string& fname, bool decl_static) const {
    assertInit();

    // Add standard math
    g.addInclude("math.h");

    // Add auxiliaries. TODO: Only add the auxiliaries that are actually used
    g.addAuxiliary(CodeGenerator::AUX_SQ);
    g.addAuxiliary(CodeGenerator::AUX_SIGN);

    // Generate declarations
    generateDeclarations(g);

    // Define function
    g.body << "/* " << getSanitizedName() << " */" << endl;
    if (decl_static) {
      g.body << "static ";
    } else if (g.cpp) {
      g.body << "extern \"C\" ";
    }
    g.body << "int " << fname << "(const real_t** arg, real_t** res, int* iw, real_t* w) {"
      << endl;

    // Insert the function body
    generateBody(g);

    // Finalize the function
    g.body
      << "  return 0;" << endl
      << "}" << endl
      << endl;
  }

  void FunctionInternal::generate(CodeGenerator& g,
                                  const vector<int>& arg, const vector<int>& res) const {
    if (simplifiedCall()) {

      // Collect input arguments
      for (int i=0; i<arg.size(); ++i) {
        g.body << "  w[" << i << "]=" << g.workel(arg[i]) << ";" << endl;
      }

      // Call function
      g.body << "  " << generateCall(g, "w", "w+"+g.to_string(arg.size())) << ";" << endl;

      // Collect output arguments
      for (int i=0; i<res.size(); ++i) {
        if (res[i]>=0) {
          g.body << "  " << g.workel(res[i]) << "=w[" << (arg.size()+i) << "];" << endl;
        }
      }
    } else {

      // Collect input arguments
      for (int i=0; i<arg.size(); ++i) {
        g.body << "  arg1[" << i << "]=" << g.work(arg[i], input(i).nnz()) << ";" << endl;
      }

      // Collect output arguments
      for (int i=0; i<res.size(); ++i) {
        g.body << "  res1[" << i << "]=" << g.work(res[i], output(i).nnz()) << ";" << endl;
      }

      // Call function
      g.body << "  if (" << generateCall(g, "arg1", "res1", "iw", "w") << ") return 1;" << endl;
    }
  }

  void FunctionInternal::generateMeta(CodeGenerator& g, const std::string& fname) const {
    // Short-hands
    int n_in = this->n_in();
    int n_out = this->n_out();
    stringstream &s = g.body;

    // Function that returns the number of inputs and outputs
    string tmp = "int " + fname + "_init(int *f_type, int *n_in, int *n_out, "
      "int *sz_arg, int* sz_res)";
    if (g.cpp) {
      tmp = "extern \"C\" " + tmp;  // C linkage
    }
    if (g.with_header) {
      g.header << tmp << ";" << endl;
    }
    s << tmp << " {" << endl
      << "  *f_type = 1;" << endl
      << "  *n_in = " << n_in << ";" << endl
      << "  *n_out = " << n_out << ";" << endl
      << "  *sz_arg = " << sz_arg() << ";" << endl
      << "  *sz_res = " << sz_res() << ";" << endl
      << "  return 0;" << endl
      << "}" << endl
      << endl;

    // Inputs and outputs for each parsity index
    std::multimap<int, int> io_sparsity_index;

    // The number of patterns needed for the inputs and outputs
    int num_io_patterns = 0;

    // Get the sparsity pattern of all inputs
    for (int i=0; i<n_in+n_out; ++i) {
      // Get the sparsity pattern
      const Sparsity& sp = i<n_in ? input(i).sparsity() : output(i-n_in).sparsity();

      // Print the sparsity pattern or retrieve the index of an existing pattern
      int ind = g.addSparsity(sp);
      num_io_patterns = std::max(num_io_patterns, ind+1);

      // Store the index
      io_sparsity_index.insert(std::pair<int, int>(ind, i));
    }

    // Function that returns the sparsity pattern
    tmp = "int " + fname + "_sparsity"
      "(int i, int *nrow, int *ncol, const int **colind, const int **row)";
    if (g.cpp) {
      tmp = "extern \"C\" " + tmp;  // C linkage
    }
    if (g.with_header) {
      g.header << tmp << ";" << endl;
    }
    s << tmp << " {" << endl;

    // Get the sparsity index using a switch
    s << "  const int* s;" << endl;
    s << "  switch (i) {" << endl;

    // Loop over all sparsity patterns
    for (int i=0; i<num_io_patterns; ++i) {
      // Get the range of matching sparsity patterns
      typedef std::multimap<int, int>::const_iterator it_type;
      std::pair<it_type, it_type> r = io_sparsity_index.equal_range(i);

      // Print the cases covered
      for (it_type it=r.first; it!=r.second; ++it) {
        s << "    case " << it->second << ":" << endl;
      }

      // Map to sparsity
      s << "      s = s" << i << "; break;" << endl;
    }

    // Finalize the switch
    s << "    default:" << endl;
    s << "      return 1;" << endl;
    s << "  }" << endl << endl;

    // Decompress the sparsity pattern
    s << "  *nrow = s[0];" << endl;
    s << "  *ncol = s[1];" << endl;
    s << "  *colind = s + 2;" << endl;
    s << "  *row = s + 2 + (*ncol + 1);" << endl;
    s << "  return 0;" << endl;
    s << "}" << endl << endl;

    // Function that returns work vector lengths
    tmp = "int " + fname + "_work(int *sz_iw, int *sz_w)";
    if (g.cpp) {
      tmp = "extern \"C\" " + tmp;  // C linkage
    }
    if (g.with_header) {
      g.header << tmp << ";" << endl;
    }
    s << tmp << " {" << endl;
    s << "  if (sz_iw) *sz_iw = " << sz_iw() << ";" << endl;
    s << "  if (sz_w) *sz_w = " << sz_w() << ";" << endl;
    s << "  return 0;" << endl;
    s << "}" << endl;
    s << endl;

    // Generate mex gateway for the function
    if (g.mex) {
      // Begin conditional compilation
      s << "#ifdef MATLAB_MEX_FILE" << endl;

      // Declare wrapper
      s << "void mex_" << fname
             << "(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {" << endl
             << "  int i, j;" << endl;

      // Check arguments
      s << "  if (argc>" << n_in << ") mexErrMsgIdAndTxt(\"Casadi:RuntimeError\","
             << "\"Evaluation of \\\"" << fname << "\\\" failed. Too many input arguments "
             << "(%d, max " << n_in << ")\", argc);" << endl;

      s << "  if (resc>" << n_out << ") mexErrMsgIdAndTxt(\"Casadi:RuntimeError\","
             << "\"Evaluation of \\\"" << fname << "\\\" failed. "
             << "Too many output arguments (%d, max " << n_out << ")\", resc);" << endl;

      // Work vectors, including input and output buffers
      int i_nnz = nnz_in(), o_nnz = nnz_out();
      size_t sz_w = this->sz_w();
      for (int i=0; i<n_in; ++i) {
        const Sparsity& s = input(i).sparsity();
        sz_w = max(sz_w, static_cast<size_t>(s.size1())); // To be able to copy a column
        sz_w = max(sz_w, static_cast<size_t>(s.size2())); // To be able to copy a row
      }
      sz_w += i_nnz + o_nnz;
      s << "  int iw[" << sz_iw() << "];" << endl;
      s << "  real_t w[" << sz_w << "];" << endl;
      string fw = "w+" + g.to_string(i_nnz + o_nnz);

      // Copy inputs to buffers
      int offset=0;
      s << "  const real_t* arg[" << n_in << "] = {0};" << endl;
      for (int i=0; i<n_in; ++i) {
        std::string p = "argv[" + g.to_string(i) + "]";
        s << "  if (--argc>=0) arg[" << i << "] = "
               << g.from_mex(p, "w", offset, input(i).sparsity(), fw) << endl;
        offset += input(i).nnz();
      }

      // Allocate output buffers
      s << "  real_t* res[" << n_out << "] = {0};" << endl;
      for (int i=0; i<n_out; ++i) {
        if (i==0) {
          // if i==0, always store output (possibly ans output)
          s << "  --resc;" << endl
            << "  ";
        } else {
          // Store output, if it exists
          s << "  if (--resc>=0) ";
        }
        // Create and get pointer
        s << "res[" << i << "] = w+" << g.to_string(offset) << ";" << endl;
        offset += output(i).nnz();
      }

      // Call the function
      s << "  i = " << fname << "(arg, res, iw, " << fw << ");" << endl;
      s << "  if (i) mexErrMsgIdAndTxt(\"Casadi:RuntimeError\",\"Evaluation of \\\"" << fname
             << "\\\" failed.\");" << endl;

      // Save results
      for (int i=0; i<n_out; ++i) {
        string res_i = "res[" + g.to_string(i) + "]";
        s << "  if (" << res_i << ") resv[" << i << "] = "
          << g.to_mex(output(i).sparsity(), res_i) << endl;
      }

      // End conditional compilation and function
      s << "}" << endl
        << "#endif" << endl << endl;
    }

    if (g.main) {
      // Declare wrapper
      s << "int main_" << fname << "(int argc, char* argv[]) {" << endl;

      // Work vectors and input and output buffers
      size_t nr = sz_w() + nnz_in() + nnz_out();
      s << "  int iw[" << sz_iw() << "];" << endl
             << "  real_t w[" << nr << "];" << endl;

      // Input buffers
      s << "  const real_t* arg[" << n_in << "] = {";
      int off=0;
      for (int i=0; i<n_in; ++i) {
        if (i!=0) s << ", ";
        s << g.work(off, input(i).nnz());
        off += input(i).nnz();
      }
      s << "};" << endl;

      // Output buffers
      s << "  real_t* res[" << n_out << "] = {";
      for (int i=0; i<n_out; ++i) {
        if (i!=0) s << ", ";
        s << g.work(off, output(i).nnz());
        off += output(i).nnz();
      }
      s << "};" << endl;

      // TODO(@jaeandersson): Read inputs from file. For now; read from stdin
      s << "  int j;" << endl
             << "  real_t* a = w;" << endl
             << "  for (j=0; j<" << nnz_in() << "; ++j) "
             << "scanf(\"%lf\", a++);" << endl;

      // Call the function
      s << "  int flag = eval(arg, res, iw, w+" << off << ");" << endl
             << "  if (flag) return flag;" << endl;

      // TODO(@jaeandersson): Write outputs to file. For now: print to stdout
      s << "  const real_t* r = w+" << nnz_in() << ";" << endl
             << "  for (j=0; j<" << nnz_out() << "; ++j) "
             << g.printf("%g ", "*r++") << endl;
      // End with newline
      s << "  " << g.printf("\\n") << endl;

      // Finalize function
      s << "  return 0;" << endl
        << "}" << endl << endl;
    }
  }

  std::string FunctionInternal::generateCall(const CodeGenerator& g,
                                             const std::string& arg, const std::string& res,
                                             const std::string& iw, const std::string& w) const {
    // Get the index of the function
    CodeGenerator::PointerMap::const_iterator it=g.added_dependencies_.find(this);
    casadi_assert(it!=g.added_dependencies_.end());
    int f = it->second;

    // Create a function call
    stringstream ss;
    ss << "f" << f << "(" << arg << ", " << res << ", " << iw << ", " << w << ")";
    return ss.str();
  }

  std::string FunctionInternal::generateCall(const CodeGenerator& g,
                                             const std::string& arg,
                                             const std::string& res) const {
    casadi_error("FunctionInternal::generateCall (simplified): not defined for class "
                 << typeid(*this).name());
  }


  void FunctionInternal::addDependency(CodeGenerator& g) const {
    // Get the current number of functions before looking for it
    size_t num_f_before = g.added_dependencies_.size();

    // Get index of the pattern
    int& ind = g.added_dependencies_[this];

    // Generate it if it does not exist
    if (g.added_dependencies_.size() > num_f_before) {
      // Add at the end
      ind = num_f_before;

      // Give it a name
      string name = "f" + CodeGenerator::to_string(ind);

      // Print to file
      generateFunction(g, "CASADI_PREFIX(" + name + ")", true);

      // Shorthand
      g.body
        << "#define " << name << "(arg, res, iw, w) "
        << "CASADI_PREFIX(" << name << ")(arg, res, iw, w)" << endl << endl;
    }
  }

  void FunctionInternal::generateDeclarations(CodeGenerator& g) const {
    // Nothing to declare
  }

  void FunctionInternal::generateBody(CodeGenerator& g) const {
    casadi_error("FunctionInternal::generateBody: generateBody not defined for class "
                 << typeid(*this).name());
  }

  Function FunctionInternal::dynamicCompilation(Function f, std::string fname, std::string fdescr,
                                                std::string compiler) {
    // Check if f is initialized
    bool f_is_init = f.isInit();
    if (!f_is_init) f.init();

    // Codegen and compile
    CodeGenerator g;
    g.add(f, fname);
    string dlname = g.compile(fname, compiler);

    // Load it
    ExternalFunction f_gen(fname, dlname);

    // Initialize it if f was initialized
    if (f_is_init) {
      f_gen.init();
      if (verbose_) {
        userOut() << "Dynamically loaded " << fdescr << " (" << fname << ")" << endl;
      }
    }
    return f_gen;
  }

  void FunctionInternal::spFwdSwitch(const bvec_t** arg, bvec_t** res,
                                     int* iw, bvec_t* w) {
    // TODO(@jaeandersson) Calculate from full-Jacobian sparsity  when necessary or more efficient
    spFwd(arg, res, iw, w);
  }

  void FunctionInternal::spFwd(const bvec_t** arg, bvec_t** res,
                               int* iw, bvec_t* w) {
    // Number inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Pass/clear forward seeds/adjoint sensitivities
    for (int i=0; i<n_in; ++i) {
      // Input vector
      bvec_t* input_i = reinterpret_cast<bvec_t*>(input(i).ptr());
      if (arg[i]==0) {
        // Set to zero if not used
        fill_n(input_i, input(i).nnz(), 0);
      } else {
        copy(arg[i], arg[i]+input(i).nnz(), input_i);
      }
    }

    // Pass/clear adjoint seeds/forward sensitivities
    for (int i=0; i<n_out; ++i) {
      // Output vector
      bvec_t* output_i = reinterpret_cast<bvec_t*>(output(i).ptr());
      if (res[i]==0) {
        // Set to zero if not used
        fill_n(output_i, output(i).nnz(), 0);
      } else {
        copy(res[i], res[i]+output(i).nnz(), output_i);
      }
    }

    // Propagate seeds
    spInit(true); // NOTE: should only be done once
    if (spCanEvaluate(true)) {
      spEvaluate(true);
    } else {
      spEvaluateViaJacSparsity(true);
    }

    // Get the sensitivities
    for (int i=0; i<n_out; ++i) {
      if (res[i]!=0) {
        const bvec_t* output_i = reinterpret_cast<const bvec_t*>(output(i).ptr());
        copy(output_i, output_i+output(i).nnz(), res[i]);
      }
    }

    // Clear seeds and sensitivities
    for (int i=0; i<n_in; ++i) input(i).set(0.);
    for (int i=0; i<n_out; ++i) output(i).set(0.);
  }

  void FunctionInternal::spAdjSwitch(bvec_t** arg, bvec_t** res,
                                     int* iw, bvec_t* w) {
    // TODO(@jaeandersson) Calculate from full-Jacobian sparsity  when necessary or more efficient
    spAdj(arg, res, iw, w);
  }

  void FunctionInternal::spAdj(bvec_t** arg, bvec_t** res,
                               int* iw, bvec_t* w) {
    // Number inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Pass/clear forward seeds/adjoint sensitivities
    for (int i=0; i<n_in; ++i) {
      // Input vector
      bvec_t* input_i = reinterpret_cast<bvec_t*>(input(i).ptr());
      if (arg[i]==0) {
        // Set to zero if not used
        fill_n(input_i, input(i).nnz(), 0);
      } else {
        copy(arg[i], arg[i]+input(i).nnz(), input_i);
      }
    }

    // Pass/clear adjoint seeds/forward sensitivities
    for (int i=0; i<n_out; ++i) {
      // Output vector
      bvec_t* output_i = reinterpret_cast<bvec_t*>(output(i).ptr());
      if (res[i]==0) {
        // Set to zero if not used
        fill_n(output_i, output(i).nnz(), 0);
      } else {
        copy(res[i], res[i]+output(i).nnz(), output_i);
      }
    }

    // Propagate seeds
    spInit(false); // NOTE: should only be done once
    if (spCanEvaluate(false)) {
      spEvaluate(false);
    } else {
      spEvaluateViaJacSparsity(false);
    }

    // Get the sensitivities
    for (int i=0; i<n_in; ++i) {
      if (arg[i]!=0) {
        int n = input(i).nnz();
        bvec_t* arg_i = arg[i];
        const bvec_t* input_i = reinterpret_cast<const bvec_t*>(input(i).ptr());
        for (int k=0; k<n; ++k) *arg_i++ |= *input_i++;
      }
    }

    // Clear seeds and sensitivities
    for (int i=0; i<n_in; ++i) input(i).set(0.);
    for (int i=0; i<n_out; ++i) output(i).set(0.);
  }

  void FunctionInternal::sz_work(size_t& sz_arg, size_t& sz_res,
                                 size_t& sz_iw, size_t& sz_w) const {
    sz_arg = sz_arg_;
    sz_res = sz_res_;
    sz_iw = sz_iw_;
    sz_w = sz_w_;
  }

  void FunctionInternal::alloc_arg(size_t sz_arg) {
    sz_arg_ = max(sz_arg + n_in(), sz_arg_);
  }

  void FunctionInternal::alloc_res(size_t sz_res) {
    sz_res_ = max(sz_res + n_out(), sz_res_);
  }

  void FunctionInternal::alloc_iw(size_t sz_iw) {
    sz_iw_ = max(sz_iw, sz_iw_);
  }

  void FunctionInternal::alloc_w(size_t sz_w) {
    sz_w_ = max(sz_w, sz_w_);
  }

  void FunctionInternal::alloc(const Function& f) {
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f.sz_work(sz_arg, sz_res, sz_iw, sz_w);
    alloc_arg(sz_arg);
    alloc_res(sz_res);
    alloc_iw(sz_iw);
    alloc_w(sz_w);
  }

  void FunctionInternal::alloc() {
    iw_tmp_.resize(sz_iw_);
    w_tmp_.resize(sz_w_);
  }

  void FunctionInternal::reportConstraints(std::ostream &stream, const Matrix<double> &v,
                                           const Matrix<double> &lb, const Matrix<double> &ub,
                                           const std::string &name, double tol) {

    casadi_assert_message(v.sparsity()==lb.sparsity(),
                          "reportConstraints:: sparsities must match");
    casadi_assert_message(ub.sparsity()==lb.sparsity(),
                          "reportConstraints:: sparsities must match");

    // Backup the formatting
    ios::fmtflags fmtflags_backup = stream.flags();
    int streamsize_backup = stream.precision();

    //stream.setf(ios::fixed, ios::floatfield);
    //stream.precision(8);

    // Check if any constraint is violated
    bool all_ok = true;
    for (int i=0; all_ok && i<v.nnz(); ++i) {
      all_ok = v.at(i) > ub.at(i) + tol || v.at(i) < lb.at(i) - tol;
    }
    if (all_ok) {
      stream << "All " << v.nnz() << " constraints on " << name << " are met: " << endl;
    } else {
      stream << "Problem with constraints on " << name << ": " << endl;
    }

    // Make a horizontal rule
    stream.width(60);
    stream.fill('-');
    stream  << "-" << endl;
    stream.fill(' ');

    // The length of the numeric fields
    int fieldlength = 10;
    // The length of the constraint visualizer strip
    int indicator_length = 15;

    // Loop over the elements of v
    for (int i=0;i<v.nnz();i++) {

      stream.width(5);
      stream << i << ". |   ";

      if (fabs(lb.at(i) - ub.at(i))<=tol) {
        stream.width(fieldlength);
        stream << lb.at(i) << " ==  ";
        stream.width(fieldlength);
        stream << v.at(i) << "      ";
        stream.width(fieldlength);
        stream << " ";
        stream << "  |   ";
      } else {
        // BEGIN  - construct the constraint visualizer strip
        std::string indicator(indicator_length+2, '-');
        indicator.at(0) = (fabs(v.at(i)-lb.at(i))<=tol)? 'X' : 'o';
        if (lb.at(i)==-std::numeric_limits<double>::infinity()) indicator.at(0)='8';

        indicator.at(indicator_length+1) = (fabs(v.at(i)-ub.at(i))<=tol)? 'X' : 'o';
        if (ub.at(i)==std::numeric_limits<double>::infinity()) indicator.at(indicator_length+1)='8';

        if (v.at(i) <= (ub.at(i) + tol) && v.at(i) >= (lb.at(i) - tol)) {
          int index = (v.at(i)-lb.at(i))/(ub.at(i)-lb.at(i))*(indicator_length-1);
          index = min(max(0, index), indicator_length-1);
          indicator.at(1+index) = '=';
        }
        // END - construct the constraint visualizer strip

        stream.width(fieldlength);
        stream << lb.at(i) << " <=  ";
        stream.width(fieldlength);
        stream << v.at(i) << " <= ";
        stream.width(fieldlength);
        stream << ub.at(i) << "    | ";
        if (v.at(i) <= (ub.at(i) + tol) && v.at(i) >= (lb.at(i) - tol)) {
          stream  << indicator;
        }
      }

      if (v.at(i) <= (ub.at(i) + tol) && v.at(i) >= (lb.at(i) - tol)) {
      } else {
        stream  << "  VIOLATED";
      }
      stream  << endl;
    }

    // Make a horizontal rule
    stream.width(60);
    stream.fill('-');
    stream  << "-" << endl;
    stream.fill(' ');

    // Restore the formatting
    stream.setf(fmtflags_backup);
    stream.precision(streamsize_backup);
  }

  // helper function for getSanitizedName
  bool isBadChar(char c) {
    return !std::isalnum(c);
  }

  std::string FunctionInternal::getSanitizedName() const {
      return sanitizeName(name_);
  }

  std::string FunctionInternal::sanitizeName(const std::string &name) {
    string sname = name;
    std::replace_if(sname.begin(), sname.end(), isBadChar, '_');
    return sname;
  }

  bool FunctionInternal::hasFullJacobian() const {
    return full_jacobian_.alive();
  }

  bool FunctionInternal::hasDerivative() const {
    return numDerForward()>0 || numDerReverse()>0 || hasFullJacobian();
  }

  bool FunctionInternal::fwdViaJac(int nfwd) {
    if (numDerForward()==0) return true;

    // Jacobian calculation penalty factor
    const double jac_penalty = getOption("jac_penalty");

    if (jac_penalty==-1) return false;

    // Heuristic 1: Jac calculated via forward mode likely cheaper
    if (jac_penalty*nnz_in()<nfwd) return true;

    // Heuristic 2: Jac calculated via reverse mode likely cheaper
    double w = adWeight();
    if (numDerReverse()>0 && jac_penalty*(1-w)*nnz_out()<w*nfwd)
      return true;

    return false;
  }

  bool FunctionInternal::adjViaJac(int nadj) {
    if (numDerReverse()==0) return true;

    // Jacobian calculation penalty factor
    const double jac_penalty = getOption("jac_penalty");

    if (jac_penalty==-1) return false;

    // Heuristic 1: Jac calculated via reverse mode likely cheaper
    if (jac_penalty*nnz_out()<nadj) return true;

    // Heuristic 2: Jac calculated via forward mode likely cheaper
    double w = adWeight();
    if (numDerForward()>0 && jac_penalty*w*nnz_in()<(1-w)*nadj)
      return true;

    return false;
  }

  void FunctionInternal::callForward(const std::vector<MX>& arg, const std::vector<MX>& res,
                                 const std::vector<std::vector<MX> >& fseed,
                                 std::vector<std::vector<MX> >& fsens,
                                 bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    casadi_assert_message(!always_inline, "Class " << typeid(*this).name() <<
                          " cannot be inlined in an MX expression");

    // Derivative information must be available
    casadi_assert(hasDerivative());

    // Number of directional derivatives
    int nfwd = fseed.size();
    fsens.resize(nfwd);

    // Quick return if no seeds
    if (nfwd==0) return;

    // Number inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Calculating full Jacobian and then multiplying
    if (fwdViaJac(nfwd)) {
      // Join forward seeds
      vector<MX> v(nfwd);
      for (int d=0; d<nfwd; ++d) {
        v[d] = veccat(fseed[d]);
      }

      // Multiply the Jacobian from the right
      MX J = fullJacobian()(arg).at(0);
      v = horzsplit(mul(J, horzcat(v)));

      // Vertical offsets
      vector<int> offset(n_out+1, 0);
      for (int i=0; i<n_out; ++i) {
        offset[i+1] = offset[i]+output(i).numel();
      }

      // Collect forward sensitivities
      for (int d=0; d<nfwd; ++d) {
        fsens[d] = vertsplit(v[d], offset);
        for (int i=0; i<n_out; ++i) {
          fsens[d][i] = reshape(fsens[d][i], output(i).size());
        }
      }
      return;
    }

    // All inputs and seeds
    vector<MX> darg;
    darg.reserve(n_in + n_out + n_in*nfwd);
    darg.insert(darg.end(), arg.begin(), arg.end());
    darg.insert(darg.end(), res.begin(), res.end());
    for (int d=0; d<nfwd; ++d) {
      darg.insert(darg.end(), fseed[d].begin(), fseed[d].end());
    }

    // Create derivative function
    Function dfcn = derForward(nfwd);

    // Create the evaluation node
    vector<MX> x = Call::create(dfcn, darg);
    vector<MX>::iterator x_it = x.begin();

    // Retrieve sensitivities
    for (int d=0; d<nfwd; ++d) {
      fsens[d].resize(n_out);
      for (int i=0; i<n_out; ++i) {
        fsens[d][i] = *x_it++;
      }
    }
    casadi_assert(x_it==x.end());
  }

  void FunctionInternal::callReverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                                 const std::vector<std::vector<MX> >& aseed,
                                 std::vector<std::vector<MX> >& asens,
                                 bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    casadi_assert_message(!always_inline, "Class " << typeid(*this).name() <<
                          " cannot be inlined in an MX expression");

    // Derivative information must be available
    casadi_assert(hasDerivative());

    // Number of directional derivatives
    int nadj = aseed.size();
    asens.resize(nadj);

    // Quick return if no seeds
    if (nadj==0) return;

    // Number inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Calculating full Jacobian and then multiplying likely cheaper
    if (adjViaJac(nadj)) {
      // Join adjoint seeds
      vector<MX> v(nadj);
      for (int d=0; d<nadj; ++d) {
        v[d] = veccat(aseed[d]);
      }

      // Multiply the transposed Jacobian from the right
      MX J = fullJacobian()(arg).at(0);
      v = horzsplit(mul(J.T(), horzcat(v)));

      // Vertical offsets
      vector<int> offset(n_in+1, 0);
      for (int i=0; i<n_in; ++i) {
        offset[i+1] = offset[i]+input(i).numel();
      }

      // Collect adjoint sensitivities
      for (int d=0; d<nadj; ++d) {
        asens[d].resize(n_in);
        vector<MX> a = vertsplit(v[d], offset);
        for (int i=0; i<n_in; ++i) {
          if (asens[d][i].isempty(true)) {
            asens[d][i] = reshape(a[i], input(i).size());
          } else {
            asens[d][i] += reshape(a[i], input(i).size());
          }
        }
      }
      return;
    }

    // All inputs and seeds
    vector<MX> darg;
    darg.reserve(n_in + n_out + n_out*nadj);
    darg.insert(darg.end(), arg.begin(), arg.end());
    darg.insert(darg.end(), res.begin(), res.end());
    for (int d=0; d<nadj; ++d) {
      darg.insert(darg.end(), aseed[d].begin(), aseed[d].end());
    }

    // Create derivative function
    Function dfcn = derReverse(nadj);

    // Create the evaluation node
    vector<MX> x = Call::create(dfcn, darg);
    vector<MX>::iterator x_it = x.begin();

    // Retrieve sensitivities
    for (int d=0; d<nadj; ++d) {
      asens[d].resize(n_in);
      for (int i=0; i<n_in; ++i) {
        if (asens[d][i].isempty(true)) {
          asens[d][i] = *x_it++;
        } else {
          asens[d][i] += *x_it++;
        }
      }
    }
    casadi_assert(x_it==x.end());
  }

  void FunctionInternal::callForward(const std::vector<SX>& arg, const std::vector<SX>& res,
                                 const std::vector<std::vector<SX> >& fseed,
                                 std::vector<std::vector<SX> >& fsens,
                                 bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    if (fseed.empty()) { // Quick return if no seeds
      fsens.clear();
      return;
    }
    casadi_error("FunctionInternal::callForward(SX) not defined for class "
                 << typeid(*this).name());
  }

  void FunctionInternal::callReverse(const std::vector<SX>& arg, const std::vector<SX>& res,
                                 const std::vector<std::vector<SX> >& aseed,
                                 std::vector<std::vector<SX> >& asens,
                                 bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    if (aseed.empty()) { // Quick return if no seeds
      asens.clear();
      return;
    }
    casadi_error("FunctionInternal::callReverse(SX) not defined for class "
                 << typeid(*this).name());
  }

  void FunctionInternal::
  callForward(const std::vector<DMatrix>& arg, const std::vector<DMatrix>& res,
          const std::vector<std::vector<DMatrix> >& fseed,
          std::vector<std::vector<DMatrix> >& fsens,
          bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    casadi_error("DMatrix does not support call-nodes");

    // TODO(@jaeandersson): Evaluate in "batches"

    // Retrieve derivative function
    int nfwd = fseed.size();
    Function dfcn = derForward(nfwd);

    // Pass inputs
    int n_in = this->n_in();
    int di=0;
    casadi_assert(arg.size()==n_in);
    for (int i=0; i<n_in; ++i) {
      dfcn.setInput(arg[i], di++);
    }

    // Pass outputs
    int n_out = this->n_out();
    casadi_assert(res.size()==n_out);
    for (int i=0; i<n_out; ++i) {
      dfcn.setInput(res[i], di++);
    }

    // Pass seeds
    for (int d=0; d<nfwd; ++d) {
      casadi_assert(fseed[d].size()==n_in);
      for (int i=0; i<n_in; ++i) {
        dfcn.setInput(fseed[d][i], di++);
      }
    }
    casadi_assert(di==dfcn.n_in()); // Consistency check

    // Calculate derivatives
    dfcn.evaluate();

    // Get sensitivities
    fsens.resize(nfwd);
    di = 0;
    for (int d=0; d<nfwd; ++d) {
      fsens[d].resize(n_out);
      for (int i=0; i<n_out; ++i) {
        fsens[d][i] = dfcn.output(di++);
      }
    }
    casadi_assert(di==dfcn.n_out()); // Consistency check
  }

  void FunctionInternal::
  callReverse(const std::vector<DMatrix>& arg, const std::vector<DMatrix>& res,
          const std::vector<std::vector<DMatrix> >& aseed,
          std::vector<std::vector<DMatrix> >& asens,
          bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    casadi_error("DMatrix does not support call-nodes");

    // TODO(@jaeandersson): Evaluate in "batches"

    // Retrieve derivative function
    int nadj = aseed.size();
    Function dfcn = derReverse(nadj);

    // Pass inputs
    int n_in = this->n_in();
    int di=0;
    casadi_assert(arg.size()==n_in);
    for (int i=0; i<n_in; ++i) {
      dfcn.setInput(arg[i], di++);
    }

    // Pass outputs
    int n_out = this->n_out();
    casadi_assert(res.size()==n_out);
    for (int i=0; i<n_out; ++i) {
      dfcn.setInput(res[i], di++);
    }

    // Pass seeds
    for (int d=0; d<nadj; ++d) {
      casadi_assert(aseed[d].size()==n_out);
      for (int i=0; i<n_out; ++i) {
        dfcn.setInput(aseed[d][i], di++);
      }
    }
    casadi_assert(di==dfcn.n_in()); // Consistency check

    // Calculate derivatives
    dfcn.evaluate();

    // Get sensitivities
    asens.resize(nadj);
    di = 0;
    for (int d=0; d<nadj; ++d) {
      asens[d].resize(n_in);
      for (int i=0; i<n_in; ++i) {
        asens[d][i] = dfcn.output(di++);
      }
    }
    casadi_assert(di==dfcn.n_out()); // Consistency check
  }

  double FunctionInternal::adWeight() {
    // If reverse mode derivatives unavailable, use forward
    if (numDerReverse()==0) return 0;

    // If forward mode derivatives unavailable, use reverse
    if (numDerForward()==0) return 1;

    // A user-provided option overrides default value
    if (hasSetOption("ad_weight")) return getOption("ad_weight");

    // By default, reverse mode is about twice as expensive as forward mode
    return 0.33;  // i.e. nf <= 2*na <=> 1/3*nf <= (1-1/3)*na, forward when tie
  }

  double FunctionInternal::adWeightSp() {
    // If reverse mode propagation unavailable, use forward
    if (!spCanEvaluate(false)) return 0;

    // If forward mode propagation unavailable, use reverse
    if (!spCanEvaluate(true)) return 1;

    // A user-provided option overrides default value
    if (hasSetOption("ad_weight_sp")) return getOption("ad_weight_sp");

    // Both modes equally expensive by default (no "taping" needed)
    return 0.49; // Forward when tie
  }

  MX FunctionInternal::grad_mx(int iind, int oind) {
    casadi_error("FunctionInternal::grad_mx: not defined for class "
                 << typeid(*this).name());
  }

  MX FunctionInternal::tang_mx(int iind, int oind) {
    casadi_error("FunctionInternal::tang_mx: not defined for class "
                 << typeid(*this).name());
  }

  MX FunctionInternal::jac_mx(int iind, int oind, bool compact, bool symmetric,
                              bool always_inline, bool never_inline) {
    casadi_error("FunctionInternal::jac_mx: not defined for class "
                 << typeid(*this).name());
  }

  SX FunctionInternal::grad_sx(int iind, int oind) {
    casadi_error("FunctionInternal::grad_sx: not defined for class "
                 << typeid(*this).name());
  }

  SX FunctionInternal::tang_sx(int iind, int oind) {
    casadi_error("FunctionInternal::tang_sx: not defined for class "
                 << typeid(*this).name());
  }

  SX FunctionInternal::jac_sx(int iind, int oind, bool compact, bool symmetric,
                              bool always_inline, bool never_inline) {
    casadi_error("FunctionInternal::jac_sx: not defined for class "
                 << typeid(*this).name());
  }

} // namespace casadi
