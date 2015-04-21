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
#include "../mx/casadi_map.hpp"
#include <typeinfo>
#include "../std_vector_tools.hpp"
#include "mx_function.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"
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
  FunctionInternal::FunctionInternal() {
    setOption("name", "unnamed_function"); // name of the function
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
    //addOption("ad_mode",                  OT_STRING,              "automatic",
    //          "Deprecated option, use \"ad_weight\" instead. Ignored.");
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
    addOption("custom_forward",  OT_DERIVATIVEGENERATOR,   GenericType(),
              "Function that returns a derivative function given a number of forward "
              "mode directional derivatives. Overrides default routines.");
    addOption("custom_reverse",  OT_DERIVATIVEGENERATOR,   GenericType(),
              "Function that returns a derivative function given a number of reverse "
              "mode directional derivatives. Overrides default routines.");
    addOption("full_jacobian",                 OT_FUNCTION,              GenericType(),
              "The Jacobian of all outputs with respect to all inputs.");

    verbose_ = false;
    user_data_ = 0;
    monitor_inputs_ = false;
    monitor_outputs_ = false;
  }


  FunctionInternal::~FunctionInternal() {
  }

  void FunctionInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    OptionsFunctionalityNode::deepCopyMembers(already_copied);
    for (vector<WeakRef>::iterator j=derivative_fwd_.begin(); j!=derivative_fwd_.end(); ++j) {
      if (!j->isNull()) *j = getcopy(j->shared(), already_copied);
    }
    for (vector<WeakRef>::iterator j=derivative_adj_.begin(); j!=derivative_adj_.end(); ++j) {
      if (!j->isNull()) *j = getcopy(j->shared(), already_copied);
    }


    full_jacobian_ = getcopy(full_jacobian_, already_copied);
  }

  void FunctionInternal::init() {
    verbose_ = getOption("verbose");
    regularity_check_ = getOption("regularity_check");

    // Warn for functions with too many inputs or outputs
    casadi_assert_warning(getNumInputs()<10000, "Function " << getOption("name")
                          << " has a large number of inputs. "
                          "Changing the problem formulation is strongly encouraged.");
    casadi_assert_warning(getNumOutputs()<10000, "Function " << getOption("name")
                          << " has a large number of outputs. "
                          "Changing the problem formulation is strongly encouraged.");

    // Resize the matrix that holds the sparsity of the Jacobian blocks
    jac_sparsity_ = jac_sparsity_compact_ =
        SparseStorage<Sparsity>(Sparsity(getNumOutputs(), getNumInputs()));
    jac_ = jac_compact_ = SparseStorage<WeakRef>(Sparsity(getNumOutputs(), getNumInputs()));

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

    monitor_inputs_ = monitored("inputs");
    monitor_outputs_ = monitored("outputs");

    gather_stats_ = getOption("gather_stats");

    inputs_check_ = getOption("inputs_check");

    // Mark the function as initialized
    is_init_ = true;
  }

  void FunctionInternal::print(ostream &stream) const {
    if (getNumInputs()==1) {
      stream << " Input: " << input().dimString() << endl;
    } else {
      if (input_.scheme.isNull()) {
        stream << " Inputs (" << getNumInputs() << "):" << endl;
        for (int i=0;i<getNumInputs();i++) {
          stream << "  " << i << ". " << input(i).dimString() << endl;
        }
      } else {
        stream << " Inputs (" << input_.scheme.name() << ": " << getNumInputs()
               << "):" << endl;
        for (int i=0;i<getNumInputs();i++) {
          stream << "  " << i  << ". (" << input_.scheme.describe(i) << ")   "
                 << input(i).dimString() << endl;
        }
      }
    }
    if (getNumOutputs()==1) {
      stream << " Output: " << output().dimString() << endl;
    } else {
      if (output_.scheme.isNull()) {
        stream << " Outputs (" << getNumOutputs() << "):" << endl;
        for (int i=0;i<getNumOutputs();i++) {
          stream << "  " << i << ". " << output(i).dimString() << endl;
        }
      } else {
        stream << " Outputs (" << output_.scheme.name() << ": "
               << getNumOutputs() << "):" << endl;
        for (int i=0;i<getNumOutputs();i++) {
          stream << "  " << i << ". (" << output_.scheme.describe(i) << ")   "
                 << output(i).dimString() << endl;
        }
      }
    }
  }

  void FunctionInternal::repr(ostream &stream) const {
    stream << "function(\"" << getOption("name") << "\")";
  }

  Function FunctionInternal::gradient(int iind, int oind) {
    // Assert scalar
    casadi_assert_message(output(oind).isScalar(),
                          "Only gradients of scalar functions allowed. Use jacobian instead.");

    // Generate gradient function
    Function ret = getGradient(iind, oind);

    // Give it a suitable name
    stringstream ss;
    ss << "gradient_" << getOption("name") << "_" << iind << "_" << oind;
    ret.setOption("name", ss.str());

    // Same input scheme
    ret.setInputScheme(input_.scheme);

    // Output names
    std::vector<std::string> ionames;
    ionames.reserve(ret.getNumOutputs());
    ionames.push_back("grad");
    for (int i=0;i<getNumOutputs();++i) {
      ionames.push_back(output_.scheme.entryLabel(i));
    }

    ret.setOutputScheme(ionames);

    return ret;
  }

  Function FunctionInternal::tangent(int iind, int oind) {
    // Assert scalar
    casadi_assert_message(input(iind).isScalar(),
                          "Only tangent of scalar input functions allowed. Use jacobian instead.");

    // Generate gradient function
    Function ret = getTangent(iind, oind);

    // Give it a suitable name
    stringstream ss;
    ss << "tangent_" << getOption("name") << "_" << iind << "_" << oind;
    ret.setOption("name", ss.str());

    // Same input scheme
    ret.setInputScheme(input_.scheme);

    // Output names
    std::vector<std::string> ionames;
    ionames.reserve(ret.getNumOutputs());
    ionames.push_back("tangent");
    for (int i=0;i<getNumOutputs();++i) {
      ionames.push_back(output_.scheme.entryLabel(i));
    }

    ret.setOutputScheme(ionames);

    return ret;
  }

  Function FunctionInternal::hessian(int iind, int oind) {
    log("FunctionInternal::hessian");

    // Assert scalar
    casadi_assert_message(output(oind).isScalar(), "Only hessians of scalar functions allowed.");

    // Generate gradient function
    Function ret = getHessian(iind, oind);

    // Give it a suitable name
    stringstream ss;
    ss << "hessian_" << getOption("name") << "_" << iind << "_" << oind;
    ret.setOption("name", ss.str());

    // Same input scheme
    ret.setInputScheme(input_.scheme);

    // Output names
    std::vector<std::string> ionames;
    ionames.reserve(ret.getNumOutputs());
    ionames.push_back("hess");
    ionames.push_back("grad");
    for (int i=0;i<getNumOutputs();++i) {
      ionames.push_back(output_.scheme.entryLabel(i));
    }

    ret.setOutputScheme(ionames);

    return ret;
  }

  Function FunctionInternal::getGradient(int iind, int oind) {
    Function f = wrapMXFunction();
    f.init();
    return f.gradient(iind, oind);
  }

  Function FunctionInternal::getTangent(int iind, int oind) {
    Function f = wrapMXFunction();
    f.init();
    return f.tangent(iind, oind);
  }

  MXFunction FunctionInternal::wrapMXFunction() {
    vector<MX> arg = symbolicInput();
    vector<MX> res = shared_from_this<Function>()(arg);

    MXFunction f = MXFunction(arg, res);
    f.setOption("name", "wrap_" + string(getOption("name")));
    f.setInputScheme(getInputScheme());
    f.setOutputScheme(getOutputScheme());
    f.setOption("ad_weight", adWeight());
    f.setOption("ad_weight_sp", adWeightSp());
    for (int i=0; i<3; ++i) {
      string n;
      switch (i) {
      case 0: n="custom_forward"; break;
      case 1: n="custom_reverse"; break;
      case 2: n="full_jacobian"; break;
      }
      if (hasSetOption(n)) f.setOption(n, getOption(n));
    }
    return f;
  }

  Function FunctionInternal::getHessian(int iind, int oind) {
    log("FunctionInternal::getHessian");

    // Create gradient function
    log("FunctionInternal::getHessian generating gradient");
    Function g = gradient(iind, oind);
    g.setOption("verbose", getOption("verbose"));
    g.setInputScheme(input_.scheme);
    g.init();

    // Return the Jacobian of the gradient, exploiting symmetry (the gradient has output index 0)
    log("FunctionInternal::getHessian generating Jacobian of gradient");
    return g.jacobian(iind, 0, false, true);
  }

  void FunctionInternal::log(const string& msg) const {
    if (verbose()) {
      cout << "CasADi log message: " << msg << endl;
    }
  }

  void FunctionInternal::log(const string& fcn, const string& msg) const {
    if (verbose()) {
      cout << "CasADi log message: In \"" << fcn << "\" --- " << msg << endl;
    }
  }

  bool FunctionInternal::verbose() const {
    return verbose_;
  }

  bool FunctionInternal::monitored(const string& mod) const {
    return monitors_.count(mod)>0;
  }

  const Dictionary & FunctionInternal::getStats() const {
    return stats_;
  }

  GenericType FunctionInternal::getStat(const string & name) const {
    // Locate the statistic
    Dictionary::const_iterator it = stats_.find(name);

    // Check if found
    if (it == stats_.end()) {
      casadi_error("Statistic: " << name << " has not been set." << endl <<
                   "Note: statistcs are only set after an evaluate call");
    }

    return GenericType(it->second);
  }

  std::vector<MX> FunctionInternal::symbolicInput() const {
    assertInit();
    vector<MX> ret(getNumInputs());
    for (int i=0; i<ret.size(); ++i) {
      stringstream name;
      name << "x_" << i;
      ret[i] = MX::sym(name.str(), input(i).sparsity());
    }
    return ret;
  }

  std::vector<MX> FunctionInternal::symbolicOutput() const {
    assertInit();
    vector<MX> ret(getNumOutputs());
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
    vector<SX> ret(getNumInputs());
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
    for (int ind=0; ind<getNumInputs(); ++ind) {
      vector<double> &v = inputNoCheck(ind).data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Clear the adjoint seeds/forward sensitivities
    for (int ind=0; ind<getNumOutputs(); ++ind) {
      vector<double> &v = outputNoCheck(ind).data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Get seeds and sensitivities
    bvec_t* input_v = get_bvec_t(inputNoCheck(iind).data());
    bvec_t* output_v = get_bvec_t(outputNoCheck(oind).data());
    bvec_t* seed_v = use_fwd ? input_v : output_v;
    bvec_t* sens_v = use_fwd ? output_v : input_v;

    // Number of sweeps needed
    int nsweep = use_fwd ? nsweep_fwd : nsweep_adj;

    // The number of zeros in the seed and sensitivity directions
    int nz_seed = use_fwd ? nz_in  : nz_out;
    int nz_sens = use_fwd ? nz_out : nz_in;

    // Print
    if (verbose()) {
      std::cout << "FunctionInternal::getJacSparsity: using "
                << (use_fwd ? "forward" : "adjoint") << " mode: ";
      std::cout << nsweep << " sweeps needed for " << nz_seed << " directions" << endl;
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
          std::cout << progress << " %"  << endl;
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
    for (int ind=0; ind<getNumInputs(); ++ind) input(ind).setZero();
    for (int ind=0; ind<getNumOutputs(); ++ind) output(ind).setZero();

    // Construct sparsity pattern
    Sparsity ret = Sparsity::triplet(nz_out, nz_in, use_fwd ? jcol : jrow, use_fwd ? jrow : jcol);

    casadi_log("Formed Jacobian sparsity pattern (dimension " << ret.shape() << ", "
               << ret.nnz() << " nonzeros, " << (100.0*ret.nnz())/ret.numel() << " % nonzeros).");
    casadi_log("FunctionInternal::getJacSparsity end ");

    // Return sparsity pattern
    return ret;
  }

  Sparsity FunctionInternal::getJacSparsityHierarchicalSymm(int iind, int oind) {
    casadi_assert(spCanEvaluate(true));

    // Number of nonzero inputs
    int nz = input(iind).nnz();

    // Clear the forward seeds/adjoint sensitivities
    for (int ind=0; ind<getNumInputs(); ++ind) {
      vector<double> &v = inputNoCheck(ind).data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Clear the adjoint seeds/forward sensitivities
    for (int ind=0; ind<getNumOutputs(); ++ind) {
      vector<double> &v = outputNoCheck(ind).data();
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
      casadi_log("Block size: " << granularity);

      // Clear the sparsity triplet acccumulator
      jcol.clear();
      jrow.clear();

      // Clear the fine block structure
      fine.clear();

      Sparsity D = r.starColoring();

      casadi_log("Star coloring on " << r.dimString() << ": " << D.size2() << " <-> " << D.size1());

      // Reset the virtual machine
      spInit(true);

      // Get seeds and sensitivities
      bvec_t* input_v = get_bvec_t(inputNoCheck(iind).data());
      bvec_t* output_v = get_bvec_t(outputNoCheck(oind).data());
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
            duplicates.makeSparse();
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
            for (int ind=0; ind<getNumInputs(); ++ind) {
              vector<double> &v = inputNoCheck(ind).data();
              if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
            }

            // Clear the adjoint seeds/forward sensitivities, ready for next bvec sweep
            for (int ind=0; ind<getNumOutputs(); ++ind) {
              vector<double> &v = outputNoCheck(ind).data();
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
      coarse = fine;
      hasrun = true;
    }

    casadi_log("Number of sweeps: " << nsweeps);
    casadi_log("Formed Jacobian sparsity pattern (dimension " << r.shape() <<
               ", " << r.nnz() << " nonzeros, " << (100.0*r.nnz())/r.numel() << " % nonzeros).");

    return r.T();
  }

  Sparsity FunctionInternal::getJacSparsityHierarchical(int iind, int oind) {
    // Number of nonzero inputs
    int nz_in = input(iind).nnz();

    // Number of nonzero outputs
    int nz_out = output(oind).nnz();

    // Clear the forward seeds/adjoint sensitivities
    for (int ind=0; ind<getNumInputs(); ++ind) {
      vector<double> &v = inputNoCheck(ind).data();
      if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
    }

    // Clear the adjoint seeds/forward sensitivities
    for (int ind=0; ind<getNumOutputs(); ++ind) {
      vector<double> &v = outputNoCheck(ind).data();
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
      casadi_log("Block size: " << granularity_col << " x " << granularity_row);

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

      casadi_log("Coloring on " << r.dimString() << " (fwd seeps: " << D1.size2() <<
                 " , adj sweeps: " << D2.size1() << ")");

      // Use forward mode?
      int fwd_cost = use_fwd ? granularity_row: granularity_col;
      int adj_cost = use_fwd ? granularity_col: granularity_row;

      // Use whatever required less colors if we tried both (with preference to forward mode)
      if ((w*D1.size2()*fwd_cost <= (1-w)*D2.size2()*adj_cost)) {
        use_fwd = true;
        casadi_log("Forward mode chosen (fwd cost: " << w*D1.size2()*fwd_cost << ", adj cost: "
                   << (1-w)*D2.size2()*adj_cost << ")");
      } else {
        use_fwd = false;
        casadi_log("Adjoint mode chosen (adj cost: " << w*D1.size2()*fwd_cost << ", adj cost: "
                   << (1-w)*D2.size2()*adj_cost << ")");
      }

      // Reset the virtual machine
      spInit(use_fwd);

      // Get seeds and sensitivities
      bvec_t* input_v = get_bvec_t(inputNoCheck(iind).data());
      bvec_t* output_v = get_bvec_t(outputNoCheck(oind).data());
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
            for (int ind=0; ind<getNumInputs(); ++ind) {
              vector<double> &v = inputNoCheck(ind).data();
              if (!v.empty()) fill_n(get_bvec_t(v), v.size(), bvec_t(0));
            }

            // Clear the adjoint seeds/forward sensitivities, ready for next bvec sweep
            for (int ind=0; ind<getNumOutputs(); ++ind) {
              vector<double> &v = outputNoCheck(ind).data();
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
    casadi_log("Number of sweeps: " << nsweeps);
    casadi_log("Formed Jacobian sparsity pattern (dimension " << r.shape() <<
               ", " << r.nnz() << " nonzeros, " << (100.0*r.nnz())/r.numel() << " % nonzeros).");

    return r.T();
  }

  Sparsity FunctionInternal::getJacSparsity(int iind, int oind, bool symmetric) {
    // Check if we are able to propagate dependencies through the function
    if (spCanEvaluate(true) || spCanEvaluate(false)) {

      if (input(iind).nnz()>3*bvec_size && output(oind).nnz()>3*bvec_size) {
        if (symmetric) {
          return getJacSparsityHierarchicalSymm(iind, oind);
        } else {
          return getJacSparsityHierarchical(iind, oind);
        }
      } else {
        return getJacSparsityPlain(iind, oind);
      }


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
      casadi_assert(hasDerForward());

      // Star coloring if symmetric
      log("FunctionInternal::getPartition starColoring");
      D1 = A.starColoring();
      casadi_log("Star coloring completed: " << D1.size2() << " directional derivatives needed ("
                 << A.size1() << " without coloring).");

    } else {
      casadi_assert(hasDerForward() || hasDerReverse());
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
            if (verbose()) cout << "Forward mode coloring interrupted (more than "
                               << max_colorings_to_test << " needed)." << endl;
          } else {
            if (verbose()) cout << "Forward mode coloring completed: "
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
            if (verbose()) cout << "Adjoint mode coloring interrupted (more than "
                               << max_colorings_to_test << " needed)." << endl;
          } else {
            if (verbose()) cout << "Adjoint mode coloring completed: "
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
    size_t ni, nr;
    nTmp(ni, nr);
    itmp_.resize(ni);
    rtmp_.resize(nr);

    // Get pointers to input arguments
    vector<cp_double> arg(getNumInputs());
    for (int i=0; i<arg.size(); ++i) arg[i]=input(i).ptr();

    // Get pointers to output arguments
    vector<p_double> res(getNumOutputs());
    for (int i=0; i<res.size(); ++i) res[i]=output(i).ptr();

    // Call memory-less
    evalD(getPtr(arg), getPtr(res), getPtr(itmp_), getPtr(rtmp_));
  }

  void FunctionInternal::evalD(cp_double* arg,
                               p_double* res, int* itmp, double* rtmp) {
    // Number of inputs and outputs
    int num_in = getNumInputs();
    int num_out = getNumOutputs();

    // Pass the inputs to the function
    for (int i=0; i<num_in; ++i) {
      if (arg[i] != 0) {
        setInput(arg[i], i);
      } else {
        setInput(0., i);
      }
    }

    // Evaluate
    evaluate();

    // Get the outputs
    for (int i=0; i<num_out; ++i) {
      if (res[i] != 0) getOutput(res[i], i);
    }
  }

  void FunctionInternal::evalSX(cp_SXElement* arg, p_SXElement* res,
                                int* itmp, SXElement* rtmp) {
    // Number of inputs and outputs
    int num_in = getNumInputs();
    int num_out = getNumOutputs();

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
    // Get the number of inputs and outputs
    int num_in = getNumInputs();
    int num_out = getNumOutputs();

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
          arg2[i] = arg2[i].setSparse(input(i).sparsity());
      return evalSX(arg2, res);
    }

    // Allocate results
    res.resize(num_out);
    for (int i=0; i<num_out; ++i)
      if (res[i].sparsity()!=output(i).sparsity())
        res[i] = SX::zeros(output(i).sparsity());

    // Allocate temporary memory if needed
    size_t ni, nr;
    nTmp(ni, nr);
    itmp_.resize(ni);
    vector<SXElement> rtmp(nr);

    // Get pointers to input arguments
    vector<cp_SXElement> argp(arg.size());
    for (int i=0; i<argp.size(); ++i) argp[i]=getPtr(arg[i]);

    // Get pointers to output arguments
    vector<p_SXElement> resp(getNumOutputs());
    for (int i=0; i<resp.size(); ++i) resp[i]=getPtr(res[i]);

    // Call memory-less
    evalSX(getPtr(argp), getPtr(resp), getPtr(itmp_), getPtr(rtmp));
  }

  void FunctionInternal::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    MXFunction f = wrapMXFunction();
    f.init();
    f.call(arg, res, true);
  }

  void FunctionInternal::spEvaluate(bool fwd) {
    // Allocate temporary memory if needed
    size_t ni, nr;
    nTmp(ni, nr);
    itmp_.resize(ni);
    rtmp_.resize(nr);
    int *iw = getPtr(itmp_);
    bvec_t *w = reinterpret_cast<bvec_t*>(getPtr(rtmp_));

    // Get pointers to output arguments
    vector<p_bvec_t> res(getNumOutputs());
    for (int i=0; i<res.size(); ++i) res[i]=reinterpret_cast<bvec_t*>(output(i).ptr());

    if (fwd) {
      // Get pointers to input arguments
      vector<cp_bvec_t> arg(getNumInputs());
      for (int i=0; i<arg.size(); ++i) arg[i]=reinterpret_cast<const bvec_t*>(input(i).ptr());

      // Call memory-less
      spFwd(getPtr(arg), getPtr(res), iw, w);
    } else {
      // Get pointers to input arguments
      vector<p_bvec_t> arg(getNumInputs());
      for (int i=0; i<arg.size(); ++i) arg[i]=reinterpret_cast<bvec_t*>(input(i).ptr());

      // Call memory-less
      spAdj(getPtr(arg), getPtr(res), iw, w);
    }
  }

  void FunctionInternal::spEvaluateViaJacSparsity(bool fwd) {
    if (fwd) {
      // Clear the outputs
      for (int oind = 0; oind < getNumOutputs(); ++oind) {
        // Get data array for output and clear it
        bvec_t *outputd = get_bvec_t(output(oind).data());
        fill_n(outputd, output(oind).nnz(), 0);
      }
    }

    // Loop over inputs
    for (int iind = 0; iind < getNumInputs(); ++iind) {
      // Skip if no seeds
      if (fwd && input(iind).isEmpty())
        continue;

      // Get data array for input
      bvec_t *inputd = get_bvec_t(input(iind).data());

      // Loop over outputs
      for (int oind = 0; oind < getNumOutputs(); ++oind) {

        // Skip if no seeds
        if (!fwd && output(oind).isEmpty())
          continue;

        // Get the sparsity of the Jacobian block
        Sparsity sp = jacSparsity(iind, oind, true, false);
        if (sp.isNull() || sp.nnz() == 0)
          continue; // Skip if zero

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
      for (int oind=0; oind < getNumOutputs(); ++oind) {
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
      // Generate a Jacobian
      Function ret = getJacobian(iind, oind, compact, symmetric);

      // Give it a suitable name
      stringstream ss;
      ss << "jacobian_" << getOption("name") << "_" << iind << "_" << oind;
      ret.setOption("name", ss.str());
      ret.setOption("verbose", getOption("verbose"));

      // Same input scheme
      ret.setInputScheme(input_.scheme);

      // Output names
      std::vector<std::string> ionames;
      ionames.reserve(ret.getNumOutputs());
      ionames.push_back("jac");
      for (int i=0; i<getNumOutputs(); ++i) {
        ionames.push_back(output_.scheme.entryLabel(i));
      }
      ret.setOutputScheme(ionames);

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

  Function FunctionInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
    return getNumericJacobian(iind, oind, compact, symmetric);
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

    // Return value
    Function ret;
    if (hasSetOption("custom_forward")) {
      /// User-provided derivative generator function
      DerivativeGenerator dergen = getOption("custom_forward");
      Function this_ = shared_from_this<Function>();
      ret = dergen(this_, nfwd, user_data_);
      // Fails for ImplicitFunction
    } else {
      casadi_assert(hasDerForward());
      ret = getDerForward(nfwd);
    }

    // Give it a suitable name
    stringstream ss;
    ss << "derForward_" << getOption("name") << "_" << nfwd;
    ret.setOption("name", ss.str());

    // Get the number of inputs and outputs
    int n_in = getNumInputs();
    int n_out = getNumOutputs();

    // Names of inputs
    std::vector<std::string> i_names;
    i_names.reserve(n_in + n_out + n_in*nfwd);

    // Nondifferentiated inputs
    for (int i=0; i<n_in; ++i) {
      i_names.push_back("der_" + input_.scheme.entryLabel(i));
    }

    // Nondifferentiated outputs (given)
    for (int i=0; i<n_out; ++i) {
      i_names.push_back("der_" + output_.scheme.entryLabel(i));
    }

    // Forward seeds
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<n_in; ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << input_.scheme.entryLabel(i);
        i_names.push_back(ss.str());
      }
    }

    // Pass to return object
    ret.setInputScheme(i_names);

    // Names of outputs
    std::vector<std::string> o_names;
    o_names.reserve(n_out*nfwd);

    // Forward sensitivities
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<n_out; ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << output_.scheme.entryLabel(i);
        o_names.push_back(ss.str());
      }
    }

    // Pass to return object
    ret.setOutputScheme(o_names);

    // Initialize it
    ret.init();

    // Consistency check for inputs
    for (int i=0; i<ret.getNumInputs(); ++i) {
      const Sparsity& sp = i<n_in ? input(i).sparsity() :
        i<n_in+n_out ? output(i-n_in).sparsity() :
        input((i-n_in-n_out) % n_in).sparsity();
      casadi_assert_message(ret.input(i).shape()==sp.shape(),
                            "Incorrect shape for " << ret << " input " << i << " \""
                            << i_names.at(i) << "\". Expected " << sp.shape()
                            << " but got " << ret.input(i).shape());
    }

    // Consistency check for outputs
    for (int i=0; i<ret.getNumOutputs(); ++i) {
      const Sparsity& sp = output(i % n_out).sparsity();
      casadi_assert_message(ret.output(i).shape()==sp.shape(),
                            "Incorrect shape for " << ret << " output " << i << " \""
                            << o_names.at(i) << "\". Expected " << sp.shape()
                            << " but got " << ret.output(i).shape());
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

    // Return value
    Function ret;
    if (hasSetOption("custom_reverse")) {
      /// User-provided derivative generator function
      DerivativeGenerator dergen = getOption("custom_reverse");
      Function this_ = shared_from_this<Function>();
      ret = dergen(this_, nadj, user_data_);
    } else {
      casadi_assert(hasDerReverse());
      ret = getDerReverse(nadj);
    }

    // Give it a suitable name
    stringstream ss;
    ss << "derReverse_" << getOption("name") << "_" << nadj;
    ret.setOption("name", ss.str());

    // Get the number of inputs and outputs
    int n_in = getNumInputs();
    int n_out = getNumOutputs();

    // Names of inputs
    std::vector<std::string> i_names;
    i_names.reserve(n_in + n_out + n_out*nadj);

    // Nondifferentiated inputs
    for (int i=0; i<n_in; ++i) {
      i_names.push_back("der_" + input_.scheme.entryLabel(i));
    }

    // Nondifferentiated outputs (given)
    for (int i=0; i<n_out; ++i) {
      i_names.push_back("der_" + output_.scheme.entryLabel(i));
    }

    // Adjoint seeds
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_out; ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << output_.scheme.entryLabel(i);
        i_names.push_back(ss.str());
      }
    }

    // Pass to return object
    ret.setInputScheme(i_names);

    // Names of outputs
    std::vector<std::string> o_names;
    o_names.reserve(n_in*nadj);

    // Adjoint sensitivities
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_in; ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << input_.scheme.entryLabel(i);
        o_names.push_back(ss.str());
      }
    }

    // Pass to return object
    ret.setOutputScheme(o_names);

    // Initialize it
    ret.init();

    // Consistency check for inputs
    for (int i=0; i<ret.getNumInputs(); ++i) {
      const Sparsity& sp = i<n_in ? input(i).sparsity() :
        i<n_in+n_out ? output(i-n_in).sparsity() :
        output((i-n_in-n_out) % n_out).sparsity();
      casadi_assert_message(ret.input(i).shape()==sp.shape(),
                            "Incorrect shape for " << ret << " input " << i << " \""
                            << i_names.at(i) << "\". Expected " << sp.shape()
                            << " but got " << ret.input(i).shape());
    }

    // Consistency check for outputs
    for (int i=0; i<ret.getNumOutputs(); ++i) {
      const Sparsity& sp = input(i % n_in).sparsity();
      casadi_assert_message(ret.output(i).shape()==sp.shape(),
                            "Incorrect shape for " << ret << " output " << i << " \""
                            << o_names.at(i) << "\". Expected " << sp.shape()
                            << " but got " << ret.output(i).shape());
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

  Function FunctionInternal::getDerForward(int nfwd) {
    // TODO(@jaeandersson): Fallback on finite differences
    casadi_error("FunctionInternal::getDerForward not defined for class "
                 << typeid(*this).name());
  }

  Function FunctionInternal::getDerReverse(int nadj) {
    casadi_error("FunctionInternal::getDerReverse not defined for class "
                 << typeid(*this).name());
  }

  int FunctionInternal::getNumInputNonzeros() const {
    int ret=0;
    for (int iind=0; iind<getNumInputs(); ++iind) {
      ret += input(iind).nnz();
    }
    return ret;
  }

  int FunctionInternal::getNumOutputNonzeros() const {
    int ret=0;
    for (int oind=0; oind<getNumOutputs(); ++oind) {
      ret += output(oind).nnz();
    }
    return ret;
  }

  int FunctionInternal::getNumInputElements() const {
    int ret=0;
    for (int iind=0; iind<getNumInputs(); ++iind) {
      ret += input(iind).numel();
    }
    return ret;
  }

  int FunctionInternal::getNumOutputElements() const {
    int ret=0;
    for (int oind=0; oind<getNumOutputs(); ++oind) {
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
        arg.size()<=getNumInputs(), "Function::call: number of passed-in dependencies ("
        << arg.size() << ") should not exceed the number of inputs of the function ("
        << getNumInputs() << ").");

      // Assumes initialized
      for (int i=0; i<arg.size(); ++i) {
        if (arg[i].isNull() || arg[i].isEmpty() || input(i).isEmpty()) continue;
        casadi_assert_message(
          arg[i].size2()==input(i).size2() && arg[i].size1()==input(i).size1(),
          "Evaluation::shapes of passed-in dependencies should match shapes of inputs of function."
          << endl << input_.scheme.describeInput(i) <<  " has shape (" << input(i).size2()
          << ", " << input(i).size1() << ") while a shape (" << arg[i].size2() << ", "
          << arg[i].size1() << ") was supplied.");
      }
      res = createCall(arg);
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
    res.resize(getNumOutputs());
    for (int i=0;i<res.size();++i) {
      res[i]=output(i);
    }
  }

  Function FunctionInternal::getNumericJacobian(int iind, int oind, bool compact, bool symmetric) {
    Function f = wrapMXFunction();
    f.init();
    return f->getNumericJacobian(iind, oind, compact, symmetric);
  }

  Function FunctionInternal::fullJacobian() {
    if (full_jacobian_.alive()) {
      // Return cached Jacobian
      return shared_cast<Function>(full_jacobian_.shared());
    } else {
      Function ret;
      if (hasSetOption("full_jacobian")) {
        /// User-provided Jacobian function
        ret = getOption("full_jacobian");
      } else {
        // Generate full Jacobian
        ret = getFullJacobian();
      }

      // Give it a suitable name
      stringstream ss;
      ss << "jacobian_" << getOption("name");
      ret.setOption("name", ss.str());
      ret.setOption("verbose", getOption("verbose"));

      // Set input and output schemes
      ret.setInputScheme(input_.scheme);
      std::vector<std::string> oscheme(1, "jac");
      ret.setOutputScheme(oscheme);

      // Initialize
      ret.init();

      // Consistency check
      casadi_assert(ret.getNumInputs()==getNumInputs());
      casadi_assert(ret.getNumOutputs()==1);

      // Cache it for reuse and return
      full_jacobian_ = ret;
      return ret;
    }
  }

  Function FunctionInternal::getFullJacobian() {
    casadi_assert(hasDerForward() || hasDerReverse());

    // Number inputs and outputs
    int n_in = getNumInputs();
    int n_out = getNumOutputs();

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
      MXFunction tmp(arg, res);
      tmp.setOption("ad_weight", adWeight());
      tmp.init();
      J = tmp.jac();
    }

    // Make sure argv is the input of J
    if (n_in!=1) {
      J = substitute(J, arg, veccat(ret_argv));
    }

    // Form an expression for the full Jacobian
    return MXFunction(ret_argv, J);
  }

  void FunctionInternal::generateFunction(CodeGenerator& gen,
                                          const std::string& fname) const {
    assertInit();

    // Add standard math
    gen.addInclude("math.h");

    // Add auxiliaries. TODO: Only add the auxiliaries that are actually used
    gen.addAuxiliary(CodeGenerator::AUX_SQ);
    gen.addAuxiliary(CodeGenerator::AUX_SIGN);

    // Generate declarations
    generateDeclarations(gen.functions, gen.real_t, gen);

    // Define function
    gen.functions
      << "/* " << getSanitizedName() << " */" << endl
      << "int " << fname << "(const " << gen.real_t << "* const* arg, " << gen.real_t
      << "* const* res, int* iii, " << gen.real_t << "* w) {" << endl;

    // Insert the function body
    generateBody(gen.functions, gen.real_t, gen);

    // Finalize the function
    gen.functions << "  return 0;" << endl
                  << "}" << endl
                  << endl;
  }

  void FunctionInternal::generateMeta(CodeGenerator& gen, const std::string& fname) const {
    // Short-hands
    int n_in = getNumInputs();
    int n_out = getNumOutputs();
    stringstream &s = gen.functions;

    // Function that returns the number of inputs and outputs
    s << "int " << fname << "_narg(int *n_in, int *n_out) {" << endl
      << "  *n_in = " << n_in << ";" << endl
      << "  *n_out = " << n_out << ";" << endl
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
      int ind = gen.addSparsity(sp);
      num_io_patterns = std::max(num_io_patterns, ind+1);

      // Store the index
      io_sparsity_index.insert(std::pair<int, int>(ind, i));
    }

    // Function that returns the sparsity pattern
    s << "int " << fname << "_sparsity"
      << "(int i, int *nrow, int *ncol, int **colind, int **row) {" << endl;

    // Get the sparsity index using a switch
    s << "  int* s;" << endl;
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
    size_t ni, nr;
    nTmp(ni, nr);
    s << "int " << fname << "_work(int *ni, int *nr) {" << endl;
    s << "  if (ni) *ni = " << ni << ";" << endl;
    s << "  if (nr) *nr = " << nr << ";" << endl;
    s << "  return 0;" << endl;
    s << "}" << endl;
    s << endl;

    // Generate mex gateway for the function
    if (gen.mex) {
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

      // Work vectors, including input buffers
      size_t ni, nr;
      nTmp(ni, nr);
      int i_nnz = 0;
      for (int i=0; i<n_in; ++i) {
        const Sparsity& s = input(i).sparsity();
        i_nnz += s.nnz();
        nr = max(nr, static_cast<size_t>(s.size1())); // To be able to copy a column
        nr = max(nr, static_cast<size_t>(s.size2())); // To be able to copy a row
      }
      nr += i_nnz;
      s << "  int iw[" << ni << "];" << endl;
      s << "  d w[" << nr << "];" << endl;
      string fw = "w+" + gen.numToString(i_nnz);

      // Copy inputs to buffers
      int offset=0;
      s << "  const d* arg[" << n_in << "] = {0};" << endl;
      for (int i=0; i<n_in; ++i) {
        std::string p = "argv[" + gen.numToString(i) + "]";
        s << "  if (--argc>=0) arg[" << i << "] = "
               << gen.from_mex(p, "w", offset, input(i).sparsity(), fw) << endl;
        offset += input(i).nnz();
      }

      // Output sparsities
      s << "  d* res[" << n_out << "] = {0};" << endl;
      for (int i=0; i<n_out; ++i) {
        s << "  if (--resc>=0) resv[" << i << "] = "
               << gen.to_mex(output(i).sparsity(), "&res["+gen.numToString(i)+"]") << endl;
      }

      // Call the function
      s << "  i = " << fname << "(arg, res, iw, " << fw << ");" << endl;
      s << "  if (i) mexErrMsgIdAndTxt(\"Casadi:RuntimeError\",\"Evaluation of \\\"" << fname
             << "\\\" failed.\");" << endl;

      // Finalize mex gateway
      s << "}" << endl << endl;
    }
    if (gen.main) {
      // Declare wrapper
      s << "int main_" << fname << "(int argc, char* argv[]) {" << endl;

      // Work vectors and input and output buffers
      size_t ni, nr;
      nTmp(ni, nr);
      nr += getNumInputNonzeros() + getNumOutputNonzeros();
      s << "  int iw[" << ni << "];" << endl
             << "  d w[" << nr << "];" << endl;

      // Input buffers
      s << "  const d* arg[" << n_in << "] = {";
      int off=0;
      for (int i=0; i<n_in; ++i) {
        if (i!=0) s << ", ";
        s << gen.work(off);
        off += input(i).nnz();
      }
      s << "};" << endl;

      // Output buffers
      s << "  d* res[" << n_out << "] = {";
      for (int i=0; i<n_out; ++i) {
        if (i!=0) s << ", ";
        s << gen.work(off);
        off += output(i).nnz();
      }
      s << "};" << endl;

      // TODO(@jaeandersson): Read inputs from file. For now; read from stdin
      s << "  int j;" << endl
             << "  d* a = w;" << endl
             << "  for (j=0; j<" << getNumInputNonzeros() << "; ++j) "
             << "scanf(\"%lf\", a++);" << endl;

      // Call the function
      s << "  int flag = eval(arg, res, iw, " << gen.work(off) << ");" << endl
             << "  if (flag) return flag;" << endl;

      // TODO(@jaeandersson): Write outputs to file. For now: print to stdout
      s << "  const d* r = w+" << getNumInputNonzeros() << ";" << endl
             << "  for (j=0; j<" << getNumOutputNonzeros() << "; ++j) "
             << gen.printf("%g ", "*r++") << endl;
      // End with newline
      s << "  " << gen.printf("\\n") << endl;

      // Finalize function
      s << "  return 0;" << endl
        << "}" << endl << endl;
    }
  }

  void FunctionInternal::generateDeclarations(std::ostream &stream, const std::string& type,
                                              CodeGenerator& gen) const {
    // Nothing to declare
  }

  void FunctionInternal::generateBody(std::ostream &stream, const std::string& type,
                                      CodeGenerator& gen) const {
    casadi_error("FunctionInternal::generateBody: generateBody not defined for class "
                 << typeid(*this).name());
  }

  Function FunctionInternal::dynamicCompilation(Function f, std::string fname, std::string fdescr,
                                                std::string compiler) {
    // Check if f is initialized
    bool f_is_init = f.isInit();
    if (!f_is_init) f.init();

    // Filenames
    string cname = fname + ".c";
    string dlname = fname + ".so";

    // Codegen and compile
    CodeGenerator gen;
    gen.addFunction(f, fname);
    gen.compile(cname, dlname, compiler);

    // Load it
    ExternalFunction f_gen("./" + dlname);
    f_gen.setOption("name", fname + "_gen");

    // Initialize it if f was initialized
    if (f_is_init) {
      f_gen.init();
      if (verbose_) {
        cout << "Dynamically loaded " << fdescr << " (" << dlname << ")" << endl;
      }
    }
    return f_gen;
  }

  std::vector<MX> FunctionInternal::createCall(const std::vector<MX> &arg) {
    return MX::createMultipleOutput(new Call(shared_from_this<Function>(), arg));
  }

  std::vector<std::vector<MX> >
  FunctionInternal::createMap(const std::vector<std::vector<MX> > &arg,
                              const std::string& parallelization) {
    int n = arg.size();
    std::vector<std::vector<MX> > ret(n);
    if (parallelization.compare("expand")==0) {
      // Bypass the Map, call the original function n times
      for (int i=0; i<n; ++i) {
        call(arg[i], ret[i], false, false);
      }
    } else {
      // Get type of parallelization
      bool omp;
      if (parallelization.compare("openmp")==0) {
        omp = true;
      } else if (parallelization.compare("serial")==0) {
        omp = false;
      } else {
        casadi_error("Unsupported parallelization \"" << parallelization
                     << "\": Available options are expand|serial|openmp");
      }

      // Call the map
      std::vector<MX> v;
      if (omp) {
        v = MX::createMultipleOutput(new OmpMap(shared_from_this<Function>(), arg));
      } else {
        v = MX::createMultipleOutput(new Map(shared_from_this<Function>(), arg));
      }

      // Collect outputs
      std::vector<MX>::const_iterator v_it = v.begin();
      int n_out = getNumOutputs();
      for (int i=0; i<n; ++i) {
        ret[i] = std::vector<MX>(v_it, v_it+n_out);
        v_it += n_out;
      }
      casadi_assert(v_it==v.end());
    }
    return ret;
  }

  void FunctionInternal::spFwdSwitch(cp_bvec_t* arg, p_bvec_t* res,
                                     int* itmp, bvec_t* rtmp) {
    // TODO(@jaeandersson) Calculate from full-Jacobian sparsity  when necessary or more efficient
    spFwd(arg, res, itmp, rtmp);
  }

  void FunctionInternal::spFwd(cp_bvec_t* arg, p_bvec_t* res,
                               int* itmp, bvec_t* rtmp) {
    // Number inputs and outputs
    int n_in = getNumInputs();
    int n_out = getNumOutputs();

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

  void FunctionInternal::spAdjSwitch(p_bvec_t* arg, p_bvec_t* res,
                                     int* itmp, bvec_t* rtmp) {
    // TODO(@jaeandersson) Calculate from full-Jacobian sparsity  when necessary or more efficient
    spAdj(arg, res, itmp, rtmp);
  }

  void FunctionInternal::spAdj(p_bvec_t* arg, p_bvec_t* res,
                               int* itmp, bvec_t* rtmp) {
    // Number inputs and outputs
    int n_in = getNumInputs();
    int n_out = getNumOutputs();

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

  void FunctionInternal::generate(std::ostream &stream, const std::vector<int>& arg,
                                  const std::vector<int>& res, CodeGenerator& gen) const {
    // Number inputs and outputs
    int n_in = getNumInputs();
    int n_out = getNumOutputs();

    // Put in a separate scope to avoid name collisions
    stream << "  {" << endl;

    // Collect input arguments
    stream << "    const d* arg1[] = {";
    for (int i=0; i<n_in; ++i) {
      if (i!=0) stream << ", ";
      stream << gen.work(arg.at(i));
    }
    stream << "};" << endl;

    // Collect output arguments
    stream << "    d* res1[] = {";
    for (int i=0; i<n_out; ++i) {
      if (i!=0) stream << ", ";
      stream << gen.work(res.at(i));
    }
    stream << "};" << endl;

    // Get the index of the function
    int f = gen.getDependency(shared_from_this<Function>());

    // Call function
    stream << "    i=f" << f << "(arg1, res1, iii, w);" << endl;
    stream << "    if (i) return i;" << endl;

    // Finalize the function call
    stream << "  }" << endl;
  }

  void FunctionInternal::nTmp(size_t& ni, size_t& nr) const {
    ni=itmp_.size();
    nr=rtmp_.size();
  }

  void FunctionInternal::printPart(const MXNode* node, std::ostream &stream, int part) const {
    if (part == 0) {
      repr(stream);
      stream << ".call([";
    } else if (part == getNumInputs()) {
      stream << "])";
    } else {
      stream << ", ";
    }
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
      casadi_assert(hasSetOption("name"));
      string name = getOption("name");
      std::replace_if(name.begin(), name.end(), isBadChar, '_');
      return name;
  }

  bool FunctionInternal::hasDerivative() const {
    return hasDerForward() || hasDerReverse() ||
      full_jacobian_.alive() || hasSetOption("full_jacobian");
  }

  bool FunctionInternal::fwdViaJac(int nfwd) {
    if (!hasDerForward()) return true;

    // Jacobian calculation penalty factor
    const int jac_penalty = 2;

    // Heuristic 1: Jac calculated via forward mode likely cheaper
    if (jac_penalty*getNumInputNonzeros()<nfwd) return true;

    // Heuristic 2: Jac calculated via reverse mode likely cheaper
    double w = adWeight();
    if (hasDerReverse() && jac_penalty*(1-w)*getNumOutputNonzeros()<w*nfwd)
      return true;

    return false;
  }

  bool FunctionInternal::adjViaJac(int nadj) {
    if (!hasDerReverse()) return true;

    // Jacobian calculation penalty factor
    const int jac_penalty = 2;

    // Heuristic 1: Jac calculated via reverse mode likely cheaper
    if (jac_penalty*getNumOutputNonzeros()<nadj) return true;

    // Heuristic 2: Jac calculated via forward mode likely cheaper
    double w = adWeight();
    if (hasDerForward() && jac_penalty*w*getNumInputNonzeros()<(1-w)*nadj)
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
    int n_in = getNumInputs();
    int n_out = getNumOutputs();

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
          fsens[d][i] = reshape(fsens[d][i], output(i).shape());
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
    vector<MX> x = MX::createMultipleOutput(new Call(dfcn, darg));
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
    int n_in = getNumInputs();
    int n_out = getNumOutputs();

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
          if (asens[d][i].isEmpty(true)) {
            asens[d][i] = reshape(a[i], input(i).shape());
          } else {
            asens[d][i] += reshape(a[i], input(i).shape());
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
    vector<MX> x = MX::createMultipleOutput(new Call(dfcn, darg));
    vector<MX>::iterator x_it = x.begin();

    // Retrieve sensitivities
    for (int d=0; d<nadj; ++d) {
      asens[d].resize(n_in);
      for (int i=0; i<n_in; ++i) {
        if (asens[d][i].isEmpty(true)) {
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
    int n_in = getNumInputs();
    int di=0;
    casadi_assert(arg.size()==n_in);
    for (int i=0; i<n_in; ++i) {
      dfcn.setInput(arg[i], di++);
    }

    // Pass outputs
    int n_out = getNumOutputs();
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
    casadi_assert(di==dfcn.getNumInputs()); // Consistency check

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
    casadi_assert(di==dfcn.getNumOutputs()); // Consistency check
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
    int n_in = getNumInputs();
    int di=0;
    casadi_assert(arg.size()==n_in);
    for (int i=0; i<n_in; ++i) {
      dfcn.setInput(arg[i], di++);
    }

    // Pass outputs
    int n_out = getNumOutputs();
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
    casadi_assert(di==dfcn.getNumInputs()); // Consistency check

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
    casadi_assert(di==dfcn.getNumOutputs()); // Consistency check
  }

  double FunctionInternal::adWeight() {
    // If reverse mode derivatives unavailable, use forward
    if (!hasDerReverse()) return 0;

    // If forward mode derivatives unavailable, use reverse
    if (!hasDerForward()) return 1;

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

} // namespace casadi
