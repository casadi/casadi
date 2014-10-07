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
#include "../mx/call_function.hpp"
#include <typeinfo>
#include "../std_vector_tools.hpp"
#include "mx_function.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../matrix/sparsity_tools.hpp"
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
    addOption("ad_mode",                  OT_STRING,              "automatic",
              "How to calculate the Jacobians.",
              "forward: only forward mode|reverse: only adjoint mode|automatic: "
              "a heuristic decides which is more appropriate");
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
    addOption("derivative_generator",     OT_DERIVATIVEGENERATOR,   GenericType(),
              "Function that returns a derivative function given a number of forward "
              "and reverse directional derivative, overrides internal routines. "
              "Check documentation of DerivativeGenerator.");

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
    for (vector<vector<WeakRef> >::iterator i=derivative_fcn_.begin();
        i!=derivative_fcn_.end(); ++i) {
      for (vector<WeakRef>::iterator j=i->begin(); j!=i->end(); ++j) {
        if (!j->isNull()) {
          *j = getcopy(j->shared(), already_copied);
        }
      }
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
        SparseStorage<Sparsity>(Sparsity::sparse(getNumOutputs(), getNumInputs()));
    jac_ = jac_compact_ = SparseStorage<WeakRef>(Sparsity::sparse(getNumOutputs(), getNumInputs()));

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
        stream << " Inputs (" << getNumInputs() << "):" << std::endl;
        for (int i=0;i<getNumInputs();i++) {
          stream << "  " << i << ". " << input(i).dimString() << std::endl;
        }
      } else {
        stream << " Inputs (" << input_.scheme.name() << ": " << getNumInputs()
               << "):" << std::endl;
        for (int i=0;i<getNumInputs();i++) {
          stream << "  " << i  << ". (" << input_.scheme.describe(i) << ")   "
                 << input(i).dimString() << std::endl;
        }
      }
    }
    if (getNumOutputs()==1) {
      stream << " Output: " << output().dimString() << endl;
    } else {
      if (output_.scheme.isNull()) {
        stream << " Outputs (" << getNumOutputs() << "):" << std::endl;
        for (int i=0;i<getNumOutputs();i++) {
          stream << "  " << i << ". " << output(i).dimString() << std::endl;
        }
      } else {
        stream << " Outputs (" << output_.scheme.name() << ": "
               << getNumOutputs() << "):" << std::endl;
        for (int i=0;i<getNumOutputs();i++) {
          stream << "  " << i << ". (" << output_.scheme.describe(i) << ")   "
                 << output(i).dimString() << std::endl;
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
    vector<MX> res = shared_from_this<Function>().call(arg);

    MXFunction f = MXFunction(arg, res);
    f.setOption("name", "wrap_" + string(getOption("name")));
    f.setInputScheme(getInputScheme());
    f.setOutputScheme(getOutputScheme());
    f.setOption("ad_mode", getOption("ad_mode")); // Why?
    if (hasSetOption("derivative_generator"))
        f.setOption("derivative_generator", getOption("derivative_generator"));

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
    vector<MX> ret(getNumInputs());
    assertInit();
    for (int i=0; i<ret.size(); ++i) {
      stringstream name;
      name << "x_" << i;
      ret[i] = MX::sym(name.str(), input(i).sparsity());
    }
    return ret;
  }

  std::vector<MX> FunctionInternal::symbolicOutput(const std::vector<MX>& arg) {
    return shared_from_this<Function>().call(arg);
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
    int nz_in = input(iind).size();

    // Number of nonzero outputs
    int nz_out = output(oind).size();

    // Number of forward sweeps we must make
    int nsweep_fwd = nz_in/bvec_size;
    if (nz_in%bvec_size>0) nsweep_fwd++;

    // Number of adjoint sweeps we must make
    int nsweep_adj = nz_out/bvec_size;
    if (nz_out%bvec_size>0) nsweep_adj++;

    // Use forward mode?
    bool use_fwd = spCanEvaluate(true) && nsweep_fwd <= nsweep_adj;

    // Override default behavior?
    if (getOption("ad_mode") == "forward") {
      use_fwd = true;
    } else if (getOption("ad_mode") == "reverse") {
      use_fwd = false;
    }

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
               << ret.size() << " nonzeros, " << (100.0*ret.size())/ret.numel() << " % nonzeros).");
    casadi_log("FunctionInternal::getJacSparsity end ");

    // Return sparsity pattern
    return ret;
  }

  Sparsity FunctionInternal::getJacSparsityHierarchicalSymm(int iind, int oind) {
    casadi_assert(spCanEvaluate(true));

    // Number of nonzero inputs
    int nz = input(iind).size();

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
          for (int k=D.colind()[csd]; k<D.colind()[csd+1]; ++k) {
            int cci = D.row()[k];

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
              for (int cri=r.colind()[cci];cri<r.colind()[cci+1];++cri) {
                lookup_col.push_back(r.row()[cri]);
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
            duplicates.sparsify();
            // NOTE: Not intended use of SubMatrix:
            SubMatrix<Matrix<int>, Sparsity, int> temp(lookup, duplicates.sparsity(), 0);
            temp = -bvec_size;

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
               ", " << r.size() << " nonzeros, " << (100.0*r.size())/r.numel() << " % nonzeros).");

    return r.T();
  }

  Sparsity FunctionInternal::getJacSparsityHierarchical(int iind, int oind) {
    // Number of nonzero inputs
    int nz_in = input(iind).size();

    // Number of nonzero outputs
    int nz_out = output(oind).size();

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

      // Adjoint mode penalty factor (adjoint mode is usually more expensive to calculate)
      int adj_penalty = 2;

      int fwd_cost = use_fwd ? granularity_row: granularity_col;
      int adj_cost = use_fwd ? granularity_col: granularity_row;

      // Use whatever required less colors if we tried both (with preference to forward mode)
      if ((D1.size2()*fwd_cost <= adj_penalty*D2.size2()*adj_cost)) {
        use_fwd = true;
        casadi_log("Forward mode chosen (fwd cost: " << D1.size2()*fwd_cost << ", adj cost: "
                   << adj_penalty*D2.size2()*adj_cost << ")");
      } else {
        use_fwd = false;
        casadi_log("Adjoint mode chosen (adj cost: " << D1.size2()*fwd_cost << ", adj cost: "
                   << adj_penalty*D2.size2()*adj_cost << ")");
      }

      use_fwd = spCanEvaluate(true) && use_fwd;

      // Override default behavior?
      if (getOption("ad_mode") == "forward") {
        use_fwd = true;
      } else if (getOption("ad_mode") == "reverse") {
        use_fwd = false;
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
          for (int k=D.colind()[csd]; k<D.colind()[csd+1]; ++k) {
            int cci = D.row()[k];

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
              for (int cri=rT.colind()[cci];cri<rT.colind()[cci+1];++cri) {
                lookup_col.push_back(rT.row()[cri]);
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
               ", " << r.size() << " nonzeros, " << (100.0*r.size())/r.numel() << " % nonzeros).");

    return r.T();
  }

  Sparsity FunctionInternal::getJacSparsity(int iind, int oind, bool symmetric) {
    // Check if we are able to propagate dependencies through the function
    if (spCanEvaluate(true) || spCanEvaluate(false)) {

      if (input(iind).size()>3*bvec_size && output(oind).size()>3*bvec_size) {
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
      return Sparsity::dense(output(oind).size(), input(iind).size());
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
          casadi_assert(sp.size1()==output(oind).size());

          // New row for each old row
          vector<int> row_map = output(oind).sparsity().getElements();

          // Insert rows
          sp.enlargeRows(output(oind).numel(), row_map);
        }

        // Enlarge if sparse input
        if (input(iind).numel()!=sp.size2()) {
          casadi_assert(sp.size2()==input(iind).size());

          // New column for each old column
          vector<int> col_map = input(iind).sparsity().getElements();

          // Insert columns
          sp.enlargeColumns(input(iind).numel(), col_map);
        }

        // Save
        jsp = sp;
      }
    }

    // If still null, not dependent
    if (jsp.isNull()) {
      jsp = Sparsity::sparse(output(oind).size(), input(iind).size());
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

    // Which AD mode?
    bool test_ad_fwd=true, test_ad_adj=true;
    if (getOption("ad_mode") == "forward") {
      test_ad_adj = false;
    } else if (getOption("ad_mode") == "reverse") {
      test_ad_fwd = false;
    } else if (getOption("ad_mode") != "automatic") {
      casadi_error("FunctionInternal::jac: Unknown ad_mode \"" << getOption("ad_mode")
                   << "\". Possible values are \"forward\", \"reverse\" and \"automatic\".");
    }

    // Get seed matrices by graph coloring
    if (symmetric) {

      // Star coloring if symmetric
      log("FunctionInternal::getPartition starColoring");
      D1 = A.starColoring();
      casadi_log("Star coloring completed: " << D1.size2() << " directional derivatives needed ("
                 << A.size1() << " without coloring).");

    } else {

      // Adjoint mode penalty factor (adjoint mode is usually more expensive to calculate)
      int adj_penalty = 2;

      // Best coloring encountered so far
      int best_coloring = numeric_limits<int>::max();

      // Test forward mode first?
      bool test_fwd_first = A.size1() <= adj_penalty*A.size2();
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
          D1 = AT.unidirectionalColoring(A, best_coloring);
          if (D1.isNull()) {
            if (verbose()) cout << "Forward mode coloring interrupted (more than "
                               << best_coloring << " needed)." << endl;
          } else {
            if (verbose()) cout << "Forward mode coloring completed: "
                               << D1.size2() << " directional derivatives needed ("
                               << A.size1() << " without coloring)." << endl;
            D2 = Sparsity();
            best_coloring = D1.size2();
          }
        } else {
          log("FunctionInternal::getPartition unidirectional coloring (adjoint mode)");
          int max_colorings_to_test = best_coloring/adj_penalty;
          D2 = A.unidirectionalColoring(AT, max_colorings_to_test);
          if (D2.isNull()) {
            if (verbose()) cout << "Adjoint mode coloring interrupted (more than "
                               << max_colorings_to_test << " needed)." << endl;
          } else {
            if (verbose()) cout << "Adjoint mode coloring completed: "
                               << D2.size2() << " directional derivatives needed ("
                               << A.size2() << " without coloring)." << endl;
            D1 = Sparsity();
            best_coloring = D2.size2();
          }
        }
      }

    }
    log("FunctionInternal::getPartition end");
  }

  void FunctionInternal::evalSX(const std::vector<SX>& arg, std::vector<SX>& res,
                                const std::vector<std::vector<SX> >& fseed,
                                std::vector<std::vector<SX> >& fsens,
                                const std::vector<std::vector<SX> >& aseed,
                                std::vector<std::vector<SX> >& asens) {
    // Make sure initialized
    assertInit();

    // Assert number of inputs
    casadi_assert_message(getNumInputs() == arg.size(),
                          "Wrong number of inputs. Expecting "
                          << getNumInputs() << ", got " << arg.size());

    // Assert number of forward seeds
    int nfdir = fseed.size();
    for (int dir=0; dir<nfdir; ++dir) {
      casadi_assert_message(getNumInputs() == fseed[dir].size(),
                            "Wrong number of forward seeds in direction " << dir
                            << ". Expecting " << getNumInputs() << ", got " << fseed[dir].size());
    }

    // Assert number of adjoint seeds
    int nadir = aseed.size();
    for (int dir=0; dir<nadir; ++dir) {
      casadi_assert_message(getNumOutputs() == aseed[dir].size(),
                            "Wrong number of adjoint seeds in direction " << dir
                            << ". Expecting " << getNumOutputs() << ", got " << aseed[dir].size());
    }

    // Check if input sparsity pattern match (quick if sparsity matches)
    bool sparsity_matches = true;
    for (int i=0; i<getNumInputs() && sparsity_matches; ++i) {
      sparsity_matches = arg[i].sparsity()==input(i).sparsity();
    }
    if (!sparsity_matches) {
      vector<SX> arg_new(arg.size());
      for (int i=0; i<arg.size(); ++i) {
        try {
          arg_new[i] = SX(input(i).sparsity());
          arg_new[i].set(arg[i]);
        } catch(exception& ex) {
          stringstream ss;
          ss << "SXFunctionInternal::evalSX: Failed to set " << input_.scheme.describeInput(i)
             << ": " << ex.what();
          throw CasadiException(ss.str());
        }
      }
      evalSX(arg_new, res, fseed, fsens, aseed, asens);
      return;
    }

    // Check if forward seed sparsity pattern match (quick if sparsity matches)
    for (int dir=0; dir<nfdir && sparsity_matches; ++dir) {
      for (int i=0; i<getNumInputs() && sparsity_matches; ++i) {
        sparsity_matches = fseed[dir][i].sparsity()==input(i).sparsity();
      }
    }
    if (!sparsity_matches) {
      vector<vector<SX> > fseed_new(nfdir);
      for (int dir=0; dir<nfdir; ++dir) {
        fseed_new[dir].resize(getNumInputs());
        for (int i=0; i<getNumInputs(); ++i) {
          try {
            fseed_new[dir][i] = SX(input(i).sparsity());
            fseed_new[dir][i].set(fseed[dir][i]);
          } catch(exception& ex) {
            stringstream ss;
            ss << "SXFunctionInternal::evalSX: Failed to set forward seed of "
               << input_.scheme.describeInput(i) << ", direction " << dir << ": " << ex.what();
            throw CasadiException(ss.str());
          }
        }
      }
      evalSX(arg, res, fseed_new, fsens, aseed, asens);
      return;
    }

    // Check if adjoint seed sparsity pattern match (quick if sparsity matches)
    for (int dir=0; dir<nadir && sparsity_matches; ++dir) {
      for (int i=0; i<getNumOutputs() && sparsity_matches; ++i) {
        sparsity_matches = aseed[dir][i].sparsity()==output(i).sparsity();
      }
    }
    if (!sparsity_matches) {
      vector<vector<SX> > aseed_new(nadir);
      for (int dir=0; dir<nadir; ++dir) {
        aseed_new[dir].resize(getNumOutputs());
        for (int i=0; i<getNumOutputs(); ++i) {
          try {
            aseed_new[dir][i] = SX(output(i).sparsity());
            aseed_new[dir][i].set(aseed[dir][i]);
          } catch(exception& ex) {
            stringstream ss;
            ss << "SXFunctionInternal::evalSX: Failed to set adjoint seed of "
               << output_.scheme.describeOutput(i) << ", direction " << dir << ": " << ex.what();
            throw CasadiException(ss.str());
          }
        }
      }
      evalSX(arg, res, fseed, fsens, aseed_new, asens);
      return;
    }

    // Resize (if needed) the number of outputs and make sure that the sparsity
    // pattern is correct (cheap if already ok)
    res.resize(getNumOutputs());
    for (int i=0; i<getNumOutputs(); ++i) {
      if (res[i].sparsity()!=output(i).sparsity()) {
        res[i] = SX(output(i).sparsity());
      }
    }

    // Resize (if needed) the number of forward sensitivities and make sure that
    // the sparsity pattern is correct (cheap if already ok)
    fsens.resize(nfdir);
    for (int dir=0; dir<nfdir; ++dir) {
      fsens[dir].resize(getNumOutputs());
      for (int i=0; i<getNumOutputs(); ++i) {
        if (fsens[dir][i].sparsity()!=output(i).sparsity()) {
          fsens[dir][i] = SX(output(i).sparsity());
        }
      }
    }

    // Resize (if needed) the number of adjoint sensitivities and make sure that
    // the sparsity pattern is correct (cheap if already ok)
    asens.resize(nadir);
    for (int dir=0; dir<nadir; ++dir) {
      asens[dir].resize(getNumInputs());
      for (int i=0; i<getNumInputs(); ++i) {
        if (asens[dir][i].sparsity()!=input(i).sparsity()) {
          asens[dir][i] = SX(input(i).sparsity());
        }
      }
    }

    // Call the sparse version
    evalSXsparse(arg, res, fseed, fsens, aseed, asens);
  }

  void FunctionInternal::evalSXsparse(const std::vector<SX>& arg, std::vector<SX>& res,
                                      const std::vector<std::vector<SX> >& fseed,
                                      std::vector<std::vector<SX> >& fsens,
                                      const std::vector<std::vector<SX> >& aseed,
                                      std::vector<std::vector<SX> >& asens) {
    casadi_error("FunctionInternal::evalSXsparse not defined for class " << typeid(*this).name());
  }

  void FunctionInternal::evalMX(const std::vector<MX>& arg, std::vector<MX>& res,
                                const std::vector<std::vector<MX> >& fseed,
                                std::vector<std::vector<MX> >& fsens,
                                const std::vector<std::vector<MX> >& aseed,
                                std::vector<std::vector<MX> >& asens) {
    MXFunction f = wrapMXFunction();
    f.init();
    f.callDerivative(arg, res, fseed, fsens, aseed, asens, true);
  }

  void FunctionInternal::spEvaluate(bool fwd) {
    // By default, everything is assumed to depend on everything

    // Variable which depends on all everything
    bvec_t all_depend(0);
    if (fwd) {
      // Get dependency on all inputs
      for (int iind=0; iind<getNumInputs(); ++iind) {
        const DMatrix& m = inputNoCheck(iind);
        const bvec_t* v = reinterpret_cast<const bvec_t*>(m.ptr());
        for (int i=0; i<m.size(); ++i) {
          all_depend |= v[i];
        }
      }

      // Propagate to all outputs
      for (int oind=0; oind<getNumOutputs(); ++oind) {
        DMatrix& m = outputNoCheck(oind);
        bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
        for (int i=0; i<m.size(); ++i) {
          v[i] = all_depend;
        }
      }

    } else {

      // Get dependency on all outputs
      for (int oind=0; oind<getNumOutputs(); ++oind) {
        const DMatrix& m = outputNoCheck(oind);
        const bvec_t* v = reinterpret_cast<const bvec_t*>(m.ptr());
        for (int i=0; i<m.size(); ++i) {
          all_depend |= v[i];
        }
      }

      // Propagate to all inputs
      for (int iind=0; iind<getNumInputs(); ++iind) {
        DMatrix& m = inputNoCheck(iind);
        bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
        for (int i=0; i<m.size(); ++i) {
          v[i] |= all_depend;
        }
      }
    }
  }

  void FunctionInternal::spEvaluateViaJacSparsity(bool fwd) {
    if (fwd) {
      // Clear the outputs
      for (int oind = 0; oind < getNumOutputs(); ++oind) {
        // Get data array for output and clear it
        bvec_t *outputd = get_bvec_t(output(oind).data());
        fill_n(outputd, output(oind).size(), 0);
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
        if (sp.isNull() || sp.size() == 0)
          continue; // Skip if zero

        const int d1 = sp.size2();
        //const int d2 = sp.size1();
        const vector<int>& colind = sp.colind();
        const vector<int>& row = sp.row();

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

  Function FunctionInternal::derivative(int nfwd, int nadj) {
    // Quick return if 0x0
    if (nfwd==0 && nadj==0) return shared_from_this<Function>();

    // Check if there are enough forward directions allocated
    if (nfwd>=derivative_fcn_.size()) {
      derivative_fcn_.resize(nfwd+1);
    }

    // Check if there are enough adjoint directions allocated
    if (nadj>=derivative_fcn_[nfwd].size()) {
      derivative_fcn_[nfwd].resize(nadj+1);
    }

    // Quick return if already cached
    if (derivative_fcn_[nfwd][nadj].alive()) {
      return shared_cast<Function>(derivative_fcn_[nfwd][nadj].shared());
    }

    // Return value
    Function ret;

    // Generate if not already cached

    // Get the number of scalar inputs and outputs
    int num_in_scalar = getNumInputNonzeros();
    int num_out_scalar = getNumOutputNonzeros();

    // Adjoint mode penalty factor (adjoint mode is usually more expensive to calculate)
    int adj_penalty = 2;

    // Crude estimate of the cost of calculating the full Jacobian
    int full_jac_cost = std::min(num_in_scalar, adj_penalty*num_out_scalar);

    // Crude estimate of the cost of calculating the directional derivatives
    int der_dir_cost = nfwd + adj_penalty*nadj;

    // Check if it is cheaper to calculate the full Jacobian and then multiply
    if ((getOption("ad_mode")=="forward" && nadj>0) ||
        (getOption("ad_mode")=="reverse" && nfwd>0)) {
      ret = getDerivativeViaJac(nfwd, nadj);
    } else if (hasSetOption("derivative_generator")) {
      /// User-provided derivative generator function
      DerivativeGenerator dergen = getOption("derivative_generator");
      Function this_ = shared_from_this<Function>();
      ret = dergen(this_, nfwd, nadj, user_data_);
    } else if (2*full_jac_cost < der_dir_cost) {
      // Generate the Jacobian and then multiply to get the derivative
      //ret = getDerivativeViaJac(nfwd, nadj); // NOTE: Uncomment this line
                                              // (and remove the next line)
                                              // to enable this feature
      ret = getDerivative(nfwd, nadj);
    } else {
      // Generate a new function
      ret = getDerivative(nfwd, nadj);
    }

    // Give it a suitable name
    stringstream ss;
    ss << "derivative_" << getOption("name") << "_" << nfwd << "_" << nadj;
    ret.setOption("name", ss.str());

    // Names of inputs
    std::vector<std::string> i_names;
    i_names.reserve(getNumInputs()*(1+nfwd)+getNumOutputs()*nadj);

    // Nondifferentiated inputs
    for (int i=0; i<getNumInputs(); ++i) {
      i_names.push_back("der_" + input_.scheme.entryLabel(i));
    }

    // Forward seeds
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<getNumInputs(); ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << input_.scheme.entryLabel(i);
        i_names.push_back(ss.str());
      }
    }

    // Adjoint seeds
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<getNumOutputs(); ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << output_.scheme.entryLabel(i);
        i_names.push_back(ss.str());
      }
    }

    // Pass to return object
    ret.setInputScheme(i_names);

    // Names of outputs
    std::vector<std::string> o_names;
    o_names.reserve(getNumOutputs()*(1+nfwd)+getNumInputs()*nadj);

    // Nondifferentiated inputs
    for (int i=0; i<getNumOutputs(); ++i) {
      o_names.push_back("der_" + output_.scheme.entryLabel(i));
    }

    // Forward sensitivities
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<getNumOutputs(); ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << output_.scheme.entryLabel(i);
        o_names.push_back(ss.str());
      }
    }

    // Adjoint sensitivities
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<getNumInputs(); ++i) {
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
    int ind=0;
    for (int d=-1; d<nfwd; ++d) {
      for (int i=0; i<getNumInputs(); ++i, ++ind) {
        if (ret.input(ind).size()!=0 && ret.input(ind).sparsity()!=input(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " input " << ind << " \""
                       << i_names.at(ind) << "\". Expected " << input(i).dimString()
                       << " but got " << ret.input(ind).dimString());
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<getNumOutputs(); ++i, ++ind) {
        if (ret.input(ind).size()!=0 && ret.input(ind).sparsity()!=output(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " input " << ind <<
                       " \"" << i_names.at(ind) << "\". Expected " << output(i).dimString()
                       << " but got " << ret.input(ind).dimString());
        }
      }
    }

    // Consistency check for outputs
    ind=0;
    for (int d=-1; d<nfwd; ++d) {
      for (int i=0; i<getNumOutputs(); ++i, ++ind) {
        if (ret.output(ind).size()!=0 && ret.output(ind).sparsity()!=output(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " output " << ind <<
                       " \"" <<  o_names.at(ind) << "\". Expected " << output(i).dimString()
                       << " but got " << ret.output(ind).dimString());
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<getNumInputs(); ++i, ++ind) {
        if (ret.output(ind).size()!=0 && ret.output(ind).sparsity()!=input(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " output " << ind << " \""
                       << o_names.at(ind) << "\". Expected " << input(i).dimString()
                       << " but got " << ret.output(ind).dimString());
        }
      }
    }

    // Save to cache
    derivative_fcn_[nfwd][nadj] = ret;

    // Return generated function
    return ret;
  }

  void FunctionInternal::setDerivative(const Function& fcn, int nfwd, int nadj) {

    // Check if there are enough forward directions allocated
    if (nfwd>=derivative_fcn_.size()) {
      derivative_fcn_.resize(nfwd+1);
    }

    // Check if there are enough adjoint directions allocated
    if (nadj>=derivative_fcn_[nfwd].size()) {
      derivative_fcn_[nfwd].resize(nadj+1);
    }

    // Save to cache
    derivative_fcn_[nfwd][nadj] = fcn;
  }

  Function FunctionInternal::getDerivative(int nfwd, int nadj) {
    if (full_jacobian_.alive()) {
      return getDerivativeViaJac(nfwd, nadj);
    }

    casadi_error("FunctionInternal::getDerivative not defined for class " << typeid(*this).name());
  }

  Function FunctionInternal::getDerivativeViaJac(int nfwd, int nadj) {
    // Get, possibly generate, a full Jacobian function
    Function jfcn = fullJacobian();

    // Number of inputs and outputs
    const int n_in = getNumInputs();
    const int n_out = getNumOutputs();

    // Get an expression for the full Jacobian
    vector<MX> arg = symbolicInput();
    vector<MX> res = jfcn.call(arg);
    MX J = res.front().T();
    res.erase(res.begin());

    // Make room for the derivatives
    arg.reserve(n_in*(1+nfwd)+n_out*nadj);
    res.reserve(n_out*(1+nfwd)+n_in*nadj);

    // Temporary string
    stringstream ss;
    vector<MX> d;

    // Forward derivatives
    if (nfwd>0) {
      // Get col offsets for the horzsplit
      vector<int> offset(1, 0);
      for (int i=0; i<n_out; ++i) {
        offset.push_back(offset.back()+res[i].numel());
      }

      // Calculate one derivative at a time (alternative: all at once using horzcat/horzsplit)
      for (int dir=0; dir<nfwd; ++dir) {
        // Assemble the right hand side
        d.clear();
        for (int i=0; i<n_in; ++i) {
          // Create a forward seed with a suitable name and add to list of inputs
          ss.str("fwd");
          ss << dir << "_" << arg[i];
          arg.push_back(MX::sym(ss.str(), arg[i].sparsity()));

          // Add to the right-hand-side under construction
          d.push_back(transpose(vec(arg.back())));
        }
        MX d_all = horzcat(d);

        // Calculate the derivatives using a matrix multiplication with the Jacobian
        d_all = mul(d_all, J);

        // Split up the left hand sides
        d = horzsplit(d_all, offset);
        for (int i=0; i<n_out; ++i) {
          res.push_back(reshape(d[i], res[i].shape()));
          res.back() = res.back().setSparse(res[i].sparsity()+res.back().sparsity());
        }
      }
    }

    // Adjoint derivatives
    if (nadj>0) {
      // Transpose of J
      MX JT = J.T();

      // Get col offsets for the horzsplit
      vector<int> offset(1, 0);
      for (int i=0; i<n_in; ++i) {
        offset.push_back(offset.back()+arg[i].numel());
      }

      // Calculate one derivative at a time (alternative: all at once using horzcat/horzsplit)
      for (int dir=0; dir<nadj; ++dir) {
        // Assemble the right hand side
        d.clear();
        for (int i=0; i<n_out; ++i) {
          // Create a forward seed with a suitable name and add to list of inputs
          ss.str("adj");
          ss << dir << "_" << res[i];
          arg.push_back(MX::sym(ss.str(), res[i].sparsity()));

          // Add to the right-hand-side under construction
          d.push_back(transpose(vec(arg.back())));
        }
        MX d_all = horzcat(d);

        // Calculate the derivatives using a matrix multiplication with the Jacobian
        d_all = mul(d_all, JT);

        // Split up the left hand sides
        d = horzsplit(d_all, offset);
        for (int i=0; i<n_in; ++i) {
          res.push_back(reshape(d[i], arg[i].shape()));
          res.back() = res.back().setSparse(arg[i].sparsity()+res.back().sparsity());
        }
      }
    }

    // Assemble the derivative function
    MXFunction ret(arg, res);
    return ret;
  }

  int FunctionInternal::getNumInputNonzeros() const {
    int ret=0;
    for (int iind=0; iind<getNumInputs(); ++iind) {
      ret += input(iind).size();
    }
    return ret;
  }

  int FunctionInternal::getNumOutputNonzeros() const {
    int ret=0;
    for (int oind=0; oind<getNumOutputs(); ++oind) {
      ret += output(oind).size();
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
                              const MXVectorVector& fseed, MXVectorVector& fsens,
                              const MXVectorVector& aseed, MXVectorVector& asens,
                              bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");

    // Lo logic for inlining yet
    bool inline_function = always_inline;

    if (inline_function) {
      // Evaluate the function symbolically
      evalMX(arg, res, fseed, fsens, aseed, asens);

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
          << std::endl << input_.scheme.describeInput(i) <<  " has shape (" << input(i).size2()
          << ", " << input(i).size1() << ") while a shape (" << arg[i].size2() << ", "
          << arg[i].size1() << ") was supplied.");
      }
      createCall(arg, res, fseed, fsens, aseed, asens);
    }
  }

  void FunctionInternal::call(const SXVector& arg, SXVector& res,
                        const SXVectorVector& fseed, SXVectorVector& fsens,
                        const SXVectorVector& aseed, SXVectorVector& asens,
                        bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    casadi_assert_message(!never_inline, "SX expressions do not support call-nodes");
    evalSX(arg, res, fseed, fsens, aseed, asens);
  }

  void FunctionInternal::call(const DMatrixVector& arg, DMatrixVector& res,
                        const DMatrixVectorVector& fseed, DMatrixVectorVector& fsens,
                        const DMatrixVectorVector& aseed, DMatrixVectorVector& asens,
                        bool always_inline, bool never_inline) {
    if (fseed.size()==0 || aseed.size()==0) {
      for (int i=0;i<arg.size();++i) {
        setInput(arg[i], i);
      }
      evaluate();
      res.resize(getNumOutputs());
      for (int i=0;i<res.size();++i) {
        res[i]=output(i);
      }
    } else {
      casadi_error("Not implemented");
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
      // Generate a new Jacobian
      Function ret;
      if (getNumInputs()==1 && getNumOutputs()==1) {
        ret = jacobian(0, 0, true, false);
      } else {
        ret = getFullJacobian();

        // Give it a suitable name
        stringstream ss;
        ss << "jacobian_" << getOption("name");
        ret.setOption("name", ss.str());
        ret.setOption("verbose", getOption("verbose"));

        // Same input scheme
        ret.setInputScheme(input_.scheme);

        // Construct output scheme
        std::vector<std::string> oscheme;
        oscheme.reserve(ret.getNumOutputs());
        oscheme.push_back("jac");
        for (int i=0; i<getNumOutputs(); ++i) {
          oscheme.push_back(output_.scheme.entryLabel(i));
        }
        ret.setOutputScheme(oscheme);
      }


      // Initialize
      ret.init();

      // Return and cache it for reuse
      full_jacobian_ = ret;
      return ret;
    }
  }

  Function FunctionInternal::getFullJacobian() {
    // Symbolic inputs and outputs of the full Jacobian function under construction
    vector<MX> jac_argv = symbolicInput(), jac_resv;

    // Symbolic inputs and outputs of the SISO function formed for generating the Jacobian
    MX arg, res;

    // Arguments and results when calling "this" with the SISO inputs
    vector<MX> argv, resv;

    // Check if the function is already single input
    if (getNumInputs()==1) {
      // Reuse the same symbolic primitives
      argv = jac_argv;
      arg = argv.front();
      resv = symbolicOutput(argv);
    } else {
      // Need to create a new symbolic primitive

      // Sparsity pattern for the symbolic input
      Sparsity sp_arg = Sparsity::sparse(1, 0);

      // Append the sparsity patterns, keep track of col offsets
      vector<int> col_offset(1, 0);
      for (int i=0; i<getNumInputs(); ++i) {
        sp_arg.appendColumns(input(i).sparsity().reshape(1, input(i).numel()));
        col_offset.push_back(sp_arg.numel());
      }

      // Create symbolic primitive
      arg = MX::sym("x", sp_arg);

      // Split up and fix shape
      argv = horzsplit(arg, col_offset);
      for (int i=0; i<getNumInputs(); ++i) {
        argv[i] = reshape(argv[i], input(i).sparsity());
      }

      // Evaluate symbolically
      resv = shared_from_this<Function>().call(argv);
    }

    // Check if the function is already single output
    if (getNumOutputs()==1) {
      // Reuse the same output
      res = resv.front();
    } else {
      // Concatenate the outputs
      vector<MX> tmp = resv;
      for (vector<MX>::iterator i=tmp.begin(); i!=tmp.end(); ++i) {
        *i = reshape(*i, 1, i->numel());
      }
      res = horzcat(tmp);
    }

    // Create a SISO function
    MXFunction f(arg, res);
    f.setOption("ad_mode", getOption("ad_mode"));
    f.init();

    // Form an expression for the full Jacobian
    MX J = f.jac();

    // Append to list of outputs
    resv.insert(resv.begin(), J);

    // Wrap in an MXFunction to get the correct input, if necessary
    if (getNumInputs()==1) {
      // No need to wrap
      jac_resv = resv;
    } else {
      // Form a SIMO function
      MXFunction J_simo(arg, resv);
      J_simo.setOption("name", "J_simo");
      J_simo.init();

      // The inputs of J_simo in terms of jac_argv
      vector<MX> tmp = jac_argv;
      for (vector<MX>::iterator i=tmp.begin(); i!=tmp.end(); ++i) {
        *i = reshape(*i, 1, i->numel());
      }
      MX J_simo_arg = horzcat(tmp);

      // Evaluate symbolically
      jac_resv = J_simo(J_simo_arg);
    }

    // We are now ready to form the full Jacobian
    MXFunction ret(jac_argv, jac_resv);
    return ret;

  }

  void FunctionInternal::generateCode(std::ostream &cfile, bool generate_main) {
    assertInit();

    // set stream parameters
    cfile.precision(std::numeric_limits<double>::digits10+2);
    // This is really only to force a decimal dot, would be better if it can be avoided
    cfile << std::scientific;

    // Print header
    cfile << "/* This function was automatically generated by CasADi */" << std::endl;

    // Create a code generator object
    CodeGenerator gen;

    // Add standard math
    gen.addInclude("math.h");

    // Generate function inputs and outputs information
    generateIO(gen);

    // Generate the actual function
    generateFunction(gen.function_, "evaluate", "const d*", "d*", "d", gen);

    // Flush the code generator
    gen.flush(cfile);

    // Define wrapper function
    cfile << "int evaluateWrap(const d** x, d** r) {" << std::endl;
    cfile << "  evaluate(";

    // Number of inputs/outputs
    int n_i = input_.data.size();
    int n_o = output_.data.size();

    // Pass inputs
    for (int i=0; i<n_i; ++i) {
      if (i!=0) cfile << ", ";
      cfile << "x[" << i << "]";
    }

    // Pass outputs
    for (int i=0; i<n_o; ++i) {
      if (i+n_i!= 0) cfile << ", ";
      cfile << "r[" << i << "]";
    }

    cfile << "); " << std::endl;
    cfile << "  return 0;" << std::endl;
    cfile << "}" << std::endl << std::endl;

    // Create a main for debugging and profiling: TODO: Cleanup and expose to user, see #617
    if (generate_main) {
      int n_in = getNumInputs();
      int n_out = getNumOutputs();

      cfile << "#include <stdio.h>" << std::endl;
      cfile << "int main() {" << std::endl;
      cfile << "  int i, j;" << std::endl;

      // Declare input buffers
      for (int i=0; i<n_in; ++i) {
        cfile << "  d t_x" << i << "[" << input(i).sparsity().size() << "];" << std::endl;
      }

      // Declare output buffers
      for (int i=0; i<n_out; ++i) {
        cfile << "  d t_r" << i << "[" << output(i).sparsity().size() << "];" << std::endl;
      }

      // Repeat 10 times
      cfile << "  for (j=0; j<10; ++j) {" << std::endl;

      // Dummy input values
      for (int i=0; i<n_in; ++i) {
        cfile << "    for (i=0; i<" << input(i).sparsity().size() << "; ++i) t_x"
              << i << "[i] = sin(2.2*i+sqrt(4.3/(j+1)));" << std::endl;
      }

      // Pass inputs
      cfile << "    evaluate(";
      for (int i=0; i<n_in; ++i) {
        cfile << "t_x" << i;
        if (i+1<n_in+n_out)
          cfile << ", ";
      }

      // Pass output buffers
      for (int i=0; i<n_out; ++i) {
        cfile << "t_r" << i;
        if (i+1<n_out)
          cfile << ", ";
      }
      cfile << "); " << std::endl;


      // Dummy printout
      for (int i=0; i<n_out; ++i) {
        int n = output(i).sparsity().size();
        for (int j=0; j<n && j<5; ++j) {
          cfile << "    printf(\"%g \", t_r" << i << "[" << j << "]);" << std::endl;
        }
        cfile << "    printf(\"\\n\");" << std::endl;
      }

      // End repeat
      cfile << "  }" << std::endl;

      cfile << "  return 0;" << std::endl;
      cfile << "}" << std::endl << std::endl;
    }
  }

  void FunctionInternal::generateFunction(
      std::ostream &stream, const std::string& fname, const std::string& input_type,
      const std::string& output_type, const std::string& type, CodeGenerator& gen) const {

    // Generate declarations
    generateDeclarations(stream, type, gen);

    // Number of inpus and outputs
    int n_in = getNumInputs();
    int n_out = getNumOutputs();

    // Define function
    stream << "/* " << getSanitizedName() << " */" << std::endl;
    stream << "void " << fname << "(";

    // Declare inputs
    for (int i=0; i<n_in; ++i) {
      stream << input_type << " x" << i;
      if (i+1<n_in+n_out)
        stream << ", ";
    }

    // Declare outputs
    for (int i=0; i<n_out; ++i) {
      stream << output_type << " r" << i;
      if (i+1<n_out)
        stream << ", ";
    }
    stream << ") { " << std::endl;

    // Insert the function body
    generateBody(stream, type, gen);

    // Finalize the function
    stream << "}" << std::endl;
    stream << std::endl;
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

  void FunctionInternal::generateIO(CodeGenerator& gen) {
    // Short-hands
    int n_i = input_.data.size();
    int n_o = output_.data.size();
    int n_io = n_i + n_o;
    stringstream &s = gen.function_;

    // Function that returns the number of inputs and outputs
    s << "int init(int *n_in, int *n_out) {" << endl;
    s << "  *n_in = " << n_i << ";" << endl;
    s << "  *n_out = " << n_o << ";" << endl;
    s << "  return 0;" << endl;
    s << "}" << endl;
    s << endl;

    // Inputs and outputs for each parsity index
    std::multimap<int, int> io_sparsity_index;

    // The number of patterns needed for the inputs and outputs
    int num_io_patterns = 0;

    // Get the sparsity pattern of all inputs
    for (int i=0; i<n_io; ++i) {
      // Get the sparsity pattern
      const Sparsity& sp = i<n_i ? input(i).sparsity() : output(i-n_i).sparsity();

      // Print the sparsity pattern or retrieve the index of an existing pattern
      int ind = gen.addSparsity(sp);
      num_io_patterns = std::max(num_io_patterns, ind+1);

      // Store the index
      io_sparsity_index.insert(std::pair<int, int>(ind, i));
    }

    // Function that returns the sparsity pattern
    s << "int getSparsity(int i, int *nrow, int *ncol, int **colind, int **row) {" << endl;

    // Get the sparsity index using a switch
    s << "  int* sp;" << endl;
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
      s << "      sp = s" << i << "; break;" << endl;
    }

    // Finalize the switch
    s << "    default:" << endl;
    s << "      return 1;" << endl;
    s << "  }" << endl << endl;

    // Decompress the sparsity pattern
    s << "  *nrow = sp[0];" << endl;
    s << "  *ncol = sp[1];" << endl;
    s << "  *colind = sp + 2;" << endl;
    s << "  *row = sp + 2 + (*ncol + 1);" << endl;
    s << "  return 0;" << endl;
    s << "}" << endl << endl;
  }

  void FunctionInternal::assignIgnore(MX& y, const MX& x, const std::vector<int>& nz) {
    y[nz] = x;
  }

  void FunctionInternal::assignIgnore(SX& y, const SX& x, const std::vector<int>& nz) {
    vector<SXElement>& y_data = y.data();
    const vector<SXElement>& x_data = x.data();
    casadi_assert(nz.size()==x_data.size());
    for (int k=0; k<nz.size(); ++k) {
      if (nz[k]>=0) {
        y_data.at(nz[k]) = x_data.at(k);
      }
    }
  }

  Function FunctionInternal::dynamicCompilation(Function f, std::string fname, std::string fdescr,
                                                std::string compiler) {
#ifdef WITH_DL

    // Flag to get a DLL
#ifdef __APPLE__
    string dlflag = " -dynamiclib";
#else // __APPLE__
    string dlflag = " -shared";
#endif // __APPLE__

    // Check if f is initialized
    bool f_is_init = f.isInit();
    if (!f_is_init) f.init();

    // Filenames
    string cname = fname + ".c";
    string dlname = fname + ".so";

    // Remove existing files, if any
    string rm_command = "rm -rf " + cname + " " + dlname;
    int flag = system(rm_command.c_str());
    casadi_assert_message(flag==0, "Failed to remove old source");

    // Codegen it
    f.generateCode(cname);
    if (verbose_) {
      cout << "Generated c-code for " << fdescr << " (" << cname << ")" << endl;
    }

    // Compile it
    string compile_command = compiler + " " + dlflag + " " + cname + " -o " + dlname;
    if (verbose_) {
      cout << "Compiling " << fdescr <<  " using \"" << compile_command << "\"" << endl;
    }

    time_t time1 = time(0);
    flag = system(compile_command.c_str());
    time_t time2 = time(0);
    double comp_time = difftime(time2, time1);
    casadi_assert_message(flag==0, "Compilation failed");
    if (verbose_) {
      cout << "Compiled " << fdescr << " (" << dlname << ") in " << comp_time << " s."  << endl;
    }

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
#else // WITH_DL
    casadi_error("Codegen in SCPgen requires CasADi to be compiled "
                 "with option \"WITH_DL\" enabled");
#endif // WITH_DL
  }

  void FunctionInternal::createCall(const std::vector<MX> &arg,
                          std::vector<MX> &res, const std::vector<std::vector<MX> > &fseed,
                          std::vector<std::vector<MX> > &fsens,
                          const std::vector<std::vector<MX> > &aseed,
                          std::vector<std::vector<MX> > &asens) {

    if (fseed.empty() && aseed.empty()) {
      // Create the evaluation node
      res = callSelf(arg);
    } else {
      // Create derivative node
      createCallDerivative(arg, res, fseed, fsens, aseed, asens);
    }
  }

  std::vector<MX> FunctionInternal::callSelf(const std::vector<MX> &arg) {
    return MX::createMultipleOutput(new CallFunction(shared_from_this<Function>(), arg));
  }

  void FunctionInternal::createCallDerivative(
      const std::vector<MX>& arg, std::vector<MX>& res,
      const std::vector<std::vector<MX> >& fseed, std::vector<std::vector<MX> >& fsens,
      const std::vector<std::vector<MX> >& aseed, std::vector<std::vector<MX> >& asens) {


    // Number of directional derivatives
    int nfdir = fseed.size();
    int nadir = aseed.size();

    // Number inputs and outputs
    int num_in = getNumInputs();
    int num_out = getNumOutputs();

    // Create derivative function
    Function dfcn = derivative(nfdir, nadir);

    // All inputs
    vector<MX> darg;
    darg.reserve(num_in*(1+nfdir) + num_out*nadir);
    darg.insert(darg.end(), arg.begin(), arg.end());

    // Forward seeds
    for (int dir=0; dir<nfdir; ++dir) {
      darg.insert(darg.end(), fseed[dir].begin(), fseed[dir].end());
    }

    // Adjoint seeds
    for (int dir=0; dir<nadir; ++dir) {
      darg.insert(darg.end(), aseed[dir].begin(), aseed[dir].end());
    }

    // Create the evaluation node
    vector<MX> x = MX::createMultipleOutput(new CallFunction(dfcn, darg));
    vector<MX>::iterator x_it = x.begin();

    // Create the output nodes corresponding to the nondifferented function
    res.resize(num_out);
    for (int i = 0; i < num_out; ++i) {
      res[i] = *x_it++;
    }

    // Forward sensitivities
    fsens.resize(nfdir);
    for (int dir = 0; dir < nfdir; ++dir) {
      fsens[dir].resize(num_out);
      for (int i = 0; i < num_out; ++i) {
        fsens[dir][i] = *x_it++;
      }
    }

    // Adjoint sensitivities
    asens.resize(nadir);
    for (int dir = 0; dir < nadir; ++dir) {
      asens[dir].resize(num_in);
      for (int i = 0; i < num_in; ++i) {
        asens[dir][i] = *x_it++;
      }
    }

  }

  void FunctionInternal::evaluateMX(
      MXNode* node, const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens,
      const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given) {
    // Collect inputs and seeds
    vector<MX> arg = MXNode::getVector(input);
    vector<vector<MX> > fseed = MXNode::getVector(fwdSeed);
    vector<vector<MX> > aseed = MXNode::getVector(adjSeed);

    // Free adjoint seeds
    MXNode::clearVector(adjSeed);

    // Evaluate symbolically
    vector<MX> res;
    vector<vector<MX> > fsens, asens;

    if (fwdSens.size()==0 && adjSens.size()==0) {
      res = callSelf(arg);
    } else {
      createCallDerivative(arg, res, fseed, fsens, aseed, asens);
    }

    // Store the non-differentiated results
    if (!output_given) {
      for (int i=0; i<res.size(); ++i) {
        if (output[i]!=0) {
          *output[i] = res[i];
        }
      }
    }

    // Store the forward sensitivities
    for (int d=0; d<fwdSens.size(); ++d) {
      for (int i=0; i<fwdSens[d].size(); ++i) {
        if (fwdSens[d][i]!=0) {
          *fwdSens[d][i] = fsens[d][i];
        }
      }
    }

    // Store the adjoint sensitivities
    for (int d=0; d<adjSens.size(); ++d) {
      for (int i=0; i<adjSens[d].size(); ++i) {
        if (adjSens[d][i]!=0 && !asens[d][i].isEmpty(true) && !(*adjSens[d][i]).isEmpty(true)) {
          adjSens[d][i]->addToSum(asens[d][i]);
        }
      }
    }
  }

  void FunctionInternal::propagateSparsity(MXNode* node, DMatrixPtrV& arg, DMatrixPtrV& res,
                                           std::vector<int>& itmp,
                                           std::vector<double>& rtmp, bool use_fwd) {
    // Pass/clear forward seeds/adjoint sensitivities
    for (int iind = 0; iind < getNumInputs(); ++iind) {
      // Input vector
      vector<double> &v = input(iind).data();
      if (v.empty()) continue; // FIXME: remove?

      if (arg[iind] == 0) {
        // Set to zero if not used
        fill_n(get_bvec_t(v), v.size(), bvec_t(0));
      } else {
        // Copy output
        input(iind).sparsity().set(
          get_bvec_t(input(iind).data()), get_bvec_t(arg[iind]->data()), arg[iind]->sparsity());
      }
    }

    // Pass/clear adjoint seeds/forward sensitivities
    for (int oind = 0; oind < getNumOutputs(); ++oind) {
      // Output vector
      vector<double> &v = output(oind).data();
      if (v.empty()) continue; // FIXME: remove?
      if (res[oind] == 0) {
        // Set to zero if not used
        fill_n(get_bvec_t(v), v.size(), bvec_t(0));
      } else {
        // Copy output
        output(oind).sparsity().set(
          get_bvec_t(output(oind).data()), get_bvec_t(res[oind]->data()), res[oind]->sparsity());
        if (!use_fwd) fill_n(get_bvec_t(res[oind]->data()), res[oind]->size(), bvec_t(0));
      }
    }

    // Propagate seeds
    spInit(use_fwd); // NOTE: should only be done once
    if (spCanEvaluate(use_fwd)) {
      spEvaluate(use_fwd);
    } else {
      spEvaluateViaJacSparsity(use_fwd);
    }

    // Get the sensitivities
    if (use_fwd) {
      for (int oind = 0; oind < res.size(); ++oind) {
        if (res[oind] != 0) {
          res[oind]->sparsity().set(
            get_bvec_t(res[oind]->data()),
            get_bvec_t(output(oind).data()),
            output(oind).sparsity());
        }
      }
    } else {
      for (int iind = 0; iind < arg.size(); ++iind) {
        if (arg[iind] != 0) {
          arg[iind]->sparsity().bor(
            get_bvec_t(arg[iind]->data()), get_bvec_t(input(iind).data()), input(iind).sparsity());
        }
      }
    }

    // Clear seeds and sensitivities
    for (int iind = 0; iind < arg.size(); ++iind) {
      vector<double> &v = input(iind).data();
      fill(v.begin(), v.end(), 0);
    }
    for (int oind = 0; oind < res.size(); ++oind) {
      vector<double> &v = output(oind).data();
      fill(v.begin(), v.end(), 0);
    }
  }

  void FunctionInternal::generateOperation(const MXNode* node, std::ostream &stream,
                                           const std::vector<std::string>& arg,
                                           const std::vector<std::string>& res,
                                           CodeGenerator& gen) const {

    // Running index of the temporary used
    int nr=0;

    // Copy arguments with nonmatching sparsities to the temp vector
    vector<string> arg_mod = arg;
    for (int i=0; i<getNumInputs(); ++i) {
      if (node->dep(i).sparsity()!=input(i).sparsity()) {
        arg_mod[i] = "rrr+" + CodeGenerator::numToString(nr);
        nr += input(i).size();

        // Codegen "copy sparse"
        gen.addAuxiliary(CodeGenerator::AUX_COPY_SPARSE);

        int sp_arg = gen.getSparsity(node->dep(i).sparsity());
        int sp_input = gen.addSparsity(input(i).sparsity());
        stream << "  casadi_copy_sparse(" << arg[i] << ", s" << sp_arg << ", " << arg_mod[i]
               << ", s" << sp_input << ");" << std::endl;
      }
    }

    // Get the index of the function
    int f = gen.getDependency(shared_from_this<Function>());
    stream << "  f" << f << "(";

    // Pass inputs to the function input buffers
    for (int i=0; i<arg.size(); ++i) {
      stream << arg_mod.at(i);
      if (i+1<arg.size()+res.size()) stream << ", ";
    }

    // Separate arguments and results with an extra space
    stream << " ";

    // Pass results to the function input buffers
    for (int i=0; i<res.size(); ++i) {
      stream << res.at(i);
      if (i+1<res.size()) stream << ", ";
    }

    // Finalize the function call
    stream << ");" << endl;
  }

  void FunctionInternal::nTmp(MXNode* node, size_t& ni, size_t& nr) {
    // Start with no extra memory
    ni=0;
    nr=0;

    // Add memory for all inputs with nonmatching sparsity
    for (int i=0; i<getNumInputs(); ++i) {
      if (node->dep(i).isNull() || node->dep(i).sparsity()!=input(i).sparsity()) {
        nr += input(i).size();
      }
    }
  }

  void FunctionInternal::evaluateSX(MXNode* node, const SXPtrV& arg, SXPtrV& res,
                                    std::vector<int>& itmp, std::vector<SXElement>& rtmp) {

    // Create input arguments
    vector<SX> argv(arg.size());
    for (int i=0; i<arg.size(); ++i) {
      argv[i] = SX(input(i).sparsity(), 0.);
      if (arg[i] != 0)
        argv[i].set(*arg[i]);
    }

    // Evaluate symbolically
    vector<SX> resv;
    vector<vector<SX> > dummy;
    evalSX(argv, resv, dummy, dummy, dummy, dummy);

    // Collect the result
    for (int i = 0; i < res.size(); ++i) {
      if (res[i] != 0)
        *res[i] = resv[i];
    }
  }

  void FunctionInternal::evaluateD(MXNode* node, const DMatrixPtrV& arg, DMatrixPtrV& res,
                                   std::vector<int>& itmp, std::vector<double>& rtmp) {

    // Set up timers for profiling
    double time_zero=0;
    double time_start=0;
    double time_stop=0;
    double time_offset=0;

    if (CasadiOptions::profiling) {
      time_start = getRealTime(); // Start timer
    }

    // Number of inputs and outputs
    int num_in = getNumInputs();
    int num_out = getNumOutputs();

    // Pass the inputs to the function
    for (int i = 0; i < num_in; ++i) {
      DMatrix *a = arg[i];
      if (a != 0) {
        setInput(*a, i);
      } else {
        setInput(0., i);
      }
    }

    if (CasadiOptions::profiling) {
      time_zero = getRealTime();
    }

    // Evaluate
    evaluate();
    if (CasadiOptions::profiling) {
      time_offset += getRealTime() - time_zero;
    }

    // Get the outputs
    for (int i = 0; i < num_out; ++i) {
      if (res[i] != 0) getOutput(*res[i], i);
    }

    // Write out profiling information
    if (CasadiOptions::profiling) {
      time_stop = getRealTime();
      if (CasadiOptions::profilingBinary) {

      } else {
        CasadiOptions::profilingLog
            << "overhead " << this << ":" <<getOption("name") << "|"
            << (time_stop-time_start-time_offset)*1e6 << " ns" << std::endl;
      }
    }
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
    if (all(v <= ub + tol) && all(v >= lb - tol)) {
      stream << "All " << v.size() << " constraints on " << name << " are met: " << endl;
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
    for (int i=0;i<v.size();i++) {

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


} // namespace casadi
