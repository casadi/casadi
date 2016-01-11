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
#include "../std_vector_tools.hpp"

#include <typeinfo>
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
    eval_ = 0;
    simple_ = 0;
    user_data_ = 0;
    monitor_inputs_ = false;
    monitor_outputs_ = false;

    sz_arg_tmp_ = 0;
    sz_res_tmp_ = 0;
    sz_iw_tmp_ = 0;
    sz_w_tmp_ = 0;

    sz_arg_per_ = 0;
    sz_res_per_ = 0;
    sz_iw_per_ = 0;
    sz_w_per_ = 0;
  }

  FunctionInternal::~FunctionInternal() {
    for (auto&& i : mem_) delete i;
    mem_.clear();
  }

  void FunctionInternal::init() {
    setDefaultOptions();

    // Free existing memory object, if any
    for (auto&& i : mem_) delete i;
    mem_.clear();

    // Get the number of inputs and outputs
    isp_.resize(get_n_in());
    osp_.resize(get_n_out());
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Get the input and output sparsities
    for (int i=0; i<n_in; ++i) isp_[i] = get_sparsity_in(i);
    for (int i=0; i<n_out; ++i) osp_[i] = get_sparsity_out(i);

    // Allocate memory for function inputs and outputs
    sz_arg_per_ += n_in;
    sz_res_per_ += n_out;

    verbose_ = option("verbose");
    jit_ = option("jit");
    compilerplugin_ = option("compiler").to_string();
    if (hasSetOption("jit_options")) jit_options_ = option("jit_options");
    regularity_check_ = option("regularity_check");

    // Warn for functions with too many inputs or outputs
    casadi_assert_warning(n_in<10000, "Function " << name_
                          << " has a large number of inputs (" << n_in << "). "
                          "Changing the problem formulation is strongly encouraged.");
    casadi_assert_warning(n_out<10000, "Function " << name_
                          << " has a large number of outputs (" << n_out << "). "
                          "Changing the problem formulation is strongly encouraged.");

    // Resize the matrix that holds the sparsity of the Jacobian blocks
    jac_sparsity_ = jac_sparsity_compact_ =
        SparseStorage<Sparsity>(Sparsity(n_out, n_in));
    jac_ = jac_compact_ = SparseStorage<WeakRef>(Sparsity(n_out, n_in));

    if (hasSetOption("user_data")) {
      user_data_ = option("user_data").to_void_pointer();
    }

    // Pass monitors
    if (hasSetOption("monitor")) {
      const std::vector<std::string>& monitors = option("monitor");
      for (std::vector<std::string>::const_iterator it=monitors.begin();it!=monitors.end();it++) {
        monitors_.insert(*it);
      }
    }

    // Custom input scheme
    if (hasSetOption("input_scheme")) {
      ischeme_ = option("input_scheme");
    }

    // If input scheme empty, provide default names
    if (ischeme_.empty()) {
      ischeme_.resize(n_in);
      for (size_t i=0; i!=ischeme_.size(); ++i) {
        ischeme_[i] = "i" + CodeGenerator::to_string(i);
      }
    }

    // Custom output scheme
    if (hasSetOption("output_scheme")) {
      oscheme_ = option("output_scheme");
    }

    // If output scheme null, provide default names
    if (oscheme_.empty()) {
      oscheme_.resize(n_out);
      for (size_t i=0; i!=oscheme_.size(); ++i) {
        oscheme_[i] = "o" + CodeGenerator::to_string(i);
      }
    }

    monitor_inputs_ = monitored("inputs");
    monitor_outputs_ = monitored("outputs");
    gather_stats_ = option("gather_stats");
    inputs_check_ = option("inputs_check");
    alloc_arg(0);
    alloc_res(0);
  }

  void FunctionInternal::finalize() {
    if (jit_) {
      CodeGenerator gen;
      gen.add(shared_from_this<Function>(), "jit_tmp");
      gen.generate("jit_tmp");
      compiler_ = Compiler("jit_tmp.c", compilerplugin_, jit_options_);

      // Try to load with simplified syntax
      simple_ = (simple_t)compiler_.getFunction("jit_tmp_simple");
      // If not succesful, try generic syntax
      if (simple_==0) {
        eval_ = (eval_t)compiler_.getFunction("jit_tmp");
        casadi_assert_message(eval_!=0, "Cannot load JIT'ed function.");
      }
    }

    // Create memory object
    casadi_assert(mem_.empty());
    Memory* m = memory();
    casadi_assert(m!=0);
    init_memory(*m);
    mem_.push_back(m);
  }

  void FunctionInternal::_eval(const double** arg, double** res, int* iw, double* w, int mem) {
    if (simplifiedCall()) {
      // Copy arguments to input buffers
      const double* arg1=w;
      for (int i=0; i<this->n_in(); ++i) {
        *w++ = arg[i] ? *arg[i] : 0;
      }

      // Evaluate
      if (simple_) {
        simple_(arg1, w);
      } else {
        simple(arg1, w);
      }

      // Get outputs
      for (int i=0; i<this->n_out(); ++i) {
        if (res[i]) *res[i] = *w;
        ++w;
      }
    } else {
      if (eval_) {
        eval_(arg, res, iw, w, mem);
      } else {
        eval(arg, res, iw, w, mem);
      }
    }
  }

  void FunctionInternal::_eval(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    eval_sx(arg, res, iw, w, mem);
  }

  void FunctionInternal::_eval(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    spFwd(arg, res, iw, w, mem);
  }

  void FunctionInternal::printDimensions(ostream &stream) const {
    stream << " Number of inputs: " << n_in() << endl;
    for (int i=0; i<n_in(); ++i) {
      stream << "  Input " << i  << ", a.k.a. \"" << name_in(i) << "\", "
             << size_in(i) << ", " << description_in(i) << endl;
    }
    stream << " Number of outputs: " << n_out() << endl;
    for (int i=0; i<n_out(); ++i) {
      stream << "  Output " << i  << ", a.k.a. \"" << name_out(i) << "\", "
             << size_out(i) << ", " << description_out(i) << endl;
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
    casadi_assert_message(sparsity_out(oind).is_scalar(),
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
    Dict opts = {{"input_scheme", ischeme_},
                 {"output_scheme", ionames},
                 {"jit", jit_},
                 {"compiler", compilerplugin_},
                 {"jit_options", jit_options_}};
    return getGradient(ss.str(), iind, oind, opts);
  }

  Function FunctionInternal::tangent(int iind, int oind) {
    // Assert scalar
    casadi_assert_message(sparsity_in(iind).is_scalar(),
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
    Dict opts = {{"input_scheme", ischeme_},
                 {"output_scheme", ionames},
                 {"jit", jit_},
                 {"compiler", compilerplugin_},
                 {"jit_options", jit_options_}};
    return getTangent(ss.str(), iind, oind, opts);
  }

  Function FunctionInternal::hessian(int iind, int oind) {
    log("FunctionInternal::hessian");

    // Assert scalar
    casadi_assert_message(sparsity_out(oind).is_scalar(),
                          "Only hessians of scalar functions allowed.");

    // Generate gradient function
    return getHessian(iind, oind);
  }

  Function FunctionInternal::getGradient(const std::string& name, int iind, int oind,
                                         const Dict& opts) {
    Function f = wrapMXFunction();
    return f.gradient(iind, oind);
  }

  Function FunctionInternal::getTangent(const std::string& name, int iind, int oind,
                                        const Dict& opts) {
    Function f = wrapMXFunction();
    return f.tangent(iind, oind);
  }

  Function FunctionInternal::wrapMXFunction() {
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
    vector<MX> arg = mx_in();
    vector<MX> res = shared_from_this<Function>()(arg);
    return Function("wrap_" + name_, arg, res, opts);
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

  std::vector<MX> FunctionInternal::symbolicOutput(const std::vector<MX>& arg) {
    return shared_from_this<Function>()(arg);
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

  // Traits
  template<bool fwd> struct JacSparsityTraits {};
  template<> struct JacSparsityTraits<true> {
    typedef const bvec_t* arg_t;
    static inline void sp(FunctionInternal *f, const bvec_t** arg, bvec_t** res,
                          int* iw, bvec_t* w, int mem) {
      f->spFwd(arg, res, iw, w, mem);
    }
  };
  template<> struct JacSparsityTraits<false> {
    typedef bvec_t* arg_t;
    static inline void sp(FunctionInternal *f, bvec_t** arg, bvec_t** res,
                          int* iw, bvec_t* w, int mem) {
      f->spAdj(arg, res, iw, w, mem);
    }
  };

  template<bool fwd>
  Sparsity FunctionInternal::getJacSparsityGen(int iind, int oind,
                                               bool symmetric, int gr_i, int gr_o) {
    // Number of nonzero inputs and outputs
    int nz_in = nnz_in(iind);
    int nz_out = nnz_out(oind);

    // Evaluation buffers
    vector<typename JacSparsityTraits<fwd>::arg_t> arg(sz_arg(), 0);
    vector<bvec_t*> res(sz_res(), 0);
    vector<int> iw(sz_iw());
    vector<bvec_t> w(sz_w(), 0);

    // Seeds and sensitivities
    vector<bvec_t> seed(nz_in, 0);
    arg[iind] = get_ptr(seed);
    vector<bvec_t> sens(nz_out, 0);
    res[oind] = get_ptr(sens);
    if (!fwd) std::swap(seed, sens);

    // Number of forward sweeps we must make
    int nsweep = seed.size() / bvec_size;
    if (seed.size() % bvec_size) nsweep++;

    // Print
    if (verbose()) {
      userOut() << "FunctionInternal::getJacSparsityGen<" << fwd << ">: "
                << nsweep << " sweeps needed for " << seed.size() << " directions" << endl;
    }

    // Progress
    int progress = -10;

    // Temporary vectors
    std::vector<int> jcol, jrow;

    // Loop over the variables, bvec_size variables at a time
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
      int ndir_local = seed.size()-offset;
      ndir_local = std::min(bvec_size, ndir_local);

      for (int i=0; i<ndir_local; ++i) {
        seed[offset+i] |= bvec_t(1)<<i;
      }

      // Propagate the dependencies
      JacSparsityTraits<fwd>::sp(this, get_ptr(arg), get_ptr(res), get_ptr(iw), get_ptr(w), 0);

      // Loop over the nonzeros of the output
      for (int el=0; el<sens.size(); ++el) {

        // Get the sparsity sensitivity
        bvec_t spsens = sens[el];

        if (!fwd) {
          // Clear the sensitivities for the next sweep
          sens[el] = 0;
        }

        // If there is a dependency in any of the directions
        if (spsens!=0) {

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
      for (int i=0; i<ndir_local; ++i) {
        seed[offset+i] = 0;
      }
    }

    // Construct sparsity pattern and return
    if (!fwd) swap(jrow, jcol);
    Sparsity ret = Sparsity::triplet(nz_out, nz_in, jcol, jrow);
    casadi_msg("Formed Jacobian sparsity pattern (dimension " << ret.size() << ", "
               << ret.nnz() << " nonzeros, " << (100.0*ret.nnz())/ret.numel() << " % nonzeros).");
    casadi_msg("FunctionInternal::getJacSparsity end ");
    return ret;
  }

  Sparsity FunctionInternal::getJacSparsityHierarchicalSymm(int iind, int oind) {
    casadi_assert(spCanEvaluate(true));

    // Number of nonzero inputs
    int nz = nnz_in(iind);
    casadi_assert(nz==nnz_out(oind));

    // Evaluation buffers
    vector<const bvec_t*> arg(sz_arg(), 0);
    vector<bvec_t*> res(sz_res(), 0);
    vector<int> iw(sz_iw());
    vector<bvec_t> w(sz_w());

    // Seeds
    vector<bvec_t> seed(nz, 0);
    arg[iind] = get_ptr(seed);

    // Sensitivities
    vector<bvec_t> sens(nz, 0);
    res[oind] = get_ptr(sens);

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

      Sparsity D = r.star_coloring();

      casadi_msg("Star coloring on " << r.dim() << ": " << D.size2() << " <-> " << D.size1());

      // Clear the seeds
      fill(seed.begin(), seed.end(), 0);

      // Subdivide the coarse block
      for (int k=0; k<coarse.size()-1; ++k) {
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
              bvec_toggle(get_ptr(seed), fine[fci+fci_start], fine[fci+fci_start+1],
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
            IM lookup = IM::triplet(lookup_row, lookup_col, lookup_value,
                                    bvec_size, coarse.size());

            std::reverse(lookup_col.begin(), lookup_col.end());
            std::reverse(lookup_row.begin(), lookup_row.end());
            std::reverse(lookup_value.begin(), lookup_value.end());
            IM duplicates =
              IM::triplet(lookup_row, lookup_col, lookup_value, bvec_size, coarse.size())
              - lookup;
            duplicates = sparsify(duplicates);
            lookup(duplicates.sparsity()) = -bvec_size;

            // Propagate the dependencies
            spFwd(get_ptr(arg), get_ptr(res), get_ptr(iw), get_ptr(w), 0);

            // Temporary bit work vector
            bvec_t spsens;

            // Loop over the cols of coarse blocks
            for (int cri=0; cri<coarse.size()-1; ++cri) {

              // Loop over the cols of fine blocks within the current coarse block
              for (int fri=fine_lookup[coarse[cri]];fri<fine_lookup[coarse[cri+1]];++fri) {
                // Lump individual sensitivities together into fine block
                bvec_or(get_ptr(sens), spsens, fine[fri], fine[fri+1]);

                // Loop over all bvec_bits
                for (int bvec_i=0;bvec_i<bvec_size;++bvec_i) {
                  if (spsens & (bvec_t(1) << bvec_i)) {
                    // if dependency is found, add it to the new sparsity pattern
                    int ind = lookup.sparsity().get_nz(bvec_i, cri);
                    casadi_assert(ind!=-1);
                    int lk = lookup->at(ind);
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
            fill(seed.begin(), seed.end(), 0);

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
    int nz_in = nnz_in(iind);

    // Number of nonzero outputs
    int nz_out = nnz_out(oind);

    // Seeds and sensitivities
    vector<bvec_t> s_in(nz_in, 0);
    vector<bvec_t> s_out(nz_out, 0);

    // Evaluation buffers
    vector<const bvec_t*> arg_fwd(sz_arg(), 0);
    vector<bvec_t*> arg_adj(sz_arg(), 0);
    arg_fwd[iind] = arg_adj[iind] = get_ptr(s_in);
    vector<bvec_t*> res(sz_res(), 0);
    res[oind] = get_ptr(s_out);
    vector<int> iw(sz_iw());
    vector<bvec_t> w(sz_w());

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
    double sp_w = adWeightSp();

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
      Sparsity D1 = rT.uni_coloring(r);
      // Adjoint mode
      Sparsity D2 = r.uni_coloring(rT);

      casadi_msg("Coloring on " << r.dim() << " (fwd seeps: " << D1.size2() <<
                 " , adj sweeps: " << D2.size1() << ")");

      // Use forward mode?
      int fwd_cost = use_fwd ? granularity_row: granularity_col;
      int adj_cost = use_fwd ? granularity_col: granularity_row;

      // Use whatever required less colors if we tried both (with preference to forward mode)
      if ((sp_w*D1.size2()*fwd_cost <= (1-sp_w)*D2.size2()*adj_cost)) {
        use_fwd = true;
        casadi_msg("Forward mode chosen (fwd cost: " << sp_w*D1.size2()*fwd_cost << ", adj cost: "
                   << (1-sp_w)*D2.size2()*adj_cost << ")");
      } else {
        use_fwd = false;
        casadi_msg("Adjoint mode chosen (adj cost: " << sp_w*D1.size2()*fwd_cost << ", adj cost: "
                   << (1-sp_w)*D2.size2()*adj_cost << ")");
      }

      // Get seeds and sensitivities
      bvec_t* seed_v = use_fwd ? get_ptr(s_in) : get_ptr(s_out);
      bvec_t* sens_v = use_fwd ? get_ptr(s_out) : get_ptr(s_in);

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
            for (int fci = fci_offset; fci<min(fci_end-fci_start, fci_cap); ++fci) {

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
            IM lookup = IM::triplet(lookup_row, lookup_col, lookup_value, bvec_size,
                                    coarse_col.size());

            // Propagate the dependencies
            if (use_fwd) {
              spFwd(get_ptr(arg_fwd), get_ptr(res), get_ptr(iw), get_ptr(w), 0);
            } else {
              fill(w.begin(), w.end(), 0);
              spAdj(get_ptr(arg_adj), get_ptr(res), get_ptr(iw), get_ptr(w), 0);
            }

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
                    int ind = lookup.sparsity().get_nz(bvec_i, cri);
                    casadi_assert(ind!=-1);
                    jrow.push_back(bvec_i+lookup->at(ind));
                    jcol.push_back(fri);
                  }
                }
              }
            }

            // Clear the forward seeds/adjoint sensitivities, ready for next bvec sweep
            fill(s_in.begin(), s_in.end(), 0);

            // Clear the adjoint seeds/forward sensitivities, ready for next bvec sweep
            fill(s_out.begin(), s_out.end(), 0);

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
      if (nnz_in(iind)>3*bvec_size && nnz_out(oind)>3*bvec_size) {
        if (symmetric) {
          sp = getJacSparsityHierarchicalSymm(iind, oind);
        } else {
          sp = getJacSparsityHierarchical(iind, oind);
        }
      } else {
        // Number of nonzero inputs and outputs
        int nz_in = nnz_in(iind);
        int nz_out = nnz_out(oind);

        // Number of forward sweeps we must make
        int nsweep_fwd = nz_in/bvec_size;
        if (nz_in%bvec_size) nsweep_fwd++;

        // Number of adjoint sweeps we must make
        int nsweep_adj = nz_out/bvec_size;
        if (nz_out%bvec_size) nsweep_adj++;

        // Get weighting factor
        double w = adWeightSp();

        // Use forward mode?
        if (w*nsweep_fwd <= (1-w)*nsweep_adj) {
          sp = getJacSparsityGen<true>(iind, oind, false);
        } else {
          sp = getJacSparsityGen<false>(iind, oind, false);
        }
      }
      // There may be false positives here that are not present
      // in the reverse mode that precedes it.
      // This can lead to an assymetrical result
      //  cf. #1522
      if (symmetric) sp=sp*sp.T();
      return sp;
    } else {
      // Dense sparsity by default
      return Sparsity::dense(nnz_out(oind), nnz_in(iind));
    }
  }

  void FunctionInternal::set_jac_sparsity(const Sparsity& sp, int iind, int oind, bool compact) {
    if (compact) {
      jac_sparsity_compact_.elem(oind, iind) = sp;
    } else {
      jac_sparsity_.elem(oind, iind) = sp;
    }
  }

  Sparsity& FunctionInternal::sparsity_jac(int iind, int oind, bool compact, bool symmetric) {
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
        Sparsity sp = sparsity_jac(iind, oind, true, symmetric);

        // Enlarge if sparse output
        if (numel_out(oind)!=sp.size1()) {
          casadi_assert(sp.size1()==nnz_out(oind));

          // New row for each old row
          vector<int> row_map = sparsity_out(oind).find();

          // Insert rows
          sp.enlargeRows(numel_out(oind), row_map);
        }

        // Enlarge if sparse input
        if (numel_in(iind)!=sp.size2()) {
          casadi_assert(sp.size2()==nnz_in(iind));

          // New column for each old column
          vector<int> col_map = sparsity_in(iind).find();

          // Insert columns
          sp.enlargeColumns(numel_in(iind), col_map);
        }

        // Save
        jsp = sp;
      }
    }

    // If still null, not dependent
    if (jsp.isNull()) {
      jsp = Sparsity(nnz_out(oind), nnz_in(iind));
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
    Sparsity &AT = sparsity_jac(iind, oind, compact, symmetric);
    Sparsity A = symmetric ? AT : AT.T();

    // Get seed matrices by graph coloring
    if (symmetric) {
      casadi_assert(get_n_forward()>0);

      // Star coloring if symmetric
      log("FunctionInternal::getPartition star_coloring");
      D1 = A.star_coloring();
      casadi_msg("Star coloring completed: " << D1.size2() << " directional derivatives needed ("
                 << A.size1() << " without coloring).");

    } else {
      casadi_assert(get_n_forward()>0 || get_n_reverse()>0);
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
          D1 = AT.uni_coloring(A, max_colorings_to_test);
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

          D2 = A.uni_coloring(AT, max_colorings_to_test);
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

  void FunctionInternal::eval(Memory& mem,
                              const double** arg, double** res, int* iw, double* w) const {
    casadi_error("'eval' not defined for " + type_name());
  }

  void FunctionInternal::eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    eval(*mem_.at(mem), arg, res, iw, w);
  }

  void FunctionInternal::simple(const double* arg, double* res) {
    casadi_error("'simple' not defined for " + type_name());
  }

  void FunctionInternal::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    casadi_error("'eval_sx' not defined for " + type_name());
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
      Dict opts = {{"verbose", verbose_},
                   {"input_scheme", ischeme_},
                   {"output_scheme", ionames},
                   {"jit", jit_},
                   {"compiler", compilerplugin_},
                   {"jit_options", jit_options_}};
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

  Function FunctionInternal::forward(int nfwd) {
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
    Dict opts = {{"input_scheme", i_names},
                 {"output_scheme", o_names},
                 {"jit", jit_},
                 {"compiler", compilerplugin_},
                 {"jit_options", jit_options_}};

    // Return value
    casadi_assert(get_n_forward()>0);
    Function ret = get_forward(name, nfwd, opts);

    // Consistency check for inputs
    for (int i=0; i<ret.n_in(); ++i) {
      const Sparsity& sp = i<n_in ? sparsity_in(i) :
        i<n_in+n_out ? sparsity_out(i-n_in) :
        sparsity_in((i-n_in-n_out) % n_in);
      casadi_assert_message(ret.size_in(i)==sp.size(),
                            "Incorrect shape for " << ret << " input " << i << " \""
                            << i_names.at(i) << "\". Expected " << sp.size()
                            << " but got " << ret.size_in(i));
    }

    // Consistency check for outputs
    for (int i=0; i<ret.n_out(); ++i) {
      const Sparsity& sp = sparsity_out(i % n_out);
      casadi_assert_message(ret.size_out(i)==sp.size(),
                            "Incorrect shape for " << ret << " output " << i << " \""
                            << o_names.at(i) << "\". Expected " << sp.size()
                            << " but got " << ret.size_out(i));
    }

    // Save to cache
    derivative_fwd_[nfwd] = ret;

    // Return generated function
    return ret;
  }

  Function FunctionInternal::reverse(int nadj) {
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
    Dict opts = {{"input_scheme", i_names},
                 {"output_scheme", o_names},
                 {"jit", jit_},
                 {"compiler", compilerplugin_},
                 {"jit_options", jit_options_}};

    // Return value
    casadi_assert(get_n_reverse()>0);
    Function ret = get_reverse(name, nadj, opts);

    // Consistency check for inputs
    for (int i=0; i<ret.n_in(); ++i) {
      const Sparsity& sp = i<n_in ? sparsity_in(i) :
        i<n_in+n_out ? sparsity_out(i-n_in) :
        sparsity_out((i-n_in-n_out) % n_out);
      casadi_assert_message(ret.size_in(i)==sp.size(),
                            "Incorrect shape for " << ret << " input " << i << " \""
                            << i_names.at(i) << "\". Expected " << sp.size()
                            << " but got " << ret.size_in(i));
    }

    // Consistency check for outputs
    for (int i=0; i<ret.n_out(); ++i) {
      const Sparsity& sp = sparsity_in(i % n_in);
      casadi_assert_message(ret.size_out(i)==sp.size(),
                            "Incorrect shape for " << ret << " output " << i << " \""
                            << o_names.at(i) << "\". Expected " << sp.size()
                            << " but got " << ret.size_out(i));
    }

    // Save to cache
    derivative_adj_[nadj] = ret;

    // Return generated function
    return ret;
  }

  void FunctionInternal::set_forward(const Function& fcn, int nfwd) {

    // Check if there are enough forward directions allocated
    if (nfwd>=derivative_fwd_.size()) {
      derivative_fwd_.resize(nfwd+1);
    }

    // Save to cache
    derivative_fwd_[nfwd] = fcn;
  }

  void FunctionInternal::set_reverse(const Function& fcn, int nadj) {

    // Check if there are enough adjoint directions allocated
    if (nadj>=derivative_adj_.size()) {
      derivative_adj_.resize(nadj+1);
    }

    // Save to cache
    derivative_adj_[nadj] = fcn;
  }

  Function FunctionInternal::get_forward(const std::string& name, int nfwd, Dict& opts) {
    // TODO(@jaeandersson): Fallback on finite differences
    casadi_error("'get_forward' not defined for " + type_name());
  }

  Function FunctionInternal::get_reverse(const std::string& name, int nadj, Dict& opts) {
    casadi_error("'get_reverse' not defined for " + type_name());
  }

  int FunctionInternal::nnz_in() const {
    int ret=0;
    for (int iind=0; iind<n_in(); ++iind) ret += nnz_in(iind);
    return ret;
  }

  int FunctionInternal::nnz_out() const {
    int ret=0;
    for (int oind=0; oind<n_out(); ++oind) ret += nnz_out(oind);
    return ret;
  }

  int FunctionInternal::numel_in() const {
    int ret=0;
    for (int iind=0; iind<n_in(); ++iind) ret += numel_in(iind);
    return ret;
  }

  int FunctionInternal::numel_out() const {
    int ret=0;
    for (int oind=0; oind<n_out(); ++oind) ret += numel_out(oind);
    return ret;
  }

  void FunctionInternal::eval_mx(const MXVector& arg, MXVector& res,
                                 bool always_inline, bool never_inline) {
    // The code below creates a call node, to inline, wrap in an MXFunction
    if (always_inline) {
      casadi_assert_message(!never_inline, "Inconsistent options");
      Function f = wrapMXFunction();
      return f.call(arg, res, true);
    }

    // Create a call-node
    res = create_call(arg);
  }

  std::vector<MX> FunctionInternal::create_call(const std::vector<MX>& arg) {
    return Call::create(shared_from_this<Function>(), arg);
  }

  Function FunctionInternal::
  getNumericJacobian(const std::string& name, int iind, int oind, bool compact, bool symmetric,
                     const Dict& opts) {
    Function f = wrapMXFunction();
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
    casadi_assert(get_n_forward()>0 || get_n_reverse()>0);

    // Number inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Symbolic inputs of the full Jacobian function under construction
    vector<MX> ret_argv = mx_in(), argv, resv;

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
        sp_argv[i] = vec(sparsity_in(i));
        row_offset[i+1] = row_offset[i]+sp_argv[i].numel();
      }
      Sparsity sp_arg = vertcat(sp_argv);

      // Create symbolic primitive
      arg = MX::sym("x", sp_arg);

      // Split up and reshape to correct shape
      argv = vertsplit(arg, row_offset);
      for (int i=0; i<n_in; ++i) {
        argv[i] = reshape(argv[i], sparsity_in(i));
      }

      // Evaluate symbolically
      resv = shared_from_this<Function>()(argv);
    }

    // Reuse the same output, if possible
    MX res = n_out==1 ? resv.front() : veccat(resv);

    // Form Jacobian
    MX J;
    {
      Function tmp("tmp", {arg}, {res}, {{"ad_weight", adWeight()}});
      J = MX::jac(tmp);
    }

    // Make sure argv is the input of J
    if (n_in!=1) {
      J = substitute(J, arg, veccat(ret_argv));
    }

    // Form an expression for the full Jacobian
    return Function(name, ret_argv, {J}, opts);
  }

  void FunctionInternal::generateFunction(CodeGenerator& g,
                                          const std::string& fname, bool decl_static) const {
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
    g.body << signature(fname) << " {" << endl;

    // Insert the function body
    generateBody(g);

    // Finalize the function
    if (!simplifiedCall()) g.body << "  return 0;" << endl;
    g.body << "}" << endl << endl;
  }

  void FunctionInternal::generate(CodeGenerator& g, const std::string& mem,
                                  const vector<int>& arg, const vector<int>& res) const {
    if (simplifiedCall()) {

      // Collect input arguments
      for (int i=0; i<arg.size(); ++i) {
        g.body << "  w[" << i << "]=" << g.workel(arg[i]) << ";" << endl;
      }

      // Call function
      g.body << "  " << simple_call(g, "w", "w+"+g.to_string(arg.size())) << ";" << endl;

      // Collect output arguments
      for (int i=0; i<res.size(); ++i) {
        if (res[i]>=0) {
          g.body << "  " << g.workel(res[i]) << "=w[" << (arg.size()+i) << "];" << endl;
        }
      }
    } else {

      // Collect input arguments
      for (int i=0; i<arg.size(); ++i) {
        g.body << "  arg1[" << i << "]=" << g.work(arg[i], nnz_in(i)) << ";" << endl;
      }

      // Collect output arguments
      for (int i=0; i<res.size(); ++i) {
        g.body << "  res1[" << i << "]=" << g.work(res[i], nnz_out(i)) << ";" << endl;
      }

      // Call function
      g.body << "  if (" << generic_call(g, "arg1", "res1", "iw", "w", "0")
             << ") return 1;" << endl;
    }
  }

  std::string FunctionInternal::signature(const std::string& fname) const {
    if (simplifiedCall()) {
      return "void " + fname + "(const real_t* arg, real_t* res)";
    } else {
      return "int " + fname + "(const real_t** arg, real_t** res, int* iw, real_t* w, int mem)";
    }
  }

  void FunctionInternal::generateMeta(CodeGenerator& g, const std::string& fname) const {
    // Short-hands
    int n_in = this->n_in();
    int n_out = this->n_out();
    stringstream &s = g.body;

    // Function that returns the number of inputs and outputs
    string tmp = "int " + fname + "_init(int* n_in, int* n_out, int* n_mem)";
    if (g.cpp) {
      tmp = "extern \"C\" " + tmp;  // C linkage
    }
    if (g.with_header) {
      g.header << tmp << ";" << endl;
    }
    s << tmp << " {" << endl
      << "  if (n_in) *n_in = " << n_in << ";" << endl
      << "  if (n_out) *n_out = " << n_out << ";" << endl
      << "  if (n_mem) *n_mem = -1;" << endl
      << "  return 0;" << endl
      << "}" << endl
      << endl;

    // Quick return if simplified syntax
    if (simplifiedCall()) {
      return;
    }

    // Function for allocating memory
    tmp = "int " + fname + "_alloc(int mem)";
    if (g.cpp) {
      tmp = "extern \"C\" " + tmp;  // C linkage
    }
    if (g.with_header) {
      g.header << tmp << ";" << endl;
    }
    s << tmp << " {" << endl
      << "  (void)mem;" << endl
      << "  return 0;" << endl
      << "}" << endl
      << endl;

    // Function for freeing memory
    tmp = "int " + fname + "_free(int mem)";
    if (g.cpp) {
      tmp = "extern \"C\" " + tmp;  // C linkage
    }
    if (g.with_header) {
      g.header << tmp << ";" << endl;
    }
    s << tmp << " {" << endl
      << "  (void)mem;" << endl
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
      const Sparsity& sp = i<n_in ? sparsity_in(i) : sparsity_out(i-n_in);

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
      auto r = io_sparsity_index.equal_range(i);

      // Skip if unused
      if (r.first==r.second) continue;

      // Print the cases covered
      for (auto it=r.first; it!=r.second; ++it) {
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
    s << "  if (nrow) *nrow = s[0];" << endl;
    s << "  if (ncol) *ncol = s[1];" << endl;
    s << "  if (colind) *colind = s + 2;" << endl;
    s << "  if (row) *row = s + 3 + s[1];" << endl;
    s << "  return 0;" << endl;
    s << "}" << endl << endl;

    // Function that returns work vector lengths
    tmp = "int " + fname
      + "_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w)";
    if (g.cpp) {
      tmp = "extern \"C\" " + tmp;  // C linkage
    }
    if (g.with_header) {
      g.header << tmp << ";" << endl;
    }
    s << tmp << " {" << endl;
    s << "  if (sz_arg) *sz_arg = " << sz_arg() << ";" << endl;
    s << "  if (sz_res) *sz_res = " << sz_res() << ";" << endl;
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
        const Sparsity& s = sparsity_in(i);
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
          << g.from_mex(p, "w", offset, sparsity_in(i), fw) << endl;
        offset += nnz_in(i);
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
        offset += nnz_out(i);
      }

      // Call the function
      s << "  i = " << fname << "(arg, res, iw, " << fw << ", 0);" << endl;
      s << "  if (i) mexErrMsgIdAndTxt(\"Casadi:RuntimeError\",\"Evaluation of \\\"" << fname
        << "\\\" failed.\");" << endl;

      // Save results
      for (int i=0; i<n_out; ++i) {
        string res_i = "res[" + g.to_string(i) + "]";
        s << "  if (" << res_i << ") resv[" << i << "] = "
          << g.to_mex(sparsity_out(i), res_i) << endl;
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
        s << "w+" << off;
        off += nnz_in(i);
      }
      s << "};" << endl;

      // Output buffers
      s << "  real_t* res[" << n_out << "] = {";
      for (int i=0; i<n_out; ++i) {
        if (i!=0) s << ", ";
        s << "w+" << off;
        off += nnz_out(i);
      }
      s << "};" << endl;

      // TODO(@jaeandersson): Read inputs from file. For now; read from stdin
      s << "  int j;" << endl
        << "  real_t* a = w;" << endl
        << "  for (j=0; j<" << nnz_in() << "; ++j) "
        << "scanf(\"%lf\", a++);" << endl;

      // Call the function
      s << "  int flag = " << fname << "(arg, res, iw, w+" << off << ", 0);" << endl
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

  std::string FunctionInternal::generic_call(const CodeGenerator& g, const std::string& arg,
                                             const std::string& res, const std::string& iw,
                                             const std::string& w, const std::string& mem) const {
    // Get the index of the function
    auto it=g.added_dependencies_.find(this);
    casadi_assert(it!=g.added_dependencies_.end());
    int f = it->second;

    // Create a function call
    stringstream ss;
    ss << "f" << f << "(" << arg << ", " << res << ", " << iw << ", " << w << " , " << mem << ")";
    return ss.str();
  }

  std::string FunctionInternal::simple_call(const CodeGenerator& g,
                                             const std::string& arg,
                                             const std::string& res) const {
    // Get the index of the function
    auto it=g.added_dependencies_.find(this);
    casadi_assert(it!=g.added_dependencies_.end());
    int f = it->second;

    // Create a function call
    stringstream ss;
    ss << "f" << f << "(" << arg << ", " << res << ")";
    return ss.str();
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
      if (simplifiedCall()) {
        g.body
          << "#define " << name << "(arg, res) "
          << "CASADI_PREFIX(" << name << ")(arg, res)" << endl << endl;
      } else {
        g.body
          << "#define " << name << "(arg, res, iw, w, mem) "
          << "CASADI_PREFIX(" << name << ")(arg, res, iw, w, mem)" << endl << endl;
      }
    }
  }

  void FunctionInternal::generateDeclarations(CodeGenerator& g) const {
    // Nothing to declare
  }

  void FunctionInternal::generateBody(CodeGenerator& g) const {
    casadi_error("'generateBody' not defined for " + type_name());
  }

  Function FunctionInternal::dynamicCompilation(Function f, std::string fname, std::string fdescr,
                                                std::string compiler) {
    // Codegen and compile
    CodeGenerator g;
    g.add(f, fname);
    string dlname = g.compile(fname, compiler);

    // Load it
    return external(fname, dlname);
  }

  void FunctionInternal::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    // Get the number of inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Loop over outputs
    for (int oind=0; oind<n_out; ++oind) {
      // Skip if nothing to assign
      if (res[oind]==0 || nnz_out(oind)==0) continue;

      // Clear result
      casadi_fill(res[oind], nnz_out(oind), bvec_t(0));

      // Loop over inputs
      for (int iind=0; iind<n_in; ++iind) {
        // Skip if no seeds
        if (arg[iind]==0 || nnz_in(iind)==0) continue;

        // Get the sparsity of the Jacobian block
        Sparsity sp = sparsity_jac(iind, oind, true, false);
        if (sp.isNull() || sp.nnz() == 0) continue; // Skip if zero

        // Carry out the sparse matrix-vector multiplication
        int d1 = sp.size2();
        const int *colind = sp.colind(), *row = sp.row();
        for (int cc=0; cc<d1; ++cc) {
          for (int el = colind[cc]; el < colind[cc+1]; ++el) {
            res[oind][row[el]] |= arg[iind][cc];
          }
        }
      }
    }
  }

  void FunctionInternal::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    // Get the number of inputs and outputs
    int n_in = this->n_in();
    int n_out = this->n_out();

    // Loop over outputs
    for (int oind=0; oind<n_out; ++oind) {
      // Skip if nothing to assign
      if (res[oind]==0 || nnz_out(oind)==0) continue;

      // Loop over inputs
      for (int iind=0; iind<n_in; ++iind) {
        // Skip if no seeds
        if (arg[iind]==0 || nnz_in(iind)==0) continue;

        // Get the sparsity of the Jacobian block
        Sparsity sp = sparsity_jac(iind, oind, true, false);
        if (sp.isNull() || sp.nnz() == 0) continue; // Skip if zero

        // Carry out the sparse matrix-vector multiplication
        int d1 = sp.size2();
        const int *colind = sp.colind(), *row = sp.row();
        for (int cc=0; cc<d1; ++cc) {
          for (int el = colind[cc]; el < colind[cc+1]; ++el) {
            arg[iind][cc] |= res[oind][row[el]];
          }
        }
      }

      // Clear seeds
      casadi_fill(res[oind], nnz_out(oind), bvec_t(0));
    }
  }

  void FunctionInternal::sz_work(size_t& sz_arg, size_t& sz_res,
                                 size_t& sz_iw, size_t& sz_w) const {
    sz_arg = this->sz_arg();
    sz_res = this->sz_res();
    sz_iw = this->sz_iw();
    sz_w = this->sz_w();
  }

  void FunctionInternal::alloc_arg(size_t sz_arg, bool persistent) {
    if (persistent) {
      sz_arg_per_ += sz_arg;
    } else {
      sz_arg_tmp_ = max(sz_arg_tmp_, sz_arg);
    }
  }

  void FunctionInternal::alloc_res(size_t sz_res, bool persistent) {
    if (persistent) {
      sz_res_per_ += sz_res;
    } else {
      sz_res_tmp_ = max(sz_res_tmp_, sz_res);
    }
  }

  void FunctionInternal::alloc_iw(size_t sz_iw, bool persistent) {
    if (persistent) {
      sz_iw_per_ += sz_iw;
    } else {
      sz_iw_tmp_ = max(sz_iw_tmp_, sz_iw);
    }
  }

  void FunctionInternal::alloc_w(size_t sz_w, bool persistent) {
    if (persistent) {
      sz_w_per_ += sz_w;
    } else {
      sz_w_tmp_ = max(sz_w_tmp_, sz_w);
    }
  }

  void FunctionInternal::alloc(const Function& f, bool persistent) {
    if (f.isNull()) return;
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f.sz_work(sz_arg, sz_res, sz_iw, sz_w);
    alloc_arg(sz_arg, persistent);
    alloc_res(sz_res, persistent);
    alloc_iw(sz_iw, persistent);
    alloc_w(sz_w, persistent);
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
    return get_n_forward()>0 || get_n_reverse()>0 || hasFullJacobian();
  }

  bool FunctionInternal::fwdViaJac(int nfwd) {
    if (get_n_forward()==0) return true;

    // Jacobian calculation penalty factor
    const double jac_penalty = option("jac_penalty");

    if (jac_penalty==-1) return false;

    // Heuristic 1: Jac calculated via forward mode likely cheaper
    if (jac_penalty*nnz_in()<nfwd) return true;

    // Heuristic 2: Jac calculated via reverse mode likely cheaper
    double w = adWeight();
    if (get_n_reverse()>0 && jac_penalty*(1-w)*nnz_out()<w*nfwd)
      return true;

    return false;
  }

  bool FunctionInternal::adjViaJac(int nadj) {
    if (get_n_reverse()==0) return true;

    // Jacobian calculation penalty factor
    const double jac_penalty = option("jac_penalty");

    if (jac_penalty==-1) return false;

    // Heuristic 1: Jac calculated via reverse mode likely cheaper
    if (jac_penalty*nnz_out()<nadj) return true;

    // Heuristic 2: Jac calculated via forward mode likely cheaper
    double w = adWeight();
    if (get_n_forward()>0 && jac_penalty*w*nnz_in()<(1-w)*nadj)
      return true;

    return false;
  }

  void FunctionInternal::
  forward_mx(const std::vector<MX>& arg, const std::vector<MX>& res,
             const std::vector<std::vector<MX> >& fseed,
             std::vector<std::vector<MX> >& fsens,
             bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    casadi_assert_message(!always_inline, "Class " + type_name() +
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
      v = horzsplit(mtimes(J, horzcat(v)));

      // Vertical offsets
      vector<int> offset(n_out+1, 0);
      for (int i=0; i<n_out; ++i) {
        offset[i+1] = offset[i]+numel_out(i);
      }

      // Collect forward sensitivities
      for (int d=0; d<nfwd; ++d) {
        fsens[d] = vertsplit(v[d], offset);
        for (int i=0; i<n_out; ++i) {
          fsens[d][i] = reshape(fsens[d][i], size_out(i));
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
    Function dfcn = forward(nfwd);

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

  void FunctionInternal::
  reverse_mx(const std::vector<MX>& arg, const std::vector<MX>& res,
             const std::vector<std::vector<MX> >& aseed,
             std::vector<std::vector<MX> >& asens,
             bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    casadi_assert_message(!always_inline, "Class " + type_name() +
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
      v = horzsplit(mtimes(J.T(), horzcat(v)));

      // Vertical offsets
      vector<int> offset(n_in+1, 0);
      for (int i=0; i<n_in; ++i) {
        offset[i+1] = offset[i]+numel_in(i);
      }

      // Collect adjoint sensitivities
      for (int d=0; d<nadj; ++d) {
        asens[d].resize(n_in);
        vector<MX> a = vertsplit(v[d], offset);
        for (int i=0; i<n_in; ++i) {
          if (asens[d][i].is_empty(true)) {
            asens[d][i] = reshape(a[i], size_in(i));
          } else {
            asens[d][i] += reshape(a[i], size_in(i));
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
    Function dfcn = reverse(nadj);

    // Create the evaluation node
    vector<MX> x = Call::create(dfcn, darg);
    vector<MX>::iterator x_it = x.begin();

    // Retrieve sensitivities
    for (int d=0; d<nadj; ++d) {
      asens[d].resize(n_in);
      for (int i=0; i<n_in; ++i) {
        if (asens[d][i].is_empty(true)) {
          asens[d][i] = *x_it++;
        } else {
          asens[d][i] += *x_it++;
        }
      }
    }
    casadi_assert(x_it==x.end());
  }

  void FunctionInternal::
  forward_sx(const std::vector<SX>& arg, const std::vector<SX>& res,
             const std::vector<std::vector<SX> >& fseed,
             std::vector<std::vector<SX> >& fsens,
             bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    if (fseed.empty()) { // Quick return if no seeds
      fsens.clear();
      return;
    }
    casadi_error("'forward' (SX) not defined for " + type_name());
  }

  void FunctionInternal::
  reverse_sx(const std::vector<SX>& arg, const std::vector<SX>& res,
             const std::vector<std::vector<SX> >& aseed,
             std::vector<std::vector<SX> >& asens,
             bool always_inline, bool never_inline) {
    casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
    if (aseed.empty()) { // Quick return if no seeds
      asens.clear();
      return;
    }
    casadi_error("'reverse' (SX) not defined for " + type_name());
  }

  void FunctionInternal::
  forward_dm(const std::vector<DM>& arg, const std::vector<DM>& res,
          const std::vector<std::vector<DM> >& fseed,
          std::vector<std::vector<DM> >& fsens,
          bool always_inline, bool never_inline) {
    casadi_error("Not implemented");
  }

  void FunctionInternal::
  reverse_dm(const std::vector<DM>& arg, const std::vector<DM>& res,
          const std::vector<std::vector<DM> >& aseed,
          std::vector<std::vector<DM> >& asens,
          bool always_inline, bool never_inline) {
    casadi_error("Not implemented");
  }

  double FunctionInternal::adWeight() {
    // If reverse mode derivatives unavailable, use forward
    if (get_n_reverse()==0) return 0;

    // If forward mode derivatives unavailable, use reverse
    if (get_n_forward()==0) return 1;

    // A user-provided option overrides default value
    if (hasSetOption("ad_weight")) return option("ad_weight");

    // By default, reverse mode is about twice as expensive as forward mode
    return 0.33;  // i.e. nf <= 2*na <=> 1/3*nf <= (1-1/3)*na, forward when tie
  }

  double FunctionInternal::adWeightSp() {
    // If reverse mode propagation unavailable, use forward
    if (!spCanEvaluate(false)) return 0;

    // If forward mode propagation unavailable, use reverse
    if (!spCanEvaluate(true)) return 1;

    // A user-provided option overrides default value
    if (hasSetOption("ad_weight_sp")) return option("ad_weight_sp");

    // Both modes equally expensive by default (no "taping" needed)
    return 0.49; // Forward when tie
  }

  MX FunctionInternal::grad_mx(int iind, int oind) {
    casadi_error("'grad_mx' not defined for " + type_name());
  }

  MX FunctionInternal::tang_mx(int iind, int oind) {
    casadi_error("'tang_mx' not defined for " + type_name());
  }

  MX FunctionInternal::jac_mx(int iind, int oind, bool compact, bool symmetric,
                              bool always_inline, bool never_inline) {
    casadi_error("'jac_mx' not defined for " + type_name());
  }

  SX FunctionInternal::grad_sx(int iind, int oind) {
    casadi_error("'grad_sx' not defined for " + type_name());
  }

  SX FunctionInternal::tang_sx(int iind, int oind) {
    casadi_error("'tang_sx' not defined for " + type_name());
  }

  SX FunctionInternal::jac_sx(int iind, int oind, bool compact, bool symmetric,
                              bool always_inline, bool never_inline) {
    casadi_error("'jac_sx' not defined for " + type_name());
  }

  SX FunctionInternal::hess_sx(int iind, int oind) {
    casadi_error("'hess_sx' not defined for " + type_name());
  }

  const SX FunctionInternal::sx_in(int ind) const {
    return SX::sym("x_" + CodeGenerator::to_string(ind), sparsity_in(ind));
  }

  const SX FunctionInternal::sx_out(int ind) const {
    return SX::sym("r_" + CodeGenerator::to_string(ind), sparsity_out(ind));
  }

  const std::vector<SX> FunctionInternal::sx_in() const {
    vector<SX> ret(n_in());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = sx_in(i);
    }
    return ret;
  }

  const std::vector<SX> FunctionInternal::sx_out() const {
    vector<SX> ret(n_out());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = sx_out(i);
    }
    return ret;
  }

  const MX FunctionInternal::mx_in(int ind) const {
    return MX::sym("x_" + CodeGenerator::to_string(ind), sparsity_in(ind));
  }

  const MX FunctionInternal::mx_out(int ind) const {
    return MX::sym("r_" + CodeGenerator::to_string(ind), sparsity_out(ind));
  }

  const std::vector<MX> FunctionInternal::mx_in() const {
    vector<MX> ret(n_in());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = mx_in(i);
    }
    return ret;
  }

  const std::vector<MX> FunctionInternal::mx_out() const {
    vector<MX> ret(n_out());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = mx_out(i);
    }
    return ret;
  }

  std::string FunctionInternal::type_name() const {
    return typeid(*this).name();
  }

  bool FunctionInternal::is_a(const std::string& type, bool recursive) const {
    return type == "function";
  }

  std::vector<MX> FunctionInternal::free_mx() const {
    casadi_error("'free_mx' only defined for 'mxfunction'");
  }

  SX FunctionInternal::free_sx() const {
    casadi_error("'free_sx' only defined for 'sxfunction'");
  }

  void FunctionInternal::generate_lifted(Function& vdef_fcn,
                                                  Function& vinit_fcn) {
    casadi_error("'generate_lifted' only defined for 'mxfunction'");
  }

  int FunctionInternal::getAlgorithmSize() const {
    casadi_error("'getAlgorithmSize' not defined for " + type_name());
  }

  int FunctionInternal::getWorkSize() const {
    casadi_error("'getWorkSize' not defined for " + type_name());
  }

  int FunctionInternal::getAtomicOperation(int k) const {
    casadi_error("'getAtomicOperation' not defined for " + type_name());
  }

  std::pair<int, int> FunctionInternal::getAtomicInput(int k) const {
    casadi_error("'getAtomicInput' not defined for " + type_name());
  }

  double FunctionInternal::getAtomicInputReal(int k) const {
    casadi_error("'getAtomicInputReal' not defined for " + type_name());
  }

  int FunctionInternal::getAtomicOutput(int k) const {
    casadi_error("'getAtomicOutput' not defined for " + type_name());
  }

  int FunctionInternal::countNodes() const {
    casadi_error("'countNodes' not defined for " + type_name());
  }

  void FunctionInternal::linsol_factorize(Memory& mem, const double* A) const {
    casadi_error("'linsol_factorize' not defined for " + type_name());
  }

  void FunctionInternal::linsol_solve(Memory& mem, double* x, int nrhs, bool tr) const {
    casadi_error("'linsol_solve' not defined for " + type_name());
  }

  MX FunctionInternal::linsol_solve(const MX& A, const MX& B, bool tr) {
    casadi_error("'linsol_solve' not defined for " + type_name());
  }

  void FunctionInternal::linsol_spsolve(bvec_t* X, const bvec_t* B, bool tr) const {
    casadi_error("'linsol_spsolve' not defined for " + type_name());
  }

  void FunctionInternal::linsol_spsolve(DM& X, const DM& B, bool tr) const {
    casadi_error("'linsol_spsolve' not defined for " + type_name());
  }

  void FunctionInternal::linsol_solveL(Memory& mem, double* x, int nrhs, bool tr) const {
    casadi_error("'linsol_solveL' not defined for " + type_name());
  }

  Sparsity FunctionInternal::linsol_cholesky_sparsity(Memory& mem, bool tr) const {
    casadi_error("'linsol_cholesky_sparsity' not defined for " + type_name());
  }

  DM FunctionInternal::linsol_cholesky(Memory& mem, bool tr) const {
    casadi_error("'linsol_cholesky' not defined for " + type_name());
  }

  void FunctionInternal::
  linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem,
                bool tr, int nrhs) {
    casadi_error("'linsol_eval_sx' not defined for " + type_name());
  }

  void FunctionInternal::
  linsol_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& fseed,
                 std::vector<std::vector<MX> >& fsens, bool tr) {
    casadi_error("'linsol_forward' not defined for " + type_name());
  }

  void FunctionInternal::
  linsol_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& aseed,
                 std::vector<std::vector<MX> >& asens, bool tr) {
    casadi_error("'linsol_reverse' not defined for " + type_name());
  }

  void FunctionInternal::
  linsol_spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem,
               bool tr, int nrhs) {
    casadi_error("'linsol_spFwd' not defined for " + type_name());
  }

  void FunctionInternal::
  linsol_spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem,
               bool tr, int nrhs) {
    casadi_error("'linsol_spAdj' not defined for " + type_name());
  }

  std::vector<std::vector<MX> >
  FunctionInternal::map_mx(const std::vector<std::vector<MX> > &x,
                           const std::string& parallelization) {
    if (x.empty()) return x;

    // Check if arguments match
    int n = x.size();
    bool matching=true;
    for (int i=0; i<n; ++i) {
      matching = matchingArg(x[i]) && matching; // non-short circuiting
    }

    // Replace arguments if needed
    if (!matching) {
      vector<vector<MX> > x_new(n);
      for (int i=0; i<n; ++i) x_new[i] = replaceArg(x[i]);
      return map_mx(x_new, parallelization);
    }

    vector< vector<MX> > trans = swapIndices(x);
    vector< MX > x_cat(trans.size());
    for (int i=0;i<trans.size();++i) {
      x_cat[i] = horzcat(trans[i]);
    }

    // Call the internal function
    vector< MX > ret_cat = map_mx(x_cat, parallelization);
    vector< vector<MX> > ret;

    for (int i=0;i<ret_cat.size();++i) {
      ret.push_back(horzsplit(ret_cat[i], size2_out(i)));
    }
    return swapIndices(ret);
  }

  std::vector<MX>
  FunctionInternal::map_mx(const std::vector<MX > &x,
                           const std::string& parallelization) {
    if (x.empty()) return x;

    // Replace arguments if needed
    if (!matchingArg(x, true)) {
      vector< MX > x_new = replaceArg(x, true);
      return map_mx(x_new, parallelization);
    }

    int n = 1;
    for (int i=0;i<x.size();++i) {
      n = max(x[i].size2()/size2_in(i), n);
    }

    vector<bool> repeat_n;
    for (int i=0;i<x.size();++i) {
      repeat_n.push_back(x[i].size2()/size2_in(i)==n);
    }

    bool repeated = true;
    for (int i=0;i<repeat_n.size();++i) {
      repeated &= repeat_n[i];
    }

    // Call the internal function

    Dict options = {{"parallelization", parallelization}};

    Function ms;
    if (repeated) {
      ms = shared_from_this<Function>().map("map", n, options);
    } else {
      ms = shared_from_this<Function>().
        map("mapsum", n, repeat_n, vector<bool>(n_out(), true), options);
    }
    // Call the internal function
    return ms(x);
  }

  std::vector<MX>
  FunctionInternal::mapsum_mx(const std::vector<MX > &x,
                              const std::string& parallelization) {
    if (x.empty()) return x;

    // Replace arguments if needed
    if (!matchingArg(x, true)) {
      vector< MX > x_new = replaceArg(x, true);
      return mapsum_mx(x_new, parallelization);
    }

    int n = 1;
    for (int i=0;i<x.size();++i) {
      n = max(x[i].size2()/size2_in(i), n);
    }

    vector<bool> repeat_n;
    for (int i=0;i<x.size();++i) {
      repeat_n.push_back(x[i].size2()/size2_in(i)==n);
    }

    Dict options = {{"parallelization", parallelization}};
    Function ms = shared_from_this<Function>().map("mapsum", n, repeat_n,
                      vector<bool>(n_out(), false), options);

    // Call the internal function
    return ms(x);
  }

  bool FunctionInternal::checkMat(const Sparsity& arg, const Sparsity& inp, bool hcat) {
    return arg.size()==inp.size() || arg.is_empty() || arg.is_scalar() ||
      (inp.size2()==arg.size1() && inp.size1()==arg.size2()
       && (arg.is_column() || inp.is_column())) ||
      (hcat && arg.size1()==inp.size1() && arg.size2() % inp.size2()==0);
  }

  void FunctionInternal::setup(Memory& mem, const double** arg, double** res,
                               int* iw, double* w) const {
    set_work(mem, arg, res, iw, w);
    set_temp(mem, arg, res, iw, w);
  }

} // namespace casadi
