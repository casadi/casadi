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


#include "sx_function.hpp"
#include <limits>
#include <stack>
#include <deque>
#include <sstream>
#include <iomanip>
#include "sx_node.hpp"
#include "casadi_common.hpp"
#include "sparsity_internal.hpp"
#include "casadi_interrupt.hpp"
#include "serializing_stream.hpp"

namespace casadi {

  SXFunction::SXFunction(const std::string& name,
                         const std::vector<SX >& inputv,
                         const std::vector<SX >& outputv,
                         const std::vector<std::string>& name_in,
                         const std::vector<std::string>& name_out)
    : XFunction<SXFunction, SX, SXNode>(name, inputv, outputv, name_in, name_out) {

    // Default (persistent) options
    just_in_time_opencl_ = false;
    just_in_time_sparsity_ = false;
  }

  SXFunction::~SXFunction() {
    clear_mem();
  }

  int SXFunction::eval(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const {
    if (verbose_) casadi_message(name_ + "::eval");

    // Make sure no free parameters
    if (!free_vars_.empty()) {
      std::stringstream ss;
      disp(ss, false);
      casadi_error("Cannot evaluate \"" + ss.str() + "\" since variables "
                   + str(free_vars_) + " are free.");
    }

    // NOTE: The implementation of this function is very delicate. Small changes in the
    // class structure can cause large performance losses. For this reason,
    // the preprocessor macros are used below

    // Evaluate the algorithm
    for (auto&& e : algorithm_) {
      switch (e.op) {
        CASADI_MATH_FUN_BUILTIN(w[e.i1], w[e.i2], w[e.i0])

      case OP_CONST: w[e.i0] = e.d; break;
      case OP_INPUT: w[e.i0] = arg[e.i1]==nullptr ? 0 : arg[e.i1][e.i2]; break;
      case OP_OUTPUT: if (res[e.i0]!=nullptr) res[e.i0][e.i2] = w[e.i1]; break;
      default:
        casadi_error("Unknown operation" + str(e.op));
      }
    }
    return 0;
  }

  bool SXFunction::is_smooth() const {
    // Go through all nodes and check if any node is non-smooth
    for (auto&& a : algorithm_) {
      if (!operation_checker<SmoothChecker>(a.op)) {
        return false;
      }
    }
    return true;
  }

  void SXFunction::disp_more(std::ostream &stream) const {
    stream << "Algorithm:";

    // Iterator to free variables
    std::vector<SXElem>::const_iterator p_it = free_vars_.begin();

    // Normal, interpreted output
    for (auto&& a : algorithm_) {
      InterruptHandler::check();
      stream << std::endl;
      if (a.op==OP_OUTPUT) {
        stream << "output[" << a.i0 << "][" << a.i2 << "] = @" << a.i1;
      } else {
        stream << "@" << a.i0 << " = ";
        if (a.op==OP_INPUT) {
          stream << "input[" << a.i1 << "][" << a.i2 << "]";
        } else {
          if (a.op==OP_CONST) {
            stream << a.d;
          } else if (a.op==OP_PARAMETER) {
            stream << *p_it++;
          } else {
            casadi_int ndep = casadi_math<double>::ndeps(a.op);
            stream << casadi_math<double>::pre(a.op);
            for (casadi_int c=0; c<ndep; ++c) {
              if (c==0) {
                stream << "@" << a.i1;
              } else {
                stream << casadi_math<double>::sep(a.op);
                stream << "@" << a.i2;
              }

            }
            stream << casadi_math<double>::post(a.op);
          }
        }
      }
      stream << ";";
    }
  }

  void SXFunction::codegen_declarations(CodeGenerator& g) const {

    // Make sure that there are no free variables
    if (!free_vars_.empty()) {
      casadi_error("Code generation of '" + name_ + "' is not possible since variables "
                   + str(free_vars_) + " are free.");
    }
  }

  void SXFunction::codegen_body(CodeGenerator& g) const {

    // Run the algorithm
    for (auto&& a : algorithm_) {
      if (a.op==OP_OUTPUT) {
        g << "if (res[" << a.i0 << "]!=0) "
          << g.res(a.i0) << "[" << a.i2 << "]=" << g.sx_work(a.i1);
      } else {

        // Where to store the result
        g << g.sx_work(a.i0) << "=";

        // What to store
        if (a.op==OP_CONST) {
          g << g.constant(a.d);
        } else if (a.op==OP_INPUT) {
          g << g.arg(a.i1) << "? " << g.arg(a.i1) << "[" << a.i2 << "] : 0";
        } else {
          casadi_int ndep = casadi_math<double>::ndeps(a.op);
          casadi_assert_dev(ndep>0);
          if (ndep==1) g << g.print_op(a.op, g.sx_work(a.i1));
          if (ndep==2) g << g.print_op(a.op, g.sx_work(a.i1), g.sx_work(a.i2));
        }
      }
      g  << ";\n";
    }
  }

  const Options SXFunction::options_
  = {{&FunctionInternal::options_},
     {{"default_in",
       {OT_DOUBLEVECTOR,
        "Default input values"}},
      {"just_in_time_sparsity",
       {OT_BOOL,
        "Propagate sparsity patterns using just-in-time "
        "compilation to a CPU or GPU using OpenCL"}},
      {"just_in_time_opencl",
       {OT_BOOL,
        "Just-in-time compilation for numeric evaluation using OpenCL (experimental)"}},
      {"live_variables",
       {OT_BOOL,
        "Reuse variables in the work vector"}},
      {"cse",
       {OT_BOOL,
        "Perform common subexpression elimination (complexity is N*log(N) in graph size)"}},
      {"allow_free",
       {OT_BOOL,
        "Allow construction with free variables (Default: false)"}},
      {"allow_duplicate_io_names",
       {OT_BOOL,
        "Allow construction with duplicate io names (Default: false)"}}
     }
  };

  Dict SXFunction::generate_options(const std::string& target) const {
    Dict opts = FunctionInternal::generate_options(target);
    //opts["default_in"] = default_in_;
    opts["live_variables"] = live_variables_;
    opts["just_in_time_sparsity"] = just_in_time_sparsity_;
    opts["just_in_time_opencl"] = just_in_time_opencl_;
    return opts;
  }

  void SXFunction::init(const Dict& opts) {
    // Call the init function of the base class
    XFunction<SXFunction, SX, SXNode>::init(opts);
    if (verbose_) casadi_message(name_ + "::init");

    // Default (temporary) options
    live_variables_ = true;

    bool cse_opt = false;
    bool allow_free = false;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="default_in") {
        default_in_ = op.second;
      } else if (op.first=="live_variables") {
        live_variables_ = op.second;
      } else if (op.first=="just_in_time_opencl") {
        just_in_time_opencl_ = op.second;
      } else if (op.first=="just_in_time_sparsity") {
        just_in_time_sparsity_ = op.second;
      } else if (op.first=="cse") {
        cse_opt = op.second;
      } else if (op.first=="allow_free") {
        allow_free = op.second;
      }
    }

    if (cse_opt) out_ = cse(out_);

    // Check/set default inputs
    if (default_in_.empty()) {
      default_in_.resize(n_in_, 0);
    } else {
      casadi_assert(default_in_.size()==n_in_,
                            "Option 'default_in' has incorrect length");
    }

    // Stack used to sort the computational graph
    std::stack<SXNode*> s;

    // All nodes
    std::vector<SXNode*> nodes;

    // Add the list of nodes
    casadi_int ind=0;
    for (auto it = out_.begin(); it != out_.end(); ++it, ++ind) {
      casadi_int nz=0;
      for (auto itc = (*it)->begin(); itc != (*it)->end(); ++itc, ++nz) {
        // Add outputs to the list
        s.push(itc->get());
        sort_depth_first(s, nodes);

        // A null pointer means an output instruction
        nodes.push_back(static_cast<SXNode*>(nullptr));
      }
    }

    casadi_assert(nodes.size() <= std::numeric_limits<int>::max(), "Integer overflow");
    // Set the temporary variables to be the corresponding place in the sorted graph
    for (casadi_int i=0; i<nodes.size(); ++i) {
      if (nodes[i]) {
        nodes[i]->temp = static_cast<int>(i);
      }
    }

    // Sort the nodes by type
    constants_.clear();
    operations_.clear();
    for (std::vector<SXNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
      SXNode* t = *it;
      if (t) {
        if (t->is_constant())
          constants_.push_back(SXElem::create(t));
        else if (!t->is_symbolic())
          operations_.push_back(SXElem::create(t));
      }
    }

    // Input instructions
    std::vector<std::pair<int, SXNode*> > symb_loc;

    // Current output and nonzero, start with the first one
    int curr_oind, curr_nz=0;
    casadi_assert(out_.size() <= std::numeric_limits<int>::max(), "Integer overflow");
    for (curr_oind=0; curr_oind<out_.size(); ++curr_oind) {
      if (out_[curr_oind].nnz()!=0) {
        break;
      }
    }

    // Count the number of times each node is used
    std::vector<casadi_int> refcount(nodes.size(), 0);

    // Get the sequence of instructions for the virtual machine
    algorithm_.resize(0);
    algorithm_.reserve(nodes.size());
    for (std::vector<SXNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
      // Current node
      SXNode* n = *it;

      // New element in the algorithm
      AlgEl ae;

      // Get operation
      ae.op = n==nullptr ? static_cast<int>(OP_OUTPUT) : static_cast<int>(n->op());

      // Get instruction
      switch (ae.op) {
      case OP_CONST: // constant
        ae.d = n->to_double();
        ae.i0 = n->temp;
        break;
      case OP_PARAMETER: // a parameter or input
        symb_loc.push_back(std::make_pair(algorithm_.size(), n));
        ae.i0 = n->temp;
        ae.d = 0; // value not used, but set here to avoid uninitialized data in serialization
        break;
      case OP_OUTPUT: // output instruction
        ae.i0 = curr_oind;
        ae.i1 = out_[curr_oind]->at(curr_nz)->temp;
        ae.i2 = curr_nz;

        // Go to the next nonzero
        casadi_assert(curr_nz < std::numeric_limits<int>::max(), "Integer overflow");
        curr_nz++;
        if (curr_nz>=out_[curr_oind].nnz()) {
          curr_nz=0;
          casadi_assert(curr_oind < std::numeric_limits<int>::max(), "Integer overflow");
          curr_oind++;
          for (; curr_oind<out_.size(); ++curr_oind) {
            if (out_[curr_oind].nnz()!=0) {
              break;
            }
          }
        }
        break;
      default:       // Unary or binary operation
        ae.i0 = n->temp;
        ae.i1 = n->dep(0).get()->temp;
        ae.i2 = n->dep(1).get()->temp;
      }

      // Number of dependencies
      casadi_int ndeps = casadi_math<double>::ndeps(ae.op);

      // Increase count of dependencies
      for (casadi_int c=0; c<ndeps; ++c) {
        refcount.at(c==0 ? ae.i1 : ae.i2)++;
      }
      // Add to algorithm
      algorithm_.push_back(ae);
    }

    // Place in the work vector for each of the nodes in the tree (overwrites the reference counter)
    std::vector<int> place(nodes.size());

    // Stack with unused elements in the work vector
    std::stack<int> unused;

    // Work vector size
    int worksize = 0;

    // Find a place in the work vector for the operation
    for (auto&& a : algorithm_) {

      // Number of dependencies
      casadi_int ndeps = casadi_math<double>::ndeps(a.op);

      // decrease reference count of children
      // reverse order so that the first argument will end up at the top of the stack
      for (casadi_int c=ndeps-1; c>=0; --c) {
        casadi_int ch_ind = c==0 ? a.i1 : a.i2;
        casadi_int remaining = --refcount.at(ch_ind);
        if (remaining==0) unused.push(place[ch_ind]);
      }

      // Find a place to store the variable
      if (a.op!=OP_OUTPUT) {
        if (live_variables_ && !unused.empty()) {
          // Try to reuse a variable from the stack if possible (last in, first out)
          a.i0 = place[a.i0] = unused.top();
          unused.pop();
        } else {
          // Allocate a new variable
          a.i0 = place[a.i0] = worksize++;
        }
      }

      // Save the location of the children
      for (casadi_int c=0; c<ndeps; ++c) {
        if (c==0) {
          a.i1 = place[a.i1];
        } else {
          a.i2 = place[a.i2];
        }
      }

      // If binary, make sure that the second argument is the same as the first one
      // (in order to treat all operations as binary) NOTE: ugly
      if (ndeps==1 && a.op!=OP_OUTPUT) {
        a.i2 = a.i1;
      }
    }

    worksize_ = worksize;

    if (verbose_) {
      if (live_variables_) {
        casadi_message("Using live variables: work array is " + str(worksize_)
         + " instead of " + str(nodes.size()));
      } else {
        casadi_message("Live variables disabled.");
      }
    }

    // Allocate work vectors (symbolic/numeric)
    alloc_w(worksize_);

    // Reset the temporary variables
    for (casadi_int i=0; i<nodes.size(); ++i) {
      if (nodes[i]) {
        nodes[i]->temp = 0;
      }
    }

    // Now mark each input's place in the algorithm
    for (auto it=symb_loc.begin(); it!=symb_loc.end(); ++it) {
      it->second->temp = it->first+1;
    }

    // Add input instructions
    casadi_assert(in_.size() <= std::numeric_limits<int>::max(), "Integer overflow");
    for (int ind=0; ind<in_.size(); ++ind) {
      int nz=0;
      for (auto itc = in_[ind]->begin(); itc != in_[ind]->end(); ++itc, ++nz) {
        int i = itc->get_temp()-1;
        if (i>=0) {
          // Mark as input
          algorithm_[i].op = OP_INPUT;

          // Location of the input
          algorithm_[i].i1 = ind;
          algorithm_[i].i2 = nz;

          // Mark input as read
          itc->set_temp(0);
        }
      }
    }

    // Locate free variables
    free_vars_.clear();
    for (std::vector<std::pair<int, SXNode*> >::const_iterator it=symb_loc.begin();
         it!=symb_loc.end(); ++it) {
      if (it->second->temp!=0) {
        // Save to list of free parameters
        free_vars_.push_back(SXElem::create(it->second));

        // Remove marker
        it->second->temp=0;
      }
    }

    if (!allow_free && has_free()) {
      casadi_error(name_ + "::init: Initialization failed since variables [" +
      join(get_free(), "") + "] are free. These symbols occur in the output expressions "
      "but you forgot to declare these as inputs. "
      "Set option 'allow_free' to allow free variables.");
    }

    // Initialize just-in-time compilation for numeric evaluation using OpenCL
    if (just_in_time_opencl_) {
      casadi_error("OpenCL is not supported in this version of CasADi");
    }

    // Initialize just-in-time compilation for sparsity propagation using OpenCL
    if (just_in_time_sparsity_) {
      casadi_error("OpenCL is not supported in this version of CasADi");
    }

    // Print
    if (verbose_) casadi_message(str(algorithm_.size()) + " elementary operations");
  }

  SX SXFunction::instructions_sx() const {
    std::vector<SXElem> ret(algorithm_.size(), casadi_limits<SXElem>::nan);

    std::vector<SXElem>::iterator it=ret.begin();

    // Iterator to the binary operations
    std::vector<SXElem>::const_iterator b_it = operations_.begin();

    // Iterator to stack of constants
    std::vector<SXElem>::const_iterator c_it = constants_.begin();

    // Iterator to free variables
    std::vector<SXElem>::const_iterator p_it = free_vars_.begin();

    // Evaluate algorithm
    if (verbose_) casadi_message("Evaluating algorithm forward");
    for (auto&& a : algorithm_) {
      switch (a.op) {
      case OP_INPUT:
      case OP_OUTPUT:
        it++;
        break;
      case OP_CONST:
        *it++ = *c_it++;
        break;
      case OP_PARAMETER:
        *it++ = *p_it++;
        break;
      default:
        *it++ = *b_it++;
      }
    }
    casadi_assert(it==ret.end(), "Dimension mismatch");
    return ret;
  }

  int SXFunction::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w, void* mem) const {
    if (verbose_) casadi_message(name_ + "::eval_sx");

    // Iterator to the binary operations
    std::vector<SXElem>::const_iterator b_it=operations_.begin();

    // Iterator to stack of constants
    std::vector<SXElem>::const_iterator c_it = constants_.begin();

    // Iterator to free variables
    std::vector<SXElem>::const_iterator p_it = free_vars_.begin();

    // Evaluate algorithm
    if (verbose_) casadi_message("Evaluating algorithm forward");
    for (auto&& a : algorithm_) {
      switch (a.op) {
      case OP_INPUT:
        w[a.i0] = arg[a.i1]==nullptr ? 0 : arg[a.i1][a.i2];
        break;
      case OP_OUTPUT:
        if (res[a.i0]!=nullptr) res[a.i0][a.i2] = w[a.i1];
        break;
      case OP_CONST:
        w[a.i0] = *c_it++;
        break;
      case OP_PARAMETER:
        w[a.i0] = *p_it++; break;
      default:
        {
          // Evaluate the function to a temporary value
          // (as it might overwrite the children in the work vector)
          SXElem f;
          switch (a.op) {
            CASADI_MATH_FUN_BUILTIN(w[a.i1], w[a.i2], f)
          }

          // If this new expression is identical to the expression used
          // to define the algorithm, then reuse
          const casadi_int depth = 2; // NOTE: a higher depth could possibly give more savings
          f.assignIfDuplicate(*b_it++, depth);

          // Finally save the function value
          w[a.i0] = f;
        }
      }
    }
    return 0;
  }

  void SXFunction::ad_forward(const std::vector<std::vector<SX> >& fseed,
                                std::vector<std::vector<SX> >& fsens) const {
    if (verbose_) casadi_message(name_ + "::ad_forward");

    // Number of forward seeds
    casadi_int nfwd = fseed.size();
    fsens.resize(nfwd);

    // Quick return if possible
    if (nfwd==0) return;

    // Check if seeds need to have dimensions corrected
    casadi_int npar = 1;
    for (auto&& r : fseed) {
      if (!matching_arg(r, npar)) {
        casadi_assert_dev(npar==1);
        return ad_forward(replace_fseed(fseed, npar), fsens);
      }
    }

    // Make sure seeds have matching sparsity patterns
    for (auto it=fseed.begin(); it!=fseed.end(); ++it) {
      casadi_assert_dev(it->size()==n_in_);
      for (casadi_int i=0; i<n_in_; ++i) {
        if (it->at(i).sparsity()!=sparsity_in_[i]) {
          // Correct sparsity
          std::vector<std::vector<SX> > fseed2(fseed);
          for (auto&& r : fseed2) {
            for (casadi_int i=0; i<n_in_; ++i) r[i] = project(r[i], sparsity_in_[i]);
          }
          return ad_forward(fseed2, fsens);
        }
      }
    }

    // Allocate results
    for (casadi_int d=0; d<nfwd; ++d) {
      fsens[d].resize(n_out_);
      for (casadi_int i=0; i<fsens[d].size(); ++i)
        if (fsens[d][i].sparsity()!=sparsity_out_[i])
          fsens[d][i] = SX::zeros(sparsity_out_[i]);
    }

    // Iterator to the binary operations
    std::vector<SXElem>::const_iterator b_it=operations_.begin();

    // Tape
    std::vector<TapeEl<SXElem> > s_pdwork(operations_.size());
    std::vector<TapeEl<SXElem> >::iterator it1 = s_pdwork.begin();

    // Evaluate algorithm
    if (verbose_) casadi_message("Evaluating algorithm forward");
    for (auto&& e : algorithm_) {
      switch (e.op) {
      case OP_INPUT:
      case OP_OUTPUT:
      case OP_CONST:
      case OP_PARAMETER:
        break;
      default:
        {
          const SXElem& f=*b_it++;
          switch (e.op) {
            CASADI_MATH_DER_BUILTIN(f->dep(0), f->dep(1), f, it1++->d)
          }
        }
      }
    }

    // Work vector
    std::vector<SXElem> w(worksize_);

    // Calculate forward sensitivities
    if (verbose_) casadi_message("Calculating forward derivatives");
    for (casadi_int dir=0; dir<nfwd; ++dir) {
      std::vector<TapeEl<SXElem> >::const_iterator it2 = s_pdwork.begin();
      for (auto&& a : algorithm_) {
        switch (a.op) {
        case OP_INPUT:
          w[a.i0] = fseed[dir][a.i1].nonzeros()[a.i2]; break;
        case OP_OUTPUT:
          fsens[dir][a.i0].nonzeros()[a.i2] = w[a.i1]; break;
        case OP_CONST:
        case OP_PARAMETER:
          w[a.i0] = 0;
          break;
        case OP_IF_ELSE_ZERO:
          w[a.i0] = if_else_zero(it2++->d[1], w[a.i2]);
          break;
        CASADI_MATH_BINARY_BUILTIN // Binary operation
          w[a.i0] = it2->d[0] * w[a.i1] + it2->d[1] * w[a.i2];
          it2++;
          break;
        default: // Unary operation
          w[a.i0] = it2->d[0] * w[a.i1];
          it2++;
        }
      }
    }
  }

  void SXFunction::ad_reverse(const std::vector<std::vector<SX> >& aseed,
                                std::vector<std::vector<SX> >& asens) const {
    if (verbose_) casadi_message(name_ + "::ad_reverse");

    // number of adjoint seeds
    casadi_int nadj = aseed.size();
    asens.resize(nadj);

    // Quick return if possible
    if (nadj==0) return;

    // Check if seeds need to have dimensions corrected
    casadi_int npar = 1;
    for (auto&& r : aseed) {
      if (!matching_res(r, npar)) {
        casadi_assert_dev(npar==1);
        return ad_reverse(replace_aseed(aseed, npar), asens);
      }
    }

    // Make sure matching sparsity of fseed
    bool matching_sparsity = true;
    for (casadi_int d=0; d<nadj; ++d) {
      casadi_assert_dev(aseed[d].size()==n_out_);
      for (casadi_int i=0; matching_sparsity && i<n_out_; ++i)
        matching_sparsity = aseed[d][i].sparsity()==sparsity_out_[i];
    }

    // Correct sparsity if needed
    if (!matching_sparsity) {
      std::vector<std::vector<SX> > aseed2(aseed);
      for (casadi_int d=0; d<nadj; ++d)
        for (casadi_int i=0; i<n_out_; ++i)
          if (aseed2[d][i].sparsity()!=sparsity_out_[i])
            aseed2[d][i] = project(aseed2[d][i], sparsity_out_[i]);
      return ad_reverse(aseed2, asens);
    }

    // Allocate results if needed
    for (casadi_int d=0; d<nadj; ++d) {
      asens[d].resize(n_in_);
      for (casadi_int i=0; i<asens[d].size(); ++i) {
        if (asens[d][i].sparsity()!=sparsity_in_[i]) {
          asens[d][i] = SX::zeros(sparsity_in_[i]);
        } else {
          std::fill(asens[d][i]->begin(), asens[d][i]->end(), 0);
        }
      }
    }

    // Iterator to the binary operations
    std::vector<SXElem>::const_iterator b_it=operations_.begin();

    // Tape
    std::vector<TapeEl<SXElem> > s_pdwork(operations_.size());
    std::vector<TapeEl<SXElem> >::iterator it1 = s_pdwork.begin();

    // Evaluate algorithm
    if (verbose_) casadi_message("Evaluating algorithm forward");
    for (auto&& a : algorithm_) {
      switch (a.op) {
      case OP_INPUT:
      case OP_OUTPUT:
      case OP_CONST:
      case OP_PARAMETER:
        break;
      default:
        {
          const SXElem& f=*b_it++;
          switch (a.op) {
            CASADI_MATH_DER_BUILTIN(f->dep(0), f->dep(1), f, it1++->d)
          }
        }
      }
    }

    // Calculate adjoint sensitivities
    if (verbose_) casadi_message("Calculating adjoint derivatives");

    // Work vector
    std::vector<SXElem> w(worksize_, 0);

    for (casadi_int dir=0; dir<nadj; ++dir) {
      auto it2 = s_pdwork.rbegin();
      for (auto it = algorithm_.rbegin(); it!=algorithm_.rend(); ++it) {
        SXElem seed;
        switch (it->op) {
        case OP_INPUT:
          asens[dir][it->i1].nonzeros()[it->i2] = w[it->i0];
          w[it->i0] = 0;
          break;
        case OP_OUTPUT:
          w[it->i1] += aseed[dir][it->i0].nonzeros()[it->i2];
          break;
        case OP_CONST:
        case OP_PARAMETER:
          w[it->i0] = 0;
          break;
        case OP_IF_ELSE_ZERO:
          seed = w[it->i0];
          w[it->i0] = 0;
          w[it->i2] += if_else_zero(it2++->d[1], seed);
          break;
        CASADI_MATH_BINARY_BUILTIN // Binary operation
          seed = w[it->i0];
          w[it->i0] = 0;
          w[it->i1] += it2->d[0] * seed;
          w[it->i2] += it2++->d[1] * seed;
          break;
        default: // Unary operation
          seed = w[it->i0];
          w[it->i0] = 0;
          w[it->i1] += it2++->d[0] * seed;
        }
      }
    }
  }

  int SXFunction::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    // Fall back when forward mode not allowed
    if (sp_weight()==1 || sp_weight()==-1)
      return FunctionInternal::sp_forward(arg, res, iw, w, mem);
    // Propagate sparsity forward
    for (auto&& e : algorithm_) {
      switch (e.op) {
      case OP_CONST:
      case OP_PARAMETER:
        w[e.i0] = 0; break;
      case OP_INPUT:
        w[e.i0] = arg[e.i1]==nullptr ? 0 : arg[e.i1][e.i2];
        break;
      case OP_OUTPUT:
        if (res[e.i0]!=nullptr) res[e.i0][e.i2] = w[e.i1];
        break;
      default: // Unary or binary operation
        w[e.i0] = w[e.i1] | w[e.i2]; break;
      }
    }
    return 0;
  }

  int SXFunction::sp_reverse(bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    // Fall back when reverse mode not allowed
    if (sp_weight()==0 || sp_weight()==-1)
      return FunctionInternal::sp_reverse(arg, res, iw, w, mem);
    std::fill_n(w, sz_w(), 0);

    // Propagate sparsity backward
    for (auto it=algorithm_.rbegin(); it!=algorithm_.rend(); ++it) {
      // Temp seed
      bvec_t seed;

      // Propagate seeds
      switch (it->op) {
      case OP_CONST:
      case OP_PARAMETER:
        w[it->i0] = 0;
        break;
      case OP_INPUT:
        if (arg[it->i1]!=nullptr) arg[it->i1][it->i2] |= w[it->i0];
        w[it->i0] = 0;
        break;
      case OP_OUTPUT:
        if (res[it->i0]!=nullptr) {
          w[it->i1] |= res[it->i0][it->i2];
          res[it->i0][it->i2] = 0;
        }
        break;
      default: // Unary or binary operation
        seed = w[it->i0];
        w[it->i0] = 0;
        w[it->i1] |= seed;
        w[it->i2] |= seed;
      }
    }
    return 0;
  }

  const SX SXFunction::sx_in(casadi_int ind) const {
    return in_.at(ind);
  }

  const std::vector<SX> SXFunction::sx_in() const {
    return in_;
  }

  bool SXFunction::is_a(const std::string& type, bool recursive) const {
    return type=="SXFunction" || (recursive && XFunction<SXFunction,
                                  SX, SXNode>::is_a(type, recursive));
  }

  void SXFunction::export_code_body(const std::string& lang,
      std::ostream &ss, const Dict& options) const {

    // Default values for options
    casadi_int indent_level = 0;

    // Read options
    for (auto&& op : options) {
      if (op.first=="indent_level") {
        indent_level = op.second;
      } else {
        casadi_error("Unknown option '" + op.first + "'.");
      }
    }

    // Construct indent string
    std::string indent;
    for (casadi_int i=0;i<indent_level;++i) {
      indent += "  ";
    }

    // Non-cell aliases for inputs
    for (casadi_int i=0;i<n_in_;++i) {
      ss << indent << "argin_" << i <<  " = nonzeros_gen(varargin{" << i+1 << "});" << std::endl;
    }

    Function f = shared_from_this<Function>();

    for (casadi_int k=0;k<f.n_instructions();++k) {
      // Get operation
      casadi_int op = static_cast<casadi_int>(f.instruction_id(k));
      // Get input positions into workvector
      std::vector<casadi_int> o = f.instruction_output(k);
      // Get output positions into workvector
      std::vector<casadi_int> i = f.instruction_input(k);
      switch (op) {
        case OP_INPUT:
          {
            ss << indent << "w" << o[0] << " = " << "argin_" << i[0] << "(" << i[1]+1 << ");";
            ss << std::endl;
          }
          break;
        case OP_OUTPUT:
          {
            ss << indent << "argout_" << o[0] << "{" << o[1]+1 << "} = w" << i[0] << ";";
            ss << std::endl;
          }
          break;
        case OP_CONST:
          {
            std::ios_base::fmtflags fmtfl = ss.flags();
            ss << indent << "w" << o[0] << " = ";
            ss << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
            ss << f.instruction_constant(k) << ";" << std::endl;
            ss.flags(fmtfl);
          }
          break;
        case OP_SQ:
          {
            ss << indent << "w" << o[0] << " = " << "w" << i[0] << "^2;" << std::endl;
          }
          break;
        case OP_FABS:
          {
            ss << indent << "w" << o[0] << " = abs(" << "w" << i[0] << ");" << std::endl;
          }
          break;
        case OP_POW:
        case OP_CONSTPOW:
          ss << indent << "w" << o[0] << " = " << "w" << i[0] << ".^w" << i[1] << ";" << std::endl;
          break;
        case OP_NOT:
          ss << indent << "w" << o[0] << " = ~" << "w" << i[0] << ";" << std::endl;
          break;
        case OP_OR:
          ss << indent << "w" << o[0] << " = w" << i[0] << " | w" << i[1] << ";" << std::endl;
          break;
        case OP_AND:
          ss << indent << "w" << o[0] << " = w" << i[0] << " & w" << i[1] << ";" << std::endl;
          break;
        case OP_NE:
          ss << indent << "w" << o[0] << " = w" << i[0] << " ~= w" << i[1] << ";" << std::endl;
          break;
        case OP_IF_ELSE_ZERO:
          ss << indent << "w" << o[0] << " = ";
          ss << "if_else_zero_gen(w" << i[0] << ", w" << i[1] << ");" << std::endl;
          break;
        default:
          if (casadi::casadi_math<double>::ndeps(op)==2) {
            ss << indent << "w" << o[0] << " = " << casadi::casadi_math<double>::print(op,
              "w"+std::to_string(i[0]), "w"+std::to_string(i[1])) << ";" << std::endl;
          } else {
            ss << indent << "w" << o[0] << " = " << casadi::casadi_math<double>::print(op,
              "w"+std::to_string(i[0])) << ";" << std::endl;
          }
      }
    }

  }

  SXFunction::SXFunction(DeserializingStream& s) :
    XFunction<SXFunction, SX, SXNode>(s) {
    s.version("SXFunction", 1);
    size_t n_instructions;
    s.unpack("SXFunction::n_instr", n_instructions);

    s.unpack("SXFunction::worksize", worksize_);
    s.unpack("SXFunction::free_vars", free_vars_);
    s.unpack("SXFunction::operations", operations_);
    s.unpack("SXFunction::constants", constants_);
    s.unpack("SXFunction::default_in", default_in_);

    algorithm_.resize(n_instructions);
    for (casadi_int k=0;k<n_instructions;++k) {
      AlgEl& e = algorithm_[k];
      s.unpack("SXFunction::ScalarAtomic::op", e.op);
      s.unpack("SXFunction::ScalarAtomic::i0", e.i0);
      s.unpack("SXFunction::ScalarAtomic::i1", e.i1);
      s.unpack("SXFunction::ScalarAtomic::i2", e.i2);
    }

    // Default (persistent) options
    just_in_time_opencl_ = false;
    just_in_time_sparsity_ = false;

    s.unpack("SXFunction::live_variables", live_variables_);

    XFunction<SXFunction, SX, SXNode>::delayed_deserialize_members(s);
  }

  void SXFunction::serialize_body(SerializingStream &s) const {
    XFunction<SXFunction, SX, SXNode>::serialize_body(s);
    s.version("SXFunction", 1);
    s.pack("SXFunction::n_instr", algorithm_.size());

    s.pack("SXFunction::worksize", worksize_);
    s.pack("SXFunction::free_vars", free_vars_);
    s.pack("SXFunction::operations", operations_);
    s.pack("SXFunction::constants", constants_);
    s.pack("SXFunction::default_in", default_in_);

    // Loop over algorithm
    for (const auto& e : algorithm_) {
      s.pack("SXFunction::ScalarAtomic::op", e.op);
      s.pack("SXFunction::ScalarAtomic::i0", e.i0);
      s.pack("SXFunction::ScalarAtomic::i1", e.i1);
      s.pack("SXFunction::ScalarAtomic::i2", e.i2);
    }

    s.pack("SXFunction::live_variables", live_variables_);

    XFunction<SXFunction, SX, SXNode>::delayed_serialize_members(s);
  }

  ProtoFunction* SXFunction::deserialize(DeserializingStream& s) {
    return new SXFunction(s);
  }

} // namespace casadi
