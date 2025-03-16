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
#include <bitset>
#include "sx_node.hpp"
#include "output_sx.hpp"
#include "call_sx.hpp"
#include "casadi_common.hpp"
#include "sparsity_internal.hpp"
#include "casadi_interrupt.hpp"
#include "serializing_stream.hpp"
#include "global_options.hpp"

namespace casadi {

  SXFunction::ExtendedAlgEl::ExtendedAlgEl(const Function& fun) : f(fun) {
    n_dep = f.nnz_in(); n_res = f.nnz_out();
    dep.resize(n_dep); res.resize(n_res, -1);
    f_n_in = f.n_in(); f_n_out = f.n_out();
    f_nnz_in.resize(f_n_in); f_nnz_out.resize(f_n_out);
    for (casadi_int i=0;i<f_n_in;++i) f_nnz_in[i] = f.nnz_in(i);
    for (casadi_int i=0;i<f_n_out;++i) f_nnz_out[i] = f.nnz_out(i);
    copy_elision_arg.resize(f_n_in, -1);
    copy_elision_offset.resize(f_n_in, -1);
  }

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
      case OP_CALL:
        call_fwd(e, arg, res, iw, w);
      break;
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
      } else if (a.op==OP_CALL) {
        const ExtendedAlgEl& m = call_.el.at(a.i1);
        stream << "[";
        casadi_int k = 0;
        for (casadi_int i=0; i<m.f.n_out(); ++i) {
          if (m.f.nnz_out(i)>1) stream << "[";
          for (casadi_int j=0; j<m.f.nnz_out(i); ++j) {
            int el = m.res[k++];
            if (el>=0) {
              stream << "@" << el;
            } else {
              stream << "NULL";
            }
            if (j<m.f.nnz_out(i)-1) stream << ",";
          }
          if (m.f.nnz_out(i)>1) stream << "]";
          if (i<m.f.n_out()-1) stream << ",";
        }
        stream << "] = ";
        stream << m.f.name() << "(";
        k = 0;
        for (casadi_int i=0; i<m.f.n_in(); ++i) {
          if (m.f.nnz_in(i)==0) stream << "0x0";
          if (m.f.nnz_in(i)>1) stream << "[";
          for (casadi_int j=0; j<m.f.nnz_in(i); ++j) {
            stream << "@" << m.dep[k++];
            if (j<m.f.nnz_in(i)-1) stream << ",";
          }
          if (m.f.nnz_in(i)>1) stream << "]";
          if (i<m.f.n_in()-1) stream << ",";
        }
        stream << ")";
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

  size_t SXFunction::codegen_sz_w(const CodeGenerator& g) const {
    if (!g.avoid_stack()) return call_.sz_w+call_.sz_w_arg+call_.sz_w_res;
    return sz_w();
  }

  void SXFunction::codegen_declarations(CodeGenerator& g) const {

    // Make sure that there are no free variables
    if (!free_vars_.empty()) {
      casadi_error("Code generation of '" + name_ + "' is not possible since variables "
                   + str(free_vars_) + " are free.");
    }

    // Generate code for the call nodes
    for (auto&& m : call_.el) {
      g.add_dependency(m.f);
    }
  }

  void SXFunction::codegen_body(CodeGenerator& g) const {
    g.reserve_work(worksize_);

    casadi_int k=0;
    // Run the algorithm
    for (auto&& a : algorithm_) {
      if (a.op==OP_OUTPUT) {
        g << "if (res[" << a.i0 << "]!=0) "
          << g.res(a.i0) << "[" << a.i2 << "]=" << g.sx_work(a.i1) << ";\n";
      } else if (a.op==OP_CALL) {
        const ExtendedAlgEl& m = call_.el[a.i1];

        casadi_int worksize = g.avoid_stack() ? worksize_ : 0;

        // Collect input arguments
        casadi_int offset = worksize;
        for (casadi_int i=0; i<m.f_n_in; ++i) {
          if (m.copy_elision_arg[i]>=0) {
            g << "arg[" << n_in_+i << "] = "
              << "arg[" + str(m.copy_elision_arg[i]) << "]? "
              << "arg[" + str(m.copy_elision_arg[i]) << "] + "
              << str(m.copy_elision_offset[i]) << " : 0;\n";
          } else {
            if (m.f_nnz_in[i]==0) {
              g << "arg[" << n_in_+i << "]=" << 0 << ";\n";
            } else {
              g << "arg[" << n_in_+i << "]=" << "w+" + str(offset) << ";\n";
            }
          }
          offset += m.f_nnz_in[i];
        }


        casadi_int out_offset = offset;

        // Collect output arguments
        for (casadi_int i=0; i<m.f_n_out; ++i) {
          g << "res[" << n_out_+i << "]=" << "w+" + str(offset) << ";\n";
          offset += m.f_nnz_out[i];
        }
        casadi_int k=0;
        for (casadi_int i=0; i<m.f_n_in; ++i) {
          if (m.copy_elision_arg[i]==-1) {
            for (casadi_int j=0; j<m.f_nnz_in[i]; ++j) {
              g << "w["+str(k+worksize) + "] = " << g.sx_work(m.dep[k]) << ";\n";
              k++;
            }
          } else {
            k+=m.f_nnz_in[i];
          }
        }
        std::string flag =
          g(m.f, "arg+"+str(n_in_), "res+"+str(n_out_), "iw", "w+" + str(offset));
        // Call function
        g << "if (" << flag << ") return 1;\n";
        for (casadi_int i=0;i<m.n_res;++i) {
          if (m.res[i]>=0) {
            g << g.sx_work(m.res[i]) << " = ";
            g << "w[" + str(i+out_offset) + "];\n";
          }
        }
      } else if (a.op==OP_INPUT) {
          if (!copy_elision_[k]) {
            g << g.sx_work(a.i0) << "="
              << g.arg(a.i1) << "? " << g.arg(a.i1) << "[" << a.i2 << "] : 0;\n";
          }
      } else {

        // Where to store the result
        g << g.sx_work(a.i0) << "=";

        // What to store
        if (a.op==OP_CONST) {
          g << g.constant(a.d);
        } else {
          casadi_int ndep = casadi_math<double>::ndeps(a.op);
          casadi_assert_dev(ndep>0);
          if (ndep==1) g << g.print_op(a.op, g.sx_work(a.i1));
          if (ndep==2) g << g.print_op(a.op, g.sx_work(a.i1), g.sx_work(a.i2));
        }
        g  << ";\n";
      }
      k++;
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

    // Perform common subexpression elimination
    // This must be done before the lock, to avoid deadlocks
    if (cse_opt) out_ = cse(out_);

    // Check/set default inputs
    if (default_in_.empty()) {
      default_in_.resize(n_in_, 0);
    } else {
      casadi_assert(default_in_.size()==n_in_,
                            "Option 'default_in' has incorrect length");
    }

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    std::lock_guard<std::mutex> lock(SX::get_mutex_temp());
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

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
        else if (!t->is_symbolic() && t->op()>=0)
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

    // Mapping of node index (cfr. temp) to algorithm index
    std::vector<int> alg_index;
    alg_index.reserve(nodes.size());

    for (std::vector<SXNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
      // Current node
      SXNode* n = *it;

      // New element in the algorithm
      AlgEl ae;

      // Get operation
      ae.op = n==nullptr ? static_cast<int>(OP_OUTPUT) : static_cast<int>(n->op());

      // Default dependencies
      int* dep = &ae.i1;
      casadi_int ndeps = ae.op == -1 ? 1 : casadi_math<double>::ndeps(ae.op);

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
      case OP_CALL: // Call node
        {
          ae.i0 = n->temp;

          // Index into ExtentedAlgEl collection
          ae.i1 = call_.el.size();

          // Create ExtentedAlgEl instance
          // This allocates space for dep and res
          const Function& f = static_cast<const CallSX*>(n)->f_;
          call_.el.emplace_back(f);

          // Make sure we have enough space to evaluate the Function call,
          // noting that we wil only ever evaluate one call at a time.
          call_.sz_arg = std::max(call_.sz_arg, f.sz_arg());
          call_.sz_res = std::max(call_.sz_res, f.sz_res());
          call_.sz_iw  = std::max(call_.sz_iw, f.sz_iw());
          call_.sz_w   = std::max(call_.sz_w, f.sz_w());
          call_.sz_w_arg   = std::max(call_.sz_w_arg, static_cast<size_t>(f.nnz_in()));
          call_.sz_w_res   = std::max(call_.sz_w_res,  static_cast<size_t>(f.nnz_out()));

          // Set the dependency pointer to the (uninitialised) slots of the ExtendedAlgEl
          ExtendedAlgEl& m = call_.el.at(ae.i1);
          dep = get_ptr(m.dep);
          ndeps = m.n_dep;

          // Populate the dependency slots with node ids.
          for (casadi_int i=0; i<ndeps; ++i) {
            dep[i] = n->dep(i).get()->temp;
          }
        }
        break;
      case -1: // Output extraction node
        {
          dep = &algorithm_.at(alg_index.at(n->dep(0).get()->temp)).i1;
          int oind = static_cast<OutputSX*>(n)->oind_;
          casadi_assert(call_.el.at(dep[0]).res.at(oind)==-1, "Duplicate");
          call_.el.at(dep[0]).res.at(oind) = n->temp;
        }
        break;
      default:       // Unary or binary operation
        ae.i0 = n->temp;
        ae.i1 = n->dep(0).get()->temp;
        ae.i2 = n->dep(1).get()->temp;
      }

      // Increase count of dependencies
      for (casadi_int c=0; c<ndeps; ++c) {
        refcount.at(dep[c])++;
      }

      // Amend node index to algorithm index mapping
      alg_index.push_back(algorithm_.size());

      // Add to algorithm
      if (ae.op>=0) algorithm_.push_back(ae);

    }

    // Place in the work vector for each of the nodes in the tree (overwrites the reference counter)
    std::vector<int> place(nodes.size());

    // Stack with unused elements in the work vector
    std::stack<int> unused;

    // Work vector size
    int worksize = 0;

    // Find a place in the work vector for the operation
    for (auto&& a : algorithm_) {

      // Default dependencies
      int* dep = &a.i1;
      casadi_int ndeps = casadi_math<double>::ndeps(a.op);

      // Default outputs
      int* res = &a.i0;
      casadi_int nres = 1;

      // Call node overrides these defaults
      if (a.op==OP_CALL) {
        ExtendedAlgEl& e = call_.el.at(a.i1);
        ndeps = e.n_dep;
        dep = get_ptr(e.dep);
        nres = e.n_res;
        res = get_ptr(e.res);
      }

      // decrease reference count of children
      // reverse order so that the first argument will end up at the top of the stack
      for (casadi_int c=ndeps-1; c>=0; --c) {
        casadi_int ch_ind = dep[c];
        casadi_int remaining = --refcount.at(ch_ind);
        if (remaining==0) unused.push(place[ch_ind]);
      }

      // Find a place to store the variable
      if (a.op!=OP_OUTPUT) {
        for (casadi_int c=0; c<nres; ++c) {
          if (res[c]<0) continue;
          if (live_variables_ && !unused.empty()) {
            // Try to reuse a variable from the stack if possible (last in, first out)
            res[c] = place[res[c]] = unused.top();
            unused.pop();
          } else {
            // Allocate a new variable
            res[c] = place[res[c]] = worksize++;
          }
        }
      }

      // Save the location of the children
      for (casadi_int c=0; c<ndeps; ++c) {
        dep[c] = place[dep[c]];
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

    alloc_arg(call_.sz_arg, true);
    alloc_res(call_.sz_res, true);
    alloc_iw(call_.sz_iw, true);
    alloc_w(call_.sz_w+call_.sz_w_arg+call_.sz_w_res, true);

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
      join(get_free(), ", ") + "] are free. These symbols occur in the output expressions "
      "but you forgot to declare these as inputs. "
      "Set option 'allow_free' to allow free variables.");
    }

    init_copy_elision();

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

  void SXFunction::init_copy_elision() {
    if (GlobalOptions::copy_elision_min_size==-1) {
      copy_elision_.resize(algorithm_.size(), false);
      return;
    }
    // Perform copy elision (codegen-only)
    // Remove nodes that only serve to compose CALL inputs

    // For work vector elements, store the arg source (-1 for no trivial source)
    std::vector<int> arg_i(worksize_, -1);
    std::vector<int> nz_i(worksize_, -1);

    // Which algel corresponds to this source?
    std::vector<casadi_int> alg_i(worksize_, -1);

    // Is this algel to be elided?
    copy_elision_.resize(algorithm_.size(), false);

    casadi_int k=0;
    for (auto&& e : algorithm_) {
      switch (e.op) {
      case OP_INPUT:
        // Make source association
        arg_i[e.i0] = e.i1;
        nz_i[e.i0] = e.i2;
        alg_i[e.i0] = k;
        copy_elision_[k] = true;
        break;
      case OP_OUTPUT:
        if (arg_i[e.i1]>=0) {
          copy_elision_[alg_i[e.i1]] = false;
        }
        break;
      case OP_CALL:
        {
          auto& m = call_.el[e.i1];

          // Inspect input arguments
          casadi_int offset_input = 0;
          for (casadi_int i=0; i<m.f_n_in; ++i) {
            // Pattern match results
            casadi_int arg = -1;
            casadi_int offset = -1;
            for (casadi_int j=0; j<m.f_nnz_in[i]; ++j) {
              casadi_int k = offset_input+j;
              if (j==0) {
                arg = arg_i[m.dep[k]];
                offset = nz_i[m.dep[k]];
              }
              if (arg_i[m.dep[k]]==-1) {
                arg = -1;
                // Pattern match failed
                break;
              }
              if (nz_i[m.dep[k]]!=offset+j) {
                arg = -1;
                // Pattern match failed
                break;
              }
            }

            // If we cannot perform elision
            if (arg==-1) {
              // We need copies for all nonzeros of input i
              for (casadi_int j=0; j<m.f_nnz_in[i]; ++j) {
                casadi_int k = offset_input+j;
                if (arg_i[m.dep[k]]>=0) {
                  copy_elision_[alg_i[m.dep[k]]] = false;
                }
              }
            }
            // Store pattern match results
            m.copy_elision_arg[i] = arg;
            m.copy_elision_offset[i] = offset;

            offset += m.f_nnz_in[i];
            offset_input += m.f_nnz_in[i];
          }

          // Remove source association of all outputs
          for (casadi_int i=0; i<m.n_res; ++i) {
            if (m.res[i]>=0) {
              arg_i[m.res[i]] = -1;
            }
          }
        }
        break;
      case OP_CONST:
      case OP_PARAMETER:
        // Remove source association
        arg_i[e.i0] = -1;
        break;
      default:
        if (arg_i[e.i1]>=0) {
          copy_elision_[alg_i[e.i1]] = false;
        }
        if (!casadi_math<double>::is_unary(e.op)) {
          if (arg_i[e.i2]>=0) {
            copy_elision_[alg_i[e.i2]] = false;
          }
        }
        // Remove source association
        arg_i[e.i0] = -1;
      }
      k++;
    }
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
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w, void* mem,
    bool always_inline, bool never_inline) const {

    always_inline = always_inline || always_inline_;
    never_inline = never_inline || never_inline_;

    // non-inlining call is implemented in the base-class
    if (!should_inline(true, always_inline, never_inline)) {
      return FunctionInternal::eval_sx(arg, res, iw, w, mem, false, true);
    }

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
      case OP_CALL:
        {
          const ExtendedAlgEl& m = call_.el.at(a.i1);
          const SXElem& orig = *b_it++;
          std::vector<SXElem> deps(m.n_dep);
          bool identical = true;

          std::vector<SXElem> ret;
          for (casadi_int i=0;i<m.n_dep;++i) {
            identical &= SXElem::is_equal(w[m.dep.at(i)], orig->dep(i), 2);
          }
          if (identical) {
            ret = OutputSX::split(orig, m.n_res);
          } else {
            for (casadi_int i=0;i<m.n_dep;++i) deps[i] = w[m.dep[i]];
            ret = SXElem::call(m.f, deps);
          }
          for (casadi_int i=0;i<m.n_res;++i) {
            if (m.res[i]>=0) w[m.res[i]] = ret[i];
          }
        }
        break;
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


  template<typename T>
  inline void my_propagate_interval(unsigned char op, const T& L1, const T& R1,
        const T& L2, const T& R2, T& L, T& R) {
    switch (op) {
      case OP_NEG:
        L = -R1;
        R = -L1;
        break;
      case OP_TWICE:
        L = 2*L1;
        R = 2*R1;
        break;
      case OP_FMIN:
        L = fmin(L1, L2);
        R = fmin(R1, R2);
        break;
      case OP_FMAX:
        L = fmax(L1, L2);
        R = fmax(R1, R2);
        break;
      case OP_ADD:
        L = L1 + L2;
        R = R1 + R2;
        break;
      case OP_SUB:
        L = L1 - R2;
        R = R1 - L2;
        break;
      case OP_SQRT:
        L = sqrt(fmax(L1, 0));
        R = sqrt(fmax(R1, 0));
        break;
      case OP_LE:
        L = L2>=R1;
        R = R2>=L1;
        break;
      case OP_LT:
        L = L2>R1;
        R = R2>L1;
        break;
      case OP_EQ:
        L = logic_and(logic_and(L1==R1, L2==R2), L1==L2);
        R = logic_and(R1>=L2, L1<=R2);
        break;
      case OP_NE:
        {
          T singular = logic_and(logic_and(L1==R1, L2==R2), L1==L2);
          L = logic_not(logic_and(R1>=L2, L1<=R2));
          R = logic_not(singular);
          
        }
        break;
      case OP_NOT:
        {
          T contains_zero = logic_and(L1<=0, R1>=0);
          L = logic_and(L1==R1, L1==0);
          R = contains_zero;
        }
        break;
      case OP_AND:
        {
          T contains_zero1 = logic_and(L1<=0, R1>=0);
          T contains_zero2 = logic_and(L2<=0, R2>=0);
          L = logic_not(logic_or(contains_zero1, contains_zero2));
          R = logic_and(logic_or(L1!=R1, L1!=0), logic_or(L2!=R2, L2!=0));
        }
        break;
      case OP_OR:
        {
          T is_zero1 = logic_and(L1==R1, L1==0);
          T is_zero2 = logic_and(L2==R2, L2==0);
          T contains_zero1 = logic_and(L1<=0, R1>=0);
          T contains_zero2 = logic_and(L2<=0, R2>=0);
          L = logic_not(logic_and(contains_zero1, contains_zero2));
          R = logic_not(logic_and(is_zero1, is_zero2));
        }
        break;
      case OP_IF_ELSE_ZERO:
        {
          T zero_free = logic_or(L1>0, R1<0);
          T is_not_zero = logic_or(L1!=R1, L1!=0);
          L = is_not_zero*fmin(if_else_zero(zero_free, L2), L2);
          R = is_not_zero*fmax(if_else_zero(zero_free, R2), R2);
        }
        break;
      case OP_FABS:
        {
          T zero_free = logic_or(L1>0, R1<0);
          T is_not_zero = logic_or(L1!=R1, L1!=0);
          L = is_not_zero*fmin(if_else_zero(zero_free, L2), L2);
          R = is_not_zero*fmax(if_else_zero(zero_free, R2), R2);
        }
        break;
      case OP_COS:
        {
          T ge_2pi = ((R1-L1)>=2*M_PI);
          T L1_cos = cos(L1);
          T R1_cos = cos(R1);
          T lb = fmin(L1_cos, R1_cos);
          T ub = fmax(L1_cos, R1_cos);
          T contains_max = ceil(L1/(2*M_PI))*2*M_PI<=R1;
          T contains_min = ceil((L1-M_PI)/(2*M_PI))*2*M_PI+M_PI<=R1;
          L = if_else(logic_or(ge_2pi, contains_min), -1, lb);
          R = if_else(logic_or(ge_2pi, contains_max), 1, ub);
        }
        break;
      case OP_SIN:
        {
          T L1_mov = L1 - M_PI/2;
          T R1_mov = R1 - M_PI/2;
          my_propagate_interval(OP_COS, L1_mov, R1_mov, L2, R2, L, R); 
        }
        break;
      case OP_ASIN:
        {
          T out_of_domain = logic_or(R1<-1, L1>1);
          L = if_else(out_of_domain, nan, if_else(L1<=-1, -pi/2, asin(L1)));
          R = if_else(out_of_domain, nan, if_else(R1>=1, pi/2, asin(R1)));
        }
        break;
      case OP_ACOS:
        {
          T out_of_domain = logic_or(R1<-1, L1>1);
          L = if_else(out_of_domain, nan, if_else(R1>=1, 0, acos(R1)));
          R = if_else(out_of_domain, nan, if_else(L1<=-1, pi, acos(L1)));
        }
        break;
      case OP_MUL:
      {
        L = fmin(fmin(fmin(L1*L2, L1*R2), R1*L2), R1*R2);
        R = fmax(fmax(fmax(L1*L2, L1*R2), R1*L2), R1*R2);
        break;
      }
      case OP_SQ:
      {
          T zero_free = logic_or(L1>0, R1<0);
          L = if_else_zero(zero_free, if_else(L1>0, L1*L1, R1*R1));
          R = if_else(zero_free, if_else(L1>0, R1*R1, L1*L1), fmax(L1*L1, R1*R1));
        break;
      }
      default:
        casadi_warning("Not implemented: "+str(casadi_math<MX>::op_as_string(op)) );
    }
  }

  Function SXFunction::get_interval_propagator() const {
    std::vector< std::vector<SXElem> > arg_L, arg_R;
    for (casadi_int k=0;k<n_in_;++k) {
      std::vector<SXElem> in = in_[k].get_nonzeros();
      for (SXElem &e : in) {
        e = SXElem::sym(e.name() + "_L");
      }
      arg_L.push_back(in);
      in= in_[k].get_nonzeros();
      for (SXElem &e : in) {
        e = SXElem::sym(e.name() + "_R");
      }
      arg_R.push_back(in);
    }
    std::vector< std::vector<SXElem> > res_L, res_R;
    for (casadi_int k=0;k<n_out_;++k) {
      res_L.emplace_back(out_[k].nnz());
      res_R.emplace_back(out_[k].nnz());
    }

    // Iterator to stack of constants
    std::vector<SXElem>::const_iterator c_it = constants_.begin();

    // Iterator to free variables
    std::vector<SXElem>::const_iterator p_it = free_vars_.begin();
    std::vector<SXElem> w_L(worksize_);
    std::vector<SXElem> w_R(worksize_);

    // Evaluate algorithm
    if (verbose_) casadi_message("Evaluating algorithm forward");
    for (auto&& a : algorithm_) {
      switch (a.op) {
      case OP_INPUT:
        w_L[a.i0] = arg_L[a.i1][a.i2];
        w_R[a.i0] = arg_R[a.i1][a.i2];
        break;
      case OP_OUTPUT:
        res_L[a.i0][a.i2] = w_L[a.i1];
        res_R[a.i0][a.i2] = w_R[a.i1];
        break;
      case OP_CONST:
        w_L[a.i0] = *c_it++;
        w_R[a.i0] = w_L[a.i0];
        break;
      case OP_PARAMETER:
        w_L[a.i0] = *p_it++;
        w_R[a.i0] = w_L[a.i0];
        break;
      case OP_CALL:
        casadi_error("Not implemented");
        break;
      default:
        {
          // Do not store directly into w_L[a.i0] since a.i0 may be equal to a.i1 or a.i2
          SXElem L, R;
          my_propagate_interval(a.op, w_L[a.i1], w_R[a.i1],
                                                        w_L[a.i2], w_R[a.i2],
                                                        L, R);
          w_L[a.i0] = L;
          w_R[a.i0] = R;
        }
      }
    }

    std::vector<SX> args;

    for (casadi_int k=0;k<n_in_;++k) {
      args.push_back(SX(in_[k].sparsity(), arg_L[k]));
    }
    for (casadi_int k=0;k<n_in_;++k) {
      args.push_back(SX(in_[k].sparsity(), arg_R[k]));
    }
    std::vector<SX> exprs;
    for (casadi_int k=0;k<n_out_;++k) {
      exprs.push_back(SX(sparsity_out_[k], res_L[k]));
    }
    for (casadi_int k=0;k<n_out_;++k) {
      exprs.push_back(SX(sparsity_out_[k], res_R[k]));
    }


    std::vector<std::string> name_in;
    for (const std::string& e : name_in_) {
      name_in.push_back(e + "_L");
    }
    for (const std::string& e : name_in_) {
      name_in.push_back(e + "_R");
    }

    std::vector<std::string> name_out;
    for (const std::string& e : name_out_) {
      name_out.push_back(e + "_L");
    }
    for (const std::string& e : name_out_) {
      name_out.push_back(e + "_R");
    }

    return Function(name_ + "_interval_propagator", args, exprs, name_in, name_out);

  }


  void SXFunction::eval_mx(const MXVector& arg, MXVector& res,
                          bool always_inline, bool never_inline) const {
    always_inline = always_inline || always_inline_;
    never_inline = never_inline || never_inline_;

    // non-inlining call is implemented in the base-class
    if (!always_inline) {
      return FunctionInternal::eval_mx(arg, res, false, true);
    }

    if (verbose_) casadi_message(name_ + "::eval_mx");

    // Iterator to stack of constants
    std::vector<SXElem>::const_iterator c_it = constants_.begin();

    casadi_assert(!has_free(),
      "Free variables not supported in inlining call to SXFunction::eval_mx");

    // Resize the number of outputs
    casadi_assert(arg.size()==n_in_, "Wrong number of input arguments");
    res.resize(out_.size());

    // Symbolic work, non-differentiated
    std::vector<MX> w(sz_w());
    if (verbose_) casadi_message("Allocated work vector");

    // Split up inputs analogous to symbolic primitives
    std::vector<std::vector<MX> > arg_split(in_.size());
    for (casadi_int i=0; i<in_.size(); ++i) {
      // Get nonzeros of argument
      std::vector<MX> orig = arg[i].get_nonzeros();

      // Project to needed sparsity
      std::vector<MX> target(sparsity_in_[i].nnz(), 0);
      std::vector<MX> w(arg[i].size1());
      casadi_project(get_ptr(orig), arg[i].sparsity(),
                     get_ptr(target), sparsity_in_[i], get_ptr(w));

      // Store
      arg_split[i] = target;
    }

    // Allocate storage for split outputs
    std::vector<std::vector<MX> > res_split(out_.size());
    for (casadi_int i=0; i<out_.size(); ++i) res_split[i].resize(nnz_out(i));

    // Evaluate algorithm
    if (verbose_) casadi_message("Evaluating algorithm forward");
    for (auto&& a : algorithm_) {
      switch (a.op) {
      case OP_INPUT:
        w[a.i0] = arg_split[a.i1][a.i2];
        break;
      case OP_OUTPUT:
        res_split[a.i0][a.i2] = w[a.i1];
        break;
      case OP_CONST:
        w[a.i0] = static_cast<double>(*c_it++);
        break;
      case OP_CALL:
        {
          const ExtendedAlgEl& m = call_.el.at(a.i1);
          std::vector<MX> deps(m.n_dep);
          std::vector<MX> args;

          casadi_int k = 0;
          // Construct matrix-valued function arguments
          for (casadi_int i=0;i<m.f_n_in;++i) {
            std::vector<MX> arg;
            for (casadi_int j=0;j<m.f_nnz_in[i];++j) {
              arg.push_back(w[m.dep[k++]]);
            }
            args.push_back(sparsity_cast(vertcat(arg), m.f.sparsity_in(i)));
          }


          std::vector<MX> ret = m.f(args);
          std::vector<MX> res;

          // Break apart matriv-valued outputs into scalar components
          for (casadi_int i=0;i<m.f_n_out;++i) {
            std::vector<MX> nz = ret[i].get_nonzeros();
            res.insert(res.end(), nz.begin(), nz.end());
          }

          // Store into work vector
          for (casadi_int i=0;i<m.n_res;++i) {
            if (m.res[i]>=0) w[m.res[i]] = res[i];
          }
        }
        break;
      default:
        // Evaluate the function to a temporary value
        // (as it might overwrite the children in the work vector)
        MX f;
        switch (a.op) {
          CASADI_MATH_FUN_BUILTIN(w[a.i1], w[a.i2], f)
        }

        // Finally save the function value
        w[a.i0] = f;
      }
    }

    // Join split outputs
    for (casadi_int i=0; i<res.size(); ++i) {
      res[i] = sparsity_cast(vertcat(res_split[i]), sparsity_out_[i]);
    }
  }

  bool SXFunction::should_inline(bool with_sx, bool always_inline, bool never_inline) const {
    // If inlining has been specified
    casadi_assert(!(always_inline && never_inline),
      "Inconsistent options for " + definition());
    casadi_assert(!(never_inline && has_free()),
      "Must inline " + definition());
    if (always_inline) return true;
    if (never_inline) return false;
    // Functions with free variables must be inlined
    if (has_free()) return true;
    // Inlining by default
    return true;
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
            case OP_CALL:
              it1++->d[0] = f;
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
        case OP_CALL:
          {
            auto& m = call_.el.at(a.i1);

            // Construct forward sensitivity function
            Function ff = m.f.forward(1);

            // Symbolic inputs to forward sensitivity function
            std::vector<SXElem> deps;
            deps.reserve(2*m.n_dep);

            // Set nominal inputs from node
            CallSX* it2_node = static_cast<CallSX*>(it2->d[0].get());
            for (casadi_int i=0;i<m.n_dep;++i) {
              deps.push_back(it2_node->dep(i));
            }

            // Do not set nominal outputs
            for (casadi_int i=0;i<m.f_n_out;++i) {
              casadi_int nnz = ff.nnz_in(i+m.f_n_in);
              casadi_assert(nnz==0, "Not implemented");
            }

            // Read in forward seeds from work vector
            casadi_int offset = 0;
            for (casadi_int i=0;i<m.f_n_in;++i) {
              casadi_int nnz = ff.nnz_in(i+m.f_n_in+m.f_n_out);
              // nnz=0 occurs for is_diff_in[i] false
              casadi_assert(nnz==0 || nnz==m.f.nnz_in(i), "Not implemented");
              if (nnz) {
                for (casadi_int j=0;j<nnz;++j) {
                  deps.push_back(w[m.dep[offset+j]]);
                }
              }
              offset += m.f_nnz_in[i];
            }

            // Call forward sensitivity function
            std::vector<SXElem> ret = SXElem::call(ff, deps);

            // Retrieve sensitivities
            offset = 0;
            casadi_int k = 0;
            for (casadi_int i=0;i<m.f_n_out;++i) {
              casadi_int nnz = ff.nnz_out(i);
              // nnz=0 occurs for is_diff_out[i] false
              casadi_assert(nnz==0 || nnz==m.f_nnz_out[i], "Not implemented");
              if (nnz) {
                for (casadi_int j=0;j<nnz;++j) {
                  if (m.res[offset+j]>=0) w[m.res[offset+j]] = ret[k];
                  k++;
                }
              }
              offset += m.f_nnz_out[i];
            }
          }
          it2++;
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
            case OP_CALL:
              it1++->d[0] = f;
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
        case OP_CALL:
          {
            auto& m = call_.el.at(it->i1);

            // Construct reverse sensitivity function
            Function fr = m.f.reverse(1);

            // Symbolic inputs to reverse sensitivity function
            std::vector<SXElem> deps;
            deps.reserve(m.n_dep+m.n_res);

            // Set nominal inputs from node
            CallSX* it2_node = static_cast<CallSX*>(it2->d[0].get());
            for (casadi_int i=0;i<m.n_dep;++i) {
              deps.push_back(it2_node->dep(i));
            }

            // Do not set nominal outputs
            for (casadi_int i=0;i<m.f_n_out;++i) {
              casadi_int nnz = fr.nnz_in(i+m.f_n_in);
              casadi_assert(nnz==0, "Not implemented");
            }

            // Read in reverse seeds from work vector
            casadi_int offset = 0;
            for (casadi_int i=0;i<m.f_n_out;++i) {
              casadi_int nnz = fr.nnz_in(i+m.f_n_in+m.f_n_out);
              // nnz=0 occurs for is_diff_out[i] false
              casadi_assert(nnz==0 || nnz==m.f.nnz_out(i), "Not implemented");
              if (nnz) {
                for (casadi_int j=0;j<nnz;++j) {
                  deps.push_back((m.res[offset+j]>=0) ? w[m.res[offset+j]] : 0);
                }
              }
              offset += m.f.nnz_out(i);
            }

            // Call reverse sensitivity function
            std::vector<SXElem> ret = SXElem::call(fr, deps);

            // Clear out reverse seeds
            for (casadi_int i=0;i<m.n_res;++i) {
              if (m.res[i]>=0) w[m.res[i]] = 0;
            }

            // Store reverse sensitivities into work vector
            offset = 0;
            casadi_int k = 0;
            for (casadi_int i=0;i<m.f_n_in;++i) {
              casadi_int nnz = fr.nnz_out(i);
              // nnz=0 occurs for is_diff_in[i] false
              casadi_assert(nnz==0 || nnz==m.f_nnz_in[i], "Not implemented");
              if (nnz) {
                for (casadi_int j=0;j<nnz;++j) {
                  w[m.dep[offset+j]] += ret[k++];
                }
              }
              offset += m.f_nnz_in[i];
            }
          }
          it2++;
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

  template<typename T, typename CT>
  void SXFunction::call_setup(const ExtendedAlgEl& m,
    CT*** call_arg, T*** call_res, casadi_int** call_iw, T** call_w, T** nz_in, T** nz_out) const {
    *call_arg += n_in_;
    *call_res += n_out_;
    *nz_in    = *call_w + worksize_;
    *nz_out   = *call_w + worksize_ + call_.sz_w_arg;
    *call_w   = *call_w + worksize_ + call_.sz_w_arg + call_.sz_w_res;

    // Set up call_arg to point to nz_in
    T* ptr_w = *nz_in;
    for (casadi_int i=0;i<m.f_n_in;++i) {
      (*call_arg)[i] = ptr_w;
      ptr_w+=m.f_nnz_in[i];
    }

    // Set up call_res to point to nz_out
    ptr_w = *nz_out;
    for (casadi_int i=0;i<m.f_n_out;++i) {
      (*call_res)[i] = ptr_w;
      ptr_w+=m.f_nnz_out[i];
    }
  }

  template<typename T>
  void SXFunction::call_fwd(const AlgEl& e, const T** arg, T** res, casadi_int* iw, T* w) const {
    auto& m = call_.el[e.i1];
    const T** call_arg   = arg;
    T** call_res         = res;
    casadi_int* call_iw  = iw;
    T* call_w            = w;
    T* nz_in;
    T* nz_out;

    call_setup(m, &call_arg, &call_res, &call_iw, &call_w, &nz_in, &nz_out);

    // Populate nz_in from work vector
    for (casadi_int i=0;i<m.n_dep;++i) {
      nz_in[i] = w[m.dep[i]];
    }
    // Perform call nz_in -> nz_out
    m.f(call_arg, call_res, call_iw, call_w);

    // Store nz_out results back in workvector
    for (casadi_int i=0;i<m.n_res;++i) {
      // Only if the result is actually needed
      if (m.res[i]>=0) {
        w[m.res[i]] = nz_out[i];
      }
    }
  }


  template<typename T>
  void SXFunction::call_rev(const AlgEl& e, T** arg, T** res, casadi_int* iw, T* w) const {
    auto& m = call_.el[e.i1];
    bvec_t** call_arg = arg;
    bvec_t** call_res       = res;
    casadi_int* call_iw     = iw;
    bvec_t* call_w          = w;
    bvec_t* nz_in;
    bvec_t* nz_out;

    call_setup(m, &call_arg, &call_res, &call_iw, &call_w, &nz_in, &nz_out);

    std::fill_n(nz_in, m.n_dep, 0);

    // Read in reverse seeds nz_out from work vector
    for (casadi_int i=0;i<m.n_res;++i) {
      nz_out[i] = (m.res[i]>=0) ? w[m.res[i]] : 0;
    }

    // Perform reverse mode call nz_out -> nz_in
    m.f.rev(call_arg, call_res, call_iw, call_w);

    // Clear out reverse seeds
    for (casadi_int i=0;i<m.n_res;++i) {
      if (m.res[i]>=0) w[m.res[i]] = 0;
    }

    // Store reverse sensitivities into work vector
    for (casadi_int i=0;i<m.n_dep;++i) {
      w[m.dep[i]] |= nz_in[i];
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
      case OP_CALL:
        call_fwd(e, arg, res, iw, w);
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
      case OP_CALL:
        call_rev(*it, arg, res, iw, w);
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
    int version = s.version("SXFunction", 1, 2);
    size_t n_instructions;
    s.unpack("SXFunction::n_instr", n_instructions);

    s.unpack("SXFunction::worksize", worksize_);
    s.unpack("SXFunction::free_vars", free_vars_);
    s.unpack("SXFunction::operations", operations_);
    s.unpack("SXFunction::constants", constants_);
    s.unpack("SXFunction::default_in", default_in_);

    if (version>=2) {

      s.unpack("SXFunction::call_sz_arg", call_.sz_arg);
      s.unpack("SXFunction::call_sz_res", call_.sz_res);
      s.unpack("SXFunction::call_sz_iw", call_.sz_iw);
      s.unpack("SXFunction::call_sz_w", call_.sz_w);
      s.unpack("SXFunction::call_sz_arg", call_.sz_w_arg);
      s.unpack("SXFunction::call_sz_res", call_.sz_w_res);

      size_t el_size;
      s.unpack("SXFunction::call_el_size", el_size);
      call_.el.reserve(el_size);

      // Loop over nodes
      for (casadi_int k=0;k<el_size;++k) {
        Function f;
        s.unpack("SXFunction::call_el_f", f);
        call_.el.emplace_back(f);
        auto& e = call_.el[k];
        s.unpack("SXFunction::call_el_dep", e.dep);
        s.unpack("SXFunction::call_el_res", e.res);
        s.unpack("SXFunction::call_el_copy_elision_arg", e.copy_elision_arg);
        s.unpack("SXFunction::call_el_copy_elision_offset", e.copy_elision_offset);
      }

      s.unpack("SXFunction::copy_elision", copy_elision_);

    } else {
      call_.sz_arg = 0;
      call_.sz_res = 0;
      call_.sz_iw = 0;
      call_.sz_w = 0;
      call_.sz_w_arg = 0;
      call_.sz_w_res = 0;
      call_.el.clear();
      copy_elision_.resize(n_instructions, false);
    }

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
    s.version("SXFunction", 2);
    s.pack("SXFunction::n_instr", algorithm_.size());

    s.pack("SXFunction::worksize", worksize_);
    s.pack("SXFunction::free_vars", free_vars_);
    s.pack("SXFunction::operations", operations_);
    s.pack("SXFunction::constants", constants_);
    s.pack("SXFunction::default_in", default_in_);

    s.pack("SXFunction::call_sz_arg", call_.sz_arg);
    s.pack("SXFunction::call_sz_res", call_.sz_res);
    s.pack("SXFunction::call_sz_iw", call_.sz_iw);
    s.pack("SXFunction::call_sz_w", call_.sz_w);
    s.pack("SXFunction::call_sz_arg", call_.sz_w_arg);
    s.pack("SXFunction::call_sz_res", call_.sz_w_res);

    s.pack("SXFunction::call_el_size", call_.el.size());
    // Loop over ExtendedALgEl elements
    for (const auto& n : call_.el) {
      s.pack("SXFunction::call_el_f", n.f);
      s.pack("SXFunction::call_el_dep", n.dep);
      s.pack("SXFunction::call_el_res", n.res);
      s.pack("SXFunction::call_el_copy_elision_arg", n.copy_elision_arg);
      s.pack("SXFunction::call_el_copy_elision_offset", n.copy_elision_offset);
    }

    s.pack("SXFunction::copy_elision", copy_elision_);

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

  void SXFunction::find(std::map<FunctionInternal*, Function>& all_fun,
      casadi_int max_depth) const {
    for (auto&& e : algorithm_) {
      if (e.op == OP_CALL) {
        const ExtendedAlgEl& m = call_.el.at(e.i1);
        add_embedded(all_fun, m.f, max_depth);
      }
    }
  }

  std::vector<SX> SXFunction::order(const std::vector<SX>& expr) {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    std::lock_guard<std::mutex> lock(SX::get_mutex_temp());
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
    // Stack used to sort the computational graph
    std::stack<SXNode*> s;

    // All nodes
    std::vector<SXNode*> nodes;

    // Add the list of nodes
    casadi_int ind=0;
    for (auto it = expr.begin(); it != expr.end(); ++it, ++ind) {
      casadi_int nz=0;
      for (auto itc = (*it)->begin(); itc != (*it)->end(); ++itc, ++nz) {
        // Add outputs to the list
        s.push(itc->get());
        XFunction<SXFunction, SX, SXNode>::sort_depth_first(s, nodes);
      }
    }

    // Clear temporary markers
    for (casadi_int i=0; i<nodes.size(); ++i) {
      nodes[i]->temp = 0;
    }

    std::vector<SX> ret(nodes.size());
    for (casadi_int i=0; i<nodes.size(); ++i) {
      ret[i] = SXElem::create(nodes[i]);
    }

    return ret;
  }

} // namespace casadi
