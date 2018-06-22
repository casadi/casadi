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


#include "sx_function.hpp"
#include <limits>
#include <stack>
#include <deque>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <bitset>
#include "casadi_misc.hpp"
#include "sx_node.hpp"
#include "output_sx.hpp"
#include "call_sx.hpp"
#include "casadi_common.hpp"
#include "sparsity_internal.hpp"
#include "global_options.hpp"
#include "casadi_interrupt.hpp"
#include "serializer.hpp"

namespace casadi {

  using namespace std;


  SXFunction::CallInfo::Node::Node(const Function& fun) : f(fun) {
    n_dep = f.nnz_in(); n_out = f.nnz_out();
    dep.resize(n_dep); out.resize(n_out, -1); out_sx.resize(n_out, 0);
    f_n_in = f.n_in(); f_n_out = f.n_out();
    f_nnz_in.resize(f_n_in); f_nnz_out.resize(f_n_out);
    for (casadi_int i=0;i<f_n_in;++i) f_nnz_in[i] = f.nnz_in(i);
    for (casadi_int i=0;i<f_n_out;++i) f_nnz_out[i] = f.nnz_out(i);
  }

  SXFunction::SXFunction(const std::string& name,
                         const vector<SX >& inputv,
                         const vector<SX >& outputv,
                         const vector<std::string>& name_in,
                         const vector<std::string>& name_out)
    : XFunction<SXFunction, SX, SXNode>(name, inputv, outputv, name_in, name_out) {

    // Default (persistent) options
    just_in_time_opencl_ = false;
    just_in_time_sparsity_ = false;
  }

  SXFunction::~SXFunction() {
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
        {
          const double** call_arg = arg + n_in_;
          double** call_res       = res + n_out_;
          casadi_int* call_iw     = iw;
          double* call_w          = w + worksize_;
          double* call_w_arg      = call_w + call_.sz_w;
          double* call_w_res      = call_w_arg + call_.sz_w_arg;

          auto& m = call_.nodes[e.i1];
          double* ptr_w = call_w_arg;
          for (casadi_int i=0;i<m.f_n_in;++i) {
            call_arg[i] = ptr_w;
            ptr_w+=m.f_nnz_in[i];
          }
          ptr_w = call_w_res;
          for (casadi_int i=0;i<m.f_n_out;++i) {
            call_res[i] = ptr_w;
            ptr_w+=m.f_nnz_out[i];
          }
          for (casadi_int i=0;i<m.n_dep;++i) call_w_arg[i] = w[m.dep[i]];
          // TODO(jgillis): should we do a checkout upfront?
          m.f(call_arg, call_res, call_iw, call_w);
          for (casadi_int i=0;i<m.n_out;++i) {
            if (m.out[i]>=0) w[m.out[i]] = call_w_res[i];
          }
        }
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

  void SXFunction::disp_more(ostream &stream) const {
    stream << "Algorithm:";

    // Iterator to free variables
    vector<SXElem>::const_iterator p_it = free_vars_.begin();

    // Normal, interpreted output
    for (auto&& a : algorithm_) {
      InterruptHandler::check();
      stream << endl;
      if (a.op==OP_OUTPUT) {
        stream << "output[" << a.i0 << "][" << a.i2 << "] = @" << a.i1;
      } else if (a.op==OP_CALL) {
        auto& m = call_.nodes.at(a.i1);
        stream << "[";
        for (casadi_int i=0;i<m.n_out;++i) {
          stream << "@" << m.out[i];
          if (i < m.n_out-1) stream << ",";
        }
        stream << "] = ";
        stream << m.f.name() << "(";
        for (casadi_int i=0;i<m.n_dep;++i) {
          stream << "@" << m.dep[i];
          if (i < m.n_dep-1) stream << ",";
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

  void SXFunction::codegen_declarations(CodeGenerator& g) const {

    // Make sure that there are no free variables
    if (!free_vars_.empty()) {
      casadi_error("Code generation is not possible since variables "
                   + str(free_vars_) + " are free.");
    }

    // Generate code for the call nodes
    for (auto&& m : call_.nodes) {
      g.add_dependency(m.f);
    }
  }

  void SXFunction::codegen_body(CodeGenerator& g) const {

    // Run the algorithm
    for (auto&& a : algorithm_) {
      if (a.op==OP_OUTPUT) {
        g << "if (res[" << a.i0 << "]!=0) "
          << "res["<< a.i0 << "][" << a.i2 << "]=" << g.sx_work(a.i1) << ";\n";
      } else if (a.op==OP_CALL) {
        auto& m = call_.nodes[a.i1];

        // Collect input arguments
        casadi_int offset = worksize_+call_.sz_w;
        for (casadi_int i=0; i<m.f_n_in; ++i) {
          g << "arg[" << n_in_+i << "]=" << "w+" + str(i+offset) << ";\n";
          offset += m.f_nnz_in[i];
        }

        offset = worksize_+call_.sz_w+call_.sz_w_arg;
        // Collect output arguments
        for (casadi_int i=0; i<m.f_n_out; ++i) {
          g << "res[" << n_out_+i << "]=" << "w+" + str(i+offset) << ";\n";
          offset += m.f_nnz_out[i];
        }
        for (casadi_int i=0;i<m.n_dep;++i) {
          g << "w["+str(i+worksize_+call_.sz_w) + "] = " << g.sx_work(m.dep[i]) << ";\n";
        }
        // Call function
        g << "if (";
        g << g(m.f, "arg+"+str(n_in_), "res+"+str(n_out_), "iw", "w+" + str(worksize_));
        g << ") return 1;\n";
        for (casadi_int i=0;i<m.n_out;++i) {
          if (m.out[i]>=0)
            g << g.sx_work(m.out[i]) << " = ";
            g << "w[" + str(i+worksize_+call_.sz_w+call_.sz_w_arg) + "];\n";
        }
      } else {

        // Where to store the result
        g << g.sx_work(a.i0) << "=";

        // What to store
        if (a.op==OP_CONST) {
          g << g.constant(a.d);
        } else if (a.op==OP_INPUT) {
          g << "arg[" << a.i1 << "] ? arg[" << a.i1 << "][" << a.i2 << "] : 0";
        } else {
          casadi_int ndep = casadi_math<double>::ndeps(a.op);
          casadi_assert_dev(ndep>0);
          if (ndep==1) g << g.print_op(a.op, g.sx_work(a.i1));
          if (ndep==2) g << g.print_op(a.op, g.sx_work(a.i1), g.sx_work(a.i2));
        }
        g  << ";\n";
      }
    }
  }

  Options SXFunction::options_
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
        "Reuse variables in the work vector"}}
     }
  };

  void SXFunction::init(const Dict& opts) {
    // Call the init function of the base class
    XFunction<SXFunction, SX, SXNode>::init(opts);
    if (verbose_) casadi_message(name_ + "::init");

    // Default (temporary) options
    bool live_variables = true;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="default_in") {
        default_in_ = op.second;
      } else if (op.first=="live_variables") {
        live_variables = op.second;
      } else if (op.first=="just_in_time_opencl") {
        just_in_time_opencl_ = op.second;
      } else if (op.first=="just_in_time_sparsity") {
        just_in_time_sparsity_ = op.second;
      }
    }

    // Check/set default inputs
    if (default_in_.empty()) {
      default_in_.resize(n_in_, 0);
    } else {
      casadi_assert(default_in_.size()==n_in_,
                            "Option 'default_in' has incorrect length");
    }

    // Stack used to sort the computational graph
    stack<SXNode*> s;

    // All nodes
    vector<SXNode*> nodes;

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
    for (vector<SXNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
      SXNode* t = *it;
      if (t) {
        if (t->is_constant())
          constants_.push_back(SXElem::create(t));
        else if (!t->is_symbolic() && t->op()>=0)
          operations_.push_back(SXElem::create(t));
      }
    }

    // Input instructions
    vector<pair<int, SXNode*> > symb_loc;

    // Current output and nonzero, start with the first one
    int curr_oind, curr_nz=0;
    casadi_assert(out_.size() <= std::numeric_limits<int>::max(), "Integer overflow");
    for (curr_oind=0; curr_oind<out_.size(); ++curr_oind) {
      if (out_[curr_oind].nnz()!=0) {
        break;
      }
    }

    // Count the number of times each node is used
    vector<casadi_int> refcount(nodes.size(), 0);

    // Get the sequence of instructions for the virtual machine
    algorithm_.resize(0);
    algorithm_.reserve(nodes.size());

    // Mapping of node index (cfr. temp) to algorithm index
    std::vector<int> alg_index;
    alg_index.reserve(nodes.size());

    for (vector<SXNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {

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
        symb_loc.push_back(make_pair(algorithm_.size(), n));
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

          // Reserve memory for call node
          ae.i1 = call_.nodes.size();

          const Function& f = static_cast<const CallSX*>(n)->f_;
          call_.nodes.emplace_back(f);
          auto& m = call_.nodes.at(ae.i1);

          call_.sz_arg = std::max(call_.sz_arg, f.sz_arg());
          call_.sz_res = std::max(call_.sz_res, f.sz_res());
          call_.sz_iw  = std::max(call_.sz_iw, f.sz_iw());
          call_.sz_w   = std::max(call_.sz_w, f.sz_w());
          call_.sz_w_arg   = std::max(call_.sz_w_arg, static_cast<size_t>(f.nnz_in()));
          call_.sz_w_res   = std::max(call_.sz_w_res,  static_cast<size_t>(f.nnz_out()));

          dep = get_ptr(m.dep);
          ndeps = m.n_dep;
          for (casadi_int i=0; i<ndeps; ++i) {
            dep[i] = n->dep(i).get()->temp;
          }
        }
        break;
      case -1: // Output extraction node
        {
          dep = &algorithm_.at(alg_index.at(n->dep(0).get()->temp)).i1;
          int oind = static_cast<OutputSX*>(n)->oind_;
          casadi_assert(call_.nodes.at(dep[0]).out.at(oind)==-1, "Duplicate");
          call_.nodes.at(dep[0]).out.at(oind) = n->temp;
          call_.nodes.at(dep[0]).out_sx.at(oind) = SXElem(n, false);
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
    vector<int> place(nodes.size());

    // Stack with unused elements in the work vector
    stack<int> unused;

    // Work vector size
    int worksize = 0;

    // Find a place in the work vector for the operation
    for (auto&& a : algorithm_) {

      // Default dependencies
      int* dep = &a.i1;
      casadi_int ndeps = casadi_math<double>::ndeps(a.op);

      // Default outputs
      int* out = &a.i0;
      casadi_int nout = 1;

      if (a.op==OP_CALL) { // Call node overrides these defaults
        auto& e = call_.nodes.at(a.i1);
        ndeps = e.n_dep;
        dep = get_ptr(e.dep);
        nout= e.n_out;
        out = get_ptr(e.out);
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
        for (casadi_int c=0; c<nout; ++c) {
          if (out[c]<0) continue;
          if (live_variables && !unused.empty()) {
            // Try to reuse a variable from the stack if possible (last in, first out)
            out[c] = place[out[c]] = unused.top();
            unused.pop();
          } else {
            // Allocate a new variable
            out[c] = place[out[c]] = worksize++;
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
      if (live_variables) {
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
    for (vector<pair<int, SXNode*> >::const_iterator it=symb_loc.begin();
         it!=symb_loc.end(); ++it) {
      if (it->second->temp!=0) {
        // Save to list of free parameters
        free_vars_.push_back(SXElem::create(it->second));

        // Remove marker
        it->second->temp=0;
      }
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

  int SXFunction::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w, void* mem) const {
    if (verbose_) casadi_message(name_ + "::eval_sx");

    // Iterator to the binary operations
    vector<SXElem>::const_iterator b_it=operations_.begin();

    // Iterator to stack of constants
    vector<SXElem>::const_iterator c_it = constants_.begin();

    // Iterator to free variables
    vector<SXElem>::const_iterator p_it = free_vars_.begin();

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
          auto& m = call_.nodes.at(a.i1);
          const SXElem& orig = *b_it++;
          std::vector<SXElem> deps(m.n_dep);
          bool identical = true;

          std::vector<SXElem> ret;
          for (casadi_int i=0;i<m.n_dep;++i) {
            identical &= SXElem::is_equal(w[m.dep.at(i)], orig->dep(i), 2);
          }
          if (identical) {
            ret = OutputSX::split(orig, m.n_out);
            for (casadi_int i=0;i<m.n_out;++i) {
              if (!m.out_sx.at(i).is_constant()) {
                ret[i] = m.out_sx.at(i);
              }
            }
          } else {
            for (casadi_int i=0;i<m.n_dep;++i) deps[i] = w[m.dep[i]];
            ret = SXElem::call_fun(m.f, deps);
          }
          for (casadi_int i=0;i<m.n_out;++i) {
            if (m.out[i]>=0) w[m.out[i]] = ret[i];
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

  void SXFunction::ad_forward(const vector<vector<SX> >& fseed,
                                vector<vector<SX> >& fsens) const {
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
          vector<vector<SX> > fseed2(fseed);
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
    vector<SXElem>::const_iterator b_it=operations_.begin();

    // Tape
    vector<TapeEl<SXElem> > s_pdwork(operations_.size());
    vector<TapeEl<SXElem> >::iterator it1 = s_pdwork.begin();

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
    vector<SXElem> w(worksize_);

    // Calculate forward sensitivities
    if (verbose_) casadi_message("Calculating forward derivatives");
    for (casadi_int dir=0; dir<nfwd; ++dir) {
      vector<TapeEl<SXElem> >::const_iterator it2 = s_pdwork.begin();
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
        case OP_CALL:
          {
            auto& m = call_.nodes.at(a.i1);
            CallSX* e = static_cast<CallSX*>(it2->d[0].get());
            Function ff = m.f.forward(1);
            std::vector<SXElem> deps;
            // Add nominal input SXElem
            for (casadi_int i=0;i<m.n_dep;++i) {
              deps.push_back(e->dep(i));
            }
            //for (casadi_int i=0;i<m.n_out;++i) {
            //  deps.push_back(w[m.out[i]]);
            //}
            for (casadi_int i=0;i<m.n_dep;++i) {
              deps.push_back(w[m.dep[i]]);
            }
            std::vector<SXElem> ret = SXElem::call_fun(ff, deps);
            // Set resulting dot variables
            for (casadi_int i=0;i<m.n_out;++i) {
              if (m.out[i]>=0) w[m.out[i]] = ret[i];
            }
          }
          it2++;
          break;
          CASADI_MATH_BINARY_BUILTIN // Binary operation
            w[a.i0] = it2->d[0] * w[a.i1] + it2->d[1] * w[a.i2];it2++;break;
        default: // Unary operation
          w[a.i0] = it2->d[0] * w[a.i1]; it2++;
        }
      }
    }
  }

  void SXFunction::ad_reverse(const vector<vector<SX> >& aseed,
                                vector<vector<SX> >& asens) const {
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
      vector<vector<SX> > aseed2(aseed);
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
          fill(asens[d][i]->begin(), asens[d][i]->end(), 0);
        }
      }
    }

    // Iterator to the binary operations
    vector<SXElem>::const_iterator b_it=operations_.begin();

    // Tape
    vector<TapeEl<SXElem> > s_pdwork(operations_.size());
    vector<TapeEl<SXElem> >::iterator it1 = s_pdwork.begin();

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
    vector<SXElem> w(worksize_, 0);

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
        case OP_CALL:
          {
            auto& m = call_.nodes.at(it->i1);
            CallSX* e = static_cast<CallSX*>(it2->d[0].get());
            Function fr = m.f.reverse(1);
            std::vector<SXElem> deps;
            // Add nominal input SXElem
            for (casadi_int i=0;i<m.n_dep;++i) {
              deps.push_back(e->dep(i));
            }
            //for (casadi_int i=0;i<m.n_out;++i) {
            //  deps.push_back(w[m.out[i]]);
            //}
            for (casadi_int i=0;i<m.n_out;++i) {
              if (m.out[i]>=0) {
                deps.push_back(w[m.out[i]]);
                w[m.out[i]] = 0;
              } else {
                deps.push_back(0);
              }
            }
            std::vector<SXElem> ret = SXElem::call_fun(fr, deps);
            // Set resulting dot variables
            for (casadi_int i=0;i<m.n_dep;++i) {
              w[m.dep[i]] += ret[i];
            }
          }
          it2++;
          break;
          CASADI_MATH_BINARY_BUILTIN // Binary operation
            seed = w[it->i0];
          w[it->i0] = 0;
          w[it->i1] += it2->d[0] * seed;
          w[it->i2] += it2->d[1] * seed;
          it2++;
          break;
        default: // Unary operation
          seed = w[it->i0];
          w[it->i0] = 0;
          w[it->i1] += it2->d[0] * seed;
          it2++;
        }
      }
    }
  }

  int SXFunction::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
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
        {
          const bvec_t** call_arg = arg + n_in_;
          bvec_t** call_res       = res + n_out_;
          casadi_int* call_iw     = iw;
          bvec_t* call_w          = w + worksize_;
          bvec_t* call_w_arg      = call_w + call_.sz_w;
          bvec_t* call_w_res      = call_w_arg + call_.sz_w_arg;

          auto& m = call_.nodes[e.i1];
          bvec_t* ptr_w = call_w_arg;
          for (casadi_int i=0;i<m.f_n_in;++i) {
            call_arg[i] = ptr_w;
            ptr_w+=m.f_nnz_in[i];
          }
          ptr_w = call_w_res;
          for (casadi_int i=0;i<m.f_n_out;++i) {
            call_res[i] = ptr_w;
            ptr_w+=m.f_nnz_out[i];
          }
          for (casadi_int i=0;i<m.n_dep;++i) call_w_arg[i] = w[m.dep[i]];
          m.f(call_arg, call_res, call_iw, call_w);
          for (casadi_int i=0;i<m.n_out;++i) {
            if (m.out[i]>=0) w[m.out[i]] = call_w_res[i];
          }
        }
        break;
      default: // Unary or binary operation
        w[e.i0] = w[e.i1] | w[e.i2]; break;
      }
    }
    return 0;
  }

  int SXFunction::sp_reverse(bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    fill_n(w, sz_w(), 0);

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
        {
          bvec_t** call_arg = arg + n_in_;
          bvec_t** call_res       = res + n_out_;
          casadi_int* call_iw     = iw;
          bvec_t* call_w          = w + worksize_;
          bvec_t* call_w_arg      = call_w + call_.sz_w;
          bvec_t* call_w_res      = call_w_arg + call_.sz_w_arg;

          auto& m = call_.nodes[it->i1];
          bvec_t* ptr_w = call_w_arg;
          for (casadi_int i=0;i<m.f_n_in;++i) {
            call_arg[i] = ptr_w;
            ptr_w+=m.f_nnz_in[i];
          }
          ptr_w = call_w_res;
          for (casadi_int i=0;i<m.f_n_out;++i) {
            call_res[i] = ptr_w;
            ptr_w+=m.f_nnz_out[i];
          }

          fill_n(call_w_arg, m.n_dep, 0);
          for (casadi_int i=0;i<m.n_out;++i) {
            call_w_res[i] = (m.out[i]>=0) ? w[m.out[i]] : 0;
          }
          m.f.rev(call_arg, call_res, call_iw, call_w);

          for (casadi_int i=0;i<m.n_out;++i) {
            if (m.out[i]>=0) w[m.out[i]] = 0;
          }
          for (casadi_int i=0;i<m.n_dep;++i) w[m.dep[i]] |= call_w_arg[i];
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

  Function SXFunction::get_jacobian(const std::string& name,
                                       const std::vector<std::string>& inames,
                                       const std::vector<std::string>& onames,
                                       const Dict& opts) const {
    // Jacobian expression
    SX J = SX::jacobian(veccat(out_), veccat(in_));

    // All inputs of the return function
    std::vector<SX> ret_in(inames.size());
    copy(in_.begin(), in_.end(), ret_in.begin());
    for (casadi_int i=0; i<n_out_; ++i) {
      ret_in.at(n_in_+i) = SX::sym(inames[n_in_+i], Sparsity(out_.at(i).size()));
    }

    // Assemble function and return
    return Function(name, ret_in, {J}, inames, onames, opts);
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
    std::string indent = "";
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

  SXFunction::SXFunction(const Info& e) :
    XFunction<SXFunction, SX, SXNode>(e.xfunction),
    algorithm_(e.algorithm), worksize_(e.worksize),
    free_vars_(e.free_vars), operations_(e.operations),
    constants_(e.constants), default_in_(e.default_in), call_(e.call) {

    // Default (persistent) options
    just_in_time_opencl_ = false;
    just_in_time_sparsity_ = false;
  }

  void SXFunction::serialize(Serializer &s) const {
    FunctionInternal::serialize(s);
    s.pack("SXFunction::n_instr", casadi_int(algorithm_.size()));

    s.pack("SXFunction::worksize", casadi_int(worksize_));
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

    s.pack("SXFunction::call_nodes_size", call_.nodes.size());
    // Loop over nodes
    for (const auto& n : call_.nodes) {
      s.pack("SXFunction::call_nodes_f", n.f);
      s.pack("SXFunction::call_nodes_dep", n.dep);
      s.pack("SXFunction::call_nodes_out", n.out);
      s.pack("SXFunction::call_nodes_out_sx", n.out_sx);
    }

    // Loop over algorithm
    for (const auto& e : algorithm_) {
      s.pack("SXFunction::ScalarAtomic::op", e.op);
      s.pack("SXFunction::ScalarAtomic::i0", e.i0);
      s.pack("SXFunction::ScalarAtomic::i1", e.i1);
      s.pack("SXFunction::ScalarAtomic::i2", e.i2);
    }

    s.pack(in_);
    s.pack(out_);
  }

  Function SXFunction::deserialize(DeSerializer& s) {
    Info info;
    FunctionInternal::deserialize(s, info.xfunction.function);
    casadi_int n_instructions;
    s.unpack("SXFunction::n_instr", n_instructions);

    s.unpack("SXFunction::worksize", info.worksize);
    s.unpack("SXFunction::free_vars", info.free_vars);
    s.unpack("SXFunction::operations", info.operations);
    s.unpack("SXFunction::constants", info.constants);
    s.unpack("SXFunction::default_in", info.default_in);

    s.unpack("SXFunction::call_sz_arg", info.call.sz_arg);
    s.unpack("SXFunction::call_sz_res", info.call.sz_res);
    s.unpack("SXFunction::call_sz_iw", info.call.sz_iw);
    s.unpack("SXFunction::call_sz_w", info.call.sz_w);
    s.unpack("SXFunction::call_sz_arg", info.call.sz_w_arg);
    s.unpack("SXFunction::call_sz_res", info.call.sz_w_res);

    size_t nodes_size;
    s.unpack("SXFunction::call_nodes_size", nodes_size);
    info.call.nodes.reserve(nodes_size);

    // Loop over nodes
    for (casadi_int k=0;k<nodes_size;++k) {
      Function f;
      s.unpack("SXFunction::call_nodes_f", f);
      info.call.nodes.emplace_back(f);
      auto& e = info.call.nodes[k];
      s.unpack("SXFunction::call_nodes_dep", e.dep);
      s.unpack("SXFunction::call_nodes_out", e.out);
      s.unpack("SXFunction::call_nodes_out_sx", e.out_sx);
    }

    info.algorithm.resize(n_instructions);
    for (casadi_int k=0;k<n_instructions;++k) {
      AlgEl& e = info.algorithm[k];
      s.unpack("SXFunction::ScalarAtomic::op", e.op);
      s.unpack("SXFunction::ScalarAtomic::i0", e.i0);
      s.unpack("SXFunction::ScalarAtomic::i1", e.i1);
      s.unpack("SXFunction::ScalarAtomic::i2", e.i2);
    }

    s.unpack(info.xfunction.in);
    s.unpack(info.xfunction.out);

    Function ret;
    ret.own(new SXFunction(info));
    ret->finalize();
    return ret;
  }

} // namespace casadi
