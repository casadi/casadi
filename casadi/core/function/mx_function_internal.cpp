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

#include "mx_function_internal.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"

#include "../std_vector_tools.hpp"
#include "../casadi_types.hpp"

#include <stack>
#include <typeinfo>
#include "../profiling.hpp"
#include "../casadi_options.hpp"

using namespace std;

namespace casadi {

  MXFunctionInternal::MXFunctionInternal(const std::vector<MX>& inputv,
                                         const std::vector<MX>& outputv) :
    XFunctionInternal<MXFunction, MXFunctionInternal, MX, MXNode>(inputv, outputv) {

    setOption("name", "unnamed_mx_function");
  }


  MXFunctionInternal::~MXFunctionInternal() {
  }


  void MXFunctionInternal::init() {
    log("MXFunctionInternal::init begin");

    // Call the init function of the base class
    XFunctionInternal<MXFunction, MXFunctionInternal, MX, MXNode>::init();

    // Stack used to sort the computational graph
    stack<MXNode*> s;

    // All nodes
    vector<MXNode*> nodes;

    // Add the list of nodes
    int ind=0;
    for (vector<MX>::iterator it = outputv_.begin(); it != outputv_.end(); ++it, ++ind) {
      // Add outputs to the list
      s.push(static_cast<MXNode*>(it->get()));
      sort_depth_first(s, nodes);

      // A null pointer means an output instruction
      nodes.push_back(static_cast<MXNode*>(0));
    }

    // Set the temporary variables to be the corresponding place in the sorted graph
    for (int i=0; i<nodes.size(); ++i) {
      if (nodes[i]) {
        nodes[i]->temp = i;
      }
    }

    // Place in the algorithm for each node
    vector<int> place_in_alg;
    place_in_alg.reserve(nodes.size());

    // Use live variables?
    bool live_variables = getOption("live_variables");

    // Input instructions
    vector<pair<int, MXNode*> > symb_loc;

    // Current output and nonzero, start with the first one
    int curr_oind=0;

    // Count the number of times each node is used
    vector<int> refcount(nodes.size(), 0);

    // Maximum number of input and output arguments
    max_arg_ = 0;
    max_res_ = 0;

    // Get the sequence of instructions for the virtual machine
    algorithm_.resize(0);
    algorithm_.reserve(nodes.size());
    for (vector<MXNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
      // Current node
      MXNode* n = *it;

      // Get the operation
      int op = n==0 ? OP_OUTPUT : n->getOp();

      // Store location if parameter (or input)
      if (op==OP_PARAMETER) {
        symb_loc.push_back(make_pair(algorithm_.size(), n));
      }

      // If a new element in the algorithm needs to be added
      if (op>=0) {
        AlgEl ae;
        ae.op = op;
        ae.data.assignNode(n);

        // Add input and output argument
        if (op==OP_OUTPUT) {
          ae.arg.resize(1);
          ae.arg[0] = outputv_.at(curr_oind)->temp;
          ae.res.resize(1);
          ae.res[0] = curr_oind++;
        } else {
          ae.arg.resize(n->ndep());
          for (int i=0; i<n->ndep(); ++i) {
            ae.arg[i] = n->dep(i)->temp;
          }
          ae.res.resize(n->nout());
          if (n->isMultipleOutput()) {
            fill(ae.res.begin(), ae.res.end(), -1);
          } else {
            ae.res[0] = n->temp;
          }
        }

        // Update maximum size of arg and res
        max_arg_ = std::max(max_arg_, ae.arg.size());
        max_res_ = std::max(max_res_, ae.res.size());

        // Increase the reference count of the dependencies
        for (int c=0; c<ae.arg.size(); ++c) {
          if (ae.arg[c]>=0) {
            refcount[ae.arg[c]]++;
          }
        }

        // Save to algorithm
        place_in_alg.push_back(algorithm_.size());
        algorithm_.push_back(ae);

      } else { // Function output node
        // Get the output index
        int oind = n->getFunctionOutput();

        // Get the index of the parent node
        int pind = place_in_alg[n->dep(0)->temp];

        // Save location in the algorithm element corresponding to the parent node
        int& otmp = algorithm_[pind].res.at(oind);
        if (otmp<0) {
          otmp = n->temp; // First time this function output is encountered, save to algorithm
        } else {
          n->temp = otmp; // Function output is a duplicate, use the node encountered first
        }

        // Not in the algorithm
        place_in_alg.push_back(-1);
      }
    }

    // Place in the work vector for each of the nodes in the tree (overwrites the reference counter)
    vector<int>& place = place_in_alg; // Reuse memory as it is no longer needed
    place.resize(nodes.size());

    // Stack with unused elements in the work vector, sorted by sparsity pattern
    SPARSITY_MAP<int, stack<int> > unused_all;

    // Work vector size
    int worksize = 0;

    // Find a place in the work vector for the operation
    for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it) {

      // There are two tasks, allocate memory of the result and free the
      // memory off the arguments, order depends on whether inplace is possible
      int first_to_free = 0;
      int last_to_free = it->op==OP_OUTPUT ? 1 : it->data->numInplace();
      for (int task=0; task<2; ++task) {

        // Dereference or free the memory of the arguments
        for (int c=last_to_free-1; c>=first_to_free; --c) { // reverse order so that the
                                                          // first argument will end up
                                                          // at the top of the stack

          // Index of the argument
          int& ch_ind = it->arg[c];
          if (ch_ind>=0) {

            // Decrease reference count and add to the stack of
            // unused variables if the count hits zero
            int remaining = --refcount[ch_ind];

            // Free variable for reuse
            if (live_variables && remaining==0) {

              // Get a pointer to the sparsity pattern of the argument that can be freed
              int nnz = nodes[ch_ind]->sparsity().nnz();

              // Add to the stack of unused work vector elements for the current sparsity
              unused_all[nnz].push(place[ch_ind]);
            }

            // Point to the place in the work vector instead of to the place in the list of nodes
            ch_ind = place[ch_ind];
          }
        }

        // Nothing more to allocate
        if (it->op==OP_OUTPUT || task==1) break;

        // Free the rest in the next iteration
        first_to_free = last_to_free;
        last_to_free = it->arg.size();

        // Allocate/reuse memory for the results of the operation
        for (int c=0; c<it->res.size(); ++c) {
          if (it->res[c]>=0) {

            // Are reuse of variables (live variables) enabled?
            if (live_variables) {
              // Get a pointer to the sparsity pattern node
              int nnz = it->data->sparsity(c).nnz();

              // Get a reference to the stack for the current sparsity
              stack<int>& unused = unused_all[nnz];

              // Try to reuse a variable from the stack if possible (last in, first out)
              if (!unused.empty()) {
                it->res[c] = place[it->res[c]] = unused.top();
                unused.pop();
                continue; // Success, no new element needed in the work vector
              }
            }

            // Allocate a new element in the work vector
            it->res[c] = place[it->res[c]] = worksize++;
          }
        }
      }
    }

    if (verbose()) {
      if (live_variables) {
        cout << "Using live variables: work array is "
             <<  worksize << " instead of "
             << nodes.size() << endl;
      } else {
        cout << "Live variables disabled." << endl;
      }
    }

    // Allocate work vectors (numeric)
    workloc_.resize(worksize+1);
    fill(workloc_.begin(), workloc_.end(), -1);
    size_t nitmp=0, nrtmp=0;
    size_t wind=0;
    for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it) {
      if (it->op!=OP_OUTPUT) {
        for (int c=0; c<it->res.size(); ++c) {
          if (it->res[c]>=0) {
            size_t nr=0, ni=0;
            it->data->nTmp(ni, nr);
            nitmp = std::max(nitmp, ni);
            nrtmp = std::max(nrtmp, nr);
            if (workloc_[it->res[c]] < 0) {
              workloc_[it->res[c]] = wind;
              wind += it->data->sparsity(c).nnz();
            }
          }
        }
      }
    }
    workloc_.back()=wind;
    for (int i=0; i<workloc_.size(); ++i) {
      if (workloc_[i]<0) workloc_[i] = i==0 ? 0 : workloc_[i-1];
      workloc_[i] += nrtmp;
    }
    itmp_.resize(nitmp);
    rtmp_.resize(nrtmp+wind);

    // Reset the temporary variables
    for (int i=0; i<nodes.size(); ++i) {
      if (nodes[i]) {
        nodes[i]->temp = 0;
      }
    }

    // Now mark each input's place in the algorithm
    for (vector<pair<int, MXNode*> >::const_iterator it=symb_loc.begin();
         it!=symb_loc.end(); ++it) {
      it->second->temp = it->first+1;
    }

    // Add input instructions, loop over inputs
    for (int ind=0; ind<inputv_.size(); ++ind) {
      // Loop over symbolic primitives of each input
      vector<MX> prim = inputv_[ind].getPrimitives();
      int nz_offset=0;
      for (int p=0; p<prim.size(); ++p) {
        int i = prim[p].getTemp()-1;
        if (i>=0) {
          // Mark as input
          algorithm_[i].op = OP_INPUT;

          // Location of the input
          algorithm_[i].arg.resize(3);
          algorithm_[i].arg[0] = ind;
          algorithm_[i].arg[1] = p;
          algorithm_[i].arg[2] = nz_offset;

          // Mark input as read
          prim[p].setTemp(0);
        }
        nz_offset += prim[p]->nnz();
      }
    }

    // Locate free variables
    free_vars_.clear();
    for (vector<pair<int, MXNode*> >::const_iterator it=symb_loc.begin();
         it!=symb_loc.end(); ++it) {
      int i = it->second->temp-1;
      if (i>=0) {
        // Save to list of free parameters
        free_vars_.push_back(MX::create(it->second));

        // Remove marker
        it->second->temp=0;
      }
    }

    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      profileWriteName(CasadiOptions::profilingLog, this, getOption("name"),
                       ProfilingData_FunctionType_MXFunction, algorithm_.size());
      int alg_counter = 0;
      for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end();
           ++it, ++alg_counter) {
        std::stringstream ss;
        print(ss, *it);
        if (it->op == OP_CALL) {
          profileWriteSourceLineDep(CasadiOptions::profilingLog, this, alg_counter,
                            ss.str(), it->op, it->data->getFunction().operator->());
        } else {
          profileWriteSourceLine(CasadiOptions::profilingLog, this, alg_counter, ss.str(), it->op);
        }
      }
    }

    log("MXFunctionInternal::init end");
  }

  void MXFunctionInternal::evalD(cp_double* arg,
                                 p_double* res, int* itmp, double* rtmp) {
    casadi_log("MXFunctionInternal::evalD():begin "  << getOption("name"));
    // Set up timers for profiling
    double time_zero=0;
    double time_start=0;
    double time_stop=0;
    if (CasadiOptions::profiling) {
      time_zero = getRealTime();
      if (CasadiOptions::profilingBinary) {
        profileWriteEntry(CasadiOptions::profilingLog, this);
      } else {
        CasadiOptions::profilingLog  << "start " << this << ":" <<getOption("name") << std::endl;
      }
    }

    // Work vector and temporaries to hold pointers to operation input and outputs
    vector<cp_double> oarg(max_arg_);
    vector<p_double> ores(max_res_);

    // Make sure that there are no free variables
    if (!free_vars_.empty()) {
      std::stringstream ss;
      repr(ss);
      casadi_error("Cannot evaluate \"" << ss.str() << "\" since variables "
                   << free_vars_ << " are free.");
    }

    // Evaluate all of the nodes of the algorithm:
    // should only evaluate nodes that have not yet been calculated!
    int alg_counter = 0;
    for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it, ++alg_counter) {
      if (CasadiOptions::profiling) {
        time_start = getRealTime(); // Start timer
      }

      if (it->op==OP_INPUT) {
        // Pass an input
        double *w = rtmp+workloc_[it->res.front()];
        int nnz=it->data.nnz();
        int i=it->arg.at(0);
        int nz_offset=it->arg.at(2);
        if (arg[i]==0) {
          fill(w, w+nnz, 0);
        } else {
          copy(arg[i]+nz_offset, arg[i]+nz_offset+nnz, w);
        }
      } else if (it->op==OP_OUTPUT) {
        // Get an output
        double *w = rtmp+workloc_[it->arg.front()];
        int i=it->res.front();
        if (res[i]!=0) copy(w, w+output(i).nnz(), res[i]);
      } else {
        // Point pointers to the data corresponding to the element
        for (int i=0; i<it->arg.size(); ++i)
          oarg[i] = it->arg[i]>=0 ? rtmp+workloc_[it->arg[i]] : 0;
        for (int i=0; i<it->res.size(); ++i)
          ores[i] = it->res[i]>=0 ? rtmp+workloc_[it->res[i]] : 0;

        // Evaluate
        it->data->evalD(getPtr(oarg), getPtr(ores), itmp, rtmp);
      }

      // Write out profiling information
      if (CasadiOptions::profiling) {
        time_stop = getRealTime(); // Stop timer

        if (CasadiOptions::profilingBinary) {
          profileWriteTime(CasadiOptions::profilingLog, this, alg_counter,
                           time_stop-time_start, time_stop-time_zero);
        } else {
          CasadiOptions::profilingLog  << (time_stop-time_start)*1e6 << " ns | "
                                       << (time_stop-time_zero)*1e3 << " ms | "
                                       << this << ":" <<getOption("name") << ":"
                                       << alg_counter <<"|";
          if (it->op == OP_CALL) {
            Function f = it->data->getFunction();
            CasadiOptions::profilingLog << f.get() << ":" << f.getOption("name");
          }
          CasadiOptions::profilingLog << "|";
          print(CasadiOptions::profilingLog, *it);
        }

      }
    }

    if (CasadiOptions::profiling) {
      time_stop = getRealTime();
      if (CasadiOptions::profilingBinary) {
        profileWriteExit(CasadiOptions::profilingLog, this, time_stop-time_zero);
      } else {
        CasadiOptions::profilingLog  << "stop " << this << ":"
                                     <<getOption("name") << (time_stop-time_zero)*1e3
                                     << " ms" << std::endl;
      }
    }

    casadi_log("MXFunctionInternal::evalD():end "  << getOption("name"));
  }

  void MXFunctionInternal::print(ostream &stream, const AlgEl& el) const {
    if (el.op==OP_OUTPUT) {
      stream << "output[" << el.res.front() << "] = @" << el.arg.at(0);
    } else if (el.op==OP_SETNONZEROS || el.op==OP_ADDNONZEROS) {
      if (el.res.front()!=el.arg.at(0)) {
        stream << "@" << el.res.front() << " = @" << el.arg.at(0) << "; ";
      }
      stream << "@" << el.res.front();
      el.data->printPart(stream, 1);
      stream << "@" << el.arg.at(1);
    } else {
      if (el.res.size()==1) {
        stream << "@" << el.res.front() << " = ";
      } else {
        stream << "{";
        for (int i=0; i<el.res.size(); ++i) {
          if (i!=0) stream << ", ";
          if (el.res[i]>=0) {
            stream << "@" << el.res[i];
          } else {
            stream << "NULL";
          }
        }
        stream << "} = ";
      }
      if (el.op==OP_INPUT) {
        stream << "input[" << el.arg.at(0) << "][" << el.arg.at(1) << "]";
      } else {
        el.data->printPart(stream, 0);
        for (int i=0; i<el.arg.size(); ++i) {
          if (el.arg[i]>=0) {
            stream << "@" << el.arg[i];
          } else {
            stream << "NULL";
          }
          el.data->printPart(stream, i+1);
        }
      }
    }
    stream << endl;
  }

  void MXFunctionInternal::print(ostream &stream) const {
    FunctionInternal::print(stream);
    for (vector<AlgEl>::const_iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it) {
      print(stream, *it);
    }
  }

  MXFunctionInternal* MXFunctionInternal::clone() const {
    return new MXFunctionInternal(*this);
  }

  void MXFunctionInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    XFunctionInternal<MXFunction, MXFunctionInternal, MX, MXNode>::deepCopyMembers(already_copied);
    for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it) {
      switch (it->op) {
      case OP_CALL:
      case OP_SOLVE:
        it->data.makeUnique(already_copied, false);
        it->data->getFunction() = deepcopy(it->data->getFunction(), already_copied);
        break;
      default:
        break;
      }
    }
  }

  void MXFunctionInternal::spInit(bool fwd) {
    bvec_t *iwork = get_bvec_t(rtmp_);
    fill(iwork+workloc_.front(), iwork+workloc_.back(), bvec_t(0));
  }

  void MXFunctionInternal::spFwd(cp_bvec_t* arg, p_bvec_t* res,
                                 int* itmp, bvec_t* rtmp) {
    // Tmporaries to hold pointers to operation input and outputs
    vector<cp_bvec_t> oarg(max_arg_);
    vector<p_bvec_t> ores(max_res_);

    // Propagate sparsity forward
    for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); it++) {
      if (it->op==OP_INPUT) {
        // Pass input seeds
        int nnz=it->data.nnz();
        int i=it->arg.at(0);
        int nz_offset=it->arg.at(2);
        const bvec_t* iw = arg[i];
        bvec_t* w = rtmp + workloc_[it->res.front()];
        if (iw!=0) {
          copy(iw+nz_offset, iw+nz_offset+nnz, w);
        } else {
          fill_n(w, nnz, 0);
        }
      } else if (it->op==OP_OUTPUT) {
        // Get the output sensitivities
        int i=it->res.front();
        int nnz=output(i).nnz();
        bvec_t* ow = res[i];
        bvec_t* w = rtmp + workloc_[it->arg.front()];
        if (ow!=0) copy(w, w+nnz, ow);
      } else {
        // Point pointers to the data corresponding to the element
        for (int i=0; i<it->arg.size(); ++i)
          oarg[i] = it->arg[i]>=0 ? rtmp+workloc_[it->arg[i]] : 0;
        for (int i=0; i<it->res.size(); ++i)
          ores[i] = it->res[i]>=0 ? rtmp+workloc_[it->res[i]] : 0;

        // Propagate sparsity forwards
        it->data->spFwd(getPtr(oarg), getPtr(ores), itmp, rtmp);
      }
    }
  }

  void MXFunctionInternal::spAdj(p_bvec_t* arg, p_bvec_t* res,
                                 int* itmp, bvec_t* rtmp) {
    // Tmporaries to hold pointers to operation input and outputs
    vector<bvec_t*> ores(max_res_);
    vector<bvec_t*> oarg(max_arg_); // Non-const since seeds are cleared

    size_t ni, nr;
    nTmp(ni, nr);
    fill_n(rtmp, nr, 0);

    // Propagate sparsity backwards
    for (vector<AlgEl>::reverse_iterator it=algorithm_.rbegin(); it!=algorithm_.rend(); it++) {
      if (it->op==OP_INPUT) {
        // Get the input sensitivities and clear it from the work vector
        int nnz=it->data.nnz();
        int i=it->arg.at(0);
        int nz_offset=it->arg.at(2);
        bvec_t* iw = arg[i];
        bvec_t* w = rtmp + workloc_[it->res.front()];
        if (iw!=0) for (int k=0; k<nnz; ++k) iw[nz_offset+k] |= w[k];
        fill_n(w, nnz, 0);
      } else if (it->op==OP_OUTPUT) {
        // Pass output seeds
        int i=it->res.front();
        int nnz=output(i).nnz();
        bvec_t* ow = res[i];
        bvec_t* w = rtmp + workloc_[it->arg.front()];
        if (ow!=0) {
          for (int k=0; k<nnz; ++k) w[k] |= ow[k];
          fill_n(ow, nnz, 0);
        }
      } else {
        // Point pointers to the data corresponding to the element
        for (int i=0; i<it->arg.size(); ++i)
          oarg[i] = it->arg[i]>=0 ? rtmp+workloc_[it->arg[i]] : 0;
        for (int i=0; i<it->res.size(); ++i)
          ores[i] = it->res[i]>=0 ? rtmp+workloc_[it->res[i]] : 0;

        // Propagate sparsity backwards
        it->data->spAdj(getPtr(oarg), getPtr(ores), itmp, rtmp);
      }
    }
  }

  Function MXFunctionInternal::getNumericJacobian(int iind, int oind,
                                                  bool compact, bool symmetric) {
    // Create expressions for the Jacobian
    vector<MX> ret_out;
    ret_out.reserve(1+outputv_.size());
    ret_out.push_back(jac(iind, oind, compact, symmetric, false, true));
    ret_out.insert(ret_out.end(), outputv_.begin(), outputv_.end());

    MXFunction ret(inputv_, ret_out);
    ret.setInputScheme(inputScheme());
    // Return function
    return ret;
  }

  std::vector<MX> MXFunctionInternal::symbolicOutput(const std::vector<MX>& arg) {
    // Check if input is given
    const int checking_depth = 2;
    bool input_given = true;
    for (int i=0; i<arg.size() && input_given; ++i) {
      if (!isEqual(arg[i], inputv_[i], checking_depth)) {
        input_given = false;
      }
    }

    // Return output if possible, else fall back to base class
    if (input_given) {
      return outputv_;
    } else {
      return FunctionInternal::symbolicOutput(arg);
    }
  }

  void MXFunctionInternal::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    log("MXFunctionInternal::evalMX begin");
    assertInit();
    casadi_assert_message(arg.size()==getNumInputs(), "Wrong number of input arguments");

    // Resize the number of outputs
    res.resize(outputv_.size());

    // Check if arguments matches the input expressions, in which case
    // the output is known to be the output expressions
    const int checking_depth = 2;
    bool output_given = true;
    for (int i=0; i<arg.size() && output_given; ++i) {
      if (!isEqual(arg[i], inputv_[i], checking_depth)) {
        output_given = false;
      }
    }

    // Copy output if known
    if (output_given) {
      copy(outputv_.begin(), outputv_.end(), res.begin());
      return;
    }

    // Symbolic work, non-differentiated
    vector<MX> swork(workloc_.size()-1);
    log("MXFunctionInternal::evalMX allocated work vector");

    // Split up inputs analogous to symbolic primitives
    vector<vector<MX> > arg_split(arg.size());
    for (int i=0; i<arg.size(); ++i)
      arg_split[i] = inputv_[i].splitPrimitives(arg[i]);

    vector<MX> oarg, ores;
    oarg.reserve(max_arg_);
    ores.reserve(max_res_);

    // Loop over computational nodes in forward order
    int alg_counter = 0;
    for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it, ++alg_counter) {
      if (it->op == OP_INPUT) {
        swork[it->res.front()] = arg_split.at(it->arg.at(0)).at(it->arg.at(1))
          .setSparse(it->data.sparsity(), true);
      } else if (it->op==OP_OUTPUT) {
        // Collect the results
        res[it->res.front()] = swork[it->arg.front()];
      } else if (it->op==OP_PARAMETER) {
        // Fetch parameter
        swork[it->res.front()] = it->data;
      } else {
        // Arguments of the operation
        oarg.resize(it->arg.size());
        for (int i=0; i<oarg.size(); ++i) {
          int el = it->arg[i]; // index of the argument
          oarg[i] = el<0 ? MX(it->data->dep(i).shape()) : swork[el];
        }

        // Perform the operation
        ores.resize(it->res.size());
        it->data->evalMX(oarg, ores);

        // Get the result
        for (int i=0; i<ores.size(); ++i) {
          int el = it->res[i]; // index of the output
          if (el>=0) swork[el] = ores[i];
        }
      }
    }
    log("MXFunctionInternal::evalMX end");
  }

  void MXFunctionInternal::evalFwd(const std::vector<std::vector<MX> >& fseed,
                                   std::vector<std::vector<MX> >& fsens) {
    log("MXFunctionInternal::evalFwd begin");
    assertInit();

    // Allocate results
    int nfwd = fseed.size();
    fsens.resize(nfwd);
    for (int d=0; d<nfwd; ++d) {
      fsens[d].resize(getNumOutputs());
    }

    // Quick return if no directions
    if (nfwd==0) return;

    // Check if there are any zero seeds
    for (vector<vector<MX> >::const_iterator i=fseed.begin(); i!=fseed.end(); ++i) {
      // If any direction can be skipped
      if (purgable(*i)) {
        // New argument without all-zero directions
        std::vector<std::vector<MX> > fseed_purged, fsens_purged;
        fseed_purged.reserve(nfwd);
        vector<int> index_purged;
        for (int d=0; d<nfwd; ++d) {
          if (purgable(fseed[d])) {
            for (int i=0; i<fsens[d].size(); ++i) {
              fsens[d][i] = MX(output(i).shape());
            }
          } else {
            fseed_purged.push_back(fsens[d]);
            index_purged.push_back(d);
          }
        }

        // Call recursively
        evalFwd(fseed_purged, fsens_purged);

        // Fetch result
        for (int d=0; d<fseed_purged.size(); ++d) {
          fsens[index_purged[d]] = fsens_purged[d];
        }
        return;
      }
    }

    // Allocate forward sensitivities
    for (int d=0; d<nfwd; ++d) {
      fsens[d].resize(outputv_.size());
    }

    // Work vector, forward derivatives
    std::vector<std::vector<MX> > dwork(workloc_.size()-1);
    fill(dwork.begin(), dwork.end(), std::vector<MX>(nfwd));
    log("MXFunctionInternal::evalFwd allocated derivative work vector (forward mode)");

    // Split up fseed analogous to symbolic primitives
    vector<vector<vector<MX> > > fseed_split(nfwd);
    for (int d=0; d<nfwd; ++d) {
      fseed_split[d].resize(fseed[d].size());
      for (int i=0; i<fseed[d].size(); ++i) {
        fseed_split[d][i] = inputv_[i].splitPrimitives(fseed[d][i]);
      }
    }

    // Pointers to the arguments of the current operation
    vector<vector<MX> > oseed, osens;
    oseed.reserve(nfwd);
    osens.reserve(nfwd);
    vector<bool> skip(nfwd, false);

    // Loop over computational nodes in forward order
    for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it) {
      if (it->op == OP_INPUT) {
        // Fetch forward seed
        for (int d=0; d<nfwd; ++d) {
          dwork[it->res.front()][d] = fseed_split[d].at(it->arg.at(0)).at(it->arg.at(1))
            .setSparse(it->data.sparsity(), true);
        }
      } else if (it->op==OP_OUTPUT) {
        // Collect forward sensitivity
        for (int d=0; d<nfwd; ++d) {
          fsens[d][it->res.front()] = dwork[it->arg.front()][d];
        }
      } else if (it->op==OP_PARAMETER) {
        // Fetch parameter
        for (int d=0; d<nfwd; ++d) {
          dwork[it->res.front()][d] = MX();
        }
      } else {
        // Get seeds, ignoring all-zero directions
        oseed.clear();
        for (int d=0; d<nfwd; ++d) {
          // Collect seeds, skipping directions with only zeros
          vector<MX> seed(it->arg.size());
          skip[d] = true; // All seeds are zero?
          for (int i=0; i<it->arg.size(); ++i) {
            int el = it->arg[i];
            if (el<0 || dwork[el][d].isEmpty(true)) {
              seed[i] = MX(it->data->dep(i).shape());
            } else {
              seed[i] = dwork[el][d];
            }
            if (skip[d] && !seed[i].isZero()) skip[d] = false;
          }
          if (!skip[d]) oseed.push_back(seed);
        }

        // Perform the operation
        osens.resize(oseed.size());
        if (!osens.empty()) {
          fill(osens.begin(), osens.end(), vector<MX>(it->res.size()));
          it->data->evalFwd(oseed, osens);
        }

        // Store sensitivities
        int d1=0;
        for (int d=0; d<nfwd; ++d) {
          for (int i=0; i<it->res.size(); ++i) {
            int el = it->res[i];
            if (el>=0) {
              dwork[el][d] = skip[d] ? MX(it->data->sparsity(i).shape()) : osens[d1][i];
            }
          }
          if (!skip[d]) d1++;
        }
      }
    }
    log("MXFunctionInternal::evalFwd end");
  }

  void MXFunctionInternal::evalAdj(const std::vector<std::vector<MX> >& aseed,
                                   std::vector<std::vector<MX> >& asens) {
    log("MXFunctionInternal::evalAdj begin");
    assertInit();

    // Allocate results
    int nadj = aseed.size();
    asens.resize(nadj);
    for (int d=0; d<nadj; ++d) {
      asens[d].resize(getNumInputs());
    }

    // Quick return if no directions
    if (nadj==0) return;

    // Check if there are any zero seeds
    for (vector<vector<MX> >::const_iterator i=aseed.begin(); i!=aseed.end(); ++i) {
      // If any direction can be skipped
      if (purgable(*i)) {
        // New argument without all-zero directions
        std::vector<std::vector<MX> > aseed_purged, asens_purged;
        aseed_purged.reserve(nadj);
        vector<int> index_purged;
        for (int d=0; d<nadj; ++d) {
          if (purgable(aseed[d])) {
            for (int i=0; i<asens[d].size(); ++i) {
              asens[d][i] = MX(input(i).shape());
            }
          } else {
            aseed_purged.push_back(asens[d]);
            index_purged.push_back(d);
          }
        }

        // Call recursively
        evalAdj(aseed_purged, asens_purged);

        // Fetch result
        for (int d=0; d<aseed_purged.size(); ++d) {
          asens[index_purged[d]] = asens_purged[d];
        }
        return;
      }
    }

    // Allocate splited adjoint sensitivities
    vector<vector<vector<MX> > > asens_split(nadj);
    for (int d=0; d<nadj; ++d) {
      asens_split[d].resize(inputv_.size());
      for (int i=0; i<inputv_.size(); ++i) {
        asens_split[d][i].resize(inputv_[i].numPrimitives());
      }
    }

    // Pointers to the arguments of the current operation
    vector<vector<MX> > oseed, osens;
    oseed.reserve(nadj);
    osens.reserve(nadj);
    vector<bool> skip(nadj, false);

    // Work vector, adjoint derivatives
    std::vector<std::vector<MX> > dwork(workloc_.size()-1);
    fill(dwork.begin(), dwork.end(), std::vector<MX>(nadj));
    log("MXFunctionInternal::evalAdj allocated derivative work vector (adjoint mode)");

    // Loop over computational nodes in reverse order
    for (vector<AlgEl>::reverse_iterator it=algorithm_.rbegin(); it!=algorithm_.rend(); ++it) {
      if (it->op == OP_INPUT) {
        // Get the adjoint sensitivities
        for (int d=0; d<nadj; ++d) {
          asens_split[d].at(it->arg.at(0)).at(it->arg.at(1)) = dwork[it->res.front()][d];
          dwork[it->res.front()][d] = MX();
        }
      } else if (it->op==OP_OUTPUT) {
        // Pass the adjoint seeds
        for (int d=0; d<nadj; ++d) {
          MX a = aseed[d][it->res.front()].setSparse(output(it->res.front()).sparsity(), true);
          if (dwork[it->arg.front()][d].isEmpty(true)) {
            dwork[it->arg.front()][d] = a;
          } else {
            dwork[it->arg.front()][d] += a;
          }
        }
      } else if (it->op==OP_PARAMETER) {
        // Clear adjoint seeds
        for (int d=0; d<nadj; ++d) {
          dwork[it->res.front()][d] = MX();
        }
      } else {
        // Collect and reset seeds
        oseed.clear();
        for (int d=0; d<nadj; ++d) {
          // Can the direction be skipped completely?
          skip[d] = true;

          // Seeds for direction d
          vector<MX> seed(it->res.size());
          for (int i=0; i<it->res.size(); ++i) {
            // Get and clear seed
            int el = it->res[i];
            if (el>=0) {
              seed[i] = dwork[el][d];
              dwork[el][d] = MX();
            } else {
              seed[i] = MX();
            }

            // If first time encountered, reset to zero of right dimension
            if (seed[i].isEmpty(true)) seed[i] = MX(it->data->sparsity(i).shape());

            // If nonzero seeds, keep direction
            if (skip[d] && !seed[i].isZero()) skip[d] = false;
          }
          // Add to list of derivatives
          if (!skip[d]) oseed.push_back(seed);
        }

        // Get values of sensitivities before addition
        osens.resize(oseed.size());
        int d1=0;
        for (int d=0; d<nadj; ++d) {
          if (skip[d]) continue;
          osens[d1].resize(it->arg.size());
          for (int i=0; i<it->arg.size(); ++i) {
            // Pass seed and reset to avoid counting twice
            int el = it->arg[i];
            if (el>=0) {
              osens[d1][i] = dwork[el][d];
              dwork[el][d] = MX();
            } else {
              osens[d1][i] = MX();
            }

            // If first time encountered, reset to zero of right dimension
            if (osens[d1][i].isEmpty(true)) osens[d1][i] = MX(it->data->dep(i).shape());
          }
          d1++;
        }

        // Perform the operation
        if (!osens.empty()) {
          it->data->evalAdj(oseed, osens);
        }

        // Store sensitivities
        d1=0;
        for (int d=0; d<nadj; ++d) {
          if (skip[d]) continue;
          for (int i=0; i<it->arg.size(); ++i) {
            int el = it->arg[i];
            if (el>=0) {
              if (dwork[el][d].isEmpty(true)) {
                dwork[el][d] = osens[d1][i];
              } else {
                dwork[el][d] += osens[d1][i];
              }
            }
          }
          d1++;
        }
      }
    }

    // Get adjoint sensitivities
    for (int d=0; d<nadj; ++d) {
      asens[d].resize(inputv_.size());
      for (int i=0; i<inputv_.size(); ++i) {
        asens[d][i] = inputv_[i].joinPrimitives(asens_split[d][i]);
      }
    }

    log("MXFunctionInternal::evalAdj end");
  }

  void MXFunctionInternal::evalSX(cp_SXElement* arg, p_SXElement* res,
                                  int* itmp, SXElement* rtmp) {
    // Work vector and temporaries to hold pointers to operation input and outputs
    vector<cp_SXElement> argp(max_arg_);
    vector<p_SXElement> resp(max_res_);

    // Evaluate all of the nodes of the algorithm:
    // should only evaluate nodes that have not yet been calculated!
    for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); it++) {
      if (it->op==OP_INPUT) {
        // Pass an input
        SXElement *w = rtmp+workloc_[it->res.front()];
        int nnz=it->data.nnz();
        int i=it->arg.at(0);
        int nz_offset=it->arg.at(2);
        if (arg[i]==0) {
          std::fill(w, w+nnz, 0);
        } else {
          std::copy(arg[i]+nz_offset, arg[i]+nz_offset+nnz, w);
        }
      } else if (it->op==OP_OUTPUT) {
        // Get the outputs
        SXElement *w = rtmp+workloc_[it->arg.front()];
        int i=it->res.front();
        if (res[i]!=0)
          std::copy(w, w+output(i).nnz(), res[i]);
      } else if (it->op==OP_PARAMETER) {
        continue; // FIXME
      } else {
        // Point pointers to the data corresponding to the element
        for (int i=0; i<it->arg.size(); ++i)
          argp[i] = it->arg[i]>=0 ? rtmp+workloc_[it->arg[i]] : 0;
        for (int i=0; i<it->res.size(); ++i)
          resp[i] = it->res[i]>=0 ? rtmp+workloc_[it->res[i]] : 0;

        // Evaluate
        it->data->evalSX(getPtr(argp), getPtr(resp), itmp, rtmp);
      }
    }
  }

  SXFunction MXFunctionInternal::expand(const std::vector<SX>& inputvsx) {
    assertInit();

    // Create inputs with the same name and sparsity as the matrix valued symbolic inputs
    vector<SX> arg(inputv_.size());
    if (inputvsx.empty()) { // No symbolic input provided
      for (int i=0; i<arg.size(); ++i) {
        // Start with matrix with the correct sparsity
        arg[i] = SX(inputv_[i].sparsity());

        // Divide input into primitives and create corresponding SX
        vector<SXElement>::iterator ait = arg[i].begin();
        vector<MX> prim = inputv_[i].getPrimitives();
        for (vector<MX>::const_iterator pit=prim.begin(); pit!=prim.end(); ++pit) {
          SX t = SX::sym(pit->getName(), pit->sparsity());
          copy(t.begin(), t.end(), ait);
          ait += t.nnz();
        }
        casadi_assert(ait==arg[i].end());
      }
    } else { // Use provided symbolic input
      // Make sure number of inputs matches
      casadi_assert(inputvsx.size()==inputv_.size());

      // Make sure that sparsity matches
      for (int i=0; i<inputvsx.size(); ++i) {
        casadi_assert(inputvsx[i].sparsity() == inputv_[i].sparsity());
      }

      // Copy to argument vector
      copy(inputvsx.begin(), inputvsx.end(), arg.begin());
    }

    // Create output vector with correct sparsity
    vector<SX> res(outputv_.size());
    for (int i=0; i<res.size(); ++i) {
      res[i] = SX(outputv_[i].sparsity());
    }

    // Evaluate symbolically
    FunctionInternal::evalSX(arg, res);

    // Create function
    SXFunction f(arg, res);
    f.setInputScheme(getInputScheme());
    f.setOutputScheme(getOutputScheme());
    string name = getOption("name");
    f.setOption("name", "expand_" + name);
    return f;
  }

  void MXFunctionInternal::printWork(ostream &stream) {
    for (int k=0; k<workloc_.size()-1; ++k) {
      vector<double>::const_iterator start=rtmp_.begin() + workloc_[k];
      vector<double>::const_iterator stop=rtmp_.begin() + workloc_[k+1];
      stream << "work[" << k << "] = " << vector<double>(start, stop) << endl;
    }
  }

  void MXFunctionInternal::generateDeclarations(std::ostream &stream, const std::string& type,
                                                CodeGenerator& gen) const {

    // Make sure that there are no free variables
    if (!free_vars_.empty()) {
      casadi_error("Code generation is not possible since variables "
                   << free_vars_ << " are free.");
    }

    // Generate code for the embedded functions
    for (vector<AlgEl>::const_iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it) {
      switch (it->op) {
      case OP_CALL:
      case OP_SOLVE:
        gen.addDependency(it->data->getFunction());
        break;
      default:
        break;
      }
    }
  }

  void MXFunctionInternal::generateBody(std::ostream &stream, const std::string& type,
                                        CodeGenerator& gen) const {
    // Temporary variables and vectors
    stream << "  int i, j, k, *ii, *jj, *kk;" << endl;
    stream << "  d r, s, t, *rr, *ss, *tt;" << endl;
    stream << "  const d *cr, *cs, *ct;" << endl;

    // Operation number (for printing)
    int k=0;

    // Names of operation argument and results
    vector<int> arg, res;

    // Print class (for debugging)
    bool codegen_class = true;

    // Codegen the algorithm
    for (vector<AlgEl>::const_iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it) {
      // Mark the beginning of the operation
      stream << "  /* " << k++;
      if (codegen_class) {
        if (it->data.get()!=0) {
          stream << " : " << typeid(*it->data.get()).name();

          // if this is a call node, also write the name of the Function
          MX algElem = it->data;
          if (algElem.getOp() == OP_CALL) {
            stream << " (" << algElem.getFunction().getSanitizedName() << ")";
          }
        }
      }
      stream << " */" << endl;

      // Print the operation
      if (it->op==OP_OUTPUT) {
        gen.copyVector(stream, CodeGenerator::work(workloc_[it->arg.front()]),
                       output(it->res.front()).nnz(),
                       "res[" + CodeGenerator::numToString(it->res.front()) + "]", "i", true);
      } else if (it->op==OP_INPUT) {
        std::string arg = "arg[" + CodeGenerator::numToString(it->arg.at(0)) + "]";
        if (it->arg.at(2)!=0) arg += "+" + CodeGenerator::numToString(it->arg.at(2));
        gen.copyVector(stream, arg, it->data.nnz(),
                       CodeGenerator::work(workloc_[it->res.front()]), "i", false);
      } else {
        // Get the names of the operation arguments
        arg.resize(it->arg.size());
        for (int i=0; i<it->arg.size(); ++i) {
          if (it->arg.at(i)>=0) {
            arg.at(i) = workloc_.at(it->arg.at(i));
          } else {
            arg.at(i) = -1;
          }
        }

        // Get the names of the operation results
        res.resize(it->res.size());
        for (int i=0; i<it->res.size(); ++i) {
          if (it->res.at(i)>=0) {
            res.at(i) = workloc_.at(it->res.at(i));
          } else {
            res.at(i) = -1;
          }
        }

        // Generate operation
        it->data->generate(stream, arg, res, gen);
      }
    }
  }

  void MXFunctionInternal::generateLiftingFunctions(MXFunction& vdef_fcn, MXFunction& vinit_fcn) {
    assertInit();

    vector<MX> swork(workloc_.size()-1);

    vector<MX> oarg, ores;
    oarg.reserve(max_arg_);
    ores.reserve(max_res_);

    // Definition of intermediate variables
    vector<MX> y;
    vector<MX> g;
    vector<MX> f_G(getNumOutputs());

    // Initial guess for intermediate variables
    vector<MX> x_init;

    // Temporary stringstream
    stringstream ss;

    for (int algNo=0; algNo<2; ++algNo) {
      for (vector<AlgEl>::iterator it=algorithm_.begin(); it!=algorithm_.end(); ++it) {
        switch (it->op) {
        case OP_LIFT:
          {
            MX& arg = swork[it->arg.at(0)];
            MX& arg_init = swork[it->arg.at(1)];
            MX& res = swork[it->res.front()];
            switch (algNo) {
            case 0:
              ss.str(string());
              ss << "y" << y.size();
              y.push_back(MX::sym(ss.str(), arg.sparsity()));
              g.push_back(arg);
              res = y.back();
              break;
            case 1:
              x_init.push_back(arg_init);
              res = arg_init;
              break;
            }
            break;
          }
        case OP_INPUT:
        case OP_PARAMETER:
          swork[it->res.front()] = it->data;
          break;
        case OP_OUTPUT:
          if (algNo==0) {
            f_G[it->res.front()] = swork[it->arg.front()];
          }
          break;
        default:
          {
            // Arguments of the operation
            oarg.resize(it->arg.size());
            for (int i=0; i<oarg.size(); ++i) {
              int el = it->arg[i]; // index of the argument
              oarg[i] = el<0 ? MX(it->data->dep(i).shape()) : swork[el];
            }

            // Perform the operation
            ores.resize(it->res.size());
            it->data->evalMX(oarg, ores);

            // Get the result
            for (int i=0; i<ores.size(); ++i) {
              int el = it->res[i]; // index of the output
              if (el>=0) swork[el] = ores[i];
            }
          }
        }
      }
    }

    // Definition of intermediate variables
    vector<MX> f_in = inputv_;
    f_in.insert(f_in.end(), y.begin(), y.end());
    vector<MX> f_out = f_G;
    f_out.insert(f_out.end(), g.begin(), g.end());
    vdef_fcn = MXFunction(f_in, f_out);
    vdef_fcn.setOption("name", "lifting_variable_definition");

    // Initial guess of intermediate variables
    f_in = inputv_;
    f_out = x_init;
    vinit_fcn = MXFunction(f_in, f_out);
    vinit_fcn.setOption("name", "lifting_variable_guess");
  }

} // namespace casadi

