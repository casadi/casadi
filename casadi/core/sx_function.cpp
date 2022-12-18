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
#include <sstream>
#include <iomanip>
#include "sx_node.hpp"
#include "casadi_common.hpp"
#include "sparsity_internal.hpp"
#include "casadi_interrupt.hpp"
#include "serializing_stream.hpp"
#include "fun_ref.hpp"
#include "call_sx.hpp"
#include <random>

namespace casadi {

  using namespace std;


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
      case OP_INPUT: {
        w[e.i0] = arg[e.i1]==nullptr ? 0 : arg[e.i1][e.i2];
        break;
      }
      case OP_OUTPUT: {
        if (res[e.i0]!=nullptr) res[e.i0][e.i2] = w[e.i1];
        break;
      }
      case OP_CALL: {
        int id;
        double val;
        FunRef::unpack(w[e.i1], id, val);
        if (CallSX::call(functions_[id], w[e.i0], w[e.i2], val, arg+n_in_, res+n_out_, iw, w+worksize_)) return 1;
        break;
      }
      case OP_FUNREF: w[e.i0] = FunRef::pack(e.i2, w[e.i1]); break;
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

  void SXFunction::codegen_declarations(CodeGenerator& g, const Instance& inst) const {

    // Make sure that there are no free variables
    if (!free_vars_.empty()) {
      casadi_error("Code generation of '" + name_ + "' is not possible since variables "
                   + str(free_vars_) + " are free.");
    }

    // Run the algorithm
    for (auto&& f : functions_) {
      Instance local;
      local.stride_in.resize(f.n_in(), 1);
      local.stride_out.resize(f.n_out(), 1);
      local.prefer_inline = true;
      CallSX::codegen_dependency(g, f, local, shared_from_this<Function>());
      
    }
  }

  void SXFunction::codegen_body(CodeGenerator& g, const Instance& inst) const {
    std::map<const AlgEl*, casadi_int> lookup_rets, lookup_ins;
    std::vector< const AlgEl* > rets, ins;

    bool vectorize = false;
    for (auto e : inst.stride_in) {
      if (e>1) vectorize = true;
    }
    for (auto e : inst.stride_out) {
      if (e>1) vectorize = true;
    }

    const SXNode*const* node_ptr = get_ptr(nodes_);

    // Run the algorithm
    for (auto&& a : algorithm_) {
      if (a.op==OP_OUTPUT) {
        if (!inst.res_null.empty()) {
          if (!inst.res_null[a.i0]) {
            rets.push_back(&a);
          }
        }
      } else if (a.op==OP_INPUT) {
        if (!inst.arg_null.empty()) {
          if (!inst.arg_null[a.i1]) {
            ins.push_back(&a);
          }
        }
      }
    }

    struct {
      bool operator()(const AlgEl* a, const AlgEl* b) const { 
        if (a->i0 == b->i0) {
          return a->i2 < b->i2;
        } else {
          return a->i0 < b->i0;
        }
      };
    } sortme_rets;

    struct {
      bool operator()(const AlgEl* a, const AlgEl* b) const { 
        if (a->i1 == b->i1) {
          return a->i2 < b->i2;
        } else {
          return a->i1 < b->i1;
        }
      };
    } sortme_ins;

    bool intro = false;//vectorize;
    bool outro = false;//vectorize;
    intro = vectorize;
    outro = vectorize;

    bool intro2_step = vectorize;

    //intro = false;
    //outro = false;
    //intro2_step = false;

    std::map<int, Function> store;

    if (outro) {
      std::sort(rets.begin(), rets.end(), sortme_rets);
      for (casadi_int i=0;i<rets.size();++i) {
        const AlgEl& a = *rets[i];
        g.local("ret"+str(i), "casadi_real");
        lookup_rets[&a] = i;
      }
    }
    if (intro2_step) {
      for (casadi_int i=0;i<n_in_;++i) {
        g.local("inp"+str(i), "const casadi_real", "*");
        g << "inp" << i << " = " << g.arg(i, true) << ";\n";
      }
    }
    if (intro) {
      std::sort(ins.begin(), ins.end(), sortme_ins);
      for (casadi_int i=0;i<ins.size();++i) {
        const AlgEl& a = *ins[i];
        g.local("in"+str(i), "casadi_real");
        int stride = inst.stride_in.empty() ? 1 : inst.stride_in.at(a.i1);
        stride = 1;
        bool external_i = inst.stride_in.empty() ? true : inst.stride_in.at(a.i1)>0;
        g << "in" << i << " = " << ( intro2_step ? "inp" + str(a.i1): g.arg(a.i1, true)) << "[" << a.i2*stride << (external_i ? "+i": "") << "]" << ";\n";
        lookup_ins[&a] = i;
      }
    }

    // Run the algorithm
    for (auto&& a : algorithm_) {
      if (a.op==OP_OUTPUT) {
        if (inst.res_null.empty()) {
          g << "if (res[" << a.i0 << "]!=0) ";
          int stride = inst.stride_out.empty() ? 1 : inst.stride_out.at(a.i0);
          stride = 1;
          g << g.res(a.i0, true) << "[" << a.i2*abs(stride) << "]" << (stride<0? "+": "")<<  "=" << g.sx_work(a.i1);
        } else {
          if (!inst.res_null[a.i0]) {
            if (outro) {
              g << "ret" << lookup_rets[&a] << " =" << g.sx_work(a.i1);
            } else {
              int stride = inst.stride_out.empty() ? 1 : inst.stride_out.at(a.i0);
              stride = 1;
              g << g.res(a.i0, true) << "[" << a.i2*abs(stride) << "]" << " =" << g.sx_work(a.i1);
            }
          }
        }
      } else if (a.op==OP_CALL) {
        Instance local;
        local.prefer_inline = true;
        g << CallSX::codegen(g, (*node_ptr)->dep(0), local, a.i0, a.i1, a.i2, "arg+" + str(n_in_), "res+" + str(n_out_), "iw", "w+" + str(worksize_), shared_from_this<Function>()); 
      } else {

        // Where to store the result
        g << g.sx_work(a.i0) << "=";

        // What to store
        if (a.op==OP_CONST) {
          g << g.constant(a.d);
        } else if (a.op==OP_INPUT) {
          if (inst.arg_null.empty()) {
            g << g.arg(a.i1, true) << "? " << g.arg(a.i1, true) << "[" << a.i2 << "] : 0";
          } else {
            if (inst.arg_null[a.i1]) {
              g << "0";
            } else {
              if (intro) {
                g << "in" << lookup_ins[&a]; //g.arg(a.i1) << "[" << a.i2 << "]";
              } else {
                int stride = inst.stride_in.empty() ? 1 : inst.stride_in.at(a.i1);
                stride=1;
                g << g.arg(a.i1, true) << "[" << a.i2*stride << "]";
              }
            }
          }
        } else if (a.op==OP_FUNREF) {
          g << g.sx_work(a.i1);
        } else {
          casadi_int ndep = casadi_math<double>::ndeps(a.op);
          casadi_assert_dev(ndep>0);
          if (ndep==1) g << g.print_op(a.op, g.sx_work(a.i1)); // a.op OP_NEG
          if (ndep==2) g << g.print_op(a.op, g.sx_work(a.i1), g.sx_work(a.i2)); // OP_MUL OP_ADD
        }
      }
      g  << ";\n";
      if (a.op!=OP_OUTPUT) node_ptr++;
    }
    if (outro) {
      for (casadi_int i=0;i<rets.size();++i) {
        const AlgEl& a = *rets[i];
        casadi_int stride = inst.stride_out.empty() ? 1 : inst.stride_out[a.i0];
        stride = 1;
        bool external_i = inst.stride_out.empty() ? true : inst.stride_out.at(a.i0)>0;
        g << g.res(a.i0, true) << "[" << a.i2*abs(stride) << (external_i ? "+i": "") << "]" << (stride<0? "+": "")<< "= ret" << i << ";\n";
      }
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
      {"layout_in",
       {OT_LAYOUTVECTOR,
        "Layout in"}},
      {"layout_out",
       {OT_LAYOUTVECTOR,
        "Layout out"}},
      {"stride_in",
       {OT_INTVECTOR,
        "Layout in"}},
      {"stride_out",
       {OT_INTVECTOR,
        "Layout out"}}
     }
  };

  Dict SXFunction::generate_options(bool is_temp, bool keep_dim) const {
    Dict opts = FunctionInternal::generate_options(is_temp, keep_dim);
    //opts["default_in"] = default_in_;
    opts["live_variables"] = live_variables_;
    opts["just_in_time_sparsity"] = just_in_time_sparsity_;
    opts["just_in_time_opencl"] = just_in_time_opencl_;
    return opts;
  }

  typedef int Order;
  typedef int NodeIndex;

  struct Entry {
      NodeIndex next;
      Order p;
      NodeIndex prev;
  };


  class SimulatedAnnealingInstrOrder {
    public:
      int N, nc;
      std::vector<NodeIndex> i_left, i_right;
      std::vector< std::vector<NodeIndex> > incoming, outgoing;
      bool inf;
      long long int stationary_length;

      SimulatedAnnealingInstrOrder(std::vector<SXNode*>& nodes, long long int stationary_length, bool with_inf) {
        this->inf = inf;
        this->stationary_length = stationary_length;
        // Set the temporary variables to be the corresponding place in the sorted graph
        for (casadi_int i=0; i<nodes.size(); ++i) {
          nodes[i]->temp = static_cast<int>(i);
        }

        // Accumlate edges i_right->i_left
        N = nodes.size();
        for (auto n : nodes) {
          switch (n->op()) {
            case OP_CONST:
            case OP_FUNREF:
            case OP_INPUT:
            case OP_OUTPUT:
            case OP_PARAMETER:
              break;
            case OP_CALL:
              i_left.push_back(n->temp);
              i_right.push_back(n->dep(0)->temp);
              i_left.push_back(n->temp);
              i_right.push_back(n->dep(1)->temp);
              break;
            default:
              i_left.push_back(n->temp);
              i_right.push_back(n->dep(0)->temp);
              if (n->n_dep()==2 && n->dep(0)->temp!=n->dep(1)->temp) {
                i_left.push_back(n->temp);
                i_right.push_back(n->dep(1)->temp);
              }
          }
        }
        nc = i_left.size();

        // Build up index
        incoming.resize(N);
        outgoing.resize(N);
        for (int i=0;i<nc;++i) {
          incoming[i_left[i]].push_back(i_right[i]);
          outgoing[i_right[i]].push_back(i_left[i]);
        }
      }


      int max_outgoing(const std::vector<Entry>& x, int i) {
          int m = 0;
          for (int k : outgoing[i]) {
              m = std::max(m, x[k].p-x[i].p-1);
          }
          return m;
      }


      long long int compute_score(const std::vector<Entry>& x) {
          long long int sum = 0;
          if (inf) {
            for (int i=0;i<x.size();++i) {
                sum += max_outgoing(x, i);
            }
          } else {
            for (int i=0;i<nc;++i) {
                int d = x[i_left[i]].p-x[i_right[i]].p-1;
                sum += d;
            }
          }
          return sum;
      }

      void validate(const std::vector<Entry>& x) {
          for (int i=0;i<nc;++i) {
              int d = x[i_left[i]].p-x[i_right[i]].p-1;
              casadi_assert_dev(d>=0);
          }
      }


      bool accept_neighbour(int r, std::vector<Entry>& x) {
          // Considered code location
          Entry& instr_current = x[r];

          if (instr_current.p==0) return false;

          NodeIndex r_prev = instr_current.prev;
          for (auto e : incoming[r]) {
              if (e==r_prev) return false;
          }
          Entry& instr_prev = x[r_prev];
          std::swap(instr_current.p, instr_prev.p);

          if (instr_prev.prev>-1) x[instr_prev.prev].next = r;
          if (instr_current.next>-1) x[instr_current.next].prev = r_prev;

          NodeIndex instr_currrent_next = instr_current.next;

          instr_current.prev = instr_prev.prev;
          instr_current.next = r_prev;
          instr_prev.prev = r;
          instr_prev.next = instr_currrent_next;

          return true;
      }

      bool probe_neighbour(int r, const std::vector<Entry>& x, long long int &score) {
          // Considered code location
          const Entry& instr_current = x[r];

          if (instr_current.p==0) return false;

          NodeIndex r_prev = instr_current.prev;
          for (auto e : incoming[r]) {
              if (e==r_prev) return false;
          }
          const Entry& instr_prev = x[r_prev];

          if (inf) {

            // L1(Linf)
            std::vector<Entry>& x_const = const_cast< std::vector<Entry>& >(x);
            // Substract normal cost
            for (int e : incoming[r]) {
                score -= max_outgoing(x, e);
            }
            for (int e : incoming[r_prev]) {
                score -= max_outgoing(x, e);
            }
            std::swap(x_const[r].p, x_const[r_prev].p);
            for (int e : incoming[r]) {
                score += max_outgoing(x, e);
            }
            for (int e : incoming[r_prev]) {
                score += max_outgoing(x, e);
            }
            std::swap(x_const[r].p, x_const[r_prev].p);
            score += outgoing[r].size()>0;
            score -= outgoing[r_prev].size()>0;
          } else {
            // L1
            score -= incoming[r].size();
            score += outgoing[r].size();
            score += incoming[r_prev].size();
            score -= outgoing[r_prev].size();
          }

          return true;
      }
    
      std::vector<casadi_int> sort() {
        double fac = 1-1e-13;
        double T = 5;

        std::vector<Entry> order(N);
        // Staring order (supplied sequence)
        for (int i=0;i<N;++i) {
            order[i].p = i;
            order[i].prev = i-1;
            order[i].next = i==N-1? -1 : i+1;
        }

        long long int score = compute_score(order);
        long long int score_orig = score;

        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(1); //Standard mersenne_twister_engine seeded with rd()
        casadi_assert(N>1, "Got N = " + str(N));
        std::uniform_int_distribution<> distrib(1, N-1);
        std::uniform_real_distribution<> distrib_real(0, 1);

        validate(order);

        bool stop_signal = false;

        long long int lowest_score = std::numeric_limits<int>::max();
        long long int lowest_k = -1;
        long long int k = 0;

        for (k=0;k<std::numeric_limits<int>::max();++k) {
            //double rem = 1-(k+1)/kmax;
            int r = distrib(gen);
            long long int canidate_score = score;
            // Compute neighbour score
            if (!probe_neighbour(r, order, canidate_score)) continue;
            if (canidate_score < score || exp(-(canidate_score-score)/T) >= distrib_real(gen)) {
                // Accept the neighbour
                accept_neighbour(r, order);
                score = canidate_score;
                // Print best achievement so far
                if (score<lowest_score) {
                  lowest_score = score;
                  lowest_k = k;
                }
                if (stop_signal && score==lowest_score) {
                  break;
                }
            }
            if (k%100000==0) {
              casadi_assert_dev(compute_score(order)==score);
            }
            if (k>=lowest_k+stationary_length) {
              stop_signal = true;
            }
            T *= fac;
        }

        uout() << "score " << std::to_string(lowest_score) << " (orig " << score_orig << ") after " << lowest_k << " iterations" << endl;

        validate(order);
        casadi_assert_dev(compute_score(order)==score);
        
        std::vector<casadi_int> ret;
        ret.reserve(order.size());
        for (const auto& k : order) {
          ret.push_back(k.p);
        }

        return ret;

      }
  };

  std::vector<casadi_int> order_nodes(std::vector<SXNode*>& nodes, long long int stationary_length, bool inf) {
    if (nodes.empty()) return {};
    if (nodes.size()==1) return {0};
    if (nodes.size()==2) return {0, 1};
    /*std::size_t h = 0;
    Sparsity::hash_combine(h, inf);
    Sparsity::hash_combine(h, nodes.size());
    Sparsity::hash_combine(h, stationary_length);*/
    uout() << "Starting ordering computation using simulated annealing " << std::endl;
    uout() << stationary_length << " stationary_length" << std::endl;
    uout() << nodes.size() << " nodes" << std::endl;
    uout() << inf << " inf" << std::endl;
    std::string cache_name = "cache_order_inf" + str(inf) + "_nodes" + str(nodes.size()) + "_stationary_length" + str(stationary_length) + ".txt";
    ifstream in(cache_name);
    if (in.good()) {
      vector<casadi_int> order;
      std::string line;
      while (std::getline(in, line)) {
        order.push_back(stoi(line));
      }
      uout() << "Using cache file" << cache_name << std::endl;
      return order;
    }
    SimulatedAnnealingInstrOrder sa(nodes, stationary_length, inf);
    std::vector<casadi_int> ret = sa.sort();
    ofstream out(cache_name);
    for (casadi_int e : ret) {
      out << e << std::endl;
    }
    return ret;
  }

  void SXFunction::init(const Dict& opts) {
    // Call the init function of the base class
    XFunction<SXFunction, SX, SXNode>::init(opts);
    if (verbose_) casadi_message(name_ + "::init");

    // Default (temporary) options
    live_variables_ = true;

    bool cse_opt = true;

    stride_in_.resize(n_in_, 1);
    stride_out_.resize(n_out_, 1);

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
      } else if (op.first=="stride_in") {
        stride_in_ = op.second;
      } else if (op.first=="stride_out") {
        stride_out_ = op.second;
      }
    }

    if (cse_opt) {
      casadi_int n_before = SX::n_nodes(veccat(out_));
      out_ = cse(out_);
      casadi_int n_after = SX::n_nodes(veccat(out_));
      if (verbose_) casadi_message("cse pass: " + str(n_before) + " -> " + str(n_after));
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

    // Output nodes (output index, nonzero index)
    vector< std::pair<int, int> > outputs;

    // Associates node k with outputs output_indicator[k]
    std::map< SXNode*, std::vector<int> > output_indicator;

    // Add the list of nodes
    casadi_int ind=0;
    for (auto it = out_.begin(); it != out_.end(); ++it, ++ind) {
      casadi_int nz=0;
      for (auto itc = (*it)->begin(); itc != (*it)->end(); ++itc, ++nz) {
        // Add outputs to the list
        s.push(itc->get());
        sort_depth_first(s, nodes);

        // Indicate last added node to have an output association
        output_indicator[itc->get()].push_back(outputs.size());

        // Register output
        outputs.push_back(std::make_pair(ind, nz));
      }
    }

    casadi_assert(nodes.size() <= std::numeric_limits<int>::max(), "Integer overflow");

    vector<casadi_int> order = range(nodes.size());
    if (!all_equal(stride_in_, 1) && nodes.size()>10000) {
      if (GlobalOptions::getSXReordering()=="L1inf") {
        order = order_nodes(nodes, 3*nodes.size()*nodes.size(), true);
      } else if (GlobalOptions::getSXReordering()=="L1") {
        order = order_nodes(nodes, 3*nodes.size()*nodes.size(), false);
      } else if (GlobalOptions::getSXReordering()=="none") {
        // nothing to do
      } else {
        casadi_error("Unrecognised setting " + GlobalOptions::getSXReordering());
      }
    }
    casadi_assert_dev(order.size()==nodes.size());
    casadi_assert_dev(is_permutation(order));

    // Reorder nodes
    vector<SXNode*> nodes_orig = nodes;
    for (int i=0;i<nodes_orig.size();++i) {
      nodes[order[i]] = nodes_orig[i]; 
    }

    // Set the temporary variables to be the corresponding place in the sorted graph
    for (casadi_int i=0; i<nodes.size(); ++i) {
      nodes[i]->temp = static_cast<int>(i);
    }

    // Sort the nodes by type
    constants_.clear();
    operations_.clear();
    for (vector<SXNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
      SXNode* t = *it;
      if (t->is_constant())
        constants_.push_back(SXElem::create(t));
      else if (!t->is_symbolic())
        operations_.push_back(SXElem::create(t));
    }

    // Input instructions
    vector<pair<int, SXNode*> > symb_loc;

    // Count the number of times each node is used
    vector<casadi_int> refcount(nodes.size(), 0);

    size_t sz_call_arg=0, sz_call_res=0, sz_call_iw=0, sz_call_w=0;

    // Get the sequence of instructions for the virtual machine
    algorithm_.resize(0);
    algorithm_.reserve(nodes.size()+outputs.size());
    for (vector<SXNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
      // Current node
      SXNode* n = *it;

      // New element in the algorithm
      AlgEl ae;

      // Get operation
      ae.op = static_cast<int>(n->op());

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
      case OP_FUNREF:
        {
          const FunRef* nn = dynamic_cast<const FunRef*>(n);
          ae.i0 = nn->temp;
          ae.i1 = nn->dep(0).get()->temp;
          ae.i2 = functions_.size();
          functions_.push_back(nn->f_);
          size_t sz_arg, sz_res, sz_iw, sz_w;
          nn->f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);
          sz_call_arg = std::max(sz_call_arg, sz_arg);
          sz_call_res = std::max(sz_call_res, sz_res);
          sz_call_iw = std::max(sz_call_iw, sz_iw);
          sz_call_w = std::max(sz_call_w, sz_w);
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

      // Handle output association // complexity needs fixing
      if (output_indicator.find(n)!=output_indicator.end()) {
        for (int i : output_indicator[n]) {
          const auto& out = outputs[i];
          AlgEl ae;
          ae.op = static_cast<int>(OP_OUTPUT);
          ae.i0 = out.first;
          ae.i1 = n->temp;
          ae.i2 = out.second*stride_out_[ae.i0];

          // Increase count of dependency
          refcount.at(ae.i1)++;

          // Add to algorithm
          algorithm_.push_back(ae);
        }
      }

    }

    nodes_.resize(nodes.size());
    std::copy(nodes.begin(), nodes.end(), nodes_.begin());

    // Place in the work vector for each of the nodes in the tree (overwrites the reference counter)
    vector<int> place(algorithm_.size());

    // Stack with unused elements in the work vector
    stack<int> unused;

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
      if (ndeps==1 && a.op!=OP_OUTPUT && a.op!=OP_CALL && a.op!=OP_FUNREF) {
        a.i2 = a.i1;
      }
    }

    worksize_ = worksize;

    if (verbose_) {
      if (live_variables_) {
        casadi_message("Using live variables: work array is " + str(worksize_)
         + " instead of " + str(algorithm_.size()));
      } else {
        casadi_message("Live variables disabled.");
      }
    }

    // Allocate work vectors (symbolic/numeric)
    alloc_arg(sz_call_arg);
    alloc_res(sz_call_res);
    alloc_iw(sz_call_iw);
    alloc_w(worksize_ + sz_call_w);

    // Reset the temporary variables
    for (casadi_int i=0; i<nodes.size(); ++i) {
      nodes[i]->temp = 0;
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
          algorithm_[i].i2 = nz*std::abs(stride_in_[ind]);

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

    // Do any dependencies need refcount?
    for (auto&& f : functions_) {
      has_refcount_ = has_refcount_ || f->has_refcount_;
    }

    shared_from_this<Function>().save("myfun" + name_ + "stride" + str(sum(stride_in_)) + ".casadi");
  }


  SXElem register_endpoint(const SXElem& expr, std::vector<SXElem>& endpoints, std::map<SXNode*, SXElem>& endpoint_symbols) {
    auto it = endpoint_symbols.find(expr.get());
    if (it==endpoint_symbols.end()) {
      endpoints.push_back(expr);
      SXElem symbol = SXElem::sym("endpoint_"+str(endpoints.size()));
      endpoint_symbols[expr.get()] = symbol;
      return symbol;
    } else {
      return it->second;
    }
  }

  Function SXFunction::pull_out(const std::vector<casadi_int>& in, Function& periphery) const {

    uout() << "pull_out" << name_ << ":" << in << n_in_ << std::endl;
    if (in.empty()) {
      periphery = Function(name_+"_periphery", std::vector<SX>{}, std::vector<SX>{});
      return shared_from_this<Function>();
    }
    casadi_assert(!has_free(), "Not supported for Functions with free parameters");

    // Iterator to the binary operations
    vector<SXElem>::const_iterator b_it=operations_.begin();

    // Iterator to stack of constants
    vector<SXElem>::const_iterator c_it = constants_.begin();

    // Iterator to free variables
    vector<SXElem>::const_iterator p_it = free_vars_.begin();

    vector<SXElem> w(worksize_);

    std::set<casadi_int> ins(in.begin(), in.end());

    std::vector<bool> tainted(worksize_, false);

    // List of endpoints
    std::vector<SXElem> endpoints;

    std::map<SXNode*, SXElem> endpoint_symbols;

    // Tape
    vector<TapeEl<SXElem> > s_pdwork(operations_.size());
    vector<TapeEl<SXElem> >::iterator it1 = s_pdwork.begin();

    std::vector< std::vector<SXElem> > arg(n_in_);
    for (casadi_int i=0;i<n_in_;++i) {
      arg[i] = in_[i].nonzeros();
    }
    std::vector< std::vector<SXElem> > res(n_out_);
    for (casadi_int i=0;i<n_out_;++i) {
      res[i].resize(sparsity_out(i).nnz());
    }
    for (auto&& a : algorithm_) {
      switch (a.op) {
      case OP_INPUT:
        uout() << "in" << a.i1 << a.i2 << std::endl;
        w[a.i0] = arg[a.i1][a.i2/std::abs(stride_in_[a.i1])];
        tainted[a.i0] = !ins.count(a.i1);
        break;
      case OP_OUTPUT:
        {
          SXElem in1 = w[a.i1];
          if (!tainted[a.i1]) {
            in1 = register_endpoint(w[a.i1], endpoints, endpoint_symbols);
          }
          res[a.i0][a.i2] = in1;
        }
        break;
      case OP_CONST:
        w[a.i0] = *c_it++;
        tainted[a.i0] = false;
        break;
      case OP_PARAMETER:
        casadi_error("Cannot occur");
      case OP_FUNREF:
        w[a.i0] = SXElem::create(new FunRef(functions_[a.i2], w[a.i1])); break;
        tainted[a.i0] = false;
      case OP_CALL:
        {
          SXElem in2 = w[a.i2];
          if (!tainted[a.i2]) {
            in2 = register_endpoint(w[a.i2], endpoints, endpoint_symbols);
          }
          w[a.i0] = SXElem::create(new CallSX(w[a.i1], in2));
          tainted[a.i0] = true;
        }
        break;
      default:
        {
          SXElem in1 = w[a.i1];
          SXElem in2;

          bool binary = casadi_math<SXElem>::is_binary(a.op); 
          if (binary) in2 = w[a.i2];

          bool any_tainted = tainted[a.i1];
          bool all_tainted = tainted[a.i1];
          if (binary) {
            any_tainted = any_tainted || tainted[a.i2];
            all_tainted = all_tainted && tainted[a.i2];
          }
          if (any_tainted) {
            if (!tainted[a.i1]) {
              in1 = register_endpoint(w[a.i1], endpoints, endpoint_symbols);
            }
            if (binary && !tainted[a.i2]) {
              in2 = register_endpoint(w[a.i2], endpoints, endpoint_symbols);
            }
          }

          // Evaluate the function to a temporary value
          // (as it might overwrite the children in the work vector)
          SXElem f;
          switch (a.op) {
            CASADI_MATH_FUN_BUILTIN(in1, in2, f)
          }


          // If this new expression is identical to the expression used
          // to define the algorithm, then reuse
          const casadi_int depth = 2; // NOTE: a higher depth could possibly give more savings
          f.assignIfDuplicate(*b_it++, depth);

          // Finally save the function value
          w[a.i0] = f;
          tainted[a.i0] = any_tainted;

        }
      }
    }

    // pass is-diff in
    periphery = Function("periphery_" + name_, vector_slice(in_, in), {SX(endpoints)}, vector_slice(name_in_, in), {name_+"_endpoint"});
    uout() << "periphery" << std::endl;    
    uout() << "in" << vector_slice(in_, in) << std::endl;
    uout() << "res" << SX(endpoints) << std::endl;

    casadi_assert_dev(!periphery.has_free());
    std::vector<casadi_int> in_invert = complement(in, in_.size());

    std::vector<SX> in_syms = vector_slice(in_, in_invert);
    std::vector<SXElem> endpoint_symbols_vec;
    for (auto e : endpoints) {
      endpoint_symbols_vec.push_back(endpoint_symbols[e.get()]);
    }

    in_syms.push_back(SX(endpoint_symbols_vec));
    std::vector<std::string> name_in = vector_slice(name_in_, in_invert);
    name_in.push_back(name_+"_endpoint");
    std::vector<SX> resSX(n_out_);
    for (casadi_int i=0;i<n_out_;++i) {
      resSX[i] = SX(sparsity_out(i), res[i]);
    }

    uout() << "core" << std::endl;
    uout() << "in" << in_syms << std::endl;
    uout() << "res" << resSX << std::endl;
    std::vector<casadi_int> stride_in = vector_slice(stride_in_, in);
    stride_in.push_back(1);
    Dict opts;
    opts["stride_in"] = stride_in;
    opts["stride_out"] = stride_out_;
    Function ret = Function("core_" + name_, in_syms, resSX, name_in, name_out_, opts);
    uout() << ret.get_free() << std::endl;
    casadi_assert(!ret.has_free(), name_);
    return ret;
  }


  void SXFunction::codegen_incref(CodeGenerator& g, const Instance& inst) const {
    for (auto&& f : functions_) {
      if (f->has_refcount_) {
        auto i = g.incref_added_.insert(f.get());
        if (i.second) { // prevent duplicate calls
          Instance local;
          local.stride_in.resize(f.n_in(), 1);
          local.stride_out.resize(f.n_out(), 1);
          g << f->codegen_name(g) << "_incref(); // SXFunction::codegen_incref\n";
        }
      }
    }
  }

  void SXFunction::codegen_decref(CodeGenerator& g, const Instance& inst) const {
    for (auto&& f : functions_) {
      if (f->has_refcount_) {
        auto i = g.decref_added_.insert(f.get());
        if (i.second) { // prevent duplicate calls
          Instance local;
          local.stride_in.resize(f.n_in(), 1);
          local.stride_out.resize(f.n_out(), 1);
          g << f->codegen_name(g) << "_decref();\n";
        }
      }
    }
  }

  SX SXFunction::instructions_sx() const {
    std::vector<SXElem> ret(algorithm_.size(), casadi_limits<SXElem>::nan);

    vector<SXElem>::iterator it=ret.begin();

    // Iterator to the binary operations
    vector<SXElem>::const_iterator b_it = operations_.begin();

    // Iterator to stack of constants
    vector<SXElem>::const_iterator c_it = constants_.begin();

    // Iterator to free variables
    vector<SXElem>::const_iterator p_it = free_vars_.begin();

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
      case OP_FUNREF:
        w[a.i0] = SXElem::create(new FunRef(functions_[a.i2], w[a.i1])); break;
      case OP_CALL:
        w[a.i0] = SXElem::create(new CallSX(w[a.i1], w[a.i2])); break;
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
              CallSX::der(f->dep(0), f->dep(1), f, it1++->d);
              break;
            case OP_FUNREF:
              FunRef::der(f->dep(0), f->dep(1), f, it1++->d);
              break;
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
          w[a.i0] = fseed[dir][a.i1].nonzeros()[a.i2]; break; // / layout_in_[a.i1].stride()
        case OP_OUTPUT:
          fsens[dir][a.i0].nonzeros()[a.i2] = w[a.i1]; break;//  / layout_out_[a.i0].stride()
        case OP_CONST:
        case OP_PARAMETER:
          w[a.i0] = 0;
          break;
        CASADI_MATH_BINARY_BUILTIN // Binary operation
        case OP_CALL:
          w[a.i0] = it2->d[0] * w[a.i1] + it2->d[1] * w[a.i2];
          it2++;
          break;
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
              CallSX::der(f->dep(0), f->dep(1), f, it1++->d);
              break;
            case OP_FUNREF:
              FunRef::der(f->dep(0), f->dep(1), f, it1++->d);
              break;
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
        CASADI_MATH_BINARY_BUILTIN // Binary operation
        case OP_CALL:
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
    // Fall back when forward mode not allowed
    if (sp_weight()==1 || sp_weight()==-1)
      return FunctionInternal::sp_forward(arg, res, iw, w, mem);
    // Propagate sparsity forward
    for (auto&& e : algorithm_) {
      switch (e.op) {
      case OP_CONST:
      case OP_PARAMETER:
      case OP_FUNREF:
        w[e.i0] = 0; break;
      case OP_INPUT:
        w[e.i0] = arg[e.i1]==nullptr ? 0 : arg[e.i1][e.i2];
        break;
      case OP_OUTPUT:
        if (res[e.i0]!=nullptr) res[e.i0][e.i2] = w[e.i1];
        break;
      case OP_CALL:
        w[e.i0] = w[e.i2]; break;
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
    fill_n(w, sz_w(), 0);

    // Propagate sparsity backward
    for (auto it=algorithm_.rbegin(); it!=algorithm_.rend(); ++it) {
      // Temp seed
      bvec_t seed;

      // Propagate seeds
      switch (it->op) {
      case OP_CONST:
      case OP_PARAMETER:
      case OP_FUNREF:
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
        seed = w[it->i0];
        w[it->i0] = 0;
        w[it->i2] |= seed;
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

  std::vector<casadi_int> order_markowitz(const Sparsity& M, casadi_int p) {
    std::vector<casadi_int> order;
    casadi_int m = M.size1()-p;
    casadi_int n = M.size2()-p;
    casadi_int total = 0;
    Matrix<casadi_int> MM(M, 1);
    for (casadi_int ii=0;ii<p;++ii) {
      std::vector<casadi_int> mu = densify(vec(sum1(MM)(0, Slice(n, n+p)))*sum2(MM)(Slice(0, p), 0)).nonzeros();

      for (casadi_int i : order) {
        mu[i] = std::numeric_limits<casadi_int>::max();
      }
      casadi_int i = std::distance(mu.begin(), std::min_element(mu.begin(), mu.end()));

      total += mu[i];
      order.push_back(i);
      MM(kron(MM(i, Slice()), MM(Slice(), n+i)).sparsity()) = 1;
      MM(i, Slice()) = Matrix<casadi_int>(1, 1);
      MM(Slice(), n+i) = Matrix<casadi_int>(1, 1);
    }
    uout() << "total " << total << std::endl;
    return order;
  }

  template <class Ita, class Itb>
  bool is_disjoint(Ita a_begin, Ita a_end, Itb b_begin, Ita b_end) {
    Ita ita = a_begin;
    Itb itb = b_begin;


    // Depending on contents, ita or itb may be already advanced

    while (ita != a_end && itb != b_end) {
      if (*ita == *itb) return false;
      if (*ita < *itb) {
        ita++;
      } else {
        itb++;
      }
    }
    return true;
  }

  void partial_product(const std::vector<  std::vector< std::pair<casadi_int, SXElem*> > > & row_view,
                  const std::vector<  std::vector< std::pair<casadi_int, SXElem*> > > & col_view,
                  casadi_int i, casadi_int j, casadi_int r) {

    auto itrow = row_view[i].begin();
    auto itcol = col_view[j].begin();

    auto endrow = row_view[i].end();
    auto endcol = col_view[j].end();

    // Depending on contents, ita or itb may be already advanced

    SXElem prod = 0;
    bool hit = false;
    while (itrow != endrow && itcol != endcol && itrow->first < r && itcol->first < r) {
      if (itrow->first == itcol->first) {
        prod += (*itrow->second) * (*itcol->second);
        hit = true;
      }
      if (itrow->first  < itcol->first ) {
        itrow++;
      } else {
        itcol++;
      }
    }

    if (!hit) return;

    // Find location (i, j)
    itrow = std::lower_bound(itrow, endrow, j, [](const std::pair<casadi_int, const SXElem*>& v, casadi_int val)
          {
              return v.first < val;
          });

    *(itrow->second) += prod;
  }


  SX SXFunction::jac_ve(const Dict& opts) const {
    const SX& x = in_[0];
    const SX& y = out_[0];
    casadi_int m = y.nnz();
    casadi_int n = x.nnz();
    casadi_int p = algorithm_.size()-nnz_out();

    std::vector< std::set<int> > in_deps(sz_w());
    std::vector< std::set<int> > out_deps(sz_w());

    std::vector<casadi_int> c_row, c_col;
    std::vector<SXElem> c_data;

    // Iterator to the binary operations
    vector<SXElem>::const_iterator b_it=operations_.begin();

    // Evaluate algorithm
    if (verbose_) casadi_message("Evaluating algorithm forward");
    for (auto&& a : algorithm_) {
      switch (a.op) {
      case OP_INPUT:
        c_row.push_back(a.i0);
        c_col.push_back(a.i2);
        c_data.push_back(1);
        in_deps[a.i0].insert(a.i2);
        break;
      case OP_OUTPUT:
        c_row.push_back(p+a.i2);
        c_col.push_back(n+a.i1);
        c_data.push_back(1);
        out_deps[a.i1].insert(a.i2);
        break;
      case OP_CONST:
        break;
      case OP_PARAMETER:
        break;
      default:
        {
          c_row.push_back(a.i0);
          c_col.push_back(n+a.i1);
          const SXElem& f=*b_it++;
          casadi_int ndeps = casadi_math<SXElem>::ndeps(a.op);
          // Place to store partials
          c_data.emplace_back();
          c_data.emplace_back();
          SXElem* d = get_ptr(c_data)+c_data.size()-2;
          in_deps[a.i0] = in_deps[a.i1];
          if (ndeps>1 && a.i1!=a.i2) {
            c_row.push_back(a.i0);
            c_col.push_back(n+a.i2);
            in_deps[a.i0].insert(in_deps[a.i2].begin(), in_deps[a.i2].end());
          }
          switch (a.op) {
            CASADI_MATH_DER_BUILTIN(f->dep(0), f->dep(1), f, d);
            case OP_CALL:
              CallSX::der(f->dep(0), f->dep(1), f, d);
              break;
            case OP_FUNREF:
              FunRef::der(f->dep(0), f->dep(1), f, d);
              break;
          }
          if (ndeps>1 && a.i1==a.i2) {
            d[0] += d[1];
            c_data.pop_back();
          } else if (ndeps==1) {
            c_data.pop_back();
          }
        }
      }
    }

    for (auto it = algorithm_.rbegin(); it!=algorithm_.rend(); ++it) {
      switch (it->op) {
      case OP_INPUT:
      case OP_OUTPUT:
      case OP_CONST:
      case OP_PARAMETER:
      case OP_CALL:
      case OP_FUNREF:
        break;
      default:
        {
          casadi_int ndeps = casadi_math<SXElem>::ndeps(it->op);
          out_deps[it->i1].insert(out_deps[it->i0].begin(), out_deps[it->i0].end());
          if (ndeps>1) {
            out_deps[it->i2].insert(out_deps[it->i0].begin(), out_deps[it->i0].end());
          }
        }
      }
    }

    SX C = SX::triplet(c_row, c_col, c_data, m+p, n+p);

    SX B = C(Slice(0, p), Slice(0, n));
    SX L = C(Slice(0, p), Slice(n, n+p));
    SX R = C(Slice(p, m+p), Slice(0, n));
    SX T = C(Slice(p, m+p), Slice(n, n+p));


    std::vector<casadi_int> perm = order_markowitz(C.sparsity(), p);
    SX P(Sparsity::permutation(perm), 1);

    SX C_star = SX::blockcat({{mtimes(mtimes(P, L), P.T())-SX::eye(p), mtimes(P, B)}, {mtimes(T, P.T()), R}});

    const casadi_int* colind = C_star.sparsity().colind();
    const casadi_int* row = C_star.sparsity().row();

    // Alternative view on sparsity pattern (with perfect redundancy)
    std::vector< std::set<casadi_int> > row_view(C_star.size1()); // List of column indices per row
    std::vector< std::set<casadi_int> > col_view(C_star.size2()); // List of row indices per column
    // Loop over columns
    for (casadi_int c = 0; c < C_star.size2(); ++c) {
      // Loop over nonzeros for the column
      for (casadi_int k = colind[c]; k < colind[c + 1]; ++k) {
        casadi_int r = row[k];
        row_view[r].insert(c);
        col_view[c].insert(r);
      }
    }

    // Traverse all (i,j) in Crout order
    for (casadi_int k=1;k<p+std::min(n, m);++k) {
      casadi_int i=k;
      for (casadi_int j=i;j<p+n;++j) {
        casadi_int r = std::min(std::min(i, j), p);
        bool empty = is_disjoint(row_view[i].begin(), row_view[i].lower_bound(r),
                               col_view[j].begin(), col_view[j].lower_bound(r));
        if (!empty) {
          row_view[i].insert(j);
          col_view[j].insert(i);
        }
      }
      casadi_int j=k;
      for (casadi_int i=j+1;i<p+m;++i) {
        casadi_int r = std::min(std::min(i, j), p);
        bool empty = is_disjoint(row_view[i].begin(), row_view[i].lower_bound(r),
                               col_view[j].begin(), col_view[j].lower_bound(r));
        if (!empty) {
          row_view[i].insert(j);
          col_view[j].insert(i);
        }
      }
    }

    c_row.clear();
    c_col.clear();
    for (casadi_int j=0;j<col_view.size();++j) {
      for (casadi_int i : col_view[j]) {
        c_row.push_back(i);
        c_col.push_back(j);
      }
    }

    Sparsity sp = Sparsity::triplet(m+p, n+p, c_row, c_col);

    C_star = project(C_star, sp);


    {

      colind = C_star.sparsity().colind();
      row = C_star.sparsity().row();
      std::vector<SXElem> & c_data = C_star.nonzeros();
      std::vector< std::vector< std::pair<casadi_int, SXElem*> > > row_view(C_star.size1()); // List of column indices per row
      std::vector< std::vector< std::pair<casadi_int, SXElem*> > > col_view(C_star.size2()); // List of row indices per column
      // Loop over columns
      for (casadi_int c = 0; c < C_star.size2(); ++c) {
        // Loop over nonzeros for the column
        for (casadi_int k = colind[c]; k < colind[c + 1]; ++k) {
          casadi_int r = row[k];
          row_view[r].push_back(std::make_pair(c, &c_data[k]));
          col_view[c].push_back(std::make_pair(r, &c_data[k]));
        }
      }

      // Traverse all (i,j) in Crout order
      for (casadi_int k=1;k<p+std::min(n, m);++k) {
        casadi_int i=k;
        for (casadi_int j=i;j<p+n;++j) {
          casadi_int r = std::min(std::min(i, j), p);
          partial_product(row_view, col_view, i, j, r);
        }
        casadi_int j=k;
        for (casadi_int i=j+1;i<p+m;++i) {
          casadi_int r = std::min(std::min(i, j), p);
          partial_product(row_view, col_view, i, j, r);
        }
      }

    }

/**
    // Traverse all (i,j) in Crout order
    for (casadi_int k=1;k<p+std::min(n, m);++k) {
      casadi_int i=k;
      for (casadi_int j=i;j<p+n;++j) {
        Slice sk(0, std::min(std::min(i, j), p));
        SX prod = mtimes(C_star(i, sk), C_star(sk, j));
        if (prod.nnz()>0) {
          C_star(i, j) += prod;
        }
      }
      casadi_int j=k;
      for (casadi_int i=j+1;i<p+m;++i) {
        Slice sk(0, std::min(std::min(i, j), p));
        SX prod = mtimes(C_star(i, sk), C_star(sk, j));
        if (prod.nnz()>0) {
          C_star(i, j) += prod;
        }
      }
    }
*/
    casadi_assert_dev(sp==C_star.sparsity());

    SX J = sparsify(C_star(Slice(p, m+p), Slice(p, n+p)));
    return cse(J);
  }

  std::vector<SX> SXFunction::jac_alt(const Dict& opts) const {
    Function f("f", in_, out_, {{"live_variables", false}, {"cse", true}});
    const SXFunction* fsx = dynamic_cast<const SXFunction*>(f.get());
    return {fsx->jac_ve(opts)};
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
    s.unpack("SXFunction::stride_in", stride_in_);
    s.unpack("SXFunction::stride_out", stride_out_);

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

    s.unpack("SXFunction::functions", functions_);
    s.unpack("SXFunction::nodes", nodes_);
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
    s.pack("SXFunction::stride_in", stride_in_);
    s.pack("SXFunction::stride_out", stride_out_);
    // Loop over algorithm
    for (const auto& e : algorithm_) {
      s.pack("SXFunction::ScalarAtomic::op", e.op);
      s.pack("SXFunction::ScalarAtomic::i0", e.i0);
      s.pack("SXFunction::ScalarAtomic::i1", e.i1);
      s.pack("SXFunction::ScalarAtomic::i2", e.i2);
    }

    s.pack("SXFunction::live_variables", live_variables_);
    XFunction<SXFunction, SX, SXNode>::delayed_serialize_members(s);

    s.pack("SXFunction::functions", functions_);
    s.pack("SXFunction::nodes", nodes_);
  }

  ProtoFunction* SXFunction::deserialize(DeserializingStream& s) {
    return new SXFunction(s);
  }

} // namespace casadi
