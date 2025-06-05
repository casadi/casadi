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
#include <random>

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
    print_instructions_ = false;
  }

  SXFunction::~SXFunction() {
    clear_mem();
  }

  int SXFunction::eval(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const {
    if (verbose_) casadi_message(name_ + "::eval");
    setup(mem, arg, res, iw, w);

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

    if (print_instructions_) {
      int k = 0;
      // Evaluate the algorithm
      for (auto&& e : algorithm_) {
        print_arg(uout(), k, e, w);
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
        print_res(uout(), k, e, w);
        k++;
      }
    } else {
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
  std::string SXFunction::print(const ScalarAtomic& a) const {
    std::stringstream stream;
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
          stream << free_vars_[a.i1];
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
    return stream.str();
  }

  void SXFunction::disp_more(std::ostream &stream) const {
    stream << "Algorithm:";

    // Normal, interpreted output
    for (auto&& a : algorithm_) {
      InterruptHandler::check();
      stream << std::endl;
      stream << print(a);
      stream << ";";
    }
  }

  size_t SXFunction::codegen_sz_w(const CodeGenerator& g) const {
    if (!g.avoid_stack()) return call_.sz_w+call_.sz_w_arg+call_.sz_w_res;
    return sz_w();
  }

  void SXFunction::codegen_declarations(CodeGenerator& g, const Instance& inst) const {

    // Make sure that there are no free variables
    if (!free_vars_.empty()) {
      casadi_error("Code generation of '" + name_ + "' is not possible since variables "
                   + str(free_vars_) + " are free.");
    }

    // Generate code for the call nodes
    for (auto&& m : call_.el) {
      Instance local;
      const Function& f = m.f;
      local.stride_in.resize(f.n_in(), 1);
      local.stride_out.resize(f.n_out(), 1);
      local.prefer_inline = true;

      bool added = g.has_dependency(f, local);
      std::string fname = g.add_dependency(m.f, local, shared_from_this<Function>());
      if (!added) {

        std::string name = fname+"_wrap";

        //g.register_extra(this,f, inst, "code", inline);

        std::stringstream s;
        g.flush(s);
        g << "#pragma omp declare simd simdlen("<< GlobalOptions::vector_width_real << ") notinbranch\n";
        // 'const' is needed to vectorize succesfully; 'pure' is not enough in some circumstances
        // __attribute__((optimize("-O0")))
        g << "static __attribute__((noinline)) __attribute__((const)) " << g.vector_width_attribute() << " casadi_real " << name << "(casadi_int index, casadi_real a) {\n";
        g << "const casadi_real* arg[" << f.sz_arg() << "];\n";
        g << "casadi_real* res[" << f.sz_res() << "];\n";
        g << "casadi_int iw[" << f.sz_iw() << "];\n";
        g << "casadi_real w[" << f.sz_w() << "]";
        g << " __attribute__((aligned (" << 64 << ")))";
        g << ";\n";
        g << "casadi_real r;\n";
        if (f.nnz_in()==1) { 
          g << "arg[0] = &a;\n";
        } else {
          g << "arg[0] = (casadi_real*) &index;\n";
          g << "arg[1] = &a;\n";
        }
        g << "res[0] = &r;\n";
        g << g(f, "arg", "res", "iw", "w", "1", local) << ";\n";
        g << "return r;\n";
        g << "}\n";

        g.flush(s);

        g.casadi_headers << s.str();

        //g.add_extra_declarations(f, s.str());
        //g.add_extra_definitions(f, s.str());
      }

    }
  }


  void SXFunction::print_arg(std::ostream &stream, casadi_int k, const ScalarAtomic& el,
      const double* w) const {
    if (el.op==OP_INPUT || el.op==OP_OUTPUT || el.op==OP_CONST) return;
    stream << name_ << ":" << k << ": " << print(el) << " inputs:" << std::endl;

    // Default dependencies
    const int* dep = &el.i1;
    casadi_int ndeps = casadi_math<double>::ndeps(el.op);

    // Call node overrides these defaults
    if (el.op==OP_CALL) {
      const ExtendedAlgEl& e = call_.el.at(el.i1);
      ndeps = e.n_dep;
      dep = get_ptr(e.dep);
      stream << "[";
      for (size_t i = 0; i < ndeps; ++i) {
        if (i>0) stream << ", ";
        if (print_canonical_) {
          print_canonical(stream, w[dep[i]]);
        } else {
          DM::print_scalar(stream, w[dep[i]]);
        }
      }
      stream << "]";
      stream << std::endl;
      return;
    }

    for (size_t i = 0; i < ndeps; ++i) {
      stream << i << ": ";
      if (print_canonical_) {
        print_canonical(stream, w[dep[i]]);
      } else {
        DM::print_scalar(stream, w[dep[i]]);
      }
      stream << std::endl;
    }
  }

  void SXFunction::print_arg(CodeGenerator& g, casadi_int k, const ScalarAtomic& el) const {
    if (el.op==OP_INPUT || el.op==OP_OUTPUT || el.op==OP_CONST) return;
    g << g.printf(name_ + ":" + str(k) + ": " + print(el) + " inputs:\\n") << "\n";
    if (el.op==OP_CALL) {
      const ExtendedAlgEl& m = call_.el[el.i1];
      g << g.print_vector(m.f.nnz_in(), "arg[" + str(n_in_) + "]");
      g << g.printf("\\n");
    } else {
      casadi_int ndeps = casadi_math<double>::ndeps(el.op);
      if (ndeps==1) {
        g << g.printf("0: %.16e\\n", g.sx_work(el.i1));
      } else if (ndeps==2) {
        g << g.printf("0: %.16e\\n1: %.16e\\n", g.sx_work(el.i1), g.sx_work(el.i2));
      }
    }
    g << "\n";
  }

  void SXFunction::print_res(CodeGenerator& g, casadi_int k, const ScalarAtomic& el) const {
    if (el.op==OP_INPUT || el.op==OP_OUTPUT) return;
    g << g.printf(name_ + ":" + str(k) + ": " + print(el) + " outputs:\\n") << "\n";
    if (el.op==OP_CALL) {
      const ExtendedAlgEl& m = call_.el[el.i1];
      g << g.print_vector(m.f.nnz_out(), "w+" + str(m.f.nnz_in()));
      g << g.printf("\\n");
    } else {
      g << g.printf("0: %.16e\\n", g.sx_work(el.i0));
    }
    g << "\n";
  }

  void SXFunction::print_res(std::ostream &stream, casadi_int k, const ScalarAtomic& el,
      const double* w) const {
    if (el.op==OP_INPUT || el.op==OP_OUTPUT) return;
    stream << name_ << ":" << k << ": " << print(el) << " outputs:" << std::endl;

    // Default outputs
    const int* res = &el.i0;
    casadi_int nres = 1;

    // Call node overrides these defaults
    if (el.op==OP_CALL) {
      const ExtendedAlgEl& e = call_.el.at(el.i1);
      nres = e.n_res;
      res = get_ptr(e.res);
      stream << "[";
      for (size_t i = 0; i < nres; ++i) {
        if (i>0) stream << ", ";
        if (print_canonical_) {
          print_canonical(stream, w[res[i]]);
        } else {
          DM::print_scalar(stream, w[res[i]]);
        }
      }
      stream << "]";
      stream << std::endl;
      return;
    }

    for (size_t i = 0; i < nres; ++i) {
      stream << i << ": ";
      if (print_canonical_) {
        print_canonical(stream, w[res[i]]);
      } else {
        DM::print_scalar(stream, w[res[i]]);
      }
      stream << std::endl;
    }

  }

  void SXFunction::codegen_body(CodeGenerator& g, const Instance& inst) const {
    g.reserve_work(worksize_);
    std::map<const AlgEl*, casadi_int> lookup_rets, lookup_ins;
    std::vector< const AlgEl* > rets, ins;

    bool vectorize = false;
    for (auto e : inst.stride_in) {
      if (e>1) vectorize = true;
    }
    for (auto e : inst.stride_out) {
      if (e>1) vectorize = true;
    }

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

    casadi_int cnt = 0;
    // Run the algorithm
    for (auto&& a : algorithm_) {
      if (a.op==OP_OUTPUT) {
        if (inst.res_null.empty()) {
          g << "if (res[" << a.i0 << "]!=0) ";
          int stride = inst.stride_out.empty() ? 1 : inst.stride_out.at(a.i0);
          stride = 1;
          g << g.res(a.i0, true) << "[" << a.i2*abs(stride) << "]" << (stride<0? "+": "")<<  "=" << g.sx_work(a.i1)  << ";\n";
        } else {
          if (!inst.res_null[a.i0]) {
            if (outro) {
              g << "ret" << lookup_rets[&a] << " = " << g.sx_work(a.i1) << ";\n";
            } else {
              int stride = inst.stride_out.empty() ? 1 : inst.stride_out.at(a.i0);
              stride = 1;
              g << g.res(a.i0, true) << "[" << a.i2*abs(stride) << "]" << " =" << g.sx_work(a.i1) << ";\n";
            }
          }
        }
      } else if (a.op==OP_CALL) {
        const ExtendedAlgEl& m = call_.el[a.i1];

        casadi_int worksize = g.avoid_stack() ? worksize_ : 0;

        // Collect input arguments
        /*casadi_int offset = worksize;
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
        }*/

        if (print_instructions_) print_arg(g, cnt, a);

        Instance local;
        local.prefer_inline = true;
        local.stride_in.resize(m.f.n_in(), 1);
        local.stride_out.resize(m.f.n_out(), 1);
        std::string fname = g.add_dependency(m.f, local, shared_from_this<Function>());
        std::string name = fname+"_wrap";

        casadi_assert(m.f.n_in() == 2, "The SXFunction wrapper only supports functions with two inputs at the moment.");
        casadi_assert(m.f.n_out() == 1, "The SXFunction wrapper only supports functions with one output at the moment.");

        g << g.sx_work(m.res[0]) << "=" + name + "(casadi_real2int(" + g.sx_work(m.dep[0]) + ")," + g.sx_work(m.dep[1]) + ");\n";

        if (print_instructions_) print_res(g, cnt, a);

        //return g.sx_work(i0) + "=" + name + "(casadi_real2int(" + g.sx_work(i1) + ")," + g.sx_work(i2) + ")";
        /*std::string flag =
          g(m.f, "arg+"+str(n_in_), "res+"+str(n_out_), "iw", "w+" + str(offset));
        // Call function
        g << "if (" << flag << ") return 1;\n";
        if (print_instructions_) print_res(g, cnt, a);
        for (casadi_int i=0;i<m.n_res;++i) {
          if (m.res[i]>=0) {
            g << g.sx_work(m.res[i]) << " = ";
            g << "w[" + str(i+out_offset) + "];\n";
          }
        }*/
      } else if (a.op==OP_INPUT) {
        // Where to store the result
        g << g.sx_work(a.i0) << "=";
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
        g << ";\n";
      } else {
        if (print_instructions_) print_arg(g, cnt, a);

        // Where to store the result
        g << g.sx_work(a.i0) << "=";

        // What to store
        if (a.op==OP_CONST) {
          g << g.constant(a.d);
        } else {
          casadi_int ndep = casadi_math<double>::ndeps(a.op);
          casadi_assert_dev(ndep>0);
          if (ndep==1) g << g.print_op(a.op, g.sx_work(a.i1)); // a.op OP_NEG
          if (ndep==2) g << g.print_op(a.op, g.sx_work(a.i1), g.sx_work(a.i2)); // OP_MUL OP_ADD
        }

        g << ";\n";

        if (print_instructions_) print_res(g, cnt, a);
      }
      cnt++;
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
      {"allow_free",
       {OT_BOOL,
        "Allow construction with free variables (Default: false)"}},
      {"allow_duplicate_io_names",
       {OT_BOOL,
        "Allow construction with duplicate io names (Default: false)"}},
      {"print_instructions",
       {OT_BOOL,
        "Print each operation during evaluation. Influenced by print_canonical."}},
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

  Dict SXFunction::generate_options(const std::string& target, bool keep_dim) const {
    Dict opts = FunctionInternal::generate_options(target, keep_dim);
    //opts["default_in"] = default_in_;
    opts["live_variables"] = live_variables_;
    opts["just_in_time_sparsity"] = just_in_time_sparsity_;
    opts["just_in_time_opencl"] = just_in_time_opencl_;
    opts["print_instructions"] = print_instructions_;
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

        uout() << "score " << std::to_string(lowest_score) << " (orig " << score_orig << ") after " << lowest_k << " iterations" << std::endl;

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
    std::ifstream in(cache_name);
    if (in.good()) {
      std::vector<casadi_int> order;
      std::string line;
      while (std::getline(in, line)) {
        order.push_back(stoi(line));
      }
      uout() << "Using cache file" << cache_name << std::endl;
      return order;
    }
    SimulatedAnnealingInstrOrder sa(nodes, stationary_length, inf);
    std::vector<casadi_int> ret = sa.sort();
    std::ofstream out(cache_name);
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
    bool allow_free = false;

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
      } else if (op.first=="allow_free") {
        allow_free = op.second;
      } else if (op.first=="print_instructions") {
        print_instructions_ = op.second;
      } else if (op.first=="stride_in") {
        stride_in_ = op.second;
      } else if (op.first=="stride_out") {
        stride_out_ = op.second;
      }
    }

    // Perform common subexpression elimination
    // This must be done before the lock, to avoid deadlocks
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

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    std::lock_guard<std::mutex> lock(SX::get_mutex_temp());
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    // Stack used to sort the computational graph
    std::stack<SXNode*> s;

    // All nodes
    std::vector<SXNode*> nodes;

    // Output nodes (output index, nonzero index)
    std::vector< std::pair<int, int> > outputs;

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

    std::vector<casadi_int> order = range(nodes.size());
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
    std::vector<SXNode*> nodes_orig = nodes;
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

    // Count the number of times each node is used
    std::vector<casadi_int> refcount(nodes.size(), 0);

    size_t sz_call_arg=0, sz_call_res=0, sz_call_iw=0, sz_call_w=0;

    // Get the sequence of instructions for the virtual machine
    algorithm_.resize(0);
    algorithm_.reserve(nodes.size()+outputs.size());

    // Mapping of node index (cfr. temp) to algorithm index
    std::vector<int> alg_index;
    alg_index.reserve(nodes.size()+outputs.size());

    for (std::vector<SXNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
      // Current node
      SXNode* n = *it;

      // New element in the algorithm
      AlgEl ae;

      // Get operation
      ae.op = static_cast<int>(n->op());

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

      // so OP_OUTPUT has two appearances in the algorithm
      // First as fake 'Unary or binary operation' than for real.

      // Increase count of dependencies
      for (casadi_int c=0; c<ndeps; ++c) {
        refcount.at(dep[c])++;
      }

      // Amend node index to algorithm index mapping
      alg_index.push_back(algorithm_.size());

      // Add to algorithm
      if (ae.op>=0) algorithm_.push_back(ae);

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

    // Place in the work vector for each of the nodes in the tree (overwrites the reference counter)
    std::vector<int> place(std::max(algorithm_.size(), nodes.size()));

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
      // if (ndeps==1 && a.op!=OP_OUTPUT && a.op!=OP_CALL && a.op!=OP_FUNREF
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

    alloc_arg(call_.sz_arg, true);
    alloc_res(call_.sz_res, true);
    alloc_iw(call_.sz_iw, true);
    alloc_w(call_.sz_w+call_.sz_w_arg+call_.sz_w_res, true);

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
    for (std::vector<std::pair<int, SXNode*> >::const_iterator it=symb_loc.begin();
         it!=symb_loc.end(); ++it) {
      if (it->second->temp!=0) {
        // Store the index into free_vars
        algorithm_[it->first].i1 = free_vars_.size();

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

    // Do any dependencies need refcount?
    for (auto&& m : call_.el) {
      has_refcount_ = has_refcount_ || m.f->has_refcount_;
    }
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
    std::vector<SXElem>::const_iterator b_it=operations_.begin();

    // Iterator to stack of constants
    std::vector<SXElem>::const_iterator c_it = constants_.begin();

    // Iterator to free variables
    std::vector<SXElem>::const_iterator p_it = free_vars_.begin();

    std::vector<SXElem> w(worksize_);

    std::set<casadi_int> ins(in.begin(), in.end());

    std::vector<bool> tainted(worksize_, false);

    // List of endpoints
    std::vector<SXElem> endpoints;

    std::map<SXNode*, SXElem> endpoint_symbols;

    // Tape
    std::vector<TapeEl<SXElem> > s_pdwork(operations_.size());
    std::vector<TapeEl<SXElem> >::iterator it1 = s_pdwork.begin();

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
      case OP_CALL:
        casadi_error("Not implemented");
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
        // Generate code for the call nodes
    for (auto&& m : call_.el) {
      const Function& f = m.f;
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
    for (auto&& m : call_.el) {
      const Function& f = m.f;
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
          w[a.i0] = fseed[dir][a.i1].nonzeros()[a.i2]; break; // / layout_in_[a.i1].stride()
        case OP_OUTPUT:
          fsens[dir][a.i0].nonzeros()[a.i2] = w[a.i1]; break;//  / layout_out_[a.i0].stride()
        case OP_CONST:
        case OP_PARAMETER:
          w[a.i0] = 0;
          break;
        case OP_CALL:
          {
            auto& m = call_.el.at(a.i1);
            CallSX* call_node = static_cast<CallSX*>(it2->d[0].get());

            // Construct forward sensitivity function
            Function ff = m.f.forward(1);

            // Symbolic inputs to forward sensitivity function
            std::vector<SXElem> deps;
            deps.reserve(2*m.n_dep);

            // Set nominal inputs from node
            casadi_int offset = 0;
            for (casadi_int i=0;i<m.f_n_in;++i) {
              casadi_int nnz = ff.nnz_in(i);
              casadi_assert(nnz==0 || nnz==m.f.nnz_in(i), "Not implemented");
              for (casadi_int j=0;j<nnz;++j) {
                deps.push_back(call_node->dep(offset+j));
              }
              offset += m.f_nnz_in[i];
            }

            // Do not set nominal outputs
            offset = 0;
            for (casadi_int i=0;i<m.f_n_out;++i) {
              casadi_int nnz = ff.nnz_in(i+m.f_n_in);
              casadi_assert(nnz==0 || nnz==m.f.nnz_out(i), "Not implemented");
              for (casadi_int j=0;j<nnz;++j) {
                deps.push_back(call_node->get_output(offset+j));
              }
              offset += m.f_nnz_out[i];
            }

            // Read in forward seeds from work vector
            offset = 0;
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
          w[a.i0] = it2->d[0] * w[a.i1]; it2++;
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
        case OP_CALL:
          {
            auto& m = call_.el.at(it->i1);
            CallSX* call_node = static_cast<CallSX*>(it2->d[0].get());

            // Construct reverse sensitivity function
            Function fr = m.f.reverse(1);

            // Symbolic inputs to reverse sensitivity function
            std::vector<SXElem> deps;
            deps.reserve(m.n_dep+m.n_res);

            // Set nominal inputs from node
            casadi_int offset = 0;
            for (casadi_int i=0;i<m.f_n_in;++i) {
              casadi_int nnz = fr.nnz_in(i);
              casadi_assert(nnz==0 || nnz==m.f.nnz_in(i), "Not implemented");
              for (casadi_int j=0;j<nnz;++j) {
                deps.push_back(call_node->dep(offset+j));
              }
              offset += m.f_nnz_in[i];
            }

            // Do not set nominal outputs
            offset = 0;
            for (casadi_int i=0;i<m.f_n_out;++i) {
              casadi_int nnz = fr.nnz_in(i+m.f_n_in);
              casadi_assert(nnz==0 || nnz==m.f.nnz_out(i), "Not implemented");
              for (casadi_int j=0;j<nnz;++j) {
                deps.push_back(call_node->get_output(offset+j));
              }
              offset += m.f_nnz_out[i];
            }

            // Read in reverse seeds from work vector
            offset = 0;
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
    std::vector<SXElem>::const_iterator b_it=operations_.begin();

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
          const SXElem& f=*b_it++;

          // Default dependencies
          const int* dep = &a.i1;
          casadi_int ndeps = casadi_math<double>::ndeps(a.op);

          // Default outputs
          const int* res = &a.i0;
          casadi_int nres = 1;

          // Call node overrides these defaults
          if (a.op==OP_CALL) {
            const ExtendedAlgEl& e = call_.el.at(a.i1);
            ndeps = e.n_dep;
            dep = get_ptr(e.dep);
            casadi_assert(e.n_res==1, "Not implemented for n_res>1");
            nres = e.n_res;
            res = get_ptr(e.res);
          }

          // Place to store partials; note that unary operations assume storage space for 2
          std::vector<SXElem> d(std::max(ndeps, static_cast<casadi_int>(2)));

          // Set input dependencies
          for (casadi_int i=0; i<ndeps; ++i) {
            in_deps[a.i0].insert(in_deps[dep[i]].begin(), in_deps[dep[i]].end());
          }

          // Determine unique groups
          std::map<int, std::set<int> > groups;

          std::vector<int> order;
          for (casadi_int i=0; i<ndeps; ++i) {
            //if (a.op==OP_CALL) {
            //  casadi_assert_dev(groups[dep[i]].empty());
            //}
            size_t before = groups[dep[i]].size();
            groups[dep[i]].insert(i);
            size_t after = groups[dep[i]].size();
            if (after>before) {
              order.push_back(dep[i]);
            }
          }

          switch (a.op) {
            CASADI_MATH_DER_BUILTIN(f->dep(0), f->dep(1), f, get_ptr(d));
            case OP_CALL:
              {
                const ExtendedAlgEl& m = call_.el.at(a.i1);
                CallSX* call_node = static_cast<CallSX*>(f.get());

                // Construct reverse sensitivity function
                std::vector<std::string> out;
                for (casadi_int i=0;i<m.f_n_in;++i) {
                  if (m.f->is_diff_in_[i]) {
                    out.push_back("jac:"+m.f.name_out(0)+":"+m.f.name_in(i));
                  }
                }

                Function fr = m.f.factory("der_"+m.f.name(), m.f.name_in(), out);

                // Symbolic inputs to reverse sensitivity function
                std::vector<SXElem> deps;
                deps.reserve(m.n_dep);

                // Set nominal inputs from node
                casadi_int offset = 0;
                for (casadi_int i=0;i<m.f_n_in;++i) {
                  casadi_int nnz = fr.nnz_in(i);
                  casadi_assert(nnz==0 || nnz==m.f.nnz_in(i), "Not implemented");
                  for (casadi_int j=0;j<nnz;++j) {
                    deps.push_back(call_node->dep(offset+j));
                  }
                  offset += m.f_nnz_in[i];
                }

                // Call reverse sensitivity function
                std::vector<SXElem> ret = SXElem::call(fr, deps);

                // Store reverse sensitivities into work vector
                offset = 0;
                casadi_int k = 0;
                casadi_int ii = 0;
                for (casadi_int i=0;i<m.f_n_in;++i) {
                  if (m.f->is_diff_in_[i]) {
                    casadi_int nnz = fr.nnz_out(ii);
                    // nnz=0 occurs for is_diff_in[i] false
                    casadi_assert(nnz==0 || nnz==m.f_nnz_in[ii], "Not implemented");
                    if (nnz) {
                      for (casadi_int j=0;j<nnz;++j) {
                        d[offset+j] = ret[k++];
                      }
                    }
                    ii++;
                  }
                  offset += m.f_nnz_in[i];
                }
              }
          }
          for (const auto& e : order) {
            c_row.push_back(res[0]);
            c_col.push_back(n+e);
            SXElem r = 0;
            for (casadi_int i : groups[e]) {
              r += d[i];
            }
            c_data.push_back(r);
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
        break;
      default:
        {
          // Default dependencies
          const int* dep = &it->i1;
          casadi_int ndeps = casadi_math<double>::ndeps(it->op);

          // Call node overrides these defaults
          if (it->op==OP_CALL) {
            const ExtendedAlgEl& e = call_.el.at(it->i1);
            ndeps = e.n_dep;
            dep = get_ptr(e.dep);
          }

          for (casadi_int i=0; i<ndeps; ++i) {
            out_deps[dep[i]].insert(out_deps[it->i0].begin(), out_deps[it->i0].end());
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
    Function f("f", in_, out_, {{"live_variables", false}, {"allow_free", true}, {"cse", true}});
    const SXFunction* fsx = dynamic_cast<const SXFunction*>(f.get());
    return {fsx->jac_ve(opts)};
  }

  const SX SXFunction::sx_in(casadi_int ind) const {
    return in_.at(ind);
  }

  const std::vector<SX> SXFunction::sx_in() const {
    return in_;
  }

  std::vector<std::string> SXFunction::get_function() const {
    std::map<std::string, bool> flagged;
    for (auto&& a : algorithm_) {
      if (a.op==OP_CALL) {
        const auto& m = call_.el.at(a.i1);
        const Function &f = m.f;
        if (flagged.find(f.name())==flagged.end()) {
          flagged[f.name()] = true;
        }
      }
    }
    std::vector<std::string> ret;
    for (auto it : flagged) {
      ret.push_back(it.first);
    }
    return ret;
  }

  const Function& SXFunction::get_function(const std::string &name) const {
    for (auto&& a : algorithm_) {
      if (a.op==OP_CALL) {
        const auto& m = call_.el.at(a.i1);
        const Function &f = m.f;
        if (name==f.name()) return f;
      }
    }
    casadi_error("No such function '" + name + "'.");
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
    int version = s.version("SXFunction", 1, 3);
    size_t n_instructions;
    s.unpack("SXFunction::n_instr", n_instructions);

    s.unpack("SXFunction::worksize", worksize_);
    s.unpack("SXFunction::free_vars", free_vars_);
    s.unpack("SXFunction::operations", operations_);
    s.unpack("SXFunction::constants", constants_);
    s.unpack("SXFunction::default_in", default_in_);
    s.unpack("SXFunction::stride_in", stride_in_);
    s.unpack("SXFunction::stride_out", stride_out_);

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
    if (version>=3) {
      s.unpack("SXFunction::print_instructions", print_instructions_);
    } else {
      print_instructions_ = false;
    }

    XFunction<SXFunction, SX, SXNode>::delayed_deserialize_members(s);
  }

  void SXFunction::serialize_body(SerializingStream &s) const {
    XFunction<SXFunction, SX, SXNode>::serialize_body(s);
    s.version("SXFunction", 3);
    s.pack("SXFunction::n_instr", algorithm_.size());

    s.pack("SXFunction::worksize", worksize_);
    s.pack("SXFunction::free_vars", free_vars_);
    s.pack("SXFunction::operations", operations_);
    s.pack("SXFunction::constants", constants_);
    s.pack("SXFunction::default_in", default_in_);
    s.pack("SXFunction::stride_in", stride_in_);
    s.pack("SXFunction::stride_out", stride_out_);

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
    s.pack("SXFunction::print_instructions", print_instructions_);

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

  void SXFunction::change_option(const std::string& option_name,
      const GenericType& option_value) {
    if (option_name == "print_instructions") {
      print_instructions_ = option_value;
    } else {
      // Option not found - continue to base classes
      XFunction<SXFunction, SX, SXNode>::change_option(option_name, option_value);
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
