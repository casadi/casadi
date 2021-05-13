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


#include "function.hpp"
#include "call_sx.hpp"
#include "fun_ref.hpp"
#include "code_generator.hpp"
#include "global_options.hpp"

namespace casadi {

  CallSX::CallSX(const SXElem& ref, const SXElem& arg) :
    ref_(ref),
    arg_(arg) {

    }


  int CallSX::call(const Function &f, double& result, double dep, double index, const double** arg, double** res, casadi_int* iw, double* w) {

    if (f.nnz_in()==2) {
      arg[0] = &index;
      arg[1] = &dep;
    } else {
      arg[0] = &dep;
    }
    res[0] = &result;
    int mem = f.checkout();
    int ret = f(arg, res, iw, w, mem);
    f.release(mem);
    return ret;
  }

  void CallSX::der(const SXElem& x, const SXElem& y, const SXElem& f, SXElem* d) {
    const FunRef* n = dynamic_cast<const FunRef*>(f.dep(0).get());
    const Function& orig = n->f_;
    if (orig.nnz_in()==2) {
      Function ff = orig.factory("der_"+orig.name(),{orig.name_in(0),orig.name_in(1)},{"jac:"+orig.name_out(0)+":"+orig.name_in(1)});
      d[0] = 0;
      d[1] = SXElem::apply(ff, f.dep(1), f.dep(0).dep(0));
    } else {
      Function ff = orig.factory("der_"+orig.name(),{orig.name_in(0)},{"jac:"+orig.name_out(0)+":"+orig.name_in(0)});
      d[0] = 0;
      d[1] = SXElem::apply(ff, f.dep(1));
    }
  }

  std::string CallSX::codegen(CodeGenerator& g, const SXElem& funref, const Instance& inst, int i0, int i1, int i2, const std::string& arg, const std::string& res, const std::string& iw, const std::string& w) {
    const FunRef* n = dynamic_cast<const FunRef*>(funref.get());
    const Function& f = n->f_;
    Instance local = inst;
    local.stride_in.resize(f.n_in(), 1);
    local.stride_out.resize(f.n_out(), 1);
    /*g << "#pragma omp ordered simd\n";
    g << "{\n";
    g << "(" << arg << ")[0] = &" << g.sx_work(i1) << ";\n";
    g << "(" << res << ")[0] = &" << g.sx_work(i0) << ";\n";
    g << g(f_, arg, res, iw, w);
    g << ";}";*/

    std::string fname = g.add_dependency(f, local);
    std::string name = fname+"_wrap";
    return g.sx_work(i0) + "=" + name + "(" + g.sx_work(i1) + "," + g.sx_work(i2) + ")";
  }

  void CallSX::codegen_dependency(CodeGenerator& g, const Function& f, const Instance& inst) {
    bool added = g.has_dependency(f, inst);
    if (!added) {
      std::string fname = g.add_dependency(f, inst);

      std::string name = fname+"_wrap";

      g << "#pragma omp declare simd simdlen("<< GlobalOptions::vector_width_real << ")\n";
      // 'const' is needed to vectorize succesfully; 'pure' is not enough in some circumstances
      g << "static __attribute__((noinline)) __attribute__((const)) casadi_real " << name << "(casadi_real index, casadi_real a) {\n";
      g << "const casadi_real* arg[" << f.sz_arg() << "];\n";
      g << "casadi_real* res[" << f.sz_res() << "];\n";
      g << "casadi_int iw[" << f.sz_iw() << "];\n";
      g << "casadi_real w[" << f.sz_w() << "];\n";
      g << "casadi_real r;\n";
      if (f.nnz_in()==1) { 
        g << "arg[0] = &a;\n";
      } else {
        g << "arg[0] = &index;\n";
        g << "arg[1] = &a;\n";
      }
      g << "res[0] = &r;\n";
      g << g(f, "arg", "res", "iw", "w",  inst) << ";\n";
      g << "return r;\n";
      g << "}\n";
    }
  }

  std::string CallSX::print(const std::string& arg1, const std::string& arg2) const {
    return arg1 + "->call(" + arg2 +  ")";
  }

  void CallSX::serialize_node(SerializingStream& s) const {
    s.pack("CallSX::ref", ref_);
    s.pack("CallSX::arg", arg_);
  }

  SXNode* CallSX::deserialize(DeserializingStream& s) {
    SXElem ref, arg;
    s.unpack("CallSX::ref", ref);
    s.unpack("CallSX::arg", arg);
    return new CallSX(ref, arg);
  }

} // namespace casadi

