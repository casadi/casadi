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
#include "code_generator.hpp"

namespace casadi {

  CallSX::CallSX(const Function &f, const SXElem& dep) :
    dep_(dep),
    f_(f) {
      casadi_assert_dev(f.n_in()==1);
      casadi_assert_dev(f.n_out()==1);
      casadi_assert_dev(f.nnz_in()==1);
      casadi_assert_dev(f.nnz_out()==1);
    }


  int CallSX::fcn(double& result, double x, const double** arg, double** res, casadi_int* iw, double* w) const {
    arg[0] = &x;
    res[0] = &result;
    int mem = f_.checkout();
    int ret = f_(arg, res, iw, w, mem);
    f_.release(mem);
    return ret;
  }

  SXElem CallSX::fcn(const SXElem& arg) const {
    return SXElem::apply(f_, arg);
  }

  void CallSX::der(const SXElem& x, const SXElem& y, const SXElem& f, SXElem* d) {
    const CallSX* n = dynamic_cast<CallSX*>(f.get());
    const Function& orig = n->f_;
    Function ff = orig.factory("der_"+orig.name(),{orig.name_in(0)},{"jac:"+orig.name_out(0)+":"+orig.name_in(0)});
    d[0] = SXElem::apply(ff, x);
    d[1] = 0;
  }

  std::string CallSX::codegen(CodeGenerator& g, int i0, int i1, const std::string& arg, const std::string& res, const std::string& iw, const std::string& w)const {
    g << "(" << arg << ")[0] = &" << g.sx_work(i1) << ";\n";
    g << "(" << res << ")[0] = &" << g.sx_work(i0) << ";\n";
    g << g(f_, arg, res, iw, w);
    return "";
  }

  void CallSX::codegen_dependency(CodeGenerator& g) const {
    g.add_dependency(f_);
  }

  std::string CallSX::print(const std::string& arg1, const std::string& arg2) const {
    return f_.name() + "(" + arg1 + ")";
  }

  void CallSX::serialize_node(SerializingStream& s) const {
    s.pack("CallSX::dep", dep_);
    s.pack("CallSX::f", f_);
  }

  SXNode* CallSX::deserialize(DeserializingStream& s) {
    SXElem dep;
    Function f;
    s.unpack("CallSX::dep", dep);
    s.unpack("CallSX::f", f);
    return new CallSX(f, dep);
  }

} // namespace casadi

