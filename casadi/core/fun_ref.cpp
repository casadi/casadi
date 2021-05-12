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
#include "fun_ref.hpp"
#include "code_generator.hpp"
#include "global_options.hpp"

namespace casadi {

  FunRef::FunRef(const Function &f, const SXElem& dep) :
    dep_(dep),
    f_(f) {
      //REMOVE uout() << f << std::endl;
      casadi_assert_dev(f.n_out()==1);
      casadi_assert_dev(f.nnz_out()==1);
      if (f.n_in()==1) {
        casadi_assert_dev(f.nnz_in(0)==1);
        has_dep_ = false;
      } else if (f.n_in()==2) {
        has_dep_ = true;
        casadi_assert_dev(f.nnz_in(0)==1);
        casadi_assert_dev(f.nnz_in(1)==1);
      }
    }

  double FunRef::pack(int id, double index) {
    casadi_assert((char) index == (casadi_int) index, "Funref index overflow");
    casadi_assert((char) id == id, "Funref id overflow");
    double ret = 0;
    char * result_parts = reinterpret_cast<char*>(&ret);
    result_parts[0] = id;
    result_parts[1] = (char) index;
    return ret;
  }

  void FunRef::unpack(double v, int& id, double& index) {
    char * result_parts = reinterpret_cast<char*>(&v);
    id = result_parts[0];
    index = result_parts[1];
  }

  void FunRef::der(const SXElem& x, const SXElem& y, const SXElem& f, SXElem* d) {
    d[0] = 1;
    d[1] = 0;
  }

  std::string FunRef::print(const std::string& arg1, const std::string& arg2) const {
    return f_.name() + "(" + arg1 + ")";
  }

  void FunRef::serialize_node(SerializingStream& s) const {
    s.pack("FunRef::dep", dep_);
    s.pack("FunRef::f", f_);
  }

  SXNode* FunRef::deserialize(DeserializingStream& s) {
    SXElem dep;
    Function f;
    s.unpack("FunRef::dep", dep);
    s.unpack("FunRef::f", f);
    return new FunRef(f, dep);
  }

} // namespace casadi

