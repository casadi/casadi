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


#ifndef CASADI_CALL_SX_HPP
#define CASADI_CALL_SX_HPP

#include "sx_node.hpp"
#include "function.hpp"
#include <deque>
#include <stack>

/// \cond INTERNAL
namespace casadi {

class CallSX : public SXNode {
  private:

    /** \brief  Constructor is private, use "create" below */
    CallSX(const Function& f, const std::vector<SXElem>& dep) :
        f_(f), dep_(dep) {}

  public:

    /** \brief  Create a binary expression */
    inline static SXElem create(const Function& f, const std::vector<SXElem>& dep) {
      casadi_assert(f.nnz_in()==dep.size(),
        "CallSX::create(f,dep): dimension mismatch: " + str(f.nnz_in()) + " vs " + str(dep.size()));
      return SXElem::create(new CallSX(f, dep));
    }

    /** \brief Destructor
    This is a rather complex destructor which is necessary since the default destructor
    can cause stack overflow due to recursive calling.
    */
    ~CallSX() override {
      for (auto & d  : dep_)
        safe_delete(d.assignNoDelete(casadi_limits<SXElem>::nan));
    }

    // Class name
    std::string class_name() const override {return "CallSX";}

    bool is_op(casadi_int op) const override { return OP_CALL==op; }

    /** \brief  Number of dependencies */
    casadi_int n_dep() const override { return dep_.size();}

    /** \brief  get the reference of a dependency */
    const SXElem& dep(casadi_int i) const override { return dep_[i];}
    SXElem& dep(casadi_int i) override { return dep_[i];}

    /** \brief  Get the operation */
    casadi_int op() const override { return OP_CALL; }

    /** \brief  Print expression */
    std::string print(const std::string& arg1, const std::string& arg2) const override {
      return "call";
    }

    Function f_;

    /** \brief  The dependencies of the node */
    std::vector<SXElem> dep_;

    void serialize_node(SerializingStream& s) const override {
      s.pack("CallSX::f", f_);
      s.pack("CallSX::dep", dep_);
    }

    static SXNode* deserialize(DeserializingStream& s) {
      std::vector<SXElem> dep;
      Function f;
      s.unpack("CallSX::f", f);
      s.unpack("CallSX::dep", dep);
      return new CallSX(f, dep);
    }

};


} // namespace casadi
/// \endcond

#endif // CASADI_CALL_SX_HPP
