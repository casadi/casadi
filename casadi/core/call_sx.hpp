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
#include "output_sx.hpp"
#include "function.hpp"
#include <deque>
#include <stack>

// Caching of outputs requires a map
#include <unordered_map>
#define CACHING_MAP std::unordered_map

/// \cond INTERNAL
namespace casadi {

class CASADI_EXPORT CallSX : public SXNode {
  private:

    /** \brief  Constructor is private, use "create" below

        \identifier{28o} */
    CallSX(const Function& f, const std::vector<SXElem>& dep) :
        f_(f), dep_(dep) {}

    mutable WeakCache<casadi_int, SharedSXElem> cache_;
  public:

    /** \brief  Create a binary expression

        \identifier{28p} */
    inline static SXElem create(const Function& f, const std::vector<SXElem>& dep) {
      casadi_assert(f.nnz_in()==dep.size(),
        "CallSX::create(f,dep): dimension mismatch: " + str(f.nnz_in()) + " vs " + str(dep.size()));
      return SXElem::create(new CallSX(f, dep));
    }

    inline static int eval_sx(const Function& f, const SXElem** arg, SXElem** res) {
      // Collect dep from arg
      std::vector<SXElem> dep;
      for (casadi_int i=0;i<f.n_in();++i) {
        dep.insert(dep.end(), arg[i], arg[i]+f.nnz_in(i));
      }

      std::vector<SXElem> ret = SXElem::call(f, dep);

      // Copy ret to res
      casadi_int offset = 0;
      for (casadi_int i=0;i<f.n_out();++i) {
        casadi_copy(get_ptr(ret)+offset, f.nnz_out(i), res[i]);
        offset += f.nnz_out(i);
      }

      return 0;
    }

    /** \brief Destructor

    This is a rather complex destructor which is necessary since the default destructor
    can cause stack overflow due to recursive calling.

        \identifier{28q} */
    ~CallSX() override {
      for (auto & d  : dep_)
        safe_delete(d.assignNoDelete(casadi_limits<SXElem>::nan));
    }

    // Class name
    std::string class_name() const override {return "CallSX";}

    bool is_op(casadi_int op) const override { return OP_CALL==op; }

    /** \brief  Number of dependencies

        \identifier{28r} */
    casadi_int n_dep() const override { return dep_.size();}

    /** \brief  Get an output

        \identifier{28s} */
    SXElem get_output(casadi_int oind) const override {
      SharedSXElem ret;
      // If not found, add it,
      if (!cache_.incache(oind, ret)) {
        ret.own(new OutputSX(SXNode::shared_from_this(), oind));
        cache_.tocache_if_missing(oind, ret);
      }
      return SXElem::create(ret.get());
    }

    /** \brief  get the reference of a dependency

        \identifier{28t} */
    const SXElem& dep(casadi_int i) const override { return dep_[i];}
    SXElem& dep(casadi_int i) override { return dep_[i];}

    /** \brief  Get the operation

        \identifier{28u} */
    casadi_int op() const override { return OP_CALL; }

    /** \brief Check if call node

        \identifier{28v} */
    bool is_call() const override { return true; }

    /** \brief  Check if a multiple output node

        \identifier{28w} */
    bool has_output() const override { return true; }

    /** \brief  Get function

        \identifier{28x} */
    Function which_function() const override { return f_; }

    /** \brief  Print expression

        \identifier{28y} */
    std::string print(const std::string& arg1, const std::string& arg2) const override {
      return f_.name();
    }

    Function f_;

    /** \brief  The dependencies of the node

        \identifier{28z} */
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
/*#ifndef CALL_SXElem_HPP
#define CALL_SXElem_HPP

#include "function.hpp"
#include "sx_node.hpp"

/// \cond INTERNAL

namespace casadi {

class CodeGenerator;
struct Instance;

class CallSX : public SXNode {
public:
  explicit CallSX(const SXElem& ref, const SXElem& arg);
  ~CallSX() override {}

  /** \brief  Get the operation
  casadi_int op() const override { return OP_CALL;}

  bool is_op(casadi_int op) const override { return op==OP_CALL; }

  // Class name
  std::string class_name() const override {return "CallSX";}

  /** \brief  Print expression
  std::string print(const std::string& arg1, const std::string& arg2) const override;

  /** \brief  Number of dependencies
  casadi_int n_dep() const override { return 2;}

  /** \brief  get the reference of a dependency
  const SXElem& dep(casadi_int i) const override { return i==0 ? ref_ : arg_ ; }
  SXElem& dep(casadi_int i) override { return i==0 ? ref_ : arg_ ; }

  static void der(const SXElem& x, const SXElem& y, const SXElem& f, SXElem* d);

  static std::string codegen(CodeGenerator& g, const SXElem& funref, const Instance& inst, int i0, int i1, int i2, const std::string& arg, const std::string& res, const std::string& iw, const std::string& w, const Function& owner);

  static void codegen_dependency(CodeGenerator& g, const Function& f, const Instance& inst, const Function& owner);

  void serialize_node(SerializingStream& s) const override;

  /** \brief Deserialize without type information
  static SXNode* deserialize(DeserializingStream& s);

  static int call(const Function &f, double& result, double dep, double index, const double** arg, double** res, casadi_int* iw, double* w);

  SXElem ref_;
  SXElem arg_;
};

} // namespace casadi
/// \endcond
#endif // CALL_SXElem_HPP*/
