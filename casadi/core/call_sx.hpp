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


#ifndef CALL_SXElem_HPP
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

  /** \brief  Get the operation */
  casadi_int op() const override { return OP_CALL;}

  bool is_op(casadi_int op) const override { return op==OP_CALL; }

  // Class name
  std::string class_name() const override {return "CallSX";}

  /** \brief  Print expression */
  std::string print(const std::string& arg1, const std::string& arg2) const override;

  /** \brief  Number of dependencies */
  casadi_int n_dep() const override { return 2;}

  /** \brief  get the reference of a dependency */
  const SXElem& dep(casadi_int i) const override { return i==0 ? ref_ : arg_ ; }
  SXElem& dep(casadi_int i) override { return i==0 ? ref_ : arg_ ; }

  static void der(const SXElem& x, const SXElem& y, const SXElem& f, SXElem* d);

  static std::string codegen(CodeGenerator& g, const SXElem& funref, const Instance& inst, int i0, int i1, int i2, const std::string& arg, const std::string& res, const std::string& iw, const std::string& w, const Function& owner);

  static void codegen_dependency(CodeGenerator& g, const Function& f, const Instance& inst, const Function& owner);

  void serialize_node(SerializingStream& s) const override;

  /** \brief Deserialize without type information */
  static SXNode* deserialize(DeserializingStream& s);

  static int call(const Function &f, double& result, double dep, double index, const double** arg, double** res, casadi_int* iw, double* w);

  SXElem ref_;
  SXElem arg_;
};

} // namespace casadi
/// \endcond
#endif // CALL_SXElem_HPP
