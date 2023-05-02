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


#include "io_instruction.hpp"
#include "serializing_stream.hpp"

namespace casadi {
  Input::Input(const Sparsity& sp, casadi_int ind, casadi_int segment, casadi_int offset)
    : IOInstruction(ind, segment, offset) {
    set_sparsity(sp);
  }

  std::string Input::disp(const std::vector<std::string>& arg) const {
    std::stringstream s;
    s << "input[" << ind_ << "][" << segment_ << "]";
    return s.str();
  }

  void Input::generate(CodeGenerator& g,
      const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    casadi_int nnz = this->nnz();
    if (nnz==0) return; // quick return
    std::string a = g.arg(ind_);
    casadi_int i = res.front();
    if (nnz==1) {
      g << g.workel(i) << " = " << a << " ? " << a << "[" << offset_ << "] : 0;\n";
    } else if (offset_==0) {
      g << g.copy(a, nnz, g.work(i, nnz)) << "\n";
    } else {
      g << g.copy(a + " ? " + a + "+" + str(offset_) + " : 0",
                          nnz, g.work(i, nnz)) << "\n";
    }
  }

  Output::Output(const MX& x, casadi_int ind, casadi_int segment, casadi_int offset)
    : IOInstruction(ind, segment, offset) {
    set_dep(x);
  }

  std::string Output::disp(const std::vector<std::string>& arg) const {
    std::stringstream s;
    s << "output[" << ind_ << "][" << segment_ << "]";
    return s.str();
  }

  void Output::generate(CodeGenerator& g,
      const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    casadi_int nnz = dep().nnz();
    if (nnz==0) return; // quick return
    casadi_int i = arg.front();
    std::string r = g.res(ind_);
    if (nnz==1) {
      g << "if (" << r << ") " << r << "[" << offset_ << "] = " << g.workel(i) << ";\n";
    } else if (offset_==0) {
      g << g.copy(g.work(i, nnz), nnz, r) << "\n";
    } else {
      g << "if (" << r << ") "
        << g.copy(g.work(i, nnz), nnz, r + "+" + str(offset_)) << "\n";
    }

  }

  Dict IOInstruction::info() const {
    Dict ret;
    ret["ind"] = ind_;
    ret["segment"] = segment_;
    ret["offset"] = offset_;
    return ret;
  }

  void IOInstruction::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("IOInstruction::ind", ind_);
    s.pack("IOInstruction::segment", segment_);
    s.pack("IOInstruction::offset", offset_);
  }

  IOInstruction::IOInstruction(DeserializingStream& s) : MXNode(s) {
    s.unpack("IOInstruction::ind", ind_);
    s.unpack("IOInstruction::segment", segment_);
    s.unpack("IOInstruction::offset", offset_);
  }

} // namespace casadi
