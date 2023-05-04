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


#include "multiple_output.hpp"
#include "function_internal.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"

namespace casadi {

  MultipleOutput::MultipleOutput() {
  }

  MultipleOutput::~MultipleOutput() {
  }

  MX MultipleOutput::get_output(casadi_int oind) const {
    MX this_ = shared_from_this<MX>();
    // No need for an OutputNode if sparsity is fully sparse
    if (this_->sparsity(oind).nnz()==0) return MX(this_->sparsity(oind));
    return MX::create(new OutputNode(this_, oind));
  }

  OutputNode::OutputNode(const MX& parent, casadi_int oind) : oind_(oind) {
    set_dep(parent);

    // Save the sparsity pattern
    set_sparsity(dep(0)->sparsity(oind));
  }

  OutputNode::~OutputNode() {
  }

  std::string OutputNode::disp(const std::vector<std::string>& arg) const {
    return arg.at(0) + "{" + str(oind_) + "}";
  }

  void OutputNode::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("OutputNode::oind", oind_);
  }

  OutputNode::OutputNode(DeserializingStream& s) : MXNode(s) {
    s.unpack("OutputNode::oind", oind_);
  }

} // namespace casadi
