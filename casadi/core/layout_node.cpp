/*
 *  This file is part of CasADi.
 *
 *  CasADi -- A symbolic framework for dynamic optimization.
 *  Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *              K.U. Leuven. All rights reserved.
 *  Copyright (C) 2011-2014 Greg Horn
 *
 *  CasADi is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  CasADi is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with CasADi; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "layout_node.hpp"
#include "layout.hpp"
#include "casadi_misc.hpp"
#include "runtime/casadi_runtime.hpp"
#include "serializing_stream.hpp"

using namespace std;

namespace casadi {


  StridedLayout::StridedLayout(const std::vector<casadi_int>& dims)
        : StridedLayout(dims, dims) {

  }

  StridedLayout::StridedLayout(const std::vector<casadi_int>& dims,
         const std::vector<casadi_int>& outer_sizes) {
    layout_.resize(3+dims.size()*3);
    layout_[0] = 1;
    layout_[1] = dims.size();

    casadi_int nnz = 1;
    for (casadi_int i=0;i<dims.size();++i) nnz*= dims[i];
    layout_[2] = nnz;

    casadi_assert(!has_negative(dims), "Dimensions must be positive");
    casadi_assert(dims.size()==outer_sizes.size(), "Dimension mismatch");
    for (casadi_int i=0;i<dims.size();++i) {
      casadi_assert(outer_sizes[i]>= dims[i], "Error in outer_sizes");
    }

    std::copy(dims.begin(), dims.end(), this->dims());
    std::copy(outer_sizes.begin(), outer_sizes.end(), this->outer_sizes());

    uout() << "fool outer_sizes" << outer_sizes << std::endl;

    strides()[0] = 1;
    for (casadi_int j=0;j<dims.size()-1;++j) {
      strides()[j+1] = strides()[j]*outer_sizes[j];
    }
    uout() << "fool" << get_strides() << std::endl;
  }

  /** \brief Get name of public class */
  std::string StridedLayout::class_name() const { return "StridedLayout"; }

  /** \brief  Print a description */
  void StridedLayout::disp(std::ostream& stream, bool more) const {
    stream << "strided_layout(" << "dims" << get_dims() << ", outer_size" << get_outer_sizes() << ", strides" << get_strides() << ")";
  }


  StridedLayout::~StridedLayout() {

  }

  Layout StridedLayout::push_right(casadi_int n) const {
    std::vector<casadi_int> dims = join(get_dims(), {n});
    std::vector<casadi_int> outer_sizes = join(get_outer_sizes(), {n});
    return Layout(dims, outer_sizes);
  }

  Layout StridedLayout::stride(casadi_int n) const {
    casadi_error("foo");
    /*std::vector<casadi_int> strides_source = get_strides_source();
    strides_source[0] *= n;
    return Layout(get_dims(), get_perms(), strides_source);*/
  }

  casadi_int StridedLayout::stride() const {
    return strides()[0];
  }

  casadi_int StridedLayout::stride_right() const {
    return strides()[n_dims()-1];
  }

  size_t StridedLayout::size() const {
    uout() << "size query" << get_strides() << get_dims() << ":" << strides()[n_dims()-1]*dims()[n_dims()-1] << std::endl;
    return strides()[n_dims()-1]*outer_sizes()[n_dims()-1];
  }

  size_t StridedLayout::n_dims() const {
    return layout_[1];
  }

  size_t StridedLayout::nnz() const {
    return layout_[2];
  }

  bool StridedLayout::operator==(const Layout& other) const {
    if (other.is_default()) return false;
    const StridedLayout* rhs = static_cast<const StridedLayout*>(other.get());
    return get_dims()==rhs->get_dims() && get_strides()==rhs->get_strides();
  }

  void StridedLayout::assert_valid_permutation(const Layout& target) const {
    casadi_assert(nnz()==target.nnz(), "permutation cannot change total amount of data");

    const StridedLayout* t = static_cast<const StridedLayout*>(target.get());

  }

  LayoutNode::operator const casadi_int*() const {
    return get_ptr(layout_);
  }

  DefaultLayout::DefaultLayout() {
    layout_ = {0};
  }

  /** \brief Get name of public class */
  std::string DefaultLayout::class_name() const { return "DefaultLayout"; }

  /** \brief  Print a description */
  void DefaultLayout::disp(std::ostream& stream, bool more) const {
    stream << "default_layout";;
  }


  DefaultLayout::~DefaultLayout() {

  }

  size_t DefaultLayout::size() const {
    casadi_error("Unknown d");
  }
  size_t DefaultLayout::nnz() const {
    casadi_error("Unknown c");
  }
  size_t DefaultLayout::n_dims() const {
    casadi_error("Unknown b");
  }

  Layout DefaultLayout::push_right(casadi_int n) const {
    return shared_from_this<Layout>();
  }

  Layout DefaultLayout::stride(casadi_int n) const {
    return shared_from_this<Layout>();
  }

  casadi_int DefaultLayout::stride() const {
    return 1;
  }

  casadi_int DefaultLayout::stride_right() const {
    return 1;
  }

  bool DefaultLayout::operator==(const Layout& other) const {
    return other.is_default();
  }

  void LayoutNode::serialize(SerializingStream& s) const {
    s.pack("LayoutNode::type", layout_[0]);
    serialize_body(s);
  }

  void LayoutNode::serialize_body(SerializingStream& s) const {
    s.pack("LayoutNode::layout", layout_);
  }

  LayoutNode* LayoutNode::deserialize(DeserializingStream& s) {
    casadi_int type;
    s.unpack("LayoutNode::type", type);
    if (type==0) {
      return new DefaultLayout(s);
    } else if (type==1) {
      return new StridedLayout(s);
    }
  }

  LayoutNode::LayoutNode(DeserializingStream& s) {
    s.unpack("LayoutNode::layout", layout_);
  }

  Layout Layout::create(LayoutNode* node) {
    return Layout(node);
  }

  Layout::Layout(LayoutNode* node) {
    own(node);
  }

} // namespace casadi
