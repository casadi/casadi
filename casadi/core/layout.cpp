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


#include "layout_node.hpp"
#include "layout.hpp"
#include "runtime/casadi_runtime.hpp"
#include "casadi_misc.hpp"
#include "code_generator.hpp"

using namespace std;
namespace casadi {

  Layout::Layout() {
    own(new DefaultLayout());
  }
  /*Layout::Layout(const Sparsity& sp) {
    own(new StridedLayout({sp.nnz()}));
  }*/

  Layout::Layout(const std::vector<casadi_int>& dims, const std::vector<casadi_int>& outer_sizes) {
    own(new StridedLayout(dims, outer_sizes));
  }

  Layout::Layout(const std::vector<casadi_int>& dims) {
    own(new StridedLayout(dims));
  }

  Layout::~Layout() {

  }

  Relayout Relayout::push_right(casadi_int n) const {
    return Relayout(source().push_right(n), join(perms(), {source().n_dims()}), target().push_right(n));
  }

  void Layout::assert_valid_permutation(const Layout& target) const {
    return (*this)->assert_valid_permutation(target); 
  }

  Layout Layout::push_right(casadi_int n) const {
    return (*this)->push_right(n);
  }

  casadi_int Layout::stride() const {
    return (*this)->stride();
  }

  Layout Layout::stride(casadi_int n) const {
    return (*this)->stride(n);
  }

  casadi_int Layout::stride_right() const {
    return (*this)->stride_right();
  }

  bool Layout::operator==(const Layout& other) const {
    casadi_assert(!is_null(), "lhs is null");
    casadi_assert(!other.is_null(), "rhs is null");
    return (*this)->operator==(other);
  }

  LayoutNode* Layout::operator->() {
    return static_cast<LayoutNode*>(SharedObject::operator->());
  }

  const LayoutNode* Layout::operator->() const {
    return static_cast<const LayoutNode*>(SharedObject::operator->());
  }

  Layout::operator const casadi_int*() const {
    return *((*this).operator->());
  }

  bool Layout::test_cast(const SharedObjectInternal* ptr) {
    return dynamic_cast<const LayoutNode*>(ptr)!=nullptr;
  }

  size_t Layout::size() const {
    return (*this)->size();
  }

  size_t Layout::nnz() const {
    return (*this)->nnz();
  }

  bool Layout::has_padding() const {
    return (*this)->has_padding();
  }

  bool Layout::is_default() const {
    return (*this)->is_default();
  }

  size_t Layout::n_dims() const {
    return (*this)->n_dims();
  }

  std::vector<casadi_int> Layout::get_compressed() const {
    return (*this)->get_compressed();
  }


  void Layout::serialize(SerializingStream& s) const {
    return (*this)->serialize(s);
  }

  Layout Layout::deserialize(DeserializingStream& s) {
    return Layout::create(LayoutNode::deserialize(s));
  }

  Relayout::Relayout(const Layout& source, const std::vector<casadi_int>& perms, const Layout& target, const std::string& label) : source_(source), perms_(perms), target_(target), label_(label) {
    const StridedLayout* s = static_cast<const StridedLayout*>(source.get());
    const StridedLayout* t = static_cast<const StridedLayout*>(target.get());

    casadi_assert(is_permutation(perms), "Not a valid permutation");
    casadi_assert(perms.size()==s->n_dims(), "Not a valid permutation");
    casadi_assert(permute(s->get_dims(), perms)==t->get_dims(), "Inconsistent dimensions");
  }

  Relayout::Relayout() : source_(Layout()), target_(Layout()) {

  }

  Relayout Relayout::invert() const {
    return Relayout(target(), invert_permutation(perms_), source());
  }


  size_t Relayout::sz_iw() const {
    const StridedLayout* s = static_cast<const StridedLayout*>(source_.get());
    return s->n_dims()*2;
  }

  bool Relayout::absorbs(const Relayout& other) const {
    if (target()==other.source()) return true;
  }

  Relayout Relayout::absorbed(const Relayout& other) const {
    casadi_assert_dev(target()==other.source());

    std::vector<casadi_int> perm = permute(perms(), other.perms());

    return Relayout(source(), perm, other.target());
  }

  std::vector<casadi_int> Relayout::nz_ref() const {
    
  }

  bool Relayout::cancels(const Relayout& other) const {
    uout() << "cancels?" << source_ << target_ << std::endl;
    uout() << "cancels?" <<  other.source() << other.target() << std::endl;
    //if (source_==other.target() && target_==other.source() && invert_permutation(perms_)==other.perms()) return true;
    if (source_.size()!=other.target().size()) return false;
    size_t sz = max(sz_iw(), other.sz_iw());
    std::vector<casadi_int> iw(sz);

    size_t size = max(max(max(source_.size(), target_.size()), other.source().size()), other.target().size());


    const StridedLayout* s = static_cast<const StridedLayout*>(source_.get());

    std::vector<casadi_int> arg = range(size);
    std::vector<casadi_int> res(size, -1);
    Layout unpadded_source = Layout(s->get_dims());
    std::vector<casadi_int> perms = range(s->n_dims());
    const casadi_int* source = unpadded_source;
    const casadi_int* target = source_;
    casadi_relayout(get_ptr(arg), get_ptr(res), source, get_ptr(perms), target,  get_ptr(iw));
    std::vector<casadi_int> ref = res;
    uout() << "a:" << ref << std::endl;
    std::copy(res.begin(), res.end(), arg.begin());

    source = source_;
    target = target_;
    uout() << "bef:" << arg << std::endl;
    casadi_relayout(get_ptr(arg), get_ptr(res), source, get_ptr(perms_), target,  get_ptr(iw));
    std::copy(res.begin(), res.end(), arg.begin());
    std::fill(res.begin(), res.end(), -1);
    source = other.source();
    target = other.target();
    casadi_relayout(get_ptr(arg), get_ptr(res), source, get_ptr(other.perms()), target, get_ptr(iw));
    uout() << "arg:" << arg << std::endl;
    uout() << "res:" << res << std::endl;
    for (casadi_int i=0;i<source_.size();++i) {
      if (ref[i]!=res[i]) {
        uout() << "not cancelled" << std::endl;
        return false;
      }
    }
    uout() << "cancelled" << std::endl;
    return true;
  }


  void Relayout::generate(CodeGenerator& g,
                        const std::string& arg,
                        const std::string& res) const {

    const StridedLayout* s = static_cast<const StridedLayout*>(source_.get());
    const StridedLayout* t = static_cast<const StridedLayout*>(target_.get());

    casadi_int n_dims = s->n_dims();
    
    std::vector<casadi_int> target_strides_perm(n_dims);
    for (casadi_int j=0;j<n_dims;++j) {
      target_strides_perm[perms_[j]] = t->strides()[j];
    }

    std::vector<casadi_int> strides = s->get_strides();
    std::vector<casadi_int> dummy, order;
    sort(strides, dummy, order);

    g.comment(str(order));

    for (casadi_int i=0;i<n_dims;++i) {
      if (s->dims()[i]>1) {
        g.local("i"+str(i), "casadi_int");
        g << "for (i" << i << "=0;i" << i << "<" << s->dims()[i] << ";++i" << i << ") ";
      }
    }

    std::vector<std::string> k;
    for (casadi_int i=0;i<n_dims;++i) {
      if (s->dims()[i]>1) k.push_back(str(s->strides()[i]) + "*i"+str(i));
    }
    std::vector<std::string> kk;
    for (casadi_int i=0;i<n_dims;++i) {
      if (s->dims()[i]>1) kk.push_back(str(target_strides_perm[i]) + "*i"+str(i));
    }
    g << res << "[" << join(kk, "+") << "] = " << arg << "[" << join(k, "+") << "];\n";

  }


  bool Relayout::is_trivial() const {
    if (source_.size()!=target_.size()) return false;
    std::vector<casadi_int> iw(sz_iw());
    size_t size = source_.size();

    const StridedLayout* s = static_cast<const StridedLayout*>(source_.get());
    std::vector<casadi_int> arg = range(size);
    std::vector<casadi_int> res(size, -1);
    Layout unpadded_source = Layout(s->get_dims());
    std::vector<casadi_int> perms = range(s->n_dims());
    const casadi_int* source = unpadded_source;
    const casadi_int* target = source_;
    casadi_relayout(get_ptr(arg), get_ptr(res), source, get_ptr(perms), target,  get_ptr(iw));
    std::vector<casadi_int> ref = res;
    uout() << "a:" << ref << std::endl;
    std::copy(res.begin(), res.end(), arg.begin());
    std::fill(res.begin(), res.end(), -1);
    source = source_;
    target = target_;
    uout() << "bef:" << arg << std::endl;
    casadi_relayout(get_ptr(arg), get_ptr(res), source, get_ptr(perms_), target,  get_ptr(iw));
    uout() << "ref:" << ref << std::endl;
    uout() << "res:" << res << std::endl;
    for (casadi_int i=0;i<source_.size();++i) {
      if (ref[i]!=res[i]) {
        uout() << "not cancelled" << std::endl;
        return false;
      }
    }
    uout() << "cancelled" << std::endl;
    return true;
  }

} // namespace casadi
