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


#ifndef CASADI_LAYOUT_HPP
#define CASADI_LAYOUT_HPP
#include "shared_object.hpp"
#include "printable.hpp"
#include "sparsity.hpp"
#include <vector>

namespace casadi {

  /** \brief  Forward declaration */
  class LayoutNode;



  /** \brief Memory Layout class

      \author Joris Gillis
      \date 2021
  */
  class CASADI_EXPORT Layout :
    public SWIG_IF_ELSE(PrintableCommon, Printable<Layout>),
    public SharedObject {
  public:
    Layout();
    Layout(LayoutNode* node);
    //Layout(const Sparsity& sp);
    Layout(const std::vector<casadi_int>& dims, const std::vector<casadi_int>& outer_sizes);
    Layout(const std::vector<casadi_int>& dims);
    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    Layout push_right(casadi_int n) const;
    casadi_int stride() const;
    casadi_int stride_right() const;
    Layout stride(casadi_int n) const;
    std::vector<casadi_int> get_compressed() const;
#ifndef SWIG
    operator const casadi_int*() const;
#endif // SWIG
    void assert_valid_permutation(const Layout& target) const;

    /// Is storage size greater than number of nonzeros?
    bool has_padding() const;
    bool is_default() const;

    /// Storage size needed for layout (includes padding specified in strides)
    size_t size() const;
    /// Number of nonzeros in the tensor (padding not included)
    size_t nnz() const;

    size_t n_dims() const;

    /** \brief Serialize an object */
    void serialize(SerializingStream& s) const;

    /** \brief Deserialize with type disambiguation */
    static Layout deserialize(DeserializingStream& s);

    static Layout create(LayoutNode* node);

    ~Layout();

#ifndef SWIG
    /// \cond INTERNAL
    ///@{
    /** \brief  Access a member of the node */
    LayoutNode* operator->();

    /** \brief  Const access a member of the node */
    const LayoutNode* operator->() const;
    ///@}
    /// \endcond
#endif // SWIG

    bool operator==(const Layout& other) const;

  };

  class CASADI_EXPORT Relayout {

  public:
    Relayout();
    Relayout(const Layout& source, const std::vector<casadi_int>& perms, const Layout& target, const std::string& label="");

    bool is_trivial() const;

    bool cancels(const Relayout& other) const;
    bool absorbs(const Relayout& other) const;
    Relayout absorbed(const Relayout& other) const;

    Relayout invert() const;
    /// Print a string representation of the object to a stream
    inline friend
      std::ostream& operator<<(std::ostream &stream, const Relayout& obj) {
      stream << "(" << obj.label_ << ":" << obj.source() << "->" << obj.target() << ")";
      return stream;
    }
    Relayout push_right(casadi_int n) const;
    inline const Layout& source() const { return source_; }
    inline const Layout& target() const { return target_; }
    inline const std::vector<casadi_int>& perms() const { return perms_; }
    size_t sz_iw() const;

    std::vector<casadi_int> nz_ref() const;

    void generate(CodeGenerator& g,
                  const std::string& arg,
                  const std::string& res) const;

    Layout source_;
    Layout target_;
    std::vector<casadi_int> perms_;
    std::string label_;

  };

} // namespace casadi

#endif // CASADI_LAYOUT_HPP
