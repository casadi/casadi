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


#ifndef CASADI_LAYOUT_NODE_HPP
#define CASADI_LAYOUT_NODE_HPP

#include "shared_object_internal.hpp"
#include "shared_object_internal.hpp"

namespace casadi {

  class Layout;
  class SerializingStream;
  class DeserializingStream;

  /** \brief memory layout class
      \author Joris Gillis
      \date 2021
      Internal class.
  */
  class CASADI_EXPORT LayoutNode : public SharedObjectInternal {
    friend class Layout;

  public:
    LayoutNode() = default;
    virtual Layout push_right(casadi_int n) const=0;
    virtual casadi_int stride() const=0;
    virtual casadi_int stride_right() const=0;
    virtual Layout stride(casadi_int n) const=0;
    virtual bool operator==(const Layout& other) const=0;
    virtual size_t size() const=0;
    virtual size_t nnz() const=0;
    virtual size_t n_dims() const=0;
    virtual void assert_valid_permutation(const Layout& target) const {};

    /// Is storage size greater than number of nonzeros?
    virtual bool has_padding() const = 0;
    virtual bool is_default() const { return false; }
    std::vector<casadi_int> get_compressed() const { return layout_; }
    operator const casadi_int*() const;
    std::vector<casadi_int> layout_;

    /** \brief Serialize an object */
    void serialize(SerializingStream& s) const;

    /** \brief Serialize an object without type information */
    virtual void serialize_body(SerializingStream& s) const;

    /** \brief Deserialize with type disambiguation
     *
     * Uses the information encoded with serialize_type to unambiguously find the (lowest level sub)class,
     * and calls its deserializing constructor.
     */
    static LayoutNode* deserialize(DeserializingStream& s);

  protected:
    /** \brief Deserializing constructor */
    explicit LayoutNode(DeserializingStream& s);
  };

  /** \brief memory layout class
      \author Joris Gillis
      \date 2021
      Internal class.
  */
  class CASADI_EXPORT DefaultLayout : public LayoutNode {

  public:
    /// Constructor
    DefaultLayout();

    /** \brief  Destructor */
    ~DefaultLayout() override;

    /** \brief Get name of public class */
    std::string class_name() const override;

    /** \brief  Print a description */
    void disp(std::ostream& stream, bool more) const override;

    Layout push_right(casadi_int n) const override;
    casadi_int stride() const override;
    casadi_int stride_right() const override;
    Layout stride(casadi_int n) const override;
    bool operator==(const Layout& other) const override;
    size_t size() const override;
    size_t nnz() const override;
    size_t n_dims() const override;

    bool is_default() const override { return true; }
    bool has_padding() const override { return false; }

    /** \brief Deserializing constructor */
    explicit DefaultLayout(DeserializingStream& s) : LayoutNode(s) {}
  };

  /** \brief memory layout class
   * 
   * 
   * 
   * A layout tells you how to interpret the (flat) data
   * 
   * e.g. data = [0,1,2,3,4,5]
   * 
   * dims  -> draw that matrix
   * 
   *                 _axis1__  
   *        axis 0  |       |
   *                |_______|
   * perms -> while going through data, perms tells you what axes to fill in first: first axis perms[0], perms[1], ...
   *           you should takes steps strides[0] when jumping along axis 0, strides[1] when jumping along axis 1, etc
   * 
   *             e.g. perms [0,1] means column-major, column by column left to right
   *                        [1,0] means row-major, row by row top to bottom
   * 
   * Example A (unittested)
   * 
   * Layout([2,3],[0,1])  interpreted as 2-by-3 [0 2 4
   *                                             1 3 5]
   *
   * Layout([3,2],[0,1])  interpreted as 3-by-2 [0 3
   *                                             1 4 
   *                                             2 5]
   * 
   * Layout([2,3],[1,0])  interpreted as 2-by-3 [0 1 2
   *                                             3 4 5]
   * 
   * Layout([3,2],[1,0])  interpreted as 3-by-2 [0 1
   *                                             2 3 
   *                                             4 5]
   * Interchangeable valid viewpoints:
   *  - <dims> tensor, but peculiar traversal order when flattening
   *  - a linear slab of memory, but you can interpret me as a <dims> tensor if you observe this recipe
   * 
   * Invalid viewpoints:
   *   <dims[perms]> tensor, natural ordering
   *   Example A contradicts this interpretation:  A2 and A3 have same dims[perms] but different outcome
   *       Note that even the strides of A2 and A3 are the same: both [1,3]
   * 
   * Examples:
   * 
   * An MX annotated with Layout([6 4 5],[1,0,2])
   *   x/6-by-y/4-by-z/5 tensor, but to flatten, first go along y axis, then x axis, tehn z axis
   * 
   * SX Function with an input strided_layout(dims[6], perms[0], strides_source[4], strides[4]):
   *   when you pass an array w, the vector [ w[0],w[4],w[8],w[12],w[16],w[20] ] will enter the virtual machine
   * 
   * MX Function with an input strided_layout(dims[6, 4], perms[1, 0], strides_source[1, 4], strides[1, 4]):
   *   it's a 6-by-4 matrix, but stored row-major
   * 
   * 
   * 
   * strides_source default: [1 dims[perms[0]], dims[perms[1]], ... ]  (excludes the last dimensions)
   *    actual strides_source must be greater than or equal to the default, element-wise
   * 
   * \author Joris Gillis
   * \date 2021
   * Internal class.
  */
  class CASADI_EXPORT StridedLayout : public LayoutNode {
  public:
    /// Constructor
    StridedLayout(const std::vector<casadi_int>& dims, const std::vector<casadi_int>& outer_sizes);
    /// Constructor
    StridedLayout(const std::vector<casadi_int>& dims);

    Layout push_right(casadi_int n) const override;
    casadi_int stride() const override;
    casadi_int stride_right() const override;
    Layout stride(casadi_int n) const override;
    size_t size() const override;
    size_t nnz() const override;
    size_t n_dims() const override;
    void assert_valid_permutation(const Layout& target) const override;

    bool operator==(const Layout& other) const override;

    /** \brief  Destructor */
    ~StridedLayout() override;

    /** \brief Get name of public class */
    std::string class_name() const override;

    /** \brief  Print a description */
    void disp(std::ostream& stream, bool more) const override;

    /** \brief dimensions */
    inline const casadi_int* dims() const { return &layout_.front()+3; }
    inline casadi_int* dims() { return &layout_.front()+3; }
    inline std::vector<casadi_int> get_dims() const {
      return std::vector<casadi_int>(dims(), dims()+n_dims());
    }

    /** \brief Order of traversal */
    inline const casadi_int* outer_sizes() const { return dims()+n_dims(); }
    inline casadi_int* outer_sizes() { return dims()+n_dims(); }
    inline std::vector<casadi_int> get_outer_sizes() const {
      return std::vector<casadi_int>(outer_sizes(), outer_sizes()+n_dims());
    }

    /** \brief Order of traversal */
    inline const casadi_int* strides() const { return outer_sizes()+n_dims(); }
    inline casadi_int* strides() { return outer_sizes()+n_dims(); }
    inline std::vector<casadi_int> get_strides() const {
      return std::vector<casadi_int>(strides(), strides()+n_dims());
    }

    bool has_padding() const override { return nnz()!=size(); }

    /** \brief Deserializing constructor */
    explicit StridedLayout(DeserializingStream& s) : LayoutNode(s) {}
  };

  /// \endcond
} // namespace casadi

#endif // CASADI_LAYOUT_NODE_HPP
